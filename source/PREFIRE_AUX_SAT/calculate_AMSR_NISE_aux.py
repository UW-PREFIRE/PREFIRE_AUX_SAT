import argparse
from datetime import datetime
import os
import netCDF4 as n
import xarray as xr
import pyhdf.SD
import numpy as np
import esmpy
import pyproj
import re
import json
import warnings

from PREFIRE_tools.utils.aux import construct_tirs_geo


FillValue = -9999.0

hemis = ['nh','sh']

amsr_grid_names = {
    'nh':'NpPolarGrid12km',
    'sh':'SpPolarGrid12km'
    }
amsr_sic_varnames = {
    'nh':'SI_12km_NH_ICECON_DAY',
    'sh':'SI_12km_SH_ICECON_DAY'
    }

nise_lat0 = {
    'nh':'90',
    'sh':'-90'
    }

amsr_product_name = 'AMSR_U2_L3_SeaIce12km_B04'
nise_product_name = 'NISE_SSMISF18'

amsr_aux_vars = [
    'AMSR_L3_seaice_concentration'
    ]

nise_aux_vars = [
    'NISE_L3_snow_fraction',
    'NISE_L3_seaice_concentration'
    ]


def compile_file_lookup(input_products_dir, tirs_times):
    """
    Create dictionary containing the file name of each daily L3 AMSR and NISE
    product file during the time range covered by tirs_times.

    Parameters
    ----------
    input_products_dir : str
        Path to parent directory of input product files used for Aux-Sat.
    tirs_times : list
        Time of each TIRS scene (type: datetime.datetime).

    Returns
    -------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.

    """
    
    tirs_tstrs = [t.strftime('%Y%m%d') for t in tirs_times]
        
    file_lookup = {
        f'{amsr_product_name}':{},
        f'{nise_product_name}':{}
        }
    
    for tstr in set(tirs_tstrs):
        # Check for existence of AMSR product file for the date
        amsr_fpath = os.path.join(input_products_dir, amsr_product_name,
                                  f"{amsr_product_name}_{tstr}.he5")
        if os.path.exists(amsr_fpath):
            file_lookup[f'{amsr_product_name}'][tstr] = amsr_fpath
        else:
            amsr_fname = os.path.basename(amsr_fpath)
            warnings.warn(f'Missing {amsr_fpath}')

        # Check for existence of NISE product file for the date
        nise_fpath = os.path.join(input_products_dir, nise_product_name,
                                  f"{nise_product_name}_{tstr}.HDFEOS")
        if os.path.exists(nise_fpath):
            file_lookup[f'{nise_product_name}'][tstr] = nise_fpath
        else:
            nise_fname = os.path.basename(nise_fpath)
            warnings.warn(f'Missing {nise_fpath}')
            
    return file_lookup


def get_AMSR_NISE_attrs(file_lookup):
    """
    Get AMSR and NISE product metadata to be assigned to Aux-Sat group in output
    file.

    Parameters
    ----------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.

    Returns
    -------
    amsr_nise_attrs : dict
        AMSR and NISE product metadata extracted from product files.

    """
    
    # AMSR: get maturity code, version, and DOI
    # - maturity code is hard-coded into AMSR product name in this script for now
    try:
        amsr_ex_fpath = \
            file_lookup[f'{amsr_product_name}'][list(file_lookup[f'{amsr_product_name}'].keys())[0]]
        with n.Dataset(amsr_ex_fpath) as amsr_nc:
            amsr_version = amsr_nc.getncattr('references').split('Version ')[1].split('.')[0]
            amsr_doi = amsr_nc.variables['DOI'][:]
        amsr_mat_code = amsr_product_name.split('AMSR_U2_L3_SeaIce12km_')[1]
    # Handle cases where no AMSR product files are available
    except:
        amsr_mat_code = 'unavailable'
        amsr_version = 'unavailable'
        amsr_doi = 'unavailable'
        
        
    # NISE: get version and PGE version
    try:
        nise_ex_fpath = \
            file_lookup[f'{nise_product_name}'][list(file_lookup[f'{nise_product_name}'].keys())[0]]
        nise_hdf = pyhdf.SD.SD(nise_ex_fpath, pyhdf.SD.SDC.READ)
        nise_meta = nise_hdf.attributes(full=1)['coremetadata.0'][0].replace(" ", "")
        nise_version = nise_meta.\
            split('OBJECT=VERSIONID')[1].split('END_OBJECT')[0].\
            split('VALUE=')[1].split('\n')[0]
        nise_pge_version = nise_meta.\
            split('OBJECT=PGEVERSION')[1].split('END_OBJECT')[0].\
            split('VALUE="v')[1].split('"')[0]
    # Handle cases where no NISE product files are available
    except:
        nise_version = 'unavailable'
        nise_pge_version = 'unavailable'
        
    amsr_nise_attrs = {
        'AMSR_L3_seaice_maturity_code': amsr_mat_code,
        'AMSR_L3_seaice_version': amsr_version,
        'AMSR_L3_seaice_DOI': amsr_doi,
        'NISE_L3_seaice_snowextent_version': nise_version,
        'NISE_L3_seaice_snowextent_PGE_version': nise_pge_version
        }

    return amsr_nise_attrs
    

def read_geo_data(file_lookup):
    """
    Read AMSR and NISE lat/lon coordinates. Coordinates are separated into
    different Northern Hemisphere and Southern Hemisphere grids.

    Parameters
    ----------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.

    Returns
    -------
    geo_data : dict
        Lat/lon coordinates of the AMSR and NISE grids. Contains separate grid
        coordinates for the Northern and Southern Hemisphere.

    """    
    
    geo_data = {
        f'{amsr_product_name}':{},
        f'{nise_product_name}':{}
        }
    
    # NISE projection x/y points are the same for both the NH and SH grids
    nise_target_crs = 'epsg:4326'
    
    nise_ul_x_m = -9036842.762500
    nise_ul_y_m = 9036842.762500
    nise_lr_x_m = 9036842.762500
    nise_lr_y_m = -9036842.762500
    
    nrows = ncols = 721
    nise_x = np.linspace(nise_ul_x_m, nise_lr_x_m, nrows)
    nise_y = np.linspace(nise_ul_y_m, nise_lr_y_m, ncols)
    nise_X, nise_Y = np.meshgrid(nise_x, nise_y)
    
    for hemi in hemis:
        # --- AMSR ---
        # Read AMSR geolocation data from first file in file_lookup
        # (lats/lons are provided in files)
        amsr_grid_name = amsr_grid_names[hemi]
        
        try:
            amsr_geo_ds = xr.open_dataset(
                file_lookup[f'{amsr_product_name}']\
                [list(file_lookup[f'{amsr_product_name}'].keys())[0]],
                group=f'/HDFEOS/GRIDS/{amsr_grid_name}/'
                )
            amsr_lon = amsr_geo_ds.lon.data
            amsr_lat = amsr_geo_ds.lat.data
            
            geo_data[f'{amsr_product_name}'][f'lon_{hemi}'] = amsr_lon
            geo_data[f'{amsr_product_name}'][f'lat_{hemi}'] = amsr_lat
        except:
            pass


        # --- NISE ---
        # NISE geolocation: create arrays of lons and lats, using the projection
        # information from NSIDC user's guide

        # Define NISE projection for the given hemisphere using proj4 string
        lat0 = nise_lat0[hemi]
        nise_proj4 = f'+proj=laea +lat_0={lat0} +lon_0=0 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs'
        nise_crs = pyproj.CRS.from_proj4(nise_proj4)
        nise_proj_to_latlon = pyproj.Transformer.from_crs(nise_crs, nise_target_crs, always_xy=True)
        nise_lon, nise_lat = nise_proj_to_latlon.transform(nise_X, nise_Y)
        
        # Change infinite values to nan
        nise_lon[nise_lon == float('+inf')] = np.nan
        nise_lat[nise_lat == float('+inf')] = np.nan
        
        geo_data[f'{nise_product_name}'][f'lon_{hemi}'] = nise_lon.astype('float64')
        geo_data[f'{nise_product_name}'][f'lat_{hemi}'] = nise_lat.astype('float64')

    return geo_data


def read_product_data(file_lookup):
    """
    Read data arrays from AMSR and NISE product files. Also change data values
    and create separate "interpolation" and "output" masks to fine tune the 
    details of interpolation to TIRS scenes and masking in output data.

    Parameters
    ----------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.

    Returns
    -------
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.

    """
    
    product_data = {
        f'{amsr_product_name}':{
            'seaice_concentration_nh':{},
            'seaice_concentration_sh':{},
            'seaice_concentration_interp_mask_nh':{},
            'seaice_concentration_interp_mask_sh':{},
            'seaice_concentration_out_mask_nh':{},
            'seaice_concentration_out_mask_sh':{}
            },
        f'{nise_product_name}':{
            'seaice_concentration_nh':{},
            'seaice_concentration_sh':{},
            'snow_fraction_nh':{},
            'snow_fraction_sh':{},
            'seaice_concentration_interp_mask_nh':{},
            'seaice_concentration_interp_mask_sh':{},
            'snow_fraction_interp_mask_nh':{},
            'snow_fraction_interp_mask_sh':{},
            'seaice_concentration_out_mask_nh':{},
            'seaice_concentration_out_mask_sh':{},
            'snow_fraction_out_mask_nh':{},
            'snow_fraction_out_mask_sh':{}
            }
        }
    
    # Read AMSR and NISE product data for each day in file_lookup
    # Skip days with missing data
    if len(file_lookup[f'{amsr_product_name}'].keys()) > 0:
        for tstr in list(file_lookup[f'{amsr_product_name}'].keys()):
            # Separate product arrays for each hemisphere
            for hemi in hemis:
                product_data = _read_AMSR_product_data(
                    file_lookup, tstr, hemi, product_data)
    
    if len(file_lookup[f'{nise_product_name}'].keys()) > 0:
        for tstr in list(file_lookup[f'{nise_product_name}'].keys()):
            # NISE HDF parent dataset is the same for both hemispheres
            nise_hdf = pyhdf.SD.SD(
                file_lookup[f'{nise_product_name}'][tstr],
                pyhdf.SD.SDC.READ
                )
            # Separate product arrays for each hemisphere
            for hemi in hemis:
                product_data = _read_NISE_product_data(
                    file_lookup, tstr, hemi, product_data, nise_hdf)

    
    return product_data

def _read_AMSR_product_data(file_lookup, tstr, hemi, product_data):
    """
    Helper to `read_product_data` that reads from the AMSR sea ice product files,
    edits product values, and creates AMSR interpolation and output masks.

    Parameters
    ----------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.
    tstr : str
        Product date in string format '%Y%m%d'.
    hemi : str
        Hemisphere ('nh' or 'sh').
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.

    Returns
    -------
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.

    """
    
    # AMSR product data
    amsr_grid_name = amsr_grid_names[hemi]
    amsr_sic_varname = amsr_sic_varnames[hemi]
    
    amsr_ds = xr.open_dataset(
        file_lookup[f'{amsr_product_name}'][tstr],
        group=f'/HDFEOS/GRIDS/{amsr_grid_name}/Data Fields/'
        )
    amsr_sic = amsr_ds[f'{amsr_sic_varname}'].data.astype('float64')
    
    amsr_interp_mask = np.ones(shape=amsr_sic.shape)
    # Interpolation mask: only mask missing values (land points will be handled
    # by the output mask)
    amsr_interp_mask[amsr_sic == 110] = 0
    product_data[f'{amsr_product_name}'][f'seaice_concentration_interp_mask_{hemi}'][tstr] = amsr_interp_mask
    
    amsr_out_mask = np.ones(shape=amsr_sic.shape)
    # Output mask: mask missing values (110) and land (120)
    amsr_out_mask[np.isin(amsr_sic, [110, 120])] = 0
    product_data[f'{amsr_product_name}'][f'seaice_concentration_out_mask_{hemi}'][tstr] = amsr_out_mask
    
    # Product array: change land values to 0. Open water also has a value of 0.
    # - This is to allow > 0 SIC values in output when a TIRS scene is located
    #   within a mixture of water and land AMSR grid cells.
    amsr_sic[amsr_sic == 120] = 0
    # Convert AMSR sea ice data to fractional units from percentage
    product_data[f'{amsr_product_name}'][f'seaice_concentration_{hemi}'][tstr] = amsr_sic*0.01
    
    return product_data

def _read_NISE_product_data(file_lookup, tstr, hemi, product_data, nise_hdf):
    """
    Helper to `read_product_data` that reads from the NISE sea ice and snow cover
    product files, edits product values, and creates NISE interpolation and
    output masks.

    Parameters
    ----------
    file_lookup : dict
        Date and file path of each daily L3 AMSR and NISE product file within
        the time range covered by tirs_times.
    tstr : str
        Product date in string format '%Y%m%d'.
    hemi : str
        Hemisphere ('nh' or 'sh').
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.
    nise_hdf : pyhdf.SD.SD
        HDF dataset object for the NISE data for the given day.

    Returns
    -------
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.

    """
    
    """
    NISE product values (from NSIDC user's guide and documentation in product files):
    - 0 = snow-free land
    - 1-100 = SIC (%)
    - 101 = permanent ice (Greenland, Antarctica)
    - 102 = not used
    - 103 = dry snow
    - 104 = wet snow
    - 105-251 = not used
    - 252 = mixed pixels at coastlines ("unable to reliably apply microwave algorithms")
    - 253 = suspect ice value ("pixel suspected of having ice")
    - 254 = corner points (undefined)
    - 255 = ocean
    """

    # See notebook: "AMSR and NISE data development notes + ice shelves" for details
    # on how it was determined that NH data is index 0 and SH data is index 2
    if hemi == 'nh':
        nise_extent_data = nise_hdf.select(0)
    elif hemi == 'sh':
        nise_extent_data = nise_hdf.select(2)
        
    nise_sic_snow = nise_extent_data[:,:].astype(np.float64)

    # NISE interpolation mask: only mask missing data (suspect ice, corner points)
    # - Same interpolation mask can be used for both SIC and snow cover
    nise_interp_mask = np.ones(shape=nise_sic_snow.shape)
    nise_interp_mask[np.isin(nise_sic_snow, [253,254])] = 0
    product_data[f'{nise_product_name}'][f'seaice_concentration_interp_mask_{hemi}'][tstr] = nise_interp_mask
    product_data[f'{nise_product_name}'][f'snow_fraction_interp_mask_{hemi}'][tstr] = nise_interp_mask

    # Split the product and mask arrays out into separate sea ice and snow cover arrays
    
    # --- Sea ice ---
    nise_sic = nise_sic_snow.copy()

    # SIC output mask: mask land areas and mixed coastline pixels, suspect ice
    # values, and corner points
    nise_sic_out_mask = np.ones(shape=nise_sic.shape)
    nise_sic_out_mask[np.isin(nise_sic, [0,101,103,104,252,253,254])] = 0
    product_data[f'{nise_product_name}'][f'seaice_concentration_out_mask_{hemi}'][tstr] = nise_sic_out_mask
    
    # SIC product array: set land, mixed coastline, and ocean pixels to 0
    # - This allows TIRS scenes with mixed land and water surroundings to have a 
    #   value > 0.
    # - Snow-free land already = 0. Change the permanent ice and the snow
    #   categories to also = 0.
    nise_sic[np.isin(nise_sic, [101,103,104,252,255])] = 0
    # Convert NISE sea ice data to fractional units from percentage
    product_data[f'{nise_product_name}'][f'seaice_concentration_{hemi}'][tstr] = nise_sic*0.01
    
    # --- Snow cover ---
    nise_snow_cover = nise_sic_snow.copy()
    # Change NISE SIC values of 1% to 2% - since the value of 1 will be used to
    # indicate snow cover (see below).
    nise_snow_cover[nise_snow_cover == 1] = 2

    # Snow fraction output mask: mask water areas (including SIC=1-100%), mixed
    # coastline pixels, suspect ice values, and corner points
    # - Note that SIC=1% has been changed to 2% - see above
    nise_snow_fraction_out_mask = np.ones(shape=nise_snow_cover.shape)
    nise_snow_fraction_out_mask[np.isin(nise_snow_cover, \
        np.concatenate((np.arange(2,101,1), np.array([252,253,254,255]))))] = 0
    product_data[f'{nise_product_name}'][f'snow_fraction_out_mask_{hemi}'][tstr] = \
        nise_snow_fraction_out_mask
    
    # Snow cover product array: classify the following NISE cell types as snow
    # covered: permanent ice, dry snow, and wet snow
    nise_snow_cover[np.isin(nise_snow_cover, [101,103,104])] = 1
    # Set sea ice, mixed coastline, and ocean points = 0 in product array, so 
    # that TIRS scenes with mixed land and water surroundings can have a snow
    # cover value > 0.
    nise_snow_cover[np.isin(nise_snow_cover, \
        np.concatenate((np.arange(2,101,1), np.array([252,255]))))] = 0

    product_data[f'{nise_product_name}'][f'snow_fraction_{hemi}'][tstr] = nise_snow_cover
    
    return product_data


def interp_AMSR_NISE_to_TIRS(tirs_times, tirs_lons, tirs_lats, tirs_shape, 
                             geo_data, product_data,
                             amsr_aux_vars, nise_aux_vars):
    """
    Interpolate AMSR sea ice and NISE sea ice and snow cover data to TIRS scene
    center points.

    Parameters
    ----------
    tirs_times : list
        Time of each TIRS scene (type: datetime.datetime).
    tirs_lons : numpy.ndarray
        Longitude of TIRS scene center points.
    tirs_lats : numpy.ndarray
        Latitude of TIRS scene center points.
    tirs_shape : tuple
        Shape of TIRS scenes in PREFIRE granule (atrack, xtrack).
    geo_data : dict
        Lat/lon coordinates of the AMSR and NISE grids. Contains separate grid
        coordinates for the Northern and Southern Hemisphere.
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.
    amsr_aux_vars : list
        List of AMSR Aux product output variables.
    nise_aux_vars : list
        List of NISE Aux product output variables.

    Returns
    -------
    output_data : dict
        AMSR sea ice and NISE sea ice and snow fraction data interpolated to TIRS
        scene center points.

    """
    
    tirs_tstrs = np.array([t.strftime('%Y%m%d') for t in tirs_times], dtype='str')
    tirs_lons_1d = np.ravel(tirs_lons)
    tirs_lats_1d = np.ravel(tirs_lats)
    
    output_data_1d = {
        v:np.full(len(tirs_lats_1d), FillValue) for v in amsr_aux_vars
        }
    
    for v in nise_aux_vars:
        output_data_1d[v] = np.full(len(tirs_lats_1d), FillValue)
    
    # Loop through the unique days on which TIRS scenes occur and determine which
    # scenes are in the NH / SH. This narrows the AMSR & NISE product data
    # to a single gridded array for the given day and hemisphere that can then
    # be passed to the ESMF interpolation function. It is analogous to finding
    # the specific ATMS granule to use for ESMF interpolation to TIRS scenes 
    # in the ATMS Anc code.
    
    # Handle missing AMSR and NISE product files
    # Check if either AMSR *or* NISE product files exist for any date in
    #  PREFIRE granule before proceeding
    tstrs = list(product_data[f'{amsr_product_name}']['seaice_concentration_nh'].keys())
    tstrs.extend(
        list(product_data[f'{nise_product_name}']['seaice_concentration_nh'].keys())
        )
    tstrs_unique = set(tstrs)
    
    if len(tstrs_unique) > 0:
        for tstr in tstrs:
            for hemi in hemis:
                # First, create ESMF Field of TIRS points to which data will be
                # interpolated (separate Field for NH and SH). These destination
                # Fields will be used for both the AMSR and NISE data interpolation.
                locstream_tirs, tirs_ixs_day = _create_locstream_tirs(tstr, hemi,
                    tirs_tstrs, tirs_lons_1d, tirs_lats_1d)
                if locstream_tirs is None:
                    continue
                
                # Process AMSR data if file exists for the date
                if tstr in product_data[f'{amsr_product_name}']['seaice_concentration_nh']:
                    for var in amsr_aux_vars:
                        output_data_1d = _interp_data(geo_data, product_data, tstr, hemi,
                                                      locstream_tirs, tirs_ixs_day,
                                                      var, output_data_1d)
                
                # Process NISE data if file exists for the date
                if tstr in product_data[f'{nise_product_name}']['seaice_concentration_nh']:
                    for var in nise_aux_vars:
                        output_data_1d = _interp_data(geo_data, product_data, tstr, hemi,
                                                      locstream_tirs, tirs_ixs_day,
                                                      var, output_data_1d)
                
    # Reshape output data to 2D
    output_data = {}
    for v in output_data_1d.keys():
        output_data[v] = np.reshape(output_data_1d[v], tirs_lons.shape)
    
    return output_data

def _create_locstream_tirs(tstr, hemi, tirs_tstrs, tirs_lons_1d, tirs_lats_1d):
    """
    Helper to `interp_AMSR_NISE_to_TIRS` that creates an ESMF locstream of TIRS
    scene center points.

    Parameters
    ----------
    tstr : str
        Product date in string format '%Y%m%d'.
    hemi : str
        Hemisphere ('nh' or 'sh').
    tirs_tstrs : numpy.ndarray
        1D array of TIRS times in string format '%Y%m%d'.
    tirs_lons_1d : numpy.ndarray
        1D array of longitude of TIRS scene center points.
    tirs_lats_1d : numpy.ndarray
        1D array of latitude of TIRS scene center points.

    Returns
    -------
    locstream_tirs : esmpy.LocStream
        ESMF locstream of TIRS scene center points.
    tirs_ixs_day : numpy.ndarray
        1D indices of TIRS scenes that match to the same day as `tstr` and the same
        hemisphere as `hemi`.

    """
    
    if hemi == 'nh':
        tirs_ixs_day = np.where(np.logical_and(tirs_tstrs==tstr,
                                               tirs_lats_1d > 0))[0]
    elif hemi == 'sh':
        tirs_ixs_day = np.where(np.logical_and(tirs_tstrs==tstr,
                                               tirs_lats_1d < 0))[0]

    if len(tirs_ixs_day) == 0:
        return (None, None)  # No PREFIRE FOVs for this hemisphere
    else:
        npts_tirs = len(tirs_ixs_day)
        locstream_tirs = esmpy.LocStream(npts_tirs,
                                         coord_sys=esmpy.CoordSys.SPH_DEG)
        locstream_tirs['ESMF:Lon'] = tirs_lons_1d[tirs_ixs_day].astype(
                                                                     'float64')
        locstream_tirs['ESMF:Lat'] = tirs_lats_1d[tirs_ixs_day].astype(
                                                                     'float64')
        return (locstream_tirs, tirs_ixs_day)


def _interp_data(geo_data, product_data, tstr, hemi, locstream_tirs, tirs_ixs_day, 
                 var, output_data_1d):
    """
    Helper to `interp_AMSR_NISE_to_TIRS` that interpolates data for the given
    variable to TIRS scenes that match to the same day as `tstr` and the same
    hemisphere as `hemi`.

    Parameters
    ----------
    geo_data : dict
        Lat/lon coordinates of the AMSR and NISE grids. Contains separate grid
        coordinates for the Northern and Southern Hemisphere.
    product_data : dict
        Arrays of AMSR and NISE product data, along with separate mask arrays
        used for ESMF interpolation and for masking Aux output data.
    tstr : str
        Product date in string format '%Y%m%d'.
    hemi : str
        Hemisphere ('nh' or 'sh').
    locstream_tirs : esmpy.LocStream
        ESMF locstream of TIRS scene center points.
    tirs_ixs_day : numpy.ndarray
        1D indices of TIRS scenes that match to the same day as `tstr` and the same
        hemisphere as `hemi`.
    var : str
        AMSR or NISE Aux output product variable name.
    output_data_1d : dict
        Contains 1D arrays of output product data interpolated to TIRS scenes.

    Returns
    -------
    output_data_1d : dict
        Contains 1D arrays of output product data interpolated to TIRS scenes.

    """
        
    sat = var[0:4]
    varname = var.split('L3_')[1]
    
    if sat == 'AMSR':
        product_name = amsr_product_name
    elif sat == 'NISE':
        product_name = nise_product_name
        
    # Source Grid coordinates and Field data: transpose arrays for efficiency.
    # See Nathan Wendt ESMPy tutorial notebook for details.
    
    lon_T = geo_data[f'{product_name}'][f'lon_{hemi}'][:].T
    lat_T = geo_data[f'{product_name}'][f'lat_{hemi}'][:].T

    src_grid = esmpy.Grid(
        np.array(lon_T.shape),
        staggerloc=esmpy.StaggerLoc.CENTER,
        coord_sys=esmpy.CoordSys.SPH_DEG
        )
    
    # Assign coordinates to grid
    src_grid_lon = src_grid.get_coords(0)
    src_grid_lon[...] = lon_T
    src_grid_lat = src_grid.get_coords(1)
    src_grid_lat[...] = lat_T
    
    # Assign mask to grid in order to handle missing data
    interp_mask = product_data[f'{product_name}'][f'{varname}_interp_mask_{hemi}'][tstr].T
    src_grid.add_item(
        esmpy.GridItem.MASK,
        staggerloc=esmpy.StaggerLoc.CENTER
        )
    src_grid.mask[0][...] = interp_mask
    
    # Interpolate data
    src_field = esmpy.Field(
        src_grid, typekind=esmpy.TypeKind.R8
        )
    src_field.data[:] = product_data[product_name][f'{varname}_{hemi}'][tstr].T
        
    dst_field = esmpy.Field(
        locstream_tirs, typekind=esmpy.TypeKind.R8
        )
    
    # Mask for product data
    regrid = esmpy.Regrid(
        src_field, dst_field,
        regrid_method=esmpy.RegridMethod.BILINEAR,
        unmapped_action=esmpy.UnmappedAction.IGNORE,
        src_mask_values=[0]
        )

    regrid(src_field, dst_field)
    output_array = dst_field.data
    
    # Now interpolate the mask variable. Any TIRS scenes with a value of 0 in
    # the interpolated mask variable will be set to FillValue. These are scenes
    # for which *all* of the source grid cells that contributed interpolation
    # weights were masked.
    src_field_mask = esmpy.Field(
        src_grid, typekind=esmpy.TypeKind.R8
        )
    src_field_mask.data[:] = product_data[product_name][f'{varname}_out_mask_{hemi}'][tstr].T
        
    dst_field_mask = esmpy.Field(
        locstream_tirs, typekind=esmpy.TypeKind.R8
        )
    regrid_mask = esmpy.Regrid(
        src_field_mask, dst_field_mask,
        regrid_method=esmpy.RegridMethod.BILINEAR,
        unmapped_action=esmpy.UnmappedAction.IGNORE
        )

    regrid_mask(src_field_mask, dst_field_mask)
    output_array[np.where(dst_field_mask.data == 0)] = FillValue
        
    output_data_1d[var][tirs_ixs_day] = output_array
    
    return output_data_1d


def write_AMSR_NISE_aux_file(output_data, output_fpath, product_filespecs,
                             amsr_nise_attrs, tirs_shape):
    """
    Write temporary file containing the Aux product data derived from AMSR and
    NISE products.

    Parameters
    ----------
    output_data : dict
        AMSR sea ice and NISE sea ice and snow fraction data interpolated to
        TIRS scene center points.
    output_fpath : str
        Path to output file.
    product_filespecs : dict
        File specs for Aux-Sat product.
    amsr_nise_attrs : dict
        AMSR and NISE product metadata extracted from product files.
    tirs_shape : tuple
        Shape of TIRS scenes in PREFIRE granule (atrack, xtrack).

    Returns
    -------
    None.

    """
    
    out_f = n.Dataset(output_fpath, 'w')
    dim_atrack = out_f.createDimension('atrack', tirs_shape[0])
    dim_xtrack = out_f.createDimension('xtrack', tirs_shape[1])
    
    for variable in output_data:
        var_data = output_data[variable]
        dtype = product_filespecs['Aux-Sat'][variable]['np_dtype']
        file_variable = out_f.createVariable(
            variable, dtype, ('atrack','xtrack'), fill_value=FillValue
        )
        file_variable[:] = var_data
        
    # Set attributes with AMSR and NISE product metadata
    for attr in amsr_nise_attrs.keys():
        out_f.setncattr(attr, amsr_nise_attrs[attr])

    out_f.close()
    

def parse_args():
    """
    Parse command line arguments passed to script at runtime.

    Returns
    -------
    args.tirs_fpath : str
        Path to a PREFIRE file containing a TIRS Geometry group.
    args.product_specs_fpath : str
        Path to product file specs.
    args.input_products_dir : str
        Parent directory of input product files.
    args.tmpfiles_dir : str
        Path to temporary directory in which intermediate processing files are stored.
    artp : 3-tuple (str, int, int)
        atrack_range_to_process, a 3-tuple containing the dimension name to
        subset, and start and stop indices (NumPy indexing convention) in the
        given granule
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('tirs_fpath', 
                        help='path to PREFIRE file with a TIRS Geometry group')
    parser.add_argument('product_specs_fpath',
                        help='path to product file specs')
    parser.add_argument('input_products_dir', 
                        help='path to AMSR and NISE input products parent directory')
    parser.add_argument('tmpfiles_dir', 
                        help='path to temporary directory in which intermediate processing files are stored')
    parser.add_argument('artp_strrep', 
                        help=("colon-delimited string containing the dimension "
                              "name to subset, and start and stop indices "
                              "(NumPy indexing convention) in the given "
                              "granule"))
    args = parser.parse_args()

    # Parse artp_strrep:
    tokens = args.artp_strrep.split(':')
    if tokens[2] == "END":
        artp = (tokens[0], int(tokens[1]), None)
    else:
        artp = (tokens[0], int(tokens[1]), int(tokens[2]))
        
    return (args.tirs_fpath, args.product_specs_fpath, args.input_products_dir,
            args.tmpfiles_dir, artp)


def main(tirs_fpath, product_specs_fpath, input_products_dir, tmpfiles_dir,
         artp):
    
    tirs_fname = os.path.basename(tirs_fpath)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Started AMSR and NISE processing of {tirs_fname} at {now}')
    
    # Get TIRS geo data
    tirs_times, tirs_shape, tirs_lons, tirs_lats = construct_tirs_geo(
                                              tirs_fpath, artp, polygons=False)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Constructed TIRS geo objects at {now}')

    # Compile AMSR / NISE file lookup and get product metadata
    file_lookup = compile_file_lookup(input_products_dir, tirs_times)
    amsr_nise_attrs = get_AMSR_NISE_attrs(file_lookup)
    
    # Read AMSR / NISE geolocation and product data
    geo_data = read_geo_data(file_lookup)
    product_data = read_product_data(file_lookup)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Read geolocation and product data at {now}')

    # Interpolate AMSR and NISE data to TIRS scene center points
    output_data = interp_AMSR_NISE_to_TIRS(
        tirs_times, tirs_lons, tirs_lats, tirs_shape,
        geo_data, product_data,
        amsr_aux_vars, nise_aux_vars
        )
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Interpolated data to TIRS scenes at {now}')

    # Write output
    tokens = os.path.basename(tirs_fpath).split('_')
    output_fpath = os.path.join(tmpfiles_dir, '_'.join(tokens[0:2])+
                                        "_Aux-AMSR_NISE_"+'_'.join(tokens[3:]))
    with open(product_specs_fpath, 'r') as f:
        product_filespecs = json.load(f)
    write_AMSR_NISE_aux_file(
        output_data, output_fpath, product_filespecs, amsr_nise_attrs, tirs_shape
        )
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Wrote AMSR and NISE Aux data at {now}')

    
if __name__ == '__main__':
    main(*parse_args())
