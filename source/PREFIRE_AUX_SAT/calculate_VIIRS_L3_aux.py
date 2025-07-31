import argparse
import json
import re
import os
import glob
from datetime import datetime
import numpy as np
import netCDF4 as n
import shapely
import warnings

from PREFIRE_tools.utils.aux import construct_tirs_geo

FillValue = -9999.0

viirs_l3_prod_names = {
    'snpp':'VNP10C1',
    'noaa20':'VJ110C1'
    }


def calc_VIIRS_L3_snow_cloud_fraction(input_products_dir, jpss_sat,
                                      product_filespecs,
                                      tirs_shape, tirs_times, tirs_polylist):
    """
    Calculate VIIRS L3 snow and cloud fraction output for all scenes in a
    PREFIRE granule.

    Parameters
    ----------
    input_products_dir : str
        Path to parent directory of input product files used for Aux-Sat.
    jpss_sat : str
        "snpp" or "noaa20".
    product_filespecs : dict
        File specs for Aux-Sat product.
    tirs_shape : tuple
        Shape of TIRS scenes in PREFIRE granule (atrack, xtrack).
    tirs_times : list
        Time of each TIRS scene (type: datetime.datetime).
    tirs_polylist : list
        Polygon outline(s) for each TIRS scene. Each list item is a list of
        polygons that has length 2 if the TIRS polygon crosses the 

    Returns
    -------
    output_data : dict
        Contains VIIRS L3 Aux-Sat output data for each variable for the given
        JPSS satellite.
    VIIRS_CMG_product_attrs : dict
        VIIRS CMG product metadata extracted from an example product file.

    """
    
    tirs_times_daily = set([datetime(t.year, t.month, t.day) for t in tirs_times])

    # Read VIIRS L3 CMG data for each unique date in PREFIRE granule
    product_name = viirs_l3_prod_names[jpss_sat]
    product_data = {}
    
    i = 0
    for t_daily in list(tirs_times_daily):
        date_str = t_daily.strftime('%Y%j')
        
        try:
            product_fpath = glob.glob(os.path.join(input_products_dir,
                                          product_name, f"*.A{date_str}.*"))[0]
            # Get product coords and attrs from the first date's product file
            if i == 0:
                product_data[date_str], VIIRS_CMG_product_attrs, cmg_lons, cmg_lats = \
                    _read_product_data(product_fpath,
                                       product_name=product_name, attrs_and_coords=True)
            else:
                product_data[date_str] = _read_product_data(product_fpath)
            i += 1
            
        # Skip missing product files with a warning
        except IndexError:
            tstr = t_daily.strftime('%Y%m%d')
            warnings.warn(f'No {product_name} file found for {tstr}')
    
    # If no product files available, manually create CMG lats/lons so that
    # the VIIRS L3 number of points and lat/lon output variables can be
    # populated, and fill product attributes with "missing"
    if len(product_data.keys()) == 0:
        cmg_lons = np.linspace(-179.975, 179.975, 7200)
        cmg_lats = np.linspace(89.975, -89.975, 3600)
        
        VIIRS_CMG_product_attrs = {
                f'{product_name}_version':'unavailable',
                f'{product_name}_PGE_name':'unavailable',
                f'{product_name}_PGE_version':'unavailable',
                f'{product_name}_DOI':'unavailable'
            }
        
    # Create 2D array of all possible pairs of CMG lon/lat coordinates, then
    # flatten to 1D to prepare for STRtree creation
    cmg_lons_2d, cmg_lats_2d = np.meshgrid(cmg_lons, cmg_lats)
    cmg_lons_flat_1d = np.ravel(cmg_lons_2d)
    cmg_lats_flat_1d = np.ravel(cmg_lats_2d)
    
    # Build global STRtree of CMG points
    cmg_tree = shapely.STRtree(shapely.points(cmg_lons_flat_1d, y=cmg_lats_flat_1d))
    
    # Loop through TIRS scenes and fill in output data
    list_data = {
        'VIIRS_L3_num_pts_total':[],
        'VIIRS_L3_snow_fraction_num_pts_used':[],
        'VIIRS_L3_snow_fraction_mean':[],
        'VIIRS_L3_snow_fraction_stdev':[],
        'VIIRS_L3_cloud_fraction_num_pts_used':[],
        'VIIRS_L3_cloud_fraction_mean':[],
        'VIIRS_L3_cloud_fraction_stdev':[]
        }
    
    for i, tirs_polys in enumerate(tirs_polylist):
        # Find indices of TIRS polygon intersections with CMG points STRtree
        # (including handling of split polygons along date line)
        if len(tirs_polys) == 2:
            tirs_intersect_ixs = []
            for split_poly in tirs_polys:
                ixs_split_scene = cmg_tree.query(
                    split_poly, predicate='contains'
                    )
                tirs_intersect_ixs.extend(ixs_split_scene)
        elif len(tirs_polys) == 1:
            tirs_intersect_ixs = cmg_tree.query(
                tirs_polys[0], predicate='contains'
                )
        
        # Calculate output data (for the correct date, if it exists) at CMG
        # indices within each TIRS scene
        date_str = tirs_times[i].strftime('%Y%j')
        try:
            qa_scene = product_data[date_str]['Basic_QA'][tirs_intersect_ixs]
            cc_scene = product_data[date_str]['Cloud_Cover'][tirs_intersect_ixs]
            sc_scene = product_data[date_str]['Snow_Cover'][tirs_intersect_ixs]

            # "Valid" snow cover points must meet the following requirements:
            # - Not masked by VIIRS CMG product FillValue
            #   - (e.g. 239 = water, 243 = Antarctica, etc.)
            # - Basic_QA != 2
            # - Cloud cover < 0.10 (handles product quality flags, which are
            #   all values of > 200)
            valid_ixs_sc = np.where(np.logical_and(
                qa_scene != 2,
                cc_scene < 0.10
                ))
            list_data['VIIRS_L3_snow_fraction_num_pts_used'].append(len(valid_ixs_sc[0]))
                        
            # "Valid" cloud cover points must meet the following requirements:
            # - Not masked by VIIRS CMG product FillValue
            #   - (e.g. 239 = water, 243 = Antarctica, etc.)
            # - Basic_QA != 2
            # - Cloud cover is < 1, accounting for floating point imprecision 
            #   possibly introduced by conversion from percent to fraction 
            #   (handles product quality flags, which are all values of > 200)
            valid_ixs_cc = np.where(np.logical_and(
                qa_scene != 2,
                cc_scene < 1.0001
                ))
            list_data['VIIRS_L3_cloud_fraction_num_pts_used'].append(len(valid_ixs_cc[0]))
            
            # Require that at least 2 CMG points within scene have valid data
            # for snow fraction and cloud fraction output to be populated
            if len(valid_ixs_sc[0]) > 1:
                list_data['VIIRS_L3_snow_fraction_mean'].append(
                    np.mean(sc_scene[valid_ixs_sc]))
                list_data['VIIRS_L3_snow_fraction_stdev'].append(
                    np.std(sc_scene[valid_ixs_sc]))
            else:
                list_data['VIIRS_L3_snow_fraction_mean'].append(FillValue)
                list_data['VIIRS_L3_snow_fraction_stdev'].append(FillValue)

            if len(valid_ixs_cc[0]) > 1:
                list_data['VIIRS_L3_cloud_fraction_mean'].append(
                    np.mean(cc_scene[valid_ixs_cc]))
                list_data['VIIRS_L3_cloud_fraction_stdev'].append(
                    np.std(cc_scene[valid_ixs_cc]))
            else:
                list_data['VIIRS_L3_cloud_fraction_mean'].append(FillValue)
                list_data['VIIRS_L3_cloud_fraction_stdev'].append(FillValue)

        # Skip missing VIIRS product files
        # - (e.g., SNPP file for 2021-07-21 is missing from NSIDC DAAC archive),
        #   but a PREFIRE granule may continue into 2021-07-22, a date on which
        #   SNPP file is available
        except KeyError:
            list_data['VIIRS_L3_snow_fraction_num_pts_used'].append(0)
            list_data['VIIRS_L3_cloud_fraction_num_pts_used'].append(0)
            list_data['VIIRS_L3_snow_fraction_mean'].append(FillValue)
            list_data['VIIRS_L3_snow_fraction_stdev'].append(FillValue)
            list_data['VIIRS_L3_cloud_fraction_mean'].append(FillValue)
            list_data['VIIRS_L3_cloud_fraction_stdev'].append(FillValue)
    
        # Fill in VIIRS number of points data
        list_data['VIIRS_L3_num_pts_total'].append(len(tirs_intersect_ixs))
                
    # Reshape output data to 2D
    output_data = {}
    for k in list_data.keys():
        output_data[k] = np.reshape(np.array(list_data[k]), tirs_shape)
        
    return output_data, VIIRS_CMG_product_attrs


def _read_product_data(product_fpath, product_name=None, attrs_and_coords=False):
    """
    Helper to calc_VIIRS_L3_snow_cover that reads arrays of data from a single
    VIIRS 0.05deg CMG product file.
    
    Optionally, arrays of VIIRS L3 CMG 0.05 deg lat/lon coordinates are
    returned, along with product metadata read from the file.

    Parameters
    ----------
    product_fpath : str
        Path to VIIRS L3 CMG 0.05 deg snow cover product file.
    product_name : str, optional
        Shortname for VIIRS L3 CMG 0.05 deg snow cover product for the given
        JPSS satellite. The default is None.
    attrs_and_coords : bool, optional
        Whether to return VIIRS product attributes and lat/lon coordinates.
        The default is False.

    Returns
    -------
    product_data_date : dict
        Arrays of QA, cloud cover, and snow cover read from the product file.
    product_attrs_date : dict, optional
        VIIRS CMG product metadata extracted from product file.
    cmg_lons : numpy.ndarray, optional
        CMG (0.05 deg) longitudes.
    cmg_lats : numpy.ndarray, optional
        CMG (0.05 deg) latitudes.

    """
    
    product_data_date = {}
    
    with n.Dataset(product_fpath) as nc:
        product_data_date['Basic_QA'] = \
            nc['/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields/Basic_QA'][:].ravel()
        # Convert cloud cover from percentage to fraction
        product_data_date['Cloud_Cover'] = \
            nc['/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields/Cloud_Cover'][:].ravel()*0.01
        # Convert snow cover from percentage to fraction
        product_data_date['Snow_Cover'] = \
            nc['/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields/Snow_Cover'][:].ravel()*0.01
        
        if attrs_and_coords:
            # CMG dataset provides coordinates of upper left corners of grid boxes.
            # Shift coordinates by 0.025 degrees to calculate the *center* of each
            # CMG grid cell.
            cmg_lons = \
                nc['/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields/longitude'][:] + 0.025
            cmg_lats = \
                nc['/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields/latitude'][:] - 0.025
            
            product_attrs_date = {}
            product_attrs_date[f'{product_name}_version'] = \
                nc.getncattr('VersionID')
            product_attrs_date[f'{product_name}_PGE_name'] = \
                nc.getncattr('PGE_Name')
            product_attrs_date[f'{product_name}_PGE_version'] = \
                nc.getncattr('PGEVersion')
            product_attrs_date[f'{product_name}_DOI'] = \
                nc.getncattr('identifier_product_doi')
        
    if attrs_and_coords:
        return product_data_date, product_attrs_date, cmg_lons, cmg_lats
    else:
        return product_data_date


def write_output(output_fpath, output_data, product_filespecs,
                 VIIRS_CMG_product_attrs, tirs_shape):
    """
    Write out VIIRS L3 aux file for the given JPSS satellite. This contains
    all the variables previously defined, whether empty or full of data.
    Compression is included.

    Parameters
    ----------
    output_fpath : str
        Path to output VIIRS L3 aux file for the given JPSS satellite.
    output_data : dict
        Contains VIIRS L3 Aux-Sat output data for each variable for the given
        JPSS satellite.
    product_filespecs : dict
        File specs for Aux-Sat product.
    VIIRS_CMG_product_attrs : dict
        VIIRS CMG product metadata extracted from an example product file.
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
        if variable in product_filespecs['Aux-Sat']:
            dtype = product_filespecs['Aux-Sat'][variable]['np_dtype']
            file_variable = out_f.createVariable(
                variable, dtype, ('atrack','xtrack'),
                fill_value=product_filespecs['Aux-Sat'][variable]['fill_value']
            )
            file_variable[:] = output_data[variable]
            
    # Set attributes with VIIRS CMG product metadata
    # Skip if attributes are missing
    if len(VIIRS_CMG_product_attrs.keys()) > 0:
        for attr in VIIRS_CMG_product_attrs.keys():
            out_f.setncattr(attr, VIIRS_CMG_product_attrs[attr])
    
    out_f.close()
    
    
def parse_args():
    """
    Parse command line arguments passed to script at runtime.

    Returns
    -------
    args.jpss_sat : str
        "snpp" or "noaa20". 
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
    parser.add_argument('jpss_sat',
                        help='jpss satellite name ("snpp" or "noaa-20")')
    parser.add_argument('tirs_fpath',
                        help='path to PREFIRE file with a TIRS Geometry group')
    parser.add_argument('product_specs_fpath',
                        help='path to product file specs')
    parser.add_argument('input_products_dir',
                        help='path to input products parent directory')
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

    return (args.jpss_sat, args.tirs_fpath, args.product_specs_fpath,
            args.input_products_dir, args.tmpfiles_dir, artp)


def main(jpss_sat, tirs_fpath, product_specs_fpath, input_products_dir,
         tmpfiles_dir, artp):
    """
    Main block that controls the VIIRS L3 snow cover processing for the given
    JPSS satellite.

    Parameters
    ----------
    jpss_sat : str
        "snpp" or "noaa20".
    tirs_fpath : str
        Path to a PREFIRE file containing a TIRS Geometry group.
    product_specs_fpath : str
        Path to product file specs.
    input_products_dir : str
        Parent directory of input product files.
    tmpfiles_dir : str
        Path to temporary directory in which intermediate processing files are
        stored.
    artp : 3-tuple (str, int, int)
        atrack_range_to_process, a 3-tuple containing the dimension name to
        subset, and start and stop indices (NumPy indexing convention) in the
        given granule

    Returns
    -------
    None.
    """
    
    # Get output file name
    tirs_fname = os.path.basename(tirs_fpath)
    tokens = tirs_fname.split('_')
    output_fpath = os.path.join(tmpfiles_dir, '_'.join(tokens[0:2])+
                                f"_Aux-VIIRS-{jpss_sat}_"+'_'.join(tokens[3:]))
    
    print(f'Started {jpss_sat} VIIRS L3 snow cover processing of {tirs_fname} at '+
          datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    
    # Build TIRS geolocation data
    tirs_times, tirs_shape, _, _, tirs_polylist = construct_tirs_geo(
                                               tirs_fpath, artp, polygons=True)

    # Get Aux-Sat product file specs
    with open(product_specs_fpath, 'r') as f:
        product_filespecs = json.load(f)
    
    # Populate output L3 snow cover data for the full PREFIRE granule
    output_data, VIIRS_CMG_product_attrs = calc_VIIRS_L3_snow_cloud_fraction(
        input_products_dir, jpss_sat,
        product_filespecs,
        tirs_shape, tirs_times, tirs_polylist
        )
    
    # Write output file
    write_output(
        output_fpath, output_data, product_filespecs, 
        VIIRS_CMG_product_attrs, tirs_shape
        )

    print(f'Finished {jpss_sat} VIIRS L3 snow cover processing of {tirs_fname} at '+
          datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == '__main__':
    main(*parse_args())
