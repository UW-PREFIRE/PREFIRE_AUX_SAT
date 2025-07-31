import argparse
import numpy as np
import netCDF4 as n
import re
import os
import glob
import json

from PREFIRE_AUX_MET.utils.classify_sfc_type_prelim import \
    classify_land_type, classify_water_type
from PREFIRE_PRD_GEN.file_read import get_PREFIRE_Lx_field


FillValue = -9999.0


def read_input_data(tmpfiles_dir, auxmet_fpath, jpss_sats, artp):
    """
    Read input data from Aux-Met file and Aux-Sat temporary files.

    Parameters
    ----------
    tmpfiles_dir : str
        Path to temporary directory in which intermediate processing files are stored.
    auxmet_fpath : str
        Path to Aux-Met file for the given PREFIRE granule.
    jpss_sats : list
        JPSS satellites with data to process.

    Returns
    -------
    input_data : dict
        Input L1B Geometry, Aux-Met, and Aux-Sat data for final surface type
        classification.

    """
    
    input_data = {}
    
    # L1B Geometry (via Aux-Met) and Aux-Met data
    
    # FillValues in input data are automatically masked upon reading in netCDF
    # file, so tests in _choose_snow_data_source and _choose_seaice_data_source
    # can check for mask values when iterating through the hierarchy of data
    # sources.
    with n.Dataset(auxmet_fpath, 'r') as nc:
        tmpA = get_PREFIRE_Lx_field(nc, "Geometry", "latitude", artp)
        tirs_shape = np.shape(tmpA)
        input_data['latitude'] = np.ravel(tmpA)

        input_data['land_fraction'] = np.ravel(
                   get_PREFIRE_Lx_field(nc, "Geometry", "land_fraction", artp))
        input_data['antarctic_ice_shelf_fraction'] = np.ravel(
                    get_PREFIRE_Lx_field(nc, "Aux-Met",
                                         "antarctic_ice_shelf_fraction", artp))
        input_data['antarctic_land_fraction'] = np.ravel(
                         get_PREFIRE_Lx_field(nc, "Aux-Met",
                                              "antarctic_land_fraction", artp))
        input_data['VIIRS_surface_type'] = np.reshape(
               get_PREFIRE_Lx_field(nc, "Aux-Met", "VIIRS_surface_type", artp),
                           (-1, nc['Aux-Met']['VIIRS_surface_type'].shape[-1]))
        input_data['auxmet_seaice'] = np.ravel(
             get_PREFIRE_Lx_field(nc, "Aux-Met", "seaice_concentration", artp))
        input_data['auxmet_snow'] = np.ravel(
             get_PREFIRE_Lx_field(nc, "Aux-Met", "snow_cover", artp))
    
    # VIIRS L3 snow cover
    for jpss_sat in jpss_sats:
        try:
            aux_viirs_fpath = glob.glob(os.path.join(tmpfiles_dir,
                                       f"PREFIRE*Aux-VIIRS-{jpss_sat}*.nc"))[0]
            with n.Dataset(aux_viirs_fpath, 'r') as nc:
                input_data[f'{jpss_sat}_VIIRS_L3_snow_fraction_mean'] = \
                                   nc['VIIRS_L3_snow_fraction_mean'][:].ravel()
        except:
            input_data[f"{jpss_sat}_VIIRS_L3_snow_fraction_mean"] = (
                             np.ma.masked_all(input_data["auxmet_snow"].shape))

    # AMSR sea ice and NISE sea ice + snow fraction
    aux_amsr_nise_fpath = glob.glob(os.path.join(tmpfiles_dir,
                                               "PREFIRE*Aux-AMSR_NISE*.nc"))[0]
    with n.Dataset(aux_amsr_nise_fpath, 'r') as nc:
        input_data['AMSR_L3_seaice_concentration'] = nc['AMSR_L3_seaice_concentration'][:].ravel()
        input_data['NISE_L3_seaice_concentration'] = nc['NISE_L3_seaice_concentration'][:].ravel()
        input_data['NISE_L3_snow_fraction'] = nc['NISE_L3_snow_fraction'][:].ravel()
    
    return input_data, tirs_shape


def calc_sfc_type_final(input_data, tirs_shape):
    """
    Calculate the "final" PREFIRE surface type classification, which uses
    L1B Geometry, Aux-Met, and Aux-Sat data.

    Parameters
    ----------
    input_data : dict
        Input L1B Geometry, Aux-Met, and Aux-Sat data for final surface type
        classification.

    Returns
    -------
    output_data : dict
        Contains "final" surface type classification for each TIRS scene, along
        with variables tracking which data sources were used to determine each
        scene's land fraction, sea ice, and snow cover.
    tirs_shape : tuple
        Shape of TIRS scenes in PREFIRE granule (atrack, xtrack).

    """
        
    # Initialize output data
    output_data = {}
    
    sfc_type_1d = []
    seaice_data_source_1d = []
    snow_data_source_1d = []
    
    # Loop through TIRS scenes and classify surface type
    for i, lat in enumerate(input_data['latitude']):
        # Antarctic scenes: use PREFIRE BAS Antarctic coastline data for land
        # and ice shelf fraction
        if lat < -60:
            lf_scene = input_data['antarctic_land_fraction'][i]
            isf_scene = input_data['antarctic_ice_shelf_fraction'][i]
            
            # Ice shelf
            if ((lf_scene + isf_scene) > 0.5) and (isf_scene > lf_scene):
                # PREFIRE surface type 5 = Antarctic ice shelf
                sfc_type_scene = 5
                seaice_data_source_scene = FillValue
                snow_data_source_scene = FillValue
            # Land
            elif ((lf_scene + isf_scene) > 0.5) and (lf_scene > isf_scene):
                seaice_data_source_scene = FillValue
                snow_data_source_scene, snow_cover_scene = \
                    _choose_snow_data_source(i, input_data)
                sfc_type_scene = classify_land_type(
                    input_data['VIIRS_surface_type'][i], snow_cover_scene
                    )
            # Water
            else:
                snow_data_source_scene = FillValue
                seaice_data_source_scene, seaice_scene = \
                    _choose_seaice_data_source(i, input_data)
                sfc_type_scene = classify_water_type(seaice_scene)
        # Non-Antarctic scenes: use Copernicus DEM (from L1B) for land fraction
        else:
            lf_scene = input_data['land_fraction'][i]
            
            if lf_scene >= 0.5:
                seaice_data_source_scene = FillValue
                snow_data_source_scene, snow_cover_scene = \
                    _choose_snow_data_source(i, input_data)
                sfc_type_scene = classify_land_type(
                    input_data['VIIRS_surface_type'][i], snow_cover_scene
                    )
            else:
                snow_data_source_scene = FillValue
                seaice_data_source_scene, seaice_scene = \
                    _choose_seaice_data_source(i, input_data)
                sfc_type_scene = classify_water_type(seaice_scene)

        sfc_type_1d.append(sfc_type_scene)
        seaice_data_source_1d.append(seaice_data_source_scene)
        snow_data_source_1d.append(snow_data_source_scene)

    # Reshape to 2D and populate output data
    output_data['merged_surface_type_final'] = \
        np.reshape(np.array(sfc_type_1d), tirs_shape)
    output_data['merged_seaice_final_data_source'] = \
        np.reshape(np.array(seaice_data_source_1d), tirs_shape)
    output_data['merged_snow_final_data_source'] = \
        np.reshape(np.array(snow_data_source_1d), tirs_shape)
        
    return output_data

def _choose_snow_data_source(i, input_data):
    """
    Helper to calc_sfc_type_final. Choose snow cover data source and extract
    snow cover value for a single TIRS scene.

    Parameters
    ----------
    i : int
        TIRS scene 1D index.
    input_data : dict
        Input L1B Geometry, Aux-Met, and Aux-Sat data for final surface type
        classification.

    Returns
    -------
    snow_data_source_scene : int
        Integer type code for snow data source.
    snow_cover_scene : float
        Snow cover for a single TIRS scene.

    """

    # ** As of 2023-10-26, SNPP and NOAA-20 L3 snow cover data are available.
    # Change this data hierarchy if VIIRS product availability changes.    
    if not np.ma.is_masked(input_data['noaa20_VIIRS_L3_snow_fraction_mean'][i]):
        # NOAA-20 VIIRS L3 snow data source = 3
        snow_data_source_scene = 3
        snow_cover_scene = input_data['noaa20_VIIRS_L3_snow_fraction_mean'][i]
    elif not np.ma.is_masked(input_data['snpp_VIIRS_L3_snow_fraction_mean'][i]):    
        # SNPP VIIRS L3 snow data source = 4
        snow_data_source_scene = 4
        snow_cover_scene = input_data['snpp_VIIRS_L3_snow_fraction_mean'][i]
    elif not np.ma.is_masked(input_data['NISE_L3_snow_fraction'][i]):
        # NISE L3 snow data source = 6
        snow_data_source_scene = 6
        snow_cover_scene = input_data['NISE_L3_snow_fraction'][i]
    elif not np.ma.is_masked(input_data['auxmet_snow'][i]):
        # GEOS-IT snow data source = 7
        snow_data_source_scene = 7
        snow_cover_scene = input_data['auxmet_snow'][i]
    else:
        snow_data_source_scene = 0
        snow_cover_scene = FillValue
        
    return snow_data_source_scene, snow_cover_scene

def _choose_seaice_data_source(i, input_data):
    """
    Helper to calc_sfc_type_final. Choose sea ice data source and extract
    sea ice value for a single TIRS scene.

    Parameters
    ----------
    i : int
        TIRS scene 1D index.
    input_data : dict
        Input L1B Geometry, Aux-Met, and Aux-Sat data for final surface type
        classification.

    Returns
    -------
    seaice_data_source_scene : int
        Integer type code for snow data source.
    seaice_scene : float
        Sea ice concentration for a single TIRS scene.

    """
    
    if not np.ma.is_masked(input_data['AMSR_L3_seaice_concentration'][i]):
        # AMSR L3 snow data source = 1
        seaice_data_source_scene = 1
        seaice_scene = input_data['AMSR_L3_seaice_concentration'][i]
    elif not np.ma.is_masked(input_data['NISE_L3_seaice_concentration'][i]):
        # NISE L3 sea ice source = 6
        seaice_data_source_scene = 6
        seaice_scene = input_data['NISE_L3_seaice_concentration'][i]
    elif not np.ma.is_masked(input_data['auxmet_seaice'][i]):
        # GEOS-IT sea ice data source = 7
        seaice_data_source_scene = 7
        seaice_scene = input_data['auxmet_seaice'][i]
    else:
        seaice_data_source_scene = 0
        seaice_scene = FillValue
    
    return seaice_data_source_scene, seaice_scene


def write_output(output_fpath, output_data, product_filespecs, tirs_shape):
    """
    Write temporary file containing the "final" PREFIRE surface type
    classification.

    Parameters
    ----------
    output_fpath : str
        Path to output file.
    output_data : dict
        Contains "final" surface type classification for each TIRS scene, along
        with variables tracking which data sources were used to determine each
        scene's land fraction, sea ice, and snow cover.
    product_filespecs : dict
        File specs for Aux-Sat product.
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
        
    out_f.close()
    

def parse_args():
    """
    Parse command line arguments passed to script at runtime.    

    Returns
    -------
    args.tmpfiles_dir : str
        Path to temporary directory in which intermediate processing files are stored.
    args.auxmet_fpath : str
        Path to Aux-Met file for the given PREFIRE granule.
    args.product_specs_fpath : str
        Path to product file specs.
    args.jpss_sats : list
        JPSS satellites with data to process.
    artp : 3-tuple (str, int, int)
        atrack_range_to_process, a 3-tuple containing the dimension name to
        subset, and start and stop indices (NumPy indexing convention) in the
        given granule
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'tmpfiles_dir',
        help='path to temporary directory in which intermediate processing files are stored'
        )
    parser.add_argument(
        'auxmet_fpath',
        help='path to Aux-Met file for the given PREFIRE granule'
        )
    parser.add_argument(
        'product_specs_fpath',
        help='path to product file specs'
        )
    parser.add_argument('artp_strrep', 
                        help=("colon-delimited string containing the dimension "
                              "name to subset, and start and stop indices "
                              "(NumPy indexing convention) in the given "
                              "granule"))
    parser.add_argument(
        '--jpss_sats', type=str,
        nargs='+',
        help='list of JPSS satellites'
        )
    args = parser.parse_args()

    # Parse artp_strrep:
    tokens = args.artp_strrep.split(':')
    if tokens[2] == "END":
        artp = (tokens[0], int(tokens[1]), None)
    else:
        artp = (tokens[0], int(tokens[1]), int(tokens[2]))

    return (args.tmpfiles_dir, args.auxmet_fpath, args.product_specs_fpath,
            args.jpss_sats, artp)


def main(tmpfiles_dir, auxmet_fpath, product_specs_fpath, jpss_sats, artp):
    
    input_data, tirs_shape = read_input_data(tmpfiles_dir, auxmet_fpath,
                                             jpss_sats, artp)
    
    output_data = calc_sfc_type_final(input_data, tirs_shape)
    
    with open(product_specs_fpath, 'r') as f:
        product_filespecs = json.load(f)

    tokens = os.path.basename(auxmet_fpath).split('_')
    output_fpath = os.path.join(tmpfiles_dir, '_'.join(tokens[0:2])+
                                    "_Aux-sfctype-final_"+'_'.join(tokens[3:]))

    write_output(output_fpath, output_data, product_filespecs, tirs_shape)


if __name__ == '__main__':
    main(*parse_args())
