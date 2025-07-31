import argparse
import netCDF4 as n
import numpy as np
import os
import datetime
import glob

from PREFIRE_PRD_GEN.file_read import load_all_vars_of_nc4group, load_all_atts_of_nc4group
from PREFIRE_PRD_GEN.file_creation import write_data_fromspec
from PREFIRE_tools.utils.filesys import mkdir_p
import PREFIRE_AUX_SAT.filepaths as AUX_SAT_fpaths

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'out_dir', 
        help='provide path to output directory'
        )
    parser.add_argument(
        'auxmet_fpath', 
        help='path to Aux-Met file for the given PREFIRE granule'
        )
    parser.add_argument(
        'product_specs_fpath',
        help='path to product file specs'
        )
    parser.add_argument(
        'product_full_version',
        help='product full version string'
        )    
    parser.add_argument(
        'aux_viirs_snpp_input_fpath', 
        help='provide path to Aux-VIIRS-snpp input file'
        )
    parser.add_argument(
        'aux_viirs_noaa20_input_fpath', 
        help='provide path to Aux-VIIRS-noaa20 input file'
        )
    parser.add_argument(
        'aux_amsr_nise_input_fpath', 
        help='provide path to Aux-AMSR_NISE input file'
        )
    parser.add_argument(
        'aux_sfctype_final_input_fpath', 
        help='provide path to Aux-sfctype-final input file'
        )
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

    return (args.out_dir, args.auxmet_fpath, args.product_specs_fpath,
            args.product_full_version, args.aux_viirs_snpp_input_fpath,
            args.aux_viirs_noaa20_input_fpath, args.aux_amsr_nise_input_fpath,
            args.aux_sfctype_final_input_fpath, artp)


def main():
    out_dir, auxmet_fpath, product_specs_fpath, product_full_version, \
         aux_viirs_snpp_input_fpath, aux_viirs_noaa20_input_fpath, \
         aux_amsr_nise_input_fpath, aux_sfctype_final_input_fpath, \
         artp = parse_args()
    
    # Output data dictionary
    dat = {}
    
    # -- Populate Aux-Sat output data and group attributes --
    
    # Load Aux-Sat data from each set of input products
    # - This loop assumes that all of the Aux-Sat temporary file paths are
    #   passed to this script in the specific order of the input_files list above
    # - Also assumes that the exact same list of output variables is present in
    #   the Aux-VIIRS-SNPP and Aux-VIIRS-NOAA20 files
    input_files = [
        aux_viirs_snpp_input_fpath,
        aux_viirs_noaa20_input_fpath,
        aux_amsr_nise_input_fpath,
        aux_sfctype_final_input_fpath
        ]

    auxsat_data = {}
    auxsat_attrs = {}
    for ix,f in enumerate(input_files):
        result = glob.glob(f)
        if len(result) == 1:
            fn = result[0]
        elif len(result) > 1:
            raise ValueError(
                   f"Two or more files match {f}. There must only be 1 match.")
        else:
            continue  # No matches; go on to the next possible input file(s)

        with n.Dataset(fn, "r") as nc:
            variables = list(nc.variables.keys())
            shapes = np.shape(nc[variables[0]][...])
            
            # SNPP and NOAA-20 Aux data (VIIRS L3)
            if ('snpp' in fn) or ('noaa20' in fn):
                # Create SNPP and NOAA-20 output variables using the list of
                # variables from SNPP
                
                # This assumes a specific order of input files in "input_files"
                # list above: aux_viirs_snpp_input_fpath, aux_viirs_noaa20_input_fpath,
                # aux_amsr_nise_input_fpath, aux_iceshelves_input_fpath
                if ix == 0:
                    for var in variables:
                        # VIIRS L3 npts variable doesn't need third dimension
                        # because L3 grid is the same for SNPP and NOAA-20
                        if var == 'VIIRS_L3_num_pts_total':
                            auxsat_data[var] = np.full([shapes[0], shapes[1]], -9999.0)
                        else:
                            auxsat_data[var] = np.full([shapes[0], shapes[1], 2], -9999.0)
                        
                for var in variables:
                    # VIIRS L3 npts variable doesn't need third dimension
                    # because L3 grid is the same for SNPP and NOAA-20
                    if var == 'VIIRS_L3_num_pts_total':
                        # SNPP (first index of third array dimension)
                        if ix == 0:
                            auxsat_data[var] = nc[var][...]
                    else:
                        # SNPP (first index of third array dimension)
                        if ix == 0:
                            auxsat_data[var][:,:,0] = nc[var][...]
                        # NOAA-20 (second index of third array dimension)
                        elif ix == 1:
                            auxsat_data[var][:,:,1] = nc[var][...]

            # AMSR / NISE and final surface type Aux data
            else:
                for var in variables:
                    auxsat_data[var] = nc[var][...]
                    
            # Get netCDF attributes containing metadata for VIIRS, AMSR, and
            # NISE products
            # Skip if attributes are missing
            if len(nc.ncattrs()) > 0:
                for attr in nc.ncattrs():
                    auxsat_attrs[attr] = nc.getncattr(attr)
    
    dat["Aux-Sat"] = auxsat_data
    dat["Aux-Sat_Group_Attributes"] = auxsat_attrs
    
    auxmet_ds = n.Dataset(auxmet_fpath)
    
    # -- Populate Geometry group attributes --
    
    # Load "Geometry" group and its group attributes from the Aux-Met file:
    dat["Geometry"] = load_all_vars_of_nc4group("Geometry", auxmet_ds, artp)
    dat["Geometry_Group_Attributes"] = load_all_atts_of_nc4group(
        "Geometry", auxmet_ds)

    # -- Populate global attributes --
    
    global_atts = {}
    
    # Copied directly from Aux-Met file
    global_atts["granule_ID"] = auxmet_ds.granule_ID
    global_atts["spacecraft_ID"] = auxmet_ds.spacecraft_ID
    global_atts["sensor_ID"] = auxmet_ds.sensor_ID
    global_atts["ctime_coverage_start_s"] = auxmet_ds.ctime_coverage_start_s
    global_atts["ctime_coverage_end_s"] = auxmet_ds.ctime_coverage_end_s
    global_atts["UTC_coverage_start"] = auxmet_ds.UTC_coverage_start
    global_atts["UTC_coverage_end"] = auxmet_ds.UTC_coverage_end
    global_atts["orbit_sim_version"] = auxmet_ds.orbit_sim_version
    
    # Provenance of packages used
    with open(AUX_SAT_fpaths.scipkg_prdgitv_fpath, 'r') as in_f:
        line_parts = in_f.readline().split('(', maxsplit=1)
        global_atts["provenance"] = "{}{} ( {}".format(line_parts[0],
                                                      product_full_version,
                                                      line_parts[1].strip())
    # Aux-Sat algorithm version
    with open(AUX_SAT_fpaths.scipkg_version_fpath) as f:
        global_atts["processing_algorithmID"] = f.readline().strip()

    # Product and netCDF version IDs
    global_atts["full_versionID"] = product_full_version
    global_atts["archival_versionID"] = (
                            product_full_version.split('_')[0].replace('R', ''))
    global_atts["netCDF_lib_version"] = n.getlibversion().split()[0]
    
    # Input product file names
    auxmet_fn = os.path.basename(auxmet_fpath)
    in_file_l = [auxmet_fn]
    global_atts["input_product_files"] = in_file_l[0]
    
    # Generate Aux-Sat output file name:
    atdim_full = auxmet_ds.dimensions[artp[0]].size
    tokens = auxmet_fn.split('_')
    fname_tmp = "PREFIRE_SAT{}_AUX-SAT_{}_{}_{}.nc".format(
        global_atts["spacecraft_ID"][-1], global_atts["full_versionID"],
        tokens[5], global_atts["granule_ID"]
        )
    if (artp[1] == 0) and (artp[2] is None or artp[2] == atdim_full-1):
        AUXSAT_fname = "raw-"+fname_tmp
    else:
        if artp[2] is None:
            tmp_idx = atdim_full-1
        else:
            tmp_idx = artp[2]-1
        AUXSAT_fname = "raw-"+fname_tmp[:-3]+ \
              f"-{artp[0]}_{artp[1]:05d}_{tmp_idx:05d}_of_{atdim_full:05d}f.nc"
    global_atts["file_name"] = AUXSAT_fname
    
    # Populate file creation time
    now_UTC_DT = datetime.datetime.now(datetime.timezone.utc)
    global_atts["UTC_of_file_creation"] = now_UTC_DT.strftime(
        "%Y-%m-%dT%H:%M:%S.%f")
    
    dat["Global_Attributes"] = global_atts
    
    auxmet_ds.close()
        
    # -- Use generic PREFIRE product writer to produce the AUX-SAT output file --
    
    AUXSAT_fpath = os.path.join(out_dir, AUXSAT_fname)
    mkdir_p(os.path.dirname(AUXSAT_fpath))
    write_data_fromspec(dat, AUXSAT_fpath, product_specs_fpath, verbose=True)


if __name__ == '__main__':
    main()
