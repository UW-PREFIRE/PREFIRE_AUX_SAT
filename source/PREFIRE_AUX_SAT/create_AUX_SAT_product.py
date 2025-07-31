import os
import argparse
from pathlib import Path
import subprocess
import shutil

_code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

jpss_sats = ['snpp','noaa20']

def create_AUX_SAT_product(
    auxmet_fpath, ancillary_dir, input_products_dir, output_dir,
    tmpfiles_dir, product_full_version, 
    atrack_range_to_process=("atrack",0,None)
    ):
    
    # Get paths to directories:
    ancillary_dirpath = str(Path(ancillary_dir))
    input_products_dirpath = str(Path(input_products_dir))
    output_dirpath = str(Path(output_dir))
    tmpfiles_dirpath = str(Path(tmpfiles_dir))

    # Make directory for temp files created during processing of this granule
    tmp_folder = Path(auxmet_fpath).stem
    tmpdir_granule = os.path.join(tmpfiles_dirpath, tmp_folder)
    if not os.path.exists(tmpdir_granule):
        os.mkdir(tmpdir_granule)

    # Get path to Aux-Sat file specs    
    product_specs_fpath = os.path.join(ancillary_dirpath,
                                       "Aux-Sat_product_filespecs.json")

    this_environ = os.environ.copy()

    artpd, artp_b, artp_e = atrack_range_to_process  # Numpy indexing
    if artp_e is None:
        artp_strrep = f"{artpd.lower()}:{artp_b:d}:END"
    else:
        artp_strrep = f"{artpd.lower()}:{artp_b:d}:{artp_e:d}"

    # Calculate Aux interpolated AMSR and NISE data
    cmd = ["python", os.path.join(_code_dir, "calculate_AMSR_NISE_aux.py"),
           auxmet_fpath, product_specs_fpath, input_products_dirpath,
           tmpdir_granule, artp_strrep]
    subprocess.run(cmd, env=this_environ)
    
    # Calculate Aux interpolated VIIRS L3 CMG snow cover
    for jpss_sat in jpss_sats:
        cmd = ["python", os.path.join(_code_dir, "calculate_VIIRS_L3_aux.py"),
               jpss_sat, auxmet_fpath, product_specs_fpath,
               input_products_dirpath, tmpdir_granule, artp_strrep]
        subprocess.run(cmd, env=this_environ)
    
    # Classify final surface type using L1B, Aux-Met, and Aux-Sat data
    cmd = ["python", os.path.join(_code_dir, "classify_sfc_type_final.py"),
           tmpdir_granule, auxmet_fpath, product_specs_fpath, artp_strrep,
           "--jpss_sats"]
    cmd.extend(jpss_sats)
    subprocess.run(cmd, env=this_environ)

    # Merge Aux-snpp and Aux-noaa20 VIIRS and ATMS data, AMSR and NISE sea ice
    # and snow cover data, and final surface type data into a combined Aux-Sat file.
    cmd = ["python", os.path.join(_code_dir, "combine_aux.py"),
           output_dirpath, auxmet_fpath, product_specs_fpath,
           product_full_version,
           os.path.join(tmpdir_granule, "PREFIRE*Aux-VIIRS-snpp*.nc"),
           os.path.join(tmpdir_granule, "PREFIRE*Aux-VIIRS-noaa20*.nc"),
           os.path.join(tmpdir_granule, "PREFIRE*Aux-AMSR_NISE*.nc"),
           os.path.join(tmpdir_granule, "PREFIRE*Aux-sfctype-final*.nc"),
           artp_strrep]
    subprocess.run(cmd, env=this_environ)

    # Delete all temporary files created by Aux-Sat processing
    shutil.rmtree(tmpdir_granule)
    

def parse_args():
    """
    Parse command line arguments passed to script at runtime.

    Returns
    -------
    args.auxmet_fpath : str
        Path to Aux-Met file for the given PREFIRE granule.
    args.ancillary_dir : str
        Path to directory containing ancillary files.
    args.input_products_dir : str
        Parent directory of input product files.
    args.output_dir : str
        Path to directory in which Aux-Sat file will be stored.

    """
    parser = argparse.ArgumentParser('Creates the PREFIRE Aux-Sat product(s)')
    parser.add_argument(
        'auxmet_fpath', 
        help='path to Aux-Met file for the given PREFIRE granule'
        )
    parser.add_argument(
        'ancillary_dir', 
        help='path to directory containing ancillary files'
        )
    parser.add_argument(
        'input_products_dir', 
        help='path to input products parent directory'
        )
    parser.add_argument(
        'output_dir', 
        help='path to parent directory in which Aux-Sat file will be stored'
        )
    parser.add_argument(
        'tmpfiles_dir', 
        help='path to directory in which temporary files will be stored'
        )
    parser.add_argument(
        "-i", "--atrack-idx-range",
        metavar="atrack_idx_range", default="0:END",
        help="along-track index range (zero-based, inclusive) "
             "to process within granule. Valid examples: "
             "'0:8140', '456:1200', '0:END' "
             "(default: %(default)s."
         )
    args = parser.parse_args()
    
    tokens = args.atrack_idx_range.split(':')
    if tokens[1] == "END":
        atrack_idx_range_np = ("atrack", int(tokens[0]), None)  # Numpy indexing
    else:
        atrack_idx_range_np = ("atrack", int(tokens[0]), int(tokens[1])+1)  # Numpy indexing

    return args.auxmet_fpath, args.ancillary_dir, args.input_products_dir, \
        args.output_dir, args.tmpfiles_dir, \
        atrack_idx_range_np
        

def main(
    auxmet_fpath, ancillary_dir, input_products_dir,
    output_dir, tmpfiles_dir, atrack_idx_range
    ):
    
    create_aux_sat_product(
        auxmet_fpath, ancillary_dir, input_products_dir,
        output_dir, tmpfiles_dir, atrack_range_to_process=atrack_idx_range
        )


if __name__ == '__main__':
    main(*parse_args())
