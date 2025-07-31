"""
Produce AUX-SAT product using AUX-MET and other non-PREFIRE satellite products.

This program requires python version 3.6 or later, and is importable as a 
python module.
"""

  # From the Python standard library:
from pathlib import Path
import os
import sys
import argparse

  # From other external Python packages:

  # Custom utilities:


#--------------------------------------------------------------------------
def main():
    """Driver routine."""

    package_top_Path = Path(os.environ["PACKAGE_TOP_DIR"])

    src_fpath = str(package_top_Path / "source")
    sys.path.append(src_fpath)
    from PREFIRE_AUX_SAT.create_AUX_SAT_product import create_AUX_SAT_product

    # This is needed for explicit script calls later:
    if "PYTHONPATH" in os.environ:
        tmp_str = os.environ["PYTHONPATH"]
        if src_fpath not in tmp_str:
            os.environ["PYTHONPATH"] = tmp_str+f":{src_fpath}"
    else:
        os.environ["PYTHONPATH"] = src_fpath

    ancillary_dir = os.environ["ANCILLARY_DATA_DIR"]
    input_AUXMET_fpath = os.environ["AUX_MET_FILE"]
    output_dir = os.environ["OUTPUT_DIR"]
    tmpfiles_dir = os.environ["TMP_DIR"]
    input_prods_dir = os.environ["IN_PRODS_DIR"]

    # If the following subdirectories do not exist in the input satellite
    #  products directory, temporarily mimic them with a symlink to
    #  'input_prods_dir' -- so that the files can be found properly:
    input_prods_subdir_l = ["AMSR_U2_L3_SeaIce12km_B04", "NISE_SSMISF18",
                            "VJ110C1", "VNP10C1"]
    for subdir in input_prods_subdir_l:
        tmp_path = os.path.join(input_prods_dir, subdir)
        if not os.path.exists(tmp_path):
            os.symlink(input_prods_dir, tmp_path)

    atrack_idx_range_str = os.environ["ATRACK_IDX_RANGE_0BI"]
    tokens = atrack_idx_range_str.split(':')
    if tokens[2] == "END":
        atrack_np_idx_range = ("atrack", int(tokens[1]), None)  # Numpy indexing
    else:
        atrack_np_idx_range = ("atrack", int(tokens[1]),
                               int(tokens[2])+1)  # Numpy indexing

      # Default product_fullver:
    if "PRODUCT_FULLVER" not in os.environ:
        product_full_version = "R01_P00"
    elif len(os.environ["PRODUCT_FULLVER"].strip()) == 0:
        product_full_version = "R01_P00"
    else:
        product_full_version = os.environ["PRODUCT_FULLVER"]

    # Create the product data:
    create_AUX_SAT_product(input_AUXMET_fpath, ancillary_dir, input_prods_dir,
                           output_dir, tmpfiles_dir, product_full_version,
                           atrack_range_to_process=atrack_np_idx_range)


if __name__ == "__main__":
    # Determine fully-qualified filesystem location of this script:
    anchor_path = os.path.abspath(os.path.dirname(sys.argv[0]))

    # Process arguments:
    arg_description = ("Produce AUX-SAT product using AUX-MET and other "
                       "non-PREFIRE satellite products.")
    arg_parser = argparse.ArgumentParser(description=arg_description)

    args = arg_parser.parse_args()

    # Run driver:
    main()
