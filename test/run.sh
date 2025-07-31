#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run.sh    OR    bash run.sh

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run.sh; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

absfpath() {
  # Generate absolute filepath from a relative (or even an absolute) filepath.
  #
  # Based on (circa Oct 2023) https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
  # 
  # $1     : a relative (or even an absolute) filepath
  # Returns the corresponding absolute filepath.
  if [ -d "$1" ]; then
    # dir
    (cd "$1"; pwd)
  elif [ -f "$1" ]; then
    # file
    if [[ $1 = /* ]]; then
      echo "$1"
    elif [[ $1 == */* ]]; then
      echo "$(cd "${1%/*}"; pwd)/${1##*/}"
    else
      echo "$(pwd)/$1"
    fi
  fi
}

activate_conda_env () {
  . "$1"/bin/activate;
}

deactivate_conda_env () {
  . "$1"/bin/deactivate;
}

set -ve;  # Exit on the first error, and print out commands as we execute them
#set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
readonly base_dir="$(absfpath ".")";

hn=`hostname -s`;  # Hostname

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

non_SDPS_hostname="longwave";

# For tests on longwave, the input, output, and tmpfile directories are symbolic
#  links to directories in /data (to avoid placing large files in /home).

input_dir="${base_dir}/inputs";

cfg_str1="ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END,${input_dir}/PREFIRE_SAT1_AUX-MET_R01_P00_20241007075724_01877.nc";
cfg_str2="ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END,${input_dir}/PREFIRE_SAT2_AUX-MET_R01_P00_20241007071543_02040.nc";


# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
MKL_NUM_THREADS=1;
NUMEXPR_NUM_THREADS=1;
OMP_NUM_THREADS=1;
VECLIB_MAXIMUM_THREADS=1;
OPENBLAS_NUM_THREADS=1;
export MKL_NUM_THREADS NUMEXPR_NUM_THREADS OMP_NUM_THREADS;
export VECLIB_MAXIMUM_THREADS OPENBLAS_NUM_THREADS;

# Some environment vars that convey configuration info to the algorithm:

this_top_dir="$(absfpath "${base_dir}/..")";

PACKAGE_TOP_DIR="${this_top_dir}";
ANCILLARY_DATA_DIR="${this_top_dir}/dist/ancillary";

OUTPUT_DIR=${base_dir}/outputs;
TMP_DIR=${base_dir}/outputs;

  # Path to directory containing input non-PREFIRE satellite products (e.g.,
  # from VIIRS, AMSR, NISE):
IN_PRODS_DIR=${input_dir};

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
PRODUCT_FULLVER="R01_P00";
  # Special form ('R00_Syy') when processing simulated observations:
#PRODUCT_FULLVER="R00_S01";

# Make required environment vars available:
export PACKAGE_TOP_DIR ANCILLARY_DATA_DIR OUTPUT_DIR TMP_DIR IN_PRODS_DIR;
export PRODUCT_FULLVER;

# Check if output and temporary file directories exist; if not, bail:
tmpdir="${OUTPUT_DIR}";
test -d "${tmpdir}" || { echo "Output directory does not exist: ${tmpdir}"; exit 1; }
tmpdir="${TMP_DIR}";
test -d "${tmpdir}" || { echo "Temporary files directory does not exist: ${tmpdir}"; exit 1; }

# If custom conda environment files exist, activate that conda environment:
conda_env_dir="${this_top_dir}/dist/c_env_for_PREFIRE_AUX_SAT";
if [ -d "${conda_env_dir}" ]; then
   activate_conda_env "${conda_env_dir}";
fi

# Execute script that writes a new 'prdgit_version.txt', which contains
#  product moniker(s) and current (latest) git hash(es) that are part of the
#  provenance of this package's product(s).
# *** This step should not be done within the SDPS, since that file is
#     created just before delivery to the SDPS.
if [ ! -f "${this_top_dir}/dist/for_SDPS_delivery.txt" ]; then
   python "${this_top_dir}/dist/determine_prdgit.py";
fi

#for cfg_str in ${cfg_str1} ${cfg_str2}
for cfg_str in ${cfg_str1}
do
   ATRACK_IDX_RANGE_0BI=${cfg_str%,*};
   AUX_MET_FILE=${cfg_str##*,};

   export AUX_MET_FILE ATRACK_IDX_RANGE_0BI;

   # Execute primary driver:
   python "${this_top_dir}/dist/produce_AUX.py";
done

# If custom conda environment files exist, DEactivate that conda environment:
if [ -d "${conda_env_dir}" ]; then
   deactivate_conda_env "${conda_env_dir}";
fi
