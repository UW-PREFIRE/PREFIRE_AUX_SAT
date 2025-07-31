#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./refresh_test_IO_files.sh   OR   bash refresh_test_IO_files.sh

# Recreate/refresh all test input and output files for this science algorithm
#  package, EXCEPT for the test input files required to bootstrap the process. 

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

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package dist/ directory)
readonly base_dir="$(absfpath ".")";

this_top_dir="$(absfpath "${base_dir}/..")";

# Change dir to the 'test/' directory:
test_path="${this_top_dir}/test";
cd "$test_path";

# Set input path:
inputs_path="$test_path"/inputs;

# Remove old AUX-SAT output files and create new ones:
outputs_path="$test_path"/outputs;
rm -rf "$outputs_path";
mkdir -p "$outputs_path";
./run.sh;

echo "Finished refreshing test input.";
