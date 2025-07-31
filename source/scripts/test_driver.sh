# ** This is a test driver script for a previous version of Aux-JPSS in which the
# nearest VIIRS granule to each TIRS scene was determined *after* all VIIRS granule
# indices had been matched to TIRS scenes.
# ** This script will not work with the current version of Aux-JPSS, but is left
# here for reference.

#!/bin/sh

# Print every command
set -v

JPSS_START="2021-04-18T21:12:00Z"
JPSS_END="2021-04-18T21:17:59Z"

PREFIRE_FILENAME=$(sips_search --sat "prefire1" --products "1B-RAD" --start "2021-04-19T00:00:00Z" --end "2021-04-19T00:05:00Z" --download . | rev | cut -d"/" -f1 | rev)
FILE_SPECS=/home/k/projects/PREFIRE_sim_tools/PREFIRE_sim_tools/datasim/file_specs.json

# Calculate intermediate/indices file for a given IMG/MOD SNPP granule
MOD_FILENAME=$(sips_search --sat "snpp" --products "VNP03MOD" --version "3.1.0" --start $JPSS_START --end $JPSS_END --download . | rev | cut -d"/" -f1 | rev)
IMG_FILENAME=$(sips_search --sat "snpp" --products "VNP03IMG" --version "3.1.0" --start $JPSS_START --end $JPSS_END --download . | rev | cut -d"/" -f1 | rev)
python3 calculate_ind.py PREFIRE_SAT1-snpp_R_S_20210418211200_now_00001.nc $PREFIRE_FILENAME $IMG_FILENAME $MOD_FILENAME

# Determine nearest JPSS granule for each TIRS footprint that has a match (for only one index file, all matches will be _that_ granule)
python3 calculate_nearest.py PREFIRE_SAT1-snpp-nearest_now_00001.npy PREFIRE_SAT1-snpp-tdiffs_now_00001.npy $PREFIRE_FILENAME PREFIRE_SAT1-snpp*.nc 

# Download auxiliary files
mkdir cldmsk snowcover mod_rad img_rad
sips_search --sat "snpp" --products "VNP02MOD" --start $JPSS_START --end $JPSS_END --download mod_rad/
sips_search --sat "snpp" --products "VNP02IMG" --start $JPSS_START --end $JPSS_END --download img_rad/
sips_search --sat "snpp" --products "CLDMSK_L2_VIIRS_SNPP" --start $JPSS_START --end $JPSS_END --download cldmsk/
sips_search --sat "snpp" --products "VNP10" --start $JPSS_START --end $JPSS_END --download snowcover/

# With nearest file and all the inputs, assemble the aux file
python3 calculate_aux.py PREFIRE_SAT1-snpp_now_00001.nc PREFIRE_SAT1-snpp-nearest_now_00001.npy PREFIRE_SAT1-snpp-tdiffs_now_00001.npy --index_files PREFIRE_SAT1-snpp_* --mod_geo_files VNP03MOD* --cloud_mask_files cldmsk/*_SNPP* --mod_rad_files mod_rad/VNP02MOD* --img_rad_files img_rad/VNP02IMG* --surf_refl_files img_rad/VNP02IMG* --snow_cover_files snowcover/VNP10* --surf_temp_files snowcover/VNP10*

# Repeat now for a NOAA20 granule
MOD_FILENAME=$(sips_search --sat "noaa20" --products "VJ103MOD" --version "3.1.0" --start $JPSS_START --end $JPSS_END --download . | rev | cut -d"/" -f1 | rev)
IMG_FILENAME=$(sips_search --sat "noaa20" --products "VJ103IMG" --version "3.1.0" --start $JPSS_START --end $JPSS_END --download . | rev | cut -d"/" -f1 | rev)
python3 calculate_ind.py PREFIRE_SAT1-noaa20_R_S_20210418211200_now_00001.nc $PREFIRE_FILENAME $IMG_FILENAME $MOD_FILENAME

python3 calculate_nearest.py PREFIRE_SAT1-noaa20-nearest_now_00001.npy PREFIRE_SAT1-noaa20-tdiffs_now_00001.npy $PREFIRE_FILENAME PREFIRE_SAT1-noaa20*.nc

sips_search --sat "noaa20" --products "VJ102MOD" --start $JPSS_START --end $JPSS_END --download mod_rad/
sips_search --sat "noaa20" --products "VJ102IMG" --start $JPSS_START --end $JPSS_END --download img_rad/
sips_search --sat "noaa20" --products "CLDMSK_L2_VIIRS_NOAA20" --start $JPSS_START --end $JPSS_END --download cldmsk/

python3 calculate_aux.py PREFIRE_SAT1-noaa20_now_00001.nc PREFIRE_SAT1-noaa20-nearest_now_00001.npy PREFIRE_SAT1-noaa20-tdiffs_now_00001.npy --index_files PREFIRE_SAT1-noaa20_* --mod_geo_files VJ103MOD* --cloud_mask_files cldmsk/*_NOAA20* --mod_rad_files mod_rad/VJ102MOD* --img_rad_files img_rad/VJ102IMG* --surf_refl_files img_rad/VJ102IMG* --snow_cover_files snowcover/VJ110* --surf_temp_files snowcover/VJ110*

# Combine SNPP and NOAA20 aux files into final AUX-JPSS file
python3 calculate_combine.py PREFIRE_SAT1-AUX-JPSS_R_S_now_00001.nc $FILE_SPECS $PREFIRE_FILENAME PREFIRE_SAT1-snpp_now_00001.nc PREFIRE_SAT1-noaa20_now_00001.nc
