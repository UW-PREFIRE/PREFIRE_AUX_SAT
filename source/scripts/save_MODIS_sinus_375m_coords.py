"""
One-time script to transform the x/y projection coordinates of the 375m MODIS
Sinusoidal Tile Grid to geographic lat/lon coordinates. These geographic
coordinates are saved to disk, so that the Aux-JPSS processing code doesn't
need to perform the time-consuming coordinate transformation during each
production run.

For each tile on the 375m MODIS Sinusoidal Tile Grid, save two sets of coordinates:
    - 2D longitude and latitude arrays containing *all* lat/lon points in the tile.
      These are used to create the tile boundary polygons for initial screening
      of TIRS scene overlaps.
    - 1D flattened longitude and latitude arrays containing only the "non-fill"
      points, i.e. coordinates from the section of each tile that is actually
      within the global map area. These are used to create the STRtree lookup
      of L3 indices for pixel-level matching with TIRS scenes.
"""

import numpy as np
import glob
import pyproj
import netCDF4 as n
import datetime as dt
import os

L3_snow_cover_dir = '/data/users/k/VIIRS_products_allfiles_for_2021-01-01_test/VNP10A1F/'
coords_output_dir = '/data/users/k/MODIS_sinus_375m_coords/'

# MODIS Sinusoidal Tile Grid projection parameters and Earth radius are from VNP10A1F user's guide
modis_sinus_proj4 = '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +R=6371007.181000 +units=m +no_defs'
modis_sinus_crs = pyproj.CRS.from_proj4(modis_sinus_proj4)
# EPSG 4326 = global lat-lon coordinate system
target_crs = 'epsg:4326'
sinus_to_latlon = pyproj.Transformer.from_crs(modis_sinus_crs, target_crs)

# Only use L3 files from one day, since MODIS Sinusoidal Tile Grid is static
# (doesn't change from day to day)
L3_fpaths = sorted(glob.glob(L3_snow_cover_dir+'VNP10A1F.A2021001*.h5'))

for fpath in L3_fpaths:
    tile_str = os.path.basename(fpath).split('.')[2]
    
    with n.Dataset(fpath, 'r') as nc:
        xs = nc['HDFEOS/GRIDS/NPP_Grid_IMG_2D/XDim'][:]
        ys = nc['HDFEOS/GRIDS/NPP_Grid_IMG_2D/YDim'][:]
        sc = nc['HDFEOS/GRIDS/NPP_Grid_IMG_2D/Data Fields/CGF_NDSI_Snow_Cover'][:].data

    X,Y = np.meshgrid(xs,ys)
    lat2d, lon2d = sinus_to_latlon.transform(X,Y)
            
    lat_nofill = lat2d[sc != 255]
    lon_nofill = lon2d[sc != 255]
    
    np.save(f'{coords_output_dir}MODIS_sinus_375m_lat2d_{tile_str}.npy', lat2d, allow_pickle=True)
    np.save(f'{coords_output_dir}MODIS_sinus_375m_lon2d_{tile_str}.npy', lon2d, allow_pickle=True)

    np.save(f'{coords_output_dir}MODIS_sinus_375m_lat1d_validpts_{tile_str}.npy', lat_nofill, allow_pickle=True)
    np.save(f'{coords_output_dir}MODIS_sinus_375m_lon1d_validpts_{tile_str}.npy', lon_nofill, allow_pickle=True)

    now = dt.datetime.now()
    print(f'Wrote coords for {tile_str} at {now}')
