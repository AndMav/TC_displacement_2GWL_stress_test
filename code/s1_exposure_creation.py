"""
Creates the base (2005) and future (SSP3-2050) population exposure of every country from the WorldPop coarsened rasters (300 arcsec).

Steps:
    1. Clips both rasters with the country's outline
    2. Matches the grid cells of the two clips and keep the common ones (the 2050 raster is shifted by 0.00125 deg,
       so for some countries, e.g. ATG and DJI, one clip has an extra empty column)
    3. Removes the cells with zero population in both 2005 and 2050
    4. Assigns impact function 1 and save

Input:  ppp_2005_1km_Aggregated_300as.tif, SSP3_2050_300as.tif        (WorldPop, 300 arcsec)
Output: {country}_base_2005_300as.h5, {country}_SSP3_2050_300as.h5    (exposure of every country)
"""

import sys
import numpy as np
from climada.entity.exposures import Exposures
from climada.util.coordinates import get_country_geometries
from config import BASE_TIF, FUTURE_TIF, EXP_DIR, COUNTRIES

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

HAZ_TYPE = "TC"
IMPF_ID = 1                     # impact function of every exposure point
GRID_DEG = 300 / 3600           # grid step of the rasters (300 arcsec)


'''
 ----------------------------------------------------------------------------------------
     II. Exposures of one country
 ----------------------------------------------------------------------------------------
'''

def grid_cells(gdf, lon0, lat0):
    """Column and row of every cell on the grid, counted from (lon0, lat0)."""
    col = np.round((gdf.geometry.x.values - lon0) / GRID_DEG).astype(int)
    row = np.round((gdf.geometry.y.values - lat0) / GRID_DEG).astype(int)
    return list(zip(col, row))


def exposures(country):
    """Base and future exposure of one country, with the same cells in the same order."""

    # 1. Clips both rasters with the country's outline
    country_shapes = get_country_geometries([country]).geometry.to_list()
    assert len(country_shapes) > 0, f"No geometry found for country '{country}'"
    exp_base = Exposures.from_raster(BASE_TIF, geometry=country_shapes)
    exp_future = Exposures.from_raster(FUTURE_TIF, geometry=country_shapes)

    # 2. Matchs the grid cells and keep the common ones
    lon0, lat0 = exp_base.gdf.geometry.x.min(), exp_base.gdf.geometry.y.min()
    cells_base = grid_cells(exp_base.gdf, lon0, lat0)
    cells_future = grid_cells(exp_future.gdf, lon0, lat0)
    common = set(cells_base) & set(cells_future)
    exp_base.gdf.drop(index=exp_base.gdf.index[[c not in common for c in cells_base]], inplace=True)
    exp_future.gdf.drop(index=exp_future.gdf.index[[c not in common for c in cells_future]], inplace=True)
    assert [c for c in cells_base if c in common] == [c for c in cells_future if c in common], \
        f"{country}: base and future do not have the same cells in the same order"

    # 3. Removes the cells with zero population in both years
    keep = (np.nan_to_num(exp_base.gdf["value"].values) > 0) | (np.nan_to_num(exp_future.gdf["value"].values) > 0)
    exp_base.gdf.drop(index=exp_base.gdf.index[~keep], inplace=True)
    exp_future.gdf.drop(index=exp_future.gdf.index[~keep], inplace=True)

    # 4. Assigns impact function
    exp_base.gdf[f"impf_{HAZ_TYPE}"] = IMPF_ID
    exp_future.gdf[f"impf_{HAZ_TYPE}"] = IMPF_ID

    return exp_base, exp_future


'''
 ----------------------------------------------------------------------------------------
     III. Computation and saving
 ----------------------------------------------------------------------------------------
'''

def main(country, out_dir=EXP_DIR):

    country_dir = out_dir / country
    country_dir.mkdir(parents=True, exist_ok=True)

    exp_base, exp_future = exposures(country)
    exp_base.write_hdf5(country_dir / f"{country}_base_2005_300as.h5")
    exp_future.write_hdf5(country_dir / f"{country}_SSP3_2050_300as.h5")

    print(f"Done: {country}, {len(exp_base.gdf)} points -> {country_dir}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    for country in (sys.argv[1:] or COUNTRIES):
        main(country)
