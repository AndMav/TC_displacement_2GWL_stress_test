"""
Coarsening of the WorldPop population rasters from 30 arcsec (~1 km) to 300 arcsec (~9 km): every block of
10 x 10 cells is summed into one cell, so the population is conserved.

Input:  config.BASE_RAW_TIF, config.FUTURE_RAW_TIF   (30 arcsec)
Output: config.BASE_TIF, config.FUTURE_TIF           (300 arcsec)

Usage:
    python s0_exposure_coarsening.py
"""

import rioxarray  # noqa: F401  (adds the .rio accessor to xarray)
import xarray as xr

from config import BASE_RAW_TIF, FUTURE_RAW_TIF, BASE_TIF, FUTURE_TIF

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

COARSE = 10                     # 10 x 10 cells of 30 arcsec -> one cell of 300 arcsec

# Inputs: the rasters at 30 arcsec
INPUT_TIFS = {
    "base_2005": BASE_RAW_TIF,          # 2005 baseline
    "ssp3_2050": FUTURE_RAW_TIF,        # SSP3, 2050
}

# Outputs: the rasters at 300 arcsec
OUTPUT_TIFS = {
    "base_2005": BASE_TIF,              # 2005 baseline
    "ssp3_2050": FUTURE_TIF,            # SSP3, 2050
}


'''
 ----------------------------------------------------------------------------------------
     II. Coarsening and saving
 ----------------------------------------------------------------------------------------
'''

def main(input_path, output_path, coarse=COARSE):

    with xr.open_dataset(input_path, chunks="auto") as ds:
        ds = ds.squeeze()                                       # remove the band dimension
        ds = ds.coarsen({"x": coarse, "y": coarse}).sum()       # sum every block of coarse x coarse cells

        output_path.parent.mkdir(parents=True, exist_ok=True)
        ds.rio.to_raster(output_path)

    print(f"Done: {input_path.name} -> {output_path}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    for name in INPUT_TIFS:
        main(INPUT_TIFS[name], OUTPUT_TIFS[name])
