"""
Creates the tropical cyclone hazard of every country: the MIT synthetic TC event sets of the 9 GCMs for the
baseline (1995-2014) and SSP3-7.0 (2041-2060), cut to the country. The four basins are concatenated first, as the
tracks of one basin also reach into its neighbours.

Input:  TC_{basin}_0300as_MIT_{gcm}_{ssp}{cal}.hdf5       (TC event sets of the four basins)
Output: TC_{country}_{gcm}_MIT_20thcal.hdf5               (baseline hazard of every country and GCM)
        TC_{country}_{gcm}_MIT_ssp370cal.hdf5             (future hazard of every country and GCM)
        (folders defined in config.py)
"""

import sys

import numpy as np

from climada.hazard import Hazard
from climada.util.coordinates import get_country_geometries, country_to_iso, latlon_bounds

from config import BASIN_HAZ_DIR, HAZ_DIR, COUNTRIES, GCMS

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

BASINS = ("WP", "SH", "IO", "AP")      # all four basins, in this order (sets the event IDs)

# Periods: name in the output file -> (SSP, 'cal' code) in the basin file names
PERIODS = {
    "20thcal": ("", "20thcal"),          # baseline 1995-2014 (no SSP)
    "ssp370cal": ("ssp370", "cal"),      # SSP3-7.0, 2041-2060
}

ARCSECS = "300"                          # resolution of the event sets
TCM = "MIT"                              # TC model

BBOX_BUFFER_DEG = 0.2                    # buffer around the country's bounding box (degrees)
REG_ID_COUNTRIES = ["FJI", "USA", "NZL"] # selected by region ID: territory across the 180 deg meridian

REMOVE_NULL_EVENTS = True                # remove the events without any intensity in the country
RENAME_EVENTS = True                     # unique event names '{gcm}-{ssp}{cal}-{tcm}-{index}'

OVERWRITE = False                        # False: skip a file that already exists


'''
 ----------------------------------------------------------------------------------------
     II. Event sets of all basins
 ----------------------------------------------------------------------------------------
'''

def basin_path(basin, gcm, ssp, cal):
    """Path of the event set of one basin, GCM and period."""
    return BASIN_HAZ_DIR / f"TC_{basin}_0{ARCSECS}as_{TCM}_{gcm}_{ssp}{cal}.hdf5"


def all_basins(gcm, ssp, cal):
    """The event sets of the four basins for one GCM and period, concatenated (WP, SH, IO, AP)."""
    hazard = Hazard.concat([Hazard.from_hdf5(basin_path(basin, gcm, ssp, cal)) for basin in BASINS])

    # Unique event names (the basin files use the same names)
    if RENAME_EVENTS:
        hazard.event_name = [f"{gcm}-{ssp}{cal}-{TCM}-{idx}" for idx in range(hazard.size)]
    return hazard


'''
 ----------------------------------------------------------------------------------------
     III. Hazard of one country
 ----------------------------------------------------------------------------------------
'''

def country_hazard(hazard, country):
    """The hazard restricted to one country: the centroids in its bounding box (+ buffer), or for FJI / USA / NZL the
    centroids of the country by region ID; then without the events that have no intensity there."""

    if country in REG_ID_COUNTRIES:
        # Centroids of the country by region ID (territory across the 180 deg meridian)
        selected = hazard.select(reg_id=[country_to_iso(country, "numeric")])
    else:
        # Centroids in the country's bounding box, plus the buffer
        lon_min, lat_min, lon_max, lat_max = get_country_geometries([country]).geometry.total_bounds
        lon_min, lat_min, lon_max, lat_max = latlon_bounds(
            lat=np.array([lat_min, lat_max]), lon=np.array([lon_min, lon_max]), buffer=BBOX_BUFFER_DEG)
        selected = hazard.select(extent=(lon_min, lon_max, lat_min, lat_max))
    assert selected is not None, f"{country} not in hazard"

    # Removes the events without intensity in the country
    if REMOVE_NULL_EVENTS:
        nonzero_rows = np.unique(selected.intensity.nonzero()[0])
        selected = selected.select(event_id=selected.event_id[nonzero_rows].tolist())
        assert selected is not None, f"{country} has no hazard events"

    return selected


'''
 ----------------------------------------------------------------------------------------
     IV. Computation and saving
 ----------------------------------------------------------------------------------------
'''

def hazard_path(country, gcm, period):
    """Path of the hazard of one country, GCM and period (as read by s4a / s5)."""
    return HAZ_DIR / country / f"TC_{country}_{gcm}_{TCM}_{period}.hdf5"


def main(countries):
    """All files of the given countries. The basin event sets are loaded once per GCM and period, then every country
    is cut out of them. A country that fails (e.g. no hazard) is reported and skipped; returns the failed ones."""
    failed = set()
    for gcm in GCMS:
        for period, (ssp, cal) in PERIODS.items():

            todo = [c for c in countries if OVERWRITE or not hazard_path(c, gcm, period).exists()]
            if not todo:
                print(f"Skipped: {gcm} {period}, all files exist")
                continue

            hazard = all_basins(gcm, ssp, cal)
            for country in todo:
                path = hazard_path(country, gcm, period)
                try:
                    selected = country_hazard(hazard, country)
                except Exception as err:                        # report the error, continue with the next country
                    print(f"FAILED: {country} {gcm} {period} - {err}")
                    failed.add(country)
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                selected.write_hdf5(path)
                print(f"Done: {country} {gcm} {period}: {selected.size} events, "
                      f"{selected.centroids.size} centroids -> {path.name}")

    if failed:
        print(f"Failed countries: {sorted(failed)}")
    return sorted(failed)

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    countries = sys.argv[1:] or COUNTRIES
    main(countries)
