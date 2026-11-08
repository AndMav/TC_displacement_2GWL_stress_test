"""
Settings shared by all scripts: data paths, country lists, GCMs and simulation scenarios.

All data are read from and written below DATA_ROOT, a folder in CLIMADA's data directory (SYSTEM_DIR, by default
~/climada/data). To use another location, change the CLIMADA data directory in climada.conf (local_data: system).
"""

from climada.util.constants import SYSTEM_DIR

# ---- Data root ---------------------------------------------------------------------------------------------------------
DATA_ROOT = SYSTEM_DIR / "tc_displacement"

# ---- External inputs: not created by the scripts, place them in these folders before running -------------------------
# WorldPop population rasters at 30 arcsec: 2005 (WorldPop) and SSP3 2050 (Wang et al. 2022) (s0)
WORLDPOP_RAW_DIR = DATA_ROOT / "exposure" / "worldpop_30as"
BASE_RAW_TIF = WORLDPOP_RAW_DIR / "ppp_2005_1km_Aggregated.tif"                 # 2005 baseline
FUTURE_RAW_TIF = WORLDPOP_RAW_DIR / "SSP3_2050_30as.tif"                        # SSP3, 2050

# MIT synthetic TC event sets per ocean basin (WindRiskTech, on request): TC_{basin}_0300as_MIT_{gcm}_{ssp}{cal}.hdf5
# (s2, s3)
BASIN_HAZ_DIR = DATA_ROOT / "hazard" / "basin_event_sets"

# ---- Exposure ----------------------------------------------------------------------------------------------------------
# WorldPop rasters at 300 arcsec (s0 -> s1)
WORLDPOP_DIR = DATA_ROOT / "exposure" / "worldpop_300as"
BASE_TIF = WORLDPOP_DIR / "ppp_2005_1km_Aggregated_300as.tif"                   # 2005 baseline
FUTURE_TIF = WORLDPOP_DIR / "SSP3_2050_300as.tif"                               # SSP3, 2050

# Country exposures: {country}/{country}_base_2005_300as.h5, {country}_SSP3_2050_300as.h5 (s1 -> s2, s4a, s5, s6c)
EXP_DIR = DATA_ROOT / "exposure" / "countries"

# ---- Hazard ------------------------------------------------------------------------------------------------------------
# Country hazards: {country}/TC_{country}_{gcm}_MIT_{20thcal|ssp370cal}.hdf5 (s3 -> s4a, s5)
HAZ_DIR = DATA_ROOT / "hazard" / "countries"

# ---- Vulnerability -----------------------------------------------------------------------------------------------------
# Calibrated (baseline, future) v_half pairs: {country}.csv (s2 -> s4a, s5)
V_HALF_DIR = DATA_ROOT / "vulnerability"

# ---- Simulation results ------------------------------------------------------------------------------------------------
BASE_RESULTS_DIR = DATA_ROOT / "results" / "base"                               # s4a -> s4b
BASE_SUMMARIES_CSV = BASE_RESULTS_DIR / "all_countries_base_summaries.csv"      # s4b -> s6b
DELTA_RESULTS_DIR = DATA_ROOT / "results" / "delta"                             # s5 -> s6a

# ---- Decomposition -----------------------------------------------------------------------------------------------------
DECOMPOSITION_DIR = DATA_ROOT / "decomposition" / "countries"                   # s6a -> s6b, s6c, s7a
WORLD_DATA_DIR = DATA_ROOT / "decomposition" / "world"                          # s7a -> s7b
DATA_S1_CSV = DATA_ROOT / "decomposition" / "DataS1.csv"                        # s8

# ---- Figures -----------------------------------------------------------------------------------------------------------
FIGURES_DIR = DATA_ROOT / "figures"                                             # s6b, s6c, s7b

# ---- Countries (ISO3) --------------------------------------------------------------------------------------------------
# The 93 countries with calibrated impact functions (Kam et al.)
COUNTRIES = [
    'ABW', 'AFG', 'ARE', 'ASM', 'ATG', 'AUS', 'BGD', 'BHR', 'BHS', 'BLZ',
    'BMU', 'BRB', 'BTN', 'CAN', 'CHN', 'COK', 'COL', 'COM', 'CPV', 'CRI',
    'CUB', 'CYM', 'DJI', 'DMA', 'DOM', 'ERI', 'ETH', 'FJI', 'GRD', 'GTM',
    'GUY', 'HKG', 'HND', 'HTI', 'IDN', 'IND', 'IRN', 'IRQ', 'JAM', 'JPN',
    'KHM', 'KNA', 'KOR', 'LAO', 'LCA', 'LKA', 'MDG', 'MDV', 'MEX', 'MMR',
    'MOZ', 'MSR', 'MUS', 'MWI', 'MYS', 'NCL', 'NIC', 'NIU', 'NPL', 'NZL',
    'OMN', 'PAK', 'PAN', 'PHL', 'PNG', 'PRI', 'PYF', 'QAT', 'SAU', 'SLB',
    'SLV', 'SOM', 'SUR', 'SWZ', 'TCA', 'THA', 'TLS', 'TON', 'TTO', 'TWN',
    'TZA', 'USA', 'VCT', 'VEN', 'VGB', 'VIR', 'VNM', 'VUT', 'WLF', 'WSM',
    'YEM', 'ZAF', 'ZWE',
]

# Countries of the country figures (s6b, s6c), top to bottom
FIGURE_COUNTRIES = ['PHL', 'VNM', 'CUB']
COUNTRY_NAMES = {'PHL': 'Philippines', 'VNM': 'Vietnam', 'CUB': 'Cuba'}

# ---- GCMs and simulation scenarios -------------------------------------------------------------------------------------
# The 9 CMIP6 models of the MIT event sets, in the order of gc_model 1..9
GCMS = ['cesm2', 'cnrm6', 'ecearth6', 'fgoals', 'ipsl6', 'miroc6', 'mpi6', 'mri6', 'ukmo6']

# Simulation scenarios
SIM_TYPES = ['all', 'haz', 'exp', 'vuln', 'haz_exp', 'haz_vuln', 'exp_vuln']
