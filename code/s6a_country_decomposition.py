"""
Decomposes the change in TC displacement risk of a country into the contributions of its drivers (hazard, exposure,
vulnerability) and their interactions, from the results of the 7 simulation scenarios (s5), at p10, median and p90.
Capital letters: decomposition terms; small letters: scenario results (relative change).

    HAZ            = haz
    EXP            = exp
    VULN           = vuln
    HAZxEXPint     = haz_exp - haz - exp
    HAZxVULNint    = haz_vuln - haz - vuln
    EXPxVULNint    = exp_vuln - exp - vuln
    HAZxEXPxVULNint = all - haz_exp - haz_vuln - exp_vuln + haz + exp + vuln

Input:  unsequa_delta_TC_{country}_{sim_type}_MIT.hdf5    (results of the 7 scenarios)
Output: {country}_decomposition_AAD.csv                   (country decomposition of AAD)
        {country}_decomposition_rp100.csv                 (country decomposition of RP100)
        {country}_EAD_exp_point_decomposition.csv         (decomposition of EAD at every exposure point, with the
                                                               dominant term)
        (folders defined in config.py)
"""


import sys
import numpy as np
import pandas as pd
from climada.engine.unsequa import UncDeltaImpactOutput
from config import DELTA_RESULTS_DIR, DECOMPOSITION_DIR, COUNTRIES, SIM_TYPES


'''
 ----------------------------------------------------------------------------------------
     I. Constants
 ----------------------------------------------------------------------------------------
'''

METRICS = ['aai_agg', 'rp100']                  # metrics of the CLIMADA uncertainty output
METRIC_NAMES = {'aai_agg': 'AAD', 'rp100': 'rp100'}   # their names in the output files
STATS = ['median', 'p10', 'p90']        

DECOMP_TERMS = ['ALL', 'HAZ', 'EXP', 'VULN', 'HAZxEXPint', 'HAZxVULNint', 'EXPxVULNint', 'HAZxEXPxVULNint']

'''
 ----------------------------------------------------------------------------------------
     II. Helper Functions
 ----------------------------------------------------------------------------------------
'''
 
def hdf5_path(country, sim_type):
    """Path of the results of one scenario of the country (s5)."""
    
    return DELTA_RESULTS_DIR / country / sim_type / f"unsequa_delta_TC_{country}_{sim_type}_MIT.hdf5"
 
 
def load_results(country):
    """Builds dictionary with the results of each simulation for the country, loading ONE scenario at a time and keeping
    only what the decomposition needs, so that memory holds one scenario's samples at a time (large countries):
        'unc':    the uncertainty table (one row per sample; aai_agg, rp100, ...)
        'eai':    the statistics of eai_exp (EAD) at every exposure point {stat: one value per point}
        'coords': the coordinates of the exposure points (latitude, longitude)
    The statistics are the same as when computed from all scenarios loaded together."""

    results = {}
    for sim_type in SIM_TYPES:
        output = UncDeltaImpactOutput.from_hdf5(hdf5_path(country, sim_type))      # one scenario
        samples = output.get_unc_df('eai_exp').values                              # samples x exposure points
        results[sim_type] = {
            'unc': output.get_uncertainty(),
            'eai': {stat: extract_stat(samples, stat, axis=0) for stat in STATS},   # one value per point
            'coords': output.coord_df[['latitude', 'longitude']].reset_index(drop=True),
        }
        del output, samples                                                        # free memory before the next one
    return results
 
def extract_stat(values, stat, axis=None):
    """Median, 10th or 90th percentile of the samples, ignoring NaN.
    Returns: 
        - One value: samples of one metric (1D array)                  -> one number          (axis=None)
        - Multiple Values: samples x exposure points (2D array)        -> one value per point (axis=0)"""

    if stat == 'median':
        return np.nanmedian(values, axis=axis)

    elif stat == 'p10':
        return np.nanpercentile(values, 10, axis=axis)

    elif stat == 'p90':
        return np.nanpercentile(values, 90, axis=axis)
 
 
def decompose(summary):
    """Splits the total change summary['all'] into pure effects and interactions.
    summary: the value of every scenario for one statistic (a number, or an array with one value per exposure point).
    Returns: ALL and the 7 terms (decomposition formulas at the top)."""

    # Scenario results
    ALL = summary['all']
    HAZ = summary['haz']
    EXP = summary['exp']
    VULN = summary['vuln']
    HAZxEXP = summary['haz_exp']
    HAZxVULN = summary['haz_vuln']
    EXPxVULN = summary['exp_vuln']

    # Interaction terms
    HAZxEXP_int = HAZxEXP - HAZ - EXP
    HAZxVULN_int = HAZxVULN - HAZ - VULN
    EXPxVULN_int = EXPxVULN - EXP - VULN
    HAZxEXPxVULN_int = (ALL - HAZxEXP - HAZxVULN - EXPxVULN + HAZ + EXP + VULN)
    
    return ALL, HAZ, EXP, VULN, HAZxEXP_int, HAZxVULN_int, EXPxVULN_int, HAZxEXPxVULN_int


'''
 ----------------------------------------------------------------------------------------
     III. Decomposition of AAD, RP100 and EAD
 ----------------------------------------------------------------------------------------
'''


def country_decomposition(unc_outputs):
    """Decomposition of AAD and RP100 for every statistic (p10, median, p90).
       Output: one table with the columns decomp_term, delta_change, metric, stat
    """

    unc_results_dfs = {}

    # Uncertainty table of every scenario: one row per sample, one column per metric
    for sim_type in SIM_TYPES:
        unc_results_dfs[sim_type] = unc_outputs[sim_type]['unc']

    tables = []                         
    
    for stat in STATS:                           # median, p10, p90
        for metric in METRICS:                   # AAD, RP100

            stat_summary = {}

            # Statistic of every scenario, ignoring NaN
            for sim_type in SIM_TYPES:
                samples = unc_results_dfs[sim_type][metric].values
                stat_summary[sim_type] = extract_stat(samples, stat)

            # no decomposition if any scenario has no valid samples (NaN)
            if any(np.isnan(value) for value in stat_summary.values()):
                continue

            # Decomposition: ALL and the 7 terms
            decomp_values = decompose(stat_summary)                 

            df = pd.DataFrame({'decomp_term': DECOMP_TERMS, 'delta_change': decomp_values})
            df['metric'] = METRIC_NAMES[metric]
            df['stat'] = stat

            tables.append(df)

    # Nothing decomposed: empty table with the column names
    if not tables:
        return pd.DataFrame(columns=['decomp_term', 'delta_change', 'metric', 'stat'])

    return pd.concat(tables, ignore_index=True)      



def spatial_decomposition(unc_outputs):
    """Decomposition of EAD (eai_exp) at every exposure point for every statistic, with the dominant term of each point.

    The dominant driver of a point is the term with the largest absolute value (ALL does not take part),
    except where a pure effect (HAZ, EXP or VULN) is exactly -1: that driver alone takes the impact to zero
    and is the dominant driver (EXP before HAZ before VULN). Points without any valid term are NO_DATA.
    ALL and the terms are written as computed.
    All scenarios share the exposure points/coordinates of the 'all' scenario (same points, same order).
    Output: one table, one row per exposure point and statistic, with the columns
        latitude, longitude, stat, dominant_driver, dominant_value, ALL and the 7 terms
    """
 
 
    # Coordinates of the exposure points (same points, same order in every scenario)
    coords = unc_outputs['all']['coords']
 
    tables = []
 
    for stat in STATS:                           # median, p10, p90
 
        stat_summary = {}
 
        # Statistic of every scenario at every exposure point (computed in load_results)
        for sim_type in SIM_TYPES:
            stat_summary[sim_type] = unc_outputs[sim_type]['eai'][stat]
 
        # Decomposition at every exposure point: ALL and the 7 terms
        decomp_values = decompose(stat_summary)

        # Dominant term of every point: the largest absolute value (ALL does not take part); NO_DATA if all are NaN
        term_values = pd.DataFrame({
            'HAZ':             decomp_values[1],
            'EXP':             decomp_values[2],
            'VULN':            decomp_values[3],
            'HAZxEXPint':      decomp_values[4],
            'HAZxVULNint':     decomp_values[5],
            'EXPxVULNint':     decomp_values[6],
            'HAZxEXPxVULNint': decomp_values[7],
        })

        abs_term_values = term_values.abs()
        no_data = term_values.isna().all(axis=1)
        abs_term_values[no_data] = 0

        dominant_driver = abs_term_values.idxmax(axis=1)

        dominant_value = []                                      
        for row, term in dominant_driver.items():
            dominant_value.append(term_values.at[row, term])
        dominant_value = pd.Series(dominant_value)

        dominant_driver[no_data] = 'NO_DATA'

        # A pure effect of exactly -1 (-100 %) takes the risk of the point to zero alone: that driver is dominant
        for term in ['VULN', 'HAZ', 'EXP']:
            is_minus_one = term_values[term] == -1
            dominant_driver[is_minus_one] = term
            dominant_value[is_minus_one] = -1.0

        # One row per point: coordinates, stat, dominant term and value, ALL and the 7 terms
        df = coords.copy()
        df['dominant_driver'] = dominant_driver
        df['dominant_value'] = dominant_value
        df['ALL'] = decomp_values[0]
        df = pd.concat([df, term_values], axis=1)
        df.insert(2, 'stat', stat)

        tables.append(df)

    return pd.concat(tables, ignore_index=True)


'''
 ----------------------------------------------------------------------------------------
     IV. Computation and saving of the tables
 ----------------------------------------------------------------------------------------
'''

def main(country):
 
    country = str(country)
    out_dir = DECOMPOSITION_DIR / country
    out_dir.mkdir(parents=True, exist_ok=True)
 
    unc_outputs = load_results(country)
 
    national = country_decomposition(unc_outputs)
    national.insert(0, 'country', country)
    
    for metric in METRICS:
        name = METRIC_NAMES[metric]
        national[national['metric'] == name].to_csv(out_dir / f"{country}_decomposition_{name}.csv", index=False)
 
    points = spatial_decomposition(unc_outputs)
    points.insert(0, 'country', country)
    points.to_csv(out_dir / f"{country}_EAD_exp_point_decomposition.csv", index=False)
 
    print(f"Done: {country} -> {out_dir}")
 
#----------------------------------------------------------------------------------
 
if __name__ == "__main__":
    for country in (sys.argv[1:] or COUNTRIES):
        main(country)