"""
Summarises the baseline simulations of all countries: the median of AAD and RP100 across the samples of every
country (NaN / Inf ignored), in one table. It is the base value of the barplots and Data S1
(change in people = relative change x base value).

Input:  uncertainty_AAD_RP100.csv                         (AAD and RP100 per sample, every country)
Output: all_countries_base_summaries.csv                  (one row per country, metric and statistic)
        (folders defined in config.py)
"""

import numpy as np
import pandas as pd

from config import BASE_RESULTS_DIR, BASE_SUMMARIES_CSV, COUNTRIES

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

METRICS = ["AAD", "rp100"]                 # columns of uncertainty_AAD_RP100.csv
STATS = ["median"]                             # the base value of s6b and Data S1


'''
 ----------------------------------------------------------------------------------------
     II. Statistics of one country
 ----------------------------------------------------------------------------------------
'''

def summarize(values, stat):
    """Median or percentile ('pXX') of an array, ignoring NaN."""
    if stat == "median":
        return float(np.nanmedian(values))
    if stat.startswith("p"):
        return float(np.nanpercentile(values, float(stat[1:])))
    raise ValueError(f"Unrecognized stat: {stat!r}")


def country_summaries(country):
    """Statistics (STATS) of AAD and RP100 of one country (long format), or None if its output is missing."""
    path = BASE_RESULTS_DIR / country / "uncertainty_AAD_RP100.csv"
    if not path.exists():
        return None

    samples = pd.read_csv(path).replace([np.inf, -np.inf], np.nan)          # Inf treated as missing
    rows = []
    for stat in STATS:
        for metric in METRICS:
            value = summarize(samples[metric].values, stat) if metric in samples.columns else np.nan
            rows.append({"country": country, "metric": metric, "stat": stat, "value": value})
    return pd.DataFrame(rows, columns=["country", "metric", "stat", "value"])


'''
 ----------------------------------------------------------------------------------------
     III. All countries and saving
 ----------------------------------------------------------------------------------------
'''

def main(countries=COUNTRIES):
    tables, missing = [], []
    for country in countries:
        table = country_summaries(country)
        if table is None:
            missing.append(country)
            continue
        tables.append(table)

    BASE_SUMMARIES_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(BASE_SUMMARIES_CSV, index=False)

    print(f"Done: {len(tables)}/{len(countries)} countries -> {BASE_SUMMARIES_CSV}")
    if missing:
        print(f"No base output: {', '.join(missing)}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()