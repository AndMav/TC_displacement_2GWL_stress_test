"""
Creates Data S1: for every country and metric (AAD, RP100), the total change in displaced people at the median, the
dominant driver and its sign at p10, median and p90, and the relative change (%) of the total (ALL) and of every
decomposition term. The dominant driver is the term with the largest absolute value; a pure effect (HAZ, EXP, VULN) of
exactly -100 % dominates (EXP before HAZ before VULN). "-": no value.

Input:  {country}_decomposition_AAD.csv, {country}_decomposition_rp100.csv   (country decompositions)
        all_countries_base_summaries.csv                                   (base values)
Output: DataS1.csv
        (folders defined in config.py)
"""

import numpy as np
import pandas as pd
import pycountry

from config import DECOMPOSITION_DIR, BASE_SUMMARIES_CSV, DATA_S1_CSV, COUNTRIES

'''
 ----------------------------------------------------------------------------------------
     I. Setup
 ----------------------------------------------------------------------------------------
'''

METRICS = {"AAD": "AAD", "rp100": "RP100"}          # metric in the input files -> label in Data S1
STATS = ["p10", "median", "p90"]
TERMS = ["HAZ", "EXP", "VULN", "HAZxEXPint", "HAZxVULNint", "EXPxVULNint", "HAZxEXPxVULNint"]
VALUE_TERMS = ["ALL"] + TERMS
NAMES = {"HAZ": "hazard", "EXP": "exposure", "VULN": "vulnerability",
         "HAZxEXPint": "hazard-exposure interaction", "HAZxVULNint": "hazard-vulnerability interaction",
         "EXPxVULNint": "exposure-vulnerability interaction", "HAZxEXPxVULNint": "three-way interaction"}
NA = "-"

# Vulnerability regions (Table S2)
REGIONS = {
    "North Atlantic 1": "ABW ATG BHS BLZ BMU BRB COL CPV CRI CUB CYM DMA DOM GRD GTM GUY HND HTI JAM KNA LCA MEX MSR NIC PAN PRI SLV SUR TCA TTO VCT VEN VGB VIR",
    "North Atlantic 2": "CAN USA",
    "North Indian": "AFG ARE BGD BHR BTN DJI ERI ETH IND IRN IRQ LKA MDV MMR NPL OMN PAK QAT SAU SOM YEM",
    "Oceania": "ASM AUS COK FJI NCL NIU NZL PNG PYF SLB TLS TON VUT WLF WSM",
    "South Indian": "COM MDG MOZ MUS MWI SWZ TZA ZAF ZWE",
    "North Western Pacific 1": "IDN KHM LAO MYS THA VNM",
    "North Western Pacific 2": "PHL",
    "North Western Pacific 3": "CHN",
    "North Western Pacific 4": "HKG JPN KOR TWN",
}
REGION = {iso: reg for reg, isos in REGIONS.items() for iso in isos.split()}


'''
 ----------------------------------------------------------------------------------------
     II. Helper functions
 ----------------------------------------------------------------------------------------
'''

def country_name(iso):
    c = pycountry.countries.get(alpha_3=iso)
    return getattr(c, "common_name", c.name) if c else iso


def load_values(iso, metric):
    """Decomposition of one country and metric (ALL and the 7 terms); None if the file is missing."""
    path = DECOMPOSITION_DIR / iso / f"{iso}_decomposition_{metric}.csv"
    return pd.read_csv(path) if path.exists() else None


def stat_values(table, stat, terms):
    """Values of the given terms for one statistic, without NaN."""
    if table is None:
        return pd.Series(dtype=float)
    return table[(table["stat"] == stat) & (table["decomp_term"].isin(terms))].set_index("decomp_term")["delta_change"].dropna()


def dominant(values):
    """Term with the largest absolute value; a pure effect of exactly -1 (-100 %) dominates, EXP before HAZ before VULN."""
    for term in ["EXP", "HAZ", "VULN"]:
        if values.get(term) == -1:
            return term
    return values.abs().idxmax()


'''
 ----------------------------------------------------------------------------------------
     III. Data S1
 ----------------------------------------------------------------------------------------
'''

def main():
    base = pd.read_csv(BASE_SUMMARIES_CSV)
    base = base[base["stat"] == "median"]

    rows = []
    for iso in sorted(COUNTRIES):
        for metric, label in METRICS.items():
            row = {"ISO3": iso, "country": country_name(iso), "vulnerability_region": REGION[iso], "metric": label}
            b = base[base["metric"] == metric].set_index("country")["value"].get(iso, np.nan)
            table = load_values(iso, metric)

            for stat in STATS:
                # dominant driver and its sign
                t = stat_values(table, stat, TERMS)
                if t.empty:
                    row[f"dominant_{stat}"], row[f"sign_{stat}"] = NA, NA
                else:
                    d = dominant(t)
                    row[f"dominant_{stat}"] = NAMES[d]
                    row[f"sign_{stat}"] = "increase" if t[d] >= 0 else "decrease"

                # relative change (%) of ALL and of every term
                v = stat_values(table, stat, VALUE_TERMS)
                for term in VALUE_TERMS:
                    row[f"{stat}_{term}"] = round(100 * float(v[term]), 1) if term in v.index else NA

                # change in displaced people at the median: median relative change of ALL x base value (as in Fig. 2)
                if stat == "median":
                    ok = "ALL" in v.index and not pd.isna(b) and b != 0
                    row["change_displaced_people_median"] = round(float(v["ALL"]) * float(b), 1) + 0.0 if ok else NA
            rows.append(row)

    cols = (["ISO3", "country", "vulnerability_region", "metric", "change_displaced_people_median"]
            + [f"{k}_{s}" for s in STATS for k in ("dominant", "sign")]
            + [f"{s}_{term}" for s in STATS for term in VALUE_TERMS])
    DATA_S1_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=cols).to_csv(DATA_S1_CSV, index=False)

    print(f"Done: {len(rows)} rows -> {DATA_S1_CSV}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
