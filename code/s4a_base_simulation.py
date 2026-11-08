"""
Uncertainty analysis of the baseline TC displacement risk of a country: the people displaced with the baseline hazard,
the 2005 exposure and the baseline vulnerability. Gives the base values of the barplots
(change in people = relative change x base value).

Uncertain inputs:
    gc_model   GCM (9 CMIP6 models)
    HE_base    seed of the hazard event resampling (80 %, with replacement), for internal variability
    idx        row of the country's calibrated v_half values (baseline column)

Input:  TC_{country}_{gcm}_MIT_20thcal.hdf5               (baseline hazard)
        {country}_base_2005_300as.h5                      (2005 exposure)
        {country}.csv                                     (calibrated v_half pairs)
Output: unsequa_base_TC_{country}_MIT.hdf5                (CLIMADA uncertainty output)
        uncertainty_AAD_RP100.csv                         (AAD and RP100, one row per sample)
        (folders defined in config.py)
"""

import sys

import numpy as np
import pandas as pd
import scipy as sp

from climada.hazard import Hazard
from climada.entity import Exposures
from climada.entity.impact_funcs import ImpactFuncSet, ImpfTropCyclone
from climada.engine.unsequa import InputVar, CalcImpact
from climada.util.coordinates import km_to_degree

from config import HAZ_DIR, EXP_DIR, V_HALF_DIR, BASE_RESULTS_DIR, GCMS

'''
 ----------------------------------------------------------------------------------------
     I. Constants and paths
 ----------------------------------------------------------------------------------------
'''

# GCMs of the MIT TC hazard sets (the integer is the sampled value of gc_model)
h1_min, h1_max = (1, 9)
gcm_model_key = dict(enumerate(GCMS, start=1))   # {1: 'cesm2', 2: 'cnrm6', ..., 9: 'ukmo6'}

N_POWER = 10               # base sample size N = 2**N_POWER (Sobol: 8 N samples with 3 inputs)
SAMPLE_PERC = 0.8          # share of events drawn (with replacement) per hazard sample
RESCALE_FREQUENCY = True   # divide the frequencies of a sample by SAMPLE_PERC: the storm rate of the full set
CENTROID_MAX_DIST_KM = 20  # max. distance exposure point -> hazard centroid

METRICS = ["aai_agg", "rp100"]   # per-sample metrics written to uncertainty_AAD_RP100.csv (aai_agg as AAD)


'''
 ----------------------------------------------------------------------------------------
     II. Helper Functions
 ----------------------------------------------------------------------------------------
'''

def load_hazards(country, period="20thcal"):
    """
    Loads the TC hazard of every GCM for one period (here the baseline, '20thcal', 1995-2014).
    Returns: Dictionary of hazard objects, one for each GCM.
    """
    haz_dict = {}
    for h1 in range(h1_min, h1_max + 1):
        haz = Hazard.from_hdf5(HAZ_DIR / country / f"TC_{country}_{gcm_model_key[h1]}_MIT_{period}.hdf5")
        haz.check()
        haz_dict[gcm_model_key[h1]] = haz

    return haz_dict

def sample_hazard(haz, seed):
    """Internal variability: bootstrap resampling of a hazard's event set.
        Draws randomly SAMPLE_PERC = 80 % of the events with replacement, reproducibly for a given seed; the repeated
        events get unique IDs. With RESCALE_FREQUENCY, the frequencies are divided by SAMPLE_PERC (storm rate of the
        full set). Returns hazard object with the selected events"""

    rng = np.random.default_rng(seed=int(seed))
    haz_sampled = haz.select(event_id=rng.choice(haz.event_id, size=int(np.round(haz.size * SAMPLE_PERC)), replace=True).tolist())
    haz_sampled.sanitize_event_ids()                       # duplicates from the bootstrap -> unique IDs

    if RESCALE_FREQUENCY:
        haz_sampled.frequency = haz_sampled.frequency / SAMPLE_PERC

    return haz_sampled

def load_exposure(country, file_name, haz):
    """
    Loads an exposure and assigns each point to the nearest hazard centroid. All points use impact function id = 1.
    """
    exp = Exposures.from_hdf5(EXP_DIR / country / file_name)
    exp.assign_centroids(haz, threshold=km_to_degree(CENTROID_MAX_DIST_KM))
    exp.gdf["impf_TC"] = 1

    return exp

def impf_func(v_half: float = 45.6, scale: float = 1.0):
    """Returns an impact function set, that contains a single Emanuel-type TC impact function for a given v_half value"""

    imp_fun = ImpfTropCyclone.from_emanuel_usa(intensity=np.arange(0, 150, 1), v_half=v_half, scale=scale)
    imp_fun.haz_type = "TC"
    imp_fun.intensity_unit = "m/s"
    imp_fun.check()

    return ImpactFuncSet([imp_fun])

def load_impf_set_base(country):
    """
    Reads the country's calibrated v_half values (baseline column) and builds one impact function set per value.
    """
    v_half_df = pd.read_csv(V_HALF_DIR / f"{country}.csv", index_col=0)
    return [impf_func(v_half=v) for v in v_half_df["baseline"]]


'''
 -----------------------------------------------------------------------------------------
     III. CalcImpact configuration
 -----------------------------------------------------------------------------------------
'''

def build_calc_imp(country):
    """
    Loads all inputs of the country and returns the CalcImpact object of the baseline.
    """

    #----Hazard-----------------------------------------------------

    haz_base = load_hazards(country, "20thcal")

    # Resampling per GCM hazard set
    def haz_base_func(gc_model, HE_base):
        return sample_hazard(haz_base[gcm_model_key[int(gc_model)]], HE_base)

    haz_base_distr = {
        "gc_model": sp.stats.randint(low=h1_min, high=h1_max + 1),  # GCM 1..9
        "HE_base": sp.stats.randint(0, 2**32 - 1),                  # resampling seed
    }
    haz_base_iv = InputVar(haz_base_func, haz_base_distr)

    #----Exposure-----------------------------------------------------

    # all GCM hazard sets share the same centroids: assigning to one of them (ukmo6) is enough
    exp_base = load_exposure(country, f"{country}_base_2005_300as.h5", haz_base["ukmo6"])

    # ---- Vulnerability ---------------------------------------------------

    impf_set_base = load_impf_set_base(country)

    def impf_func_base(idx):
        return impf_set_base[int(idx)]

    impf_distr_base = {'idx': sp.stats.randint(low=0, high=len(impf_set_base))}   # row 0..n_impf-1
    impf_base_iv = InputVar(impf_func_base, impf_distr_base)

    # ---- CalcImpact (order: exposure, impact functions, hazard) ------------

    return CalcImpact(exp_base, impf_base_iv, haz_base_iv)


'''
 -----------------------------------------------------------------------------------------
     IV. Uncertainty Analysis Computation and saving
 -----------------------------------------------------------------------------------------
'''

def nan_diagnostic(unc_df, country):
    """Number of NaN/Inf samples of aai_agg and rp100, as text."""
    lines = [f"NaN/Inf diagnostic for {country} / base", "=" * 50, ""]
    for metric in METRICS:
        if metric not in unc_df.columns:
            lines.append(f"{metric}: column not found")
            continue
        n_bad = int((~np.isfinite(unc_df[metric].astype(float))).sum())
        lines.append(f"{metric}: {n_bad}/{len(unc_df)} samples are NaN/Inf ({100 * n_bad / len(unc_df):.1f}%)")
    return "\n".join(lines)


def main(country):

    country = str(country)
    N_samples = 2**N_POWER

    out_dir = BASE_RESULTS_DIR / country
    out_dir.mkdir(parents=True, exist_ok=True)

    calc_imp = build_calc_imp(country)

    # Samples
    output_imp = calc_imp.make_sample(N=N_samples)

    # Uncertainty analysis
    output_imp = calc_imp.uncertainty(output_imp, rp=[100], calc_eai_exp=False, processes=4, chunksize=None)
    output_imp.to_hdf5(out_dir / f"unsequa_base_TC_{country}_MIT.hdf5")

    # AAD (aai_agg) and RP100 per sample (read by s4b_base_summaries.py)
    unc_df = output_imp.get_uncertainty(metric_list=["aai_agg", "freq_curve"])
    unc_df.rename(columns={"aai_agg": "AAD"}).to_csv(out_dir / "uncertainty_AAD_RP100.csv", index=False)

    # NaN/Inf check (printed)
    print(nan_diagnostic(unc_df, country))

    print(f"Done: {country} / base -> {out_dir}")

#----------------------------------------------------------------------------------

if __name__ == "__main__":
    main(*sys.argv[1:])