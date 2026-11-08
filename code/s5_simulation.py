"""
Uncertainty analysis of the change in TC displacement risk of a country between the baseline and the future, for one
simulation scenario (sim_type): only the drivers set to 'fut' change, the others stay in their baseline state.

    sim_type   hazard  exposure  vulnerability
    all        fut     fut       fut
    haz        fut     base      base
    exp        base    fut       base
    vuln       base    base      fut
    haz_exp    fut     fut       base
    haz_vuln   fut     base      fut
    exp_vuln   base    fut       fut

Uncertain inputs:
    gc_model   GCM (9 CMIP6 models), shared by the baseline and future hazard
    HE_base    seed of the baseline hazard event resampling (80 %, with replacement), for internal variability
    HE_fut     seed of the future hazard event resampling (80 %, with replacement), for internal variability
    idx        row of the country's calibrated (baseline, future) v_half pairs, shared by the baseline and future
               impact function

Input:  TC_{country}_{gcm}_MIT_20thcal.hdf5, TC_{country}_{gcm}_MIT_ssp370cal.hdf5   (baseline and future hazard)
        {country}_base_2005_300as.h5, {country}_SSP3_2050_300as.h5                  (base and future exposure)
        {country}.csv                                                               (calibrated v_half pairs)
Output: unsequa_delta_TC_{country}_{sim_type}_MIT.hdf5                              (CLIMADA uncertainty output)
        (folders defined in config.py)
"""

import sys

import numpy as np
import pandas as pd
import scipy as sp

from climada.hazard import Hazard
from climada.entity import Exposures
from climada.entity.impact_funcs import ImpactFuncSet, ImpfTropCyclone
from climada.engine.unsequa import InputVar, CalcDeltaImpact
from climada.util.coordinates import km_to_degree

from config import HAZ_DIR, EXP_DIR, V_HALF_DIR, DELTA_RESULTS_DIR, GCMS

'''
 ----------------------------------------------------------------------------------------
     I. Constants and paths 
 ---------------------------------------------------------------------------------------- 
'''



# GCMs of the MIT TC hazard sets (the integer is the sampled value of gc_model)
h1_min, h1_max = (1, 9)
gcm_model_key = dict(enumerate(GCMS, start=1))   # {1: 'cesm2', 2: 'cnrm6', ..., 9: 'ukmo6'}

# Configuration of simulation scenarios, based on which driver is set to its future state in each scenario
SCENARIO_CONFIG = {
    'all':      {'haz': 'fut',  'exp': 'fut',  'impf': 'fut'},
    'haz':       {'haz': 'fut',  'exp': 'base', 'impf': 'base'},
    'exp':      {'haz': 'base', 'exp': 'fut',  'impf': 'base'},
    'vuln':     {'haz': 'base', 'exp': 'base', 'impf': 'fut'},
    'haz_exp':   {'haz': 'fut',  'exp': 'fut',  'impf': 'base'},
    'haz_vuln':  {'haz': 'fut',  'exp': 'base', 'impf': 'fut'},
    'exp_vuln': {'haz': 'base', 'exp': 'fut',  'impf': 'fut'},
}

N_POWER = 10               # base sample size N = 2**N_POWER
SAMPLE_PERC = 0.8          # share of events drawn (with replacement) per hazard sample
CENTROID_MAX_DIST_KM = 20  # max. distance exposure point -> hazard centroid


'''
 ----------------------------------------------------------------------------------------
     II. Helper Functions 
 ----------------------------------------------------------------------------------------
'''

def load_hazards(country, period):
    """
    Loads the TC hazard of every GCM for one period.
    period: '20thcal' (baseline, 1995-2014) or 'ssp370cal' (future, SSP3-7.0)
    Returns: Dictionary of hazard objects, one for each GCM for the selected period. 
    """
    
    haz_dict = {}
    for h1 in range(h1_min, h1_max + 1):
        haz = Hazard.from_hdf5(HAZ_DIR / country / f"TC_{country}_{gcm_model_key[h1]}_MIT_{period}.hdf5")
        haz.check()
        haz_dict[gcm_model_key[h1]] = haz
        
    return haz_dict

def sample_hazard(haz, seed):
    """Internal variability: bootstrap resampling of a hazard's event set.
        draws randomly the SAMPLE_PERC = 80% of the events with replacement, reproducibly for a given seed.
        Returns hazard object with the selected events"""
    
    rng = np.random.default_rng(seed=int(seed))
    haz_sampled = haz.select(event_id=rng.choice(haz.event_id, size=int(np.round(haz.size * SAMPLE_PERC)), replace=True).tolist())
    
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
    
def load_impf_sets(country):
    """
    Reads the country's calibrated (baseline, future) v_half pairs, builds impact function sets-objects, one per v_half.
    Each set contains a single impact function-object
    Returns: one list for each (baseline, future), where element i of both lists corresponds to row i of the v_half pair in the CSV
    """
    
    v_half_df = pd.read_csv(V_HALF_DIR / f"{country}.csv", index_col=0)
    impf_set_base = [impf_func(v_half=v) for v in v_half_df["baseline"]]
    impf_set_fut = [impf_func(v_half=v) for v in v_half_df["future"]]
    
    return impf_set_base, impf_set_fut

'''
 -----------------------------------------------------------------------------------------
     III. CalcDeltaImpact configuration for the selected simulation type 
 -----------------------------------------------------------------------------------------
'''

def build_calc_delta_imp(country, sim_type):
    """
    One function that:
    1) Loads all inputs for the selected country
    2) Returns the CalcDeltaImpact object for the selected simulation type
    """

    #----Hazard-----------------------------------------------------
    
    haz_base = load_hazards(country, "20thcal")
    haz_fut = load_hazards(country, "ssp370cal")

    # Resampling per GCM hazard set
    def haz_base_func(gc_model, HE_base):
        return sample_hazard(haz_base[gcm_model_key[gc_model]], HE_base)
 
    def haz_fut_func(gc_model, HE_fut):
        return sample_hazard(haz_fut[gcm_model_key[gc_model]], HE_fut)

    # Distributions of the uncertain parameters
    haz_base_distr = {
        "gc_model": sp.stats.randint(low=h1_min, high=h1_max + 1),  # GCM 1..9
        "HE_base": sp.stats.randint(0, 2**32 - 1),                  # resampling seed
    }
    
    haz_fut_distr = {
        "gc_model": sp.stats.randint(low=h1_min, high=h1_max + 1),  # same name: same GCM
        "HE_fut": sp.stats.randint(0, 2**32 - 1),                   # own seed: independent
    }
    
    haz_base_iv = InputVar(haz_base_func, haz_base_distr)
    haz_fut_iv = InputVar(haz_fut_func, haz_fut_distr)
    
    #----Exposure-----------------------------------------------------
    
    # Exposures (no uncertain inputs); all GCM hazard sets share the same centroids:
    # assigning to one of them (ukmo6) is enough
    exp_base = load_exposure(country, f"{country}_base_2005_300as.h5", haz_base["ukmo6"])
    exp_fut = load_exposure(country, f"{country}_SSP3_2050_300as.h5", haz_fut["ukmo6"])

    # ---- Vulnerability ---------------------------------------------------

    impf_set_base, impf_set_fut = load_impf_sets(country)

    # Impact function set of row idx
    def impf_func_base(idx):
        return impf_set_base[int(idx)]

    def impf_func_fut(idx):
        return impf_set_fut[int(idx)]

    # Distributions of the uncertain parameters
    n_impf = len(impf_set_base)
    impf_distr_base = {'idx': sp.stats.randint(low=0, high=n_impf)}  # row 0..n_impf-1
    impf_distr_fut = {'idx': sp.stats.randint(low=0, high=n_impf)}   # same name: same row

    vuln_base_iv = InputVar(impf_func_base, impf_distr_base)
    vuln_fut_iv = InputVar(impf_func_fut, impf_distr_fut)

    # ---- CalcDeltaImpact -------------------------------------------------

    config = SCENARIO_CONFIG[sim_type]

    # Future components of the scenario
    haz_final_iv = haz_fut_iv if config['haz'] == 'fut' else haz_base_iv
    exp_final = exp_fut if config['exp'] == 'fut' else exp_base
    vuln_final_iv = vuln_fut_iv if config['impf'] == 'fut' else vuln_base_iv

    return CalcDeltaImpact(exp_base, vuln_base_iv, haz_base_iv,
                           exp_final, vuln_final_iv, haz_final_iv)

'''
 -----------------------------------------------------------------------------------------
     IV. Uncertainty Analysis Computation and Object save
 -----------------------------------------------------------------------------------------
'''

def main(country, sim_type):

    country = str(country)
    sim_type = str(sim_type)
    N_samples = 2**N_POWER
    
    if sim_type not in SCENARIO_CONFIG:
        raise ValueError(f"sim_type must be one of {list(SCENARIO_CONFIG)}, got '{sim_type}'")

    out_dir = DELTA_RESULTS_DIR / country / sim_type
    out_dir.mkdir(parents=True, exist_ok=True)

    calc_delta_imp = build_calc_delta_imp(country, sim_type)

    # Samples
    output_imp = calc_delta_imp.make_sample(N=N_samples)

    # Uncertainty analysis
    output_imp = calc_delta_imp.uncertainty(output_imp, rp=[100], calc_eai_exp=True, processes=4, chunksize=None)
    output_imp.to_hdf5(out_dir / f"unsequa_delta_TC_{country}_{sim_type}_MIT.hdf5")

    print(f"Done: {country} / {sim_type} -> {out_dir}")
    
#----------------------------------------------------------------------------------
        
if __name__ == "__main__":
    main(*sys.argv[1:])