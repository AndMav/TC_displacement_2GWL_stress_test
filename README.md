# README

These scripts reproduce the results of the paper with title: **Will a 2 °C adaptation plan hold? A stress test of
future tropical cyclone displacement risk by decomposing driver contributions**

*Andreas Mavrakis(1), Chahan M. Kropf(1,2,3), Simona Meiler(4), Ni Li(5,6), Yann Quilcaille(7), Lorenzo Pierini(1,2),
David N. Bresch(1,2), Lukas Riedel(1,2)*

(1) Institute for Environmental Decisions, ETH Zurich, Zurich, Switzerland

(2) Federal Office of Meteorology and Climatology MeteoSwiss, Zurich Airport, Switzerland

(3) Institute of Geography, University of Neuchatel, Neuchatel, Switzerland

(4) Department of Civil and Environmental Engineering, Stanford University, Stanford, CA, USA

(5) Department of Water and Climate, Vrije Universiteit Brussel, Brussels, Belgium

(6) Department of Hydro Sciences, TUD Dresden University of Technology, Dresden, Germany

(7) Institute for Atmospheric and Climate Science, ETH Zurich, Zurich, Switzerland

Contact: Andreas Mavrakis (andreas.mavrakis@usys.ethz.ch)


## Introduction
The repository has two folders:

* `data/`: contains the results of the vulnerability calibration, the simulation results of every country, and the
  decomposition results of every country with the world dataset that merges them (see Data). The decomposition results
  are used to create the three figures and Data S1. The simulation results do not include the expected annual
  displacement (EAD) at every exposure point of every sample: because of their size (~485 GB), these are stored at
  ETH Zurich.
* `code/`: contains `config.py` and the scripts (see Code). The first scripts create the intermediate data (exposures,
  vulnerability, hazards, simulation results and their decomposition), the last ones the figures and Data S1 of the
  paper. The scripts run with [CLIMADA](https://github.com/CLIMADA-project/climada_python) v6.1.0. Computationally
  demanding calculations of some of the scripts were run on the Euler
  cluster of ETH Zurich.


## Data
With the outputs in `data/`, the figures and Data S1 can be reproduced. These outputs are created by the vulnerability
calibration and simulation scripts, which use the proprietary hazard sets. The subfolders are the same as in
`config.py`, so that the scripts find the files:

* `vulnerability/`: calibrated (baseline, future) v_half pairs of every country (s2)
* `results/base/`, `results/delta/`: the CLIMADA uncertainty outputs of the baseline risk (base, s4a) and of the change
  in risk for every scenario (delta, s5) of every country, without the expected annual displacement at every exposure
  point (`eai_exp`) and the coordinates of the points, i.e. the sampled parameters and the AAD and RP100 of every sample.
* `results/base/all_countries_base_summaries.csv`: median AAD and RP100 of the baseline of all countries (s4b)
* `decomposition/countries/`: for every country, the decomposition of AAD and RP100 and the decomposition of EAD at
  every exposure point (s6a)
* `decomposition/world/`: the decomposition of EAD at every exposure point of all countries, merged (s7a)
* `decomposition/DataS1.csv`: Data S1 (s8)

To use them, copy the content of `data/` to `DATA_ROOT`. With these, s4b, s6b, s6c (with the exposures of s1), s7a,
s7b and s8 can be run directly.

Additionally needed, but not included:

* **WorldPop population, 2005** (Tatem 2017, https://www.worldpop.org) (s0)
* **Population SSP3, 2050** (Wang et al. 2022) (s0)
* **MIT synthetic TC event sets** of the four ocean basins (s2, s3). The hazard sets are proprietary data of
  WindRiskTech L.L.C. and cannot be publicly archived. They can be provided by WindRiskTech L.L.C. pending scientific
  review and a completed nonredistribution agreement; requests should be submitted to info@windrisktech.com.


## Code
Run the scripts in this order. Each step uses the outputs of the previous ones:
finish a step for all countries (and all scenarios in s5) before starting the next. `{ISO3}` is the ISO3 code of a
country (e.g. PHL); the countries are listed in `config.py`. s2, s4a and s5 are computationally demanding (e.g. one
cluster job per call).

#### config.py
Settings shared by all scripts: data paths, countries, GCMs and simulation scenarios. All data are read from and
written to the folder `tc_displacement` (`DATA_ROOT`), inside CLIMADA's data directory (by default `~/climada/data`).

#### s0_exposure_coarsening.py
Coarsens the WorldPop population rasters (2005 and SSP3 2050) from 30 arcsec to 300 arcsec.

`python s0_exposure_coarsening.py` (run once, for both rasters)

#### s1_exposure_creation.py
Creates the base (2005) and future (SSP3 2050) population exposure of every country from the rasters of s0, with the
same grid cells in both.

`python s1_exposure_creation.py` (run once, for all countries)

#### s2_vulnerability_calibration.py
Calibrates the future vulnerability of a country: for every baseline impact function (Kam et al. 2024), the future
v_half that gives the same average annual displacement (AAD) under the GWL2 hazard as under the baseline hazard. The
baseline impact functions of Kam et al. (2024) are included in the script.

`python s2_vulnerability_calibration.py {ISO3}` (run for each country)

#### s3_hazard_creation.py
Creates the tropical cyclone hazard of every country (baseline 1995-2014 and SSP3-7.0 2041-2060, 9 GCMs) from the MIT
event sets of the four ocean basins (loaded once for all countries; needs enough memory).

`python s3_hazard_creation.py` (run once, for all countries)

#### s4a_base_simulation.py, s4b_base_summaries.py
Uncertainty analysis of the baseline displacement risk of a country (s4a), and the median AAD and RP100 of all
countries (s4b): the base values of the change in displaced people (Fig. 2, Data S1).

`python s4a_base_simulation.py {ISO3}` (run for each country), then `python s4b_base_summaries.py` (run once, for all countries)

#### s5_simulation.py
Uncertainty analysis of the change in displacement risk of a country for one of the 7 simulation scenarios (Table S1):
`all`, `haz`, `exp`, `vuln`, `haz_exp`, `haz_vuln`, `exp_vuln`.

`python s5_simulation.py {ISO3} {scenario}` (run for each country and scenario)

#### s6a_country_decomposition.py
Decomposition of the risk change of a country into the contributions of hazard, exposure, vulnerability and their
interactions (for AAD, RP100 and the expected annual displacement (EAD) at every exposure point), and extraction of the
dominant driver.

`python s6a_country_decomposition.py` (run once, for all countries)

#### s6b_plot_barplots_grouped.py
Reproduces Fig. 2 (country decomposition of AAD and RP100 for the Philippines, Vietnam and Cuba).

`python s6b_plot_barplots_grouped.py` (run once)

#### s6c_plot_country_maps_grouped.py
Reproduces Fig. 3 (exposure change and dominant drivers at every exposure point for the Philippines, Vietnam and Cuba).

`python s6c_plot_country_maps_grouped.py` (run once)

#### s7a_world_decomposition.py, s7b_plot_world_map_grouped.py
Merges the decomposition at every exposure point of all countries (s7a) and reproduces Fig. 4 (dominant drivers at
global scale, s7b).

`python s7a_world_decomposition.py` (run once, for all countries), then `python s7b_plot_world_map_grouped.py` (run once)

#### s8_data_S1.py
Reproduces Data S1 (decomposition of the displacement risk change for all countries).

`python s8_data_S1.py` (run once, for all countries)
