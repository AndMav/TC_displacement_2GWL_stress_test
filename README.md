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


## Content

#### config.py
Settings shared by all scripts: data paths, countries, GCMs and simulation scenarios. All data are read from and
written below `DATA_ROOT`, a folder in CLIMADA's data directory (`~/climada/data/tc_displacement` by default). To use
another location, change the CLIMADA data directory in `climada.conf`.

#### s0_exposure_coarsening.py
Coarsens the WorldPop population rasters (2005 and SSP3 2050) from 30 arcsec to 300 arcsec.

#### s1_exposure_creation.py
Creates the base (2005) and future (SSP3 2050) population exposure of every country, with the same grid cells in both.

#### s2_vulnerability_calibration.py
Calibrates the future vulnerability of a country: for every baseline impact function (Kam et al. 2024), the future
v_half that gives the same average annual displacement (AAD) under the GWL2 hazard as under the baseline hazard.

#### s3_hazard_creation.py
Creates the tropical cyclone hazard of every country (baseline 1995-2014 and SSP3-7.0 2041-2060, 9 GCMs) from the MIT
event sets of the four ocean basins.

#### s4a_base_simulation.py, s4b_base_summaries.py
Uncertainty analysis of the baseline displacement risk of a country (s4a), and the median AAD and RP100 of all
countries (s4b): the base values of the change in displaced people (Fig. 2, Data S1).

#### s5_simulation.py
Uncertainty analysis of the change in displacement risk of a country for one of the 7 simulation scenarios (Table S1).

#### s6a_country_decomposition.py
Decomposition of the risk change of a country into the contributions of hazard, exposure, vulnerability and their
interactions, for AAD and RP100 and for the expected annual displacement (EAD) at every exposure point.

#### s6b_plot_barplots_grouped.py
Reproduces Fig. 2 (country decomposition of AAD and RP100 for the Philippines, Vietnam and Cuba).

#### s6c_plot_country_maps_grouped.py
Reproduces Fig. 3 (exposure change and dominant drivers at every exposure point for the Philippines, Vietnam and Cuba).

#### s7a_world_decomposition.py, s7b_plot_world_map_grouped.py
Merges the decomposition at every exposure point of all countries (s7a) and reproduces Fig. 4 (dominant drivers at
global scale, s7b).

#### s8_data_S1.py
Reproduces Data S1 (decomposition of the displacement risk change for all countries).


## Input data
The following inputs are not created by the scripts. Place them in the folders defined in `config.py` (below
`DATA_ROOT`), with the file names given there:

* **WorldPop population, 2005** (Tatem 2017, https://www.worldpop.org): `exposure/worldpop_30as/ppp_2005_1km_Aggregated.tif`
* **Population SSP3, 2050** (Wang et al. 2022): `exposure/worldpop_30as/SSP3_2050_30as.tif`
* **MIT synthetic TC event sets** of the four ocean basins, `TC_{basin}_0300as_MIT_{gcm}_{ssp}{cal}.hdf5`:
  `hazard/basin_event_sets/`. The hazard sets are proprietary data of WindRiskTech L.L.C. and cannot be publicly
  archived. They can be provided by WindRiskTech L.L.C. pending scientific review and a completed nonredistribution
  agreement; requests should be submitted to info@windrisktech.com. Needed, each for the four basins (AP, IO, SH, WP):
  * s3 (all 9 GCMs: cesm2, cnrm6, ecearth6, fgoals, ipsl6, miroc6, mpi6, mri6, ukmo6): baseline `{gcm}_20thcal` and
    SSP3-7.0 `{gcm}_ssp370cal`
  * s2 (GWL2 models): baseline `{gcm}_20thcal` of cnrm6, mri6, fgoals, mpi6, miroc6, and future `cnrm6_ssp245cal`,
    `mri6_ssp245cal`, `fgoals_ssp370cal`, `mpi6_ssp370cal`, `miroc6_ssp585cal`, `mpi6_ssp585cal`

The calibrated displacement impact functions of Kam et al. (2024) are included in `s2_vulnerability_calibration.py`.


## How to run
Run the scripts in this order, from the folder of the code. Each step uses the outputs of the previous ones: finish a
step for all countries (and all scenarios in s5) before starting the next. `{ISO3}` is the ISO3 code of a country
(e.g. PHL); the countries are listed in `config.py`.

| Step | Script | Command | Note |
|---|---|---|---|
| 0 | `s0_exposure_coarsening.py` | `python s0_exposure_coarsening.py` | Both rasters (2005, 2050) in one run |
| 1 | `s1_exposure_creation.py` | `python s1_exposure_creation.py` | All countries in one run |
| 2 | `s2_vulnerability_calibration.py` | `python s2_vulnerability_calibration.py {ISO3}` | One country per call (computationally demanding: e.g. one cluster job per country) |
| 3 | `s3_hazard_creation.py` | `python s3_hazard_creation.py` | All countries in one run (the basin event sets are loaded once for all countries; needs enough memory) |
| 4a | `s4a_base_simulation.py` | `python s4a_base_simulation.py {ISO3}` | One country per call (computationally demanding: e.g. one cluster job per country) |
| 4b | `s4b_base_summaries.py` | `python s4b_base_summaries.py` | All countries in one run |
| 5 | `s5_simulation.py` | `python s5_simulation.py {ISO3} {scenario}` | One country and one scenario per call (computationally demanding: e.g. one cluster job per country and scenario); `{scenario}`: `all`, `haz`, `exp`, `vuln`, `haz_exp`, `haz_vuln`, `exp_vuln` |
| 6a | `s6a_country_decomposition.py` | `python s6a_country_decomposition.py` | All countries in one run |
| 6b | `s6b_plot_barplots_grouped.py` | `python s6b_plot_barplots_grouped.py` | Fig. 2 |
| 6c | `s6c_plot_country_maps_grouped.py` | `python s6c_plot_country_maps_grouped.py` | Fig. 3 |
| 7a | `s7a_world_decomposition.py` | `python s7a_world_decomposition.py` | All countries in one run |
| 7b | `s7b_plot_world_map_grouped.py` | `python s7b_plot_world_map_grouped.py` | Fig. 4 |
| 8 | `s8_data_S1.py` | `python s8_data_S1.py` | Data S1 |


## Data
The folder `data/` contains the outputs needed to reproduce the figures and Data S1 without the (proprietary) hazard
sets, in the folder structure of `config.py`:

* `vulnerability/`: calibrated (baseline, future) v_half pairs of every country (s2)
* `results/base/`, `results/delta/`: the CLIMADA uncertainty outputs of every country (s4a) and scenario (s5), without
  the expected annual displacement at every exposure point (`eai_exp`) and the coordinates of the points, i.e. the
  sampled parameters and the AAD and RP100 of every sample. They load with
  `UncImpactOutput.from_hdf5` / `UncDeltaImpactOutput.from_hdf5`. The full outputs (~485 GB, mostly `eai_exp`) exceed
  repository limits and are stored at ETH Zurich.
* `results/base/all_countries_base_summaries.csv`: median AAD and RP100 of the baseline of all countries (s4b)
* `decomposition/`: country decomposition of AAD and RP100 and decomposition of EAD at every exposure point (s6a, s7a),
  and Data S1 (s8)

To use them, copy the content of `data/` to `DATA_ROOT`. With these, s4b, s6b, s6c (with the exposures of s1), s7a,
s7b and s8 can be run directly.


## Requirements
* Python environment with _CLIMADA_ v6.1.0 (best installed with conda):
        https://github.com/CLIMADA-project/climada_python
* Additional packages: `rioxarray` (s0), `pycountry` (s8)


## ETH cluster
Computationally demanding calculations (s2, s4a, s5) were run on the [Euler cluster of ETH Zurich](https://scicomp.ethz.ch/wiki/Euler).


## License
GNU General Public License v3.0 (see `LICENSE`), as CLIMADA.


## Citation
If you use this code, please cite the paper and the archived code version (Zenodo DOI, see `CITATION.cff`).
