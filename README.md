# Testing the cessation hypothesis in ultralow-viscosity basanites, La Palma

Joshua Dennett, MEarthSci Earth and Planetary Science, University of Manchester, 2026.

This study tests whether the Standard Rheologic Approach (SRA) captures the dynamic viscosity of a lava flow during emplacement, or the terminal cessation viscosity reached when the flow locks viscously, for the 1949 San Juan and 1971 Teneguía flows on La Palma. Field morphometry was propagated through Monte Carlo error analysis (10⁴ draws per station) and compared with PyFLOWGO simulations. SRA bulk estimates (10⁴ to 10⁶ Pa s) exceed the PyFLOWGO core predictions (10² Pa s) by about three orders of magnitude. For Teneguía, bootstrap extrapolation of the SRA spatial trends predicts the PyFLOWGO cessation distance within the 95% confidence interval. San Juan fails to reconstruct distal cessation where spatial coverage is below 80% of the flow length. SRA therefore records cessation rheology rather than the dynamic viscosity of the flowing core.

![La Palma study area with the 1949 San Juan and 1971 Teneguía flow paths](Results%20Figures/map_lapalma_overview_pub.png)

![SRA viscosity trend extrapolation against PyFLOWGO cessation, 1971 Teneguía and 1949 San Juan](Results%20Figures/bootstrap_viscosity_publication_20260304_140415.png)

## Repository layout

| Path | Contents |
| --- | --- |
| `pyflowgo/` | Vendored PyFLOWGO source (unmodified upstream files) plus the calibration, run and Monte Carlo scripts written for this study. `resource/LaPalma*.json` and `resource/LaPalma*_slope_profile.txt` hold the flow configurations and slope profiles; `results_flowgo/` holds the simulation outputs. |
| `sra/` | SRA Monte Carlo and bootstrap scripts, the field data (`sheet2_out.csv`, `LaPalmaData_Expanded_DATA_COMPREHENSION.xlsx`), the Monte Carlo results, plotting scripts and intermediate figures. |
| `Results Figures/` | Publication maps and figures (PNG and PDF). |
| `make_pub_maps.py` | Publication hillshade maps of both flows and the island overview. Needs the DEM tiles. |
| `build_mapbox_map.py` | Interactive Mapbox methodology map. Needs the DEM tiles and a Mapbox token. |
| `Excel_tabel.csv`, `Extract_*.xlsx` | GIS table exports of slope against distance from the vent. |

## Setup

The code was developed and run with Python 3.14 on Windows, using numpy 2.4, scipy 1.17, matplotlib 3.10, pandas 3.0, rasterio 1.5 and pyproj 3.7.

```sh
python -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Running the analysis

The steps below follow the order of the study. Several scripts use paths relative to the working directory, so run the `pyflowgo/` scripts from inside `pyflowgo/` and the `sra/` scripts from inside `sra/`. The scripts that locate the repository from their own file path (`make_pub_maps.py`, `build_mapbox_map.py`, `monte_carlo_pyflowgo.py`, `run_multi_resolution.py`, `replot_resolution_comparison.py`, and in `sra/` the bootstrap, overlay, `plot_ys_comparison.py` and resolution-comparison scripts) run from anywhere.

All configurations, calibrated parameters, simulation outputs and Monte Carlo results are already committed, so you can start at any step.

1. Calibrate PyFLOWGO to the observed flow lengths (from `pyflowgo/`).
   ```sh
   cd pyflowgo
   python calibrate_1949.py        # San Juan, target flow length 1200 m
   python calibrate_lapalma.py     # Teneguía, target flow length 2000 m
   python find_er_6000m_1949.py    # effusion-rate sweep to ~6000 m (San Juan)
   python find_q_1971_2850m.py     # effusion-rate search to ~2850 m (Teneguía)
   ```
   The calibration scripts write their best-fit parameters back into the JSON files in `resource/` (`calibrate_lapalma.py` writes `LaPalma1971_calibrated.json`; `calibrate_1949.py` overwrites `LaPalma1949.json`), so re-running them changes the committed configurations. Scratch files are git-ignored. `extend_1971.py` appends to the 1971 slope profile in place, so it has already been applied and should not be run again. `calibrate_2021_*.py` and `plot_viscosity_yield_sept29.py` cover the 2021 eruption configurations.

2. Run PyFLOWGO with the calibrated configurations.
   ```sh
   python run_1949_6000m.py
   python run_1971_2m.py
   python run_1971_1m.py
   python run_multi_resolution.py       # 2, 10 and 30 m resolution comparison
   python monte_carlo_pyflowgo.py       # PyFLOWGO parameter Monte Carlo
   ```
   Results are written to `pyflowgo/results_flowgo/`. A single run can also be made with `python main_flowgo.py resource/LaPalma1971_calibrated.json`.

3. Run the SRA Monte Carlo (from `sra/`).
   ```sh
   cd ../sra
   python monte_carlo_sra.py            # sheet2_out.csv -> LaPalmaData_MC_Results.csv
   ```
   `find_k.py` and `test_k.py` explore the downstream loss factor k.

4. Run the bootstrap cessation analysis.
   ```sh
   python bootstrap_cessation.py
   ```
   Output figures carry a timestamp in their file names.

5. Make the figures.
   ```sh
   python plot_results.py
   python plot_comparison.py
   python plot_combined.py
   python plot_uncorrected_vs_pyflowgo.py
   python plot_ys_comparison.py
   python generate_pyflowgo_overlays.py
   python regen_resolution_comparison_1949.py
   python regen_resolution_comparison_1971.py
   ```
   `plot_sanjuan_*_pub.py` and `plot_teneguia_*_pub.py` produce the single-flow publication plots. In `pyflowgo/`, `plot_viscosity_yield.py`, `plot_viscosity_yield_1949.py` and `replot_resolution_comparison.py` redraw the PyFLOWGO profiles. The maps are made from the repository root with `python make_pub_maps.py` (see the next section).

## Digital elevation models

The two DEM tiles are third-party data and are not included in this repository (they exceed GitHub's file-size limit). They are the 2 m MDT02 REGCAN95 tiles `MDT02-REGCAN95-HU28-1085-3-4-COB2` and `MDT02-REGCAN95-HU28-1087-1-2-COB2`, available from the Centro de Descargas of Spain's Instituto Geográfico Nacional (https://centrodedescargas.cnig.es/). Download them as GeoTIFFs and place them in the repository root, next to `make_pub_maps.py`:

```text
MDT02-REGCAN95-HU28-1085-3-4-COB2.tif
MDT02-REGCAN95-HU28-1087-1-2-COB2.tif
```

`make_pub_maps.py` and `build_mapbox_map.py` read these files. Nothing else in the repository needs them. `*.tif` is git-ignored.

## Mapbox token

`build_mapbox_map.py` builds an interactive map with Mapbox GL JS and reads its access token from the `MAPBOX_TOKEN` environment variable. It stops with an error if the variable is not set. Create a public token at https://account.mapbox.com/, then:

```sh
export MAPBOX_TOKEN=pk.your-token-here      # PowerShell: $env:MAPBOX_TOKEN = "pk.your-token-here"
python build_mapbox_map.py
```

The generated `la_palma_methodology_map.html` embeds the token, so it is git-ignored and should not be committed.

## Credit and citation

Flow simulations use PyFLOWGO, developed by Magdalena Oryaëlle Chevrel and Jérémie Labroquère. The copy in `pyflowgo/` is the upstream source (https://github.com/pyflowgo/pyflowgo, commit `b4297d5`) with no upstream file changed; this study adds scripts, configurations and results alongside it. The upstream PDFs, CI configuration and IDE files were left out. If you use PyFLOWGO, please cite:

Chevrel M.O., Labroquère J., Harris A. and Rowland S. (2018). PyFLOWGO: an open-source platform for simulation of channelized lava thermo-rheological properties. *Computers & Geosciences* 111, 167–180. https://doi.org/10.1016/j.cageo.2017.11.009

## Licence

PyFLOWGO is distributed under its own licence; see `pyflowgo/LICENSE`, `pyflowgo/COPYING.gpl3`, `pyflowgo/COPYING.lgpl3` and `pyflowgo/README.md`. The scripts written for this study are released under GPL-3.0 to match; the full text is in `LICENSE`. The DEM tiles remain under the terms set by the Instituto Geográfico Nacional.
