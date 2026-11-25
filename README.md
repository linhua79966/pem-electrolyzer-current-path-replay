# C19 — double-platform response and replay reproduction software (version 1.1.0)

Accompanies *Current-path effects and residual compensation in PEM electrolyzer outlet-flow estimation*, academic manuscript V0.28. Creators: Xianhuan Wu, Suhua Li, Cenglin Yao, Lesheng Liu.

**Status: author-authorized release; archived on Zenodo on 2026-11-25.** MIT applies to eligible original software only. No upstream data license is granted. No NLR raw data, NLR fitted coefficients or NLR-derived result tables are bundled. The companion dataset contains 35 Jülich aggregate payloads.

- Software V1.1.0 (MIT): [10.5281/zenodo.24531336](https://doi.org/10.5281/zenodo.24531336).
- Companion Jülich Dataset V1.0.0 (CC BY 4.0 for eligible original contributions; upstream CC0 preserved): [10.5281/zenodo.24532302](https://doi.org/10.5281/zenodo.24532302).
- Upstream Jülich dataset: [10.26165/JUELICH-DATA/FEKKSS](https://doi.org/10.26165/JUELICH-DATA/FEKKSS).

The Zenodo Software ZIP remains the immutable reference archive (SHA256 `6f8942edb1091c005a2001bb7204f10126cec98906de61e51db9a602022a5918`). This repository incorporates post-deposit DOI metadata only; the computational files are unchanged from that archive.

V1.1 provides six separate local calculation entry points in `src/reproduce.py`. It uses original numerical functions extracted from the historical code, adapted file readers and a new orchestration layer. The new pipeline has been tested on synthetic data and the supplied aggregate audit has been run on the Jülich candidate dataset. **The six new entry points have not been executed on the original research signals; complete equivalence to the published experiment is not certified.** The package does not rebuild every historical sensitivity analysis, waveform, figure or supplementary table. It is a research diagnostic, not a controller. Historical one-use execution tokens, project markers, raw data, setpoint files and production runners are excluded.

## 1. Safe initial check (no original source data)

Extract the software and dataset into separate sibling directories. Use Python 3.13 (tested with 3.13.5 on Linux); install the pinned direct dependencies into a fresh virtual environment when a package index is available. On Windows, from the software directory in the VS Code terminal:

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -B verify_candidate.py ../C19_DATASET_V1_0_CANDIDATE ../C19_ACCEPTANCE_NEW
```

On Linux/macOS, use `python3 -m venv .venv` and `.venv/bin/python` for the last two commands. The verification entry point checks manifests, runs the synthetic tests, audits the 35 supplied aggregate payloads and recreates the two timing plots. Use a new output directory. This check never calls a raw-data entry point. The extracted-package tests passed on Linux and in a user-returned Windows run (Python 3.13.14; 16 test groups, 1,363 aggregate checks and two pixel-identical decoded figures). A separate clean-install log was not supplied. No wheels or third-party source are distributed.

## 2. Independent researcher workflow (optional; not a task for the original frozen C19 experiment)

Read `acquisition/ACQUISITION.md` and `acquisition/SOURCE_SPEC.json` before obtaining files. Obtain source data separately and lawfully. Place the explicitly listed original files into independent source directories (not the historical `C:\19paper` production directory). No network access or automatic acquisition is performed by the program. Identity checks require the exact original byte contents; if the provider supplies a different revision, stop and record the mismatch instead of editing the manifest to force a match.

Each raw stage has `--preflight-only` (hashes only, no output folder created), and `--local-reproduction` (an explicit request to perform independent local calculations, not an authorization for disclosure). Examples use `python` from the active venv. The `--out` path must be new and outside the source directory. For the original authors, continue to use the archived result-only evidence; do not execute these raw stages on the frozen study data.

```text
python -B src/reproduce.py nlr-fit --source ../nlr_characterization --out ../nlr_fit_new --preflight-only
python -B src/reproduce.py nlr-fit --source ../nlr_characterization --out ../nlr_fit_new --local-reproduction
python -B src/reproduce.py nlr-replay --source ../nlr_trajectories --models ../nlr_fit_new/MODELS.json --out ../nlr_replay_new --local-reproduction
python -B src/reproduce.py nlr-lag --source ../nlr_trajectories --out ../nlr_lag_new --local-reproduction
python -B src/reproduce.py juelich-fit --source ../juelich_development --out ../juelich_fit_new --local-reproduction
python -B src/reproduce.py juelich-timing --source ../juelich_development --models ../C19_DATASET_V1_0_CANDIDATE/data/legacy/JUELICH_FROZEN_MODELS.json --out ../juelich_timing_new --local-reproduction
python -B src/reproduce.py juelich-replay --source ../juelich_evaluation --models ../C19_DATASET_V1_0_CANDIDATE/data/legacy/JUELICH_FROZEN_MODELS.json --out ../juelich_replay_new --local-reproduction
```

Jülich original fixed models are supplied only by the companion data record. The NLR replay command uses `MODELS.json` generated in the separate characterization stage. Do not select new parameters using evaluation data. FIT/SELECT/DEV_DIAGNOSTIC remain disjoint; 2022–2023 evaluation is fixed-model only. Separate output status is `LOCAL_CALCULATION_COMPLETED_UNCOMPARED` until reports are compared. File identities are matched against historical local copies, not a fresh provider certification. The output contains only local summaries, not raw traces or per-origin predictions. NLR output must not be redistributed without verifying applicable terms.

### Method and support conventions

NLR: explicit clocks, nonnegative actual and commanded current, contiguous 1 s input segments, 300-row initialization, all four horizons (1,5,10,30 s) share origins requiring finite targets at every horizon. H holds the command after consuming actual current at the origin; C consumes future command, A consumes actual, U holds actual. M0 endpoint is diagnostic, not deployable. All one-step M1 paths coincide. The time-alignment diagnostic uses all valid current segments, common anchors for lags −10…+10 s, raw level and first increments; positive lag means command leads actual. Lag support differs from flow-replay support. The 23-trajectory manifest omits three blocked archive entries and two unrecovered archives, and excludes source measurement values. The original per-collection overlap check is retained; the characterization time interval is obtained locally when fitting. Lag alone does not use characterization. Full original preflight metadata and some duplicate checks across Jülich development/evaluation collections are not included; equality to the archived results must be checked separately.

Jülich: row 51 tag header, row 52 units, rows 53+; total measured stack current and native FI-24 m³/h. Nonnegative individual currents, consistent counter/clock and gaps no larger than 1.5 times the archived modal interval. No unit conversion or imputation; missing target is excluded only from eligible metric support. Left ZOH, 300 s initialization, minimum 30 eligible records/file. At 10/30 s both horizons share origins, nearest earlier endpoint on ties within ±0.5 s. File macro differs from pooled weighting. The 126-grid comparison uses development data only and the supplied legacy model references; never choose a model using 2022–2023 results.

### Comparison and figures

Compare common keys and selected numerical fields against companion data or your legitimately held historical reports. Example:

```text
python -B src/compare_reports.py ../juelich_replay_new/REPLAY_METRICS.csv ../C19_DATASET_V1_0_CANDIDATE/data/legacy/RECOMPUTED_REPLAY_METRICS.csv --keys file_id,campaign,h --columns n,rmse_H,rmse_A,mse_R,mse_F,interaction,mse_total --output ../juelich_comparison.json
python -B src/compare_reports.py ../juelich_timing_new/CANDIDATES.csv ../C19_DATASET_V1_0_CANDIDATE/data/04_REPORTS/CANDIDATES.csv --keys index --columns family,tau_s,delay_s,beta,fit_macro_mse,select_macro_rmse,eligible,selected --output ../timing_comparison.json
python -B src/plot_local_reports.py nlr-lag ../nlr_lag_new ../nlr_lag_plots_new
python -B src/plot_local_reports.py nlr-replay ../nlr_replay_new ../nlr_replay_plots_new
python -B src/plot_local_reports.py juelich-replay ../juelich_replay_new ../juelich_plots_new
```

Compare NLR `DECOMPOSITION_METRICS.csv` against the two locally held Phase2D/2E reports, selecting common experimental rows or source-pooled rows, and `LAG_CURVES_21.csv` against your legitimately held lag reference after matching field names. No NLR numerical reference is bundled. `compare_reports.py` demands exact unique key coverage and compares counts exactly, floats with explicit tolerance (default abs/rel 1e-9); a selected-column comparison is not proof of end-to-end reproduction. Local diagnostic plots are not pixel replicas of the manuscript. Two timing plotting scripts retain the earlier plot logic. More detailed figure reconstruction and historical threshold sensitivity variants are not implemented as full raw pipelines.

## 3. Attribution and license

See `acquisition/ACQUISITION.md`, `CODE_PROVENANCE.csv`, `NEW_CODE_PROVENANCE.csv`, `LICENSE`, `RELEASE_STATUS.txt` and `THIRD_PARTY_DEPENDENCIES.md`. AI assisted development and manuscript preparation; the authors reported completing their review and accepting responsibility. The MIT license is limited to eligible original software, not upstream files, dependency licenses or the separate dataset. Cite the four creators, title, version, and the version-specific software DOI 10.5281/zenodo.24531336. The companion dataset DOI is 10.5281/zenodo.24532302.
