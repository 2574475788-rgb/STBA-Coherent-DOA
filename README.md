# ML_CFR — public code and numerical-data package

This repository is the compact public package supporting the manuscript on
stability-adaptive multi-lag spatial-factor recovery (ST-BA) for two fully
coherent sources on sparse arrays.

The repository intentionally contains **only**:
- core MATLAB estimator code;
- the final experiment drivers used for the reported numerical studies;
- final MAT/CSV numerical results underlying the manuscript;
- bootstrap/statistical post-processing code;
- a small number of qualification-gate files required by the supporting studies;
- documentation and integrity manifests.

Publication-layout and figure-rendering scripts are intentionally omitted.
The paper figures are derived from the included numerical result files.

## Quick start

```matlab
cd('path/to/repository')
clear functions
rehash path

setup_project
run_clean_pipeline('check')
```

### Recompute statistical post-processing from stored paired trials

```matlab
run_clean_pipeline('postprocess')
```

This does not regenerate observations or rerun the estimators.

### Fresh rerun of the final paired Monte Carlo

```matlab
run_clean_pipeline('core')
```

To protect the fixed manuscript results, an existing official result MAT is
not overwritten. For a full fresh rerun, first copy the repository and remove
the corresponding `results/submission_final/submission_final_*.mat` files in
that copy.

## Repository layout

```text
src/                         Core ST-BA, CML, TMR-VAI, FBSS-MUSIC, generators, CRB
experiments/submission_final Final paired Monte Carlo and paired-bootstrap analysis
experiments/formal           ULA-transfer and unequal-amplitude supporting studies
experiments/qualification    Pre-specified sanity gates for supporting studies
experiments/reproduction     Mechanism/mismatch supporting experiments
results/submission_final     Final trial-level MAT files and summary CSV files
results/frozen               Fixed mechanism/mismatch numerical results
results/candidate            Fixed ULA/unequal-amplitude numerical results
results/qualification        Supporting-study gate results
```

## What is deliberately not included

- publication figure-layout scripts;
- historical patch/development branches;
- failed Peng/MESA development code;
- diagnostic PNGs and obsolete plots;
- Monte Carlo checkpoints;
- old patch manifests and hot-fix notes;
- third-party development packages that were not promoted to the final paper curves.

## Data availability / archival release

For journal submission, use GitHub as the live development repository and
archive a tagged release in a persistent repository such as Zenodo.
Cite the archived DOI in the manuscript's Data Availability statement.
The GitHub URL may also be provided as the development repository.

## Integrity

See:
- `FILE_GUIDE_ZH.md` — per-file explanation;
- `DATA_FILES_ZH.md` — what each result MAT/CSV contains and why it is included;
- `CLEANUP_NOTES_ZH.md` — what was intentionally removed;
- `MANIFEST.csv` and `SHA256_MANIFEST.txt` — file hashes.
