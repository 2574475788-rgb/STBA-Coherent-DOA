# Numerical data files

The repository includes the numerical data underlying the manuscript so that
the reported performance values can be checked without rerunning long Monte
Carlo experiments.

## Final paired Monte Carlo

- `results/submission_final/submission_final_snr.mat`
  - trial-level paired results for the SNR sweep.
- `results/submission_final/submission_final_snapshot.mat`
  - trial-level paired results for the snapshot-count sweep.
- `results/submission_final/submission_final_separation.mat`
  - trial-level paired results for the two-source separation sweep.
- `results/submission_final/submission_final_keypoint.mat`
  - 5,000-trial low-SNR key operating point; source for the tail/stability audit.
- `results/submission_final/submission_final_runtime.mat`
  - matched implementation-runtime records.

Each corresponding `*_summary.csv` is a human-readable summary of the MAT
file and is included for convenient inspection.

## Statistical post-processing

- `submission_final_paired_STBA_vs_CML_CI.csv`
  - paired percentile-bootstrap intervals for ST-BA minus continuous CML.
- `submission_final_keypoint_MC_convergence.csv`
  - prefix Monte Carlo audit of the fixed 5,000-trial record.
- `submission_final_postprocess.mat`
  - MATLAB archive of the post-processing output.
- `FigS2_mc_stability_bootstrap.mat`
  - bootstrap stability statistics derived from the stored key-point trials.
    The filename is retained because it corresponds to the supplementary
    stability analysis; no figure-layout script is included in this repository.

## Supporting studies

- `results/frozen/v1g_stability_selector_results.mat`
  - mechanism/pilot-selection study.
- `results/frozen/v3c_robustness_formal_results.mat`
  - temporal/spatial mismatch study.
- `results/candidate/v8f_ula_reduced_snr_results.mat`
  - ULA SNR transfer control.
- `results/candidate/v8g_ula_reduced_snapshot_results.mat`
  - ULA snapshot transfer control.
- `results/candidate/unequal_amplitude_formal_results.mat`
  - unequal-amplitude operating-boundary study.

## Qualification records

- `results/qualification/v8e_ula_reduced_sanity_results.mat`
- `results/qualification/unequal_amplitude_sanity_results.mat`

These are pre-specified gates used by the supporting experiment scripts.

- `results/submission_final/external_baseline_qualification.csv`
  records the frozen qualification status of the independent Peng/MESA
  implementations. It is not a performance ranking.
