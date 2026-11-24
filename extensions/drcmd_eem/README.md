# DRCMD with native EEM

This extension adds `drcmd_eem`. Its exported `drcmd::drcmd()` call differs from
`drcmd` only in `eem_ind=TRUE`. The package remains 0.1.0 at commit
`76dc26abdb921025a2710b2826e96efe9b9df800`, with SL.gam throughout, k=5,
cv_folds=5, cutoff=.01, tml=FALSE and parallel=FALSE. It uses the native
`psi_hat_ate` arm difference, corresponding SE, and estimate +/- 1.96*SE.

Each EEM attempt reads the original observed dataset, DRCMD seed and native
folds from its checkpoint. It resets that exact fitting seed before calling the
package and verifies the newly captured native folds. DRCMD fits its own
nuisance models. Package-native RNG consumption is preserved; shared seeds do
not impose additional internal splits or force reuse of fitted models.

The EEM branch uses the package's centered pseudo-outcomes, Gaussian GAMs, and
weights `(R/kappa_hat - 1)^2`. No custom EEM loss is substituted. Read-only
tracing distinguishes the arm1, arm0 and unused direct-contrast regression roles.

## Run from the project root

```bash
module load conda_R/4.4
Rscript run.R eem init
Rscript run.R eem verify
# EEM-only sequential execution:
Rscript run.R eem 1       # alternatively 2, 3, or all
# Or submit a resumable cluster array:
bash scripts/eem/submit.sh all
# Alternatively submit only one scenario:
bash scripts/eem/submit.sh 1
bash scripts/eem/submit.sh 2
bash scripts/eem/submit.sh 3
# When all EEM jobs are finished:
Rscript run.R eem summarize
Rscript scripts/render_four_method.R
Rscript scripts/verify_four_method_pdf.R
```

For a progress summary or a completed single scenario, use
`Rscript run.R eem summarize --partial`. This marks the comparison **PARTIAL**,
records unattempted EEM datasets as pending, and does not count pending work as
fit failures. A full PDF can be rendered only after all 3,200 EEM attempts.

The array uses the same 10-dataset tasks and resource defaults as the original
workflow. `SIM_TASKS`, `SIM_CONCURRENT_TASKS`, `SIM_PARTITION`, `SIM_MEMORY`,
`SIM_WALLTIME`, and `SIM_EXCLUDE` are configurable. Resubmitting a task skips
its saved EEM checkpoints, including recorded fitting failures. It neither
regenerates data nor reruns the other three methods. Execution starts are logged.

## Separate outputs

All new outputs are under `results/drcmd_eem_v1/`:

- `checkpoints/`: native fit objects, warnings, data/fold hashes, seeds and audits.
- `eem_replicates.csv`: completed EEM attempts only.
- `comparison/replicates.csv`: original three method records plus EEM attempts.
- `comparison/summary.csv`: 64 planned setting-method rows; attempted, pending,
  finite point, valid CI, failure counts, bias/RMSE/coverage and Monte Carlo SEs.
- `comparison/matched_*`: four methods restricted to datasets with an EEM attempt.
- `comparison/paired_mse.csv`: DRCMD versus EEM, restricted to jointly finite
  estimates; differences are squared error(DRCMD) minus squared error(EEM), so
  positive values favor EEM. MCSE is SD of paired differences / sqrt(paired n).
- `comparison/paired_squared_errors.csv`, `failures.csv`, `warnings.csv`, and
  `identity_checks.csv`: detailed evidence, including unused contrast warnings.
- `comparison/{scenario}/`: scenario-specific four-method tables.
- `report/Kennedy_DRCMD_Simulation_Report.pdf`: separate full four-method report.
- `verification/`: fixed one-replicate-per-scenario checks and preserved values.

Bias and RMSE include all finite estimates and extremes; coverage uses only
valid intervals. With one paired replicate, a difference is available but its
Monte Carlo SE is not. No inference about method ranking follows from the three
verification fits.

The original `R/` source, source-verification manifest, checkpoints, results and
PDF are read-only inputs. The extension has a separate source manifest and
checksums; it does not update or bypass the original source verification.

The full run is complete: 3,200 EEM attempts, 3,200 finite estimates and 3,200
valid intervals. All 320 tasks in array 37296696 completed with exit code 0.
The four-method report validates 12,800 records and 64 summary rows against
the saved results; the original 9,600 method records are unchanged. The scheduler
audit is `results/drcmd_eem_v1/scheduler_accounting.txt`.
