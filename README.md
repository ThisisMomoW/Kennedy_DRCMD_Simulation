# Kennedy and DRCMD for a missing binary exposure

This simulation study compares four methods for estimating the marginal causal risk difference. The methods are Kennedy ratio GAM (gamma + beta), Kennedy direct GAM (gamma + m), and two DRCMD variants. Only exposure A is missing, and R = 1 indicates that it is observed.

## Scenarios

- [Observation heterogeneity](scenarios/01_observation_heterogeneity/README.md): s=1,3,5.
- [Outcome slope](scenarios/02_outcome_slope/README.md): b=0.5,1.5,3.
- [Quadratic outcome logit](scenarios/03_nonlinear_outcome/README.md): delta=2; s=1,5.

Every scenario uses independent Uniform(-1,1) covariates, n=500 and 1000,
80% population missingness, and 200 attempted datasets per setting. There are
3,200 datasets; the original three-method run contains 9,600 method records. Calibration is deterministic and fixed
across replications; realized missing counts remain random.

## Dependencies and commands

Run all commands from this directory. R 4.4.0 Patched, SuperLearner 2.0-40 and gam 1.22-4 reproduce the numerical environment. DRCMD 0.1.0 is the official commit `76dc26abdb921025a2710b2826e96efe9b9df800`. The actual loaded versions and paths are in [software/versions.csv](software/versions.csv); the session information is in [software/session_info.txt](software/session_info.txt). The installed DRCMD package and official source are included. A LaTeX installation, Pandoc, rmarkdown, knitr, ggplot2 and patchwork are needed for PDF rendering.

```bash
# On JHPCE:
module load conda_R/4.4
# Optional dependency installation in a separate checkout:
Rscript scripts/install.R
# Prepare the fixed settings, attempt plan, streams and software manifest:
Rscript run.R init
Rscript scripts/validate.R
Rscript run.R smoke
# Sequential execution (completed checkpoints are skipped):
Rscript run.R 1     # or 2, 3, or all
# Cluster execution, instead of sequential execution:
bash scripts/submit.sh
# After all jobs finish:
Rscript scripts/analyze.R
Rscript scripts/render.R
```

For the fourth method, reuse the saved datasets and run:

```bash
# Verify one existing replicate from each scenario:
Rscript run.R eem verify
# Submit the full EEM run:
bash scripts/eem/submit.sh all
# After all EEM jobs finish, summarize and render the four-method comparison:
Rscript run.R eem summarize
Rscript scripts/render_four_method.R
Rscript scripts/verify_four_method_pdf.R
```

## Statistical implementation

[R/kennedy_ratio.R](R/kennedy_ratio.R) fits the shared observation, exposure imputation and gamma functions, Gaussian beta regressions, and the common EIF.
[R/kennedy_direct.R](R/kennedy_direct.R) fits separate quasibinomial logit GAMs using lambda weights on every outer-training row. Both use `gam::s` smooths with two degrees of freedom. Shared regressions use single-library SuperLearner with five internal CV folds. Direct outcome GAMs have fixed defaults and no additional CV. Only kappa and gamma denominators are truncated to [0.01,0.99]; beta/gamma, ATEs and intervals are never clipped. Direct fitting uses the same upstream fits and scores. SEs use the SD of individual arm-score differences divided by sqrt(n).

[R/drcmd.R](R/drcmd.R) calls the exported package procedure on observed variables with SL.gam throughout, k=5, cv_folds=5, cutoff=.01, `eem_ind = FALSE` and `tml = FALSE`. The fourth method, `drcmd_eem`, changes only `eem_ind` to `TRUE`, reusing the saved datasets, original DRCMD fitting seeds and outer folds. Its native `psi_hat_ate` arm difference and corresponding SE give the normal interval. DRCMD's native outer folds are captured and reused by both Kennedy versions. DRCMD fits its own weighted nuisance regressions. Thus this compares implemented procedures with different nuisance targets and fitting steps.

The seed is 20261123. 

## Results and report

Each `scenarios/*/results/` contains replicate CSVs, summaries, paired MSE comparisons, diagnostics, checkpoints and figures. Combined files are in [results/](results/). All finite extreme estimates remain in bias/RMSE. Coverage is calculated among valid intervals, with its denominator displayed separately from attempted datasets and finite points. Failures and warnings are explicit.

The [EEM extension README](extensions/drcmd_eem/README.md) documents the fourth method and its workflow. The four-method comparison outputs and PDF are saved separately under [results/drcmd_eem_v1/](results/drcmd_eem_v1/).


The current four-method report is [Kennedy_DRCMD_Simulation_Report.Rmd](reports/Kennedy_DRCMD_Simulation_Report.Rmd) and [its PDF](reports/Kennedy_DRCMD_Simulation_Report.pdf). Render it from saved results with `Rscript scripts/render_four_method.R`.

The completed four-method comparison contains 12,800 records on all 3,200 datasets. All 3,200 EEM fits have finite estimates and valid intervals. See the [64-row summary](results/drcmd_eem_v1/comparison/summary.csv) and [paired DRCMD/EEM MSE comparisons](results/drcmd_eem_v1/comparison/paired_mse.csv). Original three-method results remain unchanged.
