# Kennedy and DRCMD for a missing binary exposure

This standalone R study compares exactly three procedures: Kennedy ratio GAM
(gamma + beta), Kennedy direct GAM (gamma + m), and DR-CMD. The target is the
marginal causal risk difference. Only exposure A is missing; R=1 means observed.

## Scenarios

- [Observation heterogeneity](scenarios/01_observation_heterogeneity/README.md): s=1,3,5.
- [Outcome slope](scenarios/02_outcome_slope/README.md): b=0.5,1.5,3.
- [Quadratic outcome logit](scenarios/03_nonlinear_outcome/README.md): delta=2; s=1,5.

Every scenario uses independent Uniform(-1,1) covariates, n=500 and 1000,
80% population missingness, and 200 attempted datasets per setting. There are
3,200 datasets and 9,600 method records. Calibration is deterministic and fixed
across replications; realized missing counts remain random.

## Dependencies and commands

Run all commands from this directory. R 4.4.0 Patched, SuperLearner 2.0-40 and gam
1.22-4 reproduce the numerical environment. DR-CMD 0.1.0 is the official commit
`76dc26abdb921025a2710b2826e96efe9b9df800`. The actual loaded versions and paths
are in [software/versions.csv](software/versions.csv); the session information
is in [software/session_info.txt](software/session_info.txt). The installed
DR-CMD package and official source are included. A LaTeX installation, Pandoc,
rmarkdown, knitr, ggplot2 and patchwork are needed for PDF rendering.

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

## Statistical implementation

[R/kennedy_ratio.R](R/kennedy_ratio.R) fits the shared observation, exposure imputation and gamma functions, Gaussian beta regressions, and the common EIF.
[R/kennedy_direct.R](R/kennedy_direct.R) fits separate quasibinomial logit GAMs using lambda weights on every outer-training row. Both use `gam::s` smooths with two degrees of freedom. Shared regressions use single-library SuperLearner with five internal CV folds. Direct outcome GAMs have fixed defaults and no additional CV. Only kappa and gamma denominators are truncated to [0.01,0.99]; beta/gamma, ATEs and intervals are never clipped. Direct fitting uses the same upstream fits and scores. SEs use the SD of individual arm-score differences divided by sqrt(n).

[R/drcmd.R](R/drcmd.R) calls the exported package procedure on observed variables with SL.gam throughout, k=5, cv_folds=5, cutoff=.01, EEM and TMLE disabled. Its native `psi_hat_ate` arm difference and corresponding SE give the normal interval. DRCMD's native outer folds are captured and reused by both Kennedy versions. DR-CMD fits its own weighted nuisance regressions. Thus this compares implemented procedures with different nuisance targets and fitting steps.

The seed is 20261123. 

## Results and report

Each `scenarios/*/results/` contains replicate CSVs, summaries, paired MSE
comparisons, diagnostics, checkpoints and figures. Combined files are in
[results/](results/). All finite extreme estimates remain in bias/RMSE. Coverage is calculated among valid intervals, with its denominator displayed separately from attempted datasets and finite points. Failures and warnings are explicit.

