# Outcome slope

X1 and X2 are independent Uniform(-1,1). A, Y and R are binary, with R=1
meaning A is observed. Only A is missing.

```text
logit g = 0.4 X1
logit m_a = -0.5 + 0.5 a + b X1
logit kappa = c_b + 0.4 X1 + 0.2 Y
```

Settings: b=0.5,1.5,3; n=500,1000; 200 attempted datasets per setting;
80% target population missingness. True marginal risk difference: 0.120114506958277, 0.104555441174916, 0.0751434473300117, respectively.
The calibrated intercepts and numerical integration results are in [settings.csv](settings.csv).
All three estimators share each observed dataset and native five outer folds.
X2 is supplied to every estimator, including when irrelevant to generation.

From the project root:

```bash
Rscript run.R 2
Rscript scripts/analyze.R
Rscript scripts/render.R
```

The scenario generator is [dgp.R](dgp.R); common fitting code is in ../../R.
Results, including every finite extreme estimate and failed fit, are in [results/](results/).
