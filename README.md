# GD-R: A Validity-Based General Dominance Method for Predictor Relative Importance

General dominance analysis (GDA) assesses predictors' relative importance by partitioning the explained variance (R²) of a regression model based on their average incremental contributions across all possible predictor subsets (Azen & Budescu, 2003; Budescu, 1993). In personnel selection, however, predictive validity (R) is often the more relevant metric because it is directly related to expected selection utility. GD-R is a validity-based GDA method that partitions the model R, rather than R², into additive predictor-specific contributions. Each predictor's GD-R value is its average incremental validity (ΔR) across all possible predictor subsets, and the GD-R values across predictors sum to the model R.

This repository provides:

- **`gdr.R`**, the analytical engine (base R only). It conducts GD-R, traditional variance-based GDA (GD-R²), GDA based on semi-partial correlations (GD-sr), and relative weight analysis (RWA) from a correlation matrix, together with conditional incremental validity analyses, user-specified operational predictor weights, and bootstrap confidence intervals (CIs).
- **An R Shiny application** that implements the same engine through a point-and-click interface: <https://blindreview-gdr.share.connect.posit.cloud>.

## Contents

- [Methods](#methods)
- [Using the R Shiny application](#using-the-r-shiny-application)
- [Using gdr.R in R](#using-gdrr-in-r)
- [Function reference](#function-reference)

## Methods

| Code | Method | Contribution of predictor Xᵢ to a subset model S | GD values sum to |
|---|---|---|---|
| `"R"` | GD-R | incremental validity, ΔR = R(Y·Xᵢ,S) − R(Y·S) | model R |
| `"R2"` | GD-R² (traditional GDA) | incremental variance, ΔR² = R²(Y·Xᵢ,S) − R²(Y·S) | model R² |
| `"spcorr"` | GD-sr | semi-partial correlation, sr = √ΔR² | not additive |
| `"rwa"` | RWA (Johnson, 2000) | not based on predictor subsets | model R² |

For GD-R², GD-R, and GD-sr, a predictor's GD value is obtained with the two-step averaging rule of GDA: its incremental contributions are first averaged across all subset models containing the same number, k, of other predictors (k = 0, …, p − 1), and these k-level averages are then averaged equally (Azen & Budescu, 2003). This averaging rule is mathematically equivalent to the Shapley value (Shapley, 1953). Because the GD values are additive, each predictor's contribution can also be expressed as a percentage-based GD value (e.g., %GD-R = GD-R value / model R).

## Using the R Shiny application

The application is available at <https://blindreview-gdr.share.connect.posit.cloud>. Inputs are entered in the sidebar (Steps 1–6), and results are displayed on the right. No analysis is conducted until **Calculate Results** is clicked.

### Step 1. Input type and data upload

Select the input type and upload the data as a `.csv` file.

- **Correlation / covariance matrix.** The file must contain a full, square, symmetric matrix with variable names in the first row. A first column containing variable names is optional. A covariance matrix is automatically converted to a correlation matrix. An example input matrix:

  ```
  ,Y,X1,X2,X3,X4
  Y,1,0.5,0.6,0.7,0.75
  X1,0.5,1,0.4,0.4,0.75
  X2,0.6,0.4,1,0.4,0.75
  X3,0.7,0.4,0.4,1,0.75
  X4,0.75,0.75,0.75,0.75,1
  ```

  The **sample size (N)** is required only for bootstrap CIs. For a meta-analytic correlation matrix, the harmonic-mean sample size can be used.

- **Raw data (rows = cases).** The file must contain one row per case and one column per variable. Non-numeric columns (e.g., IDs) are ignored, and cases with missing values on the selected variables are excluded (listwise deletion).

The **Input** tab displays the uploaded matrix (or, for raw data, the correlations among the numeric variables) to facilitate verification of the data structure.

### Step 2. Variable selection

- **2a. DV (Y):** the criterion (dependent variable).
- **2b. IVs (X):** one or more predictors (up to 15; pairwise comparisons require at least two). The selected DV is automatically removed from the list of available predictors.

### Step 3. Metric specification

Select one or more relative importance methods: GD-R², RWA, GD-R, and GD-sr. GD-R² and GD-R are selected by default, and multiple methods can be analyzed simultaneously for comparison. **Include Std. Betas?** (selected by default) reports the standardized regression coefficients of the full model.

### Step 4. Bootstrap confidence intervals (optional)

Select **Compute bootstrap CIs** to obtain percentile CIs for all reported statistics: the raw and percentage-based GD values, the standardized regression coefficients, the model R² and R, and the incremental validity values, as well as rank intervals for the predictor ranks.

- **Replications:** 1,000, 2,000, 5,000, or 10,000.
- **Confidence level:** e.g., .95.
- **Random seed:** the same seed reproduces the same intervals.

For raw data, cases are resampled with replacement (nonparametric bootstrap). For a correlation matrix, replicate matrices are drawn from the Wishart distribution with N − 1 degrees of freedom, which is equivalent to sampling N cases from a multivariate normal population with the input correlation matrix (parametric bootstrap). Rank intervals are defined by the percentiles of the bootstrap rank distribution (1 = most important). With more than 10 predictors, the bootstrap may take several minutes; a progress bar is displayed.

### Step 5. Predictor weights (optional)

By default, each predictor subset is regression-weighted, as in traditional GDA. To evaluate predictors under an operational weighting scheme, select:

- **Equal (unit) weights**, or
- **Custom weights**, entered for each predictor. The weights apply to standardized predictor scores; only their relative sizes matter, and negative weights are allowed.

Under this option, the user-specified predictor weights remain fixed across all predictor subsets, so each predictor retains its specified weight in every subset model in which it appears. The resulting GD-R², GD-R, and incremental validity values therefore represent contributions conditional on that weighting scheme. Weights that are re-optimized separately for each predictor subset do not fit this assumption; a single fixed set of weights (e.g., full-model Pareto-optimal weights) should be entered instead. Because adding a valid predictor may decrease composite validity under alternative weighting schemes (Murphy, 2019; Sackett et al., 2017), GD values can be negative. Only GD-R² and GD-R are available under this option; RWA and GD-sr are regression-based and are therefore not computed. If RWA and/or GD-sr are the only selected methods, the application stops and asks for GD-R² and/or GD-R to be selected.

### Step 6. Baseline predictors (optional)

Select one or more baseline predictors (at least one predictor must remain outside the baseline) to conduct conditional incremental validity analyses. The results are added to the main table:

- **ΔR over baseline** and **ΔR² over baseline:** the incremental validity and incremental variance of each remaining predictor when it is added to the baseline model. The **Total** column reports the increment when all remaining predictors are added together.

The R² and R of the baseline model are reported above the table. Under user-specified predictor weights, every model is scored with those weights.

### Output

**Relative importance.** A summary line reports the criterion, N, and the model R² and R (under user-specified weights, the composite R² and R together with the regression-weighted values for comparison). The results table reports, for each selected method, the raw GD values, the percentage-based GD values (% Total), and the predictor ranks (1 = most important), with bootstrap CIs in brackets when requested. Two display options change the table without rerunning the analysis:

- **Table layout:** *Vertical view* (predictors in rows; default) or *Horizontal view* (predictors in columns).
- **Decimal places:** 2 to 5 (default: 3; percentage-based values are reported with two fewer decimal places).

**Pairwise differences** (requires bootstrap CIs). For a selected method, each pair of predictors is compared in terms of the difference between their GD values, its bootstrap CI, whether the CI excludes zero ("Differs?"), and the proportion of bootstrap replications in which the first predictor's GD value is larger.

### Data validity and error messages

Dominance analysis relies on the inversion of correlation matrices, so the input matrix must be positive definite. The application checks this before any analysis and stops with a message when it is violated.

| Message | Recommended action |
|---|---|
| Matrix is Not Positive Definite | The input matrix is not a valid correlation matrix (e.g., it was assembled from different sources or heavily rounded). Verify the data source, or apply a matrix-smoothing technique (e.g., the nearest positive definite matrix) before rerunning the analysis. |
| A matrix file must be square / must be symmetric | Upload the full matrix (not a triangular matrix), with the same variables in rows and columns. |
| The matrix contains non-numeric entries | Select **Raw data** if the file contains raw data, or remove non-numeric cells. |
| Bootstrap CIs from a matrix need the sample size | Enter N in Step 1. |
| At most 15 predictors are supported | GDA evaluates all 2^p predictor subsets; reduce the number of predictors. |
| Too few complete cases on the selected variables | Raw data only: after listwise deletion, too few cases remain. Check for missing values in the selected variables. |
| The baseline must leave at least one predictor to add | Remove at least one predictor from the baseline (Step 6). |
| Enter a numeric weight for every predictor / At least one weight must be non-zero | Enter a weight for each predictor, and make sure the weights are not all zero (Step 5). |
| RW and GDsr … not computed | Expected under user-specified predictor weights (see Step 5). |
| With unit or custom weights, select GD-R² and/or GD-R | Under user-specified predictor weights, select GD-R² and/or GD-R in Step 3 (see Step 5). |

## Using gdr.R in R

`gdr.R` uses base R only and can be loaded directly from this repository:

```r
source("https://raw.githubusercontent.com/anonymosaurus/GDR/refs/heads/main/gdr.R")
```

The examples below use the following correlation matrix (criterion `Y`; predictors `X1`–`X4`):

```r
R <- matrix(
  c(1,    .50, .60, .70, .75,
    .50,  1,   .40, .40, .75,
    .60,  .40, 1,   .40, .75,
    .70,  .40, .40, 1,   .75,
    .75,  .75, .75, .75, 1),
  nrow = 5,
  dimnames = list(c("Y", "X1", "X2", "X3", "X4"), c("Y", "X1", "X2", "X3", "X4"))
)
```

To use a matrix stored in a `.csv` file (variable names in the first row), read it with `as_cor_matrix()`, which adds the row names required by the functions below and converts a covariance matrix to a correlation matrix:

```r
R <- as_cor_matrix(read.csv("matrix.csv"))   # if the file has a first column of variable names, remove it first
```

`calc_dominance()` accepts a correlation or covariance matrix directly. `analyze()` and `run_bootstrap()` require a correlation matrix with row and column names (e.g., from `as_cor_matrix()` or `cor()`); they do not convert a covariance matrix.

### GD values for one method: `calc_dominance()`

```r
gdr <- calc_dominance("Y ~ X1 + X2 + X3 + X4", R, objective = "R")
gdr$table
#>    Variable General_Dominance
#> X1       X1         0.1367050
#> X2       X2         0.1807880
#> X3       X3         0.2386030
#> X4       X4         0.2397265
gdr$total            # = model R
#> [1] 0.7958224
```

`objective` is `"R"` (GD-R), `"R2"` (GD-R²), or `"spcorr"` (GD-sr). Percentage-based GD values are `100 * gdr$table$General_Dominance / gdr$total`.

### Multiple methods simultaneously: `analyze()`

`analyze()` computes the selected methods, the standardized regression coefficients, and (optionally) conditional incremental validity from a correlation matrix, and returns all results as a single named vector. It requires the subset structure created by `make_subsets()`, which depends only on the number of predictors and can be reused.

```r
xs <- c("X1", "X2", "X3", "X4")
S <- make_subsets(length(xs))
est <- analyze(R, y = "Y", xs = xs, S = S, metrics = c("R2", "R", "spcorr", "rwa"))

est[paste0("R|", xs)]    # GD-R values
est[paste0("R%|", xs)]   # percentage-based GD-R values (%GD-R)
est[paste0("R#|", xs)]   # GD-R ranks (1 = most important)
est[paste0("R|", TOTAL)] # sum of the GD-R values (= model R)
est[c("model|R", "model|R2")]
```

Elements of the result (`m` is the method code, `X` a predictor name):

```text
beta|X              standardized regression coefficient of predictor X
m|X                 GD value (or RW value) of predictor X for method m
m|<total>           sum of the values across predictors (TOTAL is the string "<total>")
m%|X                percentage-based value of predictor X
m#|X                rank of predictor X (1 = most important)
model|R, model|R2   R and R² of the full predictor set under the specified weighting
reg|R, reg|R2       regression-based R and R² of the full predictor set
wt|X                user-specified weight of predictor X (user-specified weights only)
inc|step|stat       conditional incremental validity
                      step: base, add:X, or all
                      stat: R, R2, dR, or dR2
```

### Conditional incremental validity

```r
est <- analyze(R, "Y", xs, S, metrics = c("R2", "R"), baseline = "X4")
est[c("inc|base|R", "inc|add:X1|dR", "inc|add:X2|dR", "inc|add:X3|dR", "inc|all|dR")]
#>    inc|base|R inc|add:X1|dR inc|add:X2|dR inc|add:X3|dR    inc|all|dR
#>   0.750000000   0.005928946   0.002139805   0.028276484   0.045822426
```

`inc|base|R` is the R of the baseline model; `inc|add:X1|dR` is the incremental validity of X1 over the baseline; `inc|all|dR` is the incremental validity of all non-baseline predictors added together.

### User-specified operational predictor weights

Fixed predictor weights (in the order of the predictors, or named) are passed as `weights` in `calc_dominance()` or `w` in `analyze()`. The weights remain fixed across all predictor subsets. Only `"R2"` and `"R"` are available under this option.

```r
calc_dominance("Y ~ X1 + X2 + X3 + X4", R, objective = "R", weights = c(1, 1, 1, 1))$table

w <- c(X1 = 1, X2 = .5, X3 = .5, X4 = 1)
est_w <- analyze(R, "Y", xs, S, metrics = c("R2", "R"), w = w)
est_w[c("model|R", "reg|R")]   # composite validity vs. regression-based validity
#>   model|R     reg|R
#> 0.7452413 0.7958224
```

### Bootstrap confidence intervals

`run_bootstrap()` replicates `analyze()`, and `summarize_boot()` converts the replications into percentile CIs.

```r
metrics <- c("R2", "R")
est <- analyze(R, "Y", xs, S, metrics)

# From a correlation matrix and its sample size (parametric bootstrap)
boot <- run_bootstrap(
  type = "matrix", R = R, X = NULL, N = 300, B = 1000, seed = 8424,
  y = "Y", xs = xs, S = S, metrics = metrics, w = NULL, baseline = character(0),
  keys = names(est)
)
ci <- summarize_boot(est, boot, conf = .95)
cbind(estimate = ci$est, lower = ci$lo, upper = ci$hi)[paste0("R|", xs), ]

# From raw data (cases resampled); X is a numeric matrix with columns Y, X1, ..., X4
# boot <- run_bootstrap("data", R = cor(X), X = X, N = nrow(X), B = 1000, seed = 8424,
#                       y = "Y", xs = xs, S = S, metrics = metrics, w = NULL,
#                       baseline = character(0), keys = names(est))
```

The same seed reproduces the same intervals. Pairwise comparisons of predictor GD values can be obtained from the replications, e.g., `quantile(boot[, "R|X4"] - boot[, "R|X3"], c(.025, .975))`.

### Contributions at each k-level

`gd_levels()` returns each predictor's average incremental contribution at each k-level (conditional dominance; Azen & Budescu, 2003), and `gd_level_shares()` returns the percentage of each GD value attributable to each k-level. The first column (k = 0) is the share attributable to the predictor's contribution when entered alone.

```r
r <- subset_validities(R, "Y", xs, S)    # R of every predictor subset
round(gd_levels(r, S, "R"), 4)           # rows: predictors; columns: k = 0, ..., p - 1
round(gd_level_shares(r, S, "R")[, 1], 1) # % of each GD-R value attributable to k = 0
#> [1] 91.4 83.0 73.3 78.2
```

## Function reference

**Main functions**

| Function | Purpose |
|---|---|
| `calc_dominance(formula, cov_mat, objective = "R2", weights = NULL)` | GD values for one method, from a formula and a correlation or covariance matrix |
| `analyze(R, y, xs, S, metrics, w = NULL, baseline = character(0))` | all selected methods, standardized regression coefficients, and conditional incremental validity (`R` must be a correlation matrix) |
| `run_bootstrap(type, R, X, N, B, seed, y, xs, S, metrics, w, baseline, keys, tick = NULL)` | bootstrap replications of `analyze()` |
| `summarize_boot(est, boot, conf)` | percentile confidence intervals |

**Building blocks**

| Function | Purpose |
|---|---|
| `make_subsets(p)` | structure of all 2^p predictor subsets (bitmasks and Shapley weights) |
| `all_r2(Rxx, rxy, S)` | regression-based R² of every predictor subset |
| `all_r(Rxx, rxy, S, w = NULL)` | R of every predictor subset (regression-based, or under fixed weights) |
| `subset_validities(R, y, xs, S, w = NULL)` | `all_r()` from a full correlation matrix and variable names |
| `composite_r(Rxx, rxy, w)` | validity of a composite under fixed weights |
| `gd_from_r(r, S, metric)`, `gd_from_r2(r2, S, metric)` | GD values from the R or R² of every predictor subset |
| `gd_levels(r, S, metric)`, `gd_level_shares(r, S, metric)` | contributions at each k-level |
| `rwa_weights(Rxx, rxy)` | RW values (Johnson, 2000) |
| `is_pd(M)` | checks whether a matrix is positive definite |
| `as_cor_matrix(M)` | converts a covariance or correlation matrix to a correlation matrix with row and column names |
| `parse_formula()`, `check_vars()`, `check_metrics()`, `check_weights()` | input checks |

Each function is documented in `gdr.R`.

## References

- Azen, R., & Budescu, D. V. (2003). The dominance analysis approach for comparing predictors in multiple regression. *Psychological Methods, 8*(2), 129–148.
- Budescu, D. V. (1993). Dominance analysis: A new approach to the problem of relative importance of predictors in multiple regression. *Psychological Bulletin, 114*(3), 542–551.
- Johnson, J. W. (2000). A heuristic method for estimating the relative weight of predictor variables in multiple regression. *Multivariate Behavioral Research, 35*(1), 1–19.
- Murphy, K. R. (2019). Understanding how and why adding valid predictors can decrease the validity of selection composites: A generalization of Sackett, Dahlke, Shewach, and Kuncel (2017). *International Journal of Selection and Assessment, 27*(3), 249–255.
- Sackett, P. R., Dahlke, J. A., Shewach, O. R., & Kuncel, N. R. (2017). Effects of predictor weighting methods on incremental validity. *Journal of Applied Psychology, 102*(10), 1421–1434.
- Shapley, L. S. (1953). A value for n-person games. In H. W. Kuhn & A. W. Tucker (Eds.), *Contributions to the theory of games II* (pp. 307–317). Princeton University Press.
