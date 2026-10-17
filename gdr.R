# ==============================================================================
# gdr.R: General dominance analysis of predictive validity (GD-R)
# ==============================================================================
#
# Relative importance of correlated predictors by general dominance analysis
# (GDA). Traditional GDA partitions the model's explained variance (GD-R²);
# GD-R partitions its predictive validity, the multiple correlation R. Also
# provided: averaged semi-partial correlations (GD-sr), relative weights (RW;
# Johnson, 2000), incremental validity over baseline predictors, and bootstrap
# confidence intervals. Works from a correlation or covariance matrix; base R
# only.
#
# Main functions
#   calc_dominance()   GD values for one metric, from a formula and a matrix
#   analyze()          all metrics, standardized betas, and incremental validity
#   run_bootstrap()    bootstrap replicates of analyze()
#   summarize_boot()   percentile confidence intervals from those replicates
#
# Metric codes
#   "R2"      GD-R²: average incremental variance (ΔR²); values sum to R²
#   "R"       GD-R:  average incremental validity (ΔR); values sum to R
#   "spcorr"  GD-sr: average semi-partial correlation, sqrt(ΔR²); not additive
#   "rwa"     RW:    relative weights (analyze() only); values sum to R²
#
# User-specified predictor weights
#   By default, every subset of predictors is regression-weighted, so its
#   validity is the multiple R. Given fixed weights w, every subset is instead
#   scored as the composite sum(w_i * z_i) of its standardized predictors, with
#   validity R_c = w'r / sqrt(w'Rw), and the GD averaging is unchanged. The
#   weights are assumed to be the same in every subset; adding a predictor can
#   then lower R_c, so GD values can be negative. GD-sr and RW are defined only
#   for regression weights.
#
# Usage
#   source("gdr.R")
#   R <- matrix(c(1, .5, .6, .7,
#                 .5, 1, .4, .4,
#                 .6, .4, 1, .4,
#                 .7, .4, .4, 1), 4,
#               dimnames = list(c("Y", "X1", "X2", "X3"), c("Y", "X1", "X2", "X3")))
#   calc_dominance("Y ~ X1 + X2 + X3", R, objective = "R")
#
# Subsets of p predictors are stored as bitmasks 0, ..., 2^p - 1 (bit j set =
# predictor j included). The validity of every subset is computed once, and all
# metrics are derived from that vector, which keeps bootstrapping fast.
# ==============================================================================


# ---- Input checks -------------------------------------------------------------

#' Key under which analyze() stores each metric's total.
#' "<" cannot occur in variable names read with read.csv(check.names = TRUE).
TOTAL <- "<total>"

#' Split a formula into criterion and predictors.
#'
#' @param formula A formula or string such as "Y ~ X1 + X2".
#' @return A list with `y` (criterion name) and `xs` (predictor names).
parse_formula <- function(formula) {
  v <- all.vars(stats::as.formula(formula))
  if (length(v) < 2) {
    stop("The formula needs a criterion and at least one predictor, e.g. \"Y ~ X1 + X2\".")
  }
  list(y = v[1], xs = v[-1])
}

#' Stop unless every variable is a column of the matrix.
#'
#' @param M A matrix with column names.
#' @param vars Variable names to look for.
check_vars <- function(M, vars) {
  if (!all(vars %in% colnames(M))) {
    stop("Variables in formula not found in matrix columns.")
  }
}

#' Convert a covariance or correlation matrix to a named correlation matrix.
#'
#' Row names are copied from the column names when missing (as for matrices
#' read from a CSV file).
#'
#' @param M A square numeric matrix or data frame.
#' @return A correlation matrix with row and column names.
as_cor_matrix <- function(M) {
  M <- as.matrix(M)
  if (is.null(rownames(M))) rownames(M) <- colnames(M)
  tryCatch(
    stats::cov2cor(M),
    error = function(e) stop("Input must be a valid covariance/correlation matrix.")
  )
}

#' Check metric codes.
#'
#' Stops on unknown codes, and on RW or GD-sr when custom weights are given
#' (both are defined only for regression weights).
#'
#' @param metrics Metric codes.
#' @param w Custom weights, or NULL for regression weights.
#' @param allowed Codes accepted by the calling function.
#' @return `metrics`, invisibly.
check_metrics <- function(metrics, w = NULL, allowed = c("R2", "R", "spcorr", "rwa")) {
  bad <- setdiff(metrics, allowed)
  if (length(bad)) {
    stop(
      "Unknown metric: ", paste(bad, collapse = ", "),
      ". Use ", paste0("\"", allowed, "\"", collapse = ", "), "."
    )
  }
  if (!is.null(w) && any(metrics %in% c("rwa", "spcorr"))) {
    stop("RW and GDsr are regression-based; with custom weights use GD-R² (\"R2\") and GD-R (\"R\").")
  }
  invisible(metrics)
}

#' Check user-specified predictor weights.
#'
#' A named vector is matched to the predictors by name; an unnamed vector must
#' be in the order of `xs`. Weights must be finite and not all zero.
#'
#' @param w Numeric weights.
#' @param xs Predictor names.
#' @return The weights, named by predictor.
check_weights <- function(w, xs) {
  if (!is.null(names(w)) && all(xs %in% names(w))) w <- w[xs]
  if (!is.numeric(w) || length(w) != length(xs)) {
    stop("Supply one numeric weight per predictor (", length(xs), " predictors).")
  }
  if (any(!is.finite(w))) stop("Every weight must be a finite number.")
  if (all(w == 0)) stop("At least one weight must be non-zero.")
  setNames(as.numeric(w), xs)
}

#' Is a matrix a valid (positive definite) correlation matrix?
#'
#' @param M A symmetric matrix.
#' @param tol Smallest eigenvalue allowed.
#' @return TRUE if all entries are finite and all eigenvalues exceed `tol`.
is_pd <- function(M, tol = 1e-8) {
  all(is.finite(M)) &&
    min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > tol
}


# ---- Subsets and their validities -------------------------------------------

#' Bookkeeping for all subsets of p predictors.
#'
#' Subsets are numbered 0, ..., 2^p - 1, where bit j is set if predictor j is
#' included; position m + 1 of a validity vector belongs to subset m. For each
#' predictor i, every subset without i is paired with the same subset plus i,
#' and each pair gets the general-dominance weight 1 / (p * choose(p - 1, k)),
#' where k is the size of the subset without i. Summing weighted increments is
#' then equivalent to averaging within each k and then across k (Azen &
#' Budescu, 2003), i.e., the Shapley value. Depends only on p, so it can be
#' built once and reused.
#'
#' @param p Number of predictors.
#' @return A list with `p`, `members` (predictor indices in each subset), and
#'   `shap` (per predictor: positions of the subsets without (`wo`) and with
#'   (`w`) the predictor, and their weights `wt`).
make_subsets <- function(p) {
  masks <- 0:(2^p - 1)
  members <- lapply(masks, function(m) which(bitwAnd(m, 2^(0:(p - 1))) > 0))
  size <- lengths(members)
  shap <- lapply(seq_len(p), function(i) {
    bit <- 2^(i - 1)
    wo <- masks[bitwAnd(masks, bit) == 0]
    list(
      wo = wo + 1,
      w = wo + bit + 1,
      wt = 1 / (p * choose(p - 1, size[wo + 1]))
    )
  })
  list(p = p, members = members, shap = shap)
}

#' R² of the regression of the criterion on every subset of predictors.
#'
#' Computed from correlations alone, R² = r'R⁻¹r. Values are clipped to [0, 1]
#' to absorb rounding error; a subset whose matrix cannot be inverted gives NA.
#'
#' @param Rxx Predictor intercorrelation matrix (p x p).
#' @param rxy Predictor validities (length p).
#' @param S Subset bookkeeping from make_subsets(p).
#' @return R² for every subset, indexed by bitmask + 1 (the empty set is 0).
all_r2 <- function(Rxx, rxy, S) {
  r2 <- numeric(length(S$members))
  for (m in 2:length(r2)) {
    s <- S$members[[m]]
    v <- tryCatch(
      sum(rxy[s] * solve(Rxx[s, s, drop = FALSE], rxy[s])),
      error = function(e) NA_real_
    )
    r2[m] <- min(max(v, 0), 1)
  }
  r2
}

#' Validity of a fixed-weight composite of standardized predictors.
#'
#' @param Rxx Intercorrelations of the predictors in the composite.
#' @param rxy Their validities.
#' @param w Their weights.
#' @return w'r / sqrt(w'Rw); 0 if the composite has no variance (all weights 0).
composite_r <- function(Rxx, rxy, w) {
  v <- drop(crossprod(w, Rxx %*% w))
  if (v <= 0) 0 else sum(w * rxy) / sqrt(v)
}

#' Validity of every subset of predictors.
#'
#' With regression weights (`w = NULL`) this is the multiple R. With custom
#' weights it is the validity of the fixed-weight composite of each subset,
#' which can be negative.
#'
#' @inheritParams all_r2
#' @param w Custom weights, or NULL for regression weights.
#' @return Validity for every subset, indexed by bitmask + 1.
all_r <- function(Rxx, rxy, S, w = NULL) {
  if (is.null(w)) return(sqrt(all_r2(Rxx, rxy, S)))
  r <- numeric(length(S$members))
  for (m in 2:length(r)) {
    s <- S$members[[m]]
    r[m] <- composite_r(Rxx[s, s, drop = FALSE], rxy[s], w[s])
  }
  r
}

#' Validity of every subset, from a full correlation matrix and variable names.
#'
#' @param R Correlation matrix containing the criterion and the predictors.
#' @param y Criterion name.
#' @param xs Predictor names (their order defines the bitmasks).
#' @param S Subset bookkeeping from make_subsets(length(xs)).
#' @param w Custom weights, or NULL for regression weights.
#' @return As all_r().
subset_validities <- function(R, y, xs, S, w = NULL) {
  all_r(R[xs, xs, drop = FALSE], R[xs, y], S, w)
}


# ---- GD values from the subset validities ----------------------------------------

#' Increment rule for each GD metric.
#'
#' @param metric "R2", "R", or "spcorr".
#' @return A function of a (validity with the predictor) and b (validity
#'   without it): a² - b² (GD-R²), a - b (GD-R), or sqrt(max(a² - b², 0))
#'   (GD-sr).
increment_fn <- function(metric) {
  switch(metric,
    R2 = function(a, b) a^2 - b^2,
    R = function(a, b) a - b,
    spcorr = function(a, b) sqrt(pmax(a^2 - b^2, 0))
  )
}

#' General dominance values from subset validities.
#'
#' @param r Validity of every subset (from all_r() or subset_validities()).
#' @param S Subset bookkeeping from make_subsets(p).
#' @param metric "R2", "R", or "spcorr".
#' @return One GD value per predictor.
gd_from_r <- function(r, S, metric) {
  f <- increment_fn(metric)
  vapply(S$shap, function(z) sum(z$wt * f(r[z$w], r[z$wo])), numeric(1))
}

#' General dominance values from subset R² values (regression weights).
#'
#' @param r2 R² of every subset (from all_r2()).
#' @inheritParams gd_from_r
#' @return One GD value per predictor.
gd_from_r2 <- function(r2, S, metric) gd_from_r(sqrt(r2), S, metric)

#' Average increment of each predictor at each subset size k (conditional
#' dominance; Azen & Budescu, 2003).
#'
#' @inheritParams gd_from_r
#' @return A p x p matrix: row i is predictor i; column k + 1 is its average
#'   increment when added to a subset of k other predictors (k = 0 is the
#'   predictor on its own). Each row's mean is the predictor's GD value.
gd_levels <- function(r, S, metric) {
  f <- increment_fn(metric)
  size <- lengths(S$members)
  t(vapply(S$shap, function(z) {
    inc <- f(r[z$w], r[z$wo])
    k <- size[z$wo]
    vapply(0:(S$p - 1), function(kk) mean(inc[k == kk]), numeric(1))
  }, numeric(S$p)))
}

#' Percentage of each GD value contributed at each subset size k.
#'
#' @inheritParams gd_from_r
#' @return gd_levels() as row percentages (rows sum to 100). The first column
#'   is the share contributed at k = 0, i.e., by the predictor's standalone
#'   validity.
gd_level_shares <- function(r, S, metric) {
  lv <- gd_levels(r, S, metric)
  100 * lv / rowSums(lv)
}

#' Relative weights (Johnson, 2000).
#'
#' Regresses the criterion on the orthogonal counterparts of the predictors
#' (the symmetric square root of Rxx) and allocates the explained variance back
#' to the predictors.
#'
#' @param Rxx Predictor intercorrelation matrix.
#' @param rxy Predictor validities.
#' @return One relative weight per predictor; the weights sum to R².
rwa_weights <- function(Rxx, rxy) {
  e <- eigen(Rxx, symmetric = TRUE)
  Lam <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(Rxx)) %*% t(e$vectors)
  beta <- solve(Lam, rxy)
  as.numeric((Lam^2) %*% (beta^2))
}


# ---- Main functions ---------------------------------------------------------------

#' General dominance values for one metric.
#'
#' @param formula A formula or string such as "Y ~ X1 + X2 + X3".
#' @param cov_mat Covariance or correlation matrix containing those variables.
#' @param objective "R2" (GD-R²), "R" (GD-R), or "spcorr" (GD-sr).
#' @param weights Optional fixed predictor weights, in formula order or named
#'   by predictor. NULL (default) = regression weights.
#' @return A list with `table` (Variable, General_Dominance) and `total`, the
#'   sum of the values: the model R² for "R2" and the model R for "R" (with
#'   custom weights, the composite's R² and R).
#' @examples
#' calc_dominance("Y ~ X1 + X2 + X3", R, objective = "R")
#' calc_dominance("Y ~ X1 + X2 + X3", R, objective = "R", weights = c(1, 1, 1))
calc_dominance <- function(formula, cov_mat, objective = "R2", weights = NULL) {
  v <- parse_formula(formula)
  check_vars(cov_mat, c(v$y, v$xs))
  if (length(objective) != 1) stop("Choose one objective: \"R2\", \"R\", or \"spcorr\".")
  check_metrics(objective, weights, allowed = c("R2", "R", "spcorr"))
  if (!is.null(weights)) weights <- check_weights(weights, v$xs)
  R <- as_cor_matrix(cov_mat)

  S <- make_subsets(length(v$xs))
  r <- subset_validities(R, v$y, v$xs, S, weights)
  gen_dom <- setNames(gd_from_r(r, S, objective), v$xs)

  list(
    table = data.frame(Variable = names(gen_dom), General_Dominance = gen_dom),
    total = sum(gen_dom)
  )
}

#' All relative-importance results for one correlation matrix.
#'
#' Returns every statistic as one named vector, so that point estimates and
#' bootstrap replicates have the same structure. Keys:
#'   beta|X              standardized regression weight
#'   wt|X                custom weight (custom weights only)
#'   <m>|X, <m>|<total>  GD value (or RW) of predictor X, and the total
#'   <m>%|X              percentage of the total
#'   <m>#|X              rank (1 = most important)
#'   model|R, model|R2   validity of the full predictor set as scored
#'                       (regression-weighted, or the custom composite)
#'   reg|R, reg|R2       regression-weighted validity of the full set
#'   inc|<step>|<stat>   incremental validity over the baseline; <step> is
#'                       "base", "add:X" (baseline + X), or "all" (all
#'                       predictors); <stat> is R, R2, dR, or dR2
#'
#' @param R Correlation matrix containing the criterion and the predictors.
#' @param y Criterion name.
#' @param xs Predictor names.
#' @param S Subset bookkeeping from make_subsets(length(xs)).
#' @param metrics Metric codes: any of "R2", "R", "spcorr", "rwa".
#' @param w Custom predictor weights, or NULL for regression weights. With
#'   custom weights, only "R2" and "R" are allowed, and every model (including
#'   the incremental-validity steps) is scored as the weighted composite.
#' @param baseline Baseline predictors for incremental validity (optional).
#' @return A named numeric vector.
analyze <- function(R, y, xs, S, metrics, w = NULL, baseline = character(0)) {
  check_metrics(metrics, w)
  if (!is.null(w)) w <- check_weights(w, xs)
  Rxx <- R[xs, xs, drop = FALSE]
  rxy <- R[xs, y]
  r <- subset_validities(R, y, xs, S, w)
  r2 <- r^2
  out <- numeric(0)

  out[paste0("beta|", xs)] <- solve(Rxx, rxy)
  if (!is.null(w)) out[paste0("wt|", xs)] <- w

  for (m in metrics) {
    v <- if (m == "rwa") rwa_weights(Rxx, rxy) else gd_from_r(r, S, m)
    tot <- sum(v)
    out[paste0(m, "|", xs)] <- v
    out[paste0(m, "|", TOTAL)] <- tot
    out[paste0(m, "%|", xs)] <- 100 * v / tot
    out[paste0(m, "#|", xs)] <- rank(-v)
  }

  full <- length(r)
  out["model|R2"] <- r2[full]
  out["model|R"] <- r[full]
  out["reg|R2"] <- if (is.null(w)) r2[full] else sum(rxy * solve(Rxx, rxy))
  out["reg|R"] <- sqrt(out[["reg|R2"]])

  if (length(baseline)) {
    mask <- function(v) sum(2^(match(v, xs) - 1))
    rb <- r[mask(baseline) + 1]
    step <- function(key, val) {
      out[paste0("inc|", key, "|R2")] <<- val^2
      out[paste0("inc|", key, "|R")] <<- val
      out[paste0("inc|", key, "|dR2")] <<- val^2 - rb^2
      out[paste0("inc|", key, "|dR")] <<- val - rb
    }
    out["inc|base|R2"] <- rb^2
    out["inc|base|R"] <- rb
    added <- setdiff(xs, baseline)
    for (j in added) step(paste0("add:", j), r[mask(c(baseline, j)) + 1])
    if (length(added) > 1) step("all", r[full])
  }

  out
}

#' Bootstrap replicates of analyze().
#'
#' Raw data: cases are resampled with replacement (nonparametric bootstrap).
#' Correlation matrix: replicate matrices are drawn from the Wishart
#' distribution with N - 1 degrees of freedom, which is the sampling
#' distribution of the covariance matrix of N multivariate-normal cases
#' (parametric bootstrap). Replicates that are not positive definite are
#' dropped. Custom weights, if given, are used in every replicate.
#'
#' @param type "data" or "matrix".
#' @param R Correlation matrix (used when type = "matrix").
#' @param X Raw data matrix (used when type = "data").
#' @param N Sample size (used when type = "matrix").
#' @param B Number of replicates.
#' @param seed Random seed.
#' @param y,xs,S,metrics,w,baseline As in analyze().
#' @param keys Names of the statistics to keep (e.g., names of the point
#'   estimates from analyze()).
#' @param tick Optional progress callback, called about 100 times with the
#'   fraction completed since the last call.
#' @return A matrix with one row per usable replicate and one column per key.
run_bootstrap <- function(type, R, X, N, B, seed, y, xs, S, metrics, w, baseline,
                          keys, tick = NULL) {
  set.seed(seed)
  vars <- c(y, xs)
  if (type == "matrix") {
    W <- stats::rWishart(B, N - 1, R[vars, vars])
  }
  out <- matrix(NA_real_, B, length(keys), dimnames = list(NULL, keys))
  every <- max(1L, B %/% 100L)

  for (b in seq_len(B)) {
    Rb <- if (type == "data") {
      suppressWarnings(stats::cor(X[sample.int(nrow(X), replace = TRUE), vars, drop = FALSE]))
    } else {
      m <- stats::cov2cor(W[, , b])
      dimnames(m) <- list(vars, vars)
      m
    }
    if (is_pd(Rb)) {
      out[b, ] <- analyze(Rb, y, xs, S, metrics, w, baseline)[keys]
    }
    if (!is.null(tick) && b %% every == 0) tick(every / B)
  }
  out[stats::complete.cases(out), , drop = FALSE]
}

#' Percentile confidence intervals from bootstrap replicates.
#'
#' Ranks use quantile type 1, so their limits are attainable ranks.
#'
#' @param est Point estimates (from analyze()).
#' @param boot Replicates (from run_bootstrap()).
#' @param conf Confidence level, e.g., .95.
#' @return A list with `est`, `lo`, and `hi`, each named like `est`.
summarize_boot <- function(est, boot, conf) {
  a <- (1 - conf) / 2
  lo <- hi <- setNames(rep(NA_real_, length(est)), names(est))
  for (k in names(est)) {
    q <- stats::quantile(
      boot[, k], c(a, 1 - a),
      names = FALSE, na.rm = TRUE,
      type = if (grepl("#|", k, fixed = TRUE)) 1 else 7
    )
    lo[k] <- q[1]
    hi[k] <- q[2]
  }
  list(est = est, lo = lo, hi = hi)
}
