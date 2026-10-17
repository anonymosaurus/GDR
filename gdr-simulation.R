# =============================================================================
# gdr-simulation.R: GD-R Monte Carlo simulation (Study 1)
# =============================================================================
#
# Compares four relative-importance methods in samples drawn from known
# populations:
#   GDR2  general dominance of R² (traditional GDA)
#   GDR   general dominance of R (GD-R)
#   GDsr  general dominance of semi-partial correlations (GD-sr)
#   RW    relative weights (Johnson, 2000)
#
# Design: sample size (N) x number of predictors (p) x mean predictor
# intercorrelation (rbar) x mean validity (vbar) x predictor reliability (rel),
# fully crossed (288 conditions), with PARAM$reps replications per condition.
#
# Pipeline (see "Main" at the end of the file)
#   1. build_populations()          error-free and attenuated population matrices
#   2. run_study()                  samples and the four methods (SimDesign)
#   3. condition_table(), importance_table(), k0_table()   results as data frames
#   4. write_importance_workbook()  manuscript-style tables, one sheet per p
#   5. study1_tables()              Study 1 tables (Tables 3, 4, B3, E1, E2)
#
# Requirements: R packages SimDesign and openxlsx
#   (install.packages(c("SimDesign", "openxlsx"))) and an internet connection:
#   the GD engine (gdr.R) is sourced from the GitHub repository.
# Run from this folder:  Rscript gdr-simulation.R
#
# Output (in OUTPUT_DIR):
#   study1-results.rds   parameters, population matrices, condition results
#   study1-results.xlsx  manuscript-style tables and k = 0 shares
#   study1-tables.xlsx   Study 1 tables
# =============================================================================

library(SimDesign)
library(openxlsx)

GDR_URL <- "https://raw.githubusercontent.com/anonymosaurus/GDR/refs/heads/main/gdr.R"
source(GDR_URL)


# ---- Settings ----------------------------------------------------------------

# Simulation factors, the within-matrix SDs of the intercorrelations and
# validities, replications per condition, and the random seed.
PARAM <- list(
  N = c(100, 200, 500, 1000),
  p = c(3, 5, 7, 9),
  rbar = c(.10, .30, .50),
  vbar = c(.10, .20, .30),
  rel = c(.90, .70),
  r_sd = .10,
  v_sd = .20,
  reps = 1000,
  seed = 2026
)

# TRUE: run the simulation and save the results.
# FALSE: rebuild the output from the results saved
# by an earlier run with the same PARAM (fast).
RUN_SIMULATION <- TRUE

# Folder and names of the output files.
OUTPUT_DIR <- "sim_results"
RESULTS_RDS <- file.path(OUTPUT_DIR, "study1-results.rds")
RESULTS_XLSX <- file.path(OUTPUT_DIR, "study1-results.xlsx")
TABLES_XLSX <- file.path(OUTPUT_DIR, "study1-tables.xlsx")

# Method codes, in the order used throughout.
METHODS <- c("GDR2", "RW", "GDR", "GDsr")


# =============================================================================
# 1. Population correlation matrices
# =============================================================================

#' k values with mean exactly m and SD exactly s (divisor k).
#'
#' Evenly spaced quantiles of the normal distribution, standardized and
#' rescaled.
#'
#' @param k Number of values.
#' @param m,s Target mean and standard deviation.
#' @return A numeric vector of length k.
exact_values <- function(k, m, s) {
  z <- qnorm(ppoints(k))
  m + s * (z - mean(z)) / sqrt(mean((z - mean(z))^2))
}

#' One error-free population correlation matrix.
#'
#' The p(p - 1)/2 intercorrelations (mean rbar, SD r_sd) are assigned to
#' predictor pairs at random; the p validities have mean vbar and SD v_sd.
#' Arrangements that are not positive definite are redrawn. The multiple R
#' is not constrained.
#'
#' @param p Number of predictors.
#' @param rbar,vbar Mean intercorrelation and mean validity.
#' @param r_sd,v_sd SDs of the intercorrelations and validities.
#' @param max_tries Number of arrangements to try before giving up.
#' @return A (p + 1) x (p + 1) correlation matrix, criterion "Y" first.
make_population <- function(p, rbar, vbar, r_sd, v_sd, max_tries = 1e5) {
  intercor <- exact_values(p * (p - 1) / 2, rbar, r_sd)
  validity <- exact_values(p, vbar, v_sd)
  for (attempt in seq_len(max_tries)) {
    Rxx <- diag(p)
    Rxx[lower.tri(Rxx)] <- sample(intercor)
    Rxx <- Rxx + t(Rxx) - diag(p)
    M <- rbind(c(1, validity), cbind(validity, Rxx))
    if (is_pd(M)) {
      dimnames(M) <- list(c("Y", paste0("X", 1:p)), c("Y", paste0("X", 1:p)))
      return(M)
    }
  }
  stop("No valid population for p = ", p, ", rbar = ", rbar, ", vbar = ", vbar)
}

#' Attenuate a population matrix for predictor unreliability.
#'
#' Every predictor has reliability rel and the criterion is error-free, so
#' intercorrelations are multiplied by rel and validities by sqrt(rel). The
#' predictors are then renumbered by their GD-R² values in the attenuated
#' matrix, so X1 is the most important predictor and Xp the least.
#'
#' @param M Error-free population matrix (criterion first).
#' @param rel Predictor reliability.
#' @param S Subset bookkeeping from make_subsets(p).
#' @return The attenuated, renumbered matrix.
attenuate <- function(M, rel, S) {
  d <- sqrt(c(1, rep(rel, ncol(M) - 1)))
  M <- M * outer(d, d)
  diag(M) <- 1
  gdr2 <- gd_from_r2(all_r2(M[-1, -1], M[-1, 1], S), S, "R2")
  keep <- c(1, 1 + order(gdr2, decreasing = TRUE))
  M <- M[keep, keep]
  labels <- c("Y", paste0("X", 1:(ncol(M) - 1)))
  dimnames(M) <- list(labels, labels)
  M
}

#' Name of a population matrix, e.g. "p5_r0.30_v0.20" or "p5_r0.30_v0.20_rel0.90".
#'
#' @param p,rbar,vbar,rel Design values (rel omitted for error-free matrices).
#' @return A character vector of names.
population_key <- function(p, rbar, vbar, rel = NULL) {
  key <- sprintf("p%d_r%.2f_v%.2f", p, rbar, vbar)
  if (is.null(rel)) key else sprintf("%s_rel%.2f", key, rel)
}

#' All population matrices of the design.
#'
#' One error-free matrix per p x rbar x vbar cell (generated after setting the
#' seed), shared by all sample sizes and both reliability levels, and its
#' attenuated, renumbered version at each reliability level.
#'
#' @param design The SimDesign design table.
#' @param param Simulation parameters (r_sd, v_sd, seed).
#' @param subsets Subset bookkeeping for each p, named by p.
#' @return A list with `error_free` and `observed`, each a named list of
#'   matrices.
build_populations <- function(design, param, subsets) {
  cells <- unique(as.data.frame(design)[, c("p", "rbar", "vbar")])
  set.seed(param$seed)
  error_free <- lapply(seq_len(nrow(cells)), function(i) {
    make_population(
      cells$p[i],
      cells$rbar[i],
      cells$vbar[i],
      param$r_sd,
      param$v_sd
    )
  })
  names(error_free) <- population_key(cells$p, cells$rbar, cells$vbar)

  obs <- unique(as.data.frame(design)[, c("p", "rbar", "vbar", "rel")])
  observed <- lapply(seq_len(nrow(obs)), function(i) {
    attenuate(
      error_free[[population_key(obs$p[i], obs$rbar[i], obs$vbar[i])]],
      obs$rel[i],
      subsets[[as.character(obs$p[i])]]
    )
  })
  names(observed) <- population_key(obs$p, obs$rbar, obs$vbar, obs$rel)

  list(error_free = error_free, observed = observed)
}

#' GD-R², RW, GD-R, and GD-sr values and the full-model R² of one matrix.
#'
#' @param M Correlation matrix with the criterion in row/column 1.
#' @param S Subset bookkeeping from make_subsets(p).
#' @return A list with one vector per method and the scalar R2.
importance <- function(M, S) {
  Rxx <- M[-1, -1]
  rxy <- M[-1, 1]
  r2 <- all_r2(Rxx, rxy, S)
  list(
    GDR2 = gd_from_r2(r2, S, "R2"),
    RW = rwa_weights(Rxx, rxy),
    GDR = gd_from_r2(r2, S, "R"),
    GDsr = gd_from_r2(r2, S, "spcorr"),
    R2 = r2[length(r2)]
  )
}


# =============================================================================
# 2. Simulation (SimDesign)
# =============================================================================
# SimDesign calls Generate, Analyze, and Summarize for every condition;
# `fixed_objects` carries the population matrices, the subset bookkeeping,
# the method codes, and the largest p.

#' The attenuated population matrix that a condition's samples come from.
condition_population <- function(condition, fixed_objects) {
  key <- population_key(
    condition$p,
    condition$rbar,
    condition$vbar,
    condition$rel
  )
  fixed_objects$observed[[key]]
}

#' Generate: one sample correlation matrix of size N.
#'
#' The sample covariance matrix of N multivariate-normal cases is drawn
#' directly from the Wishart distribution with N - 1 df and rescaled to
#' correlations.
Generate <- function(condition, fixed_objects) {
  Sigma <- condition_population(condition, fixed_objects)
  stats::cov2cor(stats::rWishart(1, condition$N - 1, Sigma)[,, 1])
}

#' Analyze: every predictor's value under each method, plus the sample R².
#'
#' @return A named vector: GDR2_1, ..., GDsr_p, R2.
Analyze <- function(condition, dat, fixed_objects) {
  w <- importance(dat, fixed_objects$S[[as.character(condition$p)]])
  c(
    unlist(lapply(fixed_objects$methods, function(m) {
      setNames(w[[m]], paste0(m, "_", 1:condition$p))
    })),
    R2 = w$R2
  )
}

#' Summarize: one row of results per condition.
#'
#' For each method and predictor: the population value, the mean and SD over
#' replications, the population and mean sample % of total (sample % computed
#' in each replication), the population and mean sample rank, and the
#' population k = 0 share (NA for RW). For each method: population and mean
#' sample totals, and the proportions of replications that identify the most
#' important predictor and that reproduce the complete population ranking.
#' Predictor columns beyond p are NA, so every condition has the same columns.
Summarize <- function(condition, results, fixed_objects) {
  p <- condition$p
  p_max <- fixed_objects$p_max
  S <- fixed_objects$S[[as.character(p)]]
  M <- condition_population(condition, fixed_objects)
  pop <- importance(M, S)
  pad <- function(x) c(x, rep(NA, p_max - p))
  r_pop <- all_r(M[-1, -1], M[-1, 1], S)
  k0 <- list(
    GDR2 = gd_level_shares(r_pop, S, "R2")[, 1],
    GDR = gd_level_shares(r_pop, S, "R")[, 1],
    GDsr = gd_level_shares(r_pop, S, "spcorr")[, 1],
    RW = rep(NA_real_, p)
  )
  out <- c(
    pop_R2 = pop$R2,
    mean_R2 = mean(results[, "R2"]),
    sd_R2 = sd(results[, "R2"])
  )
  for (m in fixed_objects$methods) {
    E <- as.matrix(results[, paste0(m, "_", 1:p)])
    ranks <- t(apply(-E, 1, rank))
    stats <- list(
      pop = pop[[m]],
      mean = colMeans(E),
      sd = apply(E, 2, sd),
      poppct = 100 * pop[[m]] / sum(pop[[m]]),
      meanpct = 100 * colMeans(E / rowSums(E)),
      poprank = rank(-pop[[m]]),
      meanrank = colMeans(ranks),
      popk0 = k0[[m]]
    )
    for (s in names(stats)) {
      out[paste0(m, "_", s, "_", 1:p_max)] <- pad(stats[[s]])
    }
    out[paste0(m, "_pop_total")] <- sum(pop[[m]])
    out[paste0(m, "_mean_total")] <- mean(rowSums(E))
    out[paste0(m, "_top")] <- mean(
      max.col(E, ties.method = "first") == which.max(pop[[m]])
    )
    out[paste0(m, "_full")] <- mean(apply(ranks, 1, function(r) {
      all(r == rank(-pop[[m]]))
    }))
  }
  out
}

#' Run the simulation for every condition of the design.
#'
#' Each condition gets its own fixed seed (seed + row number), so every
#' result can be reproduced exactly.
#'
#' @param design The SimDesign design table.
#' @param param Simulation parameters (reps, p, seed).
#' @param observed Attenuated population matrices (from build_populations()).
#' @param subsets Subset bookkeeping for each p.
#' @return A data frame with one row per condition (see Summarize).
run_study <- function(design, param, observed, subsets) {
  res <- runSimulation(
    design = design,
    replications = param$reps,
    generate = Generate,
    analyse = Analyze,
    summarise = Summarize,
    fixed_objects = list(
      observed = observed,
      S = subsets,
      methods = METHODS,
      p_max = max(param$p)
    ),
    seed = param$seed + seq_len(nrow(design)),
    save = FALSE,
    progress = FALSE,
    verbose = FALSE
  )
  as.data.frame(res)
}


# =============================================================================
# 3. Results as data frames
# =============================================================================

COND_COLS <- c("N", "p", "rbar", "vbar", "rel")

#' One row per condition: population and mean sample R², method totals, and
#' top-predictor and complete-ranking recovery.
#'
#' @param res Simulation results (from run_study()).
#' @return A data frame.
condition_table <- function(res) {
  conditions <- res[, COND_COLS]
  conditions$pop_R2 <- res$pop_R2
  conditions$mean_R2 <- res$mean_R2
  conditions$sd_R2 <- res$sd_R2
  for (m in METHODS) {
    for (s in c("pop_total", "mean_total", "top", "full")) {
      conditions[[paste0(m, "_", s)]] <- res[[paste0(m, "_", s)]]
    }
  }
  conditions
}

#' One row per condition x method x predictor.
#'
#' @param res Simulation results (from run_study()).
#' @return A data frame with population and mean sample values, % of total,
#'   ranks, and the population k = 0 share.
importance_table <- function(res) {
  do.call(
    rbind,
    lapply(seq_len(nrow(res)), function(i) {
      do.call(
        rbind,
        lapply(METHODS, function(m) {
          j <- 1:res$p[i]
          g <- function(s) unlist(res[i, paste0(m, "_", s, "_", j)])
          data.frame(
            res[i, COND_COLS],
            method = m,
            predictor = paste0("X", j),
            pop_weight = g("pop"),
            mean_weight = g("mean"),
            sd_weight = g("sd"),
            pop_pct = g("poppct"),
            mean_pct = g("meanpct"),
            pop_rank = g("poprank"),
            mean_rank = g("meanrank"),
            pop_k0_pct = g("popk0"),
            row.names = NULL
          )
        })
      )
    })
  )
}

#' k = 0 shares: one row per predictor and population matrix, one column per
#' GD method (k0_GDR2, k0_GDR, k0_GDsr).
#'
#' Population values do not depend on N, so the rows of one N are used.
#'
#' @param importance_long Output of importance_table().
#' @param n Sample size whose rows are used.
#' @return A data frame.
k0_table <- function(importance_long, n) {
  k0 <- importance_long[
    importance_long$N == n &
      importance_long$method %in% c("GDR2", "GDR", "GDsr"),
    c("p", "rbar", "vbar", "rel", "method", "predictor", "pop_k0_pct")
  ]
  k0 <- reshape(
    k0,
    idvar = c("p", "rbar", "vbar", "rel", "predictor"),
    timevar = "method",
    direction = "wide",
    sep = "_"
  )
  names(k0) <- sub("pop_k0_pct_", "k0_", names(k0))
  k0
}


# =============================================================================
# 4. Excel workbook of manuscript-style tables
# =============================================================================
# One sheet per number of predictors. Each block is one population
# (mean r x mean validity x reliability). For each method, three rows give
# the weight, the % of total, and the rank: for the population (where % of
# total and rank are Excel formulas) and for the sample averages at each N,
# followed by the % of samples that reproduce the complete population
# ranking. Rank cells are highlighted with conditional formatting:
#   population group: orange = rank differs from GD-R²'s; red = the most
#     important predictor differs from GD-R²'s
#   sample groups: blue = the order by average rank differs from the
#     method's own population order
# A last sheet gives the k = 0 shares.

LABEL <- c(GDR2 = "GD-R²", RW = "RW", GDR = "GD-R", GDsr = "GDsr")
FONT <- "Arial"

#' Format numbers APA style, without the leading zero (e.g., ".30").
apa_num <- function(x, d = 2) {
  sub("^0\\.", ".", formatC(x, format = "f", digits = d))
}

#' Rows of a data frame in predictor order (X1, X2, ..., X10).
order_predictors <- function(d) {
  d[order(as.integer(sub("X", "", d$predictor))), ]
}

#' Write a row of Excel formulas (given without the leading "=").
write_formula_row <- function(wb, sheet, formulas, row, col) {
  df <- as.data.frame(as.list(formulas), stringsAsFactors = FALSE)
  for (j in seq_along(df)) {
    class(df[[j]]) <- c(class(df[[j]]), "formula")
  }
  writeData(wb, sheet, df, startCol = col, startRow = row, colNames = FALSE)
}

#' Cell styles used in the workbook (Arial 10 unless noted).
workbook_styles <- function() {
  centered <- function(...) {
    createStyle(
      fontName = FONT,
      fontSize = 10,
      halign = "center",
      valign = "center",
      ...
    )
  }
  top_line <- list(border = "top", borderColour = "#BFBFBF")
  list(
    title = createStyle(
      fontName = FONT,
      fontSize = 12,
      textDecoration = "bold",
      halign = "left"
    ),
    note = createStyle(
      fontName = FONT,
      fontSize = 9,
      textDecoration = "italic",
      fontColour = "#595959",
      halign = "left"
    ),
    block = createStyle(
      fontName = FONT,
      fontSize = 10,
      textDecoration = "bold",
      fgFill = "#D9E1F2",
      halign = "left",
      valign = "center"
    ),
    head = centered(textDecoration = "bold", fgFill = "#F2F2F2"),
    method = do.call(centered, c(list(textDecoration = "bold"), top_line)),
    stat_top = do.call(
      createStyle,
      c(
        list(
          fontName = FONT,
          fontSize = 10,
          halign = "left",
          valign = "center"
        ),
        top_line
      )
    ),
    stat = createStyle(
      fontName = FONT,
      fontSize = 10,
      halign = "left",
      valign = "center"
    ),
    w_top = do.call(centered, c(list(numFmt = "0.0000"), top_line)),
    wtot_top = do.call(
      centered,
      c(list(numFmt = "0.0000", textDecoration = "bold"), top_line)
    ),
    gap_top = do.call(centered, top_line),
    rec_top = do.call(centered, c(list(numFmt = "0%"), top_line)),
    k0 = centered(numFmt = "0.0%"),
    k0_mean = centered(
      numFmt = "0.0%",
      textDecoration = "bold",
      border = "top",
      borderColour = "#BFBFBF"
    ),
    cell = centered(),
    pct = centered(numFmt = "0.0%"),
    rank0 = centered(numFmt = "0"),
    rank1 = centered(numFmt = "0.0"),
    orange = createStyle(bgFill = "#FCE4D6"),
    red = createStyle(
      bgFill = "#F8CBAD",
      fontColour = "#9C0006",
      textDecoration = "bold"
    ),
    blue = createStyle(bgFill = "#BDD7EE")
  )
}

#' Collects cells by style while a sheet is written, then styles each group
#' of cells in one call (much faster than styling cell by cell).
#'
#' @return A list of two functions: mark(style, rows, cols) and
#'   apply(wb, sheet, styles).
style_collector <- function() {
  styled <- list()
  list(
    mark = function(style, rows, cols) {
      styled[[style]] <<- rbind(
        styled[[style]],
        expand.grid(r = rows, c = cols)
      )
    },
    apply = function(wb, sheet, styles) {
      for (s in names(styled)) {
        addStyle(
          wb,
          sheet,
          styles[[s]],
          rows = styled[[s]]$r,
          cols = styled[[s]]$c,
          gridExpand = FALSE
        )
      }
    }
  )
}

#' Column positions of an importance sheet with p predictors.
#'
#' Columns: Method | Statistic | one group (X1..Xp, Total) for the population
#' and for each N, separated by a gap column | the recovery panel (one column
#' per N).
#'
#' @param p Number of predictors.
#' @param n_values Sample sizes.
#' @return A list of column numbers.
sheet_layout <- function(p, n_values) {
  groups <- c("pop", n_values)
  width <- p + 1
  start <- setNames(3 + (seq_along(groups) - 1) * (width + 1), groups)
  rec <- unname(start[length(groups)]) + width + 1
  data_cols <- unlist(lapply(start, function(c0) c0:(c0 + p - 1)))
  pop_cols <- start["pop"]:(start["pop"] + p - 1)
  list(
    p = p,
    groups = groups,
    width = width,
    start = start,
    rec = rec,
    rec_cols = rec:(rec + length(n_values) - 1),
    last = rec + length(n_values) - 1,
    data_cols = data_cols,
    total_cols = unname(start) + p,
    gap_cols = c(unname(start[-1]) - 1, rec - 1),
    pop_cols = pop_cols,
    samp_cols = setdiff(data_cols, pop_cols)
  )
}

#' Write the title row and the two header rows of one population block.
write_block_header <- function(wb, sheet, r, blk, cond, lay, n_values, mark) {
  title <- paste0(
    "Mean r = %s  ·  Mean validity = %s  ·  Reliability = %s  ·  ",
    "Population R² = %s (R = %s)  ·  Mean sample R² at N = %s: %s"
  )
  writeData(
    wb,
    sheet,
    sprintf(
      title,
      apa_num(blk$rbar),
      apa_num(blk$vbar),
      apa_num(blk$rel),
      apa_num(cond$pop_R2[1], 3),
      apa_num(sqrt(cond$pop_R2[1]), 3),
      paste(format(n_values, big.mark = ",", trim = TRUE), collapse = " / "),
      paste(apa_num(cond$mean_R2, 3), collapse = " / ")
    ),
    startCol = 1,
    startRow = r
  )
  mark("block", r, 1:lay$last)

  writeData(
    wb,
    sheet,
    t(c("Method", "Statistic")),
    startCol = 1,
    startRow = r + 2,
    colNames = FALSE
  )
  for (g in lay$groups) {
    c0 <- unname(lay$start[g])
    mergeCells(wb, sheet, cols = c0:(c0 + lay$width - 1), rows = r + 1)
    group_title <- if (g == "pop") {
      "Population"
    } else {
      sprintf("Sample average, N = %s", format(as.numeric(g), big.mark = ","))
    }
    writeData(wb, sheet, group_title, startCol = c0, startRow = r + 1)
    writeData(
      wb,
      sheet,
      t(c(paste0("X", 1:lay$p), "Total")),
      startCol = c0,
      startRow = r + 2,
      colNames = FALSE
    )
  }
  mergeCells(wb, sheet, cols = lay$rec_cols, rows = r + 1)
  writeData(
    wb,
    sheet,
    "Full rank order recovered (% of samples)",
    startCol = lay$rec,
    startRow = r + 1
  )
  writeData(
    wb,
    sheet,
    t(sprintf("N = %s", format(n_values, big.mark = ",", trim = TRUE))),
    startCol = lay$rec,
    startRow = r + 2,
    colNames = FALSE
  )
  mark("head", (r + 1):(r + 2), 1:lay$last)
}

#' Write the three rows (weight, % of total, rank) of one method in a block.
#'
#' Population % of total and ranks are Excel formulas; sample groups hold
#' averages over the replications. Totals are SUM formulas.
write_method_rows <- function(
  wb,
  sheet,
  w_row,
  m,
  imp,
  cond,
  lay,
  n_values,
  mark
) {
  pct_row <- w_row + 1
  rank_row <- w_row + 2
  mergeCells(wb, sheet, cols = 1, rows = w_row:rank_row)
  writeData(wb, sheet, unname(LABEL[m]), startCol = 1, startRow = w_row)
  writeData(
    wb,
    sheet,
    c("Weight", "% of total", "Rank"),
    startCol = 2,
    startRow = w_row
  )

  for (g in lay$groups) {
    c0 <- unname(lay$start[g])
    first <- int2col(c0)
    last_x <- int2col(c0 + lay$p - 1)
    total <- int2col(c0 + lay$p)
    n <- if (g == "pop") n_values[1] else as.numeric(g)
    x <- order_predictors(imp[imp$method == m & imp$N == n, ])
    if (g == "pop") {
      writeData(
        wb,
        sheet,
        t(x$pop_weight),
        startCol = c0,
        startRow = w_row,
        colNames = FALSE
      )
      cells <- int2col(c0:(c0 + lay$p - 1))
      write_formula_row(
        wb,
        sheet,
        sprintf("%s%d/$%s$%d", cells, w_row, total, w_row),
        pct_row,
        c0
      )
      write_formula_row(
        wb,
        sheet,
        sprintf(
          "RANK(%s%d,$%s$%d:$%s$%d,0)",
          cells,
          w_row,
          first,
          w_row,
          last_x,
          w_row
        ),
        rank_row,
        c0
      )
    } else {
      writeData(
        wb,
        sheet,
        t(x$mean_weight),
        startCol = c0,
        startRow = w_row,
        colNames = FALSE
      )
      writeData(
        wb,
        sheet,
        t(x$mean_pct / 100),
        startCol = c0,
        startRow = pct_row,
        colNames = FALSE
      )
      writeData(
        wb,
        sheet,
        t(x$mean_rank),
        startCol = c0,
        startRow = rank_row,
        colNames = FALSE
      )
    }
    write_formula_row(
      wb,
      sheet,
      sprintf("SUM(%s%d:%s%d)", first, w_row, last_x, w_row),
      w_row,
      c0 + lay$p
    )
    write_formula_row(
      wb,
      sheet,
      sprintf("SUM(%s%d:%s%d)", first, pct_row, last_x, pct_row),
      pct_row,
      c0 + lay$p
    )
  }

  for (j in seq_along(n_values)) {
    mergeCells(wb, sheet, cols = lay$rec + j - 1, rows = w_row:rank_row)
  }
  writeData(
    wb,
    sheet,
    t(cond[[paste0(m, "_full")]]),
    startCol = lay$rec,
    startRow = w_row,
    colNames = FALSE
  )
  mark("rec_top", w_row, lay$rec_cols)

  mark("method", w_row, 1)
  mark("stat_top", w_row, 2)
  mark("stat", c(pct_row, rank_row), 2)
  mark("w_top", w_row, lay$data_cols)
  mark("wtot_top", w_row, lay$total_cols)
  mark("gap_top", w_row, lay$gap_cols)
  mark("pct", pct_row, c(lay$data_cols, lay$total_cols))
  mark("rank0", rank_row, lay$pop_cols)
  mark("rank1", rank_row, lay$samp_cols)
}

#' Conditional formatting on the rank rows of one block (see the section
#' description for the colors).
#'
#' @param rank_rows Row number of each method's rank row, named by method.
add_rank_highlighting <- function(wb, sheet, rank_rows, lay, styles) {
  c0 <- unname(lay$start["pop"])
  pop_first <- int2col(c0)
  pop_cols <- c0:(c0 + lay$p - 1)
  ref_rank <- sprintf("%s%d", pop_first, rank_rows["GDR2"])
  for (m in c("RW", "GDR", "GDsr")) {
    this_rank <- sprintf("%s%d", pop_first, rank_rows[m])
    top_differs <- sprintf("AND(%s=1,%s<>1)", this_rank, ref_rank)
    rank_differs <- sprintf(
      "AND(%s<>%s,NOT(%s))",
      this_rank,
      ref_rank,
      top_differs
    )
    conditionalFormatting(
      wb,
      sheet,
      cols = pop_cols,
      rows = rank_rows[m],
      rule = top_differs,
      style = styles$red,
      type = "expression"
    )
    conditionalFormatting(
      wb,
      sheet,
      cols = pop_cols,
      rows = rank_rows[m],
      rule = rank_differs,
      style = styles$orange,
      type = "expression"
    )
  }
  for (g in lay$groups[-1]) {
    c0 <- unname(lay$start[g])
    first <- int2col(c0)
    last_x <- int2col(c0 + lay$p - 1)
    for (m in METHODS) {
      rk <- rank_rows[m]
      sample_position <- sprintf(
        "RANK(%s%d,$%s$%d:$%s$%d,1)",
        first,
        rk,
        first,
        rk,
        last_x,
        rk
      )
      population_rank <- sprintf("%s$%d", pop_first, rk)
      conditionalFormatting(
        wb,
        sheet,
        cols = c0:(c0 + lay$p - 1),
        rows = rk,
        rule = sprintf("%s<>%s", sample_position, population_rank),
        style = styles$blue,
        type = "expression"
      )
    }
  }
}

#' Write one importance sheet: all population blocks for p predictors.
write_importance_sheet <- function(
  wb,
  p,
  importance_long,
  conditions,
  n_values,
  styles
) {
  sheet <- sprintf("%d predictors", p)
  addWorksheet(wb, sheet, gridLines = FALSE)
  lay <- sheet_layout(p, n_values)
  styler <- style_collector()

  writeData(
    wb,
    sheet,
    sprintf(
      paste0(
        "Predictor importance with %d predictors: ",
        "population values and sample averages, by method"
      ),
      p
    ),
    startCol = 1,
    startRow = 1
  )
  writeData(
    wb,
    sheet,
    paste(
      paste0(
        "Population group: orange = rank differs from GD-R²'s; ",
        "red = most important predictor differs from GD-R²'s."
      ),
      " Sample groups: blue = order by average rank differs from the method's population order."
    ),
    startCol = 1,
    startRow = 2
  )
  addStyle(wb, sheet, styles$title, rows = 1, cols = 1)
  addStyle(wb, sheet, styles$note, rows = 2, cols = 1)

  blocks <- unique(importance_long[
    importance_long$p == p,
    c("rbar", "vbar", "rel")
  ])
  blocks <- blocks[order(blocks$rbar, blocks$vbar, -blocks$rel), ]

  r <- 4
  for (i in seq_len(nrow(blocks))) {
    blk <- blocks[i, ]
    in_block <- function(d) {
      d$p == p & d$rbar == blk$rbar & d$vbar == blk$vbar & d$rel == blk$rel
    }
    imp <- importance_long[in_block(importance_long), ]
    cond <- conditions[in_block(conditions), ]
    cond <- cond[match(n_values, cond$N), ]

    write_block_header(wb, sheet, r, blk, cond, lay, n_values, styler$mark)
    rank_rows <- setNames(r + 3 + 3 * (seq_along(METHODS) - 1) + 2, METHODS)
    for (k in seq_along(METHODS)) {
      w_row <- r + 3 + 3 * (k - 1)
      write_method_rows(
        wb,
        sheet,
        w_row,
        METHODS[k],
        imp,
        cond,
        lay,
        n_values,
        styler$mark
      )
    }
    add_rank_highlighting(wb, sheet, rank_rows, lay, styles)
    r <- r + 3 + 3 * length(METHODS) + 1
  }

  styler$apply(wb, sheet, styles)
  setColWidths(wb, sheet, cols = 1:2, widths = c(9, 11))
  setColWidths(wb, sheet, cols = 3:lay$last, widths = 8.5)
  setColWidths(wb, sheet, cols = lay$gap_cols, widths = 2)
  freezePane(wb, sheet, firstActiveRow = 4, firstActiveCol = 3)
}

#' Write the k = 0 sheet: for every population matrix, the % of each
#' predictor's population GD-R², GD-R, and GD-sr value contributed at k = 0,
#' with a mean row (Excel formulas) per number of predictors.
write_k0_sheet <- function(wb, k0_wide, p_values, styles) {
  sheet <- "k = 0 share"
  addWorksheet(wb, sheet, gridLines = FALSE)
  styler <- style_collector()
  writeData(
    wb,
    sheet,
    paste0(
      "Share of each predictor's population GD weight contributed at k = 0 ",
      "(its standalone validity)"
    ),
    startCol = 1,
    startRow = 1
  )
  writeData(
    wb,
    sheet,
    paste(
      paste0(
        "k = 0 is the predictor entering an empty model, so its increment is ",
        "r² (GD-R²) or |r| (GD-R and GDsr)."
      ),
      "X1 is the most important predictor in the population (by GD-R²), Xp the least."
    ),
    startCol = 1,
    startRow = 2
  )
  addStyle(wb, sheet, styles$title, rows = 1, cols = 1)
  addStyle(wb, sheet, styles$note, rows = 2, cols = 1)

  groups <- c(GDR2 = "GD-R²", GDR = "GD-R", GDsr = "GDsr")
  r <- 4
  k0_last <- 0
  for (p in p_values) {
    x <- k0_wide[k0_wide$p == p, ]
    x <- x[order(x$rbar, x$vbar, -x$rel), ]
    pops_p <- unique(x[, c("rbar", "vbar", "rel")])
    k_start <- 4 + (seq_along(groups) - 1) * (p + 1)
    last_col <- max(k_start) + p - 1
    k0_last <- max(k0_last, last_col)

    writeData(
      wb,
      sheet,
      sprintf("%d predictors", p),
      startCol = 1,
      startRow = r
    )
    styler$mark("block", r, 1:last_col)
    writeData(
      wb,
      sheet,
      t(c("Mean r", "Mean validity", "Reliability")),
      startCol = 1,
      startRow = r + 2,
      colNames = FALSE
    )
    for (g in seq_along(groups)) {
      mergeCells(
        wb,
        sheet,
        cols = k_start[g]:(k_start[g] + p - 1),
        rows = r + 1
      )
      writeData(wb, sheet, groups[g], startCol = k_start[g], startRow = r + 1)
      writeData(
        wb,
        sheet,
        t(paste0("X", 1:p)),
        startCol = k_start[g],
        startRow = r + 2,
        colNames = FALSE
      )
    }
    styler$mark("head", (r + 1):(r + 2), 1:last_col)

    first <- r + 3
    for (i in seq_len(nrow(pops_p))) {
      row <- first + i - 1
      b <- x[
        x$rbar == pops_p$rbar[i] &
          x$vbar == pops_p$vbar[i] &
          x$rel == pops_p$rel[i],
      ]
      b <- order_predictors(b)
      writeData(
        wb,
        sheet,
        t(c(
          apa_num(pops_p$rbar[i]),
          apa_num(pops_p$vbar[i]),
          apa_num(pops_p$rel[i])
        )),
        startCol = 1,
        startRow = row,
        colNames = FALSE
      )
      for (g in seq_along(groups)) {
        writeData(
          wb,
          sheet,
          t(b[[paste0("k0_", names(groups)[g])]] / 100),
          startCol = k_start[g],
          startRow = row,
          colNames = FALSE
        )
      }
    }
    rows_p <- first:(first + nrow(pops_p) - 1)
    mean_row <- max(rows_p) + 1
    writeData(wb, sheet, "Mean", startCol = 1, startRow = mean_row)
    for (g in seq_along(groups)) {
      cols <- int2col(k_start[g]:(k_start[g] + p - 1))
      write_formula_row(
        wb,
        sheet,
        sprintf("AVERAGE(%s%d:%s%d)", cols, min(rows_p), cols, max(rows_p)),
        mean_row,
        k_start[g]
      )
    }
    k0_cols <- unlist(lapply(k_start, function(c0) c0:(c0 + p - 1)))
    styler$mark("cell", rows_p, 1:3)
    styler$mark("k0", rows_p, k0_cols)
    styler$mark("stat_top", mean_row, 1)
    styler$mark("k0_mean", mean_row, k0_cols)
    r <- mean_row + 2
  }
  styler$apply(wb, sheet, styles)
  setColWidths(wb, sheet, cols = 1:3, widths = c(9, 12, 11))
  setColWidths(wb, sheet, cols = 4:k0_last, widths = 7.5)
  freezePane(wb, sheet, firstActiveRow = 4, firstActiveCol = 4)
}

#' Remove dangling references from a saved .xlsx file.
#'
#' openxlsx lists parts it never writes (empty drawings, shared strings) in the
#' package's relationship files and manifest. Excel tolerates this, but other
#' readers (e.g., Python's openpyxl) reject the file. This removes every
#' reference to a missing part and re-zips the file in place.
#'
#' @param path Path of the .xlsx file.
#' @return `path`, invisibly.
clean_package <- function(path) {
  tmp <- tempfile("xlsx_")
  dir.create(tmp)
  unzip(path, exdir = tmp)
  parts <- list.files(tmp, recursive = TRUE, all.files = TRUE)
  part_exists <- function(p) file.exists(file.path(tmp, p))

  for (f in parts[grepl("\\.rels$", parts)]) {
    xml <- paste(
      readLines(file.path(tmp, f), warn = FALSE, encoding = "UTF-8"),
      collapse = ""
    )
    owner <- dirname(dirname(f))
    for (rel in regmatches(xml, gregexpr("<Relationship [^>]*/>", xml))[[1]]) {
      if (grepl('TargetMode="External"', rel)) {
        next
      }
      target <- sub('.*Target="([^"]+)".*', "\\1", rel)
      path_in_zip <- if (startsWith(target, "/")) {
        sub("^/", "", target)
      } else {
        file.path(owner, target)
      }
      if (!part_exists(path_in_zip)) xml <- sub(rel, "", xml, fixed = TRUE)
    }
    writeLines(xml, file.path(tmp, f), useBytes = TRUE)
  }

  manifest <- file.path(tmp, "[Content_Types].xml")
  xml <- paste(
    readLines(manifest, warn = FALSE, encoding = "UTF-8"),
    collapse = ""
  )
  for (entry in regmatches(xml, gregexpr("<Override [^>]*/>", xml))[[1]]) {
    part <- sub("^/", "", sub('.*PartName="([^"]+)".*', "\\1", entry))
    if (!part_exists(part)) xml <- sub(entry, "", xml, fixed = TRUE)
  }
  writeLines(xml, manifest, useBytes = TRUE)

  zipped <- tempfile(fileext = ".xlsx")
  old_wd <- setwd(tmp)
  zip::zipr(
    zipped,
    files = c("[Content_Types].xml", "_rels", "docProps", "xl"),
    include_directories = FALSE
  )
  setwd(old_wd)
  file.copy(zipped, path, overwrite = TRUE)
  invisible(path)
}

#' Write the manuscript-style workbook: one sheet per number of predictors and
#' the k = 0 sheet.
#'
#' @param path Output file.
#' @param importance_long,conditions,k0_wide Results data frames.
#' @param param Simulation parameters (p and N values).
write_importance_workbook <- function(
  path,
  importance_long,
  conditions,
  k0_wide,
  param
) {
  styles <- workbook_styles()
  wb <- createWorkbook()
  for (p in param$p) {
    write_importance_sheet(wb, p, importance_long, conditions, param$N, styles)
  }
  write_k0_sheet(wb, k0_wide, param$p, styles)
  saveWorkbook(wb, path, overwrite = TRUE)
  clean_package(path)
}


# =============================================================================
# 5. Study 1 tables
# =============================================================================
# Accuracy measures (one value per condition and method):
#   relative bias   100 x (mean sample estimate - population value) /
#                   population value, for the total R² (GD-R² and RW) and
#                   the total R (GD-R)
#   share error     mean over predictors of |population % share - mean
#                   sample % share| in percentage points; sample shares are
#                   computed in each replication and then averaged
#   rank recovery   % of replications reproducing the method's complete
#                   population ranking
#   rank deviation  mean over predictors of |mean sample rank - population
#                   rank|
# Differences between GD-R and GD-R² are coded so that positive values
# favor GD-R.

# Simulation factors, labeled as in the tables.
FACTORS <- c(
  "Sample size" = "N",
  "No. predictors" = "p",
  "Mean inter-predictor r" = "rbar",
  "Mean validity" = "vbar",
  "Reliability" = "rel"
)

#' Accuracy measures for every condition and method.
#'
#' @param conditions Output of condition_table().
#' @param importance_long Output of importance_table().
#' @return `conditions` with share error, rank deviation, relative bias, rank
#'   recovery (%), and the GD-R vs. GD-R² differences added.
accuracy_by_condition <- function(conditions, importance_long) {
  acc <- aggregate(
    cbind(
      share_err = abs(importance_long$pop_pct - importance_long$mean_pct),
      rank_dev = abs(importance_long$mean_rank - importance_long$pop_rank)
    ),
    by = importance_long[c(COND_COLS, "method")],
    FUN = mean
  )
  acc <- reshape(
    acc,
    idvar = COND_COLS,
    timevar = "method",
    direction = "wide",
    sep = "_"
  )
  cond <- merge(conditions, acc, by = COND_COLS)

  cond$bias_R2 <- 100 * (cond$mean_R2 / cond$pop_R2 - 1)
  cond$bias_R <- 100 * (cond$GDR_mean_total / cond$GDR_pop_total - 1)
  for (m in METHODS) {
    cond[[paste0("rec_", m)]] <- 100 * cond[[paste0(m, "_full")]]
  }

  cond$d_bias <- abs(cond$bias_R2) - abs(cond$bias_R)
  cond$d_share <- cond$share_err_GDR2 - cond$share_err_GDR
  cond$d_rec <- cond$rec_GDR - cond$rec_GDR2
  cond$d_dev <- cond$rank_dev_GDR2 - cond$rank_dev_GDR
  cond
}

#' A table with one row per factor level and one column per variable.
#'
#' Rows: "Overall", then each level of each factor. Cells: the mean over the
#' rows of `d` at that level, rounded to `digits`.
#'
#' @param d Data frame with one row per condition (or per predictor).
#' @param columns Variables to average.
#' @param factors Factors to break down by (named: label = column).
#' @param digits Decimal places.
#' @return A data frame with columns Factor, Level, and `columns`.
by_factor <- function(d, columns, factors = FACTORS, digits = 2) {
  row_for <- function(factor, level, keep) {
    data.frame(
      Factor = factor,
      Level = level,
      t(round(colMeans(d[keep, columns, drop = FALSE]), digits)),
      check.names = FALSE
    )
  }
  out <- row_for("Overall", "All conditions", rep(TRUE, nrow(d)))
  for (f in names(factors)) {
    v <- d[[factors[[f]]]]
    for (lv in sort(unique(v))) {
      label <- if (factors[[f]] %in% c("N", "p")) {
        format(lv, big.mark = ",")
      } else {
        apa_num(lv)
      }
      out <- rbind(out, row_for(f, label, v == lv))
    }
  }
  out$Factor[duplicated(out$Factor)] <- ""
  out
}

#' The Study 1 tables.
#'
#' Table 3 (Panels A and B), Table 4 (standardized regressions of the four
#' GD-R vs. GD-R² differences on the design factors, N log10-transformed,
#' with each model's R² and R), and Appendix Tables B3, E1, and E2.
#'
#' @param cond Output of accuracy_by_condition().
#' @param k0_wide Output of k0_table().
#' @return A named list; each element has a `title` and a `data` data frame.
study1_tables <- function(cond, k0_wide) {
  table3A <- by_factor(
    cond,
    c(
      "bias_R2",
      "bias_R",
      "d_bias",
      "share_err_GDR2",
      "share_err_GDR",
      "d_share"
    )
  )
  names(table3A)[-(1:2)] <- c(
    "R² relative bias (%)",
    "R relative bias (%)",
    "Difference (pp)",
    "GD-R² share error (pp)",
    "GD-R share error (pp)",
    "Difference (pp)"
  )

  table3B <- cbind(
    by_factor(cond, c("rec_GDR2", "rec_GDR", "d_rec")),
    by_factor(cond, c("rank_dev_GDR2", "rank_dev_GDR", "d_dev"), digits = 3)[,
      -(1:2)
    ]
  )
  names(table3B)[-(1:2)] <- c(
    "GD-R² rank recovery (%)",
    "GD-R rank recovery (%)",
    "Difference (pp)",
    "GD-R² rank deviation",
    "GD-R rank deviation",
    "Difference"
  )

  Z <- function(v) (v - mean(v)) / sd(v)
  cond$log10N <- log10(cond$N)
  coefs <- sapply(c("d_bias", "d_share", "d_rec", "d_dev"), function(y) {
    fit <- lm(
      Z(cond[[y]]) ~ Z(log10N) + Z(p) + Z(vbar) + Z(rbar) + Z(rel),
      data = cond
    )
    r2 <- summary(fit)$r.squared
    round(c(coef(fit)[-1], r2, sqrt(r2)), 3)
  })
  table4 <- data.frame(
    Predictor = c(
      "log10(N)",
      "Number of predictors",
      "Mean validity",
      "Mean inter-predictor r",
      "Reliability",
      "Model R²",
      "Model R"
    ),
    coefs,
    check.names = FALSE,
    row.names = NULL
  )
  names(table4)[-1] <- c(
    "Absolute % bias difference",
    "Share-error difference",
    "Rank-recovery difference",
    "Rank-deviation difference"
  )

  k0_wide$k0_diff <- k0_wide$k0_GDR - k0_wide$k0_GDR2
  tableB3 <- by_factor(
    k0_wide,
    c("k0_GDR2", "k0_GDR", "k0_diff"),
    FACTORS[-1],
    digits = 1
  )
  names(tableB3)[-(1:2)] <- c("GD-R² (%)", "GD-R (%)", "Difference (pp)")

  tableE1 <- by_factor(
    cond,
    c("bias_R2", "share_err_RW", "share_err_GDR2", "rec_RW", "rec_GDR2")
  )
  names(tableE1)[-(1:2)] <- c(
    "R² / RW relative bias (%)",
    "RW share error (pp)",
    "GD-R² share error (pp)",
    "RW rank recovery (%)",
    "GD-R² rank recovery (%)"
  )

  cond$sr_ratio <- cond$GDsr_mean_total / cond$GDR_mean_total
  tableE2 <- cbind(
    by_factor(cond, c("GDR_mean_total", "GDsr_mean_total"), digits = 3),
    by_factor(cond, "sr_ratio")[, -(1:2), drop = FALSE]
  )
  names(tableE2)[-(1:2)] <- c(
    "GD-R total (= model R)",
    "GD-sr total",
    "GD-sr total / GD-R total"
  )

  list(
    "Table 3A" = list(
      title = "Table 3, Panel A. Model-metric bias and percentage-share error",
      data = table3A
    ),
    "Table 3B" = list(
      title = "Table 3, Panel B. Rank recovery and rank deviation",
      data = table3B
    ),
    "Table 4" = list(
      title = paste(
        "Table 4. Predicting differences in sampling performance between",
        "GD-R and GD-R² (standardized coefficients)"
      ),
      data = table4
    ),
    "Table B3" = list(
      title = "Table B3. Contribution of the k = 0 level to GD-R and GD-R² (population values)",
      data = tableB3
    ),
    "Table E1" = list(
      title = "Table E1. Comparison between GD-R² and RW",
      data = tableE1
    ),
    "Table E2" = list(
      title = "Table E2. Comparison between GD-R and GD-sr",
      data = tableE2
    )
  )
}

#' Write the Study 1 tables to a workbook, one sheet per table.
#'
#' @param path Output file.
#' @param tables Output of study1_tables().
write_study1_workbook <- function(path, tables) {
  wb <- createWorkbook()
  for (sheet in names(tables)) {
    addWorksheet(wb, sheet, gridLines = FALSE)
    writeData(wb, sheet, tables[[sheet]]$title, startRow = 1)
    writeData(
      wb,
      sheet,
      tables[[sheet]]$data,
      startRow = 3,
      headerStyle = createStyle(textDecoration = "bold")
    )
    addStyle(
      wb,
      sheet,
      createStyle(textDecoration = "bold"),
      rows = 1,
      cols = 1
    )
    setColWidths(
      wb,
      sheet,
      cols = seq_along(tables[[sheet]]$data),
      widths = "auto"
    )
  }
  saveWorkbook(wb, path, overwrite = TRUE)
  clean_package(path)
}

#' Print the Study 1 tables to the console.
print_tables <- function(tables) {
  old <- options(width = 200)
  on.exit(options(old))
  for (sheet in names(tables)) {
    cat("\n== ", tables[[sheet]]$title, " ==\n", sep = "")
    print(tables[[sheet]]$data, row.names = FALSE, right = FALSE)
  }
}


# =============================================================================
# Main
# =============================================================================

design <- createDesign(
  N = PARAM$N,
  p = PARAM$p,
  rbar = PARAM$rbar,
  vbar = PARAM$vbar,
  rel = PARAM$rel
)
subsets <- setNames(lapply(PARAM$p, make_subsets), PARAM$p)
populations <- build_populations(design, PARAM, subsets)

dir.create(OUTPUT_DIR, showWarnings = FALSE)

if (RUN_SIMULATION) {
  res <- run_study(design, PARAM, populations$observed, subsets)
  saveRDS(
    list(
      param = PARAM,
      populations = populations$error_free,
      observed = populations$observed,
      results = res
    ),
    RESULTS_RDS
  )
} else {
  res <- readRDS(RESULTS_RDS)$results
}

conditions <- condition_table(res)
importance_long <- importance_table(res)
k0_wide <- k0_table(importance_long, PARAM$N[1])

write_importance_workbook(
  RESULTS_XLSX,
  importance_long,
  conditions,
  k0_wide,
  PARAM
)

tables <- study1_tables(
  accuracy_by_condition(conditions, importance_long),
  k0_wide
)
write_study1_workbook(TABLES_XLSX, tables)

print_tables(tables)
cat(sprintf(
  "\nSaved in %s/: %s, %s, and %s\n",
  OUTPUT_DIR,
  basename(RESULTS_RDS),
  basename(RESULTS_XLSX),
  basename(TABLES_XLSX)
))
