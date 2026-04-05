#!/usr/bin/env Rscript
###############################################################################
#  G25 Bayesian MCMC Admixture Estimator
#  ──────────────────────────────────────
#  Estimates ancestry proportions from Global25 PCA coordinates using
#  Bayesian inference with a Dirichlet prior and Markov Chain Monte Carlo
#  (Metropolis-Hastings) sampling.
#
#  REQUIRED PACKAGES: ggplot2, patchwork, reshape2
#  (Auto-installed if missing)
#
#  Unlike least-squares tools (e.g. Vahaduo), this produces full posterior
#  distributions over admixture proportions, giving you both point estimates
#  AND credible intervals so you know which signals are real vs. noise.
#
#  USAGE:
#    Rscript g25_bayesian_mcmc.R \
#      --source   source_populations.csv \
#      --target   target_individual.csv \
#      --mode     population          # "population" or "sample"
#      --out      results             # output file prefix
#      --iter     50000               # MCMC iterations (default 50000)
#      --burnin   10000               # burn-in to discard  (default 10000)
#      --thin     10                  # thinning interval   (default 10)
#      --max_k    0                   # max sources (0 = auto, try 2..12)
#      --n_keep   25                  # pre-selection candidates (default 25)
#      --sigma    0.01                # likelihood noise sigma (default 0.01)
#      --alpha    1.0                 # Dirichlet concentration (1=uniform)
#      --seed     42                  # RNG seed
#
#  INPUT FORMAT:
#    CSV files with no header row.
#    Column 1 = Population:SampleName  (colon-delimited)
#    Columns 2-26 = G25 coordinates (25 dimensions)
#
#  OUTPUTS:
#    <out>_summary.csv      – point estimates + 95% credible intervals
#    <out>_posterior.csv     – full posterior samples (post burn-in/thin)
#    <out>_diagnostics.txt  – convergence diagnostics
#    <out>_plots.pdf        – posterior density + trace plots (ggplot2)
###############################################################################


###############################################################################
# 0. PACKAGE MANAGEMENT & ARGUMENT PARSING
###############################################################################

# Auto-install required packages if missing
required_pkgs <- c("ggplot2", "patchwork", "reshape2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
library(ggplot2)
library(patchwork)
library(reshape2)

# ── Custom ggplot2 theme for all plots ──
theme_g25 <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans", colour = "#2d3436"),
      plot.title       = element_text(face = "bold", size = rel(1.2),
                                      hjust = 0, margin = margin(b = 10)),
      plot.subtitle    = element_text(colour = "#636e72", size = rel(0.85),
                                      hjust = 0, margin = margin(b = 12)),
      plot.caption     = element_text(colour = "#b2bec3", size = rel(0.7),
                                      hjust = 1, margin = margin(t = 10)),
      panel.grid.major.y = element_line(colour = "#dfe6e9", linewidth = 0.4),
      panel.grid.major.x = element_line(colour = "#dfe6e9", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.title       = element_text(face = "bold", size = rel(0.85)),
      axis.text        = element_text(size = rel(0.8)),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = rel(0.8)),
      legend.text      = element_text(size = rel(0.75)),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.85)),
      plot.margin      = margin(15, 15, 15, 15)
    )
}

# Colour palette — professional, colourblind-friendly
g25_palette <- c(
  "#2166AC", "#B2182B", "#4DAF4A", "#FF7F00", "#984EA3",
  "#A65628", "#E7298A", "#66A61E", "#E6AB02", "#7570B3",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"
)


parse_arguments <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)
  
  defaults <- list(
    source = NULL,
    target = NULL,
    mode   = "population",
    out    = "results",
    iter   = 50000L,
    burnin = 10000L,
    thin   = 10L,
    max_k  = 0L,
    n_keep = 25L,
    sigma  = 0.01,
    alpha  = 1.0,
    seed   = 42L
  )
  
  args <- defaults
  i <- 1
  while (i <= length(raw_args)) {
    key <- sub("^--", "", raw_args[i])
    if (key %in% names(args) && i < length(raw_args)) {
      val <- raw_args[i + 1]
      if (is.integer(defaults[[key]])) {
        args[[key]] <- as.integer(val)
      } else if (is.double(defaults[[key]]) && !is.integer(defaults[[key]])) {
        args[[key]] <- as.double(val)
      } else {
        args[[key]] <- val
      }
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  if (is.null(args$source) || is.null(args$target)) {
    cat("Usage: Rscript g25_bayesian_mcmc.R --source <file> --target <file> [options]\n\n")
    cat("Required:\n")
    cat("  --source   CSV of source/reference G25 coordinates\n")
    cat("  --target   CSV of target individual(s) G25 coordinates\n\n")
    cat("Options:\n")
    cat("  --mode     'population' (average per pop) or 'sample' [default: population]\n")
    cat("  --out      Output file prefix                         [default: results]\n")
    cat("  --iter     Total MCMC iterations                      [default: 50000]\n")
    cat("  --burnin   Burn-in iterations to discard              [default: 10000]\n")
    cat("  --thin     Thinning interval                          [default: 10]\n")
    cat("  --max_k    Max source components (0 = auto via BIC)   [default: 0]\n")
    cat("  --n_keep   Pre-selection candidates to keep           [default: 25]\n")
    cat("  --sigma    Noise std dev in likelihood                [default: 0.01]\n")
    cat("  --alpha    Dirichlet prior concentration              [default: 1.0]\n")
    cat("  --seed     Random seed                                [default: 42]\n")
    stop("Both --source and --target are required.", call. = FALSE)
  }
  
  args
}


###############################################################################
# 1. DATA LOADING & PREPARATION
###############################################################################

read_g25 <- function(filepath) {
  raw <- read.csv(filepath, header = FALSE, stringsAsFactors = FALSE,
                  strip.white = TRUE)
  labels <- trimws(raw[, 1])
  coords <- as.matrix(raw[, 2:ncol(raw)])
  
  split_labels <- strsplit(labels, ":", fixed = TRUE)
  populations  <- sapply(split_labels, `[`, 1)
  sample_names <- sapply(split_labels, function(x) {
    if (length(x) > 1) paste(x[-1], collapse = ":") else x[1]
  })
  
  list(labels = labels, populations = populations,
       samples = sample_names, coords = coords)
}


###############################################################################
# 2. BAYESIAN MODEL
###############################################################################

log_likelihood <- function(w, target_vec, source_mat, sigma) {
  predicted <- as.vector(w %*% source_mat)
  residuals <- target_vec - predicted
  -sum(residuals^2) / (2 * sigma^2)
}

log_prior_dirichlet <- function(w, alpha) {
  if (any(w <= 0)) return(-Inf)
  (alpha - 1) * sum(log(w))
}

log_posterior <- function(w, target_vec, source_mat, sigma, alpha) {
  lp <- log_prior_dirichlet(w, alpha)
  if (is.infinite(lp)) return(-Inf)
  lp + log_likelihood(w, target_vec, source_mat, sigma)
}

propose_dirichlet <- function(w_current, concentration = 500) {
  alpha_prop <- concentration * w_current
  alpha_prop[alpha_prop < 0.01] <- 0.01
  proposed <- rgamma(length(alpha_prop), shape = alpha_prop, rate = 1)
  proposed <- proposed / sum(proposed)
  list(w = proposed, alpha_used = alpha_prop)
}

log_dirichlet_density <- function(x, alpha_vec) {
  if (any(x <= 0)) return(-Inf)
  sum((alpha_vec - 1) * log(x)) - sum(lgamma(alpha_vec)) + lgamma(sum(alpha_vec))
}


###############################################################################
# 3. PRE-SELECTION
###############################################################################

solve_nnls <- function(target_vec, source_mat, max_iter = 2000, tol = 1e-8) {
  K <- nrow(source_mat)
  D <- ncol(source_mat)
  StS <- tcrossprod(source_mat)
  Sty <- as.vector(source_mat %*% target_vec)
  w <- rep(0, K)
  
  for (iter in 1:max_iter) {
    w_old <- w
    for (k in 1:K) {
      residual_k <- Sty[k] - sum(StS[k, ] * w) + StS[k, k] * w[k]
      w[k] <- max(0, residual_k / StS[k, k])
    }
    s <- sum(w)
    if (s > 0) w <- w / s
    if (max(abs(w - w_old)) < tol) break
  }
  w
}

residual_chase <- function(target_vec, source_mat, source_names, n_rounds = 10) {
  selected <- integer(0)
  residual <- target_vec
  
  for (round in 1:n_rounds) {
    rss_reduction <- numeric(nrow(source_mat))
    for (j in 1:nrow(source_mat)) {
      if (j %in% selected) { rss_reduction[j] <- -Inf; next }
      trial_set <- c(selected, j)
      trial_mat <- source_mat[trial_set, , drop = FALSE]
      w_trial <- solve_nnls(target_vec, trial_mat, max_iter = 500)
      pred <- as.vector(w_trial %*% trial_mat)
      rss_reduction[j] <- -sum((target_vec - pred)^2)
    }
    
    best <- which.max(rss_reduction)
    if (best %in% selected) break
    selected <- c(selected, best)
    
    trial_mat <- source_mat[selected, , drop = FALSE]
    w_trial <- solve_nnls(target_vec, trial_mat, max_iter = 500)
    residual <- target_vec - as.vector(w_trial %*% trial_mat)
    
    if (sqrt(sum(residual^2)) < 1e-6) break
  }
  selected
}

preselect_sources <- function(target_vec, source_mat, source_names, n_keep = 25) {
  K <- nrow(source_mat)
  cat("    Running multi-strategy pre-selection on", K, "sources...\n")
  
  dists <- apply(source_mat, 1, function(s) sqrt(sum((target_vec - s)^2)))
  dist_top <- order(dists)[1:min(n_keep, K)]
  cat("      Distance-based:", length(dist_top), "candidates\n")
  
  w_nnls <- solve_nnls(target_vec, source_mat, max_iter = 2000)
  nnls_active <- which(w_nnls > 1e-4)
  nnls_top <- order(-w_nnls)[1:min(n_keep, K)]
  nnls_set <- unique(c(nnls_active, nnls_top))
  cat("      NNLS-based:", length(nnls_active), "active (weight > 0.01%),",
      length(nnls_set), "total candidates\n")
  
  chase_set <- residual_chase(target_vec, source_mat, source_names, n_rounds = 12)
  cat("      Residual-chase:", length(chase_set), "candidates\n")
  
  target_norm <- target_vec / (sqrt(sum(target_vec^2)) + 1e-15)
  cos_sim <- apply(source_mat, 1, function(s) {
    sum(s * target_norm) / (sqrt(sum(s^2)) + 1e-15)
  })
  dir_top <- order(-cos_sim)[1:min(ceiling(n_keep/3), K)]
  
  all_candidates <- unique(c(nnls_set, chase_set, dist_top, dir_top))
  
  if (length(all_candidates) > n_keep) {
    scores <- w_nnls[all_candidates] * 1000 + 1 / (dists[all_candidates] + 0.001)
    keep_idx <- order(-scores)[1:n_keep]
    all_candidates <- all_candidates[keep_idx]
  }
  
  cat("    => Pre-selected", length(all_candidates), "candidates from", K, "total\n")
  cat("    Candidates:", paste(source_names[all_candidates], collapse = ", "), "\n")
  
  sub_w <- solve_nnls(target_vec, source_mat[all_candidates, , drop = FALSE], max_iter = 2000)
  
  list(indices = all_candidates, initial_weights = sub_w)
}


###############################################################################
# 4. MODEL SELECTION
###############################################################################

compute_bic <- function(target_vec, source_mat, w, sigma, K_active) {
  predicted <- as.vector(w %*% source_mat)
  rss <- sum((target_vec - predicted)^2)
  n   <- length(target_vec)
  log_lik <- -n/2 * log(2 * pi * sigma^2) - rss / (2 * sigma^2)
  -2 * log_lik + (K_active - 1) * log(n)
}

forward_select <- function(target_vec, source_mat, source_names, sigma,
                           max_k = 12, min_k = 2) {
  K_avail <- nrow(source_mat)
  max_k <- min(max_k, K_avail)
  
  selected  <- integer(0)
  remaining <- 1:K_avail
  best_bic_overall <- Inf
  best_set  <- NULL
  bic_trail <- numeric(0)
  step_log  <- character(0)
  
  cat("    Forward stepwise model selection (up to K =", max_k, ")...\n")
  
  for (step in 1:max_k) {
    best_bic_step <- Inf
    best_add <- NA
    
    for (j in remaining) {
      trial <- c(selected, j)
      trial_mat <- source_mat[trial, , drop = FALSE]
      w_trial <- solve_nnls(target_vec, trial_mat, max_iter = 1000)
      bic_j <- compute_bic(target_vec, trial_mat, w_trial, sigma, length(trial))
      
      if (bic_j < best_bic_step) {
        best_bic_step <- bic_j
        best_add <- j
      }
    }
    
    selected  <- c(selected, best_add)
    remaining <- setdiff(remaining, best_add)
    bic_trail <- c(bic_trail, best_bic_step)
    
    trial_mat <- source_mat[selected, , drop = FALSE]
    w_trial <- solve_nnls(target_vec, trial_mat, max_iter = 1000)
    pred <- as.vector(w_trial %*% trial_mat)
    rmse <- sqrt(mean((target_vec - pred)^2))
    
    step_line <- sprintf("K=%2d  BIC=%8.2f  RMSE=%.6f  +%-30s (%.1f%%)",
                         step, best_bic_step, rmse, source_names[best_add],
                         w_trial[length(w_trial)] * 100)
    step_log <- c(step_log, step_line)
    cat("     ", step_line, "\n")
    
    if (best_bic_step < best_bic_overall) {
      best_bic_overall <- best_bic_step
      best_set <- selected
    }
    
    if (length(bic_trail) >= 4) {
      recent <- tail(bic_trail, 3)
      if (all(diff(recent) > 0)) {
        cat("      (BIC increasing for 3 steps, stopping)\n")
        step_log <- c(step_log, "(BIC increasing for 3 steps, stopping)")
        break
      }
    }
  }
  
  if (length(best_set) < min_k) best_set <- selected[1:min(min_k, length(selected))]
  
  cat("    => Selected K =", length(best_set), "sources\n")
  
  list(indices = best_set,
       weights = solve_nnls(target_vec,
                            source_mat[best_set, , drop = FALSE],
                            max_iter = 2000),
       step_log = step_log)
}


###############################################################################
# 5. MCMC SAMPLER
###############################################################################

run_mcmc <- function(target_vec, source_mat, source_names,
                     n_iter, burnin, thin, sigma, alpha,
                     proposal_concentration = 500) {
  
  K <- nrow(source_mat)
  w_current <- rep(1/K, K)
  
  n_save   <- floor((n_iter - burnin) / thin)
  chain    <- matrix(0, nrow = n_save, ncol = K)
  colnames(chain) <- source_names
  log_post_trace  <- numeric(n_iter)
  accept_count    <- 0
  save_idx        <- 0
  
  lp_current <- log_posterior(w_current, target_vec, source_mat, sigma, alpha)
  conc <- proposal_concentration
  
  cat("    Running MCMC:", n_iter, "iterations (burn-in:", burnin,
      ", thin:", thin, ")\n")
  
  for (i in 1:n_iter) {
    prop <- propose_dirichlet(w_current, conc)
    w_proposed <- prop$w
    
    lp_proposed <- log_posterior(w_proposed, target_vec, source_mat, sigma, alpha)
    
    log_q_fwd <- log_dirichlet_density(w_proposed, conc * w_current)
    log_q_rev <- log_dirichlet_density(w_current,  conc * w_proposed)
    
    log_alpha_mh <- (lp_proposed - lp_current) + (log_q_rev - log_q_fwd)
    
    if (log(runif(1)) < log_alpha_mh) {
      w_current  <- w_proposed
      lp_current <- lp_proposed
      accept_count <- accept_count + 1
    }
    
    log_post_trace[i] <- lp_current
    
    if (i > burnin && ((i - burnin) %% thin == 0)) {
      save_idx <- save_idx + 1
      if (save_idx <= n_save) {
        chain[save_idx, ] <- w_current
      }
    }
    
    if (i <= burnin && i %% 500 == 0) {
      recent_rate <- accept_count / i
      if (recent_rate < 0.20) conc <- conc * 0.8
      else if (recent_rate > 0.50) conc <- conc * 1.2
      conc <- max(50, min(conc, 5000))
    }
    
    if (i %% 10000 == 0) {
      cat("      Iteration", i, "/", n_iter,
          " | Accept rate:", round(accept_count/i, 3),
          " | Proposal conc:", round(conc, 1), "\n")
    }
  }
  
  cat("    Final acceptance rate:", round(accept_count / n_iter, 3), "\n")
  
  list(chain = chain[1:save_idx, , drop = FALSE],
       log_posterior_trace = log_post_trace,
       acceptance_rate = accept_count / n_iter,
       final_concentration = conc)
}


###############################################################################
# 6. SUMMARISE POSTERIORS
###############################################################################

summarise_posterior <- function(chain, source_names, threshold = 0.005) {
  means   <- colMeans(chain)
  medians <- apply(chain, 2, median)
  sds     <- apply(chain, 2, sd)
  ci_lo   <- apply(chain, 2, quantile, probs = 0.025)
  ci_hi   <- apply(chain, 2, quantile, probs = 0.975)
  
  df <- data.frame(
    Source     = source_names,
    Mean       = round(means, 6),
    Median     = round(medians, 6),
    SD         = round(sds, 6),
    CI_2.5     = round(ci_lo, 6),
    CI_97.5    = round(ci_hi, 6),
    Pct_Mean   = round(means * 100, 2),
    Pct_CI_lo  = round(ci_lo * 100, 2),
    Pct_CI_hi  = round(ci_hi * 100, 2),
    stringsAsFactors = FALSE
  )
  
  df <- df[order(-df$Mean), ]
  df$Significant <- df$CI_2.5 > threshold
  df
}


###############################################################################
# 7. CONVERGENCE DIAGNOSTICS
###############################################################################

effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 10) return(n)
  
  x_centred <- x - mean(x)
  max_lag <- min(n - 1, floor(n / 2))
  
  var_x <- sum(x_centred^2) / n
  if (var_x == 0) return(n)
  
  acf_vals <- numeric(max_lag + 1)
  for (lag in 0:max_lag) {
    acf_vals[lag + 1] <- sum(x_centred[1:(n - lag)] * x_centred[(1 + lag):n]) / (n * var_x)
  }
  
  tau <- 1
  k <- 1
  while (k + 1 <= max_lag) {
    pair_sum <- acf_vals[k + 1] + acf_vals[k + 2]
    if (pair_sum < 0) break
    tau <- tau + 2 * pair_sum
    k <- k + 2
  }
  
  max(1, floor(n / tau))
}

geweke_z <- function(x, frac1 = 0.1, frac2 = 0.5) {
  n <- length(x)
  n1 <- floor(n * frac1)
  n2 <- floor(n * frac2)
  
  x1 <- x[1:n1]
  x2 <- x[(n - n2 + 1):n]
  
  m1 <- mean(x1); m2 <- mean(x2)
  s1 <- var(x1) / max(1, effective_sample_size(x1))
  s2 <- var(x2) / max(1, effective_sample_size(x2))
  
  denom <- sqrt(s1 + s2)
  if (denom == 0) return(0)
  (m1 - m2) / denom
}

run_diagnostics <- function(chain, outfile) {
  means <- colMeans(chain)
  active <- which(means > 0.01)
  
  diag_data <- list(
    chain_length = nrow(chain),
    ess = list(),
    geweke = list(),
    active_names = character(0)
  )
  
  sink(outfile)
  cat("================================================================\n")
  cat("  MCMC Convergence Diagnostics\n")
  cat("================================================================\n\n")
  
  if (length(active) == 0) {
    cat("No sources with mean > 1%. Diagnostics skipped.\n")
    sink()
    return(invisible(diag_data))
  }
  
  diag_data$active_names <- colnames(chain)[active]
  
  cat("Effective Sample Sizes:\n")
  for (j in active) {
    ess <- effective_sample_size(chain[, j])
    diag_data$ess[[colnames(chain)[j]]] <- ess
    cat(sprintf("  %-45s  ESS = %6.0f\n", colnames(chain)[j], ess))
  }
  
  cat("\nGeweke Diagnostic (|z| < 2 suggests convergence):\n")
  for (j in active) {
    z <- geweke_z(chain[, j])
    status <- ifelse(abs(z) < 2, "  OK", "  !! CHECK")
    diag_data$geweke[[colnames(chain)[j]]] <- z
    cat(sprintf("  %-45s  z = %6.3f %s\n", colnames(chain)[j], z, status))
  }
  
  cat("\nChain length (post burn-in & thin):", nrow(chain), "samples\n")
  
  sink()
  cat("  Diagnostics written to:", outfile, "\n")
  
  invisible(diag_data)
}


###############################################################################
# 8. PLOTTING — ggplot2 enhanced multi-page PDF
###############################################################################

make_plots <- function(chain, summary_df, log_post_trace, target_label,
                       outfile, diag_data = NULL, run_info = NULL,
                       mcmc_result = NULL) {
  
  # ── Identify active sources ──
  active <- summary_df$Source[summary_df$Pct_Mean > 0.5]
  if (length(active) == 0) active <- summary_df$Source[1:min(5, nrow(summary_df))]
  n_active <- length(active)
  
  plot_df <- summary_df[summary_df$Source %in% active, ]
  # Factor ordering by mean (largest on top)
  plot_df$Source <- factor(plot_df$Source, levels = rev(plot_df$Source))
  
  # Assign consistent colours per source
  src_colours <- setNames(
    g25_palette[seq_along(levels(plot_df$Source)) %% length(g25_palette) + 1],
    levels(plot_df$Source)
  )
  # Re-assign in sorted order (largest first in palette)
  src_colours <- setNames(
    g25_palette[1:nlevels(plot_df$Source)],
    rev(levels(plot_df$Source))  # largest mean gets first colour
  )
  
  # Truncate long target label for titles
  short_label <- if (nchar(target_label) > 60) {
    paste0(substr(target_label, 1, 57), "...")
  } else {
    target_label
  }
  
  pdf(outfile, width = 11, height = 8.5)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 1: Horizontal bar chart (Vahaduo-style)
  # ════════════════════════════════════════════════════════════════════════
  p1 <- ggplot(plot_df, aes(x = Pct_Mean, y = Source, fill = Source)) +
    # CI range bar (lighter, behind)
    geom_segment(aes(x = Pct_CI_lo, xend = Pct_CI_hi,
                     y = Source, yend = Source, colour = Source),
                 linewidth = 6, alpha = 0.2, lineend = "round",
                 show.legend = FALSE) +
    # Main bar
    geom_col(width = 0.6, show.legend = FALSE) +
    # Mean percentage label
    geom_text(aes(label = sprintf("%.1f%%", Pct_Mean),
                  x = ifelse(Pct_Mean > 25, Pct_Mean / 2, Pct_Mean + 1.5),
                  colour = ifelse(Pct_Mean > 25, "white_text", Source)),
              fontface = "bold", size = 4.2, show.legend = FALSE) +
    # CI annotation
    geom_text(aes(label = sprintf("[%.1f%% \u2013 %.1f%%]", Pct_CI_lo, Pct_CI_hi),
                  x = ifelse(Pct_Mean > 25, Pct_Mean / 2, Pct_Mean + 1.5)),
              vjust = 2.3, size = 2.8, colour = "grey50",
              show.legend = FALSE) +
    scale_fill_manual(values = src_colours) +
    scale_colour_manual(values = c(src_colours, "white_text" = "white")) +
    scale_x_continuous(limits = c(0, 105), expand = c(0, 0),
                       breaks = seq(0, 100, 20)) +
    labs(title = paste("Admixture Composition \u2014", short_label),
         subtitle = "Solid bar = posterior mean  |  Shaded band = 95% credible interval",
         x = "Proportion (%)", y = NULL) +
    theme_g25(base_size = 13) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", size = 11))
  
  print(p1)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 2: Forest plot with CIs
  # ════════════════════════════════════════════════════════════════════════
  
  # Build significance labels
  plot_df$sig_label <- ifelse(plot_df$Significant, "sig", "n.s.")
  plot_df$sig_colour <- ifelse(plot_df$Significant, "#2166AC", "#B2182B")
  
  p2 <- ggplot(plot_df, aes(x = Pct_Mean, y = Source)) +
    # CI line
    geom_segment(aes(x = Pct_CI_lo, xend = Pct_CI_hi,
                     y = Source, yend = Source),
                 colour = "#2166AC", linewidth = 1.8, lineend = "round") +
    # CI endpoint caps
    geom_point(aes(x = Pct_CI_lo), shape = "|", size = 3.5, colour = "#2166AC") +
    geom_point(aes(x = Pct_CI_hi), shape = "|", size = 3.5, colour = "#2166AC") +
    # Mean point
    geom_point(size = 4, colour = "#2166AC", fill = "#2166AC", shape = 21) +
    # Mean label above
    geom_text(aes(label = sprintf("%.1f%%", Pct_Mean)),
              nudge_y = 0.3, size = 3.5, fontface = "bold",
              colour = "#2166AC") +
    # CI low/high labels
    geom_text(aes(x = Pct_CI_lo, label = sprintf("%.1f%%", Pct_CI_lo)),
              nudge_y = -0.3, size = 2.8, colour = "grey45") +
    geom_text(aes(x = Pct_CI_hi, label = sprintf("%.1f%%", Pct_CI_hi)),
              nudge_y = -0.3, size = 2.8, colour = "grey45") +
    # Significance marker at right
    geom_text(aes(x = max(plot_df$Pct_CI_hi) + 5,
                  label = sig_label, colour = sig_colour),
              fontface = "italic", size = 3.2, show.legend = FALSE) +
    scale_colour_identity() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(title = paste("Estimates + 95% Credible Intervals \u2014", short_label),
         subtitle = "sig = 95% CI excludes zero  |  n.s. = not significant",
         x = "Proportion (%)", y = NULL) +
    theme_g25(base_size = 13) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", size = 11))
  
  print(p2)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 3: Posterior density plots (ridgeline-style stacked or faceted)
  # ════════════════════════════════════════════════════════════════════════
  
  # Melt chain data for ggplot
  chain_active <- chain[, active, drop = FALSE] * 100
  chain_long <- reshape2::melt(as.data.frame(chain_active),
                               variable.name = "Source",
                               value.name = "Proportion")
  chain_long$Source <- factor(chain_long$Source, levels = rev(active))
  
  # Compute summary stats for annotation
  annot_df <- do.call(rbind, lapply(active, function(s) {
    x <- chain[, s] * 100
    data.frame(
      Source = s,
      mean_val = mean(x),
      ci_lo = quantile(x, 0.025),
      ci_hi = quantile(x, 0.975),
      stringsAsFactors = FALSE
    )
  }))
  annot_df$Source <- factor(annot_df$Source, levels = rev(active))
  annot_df$label <- sprintf("Mean: %.1f%%  [%.1f%% \u2013 %.1f%%]",
                            annot_df$mean_val, annot_df$ci_lo, annot_df$ci_hi)
  
  if (n_active <= 6) {
    # Faceted density plots — more detail visible
    p3 <- ggplot(chain_long, aes(x = Proportion)) +
      geom_density(aes(fill = Source), colour = NA, alpha = 0.5,
                   adjust = 1.2, show.legend = FALSE) +
      geom_density(aes(colour = Source), linewidth = 0.8,
                   adjust = 1.2, show.legend = FALSE) +
      geom_vline(data = annot_df, aes(xintercept = mean_val),
                 colour = "#2d3436", linewidth = 0.7, linetype = "solid") +
      geom_vline(data = annot_df, aes(xintercept = ci_lo),
                 colour = "grey50", linewidth = 0.5, linetype = "dashed") +
      geom_vline(data = annot_df, aes(xintercept = ci_hi),
                 colour = "grey50", linewidth = 0.5, linetype = "dashed") +
      geom_text(data = annot_df,
                aes(x = mean_val, y = Inf, label = label),
                vjust = 1.5, hjust = 0.5, size = 3, fontface = "bold",
                colour = "#2d3436") +
      facet_wrap(~ Source, scales = "free", ncol = 2) +
      scale_fill_manual(values = src_colours) +
      scale_colour_manual(values = src_colours) +
      labs(title = paste("Posterior Distributions \u2014", short_label),
           subtitle = "Solid line = mean  |  Dashed lines = 95% credible interval bounds",
           x = "Proportion (%)", y = "Density") +
      theme_g25(base_size = 11) +
      theme(strip.text = element_text(face = "bold", size = 11))
    
  } else {
    # For many sources, use overlaid ridgeline-style via facets
    p3 <- ggplot(chain_long, aes(x = Proportion)) +
      geom_density(aes(fill = Source), colour = NA, alpha = 0.4,
                   adjust = 1.2, show.legend = FALSE) +
      geom_density(aes(colour = Source), linewidth = 0.6,
                   adjust = 1.2, show.legend = FALSE) +
      geom_vline(data = annot_df, aes(xintercept = mean_val),
                 colour = "#2d3436", linewidth = 0.6) +
      geom_vline(data = annot_df, aes(xintercept = ci_lo),
                 colour = "grey55", linewidth = 0.4, linetype = "dashed") +
      geom_vline(data = annot_df, aes(xintercept = ci_hi),
                 colour = "grey55", linewidth = 0.4, linetype = "dashed") +
      facet_wrap(~ Source, scales = "free", ncol = 2) +
      scale_fill_manual(values = src_colours) +
      scale_colour_manual(values = src_colours) +
      labs(title = paste("Posterior Distributions \u2014", short_label),
           subtitle = "Solid line = mean  |  Dashed = 95% CI",
           x = "Proportion (%)", y = "Density") +
      theme_g25(base_size = 10)
  }
  
  print(p3)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 4: Log-posterior trace plot
  # ════════════════════════════════════════════════════════════════════════
  
  # Downsample to max 2000 points to keep plotting fast
  n_total <- length(log_post_trace)
  max_pts <- 2000
  if (n_total > max_pts) {
    keep_idx <- unique(round(seq(1, n_total, length.out = max_pts)))
    trace_df <- data.frame(
      Iteration    = keep_idx,
      LogPosterior = log_post_trace[keep_idx]
    )
  } else {
    trace_df <- data.frame(
      Iteration    = seq_along(log_post_trace),
      LogPosterior = log_post_trace
    )
  }
  
  # Pre-compute a rolling mean instead of using geom_smooth/loess
  window <- max(1, floor(nrow(trace_df) / 50))
  trace_df$Smoothed <- stats::filter(trace_df$LogPosterior,
                                     rep(1/window, window), sides = 2)
  # Fill NA edges from the filter
  trace_df$Smoothed[is.na(trace_df$Smoothed)] <- trace_df$LogPosterior[is.na(trace_df$Smoothed)]
  
  p4 <- ggplot(trace_df, aes(x = Iteration)) +
    geom_line(aes(y = LogPosterior), colour = "#B2182B",
              alpha = 0.3, linewidth = 0.3) +
    geom_line(aes(y = Smoothed), colour = "#B2182B", linewidth = 0.8) +
    labs(title = "Log-Posterior Trace (convergence check)",
         subtitle = "Red curve = rolling mean  |  Should plateau after burn-in",
         x = "Iteration", y = "Log-Posterior") +
    theme_g25(base_size = 13)
  
  # Add burn-in line if available
  if (!is.null(run_info) && !is.null(run_info$burnin)) {
    p4 <- p4 +
      geom_vline(xintercept = run_info$burnin, linetype = "dashed",
                 colour = "#2166AC", linewidth = 0.8) +
      annotate("text", x = run_info$burnin,
               y = max(trace_df$LogPosterior, na.rm = TRUE),
               label = sprintf("  burn-in = %d", run_info$burnin),
               hjust = 0, vjust = 1, size = 3.5, colour = "#2166AC",
               fontface = "bold")
  }
  
  print(p4)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 5: Component trace plots (individual panels via patchwork)
  # ════════════════════════════════════════════════════════════════════════
  
  n_trace <- min(n_active, 8)
  trace_sources <- active[1:n_trace]
  
  # Downsample chain for trace plots (max 2000 samples per source)
  n_chain <- nrow(chain)
  max_trace_pts <- 2000
  if (n_chain > max_trace_pts) {
    trace_idx <- unique(round(seq(1, n_chain, length.out = max_trace_pts)))
  } else {
    trace_idx <- 1:n_chain
  }
  
  # Build long-form data for the chain traces
  chain_trace_df <- data.frame(
    Sample = rep(trace_idx, times = n_trace),
    Source = rep(trace_sources, each = length(trace_idx)),
    Proportion = unlist(lapply(trace_sources, function(s) chain[trace_idx, s] * 100))
  )
  chain_trace_df$Source <- factor(chain_trace_df$Source, levels = trace_sources)
  
  # Compute means for horizontal reference lines
  trace_means <- data.frame(
    Source = factor(trace_sources, levels = trace_sources),
    mean_pct = sapply(trace_sources, function(s) mean(chain[, s] * 100))
  )
  
  p5 <- ggplot(chain_trace_df, aes(x = Sample, y = Proportion)) +
    geom_line(aes(colour = Source), alpha = 0.35, linewidth = 0.25,
              show.legend = FALSE) +
    geom_hline(data = trace_means, aes(yintercept = mean_pct),
               linetype = "dashed", colour = "#2d3436", linewidth = 0.5) +
    facet_wrap(~ Source, ncol = 1, scales = "free_y",
               strip.position = "left") +
    scale_colour_manual(values = src_colours) +
    labs(title = "Component Trace Plots (post burn-in)",
         subtitle = "Well-mixed chains look like stationary noise; trends indicate poor convergence",
         x = "MCMC sample index", y = NULL) +
    theme_g25(base_size = 10) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold",
                                           size = 9, hjust = 1),
          panel.spacing = unit(0.3, "lines"),
          axis.text.y = element_text(size = 7))
  
  print(p5)
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 6: Pairwise correlation heatmap (NEW — not in original)
  # ════════════════════════════════════════════════════════════════════════
  
  if (n_active >= 2) {
    cor_mat <- cor(chain[, active, drop = FALSE])
    cor_melt <- reshape2::melt(cor_mat, varnames = c("Source1", "Source2"),
                               value.name = "Correlation")
    cor_melt$Source1 <- factor(cor_melt$Source1, levels = active)
    cor_melt$Source2 <- factor(cor_melt$Source2, levels = rev(active))
    
    p6 <- ggplot(cor_melt, aes(x = Source1, y = Source2, fill = Correlation)) +
      geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.2f", Correlation)),
                size = ifelse(n_active <= 6, 4.5, 3), fontface = "bold",
                colour = ifelse(abs(cor_melt$Correlation) > 0.6, "white", "#2d3436")) +
      scale_fill_gradient2(low = "#2166AC", mid = "#f5f5f5", high = "#B2182B",
                           midpoint = 0, limits = c(-1, 1),
                           name = "Posterior\nCorrelation") +
      labs(title = "Posterior Correlation Matrix",
           subtitle = paste(
             "Strong negative correlations indicate sources the model cannot distinguish \u2014",
             short_label),
           x = NULL, y = NULL) +
      theme_g25(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(face = "bold"),
            panel.grid = element_blank(),
            legend.position = "right")
    
    print(p6)
  }
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 7: Convergence diagnostics summary (ggplot table-style)
  # ════════════════════════════════════════════════════════════════════════
  
  if (!is.null(diag_data) && length(diag_data$active_names) > 0) {
    diag_df <- data.frame(
      Source = diag_data$active_names,
      ESS    = sapply(diag_data$active_names, function(nm) diag_data$ess[[nm]]),
      Geweke = sapply(diag_data$active_names, function(nm)
        round(diag_data$geweke[[nm]], 3)),
      stringsAsFactors = FALSE
    )
    diag_df$ESS_Status <- ifelse(diag_df$ESS >= 200, "Good",
                                 ifelse(diag_df$ESS >= 100, "Acceptable", "Low"))
    diag_df$Geweke_Status <- ifelse(abs(diag_df$Geweke) < 2, "OK", "CHECK")
    diag_df$Source <- factor(diag_df$Source,
                             levels = rev(diag_df$Source))
    
    # ESS bar chart
    p7a <- ggplot(diag_df, aes(x = ESS, y = Source, fill = ESS_Status)) +
      geom_col(width = 0.6, show.legend = TRUE) +
      geom_text(aes(label = ESS), hjust = -0.2, fontface = "bold", size = 3.8) +
      geom_vline(xintercept = 200, linetype = "dashed",
                 colour = "#4DAF4A", linewidth = 0.6) +
      geom_vline(xintercept = 100, linetype = "dotted",
                 colour = "#FF7F00", linewidth = 0.6) +
      scale_fill_manual(values = c("Good" = "#4DAF4A",
                                   "Acceptable" = "#FF7F00",
                                   "Low" = "#B2182B"),
                        name = "ESS Quality") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "Effective Sample Size (ESS)",
           subtitle = "Green line = 200 (good)  |  Orange line = 100 (acceptable minimum)",
           x = "ESS", y = NULL) +
      theme_g25(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(face = "bold", size = 10))
    
    # Geweke z-score dot plot
    diag_df_gew <- diag_df
    diag_df_gew$Source <- factor(diag_df_gew$Source,
                                 levels = levels(diag_df$Source))
    
    p7b <- ggplot(diag_df_gew, aes(x = Geweke, y = Source,
                                    colour = Geweke_Status)) +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed",
                 colour = "grey60", linewidth = 0.5) +
      geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
      annotate("rect", xmin = -2, xmax = 2,
               ymin = -Inf, ymax = Inf,
               fill = "#4DAF4A", alpha = 0.06) +
      geom_point(size = 4) +
      geom_text(aes(label = sprintf("z = %.2f", Geweke)),
                hjust = -0.15, size = 3.2, show.legend = FALSE) +
      scale_colour_manual(values = c("OK" = "#2166AC", "CHECK" = "#B2182B"),
                          name = "Convergence") +
      labs(title = "Geweke Convergence Diagnostic",
           subtitle = "|z| < 2 = converged (green zone)  |  |z| \u2265 2 = still drifting",
           x = "Geweke z-score", y = NULL) +
      theme_g25(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(face = "bold", size = 10))
    
    # Combine with patchwork
    p7_combined <- p7a / p7b +
      plot_annotation(
        title = "MCMC Convergence Diagnostics",
        subtitle = sprintf(
          "Chain length: %d samples  |  Acceptance rate: %.3f  (optimal: 0.20\u20130.50)",
          diag_data$chain_length,
          ifelse(!is.null(mcmc_result), mcmc_result$acceptance_rate, NA)),
        theme = theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(colour = "#636e72", size = 11, hjust = 0.5),
          plot.background = element_rect(fill = "white", colour = NA)
        )
      )
    
    print(p7_combined)
  }
  
  
  # ════════════════════════════════════════════════════════════════════════
  # PAGE 8: Processing log (text-based page via ggplot annotation)
  # ════════════════════════════════════════════════════════════════════════
  
  if (!is.null(run_info)) {
    # Build log text
    log_lines <- c(
      "================================================================",
      "  G25 Bayesian MCMC Admixture Estimator",
      "================================================================",
      "",
      "Parameters:",
      sprintf("  Source file:   %s", run_info$source_file),
      sprintf("  Target file:   %s", run_info$target_file),
      sprintf("  Mode:          %s", run_info$mode),
      sprintf("  Total sources: %d populations", run_info$n_sources),
      sprintf("  Iterations:    %d (burn-in: %d, thin: %d)",
              run_info$iter, run_info$burnin, run_info$thin),
      sprintf("  Sigma:         %.4f", run_info$sigma),
      sprintf("  Alpha:         %.2f", run_info$alpha),
      sprintf("  Seed:          %d", run_info$seed)
    )
    
    if (!is.null(run_info$max_k) && run_info$max_k > 0) {
      log_lines <- c(log_lines, sprintf("  Fixed K:       %d", run_info$max_k))
    }
    
    log_lines <- c(log_lines, "",
                   sprintf("Target: %s", run_info$target_label), "")
    
    # Pre-selection info
    log_lines <- c(log_lines,
                   "Step 1: Pre-selection",
                   sprintf("  Candidates kept: %d from %d total",
                           run_info$n_presel, run_info$n_sources))
    
    # Wrap candidate list
    cand_str <- run_info$candidates
    max_line_len <- 80
    while (nchar(cand_str) > 0) {
      if (nchar(cand_str) <= max_line_len) {
        log_lines <- c(log_lines, paste("   ", cand_str))
        cand_str <- ""
      } else {
        cut_pos <- max_line_len
        while (cut_pos > 1 && substr(cand_str, cut_pos, cut_pos) != ",") {
          cut_pos <- cut_pos - 1
        }
        if (cut_pos <= 1) cut_pos <- max_line_len
        log_lines <- c(log_lines, paste("   ", substr(cand_str, 1, cut_pos)))
        cand_str <- trimws(substr(cand_str, cut_pos + 1, nchar(cand_str)))
      }
    }
    
    log_lines <- c(log_lines, "",
                   "Step 2: Forward Stepwise Selection")
    if (!is.null(run_info$fwd_steps)) {
      for (step_line in run_info$fwd_steps) {
        log_lines <- c(log_lines, paste("  ", step_line))
      }
    }
    log_lines <- c(log_lines,
                   sprintf("  => Selected K = %d sources", run_info$k_selected),
                   "", "Results:")
    
    for (r in 1:nrow(summary_df)) {
      sig <- ifelse(summary_df$Significant[r], " ***", "")
      log_lines <- c(log_lines,
                     sprintf("  %-35s %5.1f%%  [%5.1f%% - %5.1f%%]%s",
                             summary_df$Source[r],
                             summary_df$Pct_Mean[r],
                             summary_df$Pct_CI_lo[r],
                             summary_df$Pct_CI_hi[r], sig))
    }
    
    log_lines <- c(log_lines, "",
                   sprintf("  Fit RMSE: %.6f", run_info$rmse))
    if (!is.null(mcmc_result)) {
      log_lines <- c(log_lines,
                     sprintf("  Acceptance rate: %.3f",
                             mcmc_result$acceptance_rate))
    }
    
    # Render as a ggplot text page
    log_text <- paste(log_lines, collapse = "\n")
    n_lines <- length(log_lines)
    
    p8 <- ggplot() +
      annotate("text", x = 0, y = 1, label = log_text,
               hjust = 0, vjust = 1, size = 2.8, family = "mono",
               colour = "#2d3436") +
      labs(title = "Processing Log") +
      xlim(-0.05, 1.2) + ylim(-0.1, 1.15) +
      theme_void() +
      theme(plot.title = element_text(face = "bold", size = 16,
                                      hjust = 0.5, margin = margin(b = 10)),
            plot.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(20, 30, 20, 30))
    
    print(p8)
  }
  
  
  dev.off()
  cat("  Plots written to:", outfile, "\n")
}


###############################################################################
# 9. MAIN
###############################################################################

args <- parse_arguments()
set.seed(args$seed)

cat("================================================================\n")
cat("  G25 Bayesian MCMC Admixture Estimator\n")
cat("  (Enhanced ggplot2 graphics edition)\n")
cat("================================================================\n\n")

# Load data
cat("Loading source file:", args$source, "\n")
src <- read_g25(args$source)
cat("  ->", nrow(src$coords), "samples across",
    length(unique(src$populations)), "populations\n")
cat("  -> Dimensions:", ncol(src$coords), "\n")

cat("Loading target file:", args$target, "\n")
tgt <- read_g25(args$target)
cat("  ->", nrow(tgt$coords), "target(s)\n\n")

# Aggregate if population mode
if (tolower(args$mode) == "population") {
  cat("Mode: POPULATION (averaging samples within each population)\n\n")
  pop_names <- unique(src$populations)
  pop_coords <- t(sapply(pop_names, function(p) {
    idx <- which(src$populations == p)
    if (length(idx) == 1) return(src$coords[idx, ])
    colMeans(src$coords[idx, , drop = FALSE])
  }))
  rownames(pop_coords) <- pop_names
  source_labels <- pop_names
  source_coords <- pop_coords
} else {
  cat("Mode: SAMPLE (each sample treated individually)\n\n")
  source_labels <- src$labels
  source_coords <- src$coords
  rownames(source_coords) <- source_labels
}

K_total <- nrow(source_coords)
D       <- ncol(source_coords)
cat("Source matrix:", K_total, "sources x", D, "dimensions\n\n")

# Process each target
cat("================================================================\n")
cat("  Beginning analysis\n")
cat("================================================================\n\n")

all_summaries <- list()

for (t_idx in 1:nrow(tgt$coords)) {
  
  target_vec   <- tgt$coords[t_idx, ]
  target_label <- tgt$labels[t_idx]
  
  cat("------------------------------------------------------------\n")
  cat("  Target:", target_label, "\n")
  cat("------------------------------------------------------------\n\n")
  
  # Pre-select sources
  cat("  Step 1: Pre-selecting candidate sources...\n")
  presel <- preselect_sources(target_vec, source_coords, source_labels,
                              n_keep = args$n_keep)
  sel_coords <- source_coords[presel$indices, , drop = FALSE]
  sel_names  <- source_labels[presel$indices]
  
  # Model selection
  if (args$max_k > 0) {
    cat("\n  Step 2: Forward selection with fixed K =", args$max_k, "\n")
    fwd <- forward_select(target_vec, sel_coords, sel_names,
                          sigma = args$sigma,
                          max_k = args$max_k, min_k = args$max_k)
    sel_coords <- sel_coords[fwd$indices, , drop = FALSE]
    sel_names  <- sel_names[fwd$indices]
  } else {
    cat("\n  Step 2: Auto model selection (forward stepwise + BIC)...\n")
    fwd <- forward_select(target_vec, sel_coords, sel_names,
                          sigma = args$sigma,
                          max_k = min(12, nrow(sel_coords)), min_k = 2)
    sel_coords <- sel_coords[fwd$indices, , drop = FALSE]
    sel_names  <- sel_names[fwd$indices]
  }
  
  # Run MCMC
  cat("\n  Step 3: Running Bayesian MCMC...\n")
  mcmc_result <- run_mcmc(
    target_vec   = target_vec,
    source_mat   = sel_coords,
    source_names = sel_names,
    n_iter       = args$iter,
    burnin       = args$burnin,
    thin         = args$thin,
    sigma        = args$sigma,
    alpha        = args$alpha
  )
  
  # Summarise
  cat("\n  Step 4: Summarising posteriors...\n")
  summary_df <- summarise_posterior(mcmc_result$chain, sel_names)
  all_summaries[[target_label]] <- summary_df
  
  # Print results
  cat("\n  +------------------------------------------------------------+\n")
  cat("  |  RESULTS:", target_label, "\n")
  cat("  +------------------------------------------------------------+\n")
  for (r in 1:nrow(summary_df)) {
    sig_marker <- ifelse(summary_df$Significant[r], " ***", "")
    cat(sprintf("  |  %-35s  %5.1f%%  [%5.1f%% - %5.1f%%]%s\n",
                summary_df$Source[r],
                summary_df$Pct_Mean[r],
                summary_df$Pct_CI_lo[r],
                summary_df$Pct_CI_hi[r],
                sig_marker))
  }
  cat("  +------------------------------------------------------------+\n")
  cat("    *** = 95% CI excludes zero (signal likely real)\n\n")
  
  # Fit quality
  w_mean <- colMeans(mcmc_result$chain)
  predicted <- as.vector(w_mean %*% sel_coords)
  residuals <- target_vec - predicted
  rmse <- sqrt(mean(residuals^2))
  cat("  Fit RMSE:", round(rmse, 6), "\n")
  cat("  (Lower = better; typical good fits are < 0.01)\n\n")
  
  # Save outputs
  out_prefix <- if (nrow(tgt$coords) > 1) {
    paste0(args$out, "_target", t_idx)
  } else {
    args$out
  }
  
  summary_file <- paste0(out_prefix, "_summary.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("  Summary saved:", summary_file, "\n")
  
  posterior_file <- paste0(out_prefix, "_posterior.csv")
  write.csv(as.data.frame(mcmc_result$chain), posterior_file, row.names = FALSE)
  cat("  Posterior chain saved:", posterior_file, "\n")
  
  diag_file <- paste0(out_prefix, "_diagnostics.txt")
  diag_data <- run_diagnostics(mcmc_result$chain, diag_file)
  
  # Build run_info for PDF log page
  run_info <- list(
    source_file  = args$source,
    target_file  = args$target,
    mode         = args$mode,
    n_sources    = K_total,
    iter         = args$iter,
    burnin       = args$burnin,
    thin         = args$thin,
    sigma        = args$sigma,
    alpha        = args$alpha,
    seed         = args$seed,
    max_k        = args$max_k,
    target_label = target_label,
    n_presel     = length(presel$indices),
    candidates   = paste(source_labels[presel$indices], collapse = ", "),
    fwd_steps    = fwd$step_log,
    k_selected   = length(sel_names),
    rmse         = rmse
  )
  
  plot_file <- paste0(out_prefix, "_plots.pdf")
  make_plots(mcmc_result$chain, summary_df, mcmc_result$log_posterior_trace,
             target_label, plot_file,
             diag_data = diag_data, run_info = run_info,
             mcmc_result = mcmc_result)
  
  cat("\n")
}

# Combined summary for multiple targets
if (length(all_summaries) > 1) {
  cat("================================================================\n")
  cat("  Combined summary for all targets\n")
  cat("================================================================\n\n")
  
  combined <- do.call(rbind, lapply(names(all_summaries), function(nm) {
    df <- all_summaries[[nm]]
    df$Target <- nm
    df
  }))
  combined_file <- paste0(args$out, "_all_targets_summary.csv")
  write.csv(combined, combined_file, row.names = FALSE)
  cat("  Combined results saved:", combined_file, "\n")
}

cat("\n================================================================\n")
cat("  Analysis complete!\n")
cat("================================================================\n")
