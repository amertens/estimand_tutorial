# summarise.R -- performance with Monte Carlo SEs, from stored rows and truths.
PAIRS <- data.frame(
  method   = c("itt", "astreated", "astreated", "astreated_w", "composite", "adherers"),
  strategy = c("treatment_policy", "while_on_treatment", "sustained_use", "sustained_use",
               "composite", "principal_stratum"))

summarise_runs <- function(rows, truth) {
  # a replicate that failed as a whole (run_scenario's .failed_replicate) has no
  # estimator identity: count it, then summarise the rest
  n_fail_rep <- length(unique(rows$rep[is.na(rows$family)]))
  rows <- rows[!is.na(rows$family), ]
  if (!all(rows$fingerprint == truth$fingerprint))
    stop("truth fingerprint does not match the rows' fingerprint")
  cox <- rows[rows$family == "cox", ]
  # a Cox suite that failed as a whole has no method: count it as a failure of every analysis
  whole <- cox[is.na(cox$method), ]
  if (nrow(whole)) {
    grid <- expand.grid(method = unique(PAIRS$method), window = c("all", "1y"), stringsAsFactors = FALSE)
    exp_rows <- do.call(rbind, lapply(seq_len(nrow(whole)), function(i) {
      w <- whole[rep(i, nrow(grid)), ]; w$method <- grid$method; w$window <- grid$window; w
    }))
    cox <- rbind(cox[!is.na(cox$method), ], exp_rows)
  }
  cox <- merge(cox[, setdiff(names(cox), "strategy")], PAIRS, by = "method")
  rows <- rbind(rows[rows$family != "cox", ], cox[, names(rows)])
  st <- truth$strategies
  key <- paste(rows$family, rows$strategy, rows$method, rows$window, rows$delta)
  do.call(rbind, lapply(split(rows, key), function(z) {
    ok <- z[z$ok %in% TRUE, ]
    red <- z$strategy[1] == "reduced_leaving"
    tr <- if (red) truth$reduced[truth$reduced$delta == z$delta[1], ] else st[st$strategy == z$strategy[1], ]
    true_hr <- if (red) NA_real_ else tr$hr_marginal
    n <- nrow(ok)
    if (z$family[1] == "cox") {
      le <- log(ok$est); bias <- mean(le) - log(tr$rr)
      cover <- mean(ok$lo <= tr$rr & tr$rr <= ok$hi); mean_est <- exp(mean(le)); emp <- sd(le)
      cover_hr <- mean(ok$lo <= true_hr & true_hr <= ok$hi)
      rmse <- sqrt(mean((le - log(tr$rr))^2))
      bias_logrr <- bias; mcse_logrr <- emp / sqrt(n); cover_rr <- cover; n_rr <- n
    } else {
      bias <- mean(ok$est) - tr$rd; mean_est <- mean(ok$est); emp <- sd(ok$est)
      cover <- mean(ok$lo <= tr$rd & tr$rd <= ok$hi); cover_hr <- NA_real_
      rmse <- sqrt(mean((ok$est - tr$rd)^2))
      # the RR metrics use rows with a positive risk ratio (an unstable fit can give a
      # negative risk); n_rr says how many
      pos <- is.finite(ok$rr) & ok$rr > 0
      lrr <- log(ok$rr[pos]); n_rr <- sum(pos)
      bias_logrr <- mean(lrr) - log(tr$rr); mcse_logrr <- sd(lrr) / sqrt(n_rr)
      cover_rr <- mean(ok$rr_lo[pos] <= tr$rr & tr$rr <= ok$rr_hi[pos])
    }
    mean_share <- if ("share" %in% names(ok) && any(!is.na(ok$share))) mean(ok$share, na.rm = TRUE) else NA_real_
    data.frame(scenario = z$scenario[1], family = z$family[1], strategy = z$strategy[1],
               method = z$method[1], window = z$window[1], delta = z$delta[1], n = n,
               n_fail = sum(!(z$ok %in% TRUE)), n_maxw100 = sum(ok$max_w > 100, na.rm = TRUE),
               true_rd = tr$rd, true_rr = tr$rr, true_hr = true_hr, mean_est = mean_est,
               bias = bias, mcse_bias = emp / sqrt(n), emp_se = emp, rmse = rmse,
               cover = cover, mcse_cover = sqrt(cover * (1 - cover) / n), cover_hr = cover_hr,
               bias_logrr = bias_logrr, mcse_logrr = mcse_logrr, cover_rr = cover_rr, n_rr = n_rr,
               mean_share = mean_share, n_fail_rep = n_fail_rep)
  }))
}

# Methods whose near-zero bias the paper claims (the main SDR estimator of each
# strategy and SDR in the true stratum); the pilot sizes n_sim on these. The
# untrimmed and flexible-library sensitivity analyses, the observed-adherer
# comparator and principal-score weighting are reported, not sized on.
sizes_equivalence <- function(method) method %in% c("SDR", "SDR, true stratum")
