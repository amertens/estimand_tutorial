# truths_all.R -- per-scenario truths, tagged with the parameter fingerprint.
scenario_truth <- function(name, p, N = 1e6, chunk = 2.5e5,
                           deltas = if (name == "base") c(1, 0.75, 0.5, 0.25, 0.1, 0) else c(1, 0.5, 0)) {
  list(scenario = name, fingerprint = params_fingerprint(p),
       strategies = rwe_truth(p, N = N, chunk = chunk),
       reduced = rwe_truth_reduced(p, deltas, N = N, chunk = chunk))
}

get_truth <- function(name) {
  tr <- readRDS(rwe_path(sprintf("sandbox_rwe/results/truth_%s.rds", name)))
  if (!identical(tr$fingerprint, params_fingerprint(get_scenario_params(name))))
    stop("truth for ", name, " was computed with different parameters; rerun run/10_truths.R")
  tr
}
