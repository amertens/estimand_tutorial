# Debug concrete call — smaller N, simpler config
library(here); library(dplyr); library(data.table); library(survival); library(concrete)

source(here("DGP.R")); source(here("R", "helpers.R"))

tau <- 180
options(width = 200)

# Use N=2000 for faster debug
cat("Generating N=2000 dataset...\n"); flush.console()
set.seed(5001)
dat <- generate_hep_data(
  N = 2000, h0 = 5e-3, HR_early = 2.0, HR_late = 0.90,
  np_hazard = TRUE, complexity = FALSE, dep_censor = TRUE, switch_on = TRUE,
  lambda_sw0 = 1e-3, gamma_A = 0.80, gamma_ckd = 0.60, seed = 5001)

dat_wot <- dat %>%
  mutate(
    wot_time = pmin(event_time, switch_time, censor_admin),
    wot_status = as.integer(event_time <= switch_time &
                            event_time <= censor_admin)
  )

baseline <- c("age", "sex_male", "ckd", "cirrhosis", "nsaid",
              "diabetes", "hypertension", "heart_failure")
dt <- data.table::as.data.table(
  dat_wot[, c("wot_time", "wot_status", "treatment", baseline)]
)
setnames(dt, c("wot_time", "wot_status", "treatment"),
         c("time", "status", "trt"))

cat("Events:", sum(dt$status == 1), "/", nrow(dt), "\n"); flush.console()
cat("Max time:", max(dt$time), "\n"); flush.console()

cat("\nCall formatArguments...\n"); flush.console()
args <- concrete::formatArguments(
  DataTable = dt,
  EventTime = "time", EventType = "status", Treatment = "trt",
  TargetTime = tau, Intervention = 0:1
)
cat("formatArguments OK\n"); flush.console()

cat("\nDefault Model structure:\n"); flush.console()
str(args$Model, max.level = 3)
flush.console()

cat("\nCall doConcrete with default Model (no coxnet override)...\n"); flush.console()
est <- tryCatch(
  concrete::doConcrete(ConcreteArgs = args),
  error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(est)) {
  cat("\ndoConcrete OK. Calling getOutput...\n"); flush.console()
  out <- concrete::getOutput(est, Estimand = c("RD", "Risk"), Simultaneous = FALSE)
  cat("Columns:\n"); print(names(out))
  cat("\nUnique Estimands / Estimators / Events / Times:\n")
  cat(" Estimand:", paste(unique(out$Estimand), collapse=", "), "\n")
  cat(" Estimator:", paste(unique(out$Estimator), collapse=", "), "\n")
  cat(" Event:", paste(unique(out$Event), collapse=", "), "\n")
  cat(" Time:", paste(unique(out$Time), collapse=", "), "\n")
  cat("\nRD rows:\n")
  print(out[out$Estimand == "RD", ])
}
