#!/usr/bin/env Rscript
# 시험 수준: 제품 차이 시나리오 14개 × 500회, B0 일정 (지시서 §6). 결과: results/trials/products_*.csv
# 사용법: Rscript scripts/12_trial_products.R [n_trials] [cores] [variant]  (variant: base | struct2020)
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
n_trials <- if (length(args) >= 1) as.integer(args[1]) else design$mc$n_trials
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores())
MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "trials"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
variant <- if (length(args) >= 3) args[3] else "base"
rv <- resolve_variant(variant, design, sc); p <- rv$p
sfx <- if (variant == "base") "" else paste0("_", variant)
combos <- CJ(scenario = names(sc$scenarios), schedule = "B0")
logfile <- start_run_log(paste0("trial_products_", variant), master_seed = MASTER_SEED, run_mode = "final", extra = list(n_trials = n_trials, cores = cores))
t0 <- Sys.time()
res <- run_trials(n_trials, p, design, sc$scenarios, combos, MASTER_SEED, rv$wt_spec, jitter = TRUE, methods = c("pooled_t", "ancova_weight"), model_id = p$model_id, cores = cores)
cat("elapsed:", format(Sys.time() - t0), "\n")
fwrite(res$be, file.path(out_dir, paste0("products_be_raw", sfx, ".csv.gz")))
fwrite(res$ind[, lapply(.SD, mean, na.rm = TRUE), by = .(scenario, schedule, arm), .SDcols = is.numeric], file.path(out_dir, paste0("products_trial_individual_means", sfx, ".csv")))
pp <- postprocess_products(variant)
append_run_log(logfile, "done")
