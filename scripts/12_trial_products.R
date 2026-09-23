#!/usr/bin/env Rscript
# 시험 수준: 제품 차이 시나리오 14개 × 500회, B0 일정 (지시서 §6). 결과: results/trials/products_*.csv
# 사용법: Rscript scripts/12_trial_products.R [n_trials] [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
n_trials <- if (length(args) >= 1) as.integer(args[1]) else design$mc$n_trials
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores())
MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "trials"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rv <- resolve_variant("base", design, sc); p <- rv$p
combos <- CJ(scenario = names(sc$scenarios), schedule = "B0")
logfile <- start_run_log("trial_products", master_seed = MASTER_SEED, run_mode = "final", extra = list(n_trials = n_trials, cores = cores))
t0 <- Sys.time()
res <- run_trials(n_trials, p, design, sc$scenarios, combos, MASTER_SEED, rv$wt_spec, jitter = TRUE, methods = c("pooled_t", "ancova_weight"), model_id = p$model_id, cores = cores)
cat("elapsed:", format(Sys.time() - t0), "\n")
fwrite(res$be, file.path(out_dir, "products_be_raw.csv"))
st <- summarize_trials(res$be, method = "pooled_t"); st_anc <- summarize_trials(res$be, method = "ancova_weight")
fwrite(st$per_endpoint, file.path(out_dir, "products_per_endpoint.csv")); fwrite(st$concordance, file.path(out_dir, "products_concordance.csv"))
fwrite(st_anc$per_endpoint, file.path(out_dir, "products_per_endpoint_ancova.csv"))
fwrite(res$ind[, lapply(.SD, mean, na.rm = TRUE), by = .(scenario, schedule, arm), .SDcols = is.numeric], file.path(out_dir, "products_trial_individual_means.csv"))
print(dcast(st$per_endpoint[, .(scenario, endpoint, pass_rate = round(pass_rate, 1))], scenario ~ endpoint, value.var = "pass_rate"))
print(dcast(st$per_endpoint[, .(scenario, endpoint, GMR = round(GMR_mean, 4))], scenario ~ endpoint, value.var = "GMR"))
print(st$concordance[, .(scenario, agree_last_infrel, last_pass_infrel_fail, last_fail_infrel_pass, cor_logGMR_last_infrel = round(cor_logGMR_last_infrel, 3), cor_logGMR_last_true = round(cor_logGMR_last_true, 3))])
append_run_log(logfile, "done")
