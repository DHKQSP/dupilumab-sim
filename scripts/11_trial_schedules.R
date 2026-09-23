#!/usr/bin/env Rscript
# 시험 수준: 일정 × 시나리오(S00, F090, VM125) × 500회, arm당 117명 (지시서 §5). 판정 규칙과 한계 가치 표. 결과: results/trials/
# 사용법: Rscript scripts/11_trial_schedules.R [variant] [n_trials] [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "base"
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
n_trials <- if (length(args) >= 2) as.integer(args[2]) else design$mc$n_trials
cores <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 0L)
MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "trials"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rv <- resolve_variant(variant, design, sc); p <- rv$p
scheds <- if (variant == "base") design$schedule_analysis else c("B0", "D2")
scens <- sc$scenarios[sc$schedule_analysis]
if (variant != "base") scens <- scens[c("S00", "F090")]
combos <- CJ(scenario = names(scens), schedule = scheds)
logfile <- start_run_log(paste0("trial_schedules_", variant), master_seed = MASTER_SEED, run_mode = "final",
                         extra = list(variant = variant, n_trials = n_trials, cores = cores, schedules = scheds, scenarios = names(scens)))
cat(sprintf("variant=%s n_trials=%d cores=%d combos=%d\n", variant, n_trials, cores, nrow(combos)))
t0 <- Sys.time()
res <- run_trials(n_trials, p, design, scens, combos, MASTER_SEED, rv$wt_spec, jitter = TRUE, methods = c("pooled_t", "ancova_weight"), model_id = p$model_id, cores = cores)
cat("elapsed:", format(Sys.time() - t0), "\n")
fwrite(res$be, file.path(out_dir, sprintf("schedules_be_raw_%s.csv", variant)))
st <- summarize_trials(res$be, method = "pooled_t")
st_anc <- summarize_trials(res$be, method = "ancova_weight")
fwrite(st$per_endpoint, file.path(out_dir, sprintf("schedules_per_endpoint_%s.csv", variant)))
fwrite(st$concordance, file.path(out_dir, sprintf("schedules_concordance_%s.csv", variant)))
fwrite(st_anc$per_endpoint, file.path(out_dir, sprintf("schedules_per_endpoint_ancova_%s.csv", variant)))
ind_mean <- res$ind[, lapply(.SD, mean, na.rm = TRUE), by = .(scenario, schedule, arm), .SDcols = is.numeric]
fwrite(ind_mean, file.path(out_dir, sprintf("schedules_trial_individual_means_%s.csv", variant)))
print(st$per_endpoint[endpoint %in% c("AUClast", "Cmax", "AUCinf_reliable"), .(scenario, schedule, endpoint, pass_rate = round(pass_rate, 1), GMR_mean = round(GMR_mean, 4), width_mean_pp = round(width_mean_pp, 2), n_T_mean = round(n_T_mean, 1))])
print(st$concordance[, .(scenario, schedule, pass_both_primary, agree_last_infrel, last_pass_infrel_fail, cor_logGMR_last_infrel = round(cor_logGMR_last_infrel, 4))])
# 판정 규칙(개인 수준 결과 필요)
indf <- proj_path("results", "individual", sprintf("individual_%s.csv", variant))
if (file.exists(indf)) {
  ind <- fread(indf)
  dec <- schedule_decision(ind, st$per_endpoint, design, ref = "B0", scenario = "S00")
  fwrite(dec, file.path(out_dir, sprintf("schedule_decision_%s.csv", variant)))
  cat("\n판정 규칙(B0 대비):\n"); print(dec)
}
append_run_log(logfile, "done")
