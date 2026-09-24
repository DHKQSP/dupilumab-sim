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
# 변형별 조합(검토 의견 3차 §5):
#  base, struct2020: 6개 일정 × S00/F090/VM125/KE110
#  vmax080_both, vmax125_both(곡률 민감도): 6개 일정 × S00/KE110
#  나머지 변형: B0·D2 × S00/F090/KE110
all_scens <- sc$scenarios[sc$schedule_analysis]
if (variant %in% c("base", "struct2020")) { scheds <- design$schedule_analysis; scens <- all_scens
} else if (variant %in% c("vmax080_both", "vmax125_both")) { scheds <- design$schedule_analysis; scens <- all_scens[c("S00", "KE110")]
} else if (variant == "noresid") { scheds <- design$schedule_analysis; scens <- all_scens["S00"]      # §2 폭 확대 원인 분해
} else { scheds <- c("B0", "D2"); scens <- all_scens[c("S00", "F090", "KE110")] }
combos <- CJ(scenario = names(scens), schedule = scheds)
logfile <- start_run_log(paste0("trial_schedules_", variant), master_seed = MASTER_SEED, run_mode = "final",
                         extra = list(variant = variant, n_trials = n_trials, cores = cores, schedules = scheds, scenarios = names(scens)))
cat(sprintf("variant=%s n_trials=%d cores=%d combos=%d\n", variant, n_trials, cores, nrow(combos)))
t0 <- Sys.time()
res <- run_trials(n_trials, p, design, scens, combos, MASTER_SEED, rv$wt_spec, jitter = TRUE, methods = c("pooled_t", "ancova_weight"), model_id = p$model_id, cores = cores)
cat("elapsed:", format(Sys.time() - t0), "\n")
fwrite(res$be, file.path(out_dir, sprintf("schedules_be_raw_%s.csv.gz", variant)))
ind_mean <- res$ind[, lapply(.SD, mean, na.rm = TRUE), by = .(scenario, schedule, arm), .SDcols = is.numeric]
fwrite(ind_mean, file.path(out_dir, sprintf("schedules_trial_individual_means_%s.csv", variant)))
pp <- postprocess_schedules(variant, design)
append_run_log(logfile, "done")
