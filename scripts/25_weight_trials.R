#!/usr/bin/env Rscript
# §6-2 체중 일반화(Pillar 2): 비만 비중이 높은 모집단(체중 N(100, 20) 절단 60–150 kg, 층은 절단 범위에서 자동 생성), arm당 117명, B0.
# 모델 (a) base, (b) struct2020, (d) k2016_bmi_vc0817. 시나리오 5개 × 2,000회. 탈락자(신뢰 기준 미충족) 체중 분포 비교 포함.
# 사용법: Rscript scripts/25_weight_trials.R <a|b|d> [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); mk <- args[1]; cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
MODELS <- c(a = "base", b = "struct2020", d = "k2016_bmi_vc0817")
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); rv <- resolve_variant(MODELS[[mk]], design, sc); p <- rv$p
wt <- list(mean = 100, sd = 20, trunc = c(60, 150), height = list(mean = 170, sd = 9, trunc = c(150, 195)))
MASTER_SEED <- 20260923L; N <- 2000L; SCEN <- c("S00", "F090", "KE120", "VM125", "KM10")
out_dir <- proj_path("results", "weight_generalization"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("weight_trials_", mk), master_seed = MASTER_SEED, run_mode = "final", extra = list(model = MODELS[[mk]], n = N, cores = cores))
res <- run_trials(N, p, design, sc$scenarios[SCEN], CJ(scenario = SCEN, schedule = "B0"), MASTER_SEED, wt, methods = "pooled_t", model_id = p$model_id, cores = cores, progress_every = 250)
fwrite(res$be, file.path(out_dir, sprintf("obese_trials_be_raw_%s.csv.gz", mk)))
st <- summarize_trials_ci(res$be)
drop <- res$ind[, .(n_trials = uniqueN(trial), n_reliable_mean = mean(n_reliable), n_dropout_mean = mean(n_dropout),
                    wt_retained_mean = mean(wt_retained_mean, na.rm = TRUE), wt_dropout_mean = mean(wt_dropout_mean, na.rm = TRUE),
                    dropout_pct = 100 * sum(n_dropout) / sum(n_dropout + n_reliable)), by = .(scenario, arm)]
fwrite(st$per_endpoint[, model := mk], file.path(out_dir, sprintf("obese_trials_per_endpoint_%s.csv", mk)))
fwrite(st$concordance[, model := mk], file.path(out_dir, sprintf("obese_trials_concordance_%s.csv", mk)))
fwrite(drop[, model := mk], file.path(out_dir, sprintf("obese_trials_dropout_%s.csv", mk)))
print(st$per_endpoint[endpoint %in% c("AUClast", "AUCinf_reliable", "AUCinf_true")]); print(st$concordance); print(drop)
append_run_log(logfile, "done")
