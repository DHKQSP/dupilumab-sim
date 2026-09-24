#!/usr/bin/env Rscript
# 기준 (d) 재평가: 개체 수준 200,000명, B0 vs D3 (검토 의견 통합본 §4-2). 20,000명 × 10묶음(묶음별 시드 태그), 대응 부트스트랩 구간.
# 최종 권고 B0는 결과와 무관하게 유지(§1). "취약" 표기만 갱신한다.
# 사용법: Rscript scripts/16_criterion_d_200k.R <variant>   (base | struct2020 | vmax080_both)
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); variant <- if (length(args)) args[1] else "base"
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); rv <- resolve_variant(variant, design, sc); p <- rv$p
MASTER_SEED <- 20260923L; N_CHUNK <- 20000L; K <- 10L
out_dir <- proj_path("results", "individual200k"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("criterion_d_200k_", variant), master_seed = MASTER_SEED, run_mode = "final", extra = list(variant = variant, n = N_CHUNK * K))
parts <- lapply(seq_len(K), function(k) {
  pop <- run_individual_population(N_CHUNK, p, design, c("B0", "D3"), MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = sprintf("d200k_%s_%02d", variant, k))
  cat(sprintf("  chunk %d/%d %s\n", k, K, format(Sys.time(), "%H:%M:%S")))
  pop$nca[, .(id = id + (k - 1L) * N_CHUNK, schedule, reliable, lambda_ok, pct_extrap, pct_extrap_true, err_AUClast_vs_true_inf, AUClast, WT)]
})
nca <- rbindlist(parts)
pi <- merge(paired_individual_vs_ref(nca, "B0"), paired_bootstrap_cd(nca, "B0", B = 4000), by = "schedule")
summ <- nca[, .(n = .N, reliable_pct = 100 * mean(reliable), extrap_gt20_pct = 100 * mean(pct_extrap > 20, na.rm = TRUE), n_extrap_gt20 = sum(pct_extrap > 20, na.rm = TRUE)), by = schedule]
res <- cbind(data.table(variant = variant, n_subjects = N_CHUNK * K), pi[, .(schedule, reliable_gain_pp, c_gain_boot_lo, c_gain_boot_hi, extrap_gt20_ratio, d_ratio_boot_lo, d_ratio_boot_hi, d_n_ref, d_n_sched, extrap_gt20_mcnemar_p)])
res[, `:=`(d_abs_change_per_arm = (d_n_sched - d_n_ref) / n_subjects * design$n_per_arm,
           d_meets_point = extrap_gt20_ratio <= design$decision_rule$extrap_gt20_ratio_max,
           d_fragile = d_ratio_boot_lo <= design$decision_rule$extrap_gt20_ratio_max & d_ratio_boot_hi >= design$decision_rule$extrap_gt20_ratio_max)]
fwrite(res, file.path(out_dir, sprintf("criterion_d_200k_%s.csv", variant))); fwrite(summ, file.path(out_dir, sprintf("summary_200k_%s.csv", variant)))
print(res); print(summ)
append_run_log(logfile, "done")
