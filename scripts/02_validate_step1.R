#!/usr/bin/env Rscript
# 단계 1 검증 표 산출 (지시서 §1). 테스트(tests/testthat/test-step1-clot2021.R)와 같은 계산을 표로 남긴다. 결과: results/step1/
# 사용법: Rscript scripts/02_validate_step1.R [k2016|k2020]  (k2020 = 사전 명시 구조 민감도 모델, 결과는 results/step1_k2020/)
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) args[1] else "k2016"
p <- load_params(MODEL); assert_params_confirmed(p)
dz <- read_cfg("design_clot2021.yaml"); design <- read_cfg("trial_design.yaml")
days <- as.numeric(dz$sample_days_post_dose)
out_dir <- proj_path("results", if (MODEL == "k2016") "step1" else paste0("step1_", MODEL)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("validate_step1_", MODEL), master_seed = "derive_seed(tag)", run_mode = "final")
sr <- dz$cohort_sim$sex_ratio_male$value
bm <- as.numeric(dz$cohort_sim$weight_bounds_kg$male); bf <- as.numeric(dz$cohort_sim$weight_bounds_kg$female)

# (a) 코호트 중앙값 tlast
coh <- rbindlist(lapply(dz$groups, function(g) {
  sim <- simulate_dataset(p, g$dose_mg, g$weight_mean, g$weight_sd, days, dz$cohort_sim$n_cohorts * dz$n_per_cohort, paste0("clot_cohort_", g$dose_mg), sr, bm, bf)
  cm <- cohort_median_tlast(sim$nca, dz$n_per_cohort)
  dist <- cm[, .N, by = tlast_median][order(tlast_median)][, pct := 100 * N / sum(N)]
  fwrite(dist, file.path(out_dir, sprintf("cohort_median_tlast_dist_%dmg.csv", g$dose_mg)))
  data.table(dose_mg = g$dose_mg, weight_mean = g$weight_mean, weight_sd = g$weight_sd, obs_median = g$tlast_median_obs,
             obs_range = if (is.null(g$tlast_range_obs)) NA_character_ else paste(g$tlast_range_obs, collapse = "-"),
             sim_median_of_cohort_medians = median(cm$tlast_median), p05 = q05(cm$tlast_median), p95 = q95(cm$tlast_median),
             pct_cohorts_at_obs = 100 * mean(cm$tlast_median == g$tlast_median_obs),
             cohort_min_median = median(cm$tlast_min), cohort_max_median = median(cm$tlast_max),
             subj_tlast_median = median(sim$nca$tlast_planned, na.rm = TRUE), subj_p05 = q05(sim$nca$tlast_planned), subj_p95 = q95(sim$nca$tlast_planned),
             AUClast_geo = geo_mean(sim$nca$AUClast), AUClast_mean = mean(sim$nca$AUClast, na.rm = TRUE), AUClast_logcv = log_cv_pct(sim$nca$AUClast),
             Cmax_mean = mean(sim$nca$Cmax, na.rm = TRUE), tmax_median = median(sim$nca$tmax, na.rm = TRUE),
             pass = g$tlast_median_obs >= q05(cm$tlast_median) & g$tlast_median_obs <= q95(cm$tlast_median))
}))
fwrite(coh, file.path(out_dir, "step1a_cohort_tlast.csv")); cat("\n(a) 코호트 중앙값 tlast:\n"); print(coh)

# (b) 정량 gate
tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(dz$gate$log_cv_range_pct); nsim <- dz$gate$n_sim_per_dataset
gate <- rbindlist(lapply(dz$datasets, function(ds) {
  sched <- dz$dataset_schedules[[ds$schedule]]; sched <- as.numeric(if (is.list(sched)) sched$days else sched)
  b <- ds_weight_bounds(ds)
  sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), sched, nsim, paste0("gate_", ds$id), sr, b$male, b$female)
  gate_row(ds, sim$nca, tol, cvr)
}))
fwrite(gate, file.path(out_dir, "step1b_quant_gate.csv")); cat("\n(b) 정량 gate:\n"); print(gate[, .(id, dose_mg, wt_mean, AUClast_obs_mean, AUClast_sim_mean = round(AUClast_sim_mean), AUClast_ratio = round(AUClast_ratio, 3), AUClast_sim_logcv = round(AUClast_sim_logcv, 1), Cmax_obs_mean, Cmax_sim_mean = round(Cmax_sim_mean, 1), Cmax_ratio = round(Cmax_ratio, 3), pass_mean, pass_cv)])

# (c) 체중 기울기와 체중 범위별 CV
wr <- as.numeric(dz$weight_slope_check$weight_range_kg)
wide <- simulate_dataset(p, 300, mean(wr), 1e3, days, 20000, "wt_slope", sr, wr, wr)
fit <- lm(AUClast ~ WT, data = wide$nca); slope5 <- unname(coef(fit)[2]) * 5
fit_log <- lm(log(AUClast) ~ WT, data = wide$nca)
slope_tab <- data.table(slope_per_5kg = slope5, expected = -dz$weight_slope_check$expected_auclast_decrease_per_5kg,
                        pct_change_per_5kg = 100 * (exp(unname(coef(fit_log)[2]) * 5) - 1), r2 = summary(fit)$r.squared)
fwrite(slope_tab, file.path(out_dir, "step1c_weight_slope.csv")); cat("\n(c) 체중 기울기:\n"); print(slope_tab)
# B0 일정(시험 설계)에서 체중 범위별 AUClast log-CV — 표본 수 논의용
wideB0 <- simulate_dataset(p, 300, mean(wr), 1e3, get_schedule(design, "B0"), 20000, "wt_cv_B0", sr, wr, wr)
ranges <- list(c(50, 115), c(60, 90), c(60, 75), c(75, 90), c(70, 100), c(55, 95), c(50, 100), c(65, 85))
cvtab <- rbindlist(lapply(ranges, function(r) { x <- wideB0$nca[WT >= r[1] & WT <= r[2]]
  data.table(wt_min = r[1], wt_max = r[2], n = nrow(x), AUClast_logcv = log_cv_pct(x$AUClast), AUClast_geo = geo_mean(x$AUClast), Cmax_logcv = log_cv_pct(x$Cmax),
             AUCinf_rel_logcv = log_cv_pct(x[reliable == TRUE, AUCinf]), reliable_pct = 100 * mean(x$reliable)) }))
fwrite(cvtab, file.path(out_dir, "step1e_auclast_cv_by_weight_range.csv")); cat("\n(e) 체중 범위별 AUClast log-CV (B0 일정):\n"); print(cvtab)

# (d) Cohen 2022 log SD
dsA <- dz$datasets[[which(vapply(dz$datasets, function(d) d$id == "cohen2022_200_armA", logical(1)))]]
b <- ds_weight_bounds(dsA)
coh200 <- simulate_dataset(p, 200, dsA$weight_mean, ds_weight_sd(dsA), as.numeric(dz$dataset_schedules$cohen_d42$days), 20000, "cohen_sd", sr, b$male, b$female)
cohen_tab <- data.table(sim_log_sd = sd(log(coh200$nca$AUClast), na.rm = TRUE), obs_log_sd = dz$variability_checks$cohen2022$log_sd,
                        sim_AUClast_mean = mean(coh200$nca$AUClast, na.rm = TRUE), obs_AUClast_mean_A = dsA$auclast_mean)
fwrite(cohen_tab, file.path(out_dir, "step1d_cohen2022_logsd.csv")); cat("\n(d) Cohen 2022:\n"); print(cohen_tab)

status <- data.table(
  check = c("(a) 코호트 중앙값 tlast 5–95%", "(b) AUClast 평균 ±15%", "(b) AUClast log-CV 35–51%", "(c) 체중 기울기 부호·크기", "(d) Cohen 2022 log SD ±0.10", "Cmax (gate 제외, 한계)"),
  status = c(if (all(coh$pass)) "PASS" else paste("FAIL:", paste(coh[pass == FALSE, dose_mg], collapse = ",")),
             if (all(gate$pass_mean)) "PASS" else paste("FAIL:", paste(gate[pass_mean == FALSE, sprintf("%s %.2f", id, AUClast_ratio)], collapse = "; ")),
             if (all(gate$pass_cv)) "PASS" else paste("FAIL:", paste(gate[pass_cv == FALSE, sprintf("%s %.1f", id, AUClast_sim_logcv)], collapse = "; ")),
             if (slope5 < 0 & abs(slope5) >= 25 & abs(slope5) <= 75) "PASS" else sprintf("CHECK: %.1f", slope5),
             if (abs(cohen_tab$sim_log_sd - cohen_tab$obs_log_sd) <= 0.10) "PASS" else sprintf("CHECK: %.3f", cohen_tab$sim_log_sd),
             sprintf("모의/관측 비 %s", paste(round(gate$Cmax_ratio, 2), collapse = ", "))))
fwrite(status, file.path(out_dir, "step1_status.csv")); cat("\n단계 1 상태:\n"); print(status)
append_run_log(logfile, "step1 done; all pass = ", all(coh$pass) & all(gate$pass_mean) & all(gate$pass_cv))
