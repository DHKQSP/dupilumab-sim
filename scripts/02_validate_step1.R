#!/usr/bin/env Rscript
# 단계 1 검증 표 (지시서 §1, 검토 의견 3차 §1–2). gate_role = gate 항목만 pass/fail, external 항목은 서술 보고.
# 사용법: Rscript scripts/02_validate_step1.R [k2016|k2020]   결과: results/step1/ (k2016), results/step1_k2020/
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) args[1] else "k2016"
p <- load_params(MODEL); assert_params_confirmed(p)
dz <- read_cfg("design_clot2021.yaml"); design <- read_cfg("trial_design.yaml")
days <- as.numeric(unlist(dz$sample_days_post_dose))
out_dir <- proj_path("results", if (MODEL == "k2016") "step1" else paste0("step1_", MODEL)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("validate_step1_", MODEL), master_seed = "derive_seed(tag)", run_mode = "final")
sr <- dz$cohort_sim$sex_ratio_male$value
bm <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$male)); bf <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$female))
sched_of <- function(ds) { s <- dz$dataset_schedules[[ds$schedule]]; as.numeric(unlist(if (is.list(s) && !is.null(s$days)) s$days else s)) }

# (a) 코호트 중앙값 tlast
coh <- rbindlist(lapply(dz$groups, function(g) {
  sim <- simulate_dataset(p, g$dose_mg, g$weight_mean, g$weight_sd, days, dz$cohort_sim$n_cohorts * dz$n_per_cohort, paste0("clot_cohort_", g$dose_mg), sr, bm, bf)
  cm <- cohort_median_tlast(sim$nca, dz$n_per_cohort)
  dist <- cm[, .N, by = tlast_median][order(tlast_median)][, pct := 100 * N / sum(N)]
  fwrite(dist, file.path(out_dir, sprintf("cohort_median_tlast_dist_%dmg.csv", g$dose_mg)))
  in_rng <- g$tlast_median_obs >= q05(cm$tlast_median) & g$tlast_median_obs <= q95(cm$tlast_median)
  data.table(dose_mg = g$dose_mg, gate_role = g$gate_role, weight_mean = g$weight_mean, weight_sd = g$weight_sd, obs_median = g$tlast_median_obs,
             sim_median_of_cohort_medians = median(cm$tlast_median), p05 = q05(cm$tlast_median), p95 = q95(cm$tlast_median),
             pct_cohorts_at_obs = 100 * mean(cm$tlast_median == g$tlast_median_obs),
             subj_tlast_median = median(sim$nca$tlast_planned, na.rm = TRUE), subj_p05 = q05(sim$nca$tlast_planned), subj_p95 = q95(sim$nca$tlast_planned),
             AUClast_mean = mean(sim$nca$AUClast, na.rm = TRUE), AUClast_logcv = log_cv_pct(sim$nca$AUClast), Cmax_mean = mean(sim$nca$Cmax, na.rm = TRUE),
             tmax_median = median(sim$nca$tmax, na.rm = TRUE), obs_in_5_95 = in_rng, pass = if (g$gate_role == "gate") in_rng else NA)
}))
fwrite(coh, file.path(out_dir, "step1a_cohort_tlast.csv")); cat("\n(a) 코호트 중앙값 tlast:\n"); print(coh[, .(dose_mg, gate_role, obs_median, sim_median_of_cohort_medians, p05, p95, pct_cohorts_at_obs, obs_in_5_95, pass)])

# (b) 정량 gate + 외부 점검
tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(unlist(dz$gate$log_cv_range_pct)); nsim <- dz$gate$n_sim_per_dataset
sims <- list()
gate <- rbindlist(lapply(dz$datasets, function(ds) {
  b <- ds_weight_bounds(ds)
  sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), sched_of(ds), nsim, paste0("gate_", ds$id), sr, b$male, b$female)
  sims[[ds$id]] <<- sim$nca
  gate_row(ds, sim$nca, tol, cvr)
}))
fwrite(gate, file.path(out_dir, "step1b_quant_gate.csv"))
cat("\n(b) 정량 gate(연구 제형) + 외부 점검(200 mg 1.14 mL of 175 mg/mL):\n"); print(gate[, .(id, gate_role, n_obs, AUClast_obs_mean, AUClast_sim_mean = round(AUClast_sim_mean), AUClast_ratio = round(AUClast_ratio, 3), AUClast_z = round(AUClast_z, 2), AUClast_sim_logcv = round(AUClast_sim_logcv, 1), tmax_obs_median, tmax_sim_median, Cmax_ratio = round(Cmax_ratio, 3), pass_mean, pass_cv)])

# (c) 체중 기울기, (e) 체중 범위별 CV
wr <- as.numeric(unlist(dz$weight_slope_check$weight_range_kg))
wide <- simulate_dataset(p, 300, mean(wr), 1e3, days, 20000, "wt_slope", sr, wr, wr)
fit <- lm(AUClast ~ WT, data = wide$nca); slope5 <- unname(coef(fit)[2]) * 5; fit_log <- lm(log(AUClast) ~ WT, data = wide$nca)
slope_tab <- data.table(slope_per_5kg = slope5, expected = -dz$weight_slope_check$expected_auclast_decrease_per_5kg, pct_change_per_5kg = 100 * (exp(unname(coef(fit_log)[2]) * 5) - 1), r2 = summary(fit)$r.squared,
                        pass = slope5 < 0 & abs(slope5) >= 25 & abs(slope5) <= 75)
fwrite(slope_tab, file.path(out_dir, "step1c_weight_slope.csv")); cat("\n(c) 체중 기울기:\n"); print(slope_tab)
wideB0 <- simulate_dataset(p, 300, mean(wr), 1e3, get_schedule(design, "B0"), 20000, "wt_cv_B0", sr, wr, wr)
ranges <- list(c(50, 115), c(60, 90), c(60, 75), c(75, 90), c(70, 100), c(55, 95), c(50, 100), c(65, 85))
cvtab <- rbindlist(lapply(ranges, function(r) { x <- wideB0$nca[WT >= r[1] & WT <= r[2]]
  data.table(wt_min = r[1], wt_max = r[2], n = nrow(x), AUClast_logcv = log_cv_pct(x$AUClast), AUClast_geo = geo_mean(x$AUClast), Cmax_logcv = log_cv_pct(x$Cmax),
             AUCinf_rel_logcv = log_cv_pct(x[reliable == TRUE, AUCinf]), reliable_pct = 100 * mean(x$reliable)) }))
fwrite(cvtab, file.path(out_dir, "step1e_auclast_cv_by_weight_range.csv")); cat("\n(e) 체중 범위별 AUClast log-CV (B0):\n"); print(cvtab)

# (d) Cohen 2022 log SD — 외부 점검(서술)
coh_ids <- c("Cohen2022_200mg_AI", "Cohen2022_200mg_PFSS")
cohen_tab <- data.table(gate_role = "external", sim_log_sd = sd(log(sims[["Cohen2022_200mg_AI"]]$AUClast), na.rm = TRUE), obs_log_sd = dz$variability_checks$cohen2022$log_sd)
fwrite(cohen_tab, file.path(out_dir, "step1d_cohen2022_logsd.csv")); cat("\n(d) Cohen 2022 log SD (외부 점검):\n"); print(cohen_tab)

# (f) 300 mg 개별 arm 외부 점검 (Li 2020 Table 3 p.750)
ac <- dz$arm_checks_300mg
arms <- rbindlist(lapply(ac$arms, as.data.table))
m300 <- gate[id == ac$model_reference_dataset, AUClast_sim_mean]
arm_tab <- cbind(arms, data.table(model_mean_300mg_78kg = m300, ratio_model_to_arm = m300 / arms$auclast_mean))
arm_sum <- data.table(model_mean = m300, arm_min = min(arms$auclast_mean), arm_max = max(arms$auclast_mean), in_range = m300 >= min(arms$auclast_mean) & m300 <= max(arms$auclast_mean))
fwrite(arm_tab, file.path(out_dir, "step1f_300mg_arm_check.csv")); fwrite(arm_sum, file.path(out_dir, "step1f_300mg_arm_check_summary.csv"))
cat("\n(f) 300 mg 개별 arm 외부 점검:\n"); print(arm_tab); print(arm_sum)

g_gate <- gate[gate_role == "gate"]; g_ext <- gate[gate_role == "external"]
status <- data.table(
  role = c("gate", "gate", "gate", "gate", "gate", "external", "external", "external", "external", "limitation"),
  check = c("(a) 코호트 중앙값 tlast 5–95% (300, 600 mg)", "(b) AUClast 평균 ±15% (연구 제형 5개)", "(b) AUClast log-CV 35–51% (연구 제형 5개)", "(c) 체중 기울기 부호·크기",
            "gate 종합", "(a) 200 mg 코호트 tlast", "(b) 200 mg 1.14 mL 제형 4개 arm: 모의/관측 AUClast (z)", "(d) Cohen 2022 log SD", "(f) 300 mg 개별 arm 범위(Li 2020)", "Cmax 모의/관측"),
  status = c(if (all(coh[gate_role == "gate", pass])) "PASS" else "FAIL",
             if (all(g_gate$pass_mean)) "PASS" else paste("FAIL:", paste(g_gate[pass_mean == FALSE, sprintf("%s %.2f", id, AUClast_ratio)], collapse = "; ")),
             if (all(g_gate$pass_cv)) "PASS" else "FAIL",
             if (slope_tab$pass) "PASS" else sprintf("FAIL: %.1f", slope5),
             if (all(coh[gate_role == "gate", pass]) && all(g_gate$pass_mean) && all(g_gate$pass_cv) && slope_tab$pass) "PASS" else "FAIL",
             sprintf("관측 %g일, 모의 5–95%% [%g, %g]", coh[dose_mg == 200, obs_median], coh[dose_mg == 200, p05], coh[dose_mg == 200, p95]),
             paste(g_ext[, sprintf("%s %.2f (z %.1f)", id, AUClast_ratio, AUClast_z)], collapse = "; "),
             sprintf("모의 %.3f, 관측 %.2f", cohen_tab$sim_log_sd, cohen_tab$obs_log_sd),
             sprintf("모델 %.0f, arm 범위 [%.1f, %.0f] → %s", m300, arm_sum$arm_min, arm_sum$arm_max, if (arm_sum$in_range) "범위 안" else "범위 밖"),
             paste(gate[, sprintf("%s %.2f", id, Cmax_ratio)], collapse = "; ")))
fwrite(status, file.path(out_dir, "step1_status.csv")); cat("\n단계 1 상태:\n"); print(status)
append_run_log(logfile, "step1 done; gate pass = ", status[check == "gate 종합", status])
