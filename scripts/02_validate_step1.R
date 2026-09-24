#!/usr/bin/env Rscript
# 단계 1 검증 (지시서 §1, 검토 의견 3차 §1–2, 통합본 §0-3·§7·§8). gate_role = gate 항목만 pass/fail, external은 서술.
# 연구별 채혈 일정 사용. 모델별 내부/외부(개발 자료 포함 여부) 열과 완전 외부 데이터셋만의 gate 재판정을 함께 산출한다.
# 사용법: Rscript scripts/02_validate_step1.R [k2016|k2020]   결과: results/step1/ (k2016), results/step1_k2020/
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) args[1] else "k2016"
p <- load_params(MODEL); assert_params_confirmed(p)
dz <- read_cfg("design_clot2021.yaml"); design <- read_cfg("trial_design.yaml")
days <- schedule_days(dz, "clot")
out_dir <- proj_path("results", if (MODEL == "k2016") "step1" else paste0("step1_", MODEL)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("validate_step1_", MODEL), master_seed = "derive_seed(tag)", run_mode = "final")
sr <- dz$cohort_sim$sex_ratio_male$value
bm <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$male)); bf <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$female))
devcol <- paste0("dev_", MODEL)
tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(unlist(dz$gate$log_cv_range_pct)); nsim <- dz$gate$n_sim_per_dataset

# (a) 코호트 중앙값 tlast (Clot 일정) + 문헌 정합성 커버리지 비(§8)
cov_rows <- list()
coh <- rbindlist(lapply(dz$groups, function(g) {
  sim <- simulate_dataset(p, g$dose_mg, g$weight_mean, g$weight_sd, days, dz$cohort_sim$n_cohorts * dz$n_per_cohort, paste0("clot_cohort_", g$dose_mg), sr, bm, bf)
  cm <- cohort_median_tlast(sim$nca, dz$n_per_cohort)
  dist <- cm[, .N, by = tlast_median][order(tlast_median)][, pct := 100 * N / sum(N)]
  fwrite(dist, file.path(out_dir, sprintf("cohort_median_tlast_dist_%dmg.csv", g$dose_mg)))
  cov_rows[[as.character(g$dose_mg)]] <<- cbind(data.table(source = "Clot 2021 Table 3", dose_mg = g$dose_mg, weight_mean = g$weight_mean,
                                                            tlast_median_subj = median(sim$nca$tlast_planned, na.rm = TRUE), tlast_p05_subj = q05(sim$nca$tlast_planned), tlast_p95_subj = q95(sim$nca$tlast_planned),
                                                            tlast_cohort_median = median(cm$tlast_median), tlast_cohort_p05 = q05(cm$tlast_median), tlast_cohort_p95 = q95(cm$tlast_median)), coverage_ratios(sim$nca))
  in_rng <- g$tlast_median_obs >= q05(cm$tlast_median) & g$tlast_median_obs <= q95(cm$tlast_median)
  data.table(dose_mg = g$dose_mg, gate_role = g$gate_role, dev = "external", weight_mean = g$weight_mean, weight_sd = g$weight_sd, obs_median = g$tlast_median_obs,
             sim_median_of_cohort_medians = median(cm$tlast_median), p05 = q05(cm$tlast_median), p95 = q95(cm$tlast_median),
             pct_cohorts_at_obs = 100 * mean(cm$tlast_median == g$tlast_median_obs),
             subj_tlast_median = median(sim$nca$tlast_planned, na.rm = TRUE), subj_p05 = q05(sim$nca$tlast_planned), subj_p95 = q95(sim$nca$tlast_planned),
             AUClast_mean = mean(sim$nca$AUClast, na.rm = TRUE), AUClast_logcv = log_cv_pct(sim$nca$AUClast), Cmax_mean = mean(sim$nca$Cmax, na.rm = TRUE),
             tmax_median = median(sim$nca$tmax, na.rm = TRUE), obs_in_5_95 = in_rng, pass = if (g$gate_role == "gate") in_rng else NA)
}))
fwrite(coh, file.path(out_dir, "step1a_cohort_tlast.csv")); cat("\n(a) 코호트 중앙값 tlast:\n"); print(coh[, .(dose_mg, gate_role, obs_median, sim_median_of_cohort_medians, p05, p95, pct_cohorts_at_obs, pass)])

# (b) 정량 gate + 외부 점검 (연구별 일정)
gate <- rbindlist(lapply(dz$datasets, function(ds) {
  nca <- simulate_ds(p, ds, dz, nsim, paste0("gate_", ds$id), sr)
  cbind(gate_row(ds, nca, tol, cvr), data.table(dev = ds[[devcol]], schedule = if (!is.null(ds$components)) paste(vapply(ds$components, function(c) c$schedule, ""), collapse = "+") else ds$schedule))
}))
fwrite(gate, file.path(out_dir, "step1b_quant_gate.csv"))
cat("\n(b) 정량 gate + 외부 점검:\n"); print(gate[, .(id, gate_role, dev, schedule, AUClast_obs_mean, AUClast_sim_mean = round(AUClast_sim_mean), AUClast_ratio = round(AUClast_ratio, 3), AUClast_z = round(AUClast_z, 2), AUClast_sim_logcv = round(AUClast_sim_logcv, 1), Cmax_ratio = round(Cmax_ratio, 3), pass_mean, pass_cv)])

# 변경 전후 비교(이전 일정 결과가 있으면)
prev_f <- file.path(out_dir, if (MODEL == "k2016") "step1b_quant_gate_prev_schedule.csv" else "step1b_quant_gate_prev_schedule_iiv.csv")
if (file.exists(prev_f)) {
  prev <- fread(prev_f)
  prev[, id := sub("^pkm14271_200_nonasian$", "PKM14271_200mg_test", id)]
  cmp <- merge(prev[, .(id, AUClast_sim_prev = AUClast_sim_mean, ratio_prev = AUClast_ratio, cv_prev = AUClast_sim_logcv, Cmax_ratio_prev = Cmax_ratio)],
               gate[, .(id, schedule, AUClast_sim_new = AUClast_sim_mean, ratio_new = AUClast_ratio, cv_new = AUClast_sim_logcv, Cmax_ratio_new = Cmax_ratio, pass_mean)], by = "id")
  cmp[, ratio_change := ratio_new - ratio_prev]
  fwrite(cmp, file.path(out_dir, "step1b_schedule_change_comparison.csv")); cat("\n일정 변경 전후 비교:\n"); print(cmp)
}

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
fwrite(cvtab, file.path(out_dir, "step1e_auclast_cv_by_weight_range.csv"))

# (d) Cohen 2022 log SD — 외부 점검
cohen_ds <- Filter(function(d) d$id == "Cohen2022_200mg_AI", dz$datasets)[[1]]
cnca <- simulate_ds(p, cohen_ds, dz, nsim, paste0("gate_", cohen_ds$id), sr)
cohen_tab <- data.table(gate_role = "external", sim_log_sd = sd(log(cnca$AUClast), na.rm = TRUE), obs_log_sd = dz$variability_checks$cohen2022$log_sd)
fwrite(cohen_tab, file.path(out_dir, "step1d_cohen2022_logsd.csv"))

# (f) 300 mg 개별 arm (Li 2020 Table 3): 연구별 일정으로 모의, 체중은 확정 전 가정
ac <- dz$arm_checks_300mg; aw <- ac$weight
arm_sims <- list()
for (sch in unique(vapply(ac$arms, function(a) a$schedule, ""))) {
  sim <- simulate_dataset(p, 300, aw$mean, aw$sd, schedule_days(dz, sch), nsim, paste0("arm_", sch), sr, c(50, 100), c(40, 90))
  arm_sims[[sch]] <- sim$nca
}
arm_tab <- rbindlist(lapply(ac$arms, function(a) { x <- arm_sims[[a$schedule]]; m <- mean(x$AUClast, na.rm = TRUE)
  data.table(study = a$study, arm = a$arm, schedule = a$schedule, dev = a[[devcol]], auclast_obs = a$auclast_mean, auclast_sim = m, ratio = m / a$auclast_mean,
             sim_logcv = log_cv_pct(x$AUClast), tmax_obs = a$tmax_median, tmax_sim = median(x$tmax, na.rm = TRUE),
             pass_mean = abs(m / a$auclast_mean - 1) <= tol / 100, pass_cv = { cv <- log_cv_pct(x$AUClast); cv >= cvr[1] & cv <= cvr[2] }) }))
m300 <- gate[id == ac$model_reference_dataset, AUClast_sim_mean]
arm_sum <- data.table(model_mean_pooled = m300, arm_min = min(arm_tab$auclast_obs), arm_max = max(arm_tab$auclast_obs), in_range = m300 >= min(arm_tab$auclast_obs) & m300 <= max(arm_tab$auclast_obs))
fwrite(arm_tab, file.path(out_dir, "step1f_300mg_arm_check.csv")); fwrite(arm_sum, file.path(out_dir, "step1f_300mg_arm_check_summary.csv"))
cat("\n(f) 300 mg 개별 arm (연구별 일정):\n"); print(arm_tab); print(arm_sum)
# PKM12350 대조군 커버리지(§8): AUC0-t 500 대 AUC0-inf 521
cov_rows[["pkm12350"]] <- cbind(data.table(source = "FDA BLA 761055 Table 4.2.c (PKM12350)", dose_mg = 300, weight_mean = aw$mean,
                                           tlast_median_subj = median(arm_sims$pkm12350$tlast_planned, na.rm = TRUE), tlast_p05_subj = q05(arm_sims$pkm12350$tlast_planned), tlast_p95_subj = q95(arm_sims$pkm12350$tlast_planned),
                                           tlast_cohort_median = NA_real_, tlast_cohort_p05 = NA_real_, tlast_cohort_p95 = NA_real_), coverage_ratios(arm_sims$pkm12350))
fwrite(rbindlist(cov_rows, fill = TRUE), file.path(out_dir, "step1g_literature_coverage.csv"))

# 완전 외부 데이터셋만의 gate 재판정(§7)
eo <- dz$external_only_gate[[MODEL]]
ext_ds <- gate[id %in% unlist(eo$datasets), .(item = id, type = "dataset", ratio = AUClast_ratio, cv = AUClast_sim_logcv, pass_mean, pass_cv)]
ext_arm <- arm_tab[study %in% unlist(eo$arms), .(item = paste(study, arm), type = "Li 2020 arm", ratio, cv = sim_logcv, pass_mean, pass_cv)]
ext_coh <- coh[dose_mg %in% c(300, 600), .(item = paste0("Clot ", dose_mg, " mg 코호트 tlast"), type = "cohort tlast", ratio = NA_real_, cv = NA_real_, pass_mean = obs_in_5_95, pass_cv = NA)]
ext <- rbind(ext_ds, ext_arm, ext_coh)
fwrite(ext, file.path(out_dir, "step1h_external_only_gate.csv")); cat("\n완전 외부 데이터셋 gate:\n"); print(ext)

g_gate <- gate[gate_role == "gate"]; g_ext <- gate[gate_role == "external"]
all_gate <- all(coh[gate_role == "gate", pass]) && all(g_gate$pass_mean) && all(g_gate$pass_cv) && slope_tab$pass
ext_pass <- all(ext$pass_mean %in% TRUE) && all(ext[!is.na(pass_cv), pass_cv])
status <- data.table(
  role = c("gate", "gate", "gate", "gate", "gate", "gate(외부만)", "external", "external", "external", "external", "limitation"),
  check = c("(a) 코호트 중앙값 tlast 5–95% (300, 600 mg)", "(b) AUClast 평균 ±15% (연구 제형 5개)", "(b) AUClast log-CV 35–51% (연구 제형 5개)", "(c) 체중 기울기 부호·크기",
            "gate 종합", "완전 외부 데이터셋만의 gate", "(a) 200 mg 코호트 tlast", "(b) 200 mg 1.14 mL 제형 4개 arm: 모의/관측 AUClast (z)", "(d) Cohen 2022 log SD", "(f) 300 mg 개별 arm (Li 2020, 연구별 일정)", "Cmax 모의/관측"),
  status = c(if (all(coh[gate_role == "gate", pass])) "PASS" else "FAIL",
             if (all(g_gate$pass_mean)) "PASS" else paste("FAIL:", paste(g_gate[pass_mean == FALSE, sprintf("%s %.2f", id, AUClast_ratio)], collapse = "; ")),
             if (all(g_gate$pass_cv)) "PASS" else "FAIL",
             if (slope_tab$pass) "PASS" else sprintf("FAIL: %.1f", slope5),
             if (all_gate) "PASS" else "FAIL",
             sprintf("%s (%d개 항목: %s)", if (ext_pass) "PASS" else "FAIL", nrow(ext), paste(ext[, sprintf("%s %s", item, ifelse(is.na(ratio), ifelse(pass_mean, "범위 안", "범위 밖"), sprintf("%.2f", ratio)))], collapse = "; ")),
             sprintf("관측 %g일, 모의 5–95%% [%g, %g]", coh[dose_mg == 200, obs_median], coh[dose_mg == 200, p05], coh[dose_mg == 200, p95]),
             paste(g_ext[, sprintf("%s %.2f (z %.1f)", id, AUClast_ratio, AUClast_z)], collapse = "; "),
             sprintf("모의 %.3f, 관측 %.2f", cohen_tab$sim_log_sd, cohen_tab$obs_log_sd),
             paste(arm_tab[, sprintf("%s %s %.2f", study, arm, ratio)], collapse = "; "),
             paste(gate[, sprintf("%s %.2f", id, Cmax_ratio)], collapse = "; ")))
fwrite(status, file.path(out_dir, "step1_status.csv")); cat("\n단계 1 상태:\n"); print(status)
append_run_log(logfile, "step1 done; gate pass = ", all_gate, "; external-only = ", ext_pass)
