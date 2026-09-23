# 단계 1: Clot 2021 재현과 정량 gate (지시서 2026-09-23 §1, D-014)
skip_if_no_rxode2()
p  <- load_params("k2016")
dz <- read_cfg("design_clot2021.yaml")
days <- as.numeric(dz$sample_days_post_dose)
n_coh <- dz$cohort_sim$n_cohorts; n_per <- dz$n_per_cohort
sr <- dz$cohort_sim$sex_ratio_male$value
bm <- as.numeric(dz$cohort_sim$weight_bounds_kg$male); bf <- as.numeric(dz$cohort_sim$weight_bounds_kg$female)

test_that("(a) 용량군별 실제 체중으로 8명 코호트 2,000회: 관측 중앙값 tlast가 모의 코호트 중앙값의 5–95 백분위 안", {
  res <- rbindlist(lapply(dz$groups, function(g) {
    sim <- simulate_dataset(p, g$dose_mg, g$weight_mean, g$weight_sd, days, n_coh * n_per, paste0("clot_cohort_", g$dose_mg), sr, bm, bf)
    cm <- cohort_median_tlast(sim$nca, n_per)
    data.table(dose_mg = g$dose_mg, obs_median = g$tlast_median_obs,
               sim_median_of_medians = median(cm$tlast_median), p05 = q05(cm$tlast_median), p95 = q95(cm$tlast_median),
               subj_tlast_median = median(sim$nca$tlast_planned, na.rm = TRUE), subj_p05 = q05(sim$nca$tlast_planned), subj_p95 = q95(sim$nca$tlast_planned))
  }))
  message("\n코호트 중앙값 tlast 분포:\n", paste(capture.output(print(res)), collapse = "\n"))
  for (k in seq_len(nrow(res))) expect_true(res$obs_median[k] >= res$p05[k] & res$obs_median[k] <= res$p95[k],
                                            info = sprintf("%d mg: 관측 %g, 모의 5–95%% [%g, %g]", res$dose_mg[k], res$obs_median[k], res$p05[k], res$p95[k]))
})

test_that("(b) 정량 gate: 체중 매칭 후 모의 AUClast 평균이 관측 ±15% 이내, log-scale CV 35–51% (Cmax는 gate 제외)", {
  tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(dz$gate$log_cv_range_pct); nsim <- dz$gate$n_sim_per_dataset
  tab <- rbindlist(lapply(dz$datasets, function(ds) {
    sched <- dz$dataset_schedules[[ds$schedule]]; sched <- as.numeric(if (is.list(sched)) sched$days else sched)
    b <- ds_weight_bounds(ds)
    sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), sched, nsim, paste0("gate_", ds$id), sr, b$male, b$female)
    gate_row(ds, sim$nca, tol, cvr)
  }))
  message("\n정량 gate:\n", paste(capture.output(print(tab[, .(id, wt_mean, AUClast_obs_mean, AUClast_sim_mean = round(AUClast_sim_mean, 1), AUClast_ratio = round(AUClast_ratio, 3), AUClast_sim_logcv = round(AUClast_sim_logcv, 1), Cmax_ratio = round(Cmax_ratio, 3), pass_mean, pass_cv)])), collapse = "\n"))
  expect_true(all(tab$pass_mean), info = paste("평균 ±15% 실패:", paste(tab[pass_mean == FALSE, sprintf("%s (%.3f)", id, AUClast_ratio)], collapse = "; ")))
  expect_true(all(tab$pass_cv), info = paste("CV 35–51% 실패:", paste(tab[pass_cv == FALSE, sprintf("%s (%.1f%%)", id, AUClast_sim_logcv)], collapse = "; ")))
})

test_that("(c) 체중 기울기: 체중 5 kg 증가당 AUClast 약 50 mg·day/L 감소 (부호·크기 확인, 허용 ±50%는 확정 전 가정)", {
  wr <- as.numeric(dz$weight_slope_check$weight_range_kg)
  sim <- simulate_dataset(p, dz$weight_slope_check$dose_mg, mean(wr), 1e3, days, 20000, "wt_slope", sr, wr, wr)   # 큰 SD + 절단 = 거의 균등
  fit <- lm(AUClast ~ WT, data = sim$nca)
  slope5 <- unname(coef(fit)[2]) * 5
  message(sprintf("\n체중 기울기: %.1f mg·day/L per 5 kg (기대 약 -%g)", slope5, dz$weight_slope_check$expected_auclast_decrease_per_5kg))
  expect_lt(slope5, 0)
  expect_true(abs(slope5) >= 0.5 * dz$weight_slope_check$expected_auclast_decrease_per_5kg & abs(slope5) <= 1.5 * dz$weight_slope_check$expected_auclast_decrease_per_5kg)
})

test_that("(d) Cohen 2022: 200 mg, 70–100 kg, Day 43까지 채혈에서 log(AUClast) SD ≈ 0.49 (허용 ±0.10은 확정 전 가정)", {
  ds <- dz$datasets[[which(vapply(dz$datasets, function(d) d$id == "cohen2022_200_armA", logical(1)))]]
  sched <- as.numeric(dz$dataset_schedules$cohen_d42$days)
  b <- ds_weight_bounds(ds)
  sim <- simulate_dataset(p, 200, ds$weight_mean, ds_weight_sd(ds), sched, 20000, "cohen_sd", sr, b$male, b$female)
  sdlog <- sd(log(sim$nca$AUClast), na.rm = TRUE)
  message(sprintf("\nCohen 2022 비교: 모의 log(AUClast) SD = %.3f (관측 %.2f)", sdlog, dz$variability_checks$cohen2022$log_sd))
  expect_equal(sdlog, dz$variability_checks$cohen2022$log_sd, tolerance = 0.10 / dz$variability_checks$cohen2022$log_sd)
})
