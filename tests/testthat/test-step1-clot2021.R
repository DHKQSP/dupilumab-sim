# 단계 1 gate (지시서 §1, 검토 의견 3차 §1–2, D-014/D-023). gate_role = gate 항목만 pass/fail.
# 200 mg 1.14 mL(175 mg/mL) 제형 데이터셋과 200 mg 코호트는 외부 점검 — 판정하지 않고 값만 계산(스크립트 02가 표로 보고).
# (a)(b)(c)(f)는 판정성 테스트(통과 여부가 모의 결과에 달린 과학적 판정, 느림: results/ci/test_timing_local.csv): fast 범위(CI)에서는 건너뛰고,
# full 범위에서는 실행하되 기대값 불충족은 보고 항목으로 분류한다(helper-setup.R skip_if_judgment, scripts/run_tests.R).
# 정식 판정 기록은 scripts/02_validate_step1.R의 results/step1/ 결과 파일이다. gate 범위 검사(첫 테스트)는 config 검사라 fast에 남긴다.
skip_if_no_rxode2()
p  <- load_params("k2016")
dz <- read_cfg("design_clot2021.yaml")
days <- as.numeric(unlist(dz$sample_days_post_dose))
n_coh <- dz$cohort_sim$n_cohorts; n_per <- dz$n_per_cohort; sr <- dz$cohort_sim$sex_ratio_male$value
bm <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$male)); bf <- as.numeric(unlist(dz$cohort_sim$weight_bounds_kg$female))
sched_of <- function(ds) { s <- dz$dataset_schedules[[ds$schedule]]; as.numeric(unlist(if (is.list(s) && !is.null(s$days)) s$days else s)) }

test_that("gate 범위: 연구 제형(300 mg, 600 mg)만 gate, 200 mg 1.14 mL 제형은 external (gate_decision.yaml)", {
  g <- read_cfg("gate_decision.yaml")
  expect_equal(g$decision_option, 1); expect_equal(g$primary_model, "kovalenko2016_blq")
  roles <- vapply(dz$datasets, function(d) d$gate_role, ""); doses <- vapply(dz$datasets, function(d) d$dose_mg, 0)
  expect_true(all(roles[doses == 200] == "external")); expect_true(all(roles[doses %in% c(300, 600)] == "gate"))
  expect_equal(sum(roles == "external"), 4L)
})

test_that("[판정] (a) 연구 용량 300·600 mg: 관측 중앙값 tlast가 8명 코호트 2,000회 중앙값 분포의 5–95 백분위 안", {
  skip_if_judgment("판정성 테스트: results/step1/ 결과 파일 기반 보고 항목")
  for (g in Filter(function(x) x$gate_role == "gate", dz$groups)) {
    sim <- simulate_dataset(p, g$dose_mg, g$weight_mean, g$weight_sd, days, n_coh * n_per, paste0("clot_cohort_", g$dose_mg), sr, bm, bf)
    cm <- cohort_median_tlast(sim$nca, n_per)
    expect_true(g$tlast_median_obs >= q05(cm$tlast_median) & g$tlast_median_obs <= q95(cm$tlast_median),
                info = sprintf("%d mg: 관측 %g, 모의 [%g, %g]", g$dose_mg, g$tlast_median_obs, q05(cm$tlast_median), q95(cm$tlast_median)))
  }
})

test_that("[판정] (b) 연구 제형 데이터셋: AUClast 평균 ±15%, log-CV 35–51%", {
  skip_if_judgment("판정성 테스트: results/step1/ 결과 파일 기반 보고 항목")
  tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(unlist(dz$gate$log_cv_range_pct))
  tab <- rbindlist(lapply(Filter(function(d) d$gate_role == "gate", dz$datasets), function(ds) {
    b <- ds_weight_bounds(ds)
    sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), sched_of(ds), dz$gate$n_sim_per_dataset, paste0("gate_", ds$id), sr, b$male, b$female)
    gate_row(ds, sim$nca, tol, cvr)
  }))
  expect_equal(nrow(tab), 5L)
  expect_true(all(tab$pass_mean), info = paste(tab[pass_mean == FALSE, sprintf("%s %.3f", id, AUClast_ratio)], collapse = "; "))
  expect_true(all(tab$pass_cv), info = paste(tab[pass_cv == FALSE, sprintf("%s %.1f", id, AUClast_sim_logcv)], collapse = "; "))
})

test_that("[판정] (c) 체중 기울기: 5 kg당 AUClast 약 50 mg·day/L 감소 (허용 ±50%는 확정 전 가정)", {
  skip_if_judgment("판정성 테스트: results/step1/ 결과 파일 기반 보고 항목")
  wr <- as.numeric(unlist(dz$weight_slope_check$weight_range_kg))
  sim <- simulate_dataset(p, dz$weight_slope_check$dose_mg, mean(wr), 1e3, days, 20000, "wt_slope", sr, wr, wr)
  slope5 <- unname(coef(lm(AUClast ~ WT, data = sim$nca))[2]) * 5
  expect_lt(slope5, 0)
  expect_true(abs(slope5) >= 25 & abs(slope5) <= 75)
})

test_that("[판정] (f) 300 mg 개별 arm 외부 점검: 모델 평균(약 78 kg)이 Li 2020 arm 평균 범위 안", {
  skip_if_judgment("판정성 테스트: results/step1/ 결과 파일 기반 보고 항목")
  ds <- Filter(function(d) d$id == dz$arm_checks_300mg$model_reference_dataset, dz$datasets)[[1]]
  b <- ds_weight_bounds(ds)
  sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), sched_of(ds), dz$gate$n_sim_per_dataset, paste0("gate_", ds$id), sr, b$male, b$female)
  m <- mean(sim$nca$AUClast); rng <- range(vapply(dz$arm_checks_300mg$arms, function(a) a$auclast_mean, 0))
  expect_true(m >= rng[1] & m <= rng[2], info = sprintf("모델 %.0f, 범위 [%.1f, %.0f]", m, rng[1], rng[2]))
})
