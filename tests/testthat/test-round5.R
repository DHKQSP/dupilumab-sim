# 검토 의견 통합본(2026-09-24) 신규 기능 테스트
test_that("2020 Model 1: 자체 IIV(분산 = SD²)와 잔차, MTT IIV 포함 (D-029)", {
  p <- load_params("k2020")
  expect_equal(unname(p$omega2[c("Vc", "ke", "ka", "Vmax", "MTT")]), c(0.0369, 0.0812, 0.2247, 0.0557, 0.2756))
  expect_equal(unname(p$omega[c("Vc", "ke", "ka", "Vmax", "MTT")]), sqrt(c(0.0369, 0.0812, 0.2247, 0.0557, 0.2756)))
  expect_equal(unname(p$sigma), c(0.150, 0.03)); expect_equal(p$cov$WT_ref, 75); expect_equal(p$cov$theta_WT, 0.711)
  p16 <- load_params("k2016"); expect_equal(names(p16$omega), c("Vc", "ke", "k12", "k21", "ka", "Vmax", "Km", "F"))   # 2016 난수 순서 불변
})

test_that("MTT 개체간 변동이 ktr에 반영된다", {
  p <- load_params("k2020")
  ip <- individual_params(p, data.table(id = 1:2, WT = 75, eta_MTT = c(0, log(2))))
  expect_equal(ip$ktr[1], 4 / 0.105); expect_equal(ip$ktr[2], 4 / (0.105 * 2))
})

test_that("체중 밴드 균등분포와 키·BMI", {
  p <- load_params("k2016")
  s <- with_seed(3, make_subjects(20000, p, list(dist = "uniform", trunc = c(110, 130), height = list(mean = 170, sd = 9, trunc = c(150, 195))), 0.5, 0))
  expect_true(all(s$WT >= 110 & s$WT <= 130)); expect_equal(mean(s$WT), 120, tolerance = 0.01)
  expect_true(all(s$HT >= 150 & s$HT <= 195)); expect_equal(s$BMI, s$WT / (s$HT / 100)^2)
  s2 <- with_seed(3, make_subjects(10, p, weight_spec_from_design(read_cfg("trial_design.yaml"), "base"), 0.5, 0))
  expect_false("BMI" %in% names(s2))   # 키 분포가 없으면 추출하지 않음(기존 난수 순서 불변)
})

test_that("ke~BMI 공변량과 Vc 체중 지수 대체 (§6 c, d)", {
  sc <- load_scenarios(); design <- read_cfg("trial_design.yaml")
  pc <- resolve_variant("k2016_bmi", design, sc)$p; pd <- resolve_variant("k2016_bmi_vc0817", design, sc)$p
  e <- data.table(id = 1:3, WT = c(60, 75, 140), BMI = c(20, 26, 45))
  ic <- individual_params(pc, e); i0 <- individual_params(load_params("k2016"), e)
  expect_equal(ic$ke / i0$ke, (e$BMI / 26)^0.368)
  expect_equal(ic$Vc_i, i0$Vc_i)
  id <- individual_params(pd, e); expect_equal(id$Vc_i, 2.74 * (e$WT / 75)^0.817)
  expect_error(individual_params(pc, data.table(id = 1, WT = 75)), "BMI")
})

test_that("Wilson 구간, MC 평균 구간, 쌍대 비율 차이", {
  w <- wilson_ci(50, 100); expect_equal(w$est, 50); expect_equal(c(w$lo, w$hi), c(40.38, 59.62), tolerance = 1e-3)
  w0 <- wilson_ci(0, 500); expect_equal(w0$lo, 0); expect_gt(w0$hi, 0)
  m <- mc_mean_ci(c(1, 2, 3, 4)); expect_equal(m$est, 2.5); expect_equal(m$se, sd(1:4) / 2)
  d <- paired_prop_diff_ci(c(1, 1, 0, 1), c(1, 0, 0, 1)); expect_equal(d$est, 25)
})

test_that("trial_ids: 같은 시험 번호는 묶음과 무관하게 같은 결과(적응적 상향의 전제)", {
  skip_if_no_rxode2()
  p <- load_params("k2016"); design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); wt <- weight_spec_from_design(design, "base")
  combos <- CJ(scenario = c("S00", "F097"), schedule = "B0")
  a <- run_trials(3, p, design, sc$scenarios[c("S00", "F097")], combos, 5L, wt, progress_every = 0)
  b <- run_trials(1, p, design, sc$scenarios[c("S00", "F097")], combos, 5L, wt, progress_every = 0, trial_ids = 3L)
  expect_equal(a$be[trial == 3, GMR], b$be$GMR)
  expect_true("AUCinf_subC" %in% a$be$endpoint)
  expect_true(all(c("n_reliable", "n_dropout", "wt_dropout_mean") %in% names(a$ind)))
})
