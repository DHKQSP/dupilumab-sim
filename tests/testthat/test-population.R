design <- read_cfg("trial_design.yaml")
p <- load_params("k2016")

test_that("체중: 절단 정규분포와 성별 하한", {
  w <- with_seed(1, draw_weights(4000, list(mean = 75, sd = 9, trunc = c(60, 90)), 0.5))
  expect_true(all(w$WT >= 60 & w$WT <= 90)); expect_equal(mean(w$sex == "M"), 0.5, tolerance = 0.05)
  w2 <- with_seed(2, draw_weights(4000, list(mean = 72, sd = 10, trunc = c(50, 90), female_min = 50), 0.5))
  expect_true(all(w2$WT >= 50 & w2$WT <= 90))
  w3 <- with_seed(3, draw_weights(4000, list(mean = 60, sd = 9, male_bounds = c(50, 100), female_bounds = c(40, 90)), 0.5))
  expect_true(all(w3[sex == "M", WT] >= 50) && all(w3[sex == "F", WT] >= 40 & w3[sex == "F", WT] <= 90))
})

test_that("eta 표본 SD = omega, IIV 없는 파라미터는 0", {
  e <- with_seed(2, draw_etas(5000, p$omega))
  expect_equal(sd(e$eta_ke), sqrt(0.131), tolerance = 0.03); expect_equal(sd(e$eta_ka), sqrt(0.251), tolerance = 0.03)
  expect_true(all(e$eta_F == 0) && all(e$eta_k12 == 0))
})

test_that("개체 파라미터: 체중 공변량 지수 0.705, F는 고정 0.607, ADA 열", {
  e <- data.table(id = 1:3, WT = c(60, 75, 90), eta_Vc = 0, ada = c(0L, 1L, 0L))
  ip <- individual_params(p, e)
  expect_equal(ip$Vc_i, 2.74 * (c(60, 75, 90) / 75)^0.705)
  expect_equal(unique(ip$F), 0.607); expect_equal(ip$ada, c(0, 1, 0))
})

test_that("층화 무작위배정: 각 arm 정확히 n, 층 내 1:1(±1)", {
  s <- with_seed(4, make_subjects(234, p, weight_spec_from_design(design, "base"), 0.5, 0))
  a <- with_seed(5, assign_arms_stratified(s, design$stratification$weight_breaks, design$stratification$labels))
  expect_equal(sum(a$arm == "R"), 117); expect_equal(sum(a$arm == "T"), 117)
  tab <- a[, .N, by = .(stratum, arm)]
  for (st in unique(tab$stratum)) expect_lte(abs(diff(tab[stratum == st][order(arm), N])), 1)
  expect_true(all(a[WT <= 75, stratum] == "60-75"))
})

test_that("층화 배정: 층 경계 밖 체중(50–115)도 전원 배정된다 (D-020 회귀 테스트)", {
  s <- with_seed(6, make_subjects(234, p, list(mean = 82.5, sd = 1e3, trunc = c(50, 115)), 0.5, 0))
  expect_true(any(s$WT < 60) && any(s$WT > 90))
  a <- with_seed(7, assign_arms_stratified(s, design$stratification$weight_breaks, design$stratification$labels))
  expect_false(anyNA(a$arm)); expect_false(anyNA(a$stratum))
  expect_equal(sum(a$arm == "R"), 117); expect_equal(sum(a$arm == "T"), 117)
  expect_true(all(a[WT < 60, stratum] == "60-75")); expect_true(all(a[WT > 90, stratum] == ">75-90"))
})
