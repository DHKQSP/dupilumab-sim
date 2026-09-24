# 운용특성 분석(R/oc.R, 검토 의견 2026-09-24 §3, D-040·D-042)의 구성 요소 검사

test_that("bisect_log는 단조 함수의 목표 배율을 허용오차 안에서 찾는다", {
  f <- function(m) m^-0.6                                # 배율이 커질수록 비가 작아지는 단조 함수
  r <- bisect_log(f, target = 0.80, bracket = c(1, 10), r_bracket = c(f(1), f(10)), tol_rel = 1e-4)
  expect_lt(abs(r$ratio / 0.80 - 1), 1e-4)
  expect_equal(r$m, 0.80^(-1 / 0.6), tolerance = 2e-4)
  expect_error(bisect_log(f, 1.2, c(1, 10), c(f(1), f(10)), 1e-4))   # 목표가 구간 밖이면 중단
})

test_that(".parabola_max는 간격이 다른 세 점에서도 포물선 꼭짓점을 돌려준다", {
  q <- function(x) 5 - 2 * (x - 1.3)^2
  x <- c(1, 1.25, 1.75)
  expect_equal(.parabola_max(x, q(x)), 5, tolerance = 1e-12)
  expect_equal(.parabola_max(c(0, 1, 2), c(1, 2, 3)), 3)          # 볼록하지 않으면 격자 최댓값
})

test_that("V2 기전은 k21만 나누고(Q = k12·Vc 고정) F 배율은 F <= 1을 강제한다", {
  p <- load_params("k2016")
  p2 <- apply_multipliers(p, list(V2 = 1.25))
  expect_equal(p2$theta[["k21"]], p$theta[["k21"]] / 1.25)
  expect_equal(p2$theta[["k12"]], p$theta[["k12"]])
  expect_equal(p2$theta[["Vc"]], p$theta[["Vc"]])
  expect_error(apply_multipliers(p, list(F = 1.01 / p$theta[["F"]])), "F")
  oc <- read_cfg("oc_design.yaml")
  expect_lte(mech_range(oc, "F", p)[2] * p$theta[["F"]], 1)
})

test_that("truth_metrics의 AUC0-inf는 선형 소실만 있을 때 F·Dose/(Vc·ke)와 같다(느슨한 참값 허용오차 포함)", {
  skip_if_no_rxode2()
  p <- load_params("k2016"); p$theta["Vmax"] <- 0
  ip <- individual_params(p, data.table(id = 1:4, WT = c(60, 70, 80, 90)))
  tm <- truth_metrics(ip, 300, p$model_id, cmax = TRUE)
  expect_equal(tm$AUCinf, ip$F * 300 / (ip$Vc_i * ip$ke), tolerance = 1e-5)
  expect_true(all(tm$Cmax > 0))
})

test_that("truth_ratio는 같은 제품에서 1, 알려진 배율에서 그 배율을 준다", {
  ref <- data.table(id = 1:5, AUCinf = c(100, 200, 300, 400, 500), Cmax = 1:5)
  r <- truth_ratio(ref, copy(ref))
  expect_equal(r$auc_ratio, 1); expect_equal(r$auc_se_log, 0); expect_equal(r$cmax_ratio, 1)
  r2 <- truth_ratio(ref, ref[, .(id, AUCinf = AUCinf * 0.8, Cmax = Cmax * 1.1)])
  expect_equal(r2$auc_ratio, 0.8); expect_equal(r2$cmax_ratio, 1.1)
})

test_that("config_pass는 사전 고정 구성 정의대로 계산하고 NA는 불통과로 본다", {
  oc <- read_cfg("oc_design.yaml")
  be <- CJ(trial = 1L, scenario = c("a", "b"), endpoint = OC_ENDPOINTS)
  be[, pass := TRUE]
  be[scenario == "b" & endpoint == "AUCinf_A", pass := NA]      # 규칙 A 산출 불가
  w <- config_pass(be, oc)
  expect_true(all(w[scenario == "a", c(cfg_P2, cfg_F3A, cfg_F3B, cfg_F3C, cfg_G2)]))
  expect_true(w[scenario == "b", cfg_P2]); expect_false(w[scenario == "b", cfg_F3A]); expect_false(w[scenario == "b", cfg_G2])
  expect_true(w[scenario == "b", cfg_F3C])
})
