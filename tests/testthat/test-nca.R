# NCA 규칙 검증 (지시서 §3, config/nca_rules.yaml). 손계산 예제와 lm() 대조.
mk <- function(time, conc, id = 1L) data.table(id = id, time = time, conc = conc)

test_that("AUClast: linear-up/log-down 손계산과 일치, (0,0) 앵커 포함", {
  d <- mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 10, 20, 16, 8, 4, 1))
  r <- run_nca(d)
  seg <- c(5, 15, (20 - 16) / log(20 / 16) * 2, (16 - 8) / log(2) * 4, (8 - 4) / log(2) * 4, (4 - 1) / log(4) * 12)
  expect_equal(r$AUClast, sum(seg), tolerance = 1e-10)
  expect_equal(r$Cmax, 20); expect_equal(r$tmax, 2); expect_equal(r$tlast, 24); expect_equal(r$Clast, 1); expect_equal(r$n_quant, 6L)
})

test_that("λz: 마지막 k점 후보 중 adjusted R² 최대 창 선택, AUCinf·외삽률·신뢰 기준", {
  d <- mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 10, 20, 16, 8, 4, 1))
  r <- run_nca(d)
  fit3 <- lm(log(c(8, 4, 1)) ~ c(8, 12, 24)); fit4 <- lm(log(c(16, 8, 4, 1)) ~ c(4, 8, 12, 24))
  a3 <- summary(fit3)$adj.r.squared; a4 <- summary(fit4)$adj.r.squared
  expect_gt(a3, a4)
  expect_equal(r$lambda_z, -unname(coef(fit3)[2]), tolerance = 1e-10)
  expect_equal(r$adj_r2, a3, tolerance = 1e-10)
  expect_equal(r$n_lambda, 3L)
  expect_equal(r$AUCinf, r$AUClast + 1 / r$lambda_z)
  expect_equal(r$pct_extrap, (r$AUCinf - r$AUClast) / r$AUCinf * 100)
  expect_true(r$lambda_ok && r$reliable)
  expect_equal(r$t_half, log(2) / r$lambda_z)
})

test_that("λz 동률(정확한 단일지수)이면 점 수가 많은 창을 택한다; Cmax 점은 후보 제외", {
  t <- c(0, 1, 2, 4, 8, 12, 24, 36); cc <- c(NA, 10, 20, 20 * exp(-0.1 * (c(4, 8, 12, 24, 36) - 2)))
  r <- run_nca(mk(t, cc))
  expect_equal(r$lambda_z, 0.1, tolerance = 1e-8)
  expect_equal(r$n_lambda, 5L)      # tmax 이후 5점 전부(Cmax 제외)
})

test_that("BLQ 규칙: tmax 이전 0, 중간 결측(인접점 연결), tlast 이후 제외", {
  d <- mk(c(0, 0.5, 1, 2, 4, 8, 12, 24, 36), c(NA, NA, 10, 20, 16, NA, 4, 1, NA))
  r <- run_nca(d)
  seg <- c(0, (0 + 10) / 2 * 0.5, 15, (20 - 16) / log(20 / 16) * 2, (16 - 4) / log(4) * 8, (4 - 1) / log(4) * 12)
  expect_equal(r$AUClast, sum(seg), tolerance = 1e-10)
  expect_equal(r$tlast, 24)
  expect_equal(r$n_lambda, 3L)   # 후보: (4,16),(12,4),(24,1)
})

test_that("λz 산출 불가: tmax 이후 정량점 3개 미만 또는 기울기 ≥ 0", {
  r1 <- run_nca(mk(c(0, 1, 2, 4, 8), c(NA, 10, 20, 16, 8)))
  expect_false(r1$lambda_ok); expect_true(is.na(r1$AUCinf)); expect_false(r1$reliable)
  r2 <- run_nca(mk(c(0, 1, 2, 4, 8, 12), c(NA, 10, 20, 5, 6, 7)))
  expect_false(r2$lambda_ok)
})

test_that("정량치가 전혀 없는 개체는 NA 행으로 남고 여러 개체를 한 번에 처리한다", {
  d <- rbind(mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 10, 20, 16, 8, 4, 1), id = 1L),
             mk(c(0, 1, 2, 4), c(NA, NA, NA, NA), id = 2L),
             mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 5, 10, 8, 4, 2, 0.5), id = 3L))
  r <- run_nca(d)
  expect_equal(r$id, 1:3)
  expect_equal(r[id == 2, n_quant], 0L); expect_true(is.na(r[id == 2, AUClast]))
  expect_equal(r[id == 3, AUClast], r[id == 1, AUClast] / 2, tolerance = 1e-10)
  expect_equal(r[id == 3, lambda_z], r[id == 1, lambda_z], tolerance = 1e-10)
})

test_that("BEmaster 호환 모드: 선형 사다리꼴, BLQ=0, 마지막 3점 λz", {
  d <- mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 10, 20, 16, 8, 4, 1))
  r <- run_nca(d, mode = "bemaster_compat")
  expect_equal(r$AUClast, sum(diff(c(0, 1, 2, 4, 8, 12, 24)) * (head(c(0, 10, 20, 16, 8, 4, 1), -1) + tail(c(0, 10, 20, 16, 8, 4, 1), -1)) / 2))
  fit3 <- lm(log(c(8, 4, 1)) ~ c(8, 12, 24))
  expect_equal(r$lambda_z, -unname(coef(fit3)[2]), tolerance = 1e-10)
  expect_true(is.na(r$reliable))
})

test_that("참값 지표 부착", {
  nca <- data.table(id = 1L, AUClast = 90, tlast = 24)
  subj <- data.table(id = 1L, AUClast_true = 95, AUCinf_true = 100)
  m <- add_truth_metrics(nca, subj)
  expect_equal(m$coverage_true, 0.95); expect_equal(m$pct_extrap_true, 5)
  expect_equal(m$err_AUClast_vs_true_tlast, log(90 / 95)); expect_equal(m$err_AUClast_vs_true_inf, log(0.9))
})
