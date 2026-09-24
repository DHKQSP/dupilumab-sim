# Phoenix WinNonlin 호환 NCA 엔진 규칙 (검토 의견 2026-09-24 §1, D-039)
mk <- function(time, conc, id = 1L) data.table(id = id, time = time, conc = conc)

test_that("Phoenix 표기 출력 변수와 Clast_pred·AUCINF_pred·span ratio", {
  d <- mk(c(0, 1, 2, 4, 8, 12, 24), c(NA, 10, 20, 16, 8, 4, 1))
  r <- run_nca(d)
  expect_true(all(PHOENIX_COLS %in% names(r)))
  fit <- lm(log(c(8, 4, 1)) ~ c(8, 12, 24))
  expect_equal(r$Clast_pred, unname(exp(coef(fit)[1] + coef(fit)[2] * 24)), tolerance = 1e-10)
  expect_equal(r$AUCINF_pred, r$AUClast + r$Clast_pred / r$Lambda_z)
  expect_equal(r$`AUC_%Extrap_obs`, (1 - r$AUClast / r$AUCINF_obs) * 100)
  expect_equal(r$Rsq, summary(fit)$r.squared, tolerance = 1e-10)
  expect_equal(c(r$Lambda_z_lower, r$Lambda_z_upper), c(8, 24))
  expect_equal(r$Span_ratio, (24 - 8) / (log(2) / r$Lambda_z))
  expect_equal(r$AUCinf, r$AUCINF_obs)                       # 하위 호환 별칭
})

test_that("BLQ: 첫 정량값 이전은 0, 첫 정량 이후 tmax 이전의 BLQ는 결측(이전 엔진은 0)", {
  d <- mk(c(0, 0.5, 1, 1.5, 2, 4, 8, 12, 24), c(NA, 5, NA, 15, 20, 16, 8, 4, 1))
  r <- run_nca(d); r0 <- run_nca_legacy(d)
  # 새 엔진: (1, BLQ) 행 제거 → 0.5→1.5 선형 연결
  seg <- c((0 + 5) / 2 * 0.5, (5 + 15) / 2 * 1, (15 + 20) / 2 * 0.5, (20 - 16) / log(20 / 16) * 2, (16 - 8) / log(2) * 4, (8 - 4) / log(2) * 4, (4 - 1) / log(4) * 12)
  expect_equal(r$AUClast, sum(seg), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(r$AUClast, r0$AUClast)))    # 이전 엔진은 tmax 이전 BLQ를 0으로 넣음
})

test_that("BLQ 2회 연속 뒤의 정량값은 결측(그 이후 전부 제외), 단일 BLQ는 결측으로 연결", {
  d <- mk(c(0, 1, 2, 4, 8, 12, 16, 24, 36), c(NA, 10, 20, 16, NA, NA, 3, 1, NA))
  r <- run_nca(d)
  expect_equal(r$Tlast, 4); expect_equal(r$Clast, 16); expect_false(r$lambda_ok)     # Cmax 이후 양수 1점뿐
  d2 <- mk(c(0, 1, 2, 4, 8, 12, 16, 24), c(NA, 10, 20, 16, 8, NA, 3, 1))
  r2 <- run_nca(d2)
  expect_equal(r2$Tlast, 24)
  expect_equal(r2$AUClast, sum(c(5, 15, (20 - 16) / log(20 / 16) * 2, (16 - 8) / log(2) * 4, (8 - 3) / log(8 / 3) * 8, (3 - 1) / log(3) * 8)), tolerance = 1e-10)
  pp <- preprocess_blq(d)
  expect_equal(pp$time, c(0, 1, 2, 4)); expect_equal(pp$conc[1], 0)
})

test_that("동률 판정은 엄격 부등호(|최대 − adj R²| < 1e-4), 기울기 > 0 창은 선택 전에 제외", {
  # 말단 상승: 마지막 3점은 기울기 양수 → 제외, 나머지 창 중 최대
  d <- mk(c(0, 1, 2, 4, 8, 12, 16, 24), c(NA, 10, 20, 12, 7, 6, 6.6, 7.3))
  r <- run_nca(d)
  expect_true(r$lambda_ok); expect_equal(r$No_points_lambda_z, 5L)
  skip_if_not_installed("NonCompart")
  nc <- NonCompart::sNCA(c(0, 1, 2, 4, 8, 12, 16, 24), c(0, 10, 20, 12, 7, 6, 6.6, 7.3), dose = 1, down = "Log", R2ADJ = 0)
  expect_equal(r$Lambda_z, unname(nc["LAMZ"]), tolerance = 1e-12)
})

test_that("신뢰 플래그 세 가지와 reliable, 규칙 C 대체값", {
  # 말단 3점이 짧은 간격 → span ratio < 2
  d <- mk(c(0, 1, 2, 4, 8, 10, 11, 12), c(NA, 10, 20, 16, 8, 6.9, 6.4, 6.0))
  r <- run_nca(d)
  expect_true(r$lambda_ok); expect_true(r$flag_span); expect_false(r$reliable)
  expect_lt(r$Span_ratio, 2)
  dr <- dropout_reasons(r)
  expect_equal(dr$flag_span_pct, 100)
})

test_that("무작위 모의 프로필 50개에서 NonCompart와 λz 점·파라미터 일치", {
  skip_if_not_installed("NonCompart")
  set.seed(11)
  tt <- c(0, 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56)
  prof <- rbindlist(lapply(1:50, function(i) {
    k <- exp(rnorm(1, log(0.08), 0.3)); ka <- exp(rnorm(1, log(0.3), 0.3))
    cc <- 30 * (exp(-k * tt) - exp(-ka * tt)) * exp(rnorm(length(tt), 0, 0.15)); cc[cc < 0.5] <- NA; cc[1] <- NA
    mk(tt, cc, id = i)
  }))
  own <- run_nca(prof)
  pp <- preprocess_blq(prof)
  for (i in unique(pp$id)) {
    x <- pp[id == i]; nc <- NonCompart::sNCA(x$time, x$conc, dose = 1, down = "Log", R2ADJ = 0)
    o <- own[id == i]
    if (is.na(nc["LAMZ"])) { expect_false(o$lambda_ok); next }
    expect_equal(o$No_points_lambda_z, as.integer(nc["LAMZNPT"]))
    expect_equal(c(o$Lambda_z, o$AUClast, o$AUCINF_obs, o$AUCINF_pred, o$Rsq_adjusted), unname(nc[c("LAMZ", "AUCLST", "AUCIFO", "AUCIFP", "R2ADJ")]), tolerance = 1e-8)
  }
})
