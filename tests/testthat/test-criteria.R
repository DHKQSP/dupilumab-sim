# λz 신뢰 기준 세트(R/criteria.R, 두 번째 추가 지시 2026-09-26 §1)와 로그정규 체중(R/population.R, §2)

test_that("기준 세트 (i)·(ii)는 NCA 출력의 기존 플래그와 같고, (iii) ⊂ (i), (iv) ⊂ (iii)", {
  f <- proj_path("results", "individual", "nca_base_20000.rds"); skip_if_not(file.exists(f), "저장 NCA 없음")
  x <- readRDS(f)[schedule == "B0"]
  rel_i <- (x$lambda_ok & !x$flag_rsq & !x$flag_extrap) %in% TRUE
  expect_identical(crit_ok(x, "i"), rel_i)
  expect_identical(crit_ok(x, "ii"), x$reliable %in% TRUE)
  expect_true(all(crit_ok(x, "iii") <= crit_ok(x, "i"))); expect_true(all(crit_ok(x, "iv") <= crit_ok(x, "iii")))
  expect_lt(mean(crit_ok(x, "iii")), mean(crit_ok(x, "i")))
})

test_that("로그정규 체중: 절단 구간 안, 산술 평균·SD 매개변수화, 기존 분포의 난수열 불변", {
  spec <- list(dist = "lognormal", mean = 78, sd = 19, trunc = c(40, 180))
  w <- with_seed(5L, draw_weights(2e5, spec, 0.5))$WT
  expect_true(all(w >= 40 & w <= 180))
  expect_equal(mean(w), 78, tolerance = 0.01); expect_equal(sd(w), 19, tolerance = 0.03)
  sl <- sqrt(log(1 + (19 / 78)^2)); ml <- log(78) - sl^2 / 2
  expect_equal(plnorm(60, ml, sl), 0.1653, tolerance = 2e-3); expect_equal(1 - plnorm(90, ml, sl), 0.2370, tolerance = 2e-3); expect_equal(1 - plnorm(100, ml, sl), 0.1240, tolerance = 2e-3)
  base <- list(mean = 75, sd = 9, trunc = c(60, 90))
  a <- with_seed(9L, draw_weights(50, base, 0.5)); b <- with_seed(9L, draw_weights(50, base, 0.5)); expect_identical(a, b)
})
