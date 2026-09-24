test_that("기본 파라미터는 Kovalenko 2016 Table 2 값이고 전부 confirmed", {
  p <- load_params("k2016")
  expect_equal(unname(p$theta[c("Vc", "ke", "k12", "k21", "ka", "Vmax", "Km", "F")]), c(2.74, 0.0459, 0.0652, 0.129, 0.254, 0.968, 0.01, 0.607))
  expect_equal(p$cov$theta_WT, 0.705); expect_equal(p$cov$WT_ref, 75)
  expect_equal(unname(p$omega2[c("Vc", "ke", "ka", "Vmax")]), c(0.0225, 0.131, 0.251, 0.0428))
  expect_equal(unname(p$omega2[c("k12", "k21", "F", "Km")]), rep(0, 4))
  expect_equal(unname(p$omega), sqrt(unname(p$omega2)))
  expect_equal(unname(p$sigma), c(0.242, 0.03)); expect_equal(p$lloq, 0.078)
  expect_true(all(p$status$status == "confirmed"))
  expect_error(assert_params_confirmed(p), NA)
})

test_that("Kovalenko 2020 Model 1 파라미터와 IIV 차용", {
  p <- load_params("k2020")
  expect_equal(unname(p$theta[c("Vc", "ke", "k12", "k21", "ka", "MTT", "n_transit", "Vmax", "Km", "F")]), c(2.48, 0.0534, 0.213, 0.310, 0.256, 0.105, 3, 1.07, 0.01, 0.643))
  expect_equal(p$cov$theta_WT, 0.711)
  expect_equal(p$omega2[["ke"]], 0.0812)   # 자체 IIV(D-029)
})

test_that("시나리오 배율과 민감도 변형", {
  p <- load_params("k2016")
  q <- apply_multipliers(p, list(F = 0.9, Km = 10)); expect_equal(q$theta[["F"]], 0.607 * 0.9); expect_equal(q$theta[["Km"]], 0.1); expect_equal(q$theta[["ke"]], 0.0459)
  v <- apply_variant(p, list(omega2_multiplier = 1.5)); expect_equal(v$omega2[["ke"]], 0.131 * 1.5); expect_equal(v$omega[["ke"]], sqrt(0.131 * 1.5))
  r <- apply_variant(p, list(sigma_prop = 0.12)); expect_equal(r$sigma[["prop"]], 0.12)
  design <- read_cfg("trial_design.yaml")
  a <- apply_variant(p, list(ada = TRUE), design); expect_equal(a$ada$fraction, 0.10); expect_equal(a$ada$ke_multiplier, 2)
})

test_that("CV <-> omega 변환과 log-CV", {
  expect_equal(omega_to_cv(cv_to_omega(30)), 30)
  expect_equal(log_cv_pct(exp(rnorm(1e5, 0, 0.4))), 100 * sqrt(exp(0.16) - 1), tolerance = 0.02)
})
