skip_if_no_rxode2()
p <- load_params("dev")
design <- read_cfg("trial_design.yaml")

test_that("잔차 오차와 BLQ 검열이 규칙대로 적용된다", {
  ip <- typical_subject(p)
  obs <- data.table(id = 1L, planned = c(0, 1, 7, 28, 56, 120), time = c(0, 1, 7, 28, 56, 120))
  sig <- c(prop = 0, add = 0)
  sim <- simulate_observations(ip, obs, 300, sig, lloq = 0.078)
  expect_equal(sim$y_raw, sim$C)                       # 오차 0이면 관측 = 진농도
  expect_true(sim[time == 0, blq])                      # 투여 전 BLQ
  expect_true(sim[time == 120, blq])                    # 120일: 소실 완료
  expect_true(all(is.na(sim[blq == TRUE, conc])))
  sig2 <- c(prop = 0.2, add = 0.05)
  s1 <- with_seed(5, simulate_observations(ip, obs, 300, sig2, 0.078))
  s2 <- with_seed(5, simulate_observations(ip, obs, 300, sig2, 0.078))
  expect_identical(s1$y_raw, s2$y_raw)                  # 시드 재현
  expect_false(isTRUE(all.equal(s1$y_raw[2:4], s1$C[2:4])))
})

test_that("run_trial은 BEmaster 없이 NCA 경계까지 동작한다(nca_fn=NULL)", {
  tr <- run_trial(p, p, design, seed = 11, n_per_arm = 4, nca_fn = NULL, be_fn = NULL, keep_profiles = TRUE)
  expect_equal(nrow(tr$subjects), 8)
  expect_true(all(tr$subjects$ratio_true > 0 & tr$subjects$ratio_true <= 1))
  expect_true(all(c("R", "T") %in% tr$subjects$arm))
  expect_null(tr$nca); expect_null(tr$be)
  expect_true(all(tr$sim[time == 0, blq]))
  # 같은 시드면 같은 결과
  tr2 <- run_trial(p, p, design, seed = 11, n_per_arm = 4, nca_fn = NULL, be_fn = NULL)
  expect_equal(tr$subjects$AUCinf_true, tr2$subjects$AUCinf_true)
})

test_that("시나리오 적용: T arm만 바뀐다", {
  sc <- load_scenarios()
  pp <- apply_scenario(p, sc$main$S1b)
  expect_equal(pp$R$theta[["F"]], p$theta[["F"]])
  expect_equal(pp$T$theta[["F"]], p$theta[["F"]] * 0.9)
  expect_error(apply_scenario(p, sc$sensitivity$K1), "PENDING")
})

test_that("바깥층 theta 추출은 시드에 결정적이고 RSE 0이면 대표값 그대로", {
  th1 <- draw_outer_theta(p, 7); th2 <- draw_outer_theta(p, 7)
  expect_identical(th1, th2)
  p0 <- p; p0$rse_pct[] <- 0
  expect_equal(draw_outer_theta(p0, 7), p$theta)
})
