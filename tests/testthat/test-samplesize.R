# 표본 수 분석식(R/samplesize.R, 추가 지시 2026-09-26 §3): 단조성, M1 ≥ M0, 경계 사례, 직접 모의와의 일치

inp <- data.table(sd_cmax = 0.3336, ratio_sd_cmax_auc = 0.86, rho_total = 0.753, beta_auc = -1.18, beta_cmax = -0.86)
wm <- ss_weight_moments(75, 9, 60, 90, 75)

test_that("log 체중 분산 분해는 절단 정규 표본과 맞고 층 간 + 층 안 = 전체", {
  expect_equal(wm$var_total, wm$var_between + wm$var_within, tolerance = 1e-12)
  set.seed(3); w <- log(truncnorm::rtruncnorm(4e5, 60, 90, 75, 9))
  expect_equal(wm$var_total, var(w), tolerance = 0.01)
  s2 <- w > log(75); vb <- mean(s2) * (mean(w[s2]) - mean(w))^2 + mean(!s2) * (mean(w[!s2]) - mean(w))^2
  expect_equal(wm$var_between, vb, tolerance = 0.02)
})

test_that("P2 검정력은 n에 대해 증가하고, CV에 대해 감소하며, M1이 M0 이상이다", {
  p <- function(n, cv, am, g = 0.95) ss_p2_power(n, cv, g, g, inp, am, "fixed", wm)
  expect_true(all(diff(sapply(c(60, 100, 140, 200), p, cv = 43, am = "M1")) > 0))
  expect_true(all(diff(sapply(c(35, 43, 52), function(cv) p(117, cv, "M1"))) < 0))
  for (cv in c(35, 43, 52)) expect_gte(p(117, cv, "M1"), p(117, cv, "M0"))
  expect_equal(ss_p2_power(117, 43, 0.70, 0.70, inp, "M0", "fixed", wm), 0, tolerance = 1e-6)   # 참값이 범위 밖이면 거의 0
  expect_gt(p(400, 35, "M0", 1), 0.999)
})

test_that("분석식은 같은 자료 생성 모형의 직접 모의(M1, 4,000회)와 모의 오차 안에서 같다", {
  set.seed(11); n <- 117; cv <- 43; g <- 0.95; B <- 4000
  cc <- ss_cell_cov(inp, cv, "fixed", wm)
  w <- truncnorm::rtruncnorm(B * 2 * n, 60, 90, 75, 9)
  d <- data.table(tid = rep(seq_len(B), each = 2 * n), stratum = as.integer(w > 75), lwt = log(w), k = runif(B * 2 * n))
  d[, r := frank(k, ties.method = "first"), by = .(tid, stratum)][, ns := .N, by = .(tid, stratum)]
  d[, arm := fifelse(r <= ns %/% 2L, "R", fifelse(r <= 2L * (ns %/% 2L), "T", NA_character_))]
  d[is.na(arm), arm := { o <- order(k); a <- character(.N); a[o] <- rep(c("R", "T"), length.out = .N); a }, by = tid]
  z <- matrix(rnorm(2 * nrow(d)), ncol = 2) %*% chol(matrix(c(cc$veA, cc$ce, cc$ce, cc$veC), 2))
  tt <- as.numeric(d$arm == "T"); lc <- d$lwt - wm$mean_lwt
  dA <- d[, .(tid, arm, stratum, lwt, y = exp(log(g) * tt + inp$beta_auc * lc + z[, 1]))]; dC <- d[, .(tid, arm, stratum, lwt, y = exp(log(g) * tt + inp$beta_cmax * lc + z[, 2]))]
  pa <- be_models_fast(dA, "tid")[model == "M1"]; pc <- be_models_fast(dC, "tid")[model == "M1"]
  mc <- mean(pa$pass & pc$pass)
  an <- ss_p2_power(n, cv, g, g, inp, "M1", "fixed", wm)
  expect_lt(abs(mc - an), 4 * sqrt(an * (1 - an) / B))
})
