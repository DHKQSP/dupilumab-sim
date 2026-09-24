test_that("pooled t 90% CI는 t.test(var.equal=TRUE)와 일치하고 80–125% 판정이 맞다", {
  set.seed(3)
  y <- exp(c(rnorm(50, log(500), 0.4), rnorm(50, log(520), 0.4))); arm <- rep(c("R", "T"), each = 50)
  r <- be_pooled_t(y, arm)
  tt <- t.test(log(y[arm == "T"]), log(y[arm == "R"]), var.equal = TRUE, conf.level = 0.90)
  expect_equal(c(r$CI_lower, r$CI_upper), exp(tt$conf.int[1:2]), tolerance = 1e-12)
  expect_equal(r$GMR, exp(mean(log(y[arm == "T"])) - mean(log(y[arm == "R"]))))
  expect_equal(r$pass, r$CI_lower >= 0.8 && r$CI_upper <= 1.25)
  expect_equal(r$df, 98)
})

test_that("ANCOVA(체중) CI는 lm/confint와 일치", {
  set.seed(4)
  WT <- runif(100, 60, 90); arm <- rep(c("R", "T"), each = 50)
  y <- exp(log(500) - 0.01 * (WT - 75) + rnorm(100, 0, 0.3) + ifelse(arm == "T", 0.05, 0))
  r <- be_ancova_weight(y, arm, WT)
  fit <- lm(log(y) ~ factor(arm, levels = c("R", "T")) + WT)
  ci <- confint(fit, level = 0.90)[2, ]
  expect_equal(c(r$CI_lower, r$CI_upper), unname(exp(ci)), tolerance = 1e-12)
})

test_that("be_analyze는 평가변수별·부분집합별로 계산한다", {
  set.seed(5)
  d <- data.table(arm = rep(c("R", "T"), each = 40), WT = runif(80, 60, 90),
                  Cmax = exp(rnorm(80, 3, 0.3)), AUClast = exp(rnorm(80, 6, 0.4)), AUCinf = exp(rnorm(80, 6.1, 0.4)),
                  AUCinf_true = exp(rnorm(80, 6.05, 0.4)), lambda_ok = rep(c(TRUE, FALSE), 40), reliable = rep(c(TRUE, FALSE, FALSE, TRUE), 20))
  d[, AUCinf_subC := fifelse(reliable %in% TRUE, AUCinf, AUClast)]   # 규칙 C(통합본 §5-5): trial.R과 같은 정의
  r <- be_analyze(d, BE_ENDPOINTS, methods = c("pooled_t", "ancova_weight"))
  expect_equal(nrow(r), length(BE_ENDPOINTS) * 2)
  expect_equal(r[endpoint == "AUCinf_all" & method == "pooled_t", n_R + n_T], sum(d$lambda_ok))
  expect_equal(r[endpoint == "AUCinf_reliable" & method == "pooled_t", n_R + n_T], sum(d$reliable))
  expect_equal(r[endpoint == "AUCinf_subC" & method == "pooled_t", n_R + n_T], nrow(d))   # 규칙 C는 전원 포함
  expect_error(be_analyze(d[, !"AUCinf_subC"], BE_ENDPOINTS), "평가변수 열 없음")
})
