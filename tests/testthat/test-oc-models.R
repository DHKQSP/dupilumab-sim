# 분석 모형 M0/M1/M2(R/oc_models.R, 추가 지시 2026-09-26 §1)의 닫힌 꼴 추정을 lm/confint·be_pooled_t와 대조

mk_trial <- function(seed, n = 117, rule_a_drop = 0) {
  set.seed(seed)
  wt <- pmin(pmax(rnorm(2 * n, 75, 9), 60), 90)
  arm <- rep(c("R", "T"), each = n)
  y <- exp(6 + 0.1 * (arm == "T") - 0.012 * (wt - 75) + rnorm(2 * n, 0, 0.4))
  if (rule_a_drop > 0) y[sample(2 * n, rule_a_drop)] <- NA            # 규칙 A처럼 일부 제외(군별 n이 달라짐)
  data.table(tid = seed, arm = arm, y = y, stratum = as.integer(wt > 75), lwt = log(wt), wt = wt)
}

test_that("M1(처리 + 체중 층)과 M2(처리 + log 체중)의 GMR·90% CI·자유도가 lm/confint와 일치한다", {
  d <- rbindlist(lapply(1:6, function(s) mk_trial(s, rule_a_drop = if (s %% 2) 20 else 0)))
  r <- be_models_fast(copy(d), "tid")
  for (s in 1:6) {
    z <- d[tid == s & is.finite(y)]
    f1 <- lm(log(y) ~ factor(arm, c("R", "T")) + factor(stratum), data = z)
    f2 <- lm(log(y) ~ factor(arm, c("R", "T")) + lwt, data = z)
    for (m in list(list("M1", f1), list("M2", f2))) {
      ci <- exp(confint(m[[2]], level = 0.90)[2, ]); rr <- r[tid == s & model == m[[1]]]
      expect_equal(rr$GMR, unname(exp(coef(m[[2]])[2])), tolerance = 1e-10)
      expect_equal(c(rr$CI_lower, rr$CI_upper), unname(ci), tolerance = 1e-10)
      expect_equal(rr$df, m[[2]]$df.residual)
      expect_equal(rr$n_R + rr$n_T, nrow(z))
    }
  }
})

test_that("층별 군 인원이 같으면 M1 점추정은 M0(군 평균 차)와 같고 SE는 층 효과만큼 작다", {
  set.seed(7); n <- 60
  d <- data.table(tid = 1L, arm = rep(c("R", "T"), each = 2 * n), stratum = rep(rep(0:1, each = n), 2))
  d[, y := exp(6 + 0.05 * (arm == "T") + 0.3 * stratum + rnorm(.N, 0, 0.4))][, lwt := log(70 + 10 * stratum + rnorm(.N))]
  r1 <- be_models_fast(copy(d), "tid")[model == "M1"]; r0 <- be_m0_fast(copy(d), "tid")
  expect_equal(r1$est, r0$est, tolerance = 1e-12); expect_lt(r1$se, r0$se)
})

test_that("be_m0_fast는 be_pooled_t와 같고, 층이 비면 M1은 NA(불통과)다", {
  d <- rbindlist(lapply(11:14, function(s) mk_trial(s)))
  r0 <- be_m0_fast(copy(d), "tid")
  for (s in 11:14) {
    z <- d[tid == s]; b <- be_pooled_t(z$y, z$arm)
    expect_equal(r0[tid == s, c(GMR, CI_lower, CI_upper)], c(b$GMR, b$CI_lower, b$CI_upper), tolerance = 1e-12)
  }
  one <- mk_trial(21)[, stratum := 0L]                                  # 모든 대상자가 한 층 → M1 특이
  r1 <- be_models_fast(one, "tid")
  expect_true(is.na(r1[model == "M1", GMR])); expect_true(r1[model == "M2", is.finite(GMR)])
})

test_that("run_trial_oc_models의 M0는 run_trial_oc_ext와 같고, M1은 층 균형 시 점추정이 M0와 같으며, LLOQ 재검열은 모의 LLOQ에서 원래 결과와 같다", {
  skip_if_no_rxode2()
  skip_if_fast_scope("시험 1회 모의(시나리오 2개, LLOQ 3개), 약 10초")
  oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
  rv <- resolve_variant("base", design); p <- rv$p
  scen <- list(VM = list(code = "VM", T_multipliers = list(Vmax = 0.66)), F097 = list(code = "F097", T_multipliers = list(F = 0.97)))
  ms <- as.integer(oc$trials$master_seed)
  a <- run_trial_oc_ext(3L, p, design, scen, ms, rv$wt_spec, model_id = p$model_id)
  b <- run_trial_oc_models(3L, p, design, scen, ms, rv$wt_spec, model_id = p$model_id)
  m0 <- b$be[model == "M0"][a$be, on = c("scenario", "endpoint")]
  expect_equal(nrow(m0), nrow(a$be))
  expect_identical(m0$GMR, m0$i.GMR); expect_identical(m0$CI_lower, m0$i.CI_lower); expect_identical(m0$CI_upper, m0$i.CI_upper)
  expect_identical(m0$pass, m0$i.pass); expect_equal(m0$n_R, m0$i.n_R); expect_equal(m0$n_T, m0$i.n_T)
  expect_true(all(b$be[model == "M0", df] == b$be[model == "M0", n_R + n_T - 2]))
  expect_true(all(b$be[model != "M0" & is.finite(se), df] == b$be[model != "M0" & is.finite(se), n_R + n_T - 3]))
  # 층 안 1:1 배정: 층 인원이 홀수면 남는 1명을 층 간 교대 배정하므로 군별 층 인원 차이는 최대 1(시험 3은 차이 1 → M1 점추정이 M0와 조금 다르다)
  st <- dcast(b$strata, stratum ~ arm, value.var = "n"); expect_true(all(abs(st$R - st$T) <= 1)); expect_equal(sum(st$R), sum(st$T))
  full <- b$be[endpoint %in% c("AUClast", "Cmax", "AUCinf_true")]
  e <- dcast(full, scenario + endpoint ~ model, value.var = "est")
  expect_true(all(abs(e$M1 - e$M0) < 0.01))
  expect_true(all(dcast(full, scenario + endpoint ~ model, value.var = "se")[, M1 < M0]))   # 층 효과가 있어 잔차 SD가 작다
  # LLOQ 재검열: 모의 LLOQ 포함 3개를 한 번에 → 모의 LLOQ 행은 단일 실행과 같고, 다른 LLOQ는 같은 대상자(n 동일)로 다른 값
  L0 <- study_lloq(); expect_equal(L0, p$lloq)                         # 현재 연구 LLOQ = 모델 개발 자료 LLOQ(0.078)
  c3 <- run_trial_oc_models(3L, p, design, scen, ms, rv$wt_spec, model_id = p$model_id, lloqs = c(0.02, L0, 0.5), resid = c("fixed", "scaled"), models = c("M0", "M1"))
  same <- c3$be[lloq == L0 & resid == "fixed"][b$be[model %in% c("M0", "M1")], on = c("scenario", "endpoint", "model")]
  expect_identical(same$est, same$i.est); expect_identical(same$se, same$i.se)
  # 가산 잔차 배율 = LLOQ / p$lloq: 연구 LLOQ(배율 1)에서 scaled = fixed(재추출 잔차가 저장 y_raw를 비트 단위로 재현)
  sc1 <- c3$be[lloq == L0 & resid == "scaled"][c3$be[lloq == L0 & resid == "fixed"], on = c("scenario", "endpoint", "model")]
  expect_identical(sc1$est, sc1$i.est); expect_identical(sc1$se, sc1$i.se)
  lo <- c3$be[lloq == 0.02 & resid == "fixed" & endpoint == "AUClast" & model == "M0"]; hi <- c3$be[lloq == L0 & resid == "fixed" & endpoint == "AUClast" & model == "M0"]
  expect_identical(lo$n_R, hi$n_R); expect_false(isTRUE(all.equal(lo$est, hi$est)))
  expect_false(isTRUE(all.equal(c3$be[lloq == 0.02 & resid == "scaled", est], c3$be[lloq == 0.02 & resid == "fixed", est])))
  expect_identical(c3$be[lloq == 0.5 & resid == "fixed" & endpoint == "AUCinf_true", est], c3$be[lloq == L0 & resid == "fixed" & endpoint == "AUCinf_true", est])   # 참값은 LLOQ와 무관
})

test_that("연구 LLOQ는 config/assay.yaml 한 곳에서 읽고, 민감도 격자는 연구 LLOQ를 포함하며, 옵션으로만 고정할 수 있다", {
  a <- read_cfg("assay.yaml")
  expect_equal(study_lloq(), as.numeric(a$lloq_mg_L$value))
  expect_true(any(abs(study_lloq_grid() - study_lloq()) < 1e-12)); expect_false(is.unsorted(study_lloq_grid()))
  old <- options(dupi.study_lloq = 0.2); on.exit(options(old))
  expect_equal(study_lloq(), 0.2)
  options(old); expect_equal(study_lloq(), as.numeric(a$lloq_mg_L$value))
  expect_null(read_cfg("trial_design.yaml")$lloq_mg_L)                   # 두 번째 출처가 없어야 한다
  expect_equal(load_params("k2016")$lloq, read_cfg("design_clot2021.yaml")$lloq_mg_L)   # p$lloq = 외부 자료(Clot 2021) LLOQ
})
