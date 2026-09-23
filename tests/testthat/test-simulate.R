skip_if_no_rxode2()
p <- load_params("k2016"); design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()

test_that("잔차·BLQ 규칙과 시드 재현", {
  ip <- typical_subject(p)
  planned <- c(0, 0.25, 1, 7, 28, 56, 120)
  obs <- CJ(id = 1L, planned = planned)[, time := planned]
  eps0 <- draw_eps(1L, planned, c(prop = 0, add = 0))
  s <- simulate_observations(ip, obs, 300, p$lloq, eps0)
  expect_equal(s$obs$y_raw, s$obs$C); expect_true(s$obs[planned == 0, blq]); expect_true(s$obs[planned == 120, blq])
  expect_true(all(is.na(s$obs[blq == TRUE, conc]))); expect_true(s$truth$AUCinf_true > 0)
  e1 <- with_seed(9, draw_eps(1L, planned, p$sigma)); e2 <- with_seed(9, draw_eps(1L, planned, p$sigma)); expect_identical(e1, e2)
})

test_that("일정 부분집합과 참값 부착", {
  ip <- typical_subject(p); grid <- union_grid(design, design$schedule_analysis)
  obs <- CJ(id = 1L, planned = grid)[, time := planned]; eps <- draw_eps(1L, grid, c(prop = 0, add = 0))
  s <- simulate_observations(ip, obs, 300, p$lloq, eps)
  ob <- subset_schedule(s$obs, get_schedule(design, "B0")); expect_equal(nrow(ob), 14)
  nca <- attach_truth(run_nca(ob), ob, s$truth)
  expect_equal(nca$AUClast_true, ob[time == nca$tlast, auc]); expect_true(nca$coverage_true <= 1)
})

test_that("시험 엔진: 대조군은 시나리오 간 동일(공통 난수), 시험군은 같은 eta로 배율만 다름, 각 arm 정확히 n", {
  scen <- sc$scenarios[c("S00", "F090")]
  combos <- CJ(scenario = c("S00", "F090"), schedule = c("B0", "D1"))
  r <- run_trial(1L, p, design, scen, combos, master_seed = 7L, wt_spec = weight_spec_from_design(design, "base"), n_per_arm = NULL, keep_nca = TRUE)
  expect_equal(sort(unique(r$be$endpoint)), sort(names(BE_ENDPOINTS)))
  expect_equal(nrow(r$be), nrow(combos) * length(BE_ENDPOINTS))
  n1 <- r$nca[scenario == "S00" & schedule == "B0"]; n2 <- r$nca[scenario == "F090" & schedule == "B0"]
  expect_equal(n1[arm == "R", AUClast], n2[arm == "R", AUClast])                      # 대조군 동일
  expect_equal(sum(n1$arm == "R"), design$n_per_arm); expect_equal(sum(n1$arm == "T"), design$n_per_arm)
  ratio <- n2[arm == "T", AUClast] / n1[arm == "T", AUClast]
  ratio_true <- n2[arm == "T", AUCinf_true] / n1[arm == "T", AUCinf_true]
  # F×0.9: MM 포화 소실 때문에 AUC는 10%보다 더 줄어든다(대표 개체 참값 0.846). AUClast는 같은 개체의 참값 비를 따라야 한다.
  expect_lt(median(ratio_true), 0.9)
  expect_equal(median(ratio), median(ratio_true), tolerance = 0.03)
  expect_true(cor(log(n1[arm == "T", AUClast]), log(n2[arm == "T", AUClast])) > 0.99)   # 같은 eta·잔차
  r2 <- run_trial(1L, p, design, scen, combos, master_seed = 7L, wt_spec = weight_spec_from_design(design, "base"))
  expect_equal(r$be$GMR, r2$be$GMR)                                                    # 재현성
})
