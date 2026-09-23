test_that("개인 수준 대응 비교: McNemar 불일치 쌍과 신뢰 충족률 차이", {
  n <- 100
  base <- data.table(id = 1:n, schedule = "B0", reliable = rep(c(TRUE, FALSE), each = 50), pct_extrap = c(rep(25, 10), rep(5, 90)),
                     lambda_ok = TRUE, pct_extrap_true = 1, err_AUClast_vs_true_inf = rnorm(n, 0, 0.02))
  d1 <- copy(base)[, schedule := "D1"]
  d1[id %in% 51:60, reliable := TRUE]          # 10명 신뢰 충족으로 전환
  d1[id %in% 1:8, pct_extrap := 15]             # 8명 >20% 해소
  d1[id %in% 90, pct_extrap := 22]              # 1명 새로 >20%
  r <- paired_individual_vs_ref(rbind(base, d1), "B0")
  expect_equal(r$reliable_gain_pp, 10)
  expect_equal(r$extrap_gt20_only_ref, 8L); expect_equal(r$extrap_gt20_only_sched, 1L)
  expect_equal(r$extrap_gt20_change_pp, -7)
  expect_equal(r$extrap_gt20_mcnemar_p, binom.test(1, 9, 0.5)$p.value)
  expect_equal(r$reliable_mcnemar_p, binom.test(0, 10, 0.5)$p.value)
})

test_that("시험 수준 대응 비교: CI 폭 상대 감소와 판정 규칙", {
  be <- rbind(data.table(trial = 1:50, schedule = "B0", scenario = "S00", endpoint = "AUClast", method = "pooled_t", width = 0.20, pass = TRUE),
              data.table(trial = 1:50, schedule = "D1", scenario = "S00", endpoint = "AUClast", method = "pooled_t", width = 0.18, pass = TRUE))
  pt <- paired_trial_width_vs_ref(be, "B0", "S00")
  expect_equal(pt$width_rel_decrease, 0.1)
  design <- read_cfg("trial_design.yaml")
  pi <- data.table(schedule = "D1", reliable_gain_pp = 1, extrap_gt20_change_pp = -0.1, extrap_gt20_mcnemar_p = 0.5)
  dec <- schedule_decision(pi, pt, design, data.table(schedule = "D1", n_points = 15, added_points = 2))
  expect_true(dec$crit_a_width); expect_false(dec$crit_b_reliable); expect_false(dec$crit_c_extrap); expect_true(dec$recommend)
  expect_equal(dec$added_visits_total, 2 * 260)
})

test_that("시나리오·방법이 섞인 결과에서도 필터가 정확하다 (D-022 회귀 테스트)", {
  mk <- function(sc, sh, m, w, ps) data.table(trial = 1:20, scenario = sc, schedule = sh, endpoint = "AUClast", method = m, width = w, pass = ps, GMR = 1, CI_lower = 0.9, CI_upper = 1.1, n_R = 117, n_T = 117)
  be <- rbind(mk("S00", "B0", "pooled_t", 0.20, TRUE), mk("S00", "D1", "pooled_t", 0.18, TRUE),
              mk("S00", "B0", "ancova_weight", 0.10, FALSE), mk("S00", "D1", "ancova_weight", 0.10, FALSE),
              mk("F090", "B0", "pooled_t", 0.30, FALSE), mk("F090", "D1", "pooled_t", 0.30, FALSE))
  pt <- paired_trial_width_vs_ref(be, "B0", "S00")
  expect_equal(nrow(pt), 1L); expect_equal(pt$width_rel_decrease, 0.1); expect_equal(pt$scenario, "S00")
  st <- summarize_trials(be[endpoint == "AUClast"], method = "pooled_t")
  expect_equal(st$per_endpoint[scenario == "S00" & schedule == "B0", pass_rate], 100)       # ANCOVA 행(불통과)이 섞이면 50
  expect_equal(st$per_endpoint[scenario == "S00" & schedule == "B0", width_mean_pp], 20)
  st2 <- summarize_trials(be, method = "ancova_weight")
  expect_equal(st2$per_endpoint[scenario == "S00" & schedule == "B0", pass_rate], 0)
})
