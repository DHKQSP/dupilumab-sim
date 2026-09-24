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

test_that("시험 수준 대응 비교와 판정 규칙 (a)–(d) (검토 의견 3차 §5)", {
  be <- rbind(data.table(trial = 1:50, schedule = "B0", scenario = "S00", endpoint = "AUClast", method = "pooled_t", width = 0.20, pass = TRUE),
              data.table(trial = 1:50, schedule = "D1", scenario = "S00", endpoint = "AUClast", method = "pooled_t", width = 0.195, pass = TRUE),
              data.table(trial = 1:50, schedule = "B0", scenario = "KE110", endpoint = "AUClast", method = "pooled_t", width = 0.20, pass = rep(c(TRUE, FALSE), 25)),
              data.table(trial = 1:50, schedule = "D1", scenario = "KE110", endpoint = "AUClast", method = "pooled_t", width = 0.20, pass = c(rep(TRUE, 26), rep(FALSE, 24))))
  pt <- rbind(paired_trial_width_vs_ref(be, "B0", "S00"), paired_trial_width_vs_ref(be, "B0", "KE110"))
  expect_equal(pt[scenario == "S00", width_rel_decrease], 0.025)
  per_ep <- be[, .(pass_rate = 100 * mean(pass), width_mean_pp = 100 * mean(width)), by = .(scenario, schedule, endpoint)]
  design <- read_cfg("trial_design.yaml")
  ind <- data.table(schedule = c("B0", "D1"), reliable_pct = c(86, 87), extrap_gt20_pct = c(1.0, 0.4))
  pi <- data.table(schedule = "D1", reliable_gain_pp = 1, extrap_gt20_mcnemar_p = 0.01)
  dec <- schedule_decision(ind, pi, per_ep, pt, design, data.table(schedule = c("B0", "D1"), n_points = c(13, 15), added_points = c(0, 2)))
  expect_equal(nrow(dec), 1L)
  expect_equal(dec$a_mean_width_rel_decrease, 0.025); expect_true(dec$crit_a)     # 2.5% ≥ 2%
  expect_equal(dec$b_ke110_gain_pp, 2, tolerance = 1e-12)                        # 50% → 52% = +2%p
  expect_true(dec$crit_b)
  expect_false(dec$crit_c)                                                        # +1%p < 5%p
  expect_equal(dec$d_extrap20_ratio, 0.4); expect_true(dec$crit_d)                # 0.4 ≤ 0.5
  expect_true(dec$recommend); expect_equal(dec$added_visits_total, 2 * 260); expect_equal(dec$n_criteria_evaluable, 4L)
  # KE110 없는 변형: (b)는 NA, 나머지로 판정
  dec2 <- schedule_decision(ind, pi, per_ep[scenario == "S00"], pt[scenario == "S00"], design, data.table(schedule = c("B0", "D1"), n_points = c(13, 15), added_points = c(0, 2)))
  expect_true(is.na(dec2$crit_b)); expect_equal(dec2$n_criteria_evaluable, 3L)
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
