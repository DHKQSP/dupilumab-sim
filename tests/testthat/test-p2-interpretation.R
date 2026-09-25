# P2 경계 1종 오류 해석 도우미(R/oc_interpret.R, scripts/43_p2_interpretation.R)의 빠른 검사. 모의 없음.

test_that("classify_type1은 Wilson 구간을 5%와 비교해 보수적·명목·초과로 나누고, 끝점이 5%와 같으면 명목이다", {
  expect_identical(classify_type1(c(3.36, 4.90, 5.01, 4.00, 5.00, NA), c(4.10, 5.78, 5.90, 5.00, 5.50, NA)),
                   c("conservative", "nominal", "exceeding", "nominal", "nominal", NA))
  expect_identical(classify_type1(4, 4.99, alpha_pct = 5), "conservative")
  expect_identical(classify_type1(2.6, 3, alpha_pct = 2.5), "exceeding"); expect_identical(classify_type1(2, 3, alpha_pct = 2.5), "nominal")
  expect_error(classify_type1(5, 4), "하한 > 상한")
  expect_error(classify_type1(1:2, 3), "길이")
  expect_error(classify_type1(c(1, NA), c(2, 3)), "NA")
  # wilson_ci(R/summarize.R)와 함께: 10,000회 중 532회(5.32%)는 구간이 5%를 포함, 371회는 상한 < 5%, 600회는 하한 > 5%
  w <- rbind(wilson_ci(532, 10000), wilson_ci(371, 10000), wilson_ci(600, 10000), wilson_ci(0, 10000))
  expect_identical(classify_type1(w$lo, w$hi), c("nominal", "conservative", "exceeding", "conservative"))
})

test_that("be_boundary_type1_theory: 참값이 정확히 경계면 가까운 한계 확률이 5%이고 통과 = 5% − 반대쪽 한계 확률", {
  q <- qt(0.95, 232); W <- log(1.25) - log(0.80)
  for (se in c(0.03, 0.05, 0.0632, 0.1)) {
    lo <- be_boundary_type1_theory(se, 232, 0.80); hi <- be_boundary_type1_theory(se, 232, 1.25)
    expect_equal(lo$q, q)
    expect_equal(lo$p_near, 0.05, tolerance = 1e-12)
    expect_equal(lo$p_far, pt(W / se - q, 232, lower.tail = FALSE), tolerance = 1e-12)
    expect_equal(lo$p_pass, 0.05 - lo$p_far, tolerance = 1e-15)
    expect_gt(lo$p_far, 0); expect_lte(lo$p_pass, 0.05 + 1e-15)                  # 수학적으로 5% 미만(부동소수 오차 1e-16 수준만 허용)
    # 대칭(log 0.80 = −log 1.25): 1.25에서도 같은 값
    expect_equal(hi$p_near, 0.05, tolerance = 1e-12); expect_equal(hi$p_pass, lo$p_pass, tolerance = 1e-12); expect_equal(hi$p_far, lo$p_far, tolerance = 1e-12)
  }
  # se = 0.05에서 반대쪽 한계 확률은 무시할 만큼 작다(1e-9 미만), se = 0.1이면 눈에 띈다(약 0.3%)
  expect_lt(be_boundary_type1_theory(0.05, 232, 0.80)$p_far, 1e-9)
  expect_gt(be_boundary_type1_theory(0.10, 232, 0.80)$p_far, 1e-3)
})

test_that("be_boundary_type1_theory: 참값에 따른 단조성, 구간이 너무 넓으면 0, 벡터 입력과 입력 검사", {
  p <- be_boundary_type1_theory(0.05, 232, c(0.79, 0.80, 0.81, 1.00, 1.24, 1.25, 1.26))$p_pass
  expect_true(p[1] < p[2] && p[2] < p[3] && p[3] < p[4])            # 0.80 쪽에서 안으로 갈수록 증가
  expect_true(p[7] < p[6] && p[6] < p[5] && p[5] < p[4])            # 1.25 쪽도
  expect_equal(p[4], 2 * pt((log(1.25) - qt(0.95, 232) * 0.05) / 0.05, 232) - 1, tolerance = 1e-12)   # 참값 1: P(|T| < (log 1.25 − q·se)/se)
  expect_gt(p[4], 0.99)
  expect_equal(be_boundary_type1_theory(0.2, 232, 1)$p_pass, 0)      # 2·q·se > log(1.25/0.80): 통과 불가
  v <- be_boundary_type1_theory(c(0.04, 0.05), 232, 0.80)
  expect_equal(nrow(v), 2); expect_equal(v$p_near, c(0.05, 0.05), tolerance = 1e-12)
  expect_equal(nrow(be_boundary_type1_theory(0.05, c(200, 232), 0.80)), 2)
  expect_error(be_boundary_type1_theory(0, 232, 0.8), "se")
  expect_error(be_boundary_type1_theory(0.05, 232, -1), "true_ratio")
  expect_error(be_boundary_type1_theory(0.05, 232, 0.8, limits = c(1.25, 0.8)), "limits")
})

test_that("be_pass_prob_normal은 SD = se, df → ∞에서 이론값과 같고, 1 쪽 편향은 올리고 SD < se는 내린다", {
  se <- 0.05
  th <- be_boundary_type1_theory(se, 1e7, 0.80)$p_pass
  expect_equal(be_pass_prob_normal(log(0.80), se, se, 1e7), th, tolerance = 1e-6)
  base <- be_pass_prob_normal(log(0.80), se, se, 232)
  expect_gt(be_pass_prob_normal(log(0.80) + 0.004, se, se, 232), base)  # 0.80에서 양의 편향 = 1 쪽
  expect_lt(be_pass_prob_normal(log(0.80) - 0.004, se, se, 232), base)  # 1에서 멀어짐
  expect_lt(be_pass_prob_normal(log(0.80), 0.965 * se, se, 232), 0.05)  # 시험 간 SD가 시험 안 se보다 작으면 경계 통과 < 5%
  expect_gt(be_pass_prob_normal(log(1.25) - 0.004, se, se, 232), be_pass_prob_normal(log(1.25), se, se, 232))   # 1.25에서 음의 편향 = 1 쪽
  expect_error(be_pass_prob_normal(0, 0, 0.05, 232), "sd_log")
})

test_that("bias_direction: 참값 < 1이면 음이 1에서 멀어짐, 참값 > 1이면 양이 멀어짐, 95% 구간이 0을 포함하면 none", {
  expect_identical(bias_direction(c(-0.005, 0.004, -0.005, 0.009), c(0.8, 0.8, 1.25, 1.25)), c("away_from_1", "toward_1", "toward_1", "away_from_1"))
  expect_identical(bias_direction(0.0002, 0.8, lo = -0.0007, hi = 0.0011), "none")
  expect_identical(bias_direction(0.004, 0.8, lo = 0.003, hi = 0.005), "toward_1")
  expect_error(bias_direction(0.01, 1), "1이 아닌")
})
