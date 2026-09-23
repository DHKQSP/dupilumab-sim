design <- read_cfg("trial_design.yaml")

test_that("체중은 포함기준 안에 있고 성별 혼합이 반영된다", {
  w <- with_seed(1, draw_weights(2000, design))
  lo <- design$population$weight_kg$inclusion$min; hi <- design$population$weight_kg$inclusion$max
  expect_true(all(w$WT >= lo & w$WT <= hi))
  expect_equal(mean(w$sex == "M"), design$population$sex_ratio_male$value, tolerance = 0.05)
  expect_gt(mean(w[sex == "M", WT]), mean(w[sex == "F", WT]))
})

test_that("eta의 표본 SD는 omega와 일치하고 omega=0이면 0이다", {
  om <- c(Vc = 0.3, ke = 0.2, k12 = 0, k21 = 0, F = 0.4, Vmax = 0.1, ka = 0, Km = 0)
  e <- with_seed(2, draw_etas(5000, om))
  expect_equal(sd(e$eta_Vc), 0.3, tolerance = 0.03)
  expect_equal(sd(e$eta_F), 0.4, tolerance = 0.03)
  expect_true(all(e$eta_k12 == 0) && all(e$eta_Km == 0))
})

test_that("개체 파라미터: F는 (0,1), 체중 공변량은 지수식대로", {
  p <- load_params("dev")
  e <- data.table(id = 1:3, WT = c(50, p$cov$WT_ref, 100), eta_F = c(-3, 0, 3))
  ip <- individual_params(p, e)
  expect_true(all(ip$F > 0 & ip$F < 1))
  expect_equal(ip$F[2], p$theta[["F"]])
  expect_equal(ip$Vc_i, p$theta[["Vc"]] * (e$WT / p$cov$WT_ref)^p$cov$theta_WT)
  expect_equal(unique(ip$Km), p$theta[["Km"]])
})

test_that("체중 층 분류", {
  s <- weight_stratum(c(50, 65, 79.9, 80, 99), design)
  expect_equal(as.character(s), c("<65", "65-80", "65-80", ">=80", ">=80"))
})
