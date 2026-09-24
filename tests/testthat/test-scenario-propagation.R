# 시나리오 배율이 시험 엔진의 개체 파라미터 표에 실제 반영되는지 (검토 의견 3차 §4)
skip_if_no_rxode2()
p <- load_params("k2016"); design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
gm <- function(x) exp(mean(log(x)))

test_that("모든 제품 시나리오: 시험군 개체 파라미터 기하평균이 해당 배율만큼만 이동하고 대조군은 불변", {
  combos <- CJ(scenario = names(sc$scenarios), schedule = "B0")
  sa <- simulate_trial_arms(1L, p, design, sc$scenarios, "B0", 11L, weight_spec_from_design(design, "base"))
  pars <- c("Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "F")
  cfg_name <- c(Vc_i = "Vc", ke = "ke", k12 = "k12", k21 = "k21", Vmax = "Vmax", Km = "Km", ka = "ka", F = "F")
  base_T <- sa$arms$T$S00$ipar
  for (s_ in names(sc$scenarios)) {
    ipT <- sa$arms$T[[s_]]$ipar
    expect_identical(ipT$id, base_T$id)
    mult <- sc$scenarios[[s_]]$T_multipliers
    for (par in pars) {
      expected <- if (cfg_name[[par]] %in% names(mult)) as.numeric(mult[[cfg_name[[par]]]]) else 1
      expect_equal(gm(ipT[[par]]) / gm(base_T[[par]]), expected, tolerance = 1e-10, info = paste(s_, par))
    }
  }
  expect_identical(names(sa$arms$R), "S00")
})

test_that("양 군 공통 배율(곡률 민감도)은 대조군과 시험군 모두에 적용된다", {
  v <- resolve_variant("vmax080_both", design, sc)
  sa0 <- simulate_trial_arms(1L, p, design, sc$scenarios["S00"], "B0", 11L, weight_spec_from_design(design, "base"))
  sa1 <- simulate_trial_arms(1L, v$p, design, sc$scenarios["S00"], "B0", 11L, v$wt_spec)
  for (a in c("R", "T")) {
    expect_equal(gm(sa1$arms[[a]]$S00$ipar$Vmax) / gm(sa0$arms[[a]]$S00$ipar$Vmax), 0.8, tolerance = 1e-10)
    expect_equal(gm(sa1$arms[[a]]$S00$ipar$ke) / gm(sa0$arms[[a]]$S00$ipar$ke), 1, tolerance = 1e-10)
  }
})
