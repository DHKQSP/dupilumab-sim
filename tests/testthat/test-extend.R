# 적응적 연장(지시 2026-09-25 §3, scripts/42_oc_extend.R): 선택 규칙(R/oc.R oc_extension_select)과 시험 수 문구(oc_n_trials_text). 모의 없음.

bt_synth <- function() {
  x <- data.table(model = c("m1", "m1", "m1", "m1", "m1", "m2", "m2"),
                  scenario = c("V2_up_080", "F_down_080", "Vmax_down_125", "ke_up_080", "ka_down_080", "V2_up_080", "F_up_125"),
                  mechanism = c("V2", "F", "Vmax", "ke", "ka", "V2", "F"), direction = c("up", "down", "down", "up", "down", "up", "up"),
                  target = c(0.8, 0.8, 1.25, 0.8, 0.8, 0.8, 1.25), multiplier = c(2.28, 0.87, 0.65, 1.46, 0.42, 3.21, 1.15),
                  config = "P2", n_trials = 10000L,
                  pass_pct = c(5.32, 3.70, 5.25, 4.60, 0, 5.60, 4.50),
                  lo =       c(4.90, 3.35, 5.00, 4.20, 0, 5.10, 4.10),
                  hi =       c(5.78, 4.09, 5.50, 5.00, 0.04, 6.10, 4.95))
  # 다른 구성(G2) 행은 규칙에 쓰지 않는다: 구간이 5%를 포함해도 선택하지 않는다
  rbind(x, copy(x[1])[, `:=`(config = "G2", scenario = "F_down_080", pass_pct = 5, lo = 4.6, hi = 5.4)])
}

test_that("연장 규칙: 사전 고정 수에서 P2 Wilson 95% 구간이 5%를 포함하는(lo ≤ 5 ≤ hi, 경계 포함) 시나리오만 목표 수로 늘린다", {
  d <- oc_extension_select(bt_synth(), n_from = 10000L, n_to = 20000L)
  expect_identical(names(d), c("model", "scenario", "mechanism", "direction", "target", "multiplier", "config", "n_before", "pass_pct", "lo", "hi",
                               "threshold", "rule", "selected", "n_after"))
  expect_equal(nrow(d), 7L)                                          # P2 행만(모델 × 경계 시나리오)
  expect_true(all(d$config == "P2") && all(d$n_before == 10000L) && all(d$threshold == 5))
  sel <- d[selected == TRUE, paste(model, scenario)]
  # 포함: 4.90–5.78, 하한이 정확히 5.00, 상한이 정확히 5.00. 제외: 구간 전체가 5% 아래(3.35–4.09, 4.10–4.95, 0–0.04)나 위(5.10–6.10)
  expect_setequal(sel, c("m1 V2_up_080", "m1 Vmax_down_125", "m1 ke_up_080"))
  expect_false(d[model == "m2" & scenario == "V2_up_080", selected])   # 점추정 > 5%여도 구간이 5%를 포함하지 않으면 연장하지 않는다
  expect_identical(d[selected == TRUE, unique(n_after)], 20000L); expect_identical(d[selected == FALSE, unique(n_after)], 10000L)
  expect_identical(d$selected, d$lo <= 5 & d$hi >= 5)               # scripts/20_products_5000.R의 cover 조건과 같은 식
  expect_match(d$rule[1], "P2 .*10000 trials includes 5% .*20000 trials")
  # 판정 값은 입력 그대로(사전 고정 수의 통과율·구간)
  expect_equal(d[model == "m1" & scenario == "V2_up_080", c(pass_pct, lo, hi)], c(5.32, 4.90, 5.78))
})

test_that("연장 규칙은 시나리오를 코드에 적지 않고 표에서 정하며, 입력이 규칙의 전제와 어긋나면 중단한다", {
  b <- bt_synth()
  b2 <- copy(b)[config == "P2" & scenario == "F_down_080" & model == "m1", `:=`(lo = 4.95, hi = 5.30)]   # 다른 시나리오로 바꾸면 그것이 선택된다
  expect_true(oc_extension_select(b2, 10000L, 20000L)[model == "m1" & scenario == "F_down_080", selected])
  none <- copy(b)[config == "P2", `:=`(lo = 1, hi = 2)]
  d0 <- oc_extension_select(none, 10000L, 20000L)
  expect_false(any(d0$selected)); expect_true(all(d0$n_after == 10000L))
  expect_error(oc_extension_select(copy(b)[scenario == "V2_up_080" & model == "m1" & config == "P2", n_trials := 20000L], 10000L, 20000L), "사전 고정")   # 연장 뒤 표로 다시 판정하지 않는다
  expect_error(oc_extension_select(rbind(b, b[1]), 10000L, 20000L), "중복")
  expect_error(oc_extension_select(b, 10000L, 10000L), "목표 시험 수")
  expect_error(oc_extension_select(b[config != "P2"], 10000L, 20000L), "구성 P2")
  expect_error(oc_extension_select(copy(b)[1, lo := NA], 10000L, 20000L), "NA")
  expect_error(oc_extension_select(b[, !"hi"], 10000L, 20000L), "hi")
})

test_that("시험 수 문구는 최빈 수와 예외(적응적 연장, 진행 중)를 실제 수로 적는다", {
  dec <- oc_extension_select(bt_synth(), 10000L, 20000L)
  x <- bt_synth()[config == "P2" & model == "m1"]
  expect_identical(oc_n_trials_text(x, "en", dec), "10,000 each")
  expect_identical(oc_n_trials_text(x, "ko", dec), "각 10,000회")
  x[scenario == "V2_up_080", n_trials := 20000L]
  expect_identical(oc_n_trials_text(x, "en", dec), "10,000 each (V2 up 0.80: 20,000 by the adaptive extension rule)")
  expect_identical(oc_n_trials_text(x, "ko", dec), "각 10,000회(V2 상향 0.80: 적응적 연장 규칙으로 20,000회)")
  expect_identical(oc_n_trials_text(x, "en", dec, c(m1 = "Model 1")), "10,000 each (Model 1 V2 up 0.80: 20,000 by the adaptive extension rule)")
  x[scenario == "V2_up_080", n_trials := 13500L]
  expect_identical(oc_n_trials_text(x, "en", dec), "10,000 each (V2 up 0.80: 13,500 so far, adaptive extension toward 20,000)")
  expect_identical(oc_n_trials_text(x, "ko", dec), "각 10,000회(V2 상향 0.80: 적응적 연장 진행 중 13,500회(목표 20,000회))")
  # 연장 판정이 없거나 선택되지 않은 시나리오의 다른 수는 이유 없이 수만 적는다
  expect_identical(oc_n_trials_text(x, "en"), "10,000 each (V2 up 0.80: 13,500)")
  y <- bt_synth()[config == "P2" & model == "m1"][scenario == "F_down_080", n_trials := 9500L]
  expect_identical(oc_n_trials_text(y, "en", dec), "10,000 each (F down 0.80: 9,500)")
  # 영문 문구: 한글·em dash·en dash·U+2212 없음
  en <- c(oc_n_trials_text(x, "en", dec), oc_n_trials_text(x, "en", dec, c(m1 = "Model 1")), dec$rule)
  expect_false(any(grepl("[가-힣]|—|–|−", en)))
})

test_that("시험 수 문구의 예외만 모드", {
  dec <- oc_extension_select(bt_synth(), 10000L, 20000L)
  x <- bt_synth()[config == "P2" & model == "m1"]
  expect_identical(oc_n_trials_text(x, "en", dec, exceptions_only = TRUE), "")
  x[scenario == "V2_up_080", n_trials := 20000L]
  expect_identical(oc_n_trials_text(x, "en", dec, exceptions_only = TRUE), "V2 up 0.80: 20,000 by the adaptive extension rule")
  expect_identical(oc_n_trials_text(x, "ko", dec, exceptions_only = TRUE), "V2 상향 0.80: 적응적 연장 규칙으로 20,000회")
})
