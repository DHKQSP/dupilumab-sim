# 경계 재판정(검토 의견 W2 §2): AUCinf 규칙 × 플래그 세트 구성(R/oc.R OC_REJUDGE_CONFIGS, config_pass_ext)과
# run_trial_oc_ext(scripts/40_oc_rejudge.R)의 원래 6개 평가변수 재현 검사

test_that("재판정 구성 정의는 사전 고정 구성(P2, G2, AUCinf_only)과 같고 G2 계열은 모두 Cmax를 포함한다", {
  oc <- read_cfg("oc_design.yaml")
  expect_identical(OC_ENDPOINTS_EXT[seq_along(OC_ENDPOINTS)], OC_ENDPOINTS)
  expect_setequal(setdiff(OC_ENDPOINTS_EXT, OC_ENDPOINTS), c("AUCinf_Ai", "AUCinf_Ci"))
  cf <- OC_REJUDGE_CONFIGS
  expect_setequal(cf$P2, unlist(oc$configurations$P2$endpoints))
  expect_setequal(cf$G2_Aii, unlist(oc$configurations$G2$endpoints))
  expect_setequal(cf$AUCinf_Aii, unlist(oc$configurations$AUCinf_only$endpoints))
  g2 <- grep("^G2_", names(cf), value = TRUE)
  expect_true(all(vapply(cf[g2], function(e) "Cmax" %in% e && length(e) == 2, logical(1))))
  expect_true(all(unlist(cf) %in% OC_ENDPOINTS_EXT))
  expect_identical(cf$G2_Ai, c("AUCinf_Ai", "Cmax")); expect_identical(cf$G2_Ci, c("AUCinf_Ci", "Cmax")); expect_identical(cf$G2_B, c("AUCinf_B", "Cmax"))
})

test_that("config_pass_ext는 규칙 × 플래그 세트 조합을 평가변수 통과의 논리곱으로 계산하고, pass NA는 불통과, 행 없음은 미평가(NA)로 둔다", {
  be <- CJ(trial = 1:2, scenario = c("a", "b", "c"), endpoint = OC_ENDPOINTS_EXT)
  be[, pass := TRUE]
  be[trial == 1 & scenario == "b" & endpoint == "AUCinf_A", pass := NA]        # 세트 (ii) 규칙 A 산출 불가, 세트 (i)은 통과
  be[trial == 1 & scenario == "c" & endpoint == "Cmax", pass := FALSE]         # Cmax 불통과 → P2, G2 전부 불통과, AUCinf 단독은 통과
  be[trial == 2 & scenario == "a" & endpoint == "AUCinf_Ci", pass := FALSE]    # 세트 (i) 규칙 C만 불통과
  be[trial == 2 & scenario == "b" & endpoint == "AUCinf_Ai", pass := FALSE]    # 세트 (i) 규칙 A만 불통과
  be <- be[!(trial == 2 & scenario == "c" & endpoint %in% c("AUCinf_Ai", "AUCinf_Ci"))]   # 재판정 행 없음(시험 2, c)
  be[trial == 2 & scenario == "c" & endpoint == "Cmax", pass := FALSE]
  w <- config_pass_ext(be)
  g <- function(tr, sc, cf) w[trial == tr & scenario == sc][[paste0("cfg_", cf)]]
  # 모두 통과
  expect_true(all(unlist(w[trial == 1 & scenario == "a", paste0("cfg_", names(OC_REJUDGE_CONFIGS)), with = FALSE])))
  # (1, b): 규칙 A (ii) 판정 불가 = 불통과, (i)은 통과, B·C는 영향 없음
  expect_false(g(1, "b", "G2_Aii")); expect_false(g(1, "b", "AUCinf_Aii"))
  expect_true(g(1, "b", "G2_Ai")); expect_true(g(1, "b", "AUCinf_Ai")); expect_true(g(1, "b", "G2_B")); expect_true(g(1, "b", "G2_Cii")); expect_true(g(1, "b", "P2"))
  # (1, c): Cmax 불통과
  for (cf in c("P2", "G2_Aii", "G2_B", "G2_Cii", "G2_Ai", "G2_Ci")) expect_false(g(1, "c", cf))
  for (cf in c("AUCinf_Aii", "AUCinf_B", "AUCinf_Cii", "AUCinf_Ai", "AUCinf_Ci")) expect_true(g(1, "c", cf))
  # (2, a): 세트 (i) 규칙 C만
  expect_false(g(2, "a", "G2_Ci")); expect_false(g(2, "a", "AUCinf_Ci")); expect_true(g(2, "a", "G2_Cii")); expect_true(g(2, "a", "G2_Ai"))
  # (2, b): 세트 (i) 규칙 A만
  expect_false(g(2, "b", "G2_Ai")); expect_true(g(2, "b", "G2_Aii")); expect_true(g(2, "b", "G2_Ci"))
  # (2, c): 세트 (i) 행 없음 → Cmax가 불통과여도 NA(미평가). 원래 평가변수 구성은 정상 판정
  expect_true(is.na(g(2, "c", "G2_Ai"))); expect_true(is.na(g(2, "c", "G2_Ci"))); expect_true(is.na(g(2, "c", "AUCinf_Ai")))
  expect_false(g(2, "c", "G2_Aii")); expect_false(g(2, "c", "P2")); expect_true(g(2, "c", "AUCinf_Cii"))
  # 원래 6개 평가변수 구성은 사전 고정 config_pass와 같다
  oc <- read_cfg("oc_design.yaml")
  wc <- config_pass(be[endpoint %in% OC_ENDPOINTS], oc)
  m <- merge(w, wc[, .(trial, scenario, P2_0 = cfg_P2, G2_0 = cfg_G2, AO_0 = cfg_AUCinf_only)], by = c("trial", "scenario"))
  expect_identical(m$cfg_P2, m$P2_0); expect_identical(m$cfg_G2_Aii, m$G2_0); expect_identical(m$cfg_AUCinf_Aii, m$AO_0)
  # 재판정 열이 아예 없으면 세트 (i) 구성은 전부 NA
  w0 <- config_pass_ext(be[endpoint %in% OC_ENDPOINTS])
  expect_true(all(is.na(w0$cfg_G2_Ai)) && all(is.na(w0$cfg_AUCinf_Ci)))
  expect_identical(w0$cfg_G2_Aii, w[match(paste(w0$trial, w0$scenario), paste(trial, scenario)), cfg_G2_Aii])
  # 중복 행은 오류
  expect_error(config_pass_ext(rbind(be, be[1])), "중복")
})

test_that("run_trial_oc_ext는 원래 6개 평가변수를 run_trial_oc와 비트 단위로 재현하고 세트 (i) 포함 수는 세트 (ii) 이상이다", {
  skip_if_no_rxode2()
  skip_if_fast_scope("시험 1회 모의(시나리오 1개), 약 10초")
  oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
  rv <- resolve_variant("base", design); p <- rv$p
  scen <- list(VM = list(code = "VM", T_multipliers = list(Vmax = 0.66)))
  ms <- as.integer(oc$trials$master_seed)
  a <- run_trial_oc(3L, p, design, scen, ms, rv$wt_spec, model_id = p$model_id)
  b <- run_trial_oc_ext(3L, p, design, scen, ms, rv$wt_spec, model_id = p$model_id)
  expect_identical(b$be[endpoint %in% OC_ENDPOINTS], a$be)
  expect_identical(b$be$endpoint, OC_ENDPOINTS_EXT)
  expect_identical(b$drop[, !"n_reliable_i"], a$drop)
  e <- function(ep) b$be[endpoint == ep]
  expect_gte(e("AUCinf_Ai")$n_R, e("AUCinf_A")$n_R); expect_gte(e("AUCinf_Ai")$n_T, e("AUCinf_A")$n_T)
  expect_identical(c(e("AUCinf_Ai")$n_R, e("AUCinf_Ai")$n_T), b$drop$n_reliable_i)   # 규칙 A (i) 포함 수 = 세트 (i) 충족 수(R, T)
  expect_identical(c(e("AUCinf_Ci")$n_R, e("AUCinf_Ci")$n_T), c(e("AUClast")$n_R, e("AUClast")$n_T))   # 규칙 C는 전원
  expect_true(all(b$drop$n_reliable_i >= b$drop$n_reliable) && all(b$drop$n_reliable_i <= b$drop$n_lambda))
})
