# config 스키마 검증 (검토 의견 3차 §4, D-021/D-025): YAML 불리언 예약어 키 금지 + 주요 필드 자료형
FORBIDDEN_KEYS <- c("y", "n", "yes", "no", "on", "off", "true", "false")
cfg_files <- list.files(proj_path("config"), pattern = "\\.ya?ml$", full.names = TRUE)

test_that("모든 config에 YAML 불리언 예약어 키가 없다 (원문 검사 + 파싱 결과 검사)", {
  pat <- sprintf("(^|[{,[:space:]])(%s)[[:space:]]*:", paste(c(FORBIDDEN_KEYS, toupper(FORBIDDEN_KEYS), tools::toTitleCase(FORBIDDEN_KEYS)), collapse = "|"))
  for (f in cfg_files) {
    txt <- sub("#.*$", "", readLines(f, warn = FALSE, encoding = "UTF-8"))   # 주석 제거
    txt <- gsub('"[^"]*"', '""', txt)                                          # 따옴표 문자열 제거
    hits <- grep(pat, txt, perl = TRUE)
    expect_length(hits, 0); if (length(hits)) message(basename(f), ": ", paste(hits, collapse = ","))
    walk <- function(x) { if (is.list(x)) { expect_false(any(names(x) %in% c("TRUE", "FALSE")), info = basename(f)); lapply(x, walk) }; invisible(NULL) }
    walk(yaml::read_yaml(f, fileEncoding = "UTF-8"))
  }
})

num1 <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
test_that("파라미터 config 자료형", {
  for (f in c("params_typical.yaml", "params_k2020_model1.yaml")) {
    th <- read_cfg(f)$theta
    for (nm in names(th)) expect_true(num1(th[[nm]]$value) && th[[nm]]$value > 0, info = paste(f, nm))
  }
  v <- read_cfg("params_variability.yaml")
  for (nm in names(v$iiv_omega2)) expect_true(num1(v$iiv_omega2[[nm]]$omega2) && v$iiv_omega2[[nm]]$omega2 >= 0, info = nm)
  expect_true(num1(v$residual$sigma_prop$value) && num1(v$residual$sigma_add$value))
})

test_that("시험 설계 config 자료형", {
  d <- read_cfg("trial_design.yaml")
  expect_true(num1(d$n_per_arm) && d$n_per_arm == round(d$n_per_arm)); expect_true(num1(d$dose_mg)); expect_null(d$lloq_mg_L); expect_true(num1(.num(read_cfg("assay.yaml")$lloq_mg_L)))   # 연구 LLOQ는 assay.yaml 한 곳(D-054)
  for (nm in names(d$schedules)) { x <- unlist(d$schedules[[nm]]$days); expect_true(is.numeric(x) && all(x > 0) && !is.unsorted(x) && !anyDuplicated(x), info = nm) }
  for (w in c("base", "sensitivity")) { x <- d$weight[[w]]; expect_true(num1(x$mean) && num1(x$sd) && is.numeric(x$trunc) && length(x$trunc) == 2 && x$trunc[1] < x$trunc[2], info = w) }
  expect_true(num1(d$stratification$split_kg$value))
  dr <- d$decision_rule
  for (k in c("ci_width_mean_rel_decrease_min", "ke110_pass_gain_pp_min", "reliability_gain_pp_min", "extrap_gt20_ratio_max")) expect_true(num1(dr[[k]]), info = k)
  expect_true(all(d$schedule_analysis %in% names(d$schedules)))
})

test_that("시나리오·gate·검증 config 자료형", {
  s <- read_cfg("scenarios.yaml")
  for (nm in names(s$scenarios)) for (m in names(s$scenarios[[nm]]$T_multipliers)) expect_true(num1(s$scenarios[[nm]]$T_multipliers[[m]]) && s$scenarios[[nm]]$T_multipliers[[m]] > 0, info = paste(nm, m))
  expect_true(all(s$schedule_analysis_scenarios %in% names(s$scenarios)))
  g <- read_cfg("gate_decision.yaml")
  expect_true(num1(g$decision_option) && g$decision_option %in% 1:3)
  expect_true(g$primary_model %in% names(g$model_code_map)); expect_true(all(g$sensitivity_models %in% names(g$model_code_map)))
  dz <- read_cfg("design_clot2021.yaml")
  for (ds in dz$datasets) {
    expect_true(ds$gate_role %in% c("gate", "external"), info = ds$id)
    expect_true(num1(ds$dose_mg) && num1(ds$auclast_mean) && num1(ds$auclast_sd) && num1(ds$weight_mean), info = ds$id)
    expect_true(is.null(ds$n_subj) || num1(ds$n_subj), info = ds$id)
  }
  expect_true(all(vapply(dz$groups, function(g) g$gate_role %in% c("gate", "external"), logical(1))))
})
