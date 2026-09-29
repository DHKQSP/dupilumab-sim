# deck_facts_common.R: 여러 슬라이드가 쓰는 설계 상수. 모두 config에서 읽고, 부르는 슬라이드의 출처로 기록한다.
# 문구 파일에 숫자를 쓰지 않도록 설계값(용량, arm당 인원, 체중 범위, 채혈일, LLOQ, 동등성 한계, 기준 세트 정의)도 여기서 만든다.
num_fmt <- function(d = 0) function(x) fnum(as.numeric(x), d)
f_dose <- function() dcfg("trial_design.yaml", "dose_mg", "dose (mg)", num_fmt(0))
f_n_arm <- function() dcfg("trial_design.yaml", "n_per_arm", "evaluable subjects per arm", num_fmt(0))
f_n_rand <- function() dcfg("trial_design.yaml", "n_randomized_per_arm", "randomized subjects per arm", num_fmt(0))
f_wt_range <- function() dcfg("trial_design.yaml", c("weight", "inclusion_kg"), "body weight inclusion range (kg)", function(x) rng_fmt(x[1], x[2], 0))
f_split <- function() dcfg("trial_design.yaml", c("stratification", "split_kg", "value"), "randomization stratum split (kg)", num_fmt(0))
f_lloq <- function() dcfg("assay.yaml", c("lloq_mg_L", "value"), "study LLOQ (mg/L)", num_fmt(3))
f_lloq_grid <- function() dcfg("assay.yaml", "lloq_sensitivity_mg_L", "LLOQ sensitivity range (mg/L)", function(x) rng_fmt(min(x), max(x), 2))
f_limits <- function() dcfg("trial_design.yaml", c("be", "limits"), "equivalence limits (%)", function(x) rng_fmt(100 * x[1], 100 * x[2], 2))
f_ci_level <- function() dcfg("trial_design.yaml", c("be", "ci_level"), "confidence level of the equivalence CI (%)", function(x) fnum(100 * x, 0))
# 명목 1종 오류 = (1 - CI 수준) / 2 (두 단측 검정)
f_nominal <- function() { y <- .read("config/trial_design.yaml")$be$ci_level; p <- paste0(fnum(100 * (1 - y) / 2, 0), "%")
  dderived("nominal type I error per one-sided test (%) = (1 - be.ci_level) / 2", "config/trial_design.yaml", "be.ci_level :: (1 - x) / 2 x 100", 100 * (1 - y) / 2, p) }
# 채혈일: config는 투여 후 일(days)로 적는다. 연구일 = 투여 후 일 + 1(Day 1 = 투여일)
f_study_days <- function(schedule = "B0", which = c("all", "last", "n")) {
  which <- match.arg(which); d <- unlist(.read("config/trial_design.yaml")$schedules[[schedule]]$days); premise(length(d) > 0, paste("schedule", schedule))
  sd <- d + 1; p <- switch(which, all = paste(fnum(sd[sd == round(sd)], 0), collapse = ", "), last = fnum(max(sd), 0), n = as.character(length(d)))
  dderived(sprintf("schedule %s, %s (study day = days after dose + 1)", schedule, which), "config/trial_design.yaml", sprintf("schedules.%s.days :: %s", schedule, which), sd, p)
}
# 기준 세트 정의(prereg section4)
f_set <- function(s_, what = c("r2", "extrap", "span")) {
  what <- match.arg(what); key <- switch(what, r2 = "adj_r2_min", extrap = "extrap_max_pct", span = "span_ratio_min")
  dcfg("prereg_20260926.yaml", c("section4", "criteria_sets", s_, key), sprintf("criteria set (%s) %s", s_, key), if (what == "r2") num_fmt(2) else num_fmt(0))
}
f_reps <- function(key = c("boundary", "ext")) {
  key <- match.arg(key)
  if (key == "boundary") dcfg("oc_design.yaml", c("trials", "reps_boundary"), "trials per boundary scenario (pre-specified)", function(x) fnum(as.numeric(x), 0, big = TRUE))
  else dcfg("prereg_20260926.yaml", c("section1", "extension", "to"), "trials after extension", function(x) fnum(as.numeric(x), 0, big = TRUE))
}
