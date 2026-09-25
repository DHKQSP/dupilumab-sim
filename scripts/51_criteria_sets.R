#!/usr/bin/env Rscript
# 51_criteria_sets.R — λz 신뢰 기준 관행값 민감도 요약 (두 번째 추가 지시 2026-09-26 §1, prereg section4, D-057). 새 모의 없음.
# 입력: 저장된 20,000명 B0 NCA(results/individual/nca_base_20000.rds, nca_struct2020_20000.rds), 아토피 모집단 요약(results/atopic/atopic_individual_criteria.csv),
#       재생성 시험(results/oc_models/oc_models_be_<model>.csv.gz, oc_models_crit_be_<model>.csv.gz, 연장 파일), 재시작 전 부분 산출(run1_partial_*).
# 산출(results/criteria/):
#   restart_identity_check.csv        section1 재시작 전후 같은 시험의 행 일치
#   criteria_individual.csv           모집단(60–90 kg 두 모델, 아토피 3 변형 × 3 분포) × 세트: 미달(λz 산출 불가 / 산출되나 미달), arm당 인원, 미달 대 잔류 비교
#   criteria_g2_type1.csv             경계 16칸 × 분석 모형 × G2 변형(규칙 A·C × 세트 i–iv, 규칙 B): 통과율, Wilson, 분류
#   criteria_check_vs_v10.csv         M0의 기존 변형 5개가 v1.0 저장값(oc/g2_rules_flags.csv)과 같은지(시험 수가 같은 칸)
#   criteria_g2_power.csv             S00·F097 검정력
#   criteria_bias.csv                 칸 × 분석 모형 × AUC0-inf 평가변수: 참 AUC0-inf 비 대비 편향
#   criteria_instability.csv          칸 × 분석 모형: 규칙·세트만 바꿀 때 판정이 뒤집히는 시험 비율
#   fig_criteria_failure.png, fig_criteria_g2.png, criteria_conclusion_en.md, criteria_conclusion_ko.md
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
pr <- read_cfg("prereg_20260926.yaml"); p4 <- pr$section4; oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", "criteria"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
om_dir <- proj_path("results", "oc_models"); oc_dir <- proj_path("results", "oc")
for (s_ in names(CRIT_SETS)) { a <- CRIT_SETS[[s_]]; b <- p4$criteria_sets[[s_]]
  stopifnot(isTRUE(all.equal(a$adj_r2_min, b$adj_r2_min)), isTRUE(all.equal(a$extrap_max_pct, as.numeric(b$extrap_max_pct))), identical(is.na(a$span_ratio_min), is.null(b$span_ratio_min))) }
SETS <- names(CRIT_SETS); PK <- c("k2016", "k2020"); n_bnd <- as.integer(oc$trials$reps_boundary)
MODEL_EN <- c(k2016 = "2016 model", k2020 = "Model 1", base = "2016 model", struct2020 = "Model 1", k2016_bmi_vc0817 = "2016 model + BMI and weight covariates")
f1 <- function(x) formatC(round(x + sign(x) * 1e-9, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x + sign(x) * 1e-9, 2) + 0, format = "f", digits = 2); f3 <- function(x) formatC(round(x + sign(x) * 1e-9, 3) + 0, format = "f", digits = 3)
rg <- function(x, f = f1, u = "%", sep = " to ") sprintf("%s%s%s%s%s", f(min(x)), u, sep, f(max(x)), u)
mm <- function(x) sprintf("a median of %s%% and up to %s%%", f1(median(x)), f1(max(x))); mmk <- function(x) sprintf("%s%%/%s%%", f1(median(x)), f1(max(x)))

# ---- 0) 재시작 전후 일치(section1) ----------------------------------------------------------------------------------------
RC <- rbindlist(lapply(PK, function(m) {
  f0 <- file.path(om_dir, sprintf("run1_partial_oc_models_be_%s.csv.gz", m)); if (!file.exists(f0)) return(NULL)
  a <- fread(f0); cn <- a[, .N, by = trial]; full <- cn[N == max(N), trial]; a <- a[trial %in% full]
  b <- fread(file.path(om_dir, sprintf("oc_models_be_%s.csv.gz", m)))[trial %in% full]
  mm <- merge(a, b, by = c("trial", "scenario", "endpoint", "model"), suffixes = c(".run1", ".run2"))
  data.table(pk_model = m, trials_compared = length(full), rows_run1 = nrow(a), rows_matched = nrow(mm),
             identical_rows = sum(mm$est.run1 == mm$est.run2 & mm$se.run1 == mm$se.run2 & mm$pass.run1 == mm$pass.run2 & mm$n_R.run1 == mm$n_R.run2 & mm$n_T.run1 == mm$n_T.run2, na.rm = TRUE) +
               sum(is.na(mm$est.run1) & is.na(mm$est.run2)))
}))
if (nrow(RC)) { RC[, ok := rows_run1 == rows_matched & identical_rows == rows_matched]; fwrite(RC, file.path(out_dir, "restart_identity_check.csv")); if (!all(RC$ok)) stop("재시작 전후 section1 행 불일치") }

# ---- 1) 개인 수준: 60–90 kg(저장본)과 아토피 -----------------------------------------------------------------------------
I6 <- rbindlist(lapply(c(k2016 = "base", k2020 = "struct2020"), function(v) {
  x <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", v)))[schedule == "B0"]
  crit_summary(x)[, `:=`(population = "60-90 kg (study inclusion range)", variant = v, distribution = "study_60_90", band = "all")]
}))
AT <- fread(proj_path("results", "atopic", "atopic_individual_criteria.csv"))[band == "all"][, population := "adult atopic dermatitis"]
CI <- rbind(I6, AT, fill = TRUE)
setcolorder(CI, c("population", "variant", "distribution", "band", "set"))
fwrite(CI, file.path(out_dir, "criteria_individual.csv"))

# ---- 2) 시험 수준: G2(규칙 × 세트), 편향, 판정 불안정성 ----------------------------------------------------------------------
G2V <- c(A_i = "AUCinf_Ai", A_ii = "AUCinf_A", A_iii = "AUCinf_Aiii", A_iv = "AUCinf_Aiv", B = "AUCinf_B", C_i = "AUCinf_Ci", C_ii = "AUCinf_C", C_iii = "AUCinf_Ciii", C_iv = "AUCinf_Civ")
AMS <- unlist(p4$trial$analysis_models)
sc_meta <- rbindlist(lapply(PK, function(m) fread(file.path(oc_dir, sprintf("oc_scenarios_%s.csv", m)))[, .(pk_model = model, scenario = code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio)]))
sc_meta <- sc_meta[abs(target - 0.80) < 1e-9 | abs(target - 1.25) < 1e-9]
read_trials <- function(m) {
  f <- function(nm) { p_ <- file.path(om_dir, nm); if (file.exists(p_)) fread(p_) else NULL }
  main <- rbind(f(sprintf("oc_models_be_%s.csv.gz", m)), f(sprintf("oc_models_ext_be_%s.csv.gz", m)))[model %in% AMS & endpoint %in% c("Cmax", unname(G2V))]
  crit <- rbind(f(sprintf("oc_models_crit_be_%s.csv.gz", m)), f(sprintf("oc_models_ext_crit_be_%s.csv.gz", m)))
  x <- rbind(main, crit)[, pk_model := m]
  cnt <- x[, .(n_ep = uniqueN(endpoint)), by = .(trial, scenario, model)]
  if (any(cnt$n_ep != 1 + length(G2V))) stop(m, ": 시험·시나리오마다 Cmax와 AUC0-inf 9개가 모두 있어야 합니다")
  x
}
TR <- rbindlist(lapply(PK, read_trials))
W <- dcast(TR[, .(pk_model, trial, scenario, model, endpoint, ok = pass %in% TRUE)], pk_model + trial + scenario + model ~ endpoint, value.var = "ok")
for (g in names(G2V)) W[, (paste0("G2_", g)) := Cmax & get(G2V[[g]])]
gcols <- paste0("G2_", names(G2V))
L <- melt(W, id.vars = c("pk_model", "trial", "scenario", "model"), measure.vars = gcols, variable.name = "config", value.name = "ok")
rate <- function(x) { w <- wilson_ci(sum(x), length(x)); list(n_trials = length(x), pass_pct = w$est, lo = w$lo, hi = w$hi) }
T1 <- L[scenario %in% sc_meta$scenario, rate(ok), by = .(pk_model, scenario, analysis_model = model, config)]
T1 <- sc_meta[T1, on = c("pk_model", "scenario"), nomatch = NULL][, class := classify_type1(lo, hi, 5)]
T1[, config := as.character(config)]
T1[, `:=`(rule = sub("^G2_([ABC]).*$", "\\1", config), set = fifelse(config == "G2_B", "any", sub("^G2_[AC]_", "", config)))]
fwrite(T1, file.path(out_dir, "criteria_g2_type1.csv"))
# v1.0 저장본(oc/g2_rules_flags.csv)과 대조: M0의 규칙 A·C 세트 (i)·(ii)와 규칙 B는 같은 시험이므로, 시험 수가 같은 칸은 통과율이 같아야 한다
V10 <- fread(file.path(oc_dir, "g2_rules_flags.csv"))[region == "boundary" & config %in% c("G2_Ai", "G2_Aii", "G2_B", "G2_Ci", "G2_Cii")]
V10[, config := c(G2_Ai = "G2_A_i", G2_Aii = "G2_A_ii", G2_B = "G2_B", G2_Ci = "G2_C_i", G2_Cii = "G2_C_ii")[config]]
VC <- merge(T1[analysis_model == "M0", .(pk_model, scenario, config, n_trials, pass_pct)], V10[, .(pk_model = model, scenario, config, n_trials_v10 = n_trials, pass_pct_v10 = pass_pct)], by = c("pk_model", "scenario", "config"))
VC[, `:=`(same_n = n_trials == n_trials_v10, identical = n_trials == n_trials_v10 & abs(pass_pct - pass_pct_v10) < 1e-9)]
fwrite(VC, file.path(out_dir, "criteria_check_vs_v10.csv"))
if (nrow(VC) != 80 || !any(VC$same_n) || !all(VC[same_n == TRUE, identical])) stop("v1.0 저장 G2와 불일치: ", paste(VC[same_n == TRUE & identical == FALSE, sprintf("%s %s %s", pk_model, scenario, config)], collapse = ", "))
PWt <- L[scenario %in% c("S00", "F097"), rate(ok), by = .(pk_model, scenario, analysis_model = model, config)]
fwrite(PWt, file.path(out_dir, "criteria_g2_power.csv"))
BI <- TR[endpoint %in% G2V & is.finite(est) & scenario %in% sc_meta$scenario, .(n = .N, mean_log_gmr = mean(est), sd_log_gmr = sd(est)), by = .(pk_model, scenario, analysis_model = model, endpoint)]
BI <- sc_meta[BI, on = c("pk_model", "scenario"), nomatch = NULL]
BI[, `:=`(bias_pct = 100 * (exp(mean_log_gmr - log(auc_ratio)) - 1), bias_mc_se_pct = 100 * sd_log_gmr / sqrt(n))]
BI[, bias_dir := bias_direction(mean_log_gmr - log(auc_ratio), auc_ratio, mean_log_gmr - log(auc_ratio) - 1.96 * sd_log_gmr / sqrt(n), mean_log_gmr - log(auc_ratio) + 1.96 * sd_log_gmr / sqrt(n))]
BI[, config := names(G2V)[match(endpoint, G2V)]]
fwrite(BI, file.path(out_dir, "criteria_bias.csv"))
# 판정 불안정성: (a) 9개 변형 전체, (b) 세트 고정 규칙 A·B·C, (c) 규칙 고정 세트 i–iv
flip <- function(cols) { m <- as.matrix(W[, ..cols]); rowSums(m) > 0 & rowSums(m) < ncol(m) }
W[, inst_all := flip(gcols)]
for (s_ in SETS) W[, (paste0("inst_rules_", s_)) := flip(c(paste0("G2_A_", s_), "G2_B", paste0("G2_C_", s_)))]
W[, `:=`(inst_sets_A = flip(paste0("G2_A_", SETS)), inst_sets_C = flip(paste0("G2_C_", SETS)))]
icols <- grep("^inst_", names(W), value = TRUE)
IS <- W[, c(list(n_trials = .N), lapply(.SD, function(z) 100 * mean(z))), by = .(pk_model, scenario, analysis_model = model), .SDcols = icols]
fwrite(IS, file.path(out_dir, "criteria_instability.csv"))

# ---- 3) 그림 ----------------------------------------------------------------------------------------------------------------
fd <- CI[, .(population = fifelse(distribution == "study_60_90", "60-90 kg", sprintf("atopic: %s", distribution)), variant, set, lz = lz_fail_pct, est = est_fail_pct)]
fd <- melt(fd, id.vars = c("population", "variant", "set"), variable.name = "reason", value.name = "pct")
fd[, reason := factor(fifelse(reason == "lz", "lambda-z not estimable", "estimable, criteria not met"), c("estimable, criteria not met", "lambda-z not estimable"))]
fd[, vl := factor(MODEL_EN[variant], unique(MODEL_EN[c("base", "struct2020", "k2016_bmi_vc0817")]))]
g1 <- ggplot(fd[population %in% c("60-90 kg", "atopic: primary")], aes(set, pct, fill = reason)) + geom_col(width = 0.7) +
  facet_grid(population ~ vl) + scale_fill_manual(values = c("estimable, criteria not met" = "#2a78d6", "lambda-z not estimable" = "#eb6834"), name = NULL) +
  labs(x = "Criteria set: (i) adj. R-squared 0.80; (ii) (i) + span 2; (iii) adj. R-squared 0.90; (iv) (iii) + span 3 (all with extrapolation at most 20%)", y = "Subjects without a reliable AUC0-inf (%)",
       title = "Share of subjects without a reliable AUC0-inf by criteria set", subtitle = "Planned schedule, 300 mg, 20,000 subjects per population and model") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_criteria_failure.png"), g1, width = 10, height = 6, dpi = 150, bg = "white")
g2 <- ggplot(T1[analysis_model == "M0"], aes(factor(config, gcols, sub("G2_", "", gcols)), pass_pct)) + geom_hline(yintercept = 5, linetype = "dashed", colour = "#52514e", linewidth = 0.3) +
  geom_point(aes(colour = pk_model), position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.6) +
  scale_colour_manual(values = c(k2016 = "#2a78d6", k2020 = "#eb6834"), labels = MODEL_EN[c("k2016", "k2020")], name = NULL) +
  labs(x = "AUC0-inf rule and criteria set (G2 = AUC0-inf + Cmax)", y = "Boundary type I error (%)", title = "Boundary type I error of AUC0-inf + Cmax by handling rule and criteria set",
       subtitle = "16 boundary scenarios, pooled t-test (M0), dashed line 5%") + theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_criteria_g2.png"), g2, width = 9, height = 5, dpi = 150, bg = "white")

# ---- 4) 결론 문안 -------------------------------------------------------------------------------------------------------------
c6 <- CI[distribution == "study_60_90"]; ca <- CI[distribution == "primary"]
fs <- function(d, s_) rg(d[set == s_, fail_pct])
g2line <- function(am) paste(vapply(names(G2V), function(g) { x <- T1[analysis_model == am & config == paste0("G2_", g)]
  sprintf("%s %d/16 above 5%% (max %s%%)", sub("_", " ", g), sum(x$pass_pct > 5), f2(max(x$pass_pct))) }, ""), collapse = "; ")
ins <- IS[scenario %in% sc_meta$scenario & analysis_model == "M0"]
en <- c("# Sensitivity to the conventional lambda-z reliability criteria", "",
  "Phoenix WinNonlin does not define reliability criteria for lambda-z: the Lambda Z Acceptance Criteria on the Rules tab are optional user entries (minimum adjusted R-squared, maximum extrapolated percentage, span), and profiles that do not meet them are only flagged. The value 0.80 is a statistical analysis plan convention; 0.90, sometimes with a span of at least 3 half-lives, is common in public statistical analysis plans.", "",
  "## Subjects without a reliable AUC0-inf (planned schedule, 300 mg, 20,000 subjects per model)", "",
  sprintf("- Study population (60 to 90 kg, 2016 model and Model 1): set (i) %s; set (ii) %s; set (iii) %s; set (iv) %s. Per arm of 117 subjects: %s under set (i) and %s under set (iii).",
          fs(c6, "i"), fs(c6, "ii"), fs(c6, "iii"), fs(c6, "iv"), rg(c6[set == "i", fail_per_arm], f1, ""), rg(c6[set == "iii", fail_per_arm], f1, "")),
  sprintf("- The failure share under set (i) is therefore the lower end of the conventional range; under the commonly used 0.90 it is %s in the study population.", fs(c6, "iii")),
  sprintf("- Lambda-z itself is not estimable in %s of subjects (study population); the rest of the failures have an estimable lambda-z that does not meet the criteria.", rg(c6[set == "i", lz_fail_pct], f2)),
  sprintf("- Failing subjects differ from retained subjects (study population, set (i)): body weight %s kg higher, geometric mean AUClast ratio %s. In failing subjects the observed-to-true ratio of AUClast (observed AUClast / true AUC0-inf) has a median of %s (5th percentile %s).",
          rg(c6[set == "i", wt_diff], f2, ""), rg(c6[set == "i", auclast_gm_ratio_fail_to_retained], f3, ""), rg(100 * c6[set == "i", auclast_over_true_fail_median], f1), rg(100 * c6[set == "i", auclast_over_true_fail_p05], f1)),
  sprintf("- Robustness check only (report Appendix I; not evidence for the endpoint choice): with an adult atopic dermatitis body-weight distribution (placeholder from published summaries) and three model variants, set (i) %s; set (iii) %s; set (iv) %s.", fs(ca, "i"), fs(ca, "iii"), fs(ca, "iv")), "",
  "## Boundary type I error of AUC0-inf + Cmax by rule and criteria set (16 boundary scenarios)", "",
  sprintf("- M0: %s.", g2line("M0")), sprintf("- M1: %s.", g2line("M1")), "",
  "## Decision instability (same trials, only the rule or criteria set changed; 16 boundary scenarios, M0)", "",
  sprintf("- Across all 9 variants the G2 decision changes in %s of trials per cell; across rules A, B, C with set (i) in %s, with set (iii) in %s; across sets within rule A in %s, within rule C in %s.",
          mm(ins$inst_all), mm(ins$inst_rules_i), mm(ins$inst_rules_iii), mm(ins$inst_sets_A), mm(ins$inst_sets_C)), "",
  if (nrow(RC)) sprintf("Restart check: the section1 rows of the trials completed before the restart (%s) are identical after the restart (%s rows).", paste(sprintf("%s, %s", MODEL_EN[RC$pk_model], format(RC$trials_compared, big.mark = ",")), collapse = "; "), format(sum(RC$rows_matched), big.mark = ",")) else "")
ko <- c("# λz 신뢰 기준 관행값 민감도", "",
  "Phoenix WinNonlin은 λz 신뢰 기준을 정하지 않는다: Rules 탭의 Lambda Z Acceptance Criteria는 사용자가 넣는 선택 항목(최소 adjusted R², 최대 외삽 %, span)이며 미달 프로필은 표시만 된다. 0.80은 SAP 관행값이고, 공개 SAP에서는 0.90(span 3 이상 병기 포함)도 흔하다.", "",
  "## 신뢰할 수 있는 AUC0-inf를 얻지 못하는 대상자 (B0, 300 mg, 모델당 20,000명)", "",
  sprintf("- 연구 모집단(60–90 kg, 2016·Model 1): 세트 (i) %s; (ii) %s; (iii) %s; (iv) %s. arm당 117명 기준 (i) %s명, (iii) %s명.", fs(c6, "i"), fs(c6, "ii"), fs(c6, "iii"), fs(c6, "iv"),
          rg(c6[set == "i", fail_per_arm], f1, ""), rg(c6[set == "iii", fail_per_arm], f1, "")),
  sprintf("- 따라서 0.80 기준의 탈락률은 관행 범위의 하한이며, 통상 쓰이는 0.90에서는 연구 모집단 %s로 커진다.", fs(c6, "iii")),
  sprintf("- 미달자는 잔류자보다 체중 %s kg 높고 AUClast 기하평균비 %s(연구 모집단, 세트 (i)). 미달자에서 관측 대 참 비(관측 AUClast / 참 AUC0-inf) 중앙값 %s(5백분위 %s).",
          rg(c6[set == "i", wt_diff], f2, "", " ~ "), rg(c6[set == "i", auclast_gm_ratio_fail_to_retained], f3, "", " ~ "), rg(100 * c6[set == "i", auclast_over_true_fail_median], f1, "%", " ~ "), rg(100 * c6[set == "i", auclast_over_true_fail_p05], f1, "%", " ~ ")),
  sprintf("- 견고성 확인만(보고서 부록 I, 평가변수 선택의 근거 아님): 아토피 성인 체중 분포(공개 요약 기반 자리표시자), 모델 변형 3개에서 (i) %s; (iii) %s; (iv) %s.", fs(ca, "i"), fs(ca, "iii"), fs(ca, "iv")), "",
  "## 규칙·세트별 G2 경계 1종 오류(16칸)", "", sprintf("- M0: %s.", g2line("M0")), sprintf("- M1: %s.", g2line("M1")), "",
  sprintf("## 판정 불안정성(M0, 칸별 중앙값/최대): 9개 변형 전체 %s, 세트 (i)에서 규칙만 %s, (iii)에서 %s, 규칙 A에서 세트만 %s, 규칙 C에서 %s.", mmk(ins$inst_all), mmk(ins$inst_rules_i), mmk(ins$inst_rules_iii),
          mmk(ins$inst_sets_A), mmk(ins$inst_sets_C)))
writeLines(en, file.path(out_dir, "criteria_conclusion_en.md")); writeLines(ko, file.path(out_dir, "criteria_conclusion_ko.md"))
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", paste(en, collapse = "\n"))) stop("영문 결론 규칙 위반")
cat(en, sep = "\n")
