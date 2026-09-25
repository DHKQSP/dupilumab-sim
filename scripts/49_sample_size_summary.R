#!/usr/bin/env Rscript
# 49_sample_size_summary.R — 표본 수 표 요약 (추가 지시 2026-09-26 §3, config/prereg_20260926.yaml section3, D-055). 새 모의 없음.
# 입력: results/sample_size/(scripts/47 산출), results/oc_models/oc_models_be_<model>.csv.gz(S00·F097, n = 117), results/oc_models/truth_F097.csv,
#       results/oc/oc_scenarios_<model>.csv(F_down_090·095 참값 비).
# 산출(results/sample_size/):
#   ss_table_power.csv       1차 규칙(Cmax CV 모의값 고정, Cmax GMR = AUClast GMR): 입력 × CV × GMR × n × 분석 모형의 분석식·모의 검정력
#   ss_table_n_needed.csv    CV × GMR × 분석 모형 × 목표(90%, 85%): 필요 평가 가능 n/군과 무작위배정 n/군, 민감도 규칙·입력 병기
#   ss_pk_check.csv          PK 모델 시험(n = 117) 경험 검정력 대 같은 입력의 분석식
#   fig_ss_power.png, ss_conclusion_en.md, ss_conclusion_ko.md
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
pr <- read_cfg("prereg_20260926.yaml"); p3 <- pr$section3
design <- read_cfg("trial_design.yaml")
out_dir <- proj_path(p3$out_dir); CL <- design$be$ci_level; LIMS <- as.numeric(unlist(design$be$limits))
CVS <- as.numeric(unlist(p3$grid$auclast_cv_pct)); GMRS <- as.numeric(unlist(p3$grid$true_gmr)); NS <- as.integer(unlist(p3$grid$n_evaluable_per_arm))
BASE_CV <- as.numeric(p3$base_cv_pct); SENS_CV <- as.numeric(p3$sensitivity_cv_pct)
EVAL_FRAC <- as.numeric(design$n_per_arm) / as.numeric(design$n_randomized_per_arm)
f1 <- function(x) formatC(round(x, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2); s2 <- function(x) sprintf("%+.2f", round(x, 2) + 0)
MODEL_EN <- c(k2016 = "2016 model", k2020 = "Model 1"); MODEL_KO <- c(k2016 = "2016 모델", k2020 = "Model 1")

INP <- fread(file.path(out_dir, "ss_inputs.csv")); WM <- fread(file.path(out_dir, "ss_weight_moments.csv"))
AN <- fread(file.path(out_dir, "ss_power_analytic.csv")); NN <- fread(file.path(out_dir, "ss_n_needed.csv")); MC <- fread(file.path(out_dir, "ss_power_mc.csv"))
KM <- fread(file.path(out_dir, "ss_cmax_gmr_mechanistic.csv"))[, .(k_mech = mean(k)), by = input_model]
WM2 <- ss_weight_moments(75, 9, 60, 90, as.numeric(design$stratification$split_kg$value))
wt <- weight_spec_from_design(design, "base"); stopifnot(wt$mean == 75, wt$sd == 9, identical(as.numeric(wt$trunc), c(60, 90)))
stopifnot(isTRUE(all.equal(unlist(WM2), unlist(WM), tolerance = 1e-12)))
pw <- function(n, cv, ga, gc, m, am, rule) 100 * ss_p2_power(n, cv, ga, gc, INP[input_model == m], am, rule, WM, CL, LIMS)
cg <- function(g, rule, m) if (rule == "equal") g else exp(KM[input_model == m, k_mech] * log(g))
# 검사: R/samplesize.R의 분석식이 scripts/47 저장 격자와 같은가
chk <- AN[, .(power_pct, re = mapply(function(m, cr, gr, cv, g, n, am) pw(n, cv, g, cg(g, gr, m), m, am, cr), input_model, cmax_cv_rule, cmax_gmr_rule, cv, gmr, n, analysis_model))]
if (max(abs(chk$power_pct - chk$re)) > 1e-8) stop("R/samplesize.R 분석식이 scripts/47 격자와 다릅니다")

# ---- 표 ---------------------------------------------------------------------------------------------------------------------
TP <- merge(AN[cmax_cv_rule == "fixed" & cmax_gmr_rule == "equal", .(input_model, cv, gmr, n, analysis_model, analytic_pct = power_pct)],
            MC[, .(input_model, cv, gmr, n, analysis_model, mc_trials = n_trials, mc_pct = power_pct, mc_lo = lo, mc_hi = hi, diff_pp, analytic_in_ci)],
            by = c("input_model", "cv", "gmr", "n", "analysis_model"), all.x = TRUE)
setorder(TP, input_model, gmr, cv, n, analysis_model)
fwrite(TP, file.path(out_dir, "ss_table_power.csv"))
TN <- dcast(NN, cv + gmr + analysis_model + target_pct ~ input_model + cmax_cv_rule + cmax_gmr_rule, value.var = "n_evaluable_per_arm")
prim <- NN[input_model == "k2016" & cmax_cv_rule == "fixed" & cmax_gmr_rule == "equal", .(cv, gmr, analysis_model, target_pct, n_evaluable_per_arm, n_randomized_per_arm)]
TN <- merge(prim, TN, by = c("cv", "gmr", "analysis_model", "target_pct"))
setorder(TN, gmr, target_pct, cv, analysis_model)
fwrite(TN, file.path(out_dir, "ss_table_n_needed.csv"))

# ---- PK 모델 점검(n = 117, 모델 고유 CV) ---------------------------------------------------------------------------------------
PKC <- list()
for (m in c("k2016", "k2020")) {
  sc <- fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", m)))
  tr <- rbind(data.table(scenario = "S00", auc_ratio = 1, cmax_ratio = 1), fread(proj_path("results", "oc_models", "truth_F097.csv"))[pk_model == m, .(scenario, auc_ratio, cmax_ratio)],
              sc[code %in% c("F_down_090", "F_down_095"), .(scenario = code, auc_ratio, cmax_ratio)])
  b1 <- fread(proj_path("results", "oc_models", sprintf("oc_models_be_%s.csv.gz", m)))[scenario %in% c("S00", "F097") & endpoint %in% c("AUClast", "Cmax")]
  b2 <- fread(file.path(out_dir, sprintf("ss_pk_trials_%s.csv.gz", m)))
  b <- rbind(b1[, .(trial, scenario, endpoint, model, pass)], b2[, .(trial, scenario, endpoint, model, pass)])
  w <- dcast(b, trial + scenario + model ~ endpoint, value.var = "pass")
  e <- w[, { ci <- wilson_ci(sum(AUClast %in% TRUE & Cmax %in% TRUE), .N); .(n_trials = .N, empirical_pct = ci$est, lo = ci$lo, hi = ci$hi) }, by = .(scenario, analysis_model = model)]
  e <- merge(e, tr, by = "scenario")
  cvn <- INP[input_model == m, cv_auc_pct]
  e[, analytic_pct := mapply(function(ga, gc, am) if (am == "M2") NA_real_ else pw(117L, cvn, ga, gc, m, am, "fixed"), auc_ratio, cmax_ratio, analysis_model)]
  e[, `:=`(pk_model = m, natural_cv_auc_pct = cvn, diff_pp = empirical_pct - analytic_pct, analytic_in_ci = analytic_pct >= lo & analytic_pct <= hi)]
  PKC[[m]] <- e
}
PKC <- rbindlist(PKC); setcolorder(PKC, c("pk_model", "scenario", "auc_ratio", "cmax_ratio", "analysis_model")); setorder(PKC, pk_model, scenario, analysis_model)
fwrite(PKC, file.path(out_dir, "ss_pk_check.csv"))

# ---- 그림: 1차 입력, GMR 0.95, 분석식 곡선 + 모의 점 ---------------------------------------------------------------------------
cur <- CJ(cv = CVS, n = 60:200, analysis_model = c("M0", "M1"))
cur[, power_pct := mapply(function(cv, n, am) pw(n, cv, 0.95, 0.95, "k2016", am, "fixed"), cv, n, analysis_model)]
cur[, cvl := factor(sprintf("CV %g%%", cv), sprintf("CV %g%%", CVS))]
pts <- TP[input_model == "k2016" & gmr == 0.95][, cvl := factor(sprintf("CV %g%%", cv), sprintf("CV %g%%", CVS))]
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"); names(PAL) <- levels(cur$cvl)
g <- ggplot(cur, aes(n, power_pct, colour = cvl)) + geom_hline(yintercept = c(85, 90), linetype = "dashed", colour = "#52514e", linewidth = 0.3) + geom_vline(xintercept = 117, colour = "#8a8984", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_pointrange(data = pts, aes(y = mc_pct, ymin = mc_lo, ymax = mc_hi), size = 0.2) +
  facet_wrap(~analysis_model, labeller = labeller(analysis_model = c(M0 = "M0 (pooled t)", M1 = "M1 (ANOVA with weight stratum)"))) +
  scale_colour_manual(values = PAL, name = "AUClast CV") + coord_cartesian(ylim = c(50, 100)) +
  labs(x = "Evaluable subjects per arm", y = "P2 power, % (lines analytic; points simulation, 5,000 trials, Wilson 95%)",
       title = "Power of P2 (AUClast and Cmax) at true GMR 0.95 by AUClast CV", subtitle = "Inputs from the 2016 model (Cmax CV and correlation from simulation); dashed lines 85% and 90%; grey line n = 117") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_ss_power.png"), g, width = 9.5, height = 5.5, dpi = 150, bg = "white")

# ---- 결론 문안 --------------------------------------------------------------------------------------------------------------
nn <- function(cv_, g, am, tg, col = "n_evaluable_per_arm", m = "k2016", cr = "fixed", gr = "equal") {
  x <- NN[input_model == m & cmax_cv_rule == cr & cmax_gmr_rule == gr & cv == cv_ & gmr == g & analysis_model == am & target_pct == tg][[col]]
  stopifnot(length(x) == 1); x }
nline <- function(g, am, tg) paste(vapply(CVS, function(cv) sprintf("%g%%: %d (%d)", cv, nn(cv, g, am, tg), nn(cv, g, am, tg, "n_randomized_per_arm")), ""), collapse = "; ")
an117 <- function(cv_, g, am) { x <- TP[input_model == "k2016" & cv == cv_ & gmr == g & n == 117 & analysis_model == am, analytic_pct]; stopifnot(length(x) == 1); x }
mcs <- MC[, .(max_abs = max(abs(diff_pp)), in_ci = sum(analytic_in_ci), n = .N)]
pk_line <- function(lang) paste(vapply(seq_len(nrow(PKC[analysis_model != "M2"])), function(i) { r <- PKC[analysis_model != "M2"][i]
  sprintf("%s %s %s %s%% [%s, %s] vs %s%%", MODEL_EN[[r$pk_model]], r$scenario, r$analysis_model, f2(r$empirical_pct), f2(r$lo), f2(r$hi), f2(r$analytic_pct)) }, ""), collapse = "; ")
save_g <- function(g, tg, cv) nn(cv, g, "M0", tg) - nn(cv, g, "M1", tg)
inp <- INP[input_model == "k2016"]
en <- c("# Sample size for P2 (AUClast and Cmax) joint power", "",
  sprintf("Generated by `scripts/49_sample_size_summary.R` from `results/sample_size/` (pre-registered in `config/prereg_20260926.yaml`, section3). n is the number of evaluable subjects per arm; the randomized number per arm (in parentheses) assumes the protocol ratio %d/%d evaluable.", as.integer(design$n_per_arm), as.integer(design$n_randomized_per_arm)), "",
  sprintf("Inputs (2016 model, B0, 60 to 90 kg, 20,000 subjects): Cmax log-scale CV %s%%, correlation of log AUClast and log Cmax %s, slopes on log body weight %s (AUClast) and %s (Cmax); the weight stratum explains %s%% of the AUClast variance. Model 1 inputs (sensitivity): Cmax CV %s%%, correlation %s.",
          f1(inp$cv_cmax_pct), f2(inp$rho_total), f2(inp$beta_auc), f2(inp$beta_cmax), f1(100 * inp$r2_stratum_auc), f1(INP[input_model == "k2020", cv_cmax_pct]), f2(INP[input_model == "k2020", rho_total])), "",
  "CV sources: Li 2020 Table 3, 300 mg arms, SD/mean 35.4 to 51.4% (pooled 42.0%); Cohen 2022 about 52% (200 mg, different presentation, device comparison); protocol assumption 43%; PK models 40.3% (2016 model) and 42.8% (Model 1).", "",
  "## Evaluable n per arm for P2 power at true GMR 0.95 (randomized n per arm)", "",
  sprintf("- 90%% power, M0: %s.", nline(0.95, "M0", 90)), sprintf("- 90%% power, M1: %s.", nline(0.95, "M1", 90)),
  sprintf("- 85%% power, M0: %s.", nline(0.95, "M0", 85)), sprintf("- 85%% power, M1: %s.", nline(0.95, "M1", 85)), "",
  sprintf("At the protocol CV of %g%% with n = 117, P2 power at GMR 0.95 is %s%% (M0) and %s%% (M1); at the sensitivity CV of %g%% it is %s%% and %s%%. M1 needs %d (90%%) and %d (85%%) fewer evaluable subjects per arm than M0 at CV %g%%.",
          BASE_CV, f1(an117(BASE_CV, 0.95, "M0")), f1(an117(BASE_CV, 0.95, "M1")), SENS_CV, f1(an117(SENS_CV, 0.95, "M0")), f1(an117(SENS_CV, 0.95, "M1")), save_g(0.95, 90, BASE_CV), save_g(0.95, 85, BASE_CV), BASE_CV), "",
  sprintf("At GMR 0.90, 90%% power requires (M0 / M1): %s. At GMR 1.00: %s.",
          paste(vapply(CVS, function(cv) sprintf("CV %g%% %d / %d", cv, nn(cv, 0.90, "M0", 90), nn(cv, 0.90, "M1", 90)), ""), collapse = "; "),
          paste(vapply(CVS, function(cv) sprintf("CV %g%% %d / %d", cv, nn(cv, 1.00, "M0", 90), nn(cv, 1.00, "M1", 90)), ""), collapse = "; ")), "",
  sprintf("Sensitivity (90%% power, GMR 0.95, CV %g%%, M1): Model 1 inputs %d; Cmax CV proportional to AUClast CV %d; Cmax GMR mechanistic (log Cmax ratio = %s x log AUC ratio) %d.",
          BASE_CV, nn(BASE_CV, 0.95, "M1", 90, m = "k2020"), nn(BASE_CV, 0.95, "M1", 90, cr = "ratio"), f2(KM[input_model == "k2016", k_mech]), nn(BASE_CV, 0.95, "M1", 90, gr = "mechanistic")), "",
  "## Checks", "",
  sprintf("- Statistical simulation (5,000 trials per cell, %d cells): the analytic power lies within the Wilson 95%% interval in %d cells; maximum absolute difference %s points.", mcs$n, mcs$in_ci, f2(mcs$max_abs)),
  sprintf("- PK-model trials at n = 117 and the models' own CV (empirical [Wilson 95%%] vs analytic): %s.", pk_line("en")))
ko <- c("# P2(AUClast, Cmax) 동시 검정력 표본 수", "",
  sprintf("`scripts/49_sample_size_summary.R`가 `results/sample_size/`에서 생성(사전 등록 `config/prereg_20260926.yaml` section3). n = 군당 평가 가능 대상자, 괄호 = 군당 무작위배정(프로토콜 비 %d/%d).", as.integer(design$n_per_arm), as.integer(design$n_randomized_per_arm)), "",
  sprintf("입력(2016 모델, B0, 60–90 kg, 20,000명): Cmax log-CV %s%%, log AUClast–log Cmax 상관 %s, log 체중 기울기 %s(AUClast)·%s(Cmax), 체중 층이 AUClast 분산의 %s%% 설명. Model 1(민감도): Cmax CV %s%%, 상관 %s.",
          f1(inp$cv_cmax_pct), f2(inp$rho_total), f2(inp$beta_auc), f2(inp$beta_cmax), f1(100 * inp$r2_stratum_auc), f1(INP[input_model == "k2020", cv_cmax_pct]), f2(INP[input_model == "k2020", rho_total])), "",
  "CV 출처: Li 2020 Table 3 300 mg 6개 군 SD/평균 35.4–51.4%(합동 42.0%), Cohen 2022 약 52%(200 mg, 다른 제형, 기기 비교), 프로토콜 43%, PK 모델 40.3%(2016)·42.8%(Model 1).", "",
  "## 참 GMR 0.95에서 P2 검정력에 필요한 군당 평가 가능 n (무작위배정 n)", "",
  sprintf("- 90%%, M0: %s.", nline(0.95, "M0", 90)), sprintf("- 90%%, M1: %s.", nline(0.95, "M1", 90)),
  sprintf("- 85%%, M0: %s.", nline(0.95, "M0", 85)), sprintf("- 85%%, M1: %s.", nline(0.95, "M1", 85)), "",
  sprintf("프로토콜 CV %g%%·n = 117에서 GMR 0.95의 P2 검정력은 %s%%(M0), %s%%(M1), 민감도 CV %g%%에서 %s%%, %s%%. CV %g%%에서 M1은 M0보다 군당 %d명(90%%), %d명(85%%) 적다.",
          BASE_CV, f1(an117(BASE_CV, 0.95, "M0")), f1(an117(BASE_CV, 0.95, "M1")), SENS_CV, f1(an117(SENS_CV, 0.95, "M0")), f1(an117(SENS_CV, 0.95, "M1")), BASE_CV, save_g(0.95, 90, BASE_CV), save_g(0.95, 85, BASE_CV)), "",
  sprintf("GMR 0.90에서 90%% 검정력(M0 / M1): %s. GMR 1.00: %s.",
          paste(vapply(CVS, function(cv) sprintf("CV %g%% %d / %d", cv, nn(cv, 0.90, "M0", 90), nn(cv, 0.90, "M1", 90)), ""), collapse = "; "),
          paste(vapply(CVS, function(cv) sprintf("CV %g%% %d / %d", cv, nn(cv, 1.00, "M0", 90), nn(cv, 1.00, "M1", 90)), ""), collapse = "; ")), "",
  sprintf("민감도(90%%, GMR 0.95, CV %g%%, M1): Model 1 입력 %d; Cmax CV를 AUClast CV에 비례 %d; Cmax GMR 기전 비례(log Cmax 비 = %s × log AUC 비) %d.",
          BASE_CV, nn(BASE_CV, 0.95, "M1", 90, m = "k2020"), nn(BASE_CV, 0.95, "M1", 90, cr = "ratio"), f2(KM[input_model == "k2016", k_mech]), nn(BASE_CV, 0.95, "M1", 90, gr = "mechanistic")), "",
  "## 점검", "",
  sprintf("- 통계 모의(칸당 5,000회, %d칸): 분석식이 Wilson 95%% 구간 안인 칸 %d, 최대 절대 차이 %s%%p.", mcs$n, mcs$in_ci, f2(mcs$max_abs)),
  sprintf("- PK 모델 시험(n = 117, 모델 고유 CV; 경험 [Wilson 95%%] 대 분석식): %s.", pk_line("ko")))
writeLines(en, file.path(out_dir, "ss_conclusion_en.md")); writeLines(ko, file.path(out_dir, "ss_conclusion_ko.md"))
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", paste(en, collapse = "\n"))) stop("영문 결론에 한글·en/em dash·U+2212가 있습니다")
cat(en, sep = "\n")
