#!/usr/bin/env Rscript
# 48_lloq_summary.R — 연구 LLOQ 민감도 요약 (추가 지시 2026-09-26 §2, config/prereg_20260926.yaml section2, D-054). 새 모의 없음.
# 입력: results/lloq/(scripts/46 산출). 산출(results/lloq/):
#   lloq_individual_table.csv    모델 × 잔차 변형 × LLOQ: 사전 등록 지표 + 0.078 대비 쌍대 변화(신뢰 (i)·(ii), λz 실패, 커버리지 < 80%, NCA 외삽 > 20%)
#   lloq_tlast_distribution.csv  tlast의 계획 채혈일 분포(%)
#   lloq_cliff_table.csv         절벽 끝(LLOQ 도달) 요약 + 현행 일정의 절벽 채혈점 1·2·≥3점(명목·방문 허용창, 1·2일 정의)
#   lloq_trial_type1.csv         주 모델 경계 3칸 × LLOQ × 잔차 변형 × 분석 모형 × 구성: 통과율, Wilson, 분류, 0.078 대비 쌍대 변화
#   lloq_trial_power.csv         S00 검정력(같은 구조)
#   lloq_trial_drop.csv          군별 λz·신뢰 충족 평균 수와 tlast 중앙값의 시험 평균
#   fig_lloq_individual.png, fig_lloq_trial.png, lloq_conclusion_en.md, lloq_conclusion_ko.md
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
pr <- read_cfg("prereg_20260926.yaml")$section2
design <- read_cfg("trial_design.yaml")
out_dir <- Sys.getenv("DUPI_LLOQ_OUT", proj_path(pr$out_dir)); L0 <- as.numeric(pr$reference_lloq_mg_L)   # 환경변수는 시험용
LL <- sort(as.numeric(unlist(pr$lloq_grid_mg_L))); RES <- unlist(pr$residual_variants)
MODEL_EN <- c(k2016 = "2016 model", k2020 = "Model 1"); MODEL_KO <- c(k2016 = "2016 모델", k2020 = "Model 1")
f1 <- function(x) formatC(round(x, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2)
s1 <- function(x) sprintf("%+.1f", round(x, 1) + 0); s2 <- function(x) sprintf("%+.2f", round(x, 2) + 0)
fl <- function(x) trimws(formatC(x, format = "fg", digits = 3))
is0 <- function(x) abs(x - L0) < 1e-12

# ---- 1) 개인 수준 -----------------------------------------------------------------------------------------------------------
CK <- rbindlist(lapply(c("k2016", "k2020"), function(m) fread(file.path(out_dir, sprintf("lloq_individual_check_%s.csv", m)))))
stopifnot(all(CK$pass))
IND <- rbindlist(lapply(c("k2016", "k2020"), function(m) fread(file.path(out_dir, sprintf("lloq_individual_%s.csv", m)))))
NCA <- rbindlist(lapply(c("k2016", "k2020"), function(m) readRDS(file.path(out_dir, sprintf("lloq_individual_nca_%s.rds", m)))))
ind_cols <- c("extrap_true_median", "extrap_true_p95", "coverage_lt80_pct", "extrap_median", "extrap_p95", "extrap_gt20_pct", "reliable_no_span_pct", "reliable_pct",
              "lambda_fail_pct", "tlast_median", "tlast_p05", "tlast_p95", "quant_at_last_pct", "AUClast_logcv", "err_inf_sd")
stopifnot(all(ind_cols %in% names(IND)))
ind <- IND[, c("model", "resid", "lloq", "n", ind_cols), with = FALSE]
NCA[, `:=`(ind_rel_i = rel_i, ind_rel_ii = reliable %in% TRUE, ind_lz_fail = !(lambda_ok %in% TRUE), ind_cov80 = coverage_true < 0.80, ind_ext20 = pct_extrap > 20)]
ref <- NCA[is0(lloq), .(model, resid, id, r_i = ind_rel_i, r_ii = ind_rel_ii, r_lz = ind_lz_fail, r_cov = ind_cov80, r_ext = ind_ext20)]
pd <- merge(NCA, ref, by = c("model", "resid", "id"))
pc <- function(a, b) { d <- paired_prop_diff_ci(a %in% TRUE, b %in% TRUE); list(d$est, d$lo, d$hi) }
PD <- pd[, c(setNames(pc(ind_rel_i, r_i), c("d_rel_i_pp", "d_rel_i_lo", "d_rel_i_hi")), setNames(pc(ind_rel_ii, r_ii), c("d_rel_ii_pp", "d_rel_ii_lo", "d_rel_ii_hi")),
             setNames(pc(ind_lz_fail, r_lz), c("d_lz_fail_pp", "d_lz_fail_lo", "d_lz_fail_hi")), setNames(pc(ind_cov80, r_cov), c("d_cov80_pp", "d_cov80_lo", "d_cov80_hi")),
             setNames(pc(ind_ext20, r_ext), c("d_ext20_pp", "d_ext20_lo", "d_ext20_hi"))), by = .(model, resid, lloq)]
ind <- merge(ind, PD, by = c("model", "resid", "lloq"))
setorder(ind, model, resid, lloq)
fwrite(ind, file.path(out_dir, "lloq_individual_table.csv"))
TL <- NCA[, .N, by = .(model, resid, lloq, tlast_planned)][, pct := 100 * N / sum(N), by = .(model, resid, lloq)]
setorder(TL, model, resid, lloq, tlast_planned)
fwrite(TL, file.path(out_dir, "lloq_tlast_distribution.csv"))

# ---- 2) 절벽 -----------------------------------------------------------------------------------------------------------------
stopifnot(all(fread(file.path(out_dir, "lloq_cliff_check.csv"))$pass))
CS <- fread(file.path(out_dir, "lloq_cliff_summary.csv")); CP <- fread(file.path(out_dir, "lloq_cliff_points.csv"))
CT <- merge(CP[, .(model, lloq, timing, definition_day, n_subjects, pct_ge1, ge1_lo, ge1_hi, pct_ge2, ge2_lo, ge2_hi, pct_ge3, ge3_lo, ge3_hi)], CS, by = c("model", "lloq"))
setorder(CT, model, timing, definition_day, lloq)
fwrite(CT, file.path(out_dir, "lloq_cliff_table.csv"))

# ---- 3) 시험 수준 -------------------------------------------------------------------------------------------------------------
tp <- pr$trial
BE <- fread(file.path(out_dir, sprintf("lloq_trials_be_%s.csv.gz", tp$pk_model)))
ntr <- uniqueN(BE$trial); stopifnot(identical(sort(unique(BE$trial)), seq(as.integer(tp$trials[[1]]), as.integer(Sys.getenv("DUPI_LLOQ_TEST_N", tp$trials[[2]])))))   # 환경변수는 시험용
CFG <- lapply(read_cfg("prereg_20260926.yaml")$section1$configuration_endpoints[unlist(tp$configurations)], unlist)
W <- dcast(BE[, .(trial, lloq, resid, scenario, model, endpoint, ok = pass %in% TRUE)], trial + lloq + resid + scenario + model ~ endpoint, value.var = "ok")
for (cf in names(CFG)) W[, (cf) := Reduce(`&`, .SD), .SDcols = CFG[[cf]]]
LG <- melt(W, id.vars = c("trial", "lloq", "resid", "scenario", "model"), measure.vars = names(CFG), variable.name = "config", value.name = "ok")
LG[, config := as.character(config)]
R0 <- LG[is0(lloq), .(trial, resid, scenario, model, config, ok0 = ok)]
LG <- merge(LG, R0, by = c("trial", "resid", "scenario", "model", "config"))
TT <- LG[, { w <- wilson_ci(sum(ok), .N); d <- paired_prop_diff_ci(ok, ok0)
             .(n_trials = .N, n_pass = sum(ok), pass_pct = w$est, lo = w$lo, hi = w$hi, diff_vs_ref_pp = d$est, diff_lo = d$lo, diff_hi = d$hi) }, by = .(scenario, lloq, resid, model, config)]
sc <- fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", tp$pk_model)))[, .(scenario = code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio)]
T1 <- merge(TT[scenario != "S00"], sc, by = "scenario"); T1[, class := classify_type1(lo, hi, 5)]
setorder(T1, scenario, model, resid, config, lloq)
fwrite(T1, file.path(out_dir, "lloq_trial_type1.csv"))
PW <- TT[scenario == "S00"]; setorder(PW, model, resid, config, lloq)
fwrite(PW, file.path(out_dir, "lloq_trial_power.csv"))
DR <- fread(file.path(out_dir, sprintf("lloq_trials_drop_%s.csv.gz", tp$pk_model)))
DRs <- DR[, .(n_trials = .N, n_mean = mean(n), lambda_ok_mean = mean(n_lambda), reliable_ii_mean = mean(n_reliable), reliable_i_mean = mean(n_reliable_i),
              tlast_median_mean = mean(tlast_median)), by = .(lloq, resid, scenario)]
fwrite(DRs, file.path(out_dir, "lloq_trial_drop.csv"))

# ---- 4) 그림 --------------------------------------------------------------------------------------------------------------------
COLS <- c(k2016 = "#2a78d6", k2020 = "#eb6834")
fi <- melt(ind, id.vars = c("model", "resid", "lloq"), measure.vars = c("reliable_no_span_pct", "reliable_pct", "lambda_fail_pct", "coverage_lt80_pct", "extrap_true_median", "tlast_median"))
fi[, panel := factor(c(reliable_no_span_pct = "AUC0-inf reliable, criteria (i) (%)", reliable_pct = "AUC0-inf reliable, criteria (ii) (%)", lambda_fail_pct = "Lambda-z not estimable (%)",
                       coverage_lt80_pct = "True coverage below 80% (%)", extrap_true_median = "True extrapolated share, median (%)", tlast_median = "tlast, median (days after dose)")[as.character(variable)],
                     c("AUC0-inf reliable, criteria (i) (%)", "AUC0-inf reliable, criteria (ii) (%)", "Lambda-z not estimable (%)", "True coverage below 80% (%)", "True extrapolated share, median (%)", "tlast, median (days after dose)"))]
fi[, residual := factor(fifelse(resid == "fixed", "residual error as estimated (primary)", "additive error scaled with LLOQ"), c("residual error as estimated (primary)", "additive error scaled with LLOQ"))]
g1 <- ggplot(fi, aes(lloq, value, colour = model, linetype = residual)) + geom_vline(xintercept = L0, colour = "#8a8984", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) + facet_wrap(~panel, scales = "free_y", ncol = 3) + scale_x_log10(breaks = LL, labels = fl) +
  scale_colour_manual(values = COLS, labels = MODEL_EN, name = NULL) + labs(x = "LLOQ (mg/L, log scale); grey line = current assumption", y = NULL, linetype = NULL,
  title = "Individual-level metrics by LLOQ (B0, 60 to 90 kg, 20,000 subjects per model, same subjects and residual draws)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_lloq_individual.png"), g1, width = 10, height = 6.5, dpi = 150, bg = "white")
ft <- T1[config %in% c("P2", "G2_Aii", "G2_B") & resid == "fixed"]
ft[, cell := sprintf("%s %s (x%s)", mechanism, direction, fl(multiplier))]
ft[, cfgl := factor(c(P2 = "P2", G2_Aii = "G2-A(ii)", G2_B = "G2-B")[config], c("P2", "G2-A(ii)", "G2-B"))]
g2 <- ggplot(ft, aes(lloq, pass_pct, colour = cfgl)) + geom_hline(yintercept = 5, linetype = "dashed", colour = "#52514e", linewidth = 0.3) + geom_vline(xintercept = L0, colour = "#8a8984", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03, linewidth = 0.4) + geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_grid(model ~ cell, labeller = labeller(model = c(M0 = "M0 (pooled t)", M1 = "M1 (weight stratum)"))) + scale_x_log10(breaks = LL, labels = fl) +
  scale_colour_manual(values = c("P2" = "#2a78d6", "G2-A(ii)" = "#eb6834", "G2-B" = "#1baf7a"), name = NULL) +
  labs(x = "LLOQ (mg/L, log scale)", y = "Boundary type I error, % (Wilson 95%)", title = sprintf("Boundary type I error by LLOQ, 2016 model, %s trials per cell, residual error as estimated", format(ntr, big.mark = ","))) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "fig_lloq_trial.png"), g2, width = 10, height = 6, dpi = 150, bg = "white")

# ---- 5) 결론 문안 ---------------------------------------------------------------------------------------------------------------
gi <- function(m, rv, L, col) ind[model == m & resid == rv & abs(lloq - L) < 1e-12][[col]]
rng_l <- function(m, rv, col, f = f1, u = "%") sprintf("%s%s (LLOQ %s) to %s%s (LLOQ %s)", f(gi(m, rv, min(LL), col)), u, fl(min(LL)), f(gi(m, rv, max(LL), col)), u, fl(max(LL)))
ct <- function(m, L, tm = "nominal", dd = 1) CT[model == m & abs(lloq - L) < 1e-12 & timing == tm & definition_day == dd]
p2r <- T1[config == "P2"]; g2a <- T1[config == "G2_Aii"]; g2b <- T1[config == "G2_B"]
mx <- function(x) x[which.max(pass_pct)]
cellname <- function(r) sprintf("%s %s", r$mechanism, r$direction)
line_ind <- function(m, lang) {
  if (lang == "en") sprintf("- %s: reliability under criteria (i) %s (residual as estimated) [%s with the scaled residual]; criteria (ii) %s; lambda-z not estimable %s; true coverage below 80%% in %s of subjects; median true extrapolated share %s; median tlast %s.",
                            MODEL_EN[[m]], rng_l(m, "fixed", "reliable_no_span_pct"), rng_l(m, "scaled", "reliable_no_span_pct"), rng_l(m, "fixed", "reliable_pct"), rng_l(m, "fixed", "lambda_fail_pct", f2),
                            rng_l(m, "fixed", "coverage_lt80_pct", f2), rng_l(m, "fixed", "extrap_true_median", f2), rng_l(m, "fixed", "tlast_median", u = " days"))
  else sprintf("- %s: 신뢰 기준 (i) 충족 %s(잔차 그대로) [잔차 비례 %s]; 기준 (ii) %s; λz 산출 불가 %s; 실제 커버리지 80%% 미만 %s; 실제 외삽 중앙값 %s; tlast 중앙값 %s.",
               MODEL_KO[[m]], rng_l(m, "fixed", "reliable_no_span_pct"), rng_l(m, "scaled", "reliable_no_span_pct"), rng_l(m, "fixed", "reliable_pct"), rng_l(m, "fixed", "lambda_fail_pct", f2),
               rng_l(m, "fixed", "coverage_lt80_pct", f2), rng_l(m, "fixed", "extrap_true_median", f2), rng_l(m, "fixed", "tlast_median", u = "일"))
}
cliff_line <- function(m, lang) {
  a <- ct(m, min(LL)); b <- ct(m, L0); c_ <- ct(m, max(LL))
  if (lang == "en") sprintf("- %s: the true concentration reaches the LLOQ at a median of study Day %s (LLOQ %s mg/L), %s (%s) and %s (%s); median cliff length %s, %s and %s days; subjects with 1 or more samples in the cliff (1-day definition, nominal days): %s%%, %s%% and %s%%; 2 or more: %s%%, %s%%, %s%%.",
                            MODEL_EN[[m]], f1(a$lloq_studyday_median), fl(min(LL)), f1(b$lloq_studyday_median), fl(L0), f1(c_$lloq_studyday_median), fl(max(LL)), f2(a$len1_median), f2(b$len1_median), f2(c_$len1_median),
                            f1(a$pct_ge1), f1(b$pct_ge1), f1(c_$pct_ge1), f1(a$pct_ge2), f1(b$pct_ge2), f1(c_$pct_ge2))
  else sprintf("- %s: 참 농도가 LLOQ에 닿는 연구일 중앙값 Day %s(LLOQ %s mg/L), %s(%s), %s(%s); 절벽 길이 중앙값 %s, %s, %s일; 절벽 채혈점 1점 이상(1일 정의, 명목일) %s%%, %s%%, %s%%; 2점 이상 %s%%, %s%%, %s%%.",
               MODEL_KO[[m]], f1(a$lloq_studyday_median), fl(min(LL)), f1(b$lloq_studyday_median), fl(L0), f1(c_$lloq_studyday_median), fl(max(LL)), f2(a$len1_median), f2(b$len1_median), f2(c_$len1_median),
               f1(a$pct_ge1), f1(b$pct_ge1), f1(c_$pct_ge1), f1(a$pct_ge2), f1(b$pct_ge2), f1(c_$pct_ge2))
}
fam_line <- function(x, lab, lang) {
  x <- x[resid == "fixed"]
  by_m <- vapply(c("M0", "M1"), function(a) { y <- x[model == a]; k <- mx(y)
    if (lang == "en") sprintf("%s %s%% to %s%% (maximum %s, LLOQ %s, %s [%s, %s])", a, f2(min(y$pass_pct)), f2(max(y$pass_pct)), cellname(k), fl(k$lloq), f2(k$pass_pct), f2(k$lo), f2(k$hi))
    else sprintf("%s %s ~ %s%%(최대 %s, LLOQ %s, %s [%s, %s])", a, f2(min(y$pass_pct)), f2(max(y$pass_pct)), cellname(k), fl(k$lloq), f2(k$pass_pct), f2(k$lo), f2(k$hi)) }, "")
  nx <- vapply(c("M0", "M1"), function(a) sum(x[model == a, lo > 5]), 0L)
  if (lang == "en") sprintf("- %s: %s; cells (of 3 cells x 6 LLOQs) with Wilson lower bound above 5%%: M0 %d, M1 %d.", lab, paste(by_m, collapse = "; "), nx[1], nx[2])
  else sprintf("- %s: %s; Wilson 하한 > 5%% 칸(3칸 × LLOQ 6개 중): M0 %d, M1 %d.", lab, paste(by_m, collapse = "; "), nx[1], nx[2])
}
dmax <- T1[config == "P2" & resid == "fixed"][which.max(abs(diff_vs_ref_pp))]
pw <- PW[config == "P2" & resid == "fixed"]
en <- c("# Sensitivity to the LLOQ of the study assay", "",
  sprintf("Generated by `scripts/48_lloq_summary.R` from `results/lloq/` (pre-registered in `config/prereg_20260926.yaml`, section2; the study LLOQ is set in `config/assay.yaml`, currently %s mg/L). The same subjects, sampling times and residual draws are re-censored at each LLOQ (paired). Primary residual variant: as estimated (additive SD 0.03 mg/L); secondary: additive SD scaled by LLOQ/%s.", fl(L0), fl(L0)), "",
  "## Individual level (B0, 60 to 90 kg, 20,000 subjects per model)", "", line_ind("k2016", "en"), line_ind("k2020", "en"), "",
  "## Cliff capture (current schedule)", "", cliff_line("k2016", "en"), cliff_line("k2020", "en"), "",
  sprintf("## Boundary type I error, 2016 model (%s trials; Vmax both directions and F down; residual as estimated)", format(ntr, big.mark = ",")), "",
  fam_line(p2r, "P2", "en"), fam_line(g2a, "G2-A(ii)", "en"), fam_line(g2b, "G2-B", "en"),
  sprintf("- Largest paired change of P2 versus LLOQ %s: %s, %s, LLOQ %s: %s points [%s, %s].", fl(L0), dmax$model, cellname(dmax), fl(dmax$lloq), s2(dmax$diff_vs_ref_pp), s2(dmax$diff_lo), s2(dmax$diff_hi)),
  sprintf("- Power of P2 for identical products: M0 %s%% to %s%%, M1 %s%% to %s%% across the LLOQs.", f2(min(pw[model == "M0", pass_pct])), f2(max(pw[model == "M0", pass_pct])), f2(min(pw[model == "M1", pass_pct])), f2(max(pw[model == "M1", pass_pct]))), "",
  sprintf("Verification: at LLOQ %s with the residual as estimated, the individual-level results equal the stored results (subject level and summary), the cliff results equal the stored cliff files, and the trial-level M0 equals the stored operating-characteristic trials; the scaled residual at %s equals the unscaled one.", fl(L0), fl(L0)))
ko <- c("# 연구 분석법 LLOQ 민감도", "",
  sprintf("`scripts/48_lloq_summary.R`가 `results/lloq/`에서 생성(사전 등록 `config/prereg_20260926.yaml` section2; 연구 LLOQ는 `config/assay.yaml`, 현재 %s mg/L). 같은 대상자·채혈 시각·잔차를 LLOQ마다 다시 검열(쌍대). 1차 잔차: 추정값 그대로(가산 SD 0.03 mg/L), 2차: 가산 SD × LLOQ/%s.", fl(L0), fl(L0)), "",
  "## 개인 수준 (B0, 60–90 kg, 모델당 20,000명)", "", line_ind("k2016", "ko"), line_ind("k2020", "ko"), "",
  "## 절벽 채혈 (현행 일정)", "", cliff_line("k2016", "ko"), cliff_line("k2020", "ko"), "",
  sprintf("## 경계 1종 오류, 2016 모델 (%s회; Vmax 양방향, F 하향; 잔차 추정값 그대로)", format(ntr, big.mark = ",")), "",
  fam_line(p2r, "P2", "ko"), fam_line(g2a, "G2-A(ii)", "ko"), fam_line(g2b, "G2-B", "ko"),
  sprintf("- LLOQ %s 대비 P2의 최대 쌍대 변화: %s, %s, LLOQ %s: %s%%p [%s, %s].", fl(L0), dmax$model, cellname(dmax), fl(dmax$lloq), s2(dmax$diff_vs_ref_pp), s2(dmax$diff_lo), s2(dmax$diff_hi)),
  sprintf("- 동일 제품 P2 검정력: M0 %s ~ %s%%, M1 %s ~ %s%%(LLOQ 전체).", f2(min(pw[model == "M0", pass_pct])), f2(max(pw[model == "M0", pass_pct])), f2(min(pw[model == "M1", pass_pct])), f2(max(pw[model == "M1", pass_pct]))), "",
  sprintf("검증: LLOQ %s·잔차 그대로에서 개인 수준(대상자·요약), 절벽, 시험 수준 M0가 저장본과 같고, 잔차 비례 변형은 %s에서 원래와 같다.", fl(L0), fl(L0)))
writeLines(en, file.path(out_dir, "lloq_conclusion_en.md")); writeLines(ko, file.path(out_dir, "lloq_conclusion_ko.md"))
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", paste(en, collapse = "\n"))) stop("영문 결론에 한글·en/em dash·U+2212가 있습니다")
cat(en, sep = "\n")
