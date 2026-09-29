# S21 견고성 요약표: 결론(AUC0-last + Cmax 공동 1차, AUC0-inf 이차)이 가정을 바꿔도 유지되는가. 다섯 행: 구조 모델 두 개, 곡선 모양(Km·Vmax, 양 군),
# 분석법 LLOQ 0.02~0.5 mg/L, 비례 잔차 12%, 분석 모형(M0/M1/M2). 행마다 무엇을 계산했는지(모델, 칸, 시험 수)와 계산하지 않은 것을 적는다.
# 시험 모집단(건강인 60~90 kg, 체중 층화 배정, B0)만.
# 자료: results/oc_models/type1_models.csv(경계 1종 오류, 분석 모형별), results/curve_shape/curve_shape_B0.csv(곡선 모양, 2016 모델만, 대상자 20,000명씩),
#       results/oc/inversion_all.csv(시험약 Km 역산), results/lloq/lloq_trial_type1.csv(2016 모델, 경계 3칸, 5,000회), results/lloq/lloq_individual_table.csv,
#       results/trialpop/tp_residual_sensitivity.csv, results/trials/schedule_decision_resid12.csv(D2만), results/trials/schedules_per_endpoint_resid12.csv(B0, 500회).
# 주의: 잔차 12%에서는 경계 1종 오류(P2, AUC0-inf 구성)를 계산하지 않았다(미산출). 곡선 모양의 Vmax ×0.8(양 군)은 AUC0-last 노출 gate(±15%) 밖이지만
#       stress_test가 아니어서 보고서는 관측과 맞는 변형으로 센다(노트에 적는다; gate_note 열은 한국어라 영문 코드 variant로 거른다).
s21_mult <- function(variant, par) dcfg("scenarios.yaml", c("sensitivity_variants", variant, "theta_multipliers", par), sprintf("%s multiplier, variant %s (both arms)", par, variant),
                                         function(x) sub("\\.?0+$", "", fnum(as.numeric(x), 2)))
# 행 조건에 맞는 행 수를 "k/n"으로(둘 다 추적)
s21_kn <- function(rel, where, cond, item) sprintf("%s/%s", dcount(rel, sprintf("%s & %s", where, cond), item), dcount(rel, where, paste(item, "(denominator)")))

slide_S21 <- function() {
  T1 <- "oc_models/type1_models.csv"; CSH <- "curve_shape/curve_shape_B0.csv"; IA <- "oc/inversion_all.csv"
  LT <- "lloq/lloq_trial_type1.csv"; LI <- "lloq/lloq_individual_table.csv"; TRS <- "trialpop/tp_residual_sensitivity.csv"
  SR <- "trials/schedule_decision_resid12.csv"; PR <- "trials/schedules_per_endpoint_resid12.csv"; PB <- "trials/schedules_per_endpoint_base.csv"
  PC <- "oc_models/type1_paired_change.csv"; SS <- "oc_models/sd_se_models.csv"; PW <- "oc_models/power_models.csv"
  deck_slide("S21", tag = "sim")
  L <- DK$txt$S21
  nom <- f_nominal(); nom_num <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  lim <- unlist(.read("config/trial_design.yaml")$be$limits)
  V2C <- "pk_model=='k2020' & scenario=='V2_up_080'"

  # ---- 전제: 행마다 "결론 유지"가 기대는 사실 ----
  premise(nrow(rows(T1, "config=='P2' & pass_pct > 5 & !(pk_model=='k2020' & scenario=='V2_up_080')")) == 0, "P2 above 5% nowhere but the 2020 V2 cell, all analysis models")
  premise(nrow(rows(LT, "config=='P2' & pass_pct > 5")) == 0, "P2 at or below 5% in every LLOQ cell (both analysis models, both residual versions)")
  for (m in c("k2016", "k2020")) premise(nrow(rows(T1, sprintf("analysis_model=='M1' & config=='G2_Ai' & pk_model=='%s' & pass_pct > 5", m))) >= 4, paste("AUC0-inf rule A (i) + Cmax above 5% in several cells,", m))
  cs <- rows(CSH, "stress_test==FALSE")
  premise(setequal(cs$variant, c("base", "km05_both", "km2_both", "km5_both", "km10_both", "vmax080_both", "vmax125_both")), "non-stress curve-shape variants")
  premise(all(cs$coverage_lt80_pct == 0) && all(cs$extrap_true_p95 < 10), "no subject below 80% window coverage and small 95th-percentile extrapolation in every non-stress curve-shape variant")
  km <- c(rows(IA, "mechanism=='Km' & is.finite(end_auc_ratio)")$end_auc_ratio, rows(IA, "mechanism=='Km' & reachable==TRUE")$auc_ratio)
  premise(all(km > lim[1] & km < lim[2]), "Km inversion stays inside the equivalence limits")
  premise(nrow(rows(LT, "config %in% c('G2_Aii') & resid=='fixed' & lo <= 5")) == 0, "AUC0-inf rule A (ii) + Cmax: Wilson lower bound above 5% in every LLOQ cell, both analysis models")
  premise(all(rows(LI, "model=='k2016' & resid=='fixed'")$coverage_lt80_pct < 0.1), "window coverage below 80% under 0.1% at every LLOQ (2016 model)")
  premise(length(unique(rows(LT)$scenario)) == 3 && all(rows(LT)$n_trials == 5000) && !"pk_model" %in% names(rows(LT)), "LLOQ trials: 3 boundary cells, 5,000 trials, primary model only (file has no PK model column)")
  premise(row1(TRS, "variant=='resid12' & set=='iii'")$fail_pct > 20 && row1(TRS, "variant=='resid12' & set=='iv'")$fail_pct > 50, "with residual 12% many subjects still fail sets (iii) and (iv)")
  premise(nrow(rows(SR)) == 1 && rows(SR)$schedule == "D2" && !isTRUE(rows(SR)$recommend), "residual 12%: only D2 evaluated, not recommended")
  rf <- list.files(proj_path("results"), pattern = "resid12", recursive = TRUE)
  premise(!any(grepl("^(oc|oc_models|criteria|lloq)/|type1", rf)) && all(rows(LT)$resid %in% c("fixed", "scaled")),
          "no boundary type I error was computed with the 12% residual (resid in the LLOQ file is the additive-error version, not the proportional residual)")
  for (am in c("M0", "M1", "M2")) { r <- rows(T1, sprintf("analysis_model=='%s' & config=='P2' & class!='conservative'", am))
    premise(nrow(r) == 1 && r$pk_model == "k2020" && r$scenario == "V2_up_080", paste("the only non-conservative P2 cell is the 2020 V2 cell,", am)) }

  # ---- 제목 ----
  s12 <- dv(TRS, "variant=='resid12' & set=='iii'", "sigma_prop_pct", 0, "%", "proportional residual, variant")
  deck_kicker(tx("S21.kicker")); deck_title(tx("S21.title", list(nom = nom)))

  # ---- 표 ----
  pmax <- function(m, am = "M1") dext(T1, sprintf("analysis_model=='%s' & config=='P2' & pk_model=='%s'", am, m), "pass_pct", max, 2, "%", sprintf("largest P2 boundary type I error, %s, %s", m, am))
  gcnt <- function(m) dcount(T1, sprintf("analysis_model=='M1' & config=='G2_Ai' & pk_model=='%s' & pass_pct > 5", m), sprintf("M1 G2_Ai cells above 5%%, %s", m))
  npm <- dcount(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016'", "boundary cells per model")
  wcnt <- function(m) dcount(T1, sprintf("analysis_model=='M1' & config=='G2_Ai' & pk_model=='%s' & lo > 5", m), sprintf("M1 G2_Ai cells with Wilson lower bound above 5%%, %s", m))
  r1 <- list(npm = npm, n10 = f_reps("boundary"), n20 = f_reps("ext"), p16 = pmax("k2016"), p20 = pmax("k2020"), g16 = gcnt("k2016"), g20 = gcnt("k2020"),
             w16 = wcnt("k2016"), w20 = wcnt("k2020"), nom = nom)
  premise(identical(unique(rows(T1, "n_trials != 10000")[, paste(pk_model, scenario)]), "k2020 V2_up_080") && all(rows(T1, V2C)$n_trials == 20000) && all(rows(T1, sprintf("!(%s)", V2C))$n_trials == 10000),
          "only the 2020 V2 cell was extended (20,000 trials); all other cells 10,000 (rows 1 and 5)")
  premise(row1(T1, sprintf("analysis_model=='M1' & config=='P2' & %s", V2C))$pass_pct == max(rows(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020'")$pass_pct), "the 2020 maximum is the V2 cell (row 1)")
  mfmt <- function(x) sub("\\.?0+$", "", fnum(as.numeric(x), 2))
  kmr <- dcfg("oc_design.yaml", c("mechanisms", "Km", "range"), "Km multiplier search range, test arm", function(x) sprintf("%s~×%s", mfmt(x[1]), mfmt(x[2])))
  km_rng <- local({ r <- rows(IA, "mechanism=='Km' & is.finite(end_auc_ratio)"); r2 <- rows(IA, "mechanism=='Km' & reachable==TRUE"); x <- range(c(r$end_auc_ratio, r2$auc_ratio))
    dderived("true AUC0-inf ratio range over Km x0.01 to x100 (range ends and reached rows)", IA, "mechanism=='Km' :: end_auc_ratio, auc_ratio", x, rng_fmt(x[1], x[2], 3)) })
  thr80 <- dderived("window coverage threshold in column name coverage_lt80_pct (percent)", CSH, "column name coverage_lt80_pct", 80, "80")
  r2 <- list(k1 = s21_mult("km05_both", "Km"), k2 = s21_mult("km10_both", "Km"), v1 = s21_mult("vmax080_both", "Vmax"), v2 = s21_mult("vmax125_both", "Vmax"),
             ncs = dint(CSH, "variant=='base'", "n", "subjects per curve-shape variant"), thr = thr80,
             lt = drange(CSH, "stress_test==FALSE", "coverage_lt80_pct", 0, "%", "share below 80% window coverage, non-stress curve-shape variants"),
             p95 = drange(CSH, "stress_test==FALSE", "extrap_true_p95", 2, "%", "95th percentile of true extrapolation, non-stress curve-shape variants"),
             kmr = kmr, km = km_rng)
  LW <- "resid=='fixed'"
  r3 <- list(lg = f_lloq_grid(), nl = as.character(length(unlist(.read("config/assay.yaml")$lloq_sensitivity_mg_L))),
             nc = as.character(length(unique(rows(LT)$scenario))), nt = dint(LT, "model=='M1' & resid=='fixed' & config=='P2' & scenario=='F_down_080' & abs(lloq - 0.078) < 1e-9", "n_trials", "LLOQ trials per cell"),
             p1 = drange(LT, "model=='M1' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid M1 P2 range"),
             p0 = drange(LT, "model=='M0' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid M0 P2 range"),
             k1 = s21_kn(LT, "model=='M1' & resid=='fixed' & config=='P2'", "lo > 5", "LLOQ grid M1 P2 cells with Wilson lower bound above 5%"),
             g1 = drange(LT, "model=='M1' & resid=='fixed' & config=='G2_Aii'", "pass_pct", 2, "%", "LLOQ grid M1 G2_Aii range"),
             kg = s21_kn(LT, "model=='M1' & resid=='fixed' & config=='G2_Aii'", "lo > 5", "LLOQ grid M1 G2_Aii cells with Wilson lower bound above 5%"), nom = nom)
  # nl, nc는 config 격자 길이와 파일의 칸 수: 추적 행을 남긴다
  dderived("number of LLOQ values in the sensitivity grid", "config/assay.yaml", "length(lloq_sensitivity_mg_L)", as.integer(r3$nl), r3$nl)
  dderived("boundary cells in the LLOQ trials", LT, "length(unique(scenario))", as.integer(r3$nc), r3$nc)
  r4 <- list(s12 = s12, s24 = dv(TRS, "variant=='k2016' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model"),
             i24 = dv(TRS, "variant=='k2016' & set=='i'", "fail_pct", 1, "%", "set (i) failing, residual 24.2%"), i12 = dv(TRS, "variant=='resid12' & set=='i'", "fail_pct", 1, "%", "set (i) failing, residual 12%"),
             f24 = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, residual 24.2%"), f12 = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, residual 12%"),
             g24 = dv(TRS, "variant=='k2016' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, residual 24.2%"), g12 = dv(TRS, "variant=='resid12' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, residual 12%"),
             ntr = dint(PR, "scenario=='S00' & schedule=='B0' & endpoint=='AUClast'", "n_trials", "trials per scenario, residual 12%"))
  am_rng <- function(am) drange(T1, sprintf("analysis_model=='%s' & config=='P2'", am), "pass_pct", 2, "%", sprintf("%s P2 boundary range", am))
  am_c <- function(am) dcount(T1, sprintf("analysis_model=='%s' & config=='P2' & class=='conservative'", am), sprintf("%s P2 cells conservative", am))
  nall <- dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells")
  v2 <- function(am) dv(T1, sprintf("%s & analysis_model=='%s' & config=='P2'", V2C, am), "pass_pct", 2, "%", sprintf("P2, 2020 V2 cell, %s", am))
  r5 <- list(n = nall, r0 = am_rng("M0"), r1 = am_rng("M1"), r2 = am_rng("M2"), c0 = am_c("M0"), c1 = am_c("M1"), c2 = am_c("M2"), n10 = r1$n10, n20 = r1$n20,
             m0 = v2("M0"), m1 = v2("M1"), m2 = v2("M2"))
  r5$c <- r5$c1
  premise(r5$c0 == r5$c1 && r5$c1 == r5$c2 && as.integer(r5$c1) == as.integer(nall) - 1L, "all boundary cells but the V2 cell are conservative under M0, M1 and M2 (row 5)")
  cls <- vapply(c("M0", "M1", "M2"), function(am) row1(T1, sprintf("%s & analysis_model=='%s' & config=='P2'", V2C, am))$class, "")
  premise(identical(unname(cls), c("nominal", "exceeding", "exceeding")), "V2 cell class: nominal under M0, exceeding under M1 and M2 (row 5)")
  for (am in c("M0", "M1", "M2")) premise(nrow(rows(T1, sprintf("analysis_model=='%s' & config=='G2_Ai' & lo > 5", am))) > 0, paste("AUC0-inf rule A (i) + Cmax above 5% (Wilson lower bound) in some cells under", am, "(legend)"))
  premise(nrow(rows(LT, "config=='G2_Aii' & lo > 5")) > 0 && nrow(rows(TRS, "variant=='resid12' & fail_pct > 0")) > 0, "AUC0-inf exceedance (LLOQ) and reliability failures (residual 12%) remain (legend)")
  R <- list(r1, r2, r3, r4, r5)
  df <- data.frame(a = vapply(1:5, function(i) fill(L$table$cond[[i]], R[[i]]), ""), b = vapply(1:5, function(i) fill(L$table$scope[[i]], R[[i]]), ""),
                   c = vapply(1:5, function(i) fill(L$table$result[[i]], R[[i]]), ""), d = unlist(L$table$verdict), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S21.table.head")
  TH <- 3.66
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, GEO$CW, TH), widths = c(1.95, 3.45, 5.05, 1.78), size = 12, align_num = FALSE, label = "table_robust")

  # ---- 아래: 요점 ----
  BY <- GEO$BODY_TOP + TH + 0.3
  deck_bullets(tx("S21.bullets", list(nom = nom, s12 = s12, nc = r3$nc)), box = c(GEO$ML, BY, GEO$CW, GEO$BODY_BOTTOM - BY), size = 16, gap_pt = 6)

  # ---- 노트 ----
  nts <- list(
    nom = nom, npm = npm, n10 = r1$n10, n20 = r1$n20,
    p16c = dci(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "P2 largest, 2016 model, M1"),
    p20c = dci(T1, sprintf("analysis_model=='M1' & config=='P2' & %s", V2C), "pass_pct", "lo", "hi", 2, "%", "P2 largest, 2020 model, M1"),
    m016 = pmax("k2016", "M0"), m020 = pmax("k2020", "M0"),
    gm16 = dext(T1, "analysis_model=='M1' & config=='G2_Ai' & pk_model=='k2016'", "pass_pct", max, 2, "%", "largest G2_Ai, 2016 model, M1"),
    gm20 = dext(T1, "analysis_model=='M1' & config=='G2_Ai' & pk_model=='k2020'", "pass_pct", max, 2, "%", "largest G2_Ai, 2020 model, M1"),
    fi16 = dv("trialpop/tp_failure_by_set.csv", "pk_model=='k2016' & set=='i'", "fail_pct", 1, "%", "set (i) failing, 2016 model"),
    fi20 = dv("trialpop/tp_failure_by_set.csv", "pk_model=='k2020' & set=='i'", "fail_pct", 1, "%", "set (i) failing, 2020 model"),
    lt = r2$lt, ncs = r2$ncs, k1 = r2$k1, k2 = r2$k2, v1 = r2$v1, v2 = r2$v2, thr = thr80, p95 = r2$p95,
    p95k = dv(CSH, "variant=='km10_both'", "extrap_true_p95", 2, "%", "95th percentile true extrapolation, Km x10"),
    p95v = dv(CSH, "variant=='vmax125_both'", "extrap_true_p95", 2, "%", "95th percentile true extrapolation, Vmax x1.25"),
    sx = s21_mult("vmax050_both", "Vmax"), sp95 = dv(CSH, "variant=='vmax050_both'", "extrap_true_p95", 2, "%", "95th percentile true extrapolation, Vmax x0.5 stress test"),
    slt = dv(CSH, "variant=='vmax050_both'", "coverage_lt80_pct", 3, "%", "share below 80% coverage, Vmax x0.5 stress test"),
    gate = dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "exposure gate: AUClast mean within +/- tolerance (%)", num_fmt(0)),
    gv = dv(CSH, "variant=='vmax080_both'", "AUClast_ratio_vs_obs544", 2, "", "AUClast simulated / observed, Vmax x0.8 both arms"),
    gs = dv(CSH, "variant=='vmax050_both'", "AUClast_ratio_vs_obs544", 2, "", "AUClast simulated / observed, Vmax x0.5 stress test"),
    kmr = kmr, km = km_rng,
    lg = r3$lg, nl = r3$nl, nc = r3$nc, nt = r3$nt, lp1 = r3$p1, lp0 = r3$p0, lk1 = r3$k1,
    lk0 = s21_kn(LT, "model=='M0' & resid=='fixed' & config=='P2'", "lo > 5", "LLOQ grid M0 P2 cells with Wilson lower bound above 5%"),
    lg1 = r3$g1, lg0 = drange(LT, "model=='M0' & resid=='fixed' & config=='G2_Aii'", "pass_pct", 2, "%", "LLOQ grid M0 G2_Aii range"),
    lkg0 = s21_kn(LT, "model=='M0' & resid=='fixed' & config=='G2_Aii'", "lo > 5", "LLOQ grid M0 G2_Aii cells with Wilson lower bound above 5%"),
    lb1 = drange(LT, "model=='M1' & resid=='fixed' & config=='G2_B'", "pass_pct", 2, "%", "LLOQ grid M1 G2_B range"),
    lc = drange(LI, "model=='k2016' & resid=='fixed'", "coverage_lt80_pct", 2, "%", "LLOQ grid, share below 80% window coverage, 2016 model"),
    lr = drange(LI, "model=='k2016' & resid=='fixed'", "reliable_no_span_pct", 1, "%", "LLOQ grid, reliability set (i), 2016 model"),
    lloq = f_lloq(),
    s12 = r4$s12, s24 = r4$s24, f24 = r4$f24, f12 = r4$f12, g24 = r4$g24, g12 = r4$g12,
    i24 = r4$i24, i12 = r4$i12, g16 = r1$g16, g20 = r1$g20, w16 = r1$w16, w20 = r1$w20, m0 = r5$m0, m1 = r5$m1, m2 = r5$m2,
    gt = dderived("M1 G2_Ai cells above 5% (point estimate), both models", T1, "analysis_model=='M1' & config=='G2_Ai' & pass_pct > 5 :: count over k2016 and k2020", as.integer(r1$g16) + as.integer(r1$g20), as.character(as.integer(r1$g16) + as.integer(r1$g20))),
    wt = dderived("M1 G2_Ai cells with Wilson lower bound above 5%, both models", T1, "analysis_model=='M1' & config=='G2_Ai' & lo > 5 :: count over k2016 and k2020", as.integer(r1$w16) + as.integer(r1$w20), as.character(as.integer(r1$w16) + as.integer(r1$w20))),
    ra = dv(SR, "schedule=='D2'", "a_mean_width_rel_decrease", 2, "%", "residual 12%, D2, criterion (a) (%)", 100), rc = dv(SR, "schedule=='D2'", "c_reliable_gain_pp", 2, "", "residual 12%, D2, criterion (c) (points)"),
    rd = dv(SR, "schedule=='D2'", "d_extrap20_ratio", 2, "", "residual 12%, D2, criterion (d) ratio"),
    ntr = r4$ntr, s00 = dv(PR, "scenario=='S00' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate S00, residual 12%"),
    ke = dv(PR, "scenario=='KE110' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate KE110, residual 12%"),
    f90 = dv(PR, "scenario=='F090' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate F090, residual 12%"),
    kem = dcfg("scenarios.yaml", c("scenarios", "KE110", "T_multipliers", "ke"), "ke multiplier, scenario KE110", num_fmt(2)),
    fm = dcfg("scenarios.yaml", c("scenarios", "F090", "T_multipliers", "F"), "F multiplier, scenario F090", num_fmt(2)),
    s00b = dv(PB, "scenario=='S00' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate S00, base"),
    keb = dv(PB, "scenario=='KE110' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate KE110, base"),
    f90b = dv(PB, "scenario=='F090' & schedule=='B0' & endpoint=='AUClast'", "pass_rate", 1, "%", "AUC0-last pass rate F090, base"),
    n = nall, r0 = r5$r0, r1 = r5$r1, r2 = r5$r2, c0 = r5$c0, c1 = r5$c1, c2 = r5$c2,
    dpc = drange(PC, "config=='P2' & comparison=='M1 minus M0'", "diff_pp", 2, "", "P2, M1 minus M0, paired (points)"),
    npos = dcount(PC, "config=='P2' & comparison=='M1 minus M0' & diff_pp > 0", "cells where P2 rises under M1"),
    sd0 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M0' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE AUCinf_true M0"),
    sd1 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M1' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE AUCinf_true M1"),
    pw0 = dci(PW, "pk_model=='k2016' & scenario=='S00' & analysis_model=='M0' & config=='P2'", "pass_pct", "lo", "hi", 1, "%", "power k2016 S00 M0 P2"),
    pw1 = dci(PW, "pk_model=='k2016' & scenario=='S00' & analysis_model=='M1' & config=='P2'", "pass_pct", "lo", "hi", 1, "%", "power k2016 S00 M1 P2"))
  premise(identical(range(rows(CSH, "stress_test==FALSE & variant!='vmax080_both'")$extrap_true_p95), range(cs$extrap_true_p95)) && all(rows(CSH, "stress_test==FALSE & variant!='vmax080_both'")$coverage_lt80_pct == 0),
          "the curve-shape ranges are the same without Vmax x0.8 (notes)")
  premise(row1(CSH, "variant=='vmax080_both'")$stress_test == FALSE && abs(row1(CSH, "variant=='vmax080_both'")$AUClast_ratio_vs_obs544 - 1) * 100 > as.numeric(.read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct),
          "Vmax x0.8 (both arms) lies outside the AUClast exposure gate but is not flagged as a stress test (notes)")
  deck_notes(tx("S21.notes", nts))
  deck_end()
}
