# S24 예상 질의와 답: 예상 질의응답 문서(regulatory/src/FDA_questions.Rmd)에서 임상약리·임상개발에 가장 중요한 다섯 문항(Q1, Q7, Q19+Q11, Q4+Q3, Q12+Q17)을
# 한국어로 줄였다. 왼쪽 질문 상자 + 오른쪽 답(2열). [문헌+모의]
# M1 수치는 results/oc_models/type1_models.csv에서 읽는다(FDA_questions는 Q1, Q8, Q9에서 oc/g2_rules_flags.csv의 M0 수를 인쇄한다).
S24_T1 <- "oc_models/type1_models.csv"; S24_DE <- "oc_models/p2_decomposition_models.csv"; S24_CS <- "cliff/cliff_summary.csv"; S24_CP <- "cliff/cliff_points.csv"
s24_w <- function(am, cf, extra = "") sprintf("analysis_model=='%s' & config=='%s'%s", am, cf, extra)
S24_V2 <- function(am) sprintf("pk_model=='k2020' & mechanism=='V2' & analysis_model=='%s'", am)
S24_CPW <- function(timing) sprintf("model=='k2016' & weight=='base' & timing=='%s' & definition_day==1 & schedule=='daily_29_57'", timing)
# 매일 채혈 구간(연구일): config/oc_design.yaml cliff.schedules.daily_29_57의 "a..b"(투여 후 일) + 1
s24_daily <- function() {
  s <- .read("config/oc_design.yaml")$cliff$schedules$daily_29_57; m <- regmatches(s, regexec("([0-9]+)\\.\\.([0-9]+)", s))[[1]]
  premise(length(m) == 3, "daily schedule range a..b in config"); sd <- as.numeric(m[2:3]) + 1
  dderived("daily sampling window of the cliff check (study days = days after dose + 1)", "config/oc_design.yaml", "cliff.schedules.daily_29_57 :: regex a..b (days after dose) + 1", sd, rng_fmt(sd[1], sd[2], 0))
}
# B0 마지막 채혈 연구일(s04_study_days와 같은 계산; config 목록이 정수·실수 혼합이라 unlist)
s24_last <- function() {
  d <- unlist(.read("config/trial_design.yaml")$schedules$B0$days); premise(length(d) > 0, "schedule B0"); sd <- d + 1
  dderived("schedule B0, last (study day = days after dose + 1)", "config/trial_design.yaml", "schedules.B0.days :: last", sd, fnum(max(sd), 0))
}
s24_max <- function(am, cf, item) {
  r <- rows(S24_T1, s24_w(am, cf)); premise(nrow(r) == 16, sprintf("%s %s: 16 boundary cells", am, cf)); rr <- r[which.max(pass_pct)]
  dci(S24_T1, s24_w(am, cf, sprintf(" & pk_model=='%s' & scenario=='%s'", rr$pk_model, rr$scenario)), "pass_pct", "lo", "hi", 2, "%", item)
}

slide_S24 <- function() {
  T1 <- S24_T1; DE <- S24_DE; CS <- S24_CS; CP <- S24_CP; TPF <- "trialpop/tp_failure_by_set.csv"; PC <- "rationale/pillar1_coverage_B0.csv"
  TP <- "sample_size/ss_table_power.csv"; SS <- "oc_models/sd_se_models.csv"
  deck_slide("S24", tag = "litsim")

  # ---- 전제 ----
  p2 <- rows(T1, s24_w("M1", "P2")); top <- p2[which.max(pass_pct)]
  premise(nrow(p2) == 16 && sum(p2$class == "exceeding") == 1 && top$pk_model == "k2020" && top$mechanism == "V2", "one exceeding P2 cell under M1, the 2020 model V2 cell (Q7)")
  premise(rows(DE, S24_V2("M1"))$auclast_bias_dir == "toward_1" && rows(DE, S24_V2("M1"))$auclast_bias_pct > 0, "AUC0-last GMR biased toward 1 in the exceeding cell (Q7)")
  premise(rows(DE, S24_V2("M0"))$class == "nominal", "the same cell is nominal under M0 (Q7)")
  sd_ <- rbindlist(lapply(c("base", "struct2020"), function(k) rows(sprintf("trials/schedule_decision_%s.csv", k))[, variant := k]))
  premise(nrow(sd_) > 0 && !any(sd_$recommend %in% TRUE), "the pre-specified rule recommended no added samples in either model (Q4)")
  premise(rows(CS, "model=='k2016' & weight=='base'")$lloq_studyday_median < max(unlist(.read("config/trial_design.yaml")$schedules$B0$days)) + 1, "median LLOQ day before the last B0 sample (Q3)")

  f <- list(nom = f_nominal(), r2i = f_set("i", "r2"), r2iii = f_set("iii", "r2"), last = s24_last(),
            # Q1
            ext = dv(PC, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolation, 2016 model (%)"),
            fi = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
            ncell = dcount(T1, s24_w("M1", "G2_Ai"), "boundary cells per configuration (M1)"),
            ga = dcount(T1, s24_w("M1", "G2_Ai", " & pass_pct > 5"), "M1 G2_Ai cells with point estimate above 5%"),
            # Q7
            c1 = dcount(T1, s24_w("M1", "P2", " & class=='conservative'"), "M1 P2 cells classified conservative"),
            e1 = dcount(T1, s24_w("M1", "P2", " & class=='exceeding'"), "M1 P2 cells classified exceeding"),
            r1 = drange(T1, s24_w("M1", "P2"), "pass_pct", 2, "%", "M1 P2 range"),
            mult = dv(DE, S24_V2("M1"), "multiplier", 2, "", "2020 V2 cell multiplier"),
            p2ci = dci(DE, S24_V2("M1"), "p2_pct", "lo", "hi", 2, "%", "2020 V2 cell P2 M1"),
            bias = dv(DE, S24_V2("M1"), "auclast_bias_pct", 2, "%", "2020 V2 AUC0-last bias M1"),
            m0 = dv(DE, S24_V2("M0"), "p2_pct", 2, "%", "2020 V2 cell P2 M0"),
            # Q19, Q11
            fiii = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models"),
            # Q4, Q3
            len = dv(CS, "model=='k2016' & weight=='base'", "len1_median", 2, "", "cliff length median, 2016 model (days)"),
            dd = s24_daily(),
            k3 = dderived("sample-count threshold in column name pct_ge3", CP, "column name pct_ge3 (share of subjects with 3 or more samples in the cliff)", 3, "3"),
            pct = dv(CP, S24_CPW("nominal"), "pct_ge3", 1, "%", "daily sampling, 3 or more samples on the cliff, nominal days (%)"),
            day16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_median", 1, "", "LLOQ study day median, 2016 model"),
            d58 = dderived("study-day threshold in column name lloq_after_day58_pct", CS, "column name lloq_after_day58_pct (100 x mean(t_lloq + 1 > 58))", 58, "58"),
            a16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_after_day58_pct", 1, "%", "above LLOQ after Day 58, 2016 model (%)"),
            # Q12, Q17
            cvb = dint(TP, "input_model=='k2016' & cv==43 & gmr==0.95 & n==117 & analysis_model=='M1'", "cv", "protocol CV (%)"),
            cvs = dint(TP, "input_model=='k2016' & cv==50 & gmr==0.95 & n==117 & analysis_model=='M1'", "cv", "sensitivity CV (%)"),
            g = dv(TP, "input_model=='k2016' & cv==43 & gmr==0.95 & n==117 & analysis_model=='M1'", "gmr", 2, "", "true GMR of the power table"),
            n = f_n_arm(),
            p43 = dv(TP, "input_model=='k2016' & cv==43 & gmr==0.95 & n==117 & analysis_model=='M1'", "analytic_pct", 1, "%", "power n117 CV 43 M1"),
            p43m0 = dv(TP, "input_model=='k2016' & cv==43 & gmr==0.95 & n==117 & analysis_model=='M0'", "analytic_pct", 1, "%", "power n117 CV 43 M0"),
            p50 = dv(TP, "input_model=='k2016' & cv==50 & gmr==0.95 & n==117 & analysis_model=='M1'", "analytic_pct", 1, "%", "power n117 CV 50 M1"),
            sd0 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M0' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE M0"),
            sd1 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M1' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE M1"))
  premise(nrow(rows(TP, "input_model=='k2016' & cv==43 & gmr==0.95 & n==117 & analysis_model=='M1'")) == 1 &&
          as.numeric(.read("config/prereg_20260926.yaml")$section3$base_cv_pct) == 43 && as.numeric(.read("config/prereg_20260926.yaml")$section3$sensitivity_cv_pct) == 50 &&
          as.numeric(.read("config/trial_design.yaml")$n_per_arm) == 117, "power rows are the protocol CV, sensitivity CV and protocol evaluable n (Q12)")
  deck_kicker(tx("S24.kicker")); deck_title(tx("S24.title", f))

  # ---- 2열: 질문 상자 | 답 ----
  qs <- DK$txt$S24$qa; nq <- length(qs); premise(nq == 5, "five questions")
  gap <- 0.1; y0 <- GEO$BODY_TOP; rh <- (GEO$BODY_BOTTOM - y0 - (nq - 1) * gap) / nq; qw <- 3.1; ag <- 0.2
  for (i in seq_len(nq)) {
    y <- y0 + (i - 1) * (rh + gap); q <- qs[[i]]
    deck_text(sprintf("__%s__  %s", q$id, fill(q$q, f, sprintf("S24.qa.%d.q", i))), c(GEO$ML, y, qw, rh), size = 16, bold = FALSE, bg = PAL$tint_blue, geom = "roundRect",
              label = sprintf("q%d", i), gap_pt = 0)
    deck_text(fill(q$a, f, sprintf("S24.qa.%d.a", i)), c(GEO$ML + qw + ag, y, GEO$CW - qw - ag, rh), size = 16, label = sprintf("a%d", i), gap_pt = 0)
  }

  # ---- 노트 ----
  EXf <- "oc_models/extension_decision_models.csv"; PW <- "oc_models/power_models.csv"; PCf <- "oc_models/type1_paired_change.csv"; NN <- "sample_size/ss_table_n_needed.csv"
  nn_ <- function(cv, am, col = "n_evaluable_per_arm") dint(NN, sprintf("cv==%s & gmr==0.95 & analysis_model=='%s' & target_pct==90", cv, am), col, sprintf("%s %s CV %s target 90", col, am, cv))
  deck_notes(tx("S24.notes", c(f, list(
    p95 = dv(PC, "model=='k2016' & group=='all'", "extrap_true_p95", 2, "%", "95th percentile true extrapolation, 2016 model (%)"),
    nca = dv(PC, "model=='k2016' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolation, 2016 model (%)"),
    gaii = dcount(T1, s24_w("M1", "G2_Aii", " & pass_pct > 5"), "M1 G2_Aii cells with point estimate above 5%"),
    gaiimax = s24_max("M1", "G2_Aii", "M1 G2_Aii maximum"),
    ga0 = dcount(T1, s24_w("M0", "G2_Ai", " & pass_pct > 5"), "M0 G2_Ai cells with point estimate above 5%"),
    tgt = dv(T1, s24_w("M1", "P2", " & pk_model=='k2020' & scenario=='V2_up_080'"), "target", 2, "", "true AUC0-inf ratio target, exceeding cell"),
    m16 = dv(T1, s24_w("M1", "P2", " & pk_model=='k2016' & scenario=='V2_up_080'"), "multiplier", 2, "", "2016 model V2 multiplier for the same target"),
    ntr = dint(DE, S24_V2("M1"), "n_trials", "2020 V2 cell trials M1"),
    pre = dci(EXf, "pk_model=='k2020' & scenario=='V2_up_080' & analysis_model=='M1'", "pass_pct", "lo", "hi", 2, "%", "M1 P2 V2 cell at the pre-specified 10,000 trials"),
    ref = dv(DE, S24_V2("M1"), "ref_pct", 2, "%", "2020 V2 unbiased reference M1"),
    al = dv(DE, S24_V2("M1"), "auclast_pct", 2, "%", "2020 V2 AUC0-last alone M1"),
    cm = dv(DE, S24_V2("M1"), "p2_minus_auclast_pp", 2, "", "2020 V2 Cmax term M1 (points)"),
    m0ci = dci(DE, S24_V2("M0"), "p2_pct", "lo", "hi", 2, "%", "2020 V2 cell P2 M0"),
    fiv = drange(TPF, "set=='iv'", "fail_pct", 1, "%", "trial population, set (iv) failing, two models"),
    day20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_median", 1, "", "LLOQ study day median, 2020 model"),
    a20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_after_day58_pct", 1, "%", "above LLOQ after Day 58, 2020 model (%)"),
    pctw = dv(CP, S24_CPW("windowed"), "pct_ge3", 1, "%", "daily sampling, 3 or more samples on the cliff, windowed times (%)"),
    p50m0 = dv(TP, "input_model=='k2016' & cv==50 & gmr==0.95 & n==117 & analysis_model=='M0'", "analytic_pct", 1, "%", "power n117 CV 50 M0"),
    tg = dint(NN, "cv==43 & gmr==0.95 & analysis_model=='M1' & target_pct==90", "target_pct", "target power (%)"),
    nb1 = nn_(43, "M1"), nb0 = nn_(43, "M0"), ns1 = nn_(50, "M1"), ns0 = nn_(50, "M0"),
    pw1 = dci(PW, "pk_model=='k2016' & scenario=='S00' & analysis_model=='M1' & config=='P2'", "pass_pct", "lo", "hi", 1, "%", "power k2016 S00 M1 P2"),
    pw0 = dci(PW, "pk_model=='k2016' & scenario=='S00' & analysis_model=='M0' & config=='P2'", "pass_pct", "lo", "hi", 1, "%", "power k2016 S00 M0 P2"),
    pc = drange(PCf, "comparison=='M1 minus M0' & config=='P2'", "diff_pp", 2, "", "P2 M1 minus M0, range (points)"),
    pcm = local({ r <- rows(PCf, "comparison=='M1 minus M0' & config=='P2'"); premise(nrow(r) == 16, "16 paired changes")
      dderived("P2 M1 minus M0, mean (points)", PCf, "comparison=='M1 minus M0' & config=='P2' :: mean(diff_pp)", mean(r$diff_pp), fnum(mean(r$diff_pp), 2)) })))))
  deck_end()
}
