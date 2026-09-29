# S08 논리 구조 도식: 전제(말기 절벽) → ① → ② → ③ → 결론. 둥근 사각형 다섯 개와 오른쪽 화살표, 단계마다 대표 수치 하나(headline).
# 시험 모집단(건강인, 체중 층화, B0)만. 주분석 M1, M0은 노트에 병기.
# 글자 없는 도형: deck_box는 8pt 자리 글자를 넣어 검사 6(글자 크기)에 걸리므로 16pt 빈 글상자로 그린다
s08_panel <- function(box, fill, geom = "roundRect", label = "panel") deck_text(" ", box, size = 16, bg = fill, geom = geom, label = label)
# f_study_days()와 같은 계산(config 채혈일이 정수·실수 혼합이라 yaml이 목록으로 읽으므로 unlist; 공용 함수 수정 요청)
s08_study_days <- function(schedule = "B0", which = c("all", "last", "n")) {
  which <- match.arg(which); d <- unlist(.read("config/trial_design.yaml")$schedules[[schedule]]$days); premise(length(d) > 0, paste("schedule", schedule))
  sd <- d + 1; p <- switch(which, all = paste(fnum(sd[sd == round(sd)], 0), collapse = ", "), last = fnum(max(sd), 0), n = as.character(length(d)))
  dderived(sprintf("schedule %s, %s (study day = days after dose + 1)", schedule, which), "config/trial_design.yaml", sprintf("schedules.%s.days :: %s", schedule, which), sd, p)
}
s08_inst <- function(am, what = c("median", "max")) {
  what <- match.arg(what); CSf <- "criteria/criteria_instability.csv"
  r <- rows(CSf, sprintf("analysis_model=='%s' & !scenario %%in%% c('S00','F097')", am)); premise(nrow(r) == 16, "16 boundary cells in the instability file")
  x <- if (what == "median") median(r$inst_all) else max(r$inst_all)
  dderived(sprintf("decision instability inst_all, %s over 16 boundary cells, %s", what, am), CSf,
           sprintf("analysis_model=='%s' & !scenario in (S00, F097) :: %s of inst_all", am, what), x, paste0(fnum(x, 1), "%"))
}

# 창 포착률 최솟값: "이상"이라고 쓰므로 0.1 단위로 내림(보고서 tpcv_min과 같은 파일·조건·계산)
s08_cov_min <- function() {
  TCV <- "trialpop/tp_coverage_individual.csv"; w <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  r <- rows(TCV, w); premise(nrow(r) == 2, "window coverage: two models"); x <- min(r$min) * 100
  dderived("window coverage (true AUC0-tlast / true AUC0-inf), smallest subject-level value over both models", TCV, sprintf("%s :: min(min)", w), x, sprintf("%s%%", fnum(floor(x * 10) / 10, 1)))
}

slide_S08 <- function() {
  CS <- "cliff/cliff_summary.csv"; TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; TCH <- "trialpop/tp_characteristics.csv"
  CGf <- "criteria/criteria_g2_type1.csv"; T1f <- "oc_models/type1_models.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"; CB <- "model %in% c('k2016','k2020') & weight=='base'"
  deck_slide("S08", tag = "litsim")
  deck_kicker(tx("S08.kicker")); deck_title(tx("S08.title"))

  # 전제 검사
  premise(nrow(rows(CS, CB)) == 2, "cliff summary: two models, trial weight range only (base rows)")
  premise(all(rows(TCH, "set=='i'")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUC0-inf (text: not random)")
  premise(sum(rows(CGf, "analysis_model=='M1' & config=='G2_A_i'")$pass_pct > 5) > sum(rows(CGf, "analysis_model=='M1' & config=='G2_C_i'")$pass_pct > 5),
          "rule A (i) exceeds 5% in more boundary cells than rule C (i) under M1 (text: decision depends on the rule)")
  p2 <- rows(T1f, "analysis_model=='M1' & config=='P2'"); premise(nrow(p2) == 16 && sum(p2$class == "exceeding") == 1, "one exceeding P2 cell under M1 (caption)")
  pm <- p2[which.max(pass_pct)]; premise(pm$pk_model == "k2020" && pm$scenario == "V2_up_080", "the exceeding P2 cell under M1 is the 2020 model V2 up (notes)")

  f <- list(
    nom = f_nominal(), ns = s08_study_days("B0", "n"), last = s08_study_days("B0", "last"),
    len = headline(drange(CS, CB, "len1_median", 2, "", "cliff length (1-day definition), median, two models")),
    iii = headline(drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models")),
    i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
    a = headline(dcount(CGf, "analysis_model=='M1' & config=='G2_A_i' & pass_pct > 5", "M1 G2_A_i cells above 5% (point)")),
    ncell = dcount(CGf, "analysis_model=='M1' & config=='G2_A_i'", "boundary cells per configuration (M1)"),
    c = dcount(CGf, "analysis_model=='M1' & config=='G2_C_i' & pass_pct > 5", "M1 G2_C_i cells above 5% (point)"),
    cov = headline(drange(TCV, MW, "median", 1, "%", "window coverage, median, two models", scale = 100)),
    cmin = s08_cov_min(),
    n_cons = headline(dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='conservative'", "M1 P2 cells classified conservative")),
    n_all = dcount(T1f, "analysis_model=='M1' & config=='P2'", "M1 P2 boundary cells"),
    n_exc = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='exceeding'", "M1 P2 cells classified exceeding"),
    p2max = dv(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2 maximum"))

  # 도식: 상자 다섯 개 + 화살표 네 개
  steps <- c("premise", "a1", "a2", "a3", "concl")
  fills <- c(PAL$tint_grey, PAL$tint_blue, PAL$tint_blue, PAL$tint_blue, PAL$tint_orange)
  cols <- c(PAL$ink2, PAL$blue, PAL$blue, PAL$blue, PAL$orange)
  aw <- 0.26; ag <- 0.04; bw <- (GEO$CW - 4 * (aw + 2 * ag)) / 5; y <- GEO$BODY_TOP + 0.08; bh <- GEO$BODY_BOTTOM - y - 0.12
  for (k in seq_along(steps)) {
    s_ <- steps[k]; x <- GEO$ML + (k - 1) * (bw + aw + 2 * ag); xi <- x; wi <- bw
    s08_panel(c(x, y, bw, bh), fills[k], "roundRect", sprintf("box_%s", s_))
    deck_text(tx(sprintf("S08.steps.%s.label", s_)), c(xi, y + 0.1, wi, 0.42), size = 18, bold = TRUE, color = cols[k], label = sprintf("step_%s", s_))
    deck_text(tx(sprintf("S08.steps.%s.claim", s_), f), c(xi, y + 0.52, wi, 1.95), size = 16, label = sprintf("claim_%s", s_))
    deck_text(tx(sprintf("S08.steps.%s.value", s_), f), c(xi, y + 2.5, wi, 0.52), size = 22, bold = TRUE, color = cols[k], label = sprintf("value_%s", s_))
    deck_text(tx(sprintf("S08.steps.%s.caption", s_), f), c(xi, y + 3.08, wi, 1.2), size = 16, color = PAL$ink2, label = sprintf("caption_%s", s_))
    deck_text(tx(sprintf("S08.steps.%s.ref", s_)), c(xi, y + bh - 0.47, wi, 0.4), size = 16, color = PAL$muted, label = sprintf("ref_%s", s_))
    if (k < length(steps)) s08_panel(c(x + bw + ag, y + bh / 2 - 0.25, aw, 0.5), PAL$muted, "rightArrow", sprintf("arrow_%d", k))
  }

  deck_notes(tx("S08.notes", c(f, list(
    cst = drange(CS, CB, "c_start1_median", 2, "", "cliff start concentration (mg/L), median, two models"),
    len595 = dspan(CS, CB, "len1_p05", "len1_p95", 2, "", "cliff length, 5th to 95th percentile, two models"),
    gmr = drange(TCH, "set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf ratio failing to retained, set (i)"),
    aii = dcount(CGf, "analysis_model=='M1' & config=='G2_A_ii' & pass_pct > 5", "M1 G2_A_ii cells above 5% (point)"),
    b = dcount(CGf, "analysis_model=='M1' & config=='G2_B' & pass_pct > 5", "M1 G2_B cells above 5% (point)"),
    inst = s08_inst("M1", "median"), inst0 = s08_inst("M0", "median"),
    p2ci = dci(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M1 P2 maximum with 95% CI"),
    m0c = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='conservative'", "M0 P2 cells classified conservative"),
    m0n = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='nominal'", "M0 P2 cells classified nominal")))))
  deck_end()
}
