# S02 핵심 결론: 제안(AUC0-last + Cmax 공동 1차, AUC0-inf 이차), 한 문장 요약, 논거 ①②③과 각 대표 수치(headline). 시험 모집단만, M1 주분석(M0 병기).
# 첫 내용 슬라이드이므로 여기서 처음 나오는 약어(NCA, SAP, λz)를 풀어 쓴다.
s02_inst_m1 <- function(col = "inst_all", what = c("median", "max")) {
  what <- match.arg(what); CSf <- "criteria/criteria_instability.csv"
  r <- rows(CSf, "analysis_model=='M1' & !scenario %in% c('S00','F097')"); premise(nrow(r) == 16, "16 boundary cells in the instability file (M1)")
  x <- if (what == "median") median(r[[col]]) else max(r[[col]])
  dderived(sprintf("decision instability %s, %s over 16 boundary cells, M1", col, what), CSf,
           sprintf("analysis_model=='M1' & !scenario in (S00, F097) :: %s of %s", what, col), x, paste0(fnum(x, 1), "%"))
}
s02_inst_m0 <- function(col = "inst_all", what = c("median", "max")) {
  what <- match.arg(what); CSf <- "criteria/criteria_instability.csv"
  r <- rows(CSf, "analysis_model=='M0' & !scenario %in% c('S00','F097')"); premise(nrow(r) == 16, "16 boundary cells in the instability file (M0)")
  x <- if (what == "median") median(r[[col]]) else max(r[[col]])
  dderived(sprintf("decision instability %s, %s over 16 boundary cells, M0", col, what), CSf,
           sprintf("analysis_model=='M0' & !scenario in (S00, F097) :: %s of %s", what, col), x, paste0(fnum(x, 1), "%"))
}

# 글자 없는 도형(카드 배경): deck_box는 8pt 자리 글자를 넣어 검사 6(글자 크기)에 걸리므로 16pt 빈 글상자로 그린다
s02_panel <- function(box, fill, geom = "roundRect", label = "panel") deck_text(" ", box, size = 16, bg = fill, geom = geom, label = label)

# 창 포착률 최솟값: "이상"이라고 쓰므로 0.1 단위로 내림(보고서 tpcv_min과 같은 파일·조건·계산)
s02_cov_min <- function() {
  TCV <- "trialpop/tp_coverage_individual.csv"; w <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  r <- rows(TCV, w); premise(nrow(r) == 2, "window coverage: two models"); x <- min(r$min) * 100
  dderived("window coverage (true AUC0-tlast / true AUC0-inf), smallest subject-level value over both models", TCV, sprintf("%s :: min(min)", w), x, sprintf("%s%%", fnum(floor(x * 10) / 10, 1)))
}

slide_S02 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; TCH <- "trialpop/tp_characteristics.csv"
  CGf <- "criteria/criteria_g2_type1.csv"; T1f <- "oc_models/type1_models.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  deck_slide("S02", tag = "sim")
  deck_kicker(tx("S02.kicker")); deck_title(tx("S02.title"))

  # 전제 검사: 문장이 기대는 방향·순서
  premise(all(rows(TCH, "set=='i'")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUC0-inf than retained subjects (text: not random)")
  a_i <- rows(CGf, "analysis_model=='M1' & config=='G2_A_i'"); c_i <- rows(CGf, "analysis_model=='M1' & config=='G2_C_i'")
  premise(sum(a_i$pass_pct > 5) > sum(c_i$pass_pct > 5), "rule A (i) has more boundary cells above 5% than rule C (i) under M1 (text: decision depends on the rule)")
  premise(sum(a_i$lo > 5) == sum(a_i$class == "exceeding") && sum(c_i$lo > 5) == sum(c_i$class == "exceeding"), "Wilson lower bound above 5% equals class 'exceeding' in the criteria file (text: Wilson counts)")
  c_pt <- c_i[pass_pct > 5]; premise(nrow(c_pt) == 1 && c_pt$class == "nominal", "the single rule C (i) cell above 5% (point) is Wilson nominal under M1 (text)")
  p2 <- rows(T1f, "analysis_model=='M1' & config=='P2'"); premise(nrow(p2) == 16 && sum(p2$class == "exceeding") == 1, "one exceeding P2 cell under M1 (text)")
  m0n <- rows(T1f, "analysis_model=='M0' & config=='P2' & class=='nominal'")
  premise(nrow(m0n) == 1 && m0n$pk_model == "k2020" && m0n$scenario == "V2_up_080", "the single nominal M0 P2 cell is the same 2020 model V2 up cell (notes: same cell)")
  premise(p2[class == "exceeding"]$pk_model == "k2020", "the exceeding P2 cell under M1 is from the 2020 model (text)")
  pmax <- p2[which.max(pass_pct)]; premise(pmax$pk_model == "k2020" && pmax$scenario == "V2_up_080", "the exceeding P2 cell under M1 is the 2020 model V2 up (notes)")

  f <- list(
    wt = f_wt_range(), nom = f_nominal(),
    # 논거 ① 대표 수치: 세트 (iii) 미달(두 모델)
    iii = headline(drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models")),
    i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
    r2iii = f_set("iii", "r2"), r2i = f_set("i", "r2"),
    # 논거 ② 대표 수치: M1 규칙 A (i)에서 5% 초과 칸
    a = headline(dcount(CGf, "analysis_model=='M1' & config=='G2_A_i' & pass_pct > 5", "M1 G2_A_i cells above 5% (point)")),
    ncell = dcount(CGf, "analysis_model=='M1' & config=='G2_A_i'", "boundary cells per configuration (M1)"),
    alo = dcount(CGf, "analysis_model=='M1' & config=='G2_A_i' & lo > 5", "M1 G2_A_i cells with Wilson lower bound above 5%"),
    c = dcount(CGf, "analysis_model=='M1' & config=='G2_C_i' & pass_pct > 5", "M1 G2_C_i cells above 5% (point)"),
    gmr2 = drange(TCH, "set=='i'", "true_aucinf_gmr", 2, "", "true AUC0-inf ratio failing to retained, set (i), two models"),
    # 논거 ③ 대표 수치: 창 포착률 중앙값(두 모델)
    cov = headline(drange(TCV, MW, "median", 1, "%", "window coverage, median, two models", scale = 100)),
    cmin = s02_cov_min(),
    n_all = dcount(T1f, "analysis_model=='M1' & config=='P2'", "M1 P2 boundary cells"),
    n_cons = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='conservative'", "M1 P2 cells classified conservative"),
    n_exc = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='exceeding'", "M1 P2 cells classified exceeding"),
    p2max = dv(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2 maximum"))

  # 위: 제안 카드 + 한 문장 요약
  top <- GEO$BODY_TOP; h_top <- 1.52; pw <- 5.6
  deck_text(tx("S02.proposal"), c(GEO$ML, top, pw, h_top), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "proposal", gap_pt = 4)
  deck_text(tx("S02.summary", f), c(GEO$ML + pw + 0.25, top, GEO$CW - pw - 0.25, h_top), size = 17, label = "summary")

  # 아래: 논거 ①②③ 카드 세 개(위 대표 수치, 가운데 주장, 아래 근거)
  y0 <- top + h_top + 0.12; gw <- 0.22; cw <- (GEO$CW - 2 * gw) / 3; ch <- GEO$BODY_BOTTOM - y0; sh <- 0.96
  args <- lapply(c("a1", "a2", "a3"), function(k) list(lab = tx(sprintf("S02.args.%s.label", k), f), claim = tx(sprintf("S02.args.%s.claim", k)), det = tx(sprintf("S02.args.%s.detail", k), f)))
  vals <- c(f$iii, tx("S02.args.a2.value", f), f$cov)
  for (k in seq_along(args)) {
    a <- args[[k]]; x <- GEO$ML + (k - 1) * (cw + gw)
    s02_panel(c(x, y0, cw, ch), PAL$tint_grey, "roundRect", sprintf("card_%d", k))
    deck_stat(vals[k], a$lab, c(x, y0, cw, sh), value_size = 28)
    deck_text(a$claim, c(x + 0.04, y0 + sh + 0.02, cw - 0.08, 0.78), size = 18, bold = TRUE, label = sprintf("claim_%d", k))
    deck_text(a$det, c(x + 0.04, y0 + sh + 0.80, cw - 0.08, ch - sh - 0.80), size = 16, color = PAL$ink2, label = sprintf("detail_%d", k))
  }

  deck_notes(tx("S02.notes", c(f, list(
    n_ind = dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"), dose = f_dose(),
    a_max = dext(CGf, "analysis_model=='M1' & config=='G2_A_i'", "pass_pct", max, 2, "%", "M1 G2_A_i largest boundary pass rate"),
    aii = dcount(CGf, "analysis_model=='M1' & config=='G2_A_ii' & pass_pct > 5", "M1 G2_A_ii cells above 5% (point)"),
    clo = dcount(CGf, "analysis_model=='M1' & config=='G2_C_i' & lo > 5", "M1 G2_C_i cells with Wilson lower bound above 5%"),
    c_ci = dci(CGf, "analysis_model=='M1' & config=='G2_C_i' & pass_pct > 5", "pass_pct", "lo", "hi", 2, "%", "M1 G2_C_i cell above 5% (point) with Wilson 95% CI"),
    b = dcount(CGf, "analysis_model=='M1' & config=='G2_B' & pass_pct > 5", "M1 G2_B cells above 5% (point)"),
    inst = s02_inst_m1("inst_all", "median"), inst_max = s02_inst_m1("inst_all", "max"), inst0 = s02_inst_m0("inst_all", "median"), inst0_max = s02_inst_m0("inst_all", "max"),
    p2ci = dci(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M1 P2 maximum with 95% CI"),
    tgt = dv(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "target", 2, "", "true AUC0-inf ratio target, exceeding cell"),
    ntr = dint(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials, exceeding cell"),
    m0c = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='conservative'", "M0 P2 cells classified conservative"),
    m0n = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='nominal'", "M0 P2 cells classified nominal"),
    m0max = dv(T1f, "analysis_model=='M0' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M0 P2 maximum"),
    gmr = drange(TCH, "set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf ratio failing to retained, set (i)")))))
  deck_end()
}
