# S02 핵심 결론: 제안(AUC0-last + Cmax 공동 1차, AUC0-inf 2차), 한 문장 요약, 논거 ①②③과 각 대표 수치(headline). 시험 모집단만, M1 주분석(M0 함께 제시).
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

# 큰 수치 카드(deck_stat과 같은 모양·이름 "stat": 세로 가운데 정렬)에 수치 뒤 작은 한정어(판정 구성)를 같은 줄에 붙인다.
# 설명 줄이 한 줄에 들어가게 하려고 한정어를 수치 줄로 옮긴다(카드 세 개의 수치 띠 높이를 같게 유지).
s02_stat <- function(value, suffix, label, box, value_size = 28) {
  fit_check("stat_label", label, c(box[1], box[2], box[3], box[4] - value_size * 1.25 / 72), SZ$stat_label, gap_pt = 0, card = box[4] >= CARD_MIN_H)
  wv <- text_w(strip_markup(value), value_size, TRUE) + text_w(paste0("  ", suffix), SZ$stat_label)
  if (wv > box[3] - 2 * CARD_INS[["lr"]]) stop(sprintf("S02 stat value line %.2f in does not fit %.2f in", wv, box[3] - 2 * CARD_INS[["lr"]]), call. = FALSE)
  v <- c(runs(nobreak(value), value_size, PAL$blue, bold = TRUE), list(ftext(nobreak(paste0("  ", suffix)), ftp(SZ$stat_label, PAL$ink2))))
  p1 <- do.call(fpar, c(v, list(fp_p = fp_par(text.align = "left", padding.bottom = 2, line_spacing = 1.0))))
  p2 <- para(label, SZ$stat_label, PAL$ink, FALSE, "left", gap_pt = 0)
  DK$x <- ph_with(DK$x, block_list(p1, p2), location = loc(box, "stat", bg = PAL$tint_blue, geom = "roundRect", ln = no_line())); invisible(NULL)
}

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
  # 제목은 한 줄: 제목 상자를 한 줄 높이로 줄이고 본문을 그만큼 올린다(두 줄이 되면 fit_check가 빌드를 멈춘다)
  deck_kicker(tx("S02.kicker")); deck_title(tx("S02.title"), box = c(GEO$ML, GEO$TITLE_TOP, GEO$CW - 0.1, 0.68))

  # 전제 검사: 문장이 기대는 방향·순서
  premise(all(rows(TCH, "set=='i'")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUC0-inf than retained subjects (text: not random)")
  fs <- rows(TPF, "TRUE"); premise(all(fs[, .(ok = fail_pct[set == "i"] == min(fail_pct)), by = pk_model]$ok), "set (i) has the smallest failing share in each model (text: most lenient set)")
  a_i <- rows(CGf, "analysis_model=='M1' & config=='G2_A_i'"); c_i <- rows(CGf, "analysis_model=='M1' & config=='G2_C_i'")
  premise(sum(a_i$pass_pct > 5) > sum(c_i$pass_pct > 5), "rule A (i) has more boundary cells above 5% than rule C (i) under M1 (text: decision depends on the rule)")
  premise(sum(a_i$lo > 5) == sum(a_i$class == "exceeding") && sum(c_i$lo > 5) == sum(c_i$class == "exceeding"), "Wilson lower bound above 5% equals class 'exceeding' in the criteria file (text: Wilson counts)")
  c_pt <- c_i[pass_pct > 5]; premise(nrow(c_pt) == 1 && c_pt$class == "nominal", "the single rule C (i) cell above 5% (point) is Wilson nominal under M1 (text)")
  p2 <- rows(T1f, "analysis_model=='M1' & config=='P2'"); premise(nrow(p2) == 16 && sum(p2$class == "exceeding") == 1, "one exceeding P2 cell under M1 (text)")
  m0n <- rows(T1f, "analysis_model=='M0' & config=='P2' & class=='nominal'")
  premise(nrow(m0n) == 1 && m0n$pk_model == "k2020" && m0n$scenario == "V2_up_080", "the single nominal M0 P2 cell is the same 2020 model V2 up cell (notes: same cell)")
  premise(p2[class == "exceeding"]$pk_model == "k2020", "the exceeding P2 cell under M1 is from the 2020 model (text)")
  premise(sum(p2$class == "conservative") + sum(p2$class == "exceeding") == nrow(p2), "M1 P2 cells are conservative except the exceeding one (text: only that cell)")
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
    gmr2 = drange(TCH, "set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf ratio failing to retained, set (i), two models"),
    # 논거 ③ 대표 수치: 창 포착률 중앙값(두 모델)
    cov = headline(drange(TCV, MW, "median", 1, "%", "window coverage, median, two models", scale = 100)),
    cmin = s02_cov_min(),
    n_all = dcount(T1f, "analysis_model=='M1' & config=='P2'", "M1 P2 boundary cells"),
    n_cons = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='conservative'", "M1 P2 cells classified conservative"),
    n_exc = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='exceeding'", "M1 P2 cells classified exceeding"),
    p2max = dv(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2 maximum"))

  # 위: 제안 카드(왼쪽 '제안' 표지 + 내어 쓴 네 줄) + 한 문장 요약
  top <- 1.20; h_top <- 1.52; pw <- 5.8; lw <- 0.62; ty <- top + 0.035
  s02_panel(c(GEO$ML, top, pw, h_top), PAL$tint_orange, "roundRect", "proposal_panel")
  deck_text(tx("S02.proposal_label"), c(GEO$ML + 0.04, ty, lw, 0.42), size = 16, bold = TRUE, label = "proposal_label")
  deck_text(tx("S02.proposal"), c(GEO$ML + lw, ty, pw - lw - 0.04, h_top - 0.12), size = 16, label = "proposal", gap_pt = 3)
  deck_text(tx("S02.summary", f), c(GEO$ML + pw + 0.22, ty, GEO$CW - pw - 0.22, h_top - 0.12), size = 16, label = "summary")

  # 아래: 논거 ①②③ 카드 세 개(위 대표 수치, 가운데 주장, 아래 근거) + 맨 아래 용어 풀이 한 줄(창 포착률, Wilson 분류)
  hw <- 0.39; yw <- GEO$BODY_BOTTOM - hw
  y0 <- top + h_top + 0.24; gw <- 0.18; cw <- (GEO$CW - 2 * gw) / 3; ch <- yw - 0.08 - y0; sh <- 0.96
  args <- lapply(c("a1", "a2", "a3"), function(k) list(lab = tx(sprintf("S02.args.%s.label", k), f), claim = tx(sprintf("S02.args.%s.claim", k)), det = tx(sprintf("S02.args.%s.detail", k), f)))
  vals <- c(f$iii, tx("S02.args.a2.value", f), f$cov)
  for (k in seq_along(args)) {
    a <- args[[k]]; x <- GEO$ML + (k - 1) * (cw + gw)
    s02_panel(c(x, y0, cw, ch), PAL$tint_grey, "roundRect", sprintf("card_%d", k))
    if (k == 2) s02_stat(vals[k], tx("S02.args.a2.suffix"), a$lab, c(x, y0, cw, sh)) else deck_stat(vals[k], a$lab, c(x, y0, cw, sh), value_size = 28)
    deck_text(a$claim, c(x + 0.02, y0 + sh + 0.02, cw - 0.04, 0.78), size = 18, bold = TRUE, label = sprintf("claim_%d", k))
    deck_text(a$det, c(x + 0.02, y0 + sh + 0.80, cw - 0.04, ch - sh - 0.80), size = 16, color = PAL$ink2, label = sprintf("detail_%d", k), gap_pt = 4)
  }
  deck_text(tx("S02.gloss", f), c(GEO$ML, yw, GEO$CW, hw), size = 16, color = PAL$ink2, label = "text_gloss")

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
    m0max = dv(T1f, "analysis_model=='M0' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M0 P2 maximum")))))
  deck_end()
}
