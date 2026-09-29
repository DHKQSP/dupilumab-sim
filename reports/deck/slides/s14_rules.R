# S14 논거 ② 처리 규칙에 따라 판정이 달라진다. AUC0-inf + Cmax의 경계 1종 오류를 규칙 A/B/C x 기준 세트 (i)~(iv)(9가지 변형)로,
# AUC0-inf GMR 편향(참 AUC0-inf 비 대비), 판정 불안정(같은 시험의 판정이 변형 사이에서 갈리는 비율). M1(체중 층 포함), M0 병기.
# 자료: results/criteria/criteria_g2_type1.csv, criteria_bias.csv, criteria_instability.csv(시험 모집단 재생성 시험: 건강인, 체중 층화).
# 표 6행: 규칙 A 세트별 4행, 규칙 B 1행, 규칙 C는 SAP 제안의 세트 (i) 1행. 규칙 C 세트 (ii)~(iv)는 요점(최대)과 노트(세트별 모든 열)에 세트별로 적는다
# (세트 사이 최소·최대를 한 값으로 합치지 않는다).
# 칸 수 표기: "점추정 > 5%"(pass_pct > 5)와 "Wilson 하한 > 5%"(lo > 5, 분류 exceeding)를 따로 적는다(두 수가 다른 변형이 있다).

# 판정 불안정: 경계 16칸(S00, F097 제외)의 중앙값 또는 최대(DECISIONS D-058). 보고서는 M0만 인쇄하므로 M1은 같은 파일에서 읽는다.
s14_inst <- function(am, col, what = c("median", "max"), drop_ka = FALSE) {
  what <- match.arg(what); CSf <- "criteria/criteria_instability.csv"
  w <- sprintf("analysis_model=='%s' & !scenario %%in%% c('S00','F097')%s", am, if (drop_ka) " & !grepl('^ka_', scenario)" else "")
  r <- rows(CSf, w); premise(nrow(r) == if (drop_ka) 14 else 16, sprintf("boundary cells in the instability file (%s)", am))
  x <- if (what == "median") median(r[[col]]) else max(r[[col]])
  loc_ <- sprintf("analysis_model=='%s' & !scenario in (S00, F097)%s :: %s of %s", am, if (drop_ka) " & not ka" else "", what, col)
  dderived(sprintf("decision instability %s, %s over %d boundary cells, %s", col, what, nrow(r), am), CSf, loc_, x, paste0(fnum(x, 1), "%"))
}
# 여러 변형(config) 각각의 행 수를 세어 범위로(예: 규칙 C 네 세트의 5% 초과 칸 "1", "0~1")
s14_nrng <- function(rel, am, cfs, cond, item) {
  k <- vapply(cfs, function(cf) nrow(rows(rel, sprintf("analysis_model=='%s' & config=='%s' & %s", am, cf, cond))), 1L)
  dderived(item, rel, sprintf("analysis_model=='%s' & config in (%s) & %s :: range over configs of row counts", am, paste(cfs, collapse = ", "), cond), k, rng_fmt(min(k), max(k), 0))
}
s14_in <- function(v) sprintf("c(%s)", paste(sprintf("'%s'", v), collapse = ","))

slide_S14 <- function() {
  CG <- "criteria/criteria_g2_type1.csv"; CB <- "criteria/criteria_bias.csv"; CS <- "criteria/criteria_instability.csv"
  T1 <- "oc_models/type1_models.csv"; PW <- "criteria/criteria_g2_power.csv"
  AB <- c(A_i = "G2_A_i", A_ii = "G2_A_ii", A_iii = "G2_A_iii", A_iv = "G2_A_iv", B = "G2_B")
  CC <- c(C_i = "G2_C_i", C_ii = "G2_C_ii", C_iii = "G2_C_iii", C_iv = "G2_C_iv")
  bcf <- function(cf) sub("^G2_", "", cf)          # criteria_bias.csv의 config 표기(A_i, B, C_ii ...)
  deck_slide("S14", tag = "sim")

  # ---- 전제: 문장이 기대는 사실 ----
  for (am in c("M0", "M1")) for (cf in c(AB, CC)) {
    premise(nrow(rows(CG, sprintf("analysis_model=='%s' & config=='%s'", am, cf))) == 16, sprintf("16 boundary cells, %s %s", am, cf))
    premise(nrow(rows(CB, sprintf("analysis_model=='%s' & config=='%s'", am, bcf(cf)))) == 16, sprintf("16 bias rows, %s %s", am, cf))
  }
  p2 <- rows(T1, "analysis_model=='M1' & config=='P2'"); p2m <- p2[which.max(pass_pct)]
  premise(p2m$pk_model == "k2020" && p2m$scenario == "V2_up_080", "AUC0-last + Cmax is highest in the 2020 model V2 up cell under M1 (text)")
  for (cf in CC) { r <- rows(CG, sprintf("analysis_model=='M1' & config=='%s' & pass_pct > 5", cf))
    premise(nrow(r) == 1 && r$pk_model == p2m$pk_model && r$scenario == p2m$scenario, sprintf("%s: the only M1 cell above 5%% is the cell where P2 is highest (text)", cf)) }
  cnt <- function(rel, am, cf, cond) nrow(rows(rel, sprintf("analysis_model=='%s' & config=='%s' & %s", am, cf, cond)))
  premise(min(vapply(AB, function(cf) cnt(CG, "M1", cf, "pass_pct > 5"), 1L)) > max(vapply(CC, function(cf) cnt(CG, "M1", cf, "pass_pct > 5"), 1L)),
          "every rule A and B variant has more cells above 5% than any rule C variant under M1 (text: only rule C is close)")
  premise(min(vapply(AB, function(cf) cnt(CB, "M1", bcf(cf), "bias_dir=='toward_1'"), 1L)) > 8, "rules A and B: AUC0-inf GMR biased toward 1 in most of the 16 cells under M1 (text)")
  for (cf in CC) { r <- rows(CB, sprintf("analysis_model=='M1' & config=='%s' & bias_dir=='toward_1'", bcf(cf)))
    premise(nrow(r) == 1 && r$pk_model == "k2020" && r$scenario == "V2_up_080", sprintf("%s: the only cell biased toward 1 is the 2020 V2 up cell (notes)", cf)) }
  premise(all(rows(CG, "grepl('^ka_', scenario)")$pass_pct == 0) && all(rows(CS, "grepl('^ka_', scenario)")$inst_all == 0), "ka-down cells: every variant fails every trial, instability 0 (text)")
  premise(all(rows(T1, "grepl('^ka_', scenario)")$cmax_ratio < min(unlist(.read("config/trial_design.yaml")$be$limits))), "ka-down cells: true Cmax ratio below the lower equivalence limit (text: Cmax fails)")
  im <- rows(CS, "analysis_model=='M1' & !scenario %in% c('S00','F097')"); icols <- c("inst_all", "inst_rules_i", "inst_rules_ii", "inst_rules_iii", "inst_rules_iv", "inst_sets_A")
  premise(all(median(im$inst_sets_C) < vapply(icols, function(k) median(im[[k]]), 1)), "rule C across criteria sets has the smallest median instability under M1 (text)")
  cx_v2 <- vapply(CC, function(cf) rows(CG, sprintf("analysis_model=='M1' & config=='%s' & pk_model=='k2020' & scenario=='V2_up_080'", cf))$pass_pct, 1)
  premise(all(cx_v2[c("C_iii", "C_iv")] > p2m$pass_pct) && all(cx_v2[c("C_i", "C_ii")] <= p2m$pass_pct), "rules C (iii) and C (iv), and only these, are above AUC0-last + Cmax in the V2 cell under M1 (notes)")
  cls_v2 <- vapply(CC, function(cf) rows(CG, sprintf("analysis_model=='M1' & config=='%s' & pk_model=='k2020' & scenario=='V2_up_080'", cf))$class, "")
  premise(identical(unname(cls_v2), c("nominal", "exceeding", "exceeding", "exceeding")), "rule C in the V2 cell under M1: (i) nominal, (ii) to (iv) exceeding (text, notes)")
  for (cf in CC) premise(rows(CG, sprintf("analysis_model=='M1' & config=='%s'", cf))[which.max(pass_pct)][, pk_model == "k2020" & scenario == "V2_up_080"],
                         sprintf("%s: the largest M1 pass rate is in the 2020 V2 up cell (text: this cell)", cf))
  one_ <- dderived("reference ratio of the bias direction toward 1 (identical products)", CB, "bias_dir == 'toward_1' :: GMR moves toward the ratio 1", 1, "1")

  # ---- 제목 ----
  f <- list(nom = f_nominal(), ncell = dcount(CG, "analysis_model=='M1' & config=='G2_A_i'", "boundary cells per variant (M1)"),
            ab = s14_nrng(CG, "M1", AB, "pass_pct > 5", "M1 rules A (i) to (iv) and B: cells above 5% (point), range over variants"),
            c = s14_nrng(CG, "M1", CC, "pass_pct > 5", "M1 rule C (i) to (iv): cells above 5% (point), same count in every set"))
  premise(length(unique(vapply(CC, function(cf) cnt(CG, "M1", cf, "pass_pct > 5"), 1L))) == 1, "rule C: every set has the same number of cells above 5% under M1 (title: in every set)")
  deck_kicker(tx("S14.kicker")); deck_title(tx("S14.title", f))

  # ---- 표: 규칙 A (i)~(iv), B, C (i)(SAP 제안 세트) ----
  T <- DK$txt$S14$table
  one <- function(k, cf) {
    w1 <- sprintf("analysis_model=='M1' & config=='%s'", cf); w0 <- sprintf("analysis_model=='M0' & config=='%s'", cf); wb <- sprintf("analysis_model=='M1' & config=='%s'", bcf(cf))
    c(T$rows[[k]],
      sprintf("%s (%s)", dcount(CG, paste(w1, "& pass_pct > 5"), sprintf("M1 %s cells above 5%% (point)", cf)), dcount(CG, paste(w0, "& pass_pct > 5"), sprintf("M0 %s cells above 5%% (point)", cf))),
      sprintf("%s (%s)", dcount(CG, paste(w1, "& lo > 5"), sprintf("M1 %s cells with Wilson lower bound above 5%%", cf)), dcount(CG, paste(w0, "& lo > 5"), sprintf("M0 %s cells with Wilson lower bound above 5%%", cf))),
      dext(CG, w1, "pass_pct", max, 2, "%", sprintf("M1 %s largest boundary pass rate", cf)),
      drange(CB, wb, "bias_pct", 2, "%", sprintf("M1 AUC0-inf GMR bias against the true ratio, %s", bcf(cf))),
      dcount(CB, paste(wb, "& bias_dir=='toward_1'"), sprintf("M1 %s cells with AUC0-inf GMR biased toward 1", bcf(cf))))
  }
  RW <- c(AB, C_i = CC[["C_i"]])
  m <- do.call(rbind, lapply(names(RW), function(k) one(k, RW[[k]])))
  df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- tx("S14.table.head", f)
  LW <- 7.85; th <- 2.52
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, LW, th), widths = c(1.84, 1.21, 1.48, 0.87, 1.33, 1.12), size = 12, highlight = 6)

  # ---- 요점(표 아래) ----
  abmax <- dext(CG, sprintf("analysis_model=='M1' & config %%in%% %s", s14_in(AB)), "pass_pct", max, 2, "%", "M1 rules A and B: largest boundary pass rate")
  p2v <- dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 AUC0-last + Cmax, 2020 model V2 up cell")
  WV2 <- "analysis_model=='M1' & pk_model=='k2020' & scenario=='V2_up_080'"
  civ1 <- dv(CG, sprintf("%s & config=='G2_C_i'", WV2), "pass_pct", 2, "%", "M1 G2_C_i, 2020 model V2 up cell")
  cxr <- drange(CG, sprintf("%s & config %%in%% %s", WV2, s14_in(CC[c("C_ii", "C_iii", "C_iv")])), "pass_pct", 2, "%", "M1 G2_C_ii to G2_C_iv, 2020 model V2 up cell, range over sets (ii) to (iv)")
  yb <- GEO$BODY_TOP + th + 0.14
  deck_bullets(tx("S14.bullets", list(c = f$c, nom = f$nom, p2 = p2v, abmax = abmax, ci = civ1, cx = cxr, one = one_)), box = c(GEO$ML, yb, LW, GEO$BODY_BOTTOM - yb), size = 16, gap_pt = 6)

  # ---- 오른쪽: 판정 불안정(정의, 카드 두 개) ----
  xr <- GEO$ML + LW + 0.3; wr <- GEO$W - GEO$MR - xr
  inf <- list(ncell = f$ncell, nka = dcount(CS, "analysis_model=='M1' & grepl('^ka_', scenario)", "ka-down boundary cells (M1)"),
              zero = dext(CS, "analysis_model=='M1' & grepl('^ka_', scenario)", "inst_all", max, 0, "%", "decision instability in the ka-down cells, M1"),
              ninf = dcount(CS, "analysis_model=='M1' & !scenario %in% c('S00','F097') & !grepl('^ka_', scenario)", "informative boundary cells (M1)"))
  nvar <- local({ k <- sort(unique(rows(CG, "analysis_model=='M1'")$config)); premise(setequal(k, c(AB, CC)), "nine variants in the criteria file")
    dderived("number of AUC0-inf + Cmax variants (rule x criteria set)", CG, "analysis_model=='M1' :: count of distinct config", length(k), as.character(length(k))) })
  inf$nvar <- nvar
  dh <- 1.42
  deck_text(tx("S14.inst.def", inf), c(xr, GEO$BODY_TOP, wr, dh), size = 16, label = "inst_def", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 2)
  ch <- (GEO$BODY_BOTTOM - GEO$BODY_TOP - dh - 0.24) / 2; y1 <- GEO$BODY_TOP + dh + 0.12
  i_all <- s14_inst("M1", "inst_all", "median"); i14 <- s14_inst("M1", "inst_all", "median", drop_ka = TRUE)
  deck_stat(i_all, tx("S14.inst.card1", list(nvar = nvar, ncell = f$ncell, ninf = inf$ninf, i14 = i14, max = s14_inst("M1", "inst_all", "max"),
                                             m0 = s14_inst("M0", "inst_all", "median"), m0max = s14_inst("M0", "inst_all", "max"))),
            c(xr, y1, wr, ch), color = PAL$orange, bg = PAL$tint_orange, value_size = 36)
  deck_stat(s14_inst("M1", "inst_sets_C", "median"), tx("S14.inst.card2", list(ncell = f$ncell, max = s14_inst("M1", "inst_sets_C", "max"), a = s14_inst("M1", "inst_sets_A", "median"), amax = s14_inst("M1", "inst_sets_A", "max"))),
            c(xr, y1 + ch + 0.12, wr, ch), value_size = 36)

  # ---- 노트 ----
  ci1 <- function(cf, am = "M1", pk = "k2020", sc = "V2_up_080") dci(CG, sprintf("analysis_model=='%s' & config=='%s' & pk_model=='%s' & scenario=='%s'", am, cf, pk, sc), "pass_pct", "lo", "hi", 2, "%", sprintf("%s %s, %s %s", am, cf, pk, sc))
  cmaxcell <- function(am, cf) { r <- rows(CG, sprintf("analysis_model=='%s' & config=='%s'", am, cf))[which.max(pass_pct)]; ci1(cf, am, r$pk_model, r$scenario) }
  premise(rows(CG, "analysis_model=='M1' & config=='G2_A_ii'")[which.max(pass_pct)]$scenario == "Vmax_down_125" && rows(CG, "analysis_model=='M1' & config=='G2_A_ii'")[which.max(pass_pct)]$pk_model == "k2020", "M1 rule A (ii) highest in the 2020 Vmax down cell (notes)")
  premise(rows(CG, "analysis_model=='M1' & config=='G2_B'")[which.max(pass_pct)]$scenario == "Vmax_up_080" && rows(CG, "analysis_model=='M1' & config=='G2_B'")[which.max(pass_pct)]$pk_model == "k2016", "M1 rule B highest in the 2016 Vmax up cell (notes)")
  nk <- "!grepl('^ka_', scenario)"
  cn1 <- function(am, cf, cond, what) dcount(CG, sprintf("analysis_model=='%s' & config=='%s' & %s", am, cf, cond), sprintf("%s %s cells %s", am, cf, what))
  s00 <- rows(CS, "analysis_model=='M1' & scenario=='S00'"); premise(nrow(s00) == 2, "identical-product rows, two models")
  cset <- function(k, am, cond, what) cn1(am, CC[[k]], cond, what)
  cpair <- function(k, cond, what) sprintf("%s (%s)", cset(k, "M1", cond, what), cset(k, "M0", cond, what))
  cb <- function(k) drange(CB, sprintf("analysis_model=='M1' & config=='%s'", k), "bias_pct", 2, "%", sprintf("M1 AUC0-inf GMR bias against the true ratio, %s", k))
  ct <- function(k) dcount(CB, sprintf("analysis_model=='M1' & config=='%s' & bias_dir=='toward_1'", k), sprintf("M1 %s cells with AUC0-inf GMR biased toward 1", k))
  c0m <- function(k) dext(CG, sprintf("analysis_model=='M0' & config=='%s'", CC[[k]]), "pass_pct", max, 2, "%", sprintf("M0 %s largest boundary pass rate", CC[[k]]))
  cs_ <- list(); for (q in seq_along(CC)) { k <- names(CC)[q]
    cs_[[sprintf("c%d_pt", q)]] <- cpair(k, "pass_pct > 5", "above 5% (point)"); cs_[[sprintf("c%d_lo", q)]] <- cpair(k, "lo > 5", "with Wilson lower bound above 5%")
    cs_[[sprintf("c%d_b", q)]] <- cb(k); cs_[[sprintf("c%d_t", q)]] <- ct(k); cs_[[sprintf("c0_%d", q)]] <- c0m(k) }
  deck_notes(tx("S14.notes", c(cs_, list(nvar = nvar, one = one_, zero = inf$zero, nka = inf$nka, ncell = f$ncell,
    wt = f_wt_range(), n_arm = f_n_arm(), reps = f_reps("boundary"), nom = f$nom,
    n20 = dint(CG, "analysis_model=='M1' & config=='G2_C_ii' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials, 2020 model V2 up cell"),
    ai_pt = cn1("M1", "G2_A_i", "pass_pct > 5", "above 5% (point)"), ai_lo = cn1("M1", "G2_A_i", "lo > 5", "with Wilson lower bound above 5%"),
    b_pt = cn1("M1", "G2_B", "pass_pct > 5", "above 5% (point)"), b_lo = cn1("M1", "G2_B", "lo > 5", "with Wilson lower bound above 5%"),
    ci_pt = cn1("M1", "G2_C_i", "pass_pct > 5", "above 5% (point)"), ci_lo = cn1("M1", "G2_C_i", "lo > 5", "with Wilson lower bound above 5%"),
    ci = ci1("G2_C_i"), cii = ci1("G2_C_ii"), ciii = ci1("G2_C_iii"), civ = ci1("G2_C_iv"),
    p2 = dci(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M1 AUC0-last + Cmax, 2020 model V2 up cell"),
    aii = cmaxcell("M1", "G2_A_ii"), b = cmaxcell("M1", "G2_B"),
    aii0 = dext(CG, "analysis_model=='M0' & config=='G2_A_ii'", "pass_pct", max, 2, "%", "M0 G2_A_ii largest boundary pass rate"),
    b0 = dext(CG, "analysis_model=='M0' & config=='G2_B'", "pass_pct", max, 2, "%", "M0 G2_B largest boundary pass rate"),
    tgt = dcfg("oc_design.yaml", "boundary_targets", "boundary true AUC0-inf ratios", function(x) paste(fnum(x, 2), collapse = ", ")),
    t1 = dcfg("oc_design.yaml", "boundary_targets", "lower boundary true AUC0-inf ratio", function(x) fnum(min(x), 2)),
    t2 = dcfg("oc_design.yaml", "boundary_targets", "upper boundary true AUC0-inf ratio", function(x) fnum(max(x), 2)),
    aii_b14 = drange(CB, sprintf("analysis_model=='M1' & config=='A_ii' & %s", nk), "bias_pct", 2, "%", "M1 AUC0-inf GMR bias, A_ii, without ka cells"),
    b_b14 = drange(CB, sprintf("analysis_model=='M1' & config=='B' & %s", nk), "bias_pct", 2, "%", "M1 AUC0-inf GMR bias, B, without ka cells"),
    aii_b0 = drange(CB, "analysis_model=='M0' & config=='A_ii'", "bias_pct", 2, "%", "M0 AUC0-inf GMR bias, A_ii"),
    c_b0 = drange(CB, "analysis_model=='M0' & config=='C_i'", "bias_pct", 2, "%", "M0 AUC0-inf GMR bias, C_i"),
    r_i = s14_inst("M1", "inst_rules_i", "median"), r_ii = s14_inst("M1", "inst_rules_ii", "median"),
    r_iii = s14_inst("M1", "inst_rules_iii", "median"), r_iv = s14_inst("M1", "inst_rules_iv", "median"),
    i14 = i14, ninf = inf$ninf,
    s00 = drange(CS, "analysis_model=='M1' & scenario=='S00'", "inst_all", 1, "%", "decision instability, identical products, M1, two models"),
    s00a = drange(CS, "analysis_model=='M1' & scenario=='S00'", "inst_sets_A", 1, "%", "instability across sets within rule A, identical products, M1"),
    s00r = drange(CS, "analysis_model=='M1' & scenario=='S00'", "inst_rules_i", 1, "%", "instability across rules with set (i), identical products, M1"),
    pw4 = drange(PW, "scenario=='S00' & analysis_model=='M1' & config=='G2_A_iv'", "pass_pct", 1, "%", "power, identical products, M1, G2_A_iv, two models")))))
  deck_end()
}
