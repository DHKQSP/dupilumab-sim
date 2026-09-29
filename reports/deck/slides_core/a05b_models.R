# A5b ③ 세부(별첨): 분석 모형 M0/M1/M2와 AUClast + Cmax 초과 한 칸의 분해.
# (a) 표: 판정 구성별 경계 16칸 중 점추정 > 5% 칸 수와 최댓값, M0(층 미포함)·M1(체중 층)·M2(log 체중). AUClast + Cmax(P2), AUCinf + Cmax 규칙 A 세트 (iii)
#     (criteria_g2_type1.csv G2_A_iii, 사전 등록상 M0·M1만), 참고로 규칙 A 세트 (i)·(ii), 규칙 B, 규칙 C 세트 (i)(type1_models.csv G2_Ai, G2_Aii, G2_B, G2_Ci).
# (b) 그림: 2020 모델 말초 분포 증가 칸(M1, 시험 20,000회)에서 참 AUCinf 단독 판정(편향 없는 기준) → AUClast 단독 → AUClast + Cmax의 경계 1종 오류와 Wilson 95% 구간
#     (type1_models.csv AUCinf_true_only, AUClast_only, P2), 단계별 차이는 같은 시험 안의 분해(p2_decomposition_models.csv; D-060).
# 캡션 첫 줄: AUCinf 쪽 기전(표적 소실 증가 칸의 비구획 AUCinf GMR 편향, criteria_bias.csv M1, 규칙 A (iii)과 λz 산출 전원). (b)의 효과 글자 줄은 흰 띠로 5% 점선을 가린다.
# 자료 논리: 결과보고 덱 s19_transparency.R, s15_mechanism.R.
a5b_signed <- function(p) if (grepl("^-", p)) p else paste0("+", p)

slide_A5b <- function() {
  T1 <- "oc_models/type1_models.csv"; CG <- "criteria/criteria_g2_type1.csv"; DE <- "oc_models/p2_decomposition_models.csv"
  W <- "pk_model=='k2020' & scenario=='V2_up_080'"
  WD <- function(am) sprintf("pk_model=='k2020' & mechanism=='V2' & analysis_model=='%s'", am)       # 결과보고 덱 S19와 같은 행 조건
  deck_slide("A5b", tag = "sim")
  L <- DK$txt$A5b
  nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  AM <- c("M0", "M1", "M2")

  # ---- 전제 ----
  for (am in AM) { r <- rows(T1, sprintf("analysis_model=='%s' & config=='P2' & pass_pct > 5", am))
    premise(nrow(r) == 1 && r$pk_model == "k2020" && r$scenario == "V2_up_080", sprintf("%s: the only AUClast + Cmax cell above 5%% is the 2020 V2 up cell", am)) }
  premise(nrow(rows(CG, "analysis_model=='M2'")) == 0, "criteria-set configurations were run under M0 and M1 only (M2 cell: not computed)")
  for (am in c("M0", "M1")) for (s_ in c("i", "ii")) {
    x <- rows(T1, sprintf("analysis_model=='%s' & config=='G2_A%s'", am, s_)); y <- rows(CG, sprintf("analysis_model=='%s' & config=='G2_A_%s'", am, s_))
    premise(nrow(x) == 16 && isTRUE(all.equal(x[order(pk_model, scenario), pass_pct], y[order(pk_model, scenario), pass_pct])), sprintf("G2_A%s in type1_models equals set (%s) of rule A in the criteria file (%s)", s_, s_, am)) }
  x <- rows(T1, "analysis_model=='M1' & config=='G2_Ci'"); y <- rows(CG, "analysis_model=='M1' & config=='G2_C_i'")
  premise(isTRUE(all.equal(x[order(pk_model, scenario), pass_pct], y[order(pk_model, scenario), pass_pct])), "G2_Ci equals rule C set (i)")
  tw <- rows(DE, "auclast_bias_dir=='toward_1'")
  premise(nrow(tw) == 3 && all(tw$pk_model == "k2020" & tw$scenario == "V2_up_080") && all(tw$auclast_bias_pct > 0), "AUClast GMR biased toward 1 only in this cell, under each analysis model (body)")
  premise(nrow(rows(DE, "analysis_model=='M1'")) == 16, "16 boundary cells in the decomposition file (M1)")
  dd <- rows(DE, W); premise(nrow(dd) == 3 && all(abs(dd$ref_minus_5_pp + dd$auclast_minus_ref_pp + dd$p2_minus_auclast_pp - dd$p2_minus_5_pp) < 1e-9), "decomposition adds up to P2 minus 5")
  d1 <- rows(T1, sprintf("%s & analysis_model=='M1' & config %%in%% c('AUCinf_true_only','AUClast_only','P2')", W)); de1 <- row1(DE, WD("M1"))
  premise(nrow(d1) == 3 && abs(d1[config == "AUCinf_true_only", pass_pct] - de1$ref_pct) < 1e-9 && abs(d1[config == "AUClast_only", pass_pct] - de1$auclast_pct) < 1e-9 &&
            abs(d1[config == "P2", pass_pct] - de1$p2_pct) < 1e-9, "figure values (type1_models) equal the decomposition file")
  premise(de1$ref_pct < nomv && de1$p2_pct > nomv && de1$auclast_minus_ref_pp > de1$p2_minus_5_pp && de1$p2_minus_auclast_pp < 0, "reference below 5%, the AUClast term carries the excess, Cmax lowers it slightly (title, body)")
  premise(row1(T1, sprintf("%s & analysis_model=='M1' & config=='P2'", W))$class == "exceeding" && row1(T1, sprintf("%s & analysis_model=='M0' & config=='P2'", W))$class == "nominal" &&
            row1(T1, sprintf("%s & analysis_model=='M2' & config=='P2'", W))$class == "exceeding", "M1 and M2 exceeding, M0 nominal (notes)")

  # ---- 제목 ----
  k_iii <- vapply(c("M0", "M1"), function(am) nrow(rows(CG, sprintf("analysis_model=='%s' & config=='G2_A_iii' & pass_pct > 5", am))), 1L)
  f <- list(nom = f_nominal(), n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"),
            k = dderived("AUCinf (set iii) + Cmax cells above 5% (point), range over M0 and M1", CG, "analysis_model in (M0, M1) & config=='G2_A_iii' & pass_pct > 5 :: range over models of row counts",
                         k_iii, rng_fmt(min(k_iii), max(k_iii), 0)))
  premise(all(k_iii > 0) && all(vapply(c("M0", "M1"), function(am) nrow(rows(CG, sprintf("analysis_model=='%s' & config=='G2_A_iii'", am))), 1L) == 16), "set (iii) rule A computed for all 16 cells under M0 and M1 (title)")
  one <- vapply(AM, function(am) nrow(rows(T1, sprintf("analysis_model=='%s' & config=='P2' & pass_pct > 5", am))), 1L); premise(all(one == 1), "one AUClast + Cmax cell above 5% in every analysis model (title)")
  f$one <- dderived("AUClast + Cmax cells above 5% (point), each of M0, M1, M2", T1, "config=='P2' & pass_pct > 5 :: row count per analysis model (all equal)", one, fnum(one[1], 0))
  y0 <- core_title(tx("A5b.title", f), tx("A5b.kicker", f))

  # ---- 아래: 캡션, 본문 ----
  CB <- "criteria/criteria_bias.csv"; WB <- function(pk, ep) sprintf("analysis_model=='M1' & pk_model=='%s' & scenario=='Vmax_up_080' & endpoint=='%s'", pk, ep)
  premise(all(rows(CB, "analysis_model=='M1' & scenario=='Vmax_up_080' & endpoint %in% c('AUCinf_Aiii','AUCinf_B')")$bias_dir == "toward_1") &&
          all(rows(CB, "analysis_model=='M1' & scenario=='Vmax_up_080' & endpoint %in% c('AUCinf_Aiii','AUCinf_B')")$bias_pct > 0), "Vmax up: NCA AUCinf GMR biased toward 1 (caption)")
  tv <- unique(rows(CB, "analysis_model=='M1' & scenario=='Vmax_up_080'")$target); premise(length(tv) == 1, "one target true ratio in the Vmax up cells")
  mech <- list(tvx = dv(CB, WB("k2016", "AUCinf_Aiii"), "target", 2, "", "target true AUCinf ratio, Vmax up cells"),
               bi16 = dv(CB, WB("k2016", "AUCinf_Aiii"), "bias_pct", 1, "%", "AUCinf GMR bias toward 1, set (iii) rule A, 2016 Vmax up, M1"),
               bi20 = dv(CB, WB("k2020", "AUCinf_Aiii"), "bias_pct", 1, "%", "AUCinf GMR bias toward 1, set (iii) rule A, 2020 Vmax up, M1"),
               bb16 = dv(CB, WB("k2016", "AUCinf_B"), "bias_pct", 1, "%", "AUCinf GMR bias toward 1, all estimable (rule B), 2016 Vmax up, M1"),
               bb20 = dv(CB, WB("k2020", "AUCinf_B"), "bias_pct", 1, "%", "AUCinf GMR bias toward 1, all estimable (rule B), 2020 Vmax up, M1"))
  cap <- c(tx("A5b.mech", mech), tx("A5b.caption", list(r2i = f_set("i", "r2"), ex = f_set("i", "extrap"), span = f_set("ii", "span"))))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  b <- list(p2 = dv(DE, WD("M1"), "p2_pct", 2, "%", "AUClast + Cmax, 2020 V2 up cell, M1"),
            bias = dv(DE, WD("M1"), "auclast_bias_pct", 2, "%", "AUClast GMR bias toward 1 against the true AUCinf ratio, 2020 V2 cell, M1"),
            d = a5b_signed(dv(DE, WD("M1"), "auclast_minus_ref_pp", 2, "", "AUClast term: AUClast alone minus the true-AUCinf reference (points), M1")),
            n = dcount(DE, "analysis_model=='M1'", "boundary cells in the decomposition"))
  by <- core_body(tx("A5b.body", c(b, list(nom = f$nom))), capy - 0.08)

  # ---- 왼쪽 (a) 표 ----
  RW <- list(P2 = c(T1, "P2"), A_iii = c(CG, "G2_A_iii"), A_i = c(T1, "G2_Ai"), A_ii = c(T1, "G2_Aii"), B = c(T1, "G2_B"), C_i = c(T1, "G2_Ci"))
  cell <- function(rel, cf, am) { w <- sprintf("analysis_model=='%s' & config=='%s'", am, cf)
    if (!nrow(rows(rel, w))) return(L$table$na)
    fill(L$table$cell, list(k = dcount(rel, paste(w, "& pass_pct > 5"), sprintf("%s %s cells above 5%% (point)", am, cf)),
                            mx = dext(rel, w, "pass_pct", max, 2, "%", sprintf("%s %s largest boundary type I error", am, cf)))) }
  df <- data.frame(a = unlist(L$table$rows[names(RW)]), stringsAsFactors = FALSE)
  for (am in AM) df[[am]] <- vapply(RW, function(z) cell(z[1], z[2], am), "")
  names(df) <- tx("A5b.table.head")
  hh <- 0.40; TW <- 7.35; ty <- y0 + hh; th <- by - 0.12 - ty
  deck_text(tx("A5b.head_a", f), c(GEO$ML, y0, TW, hh), size = 16, bold = TRUE, label = "label_a", gap_pt = 0)
  deck_table(df, box = c(GEO$ML, ty, TW, th), widths = c(3.3, 1.35, 1.35, 1.35), size = 14, highlight = 1:2, highlight_fill = PAL$tint_grey, label = "table_models")
  dsrc("analysis model table", c(T1, CG), "(table)")

  # ---- 오른쪽 (b) 그림: 같은 칸의 세 판정(M1) ----
  G <- L$fig
  CF <- c(ref = "AUCinf_true_only", last = "AUClast_only", p2 = "P2")
  fd <- d1[, .(k = names(CF)[match(config, CF)], est = pass_pct, lo, hi)][, y := c(ref = 3, last = 2, p2 = 1)[k]]
  fd[, col := c(ref = PAL$ink2, last = PAL$blue, p2 = PAL$blue)[k]]
  ann <- data.table(y = c(2.5, 1.5), lab = c(fill(G$d_last, list(d = a5b_signed(fnum(de1$auclast_minus_ref_pp, 2)))), fill(G$d_cmax, list(d = a5b_signed(fnum(de1$p2_minus_auclast_pp, 2))))),
                    col = c(PAL$blue, PAL$ink2))
  xl <- c(floor(min(fd$lo) * 5) / 5 - 0.2, ceiling(max(fd$hi) * 5) / 5 + 0.7)
  p <- ggplot(fd) +
    geom_vline(xintercept = nomv, linetype = "22", colour = PAL$ink2, linewidth = 0.7) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = c(ann$y - 0.19, 0.45), ymax = c(ann$y + 0.19, 0.84), fill = "white", colour = NA) +   # 효과·구간 글자 줄: 5% 점선·눈금선이 글자를 가르지 않게 흰 바탕
    geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = col), linewidth = 1.1) +
    geom_point(aes(est, y, colour = col), size = 4) +
    geom_text(aes(hi, y, label = paste0(fnum(est, 2), "%"), colour = col), hjust = -0.25, size = PT(16), family = FONT, fontface = "bold") +
    geom_text(data = fd[k == "p2"], aes(lo, y - 0.24, label = fill(G$ci, list(lo = fnum(lo, 2), hi = fnum(hi, 2))), colour = col), hjust = 0, vjust = 1, size = PT(14), family = FONT) +
    geom_text(data = ann, aes(xl[1] + 0.02, y, label = lab, colour = col), hjust = 0, size = PT(14), family = FONT) +
    scale_colour_identity() +
    scale_y_continuous(breaks = 3:1, labels = unlist(G$rows[c("ref", "last", "p2")]), limits = c(0.45, 3.4), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = xl, breaks = seq(ceiling(xl[1] * 2) / 2, ceiling(max(fd$hi) * 2) / 2, 0.5), expand = expansion(mult = 0)) +
    labs(x = fill(G$xlab, list(v = f$nom)), y = NULL, subtitle = G$sub) + theme_core(16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 14, colour = PAL$ink, lineheight = 0.95), plot.margin = margin(4, 14, 4, 4))
  xr <- GEO$ML + TW + 0.3; wr <- GEO$W - GEO$MR - xr
  deck_text(tx("A5b.head_b"), c(xr, y0, wr, hh), size = 16, bold = TRUE, label = "label_b", gap_pt = 0)
  deck_figure(p, "a5b_decomposition", c(xr, ty, wr, th), src = c(T1, DE))

  # ---- 노트 ----
  p2ci <- function(am) dci(DE, WD(am), "p2_pct", "lo", "hi", 2, "%", sprintf("AUClast + Cmax, 2020 V2 cell, %s", am))
  p2k <- function(am) dci(T1, sprintf("%s & analysis_model=='%s' & config=='P2'", W, am), "pass_pct_10k", "lo_10k", "hi_10k", 2, "%", sprintf("AUClast + Cmax, first 10,000 trials, 2020 V2 cell, %s", am))
  dec <- function(col, am, item) a5b_signed(dv(DE, WD(am), col, 3, "", sprintf("%s, 2020 V2 cell, %s", item, am)))
  deck_notes(tx("A5b.notes", c(f, b, list(
    p0 = p2ci("M0"), p1 = p2ci("M1"), p2m = p2ci("M2"), k0 = p2k("M0"), k1 = p2k("M1"),
    ntr = dint(DE, WD("M1"), "n_trials", "trials in the 2020 V2 cell"), reps = f_reps("boundary"),
    r0 = dv(DE, WD("M0"), "ref_pct", 2, "%", "true-AUCinf reference, 2020 V2 cell, M0"), r1 = dv(DE, WD("M1"), "ref_pct", 2, "%", "true-AUCinf reference, 2020 V2 cell, M1"),
    l1 = dv(DE, WD("M1"), "auclast_pct", 2, "%", "AUClast alone, 2020 V2 cell, M1"),
    t0 = dec("ref_minus_5_pp", "M0", "reference minus 5 (points)"), t1 = dec("ref_minus_5_pp", "M1", "reference minus 5 (points)"),
    a0 = dec("auclast_minus_ref_pp", "M0", "AUClast term (points)"), a1 = dec("auclast_minus_ref_pp", "M1", "AUClast term (points)"),
    c0 = dec("p2_minus_auclast_pp", "M0", "Cmax term (points)"), c1 = dec("p2_minus_auclast_pp", "M1", "Cmax term (points)"),
    s0 = dec("p2_minus_5_pp", "M0", "AUClast + Cmax minus 5 (points)"), s1 = dec("p2_minus_5_pp", "M1", "AUClast + Cmax minus 5 (points)"),
    bci = dci(DE, WD("M1"), "auclast_bias_pct", "auclast_bias_lo_pct", "auclast_bias_hi_pct", 2, "%", "AUClast GMR bias with 95% interval, M1"),
    b0 = dv(DE, WD("M0"), "auclast_bias_pct", 2, "%", "AUClast GMR bias, M0"),
    sd0 = dv(DE, WD("M0"), "ref_sd_se", 3, "", "SD/SE of the true-AUCinf reference, 2020 V2 cell, M0"), sd1 = dv(DE, WD("M1"), "ref_sd_se", 3, "", "SD/SE of the true-AUCinf reference, 2020 V2 cell, M1"),
    m20 = dv(DE, WD("M1"), "multiplier", 2, "", "V2 multiplier, 2020 model"), tr = dv(DE, WD("M1"), "target", 2, "", "target true AUCinf ratio of the V2 cells"),
    m16 = dv(DE, "pk_model=='k2016' & mechanism=='V2' & analysis_model=='M1'", "multiplier", 2, "", "V2 multiplier, 2016 model"),
    p16 = dv(T1, "pk_model=='k2016' & scenario=='V2_up_080' & analysis_model=='M1' & config=='P2'", "pass_pct", 2, "%", "AUClast + Cmax, 2016 V2 cell, M1"),
    g2v = dv(CG, sprintf("%s & analysis_model=='M1' & config=='G2_A_iii'", W), "pass_pct", 2, "%", "AUCinf (set iii) + Cmax, 2020 V2 cell, M1"),
    lo_iii = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii' & lo > 5", "AUCinf (set iii) + Cmax cells with Wilson lower bound above 5%, M1"),
    lo_b = dcount(T1, "analysis_model=='M1' & config=='G2_B' & lo > 5", "rule B cells with Wilson lower bound above 5%, M1"),
    civ = dext(CG, "analysis_model=='M1' & config=='G2_C_iv'", "pass_pct", max, 2, "%", "rule C set (iv) largest, M1"),
    aiv = dcount(CG, "analysis_model=='M1' & config=='G2_A_iv' & pass_pct > 5", "rule A set (iv) cells above 5% (point), M1"),
    aivx = dext(CG, "analysis_model=='M1' & config=='G2_A_iv'", "pass_pct", max, 2, "%", "rule A set (iv) largest, M1"),
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), split = f_split()), mech)))
  deck_end()
}
