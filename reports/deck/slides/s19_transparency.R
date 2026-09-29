# S19 논거 ③ 투명성: M1에서 AUC0-last + Cmax(P2)의 경계 1종 오류가 5%를 넘는 유일한 칸(2020 모델, 말초 분포용적 V2 증가, 참 AUC0-inf 비 0.80)을
# 사실대로 공개한다: 기전(V2 배율, 2016 모델의 같은 참 비 배율), 크기와 구간(M0·M1, 처음 10,000회와 20,000회), 같은 시험 안의 분해(AUC0-last 편향),
# 같은 칸의 AUC0-inf 구성과의 대비. 시험 모집단(건강인 60~90 kg, 체중 층화 배정, B0)만.
# 자료: results/oc_models/type1_models.csv(분석 모형별 경계 1종 오류; *_10k = 처음 10,000회), extension_decision_models.csv(확대 판정),
#       p2_decomposition_models.csv(분해, AUC0-last 편향), sd_se_models.csv(NCA AUC0-inf 편향), criteria/criteria_g2_type1.csv(세트 (iii)·(iv) 구성).
# 주의: 확대는 M0의 Wilson 구간이 5%를 포함해 정해졌다. M1은 처음 10,000회에서 이미 하한 > 5%라 스스로 확대 조건을 만들지 않았고,
#       등록 규칙(prereg section1, D-048 규칙을 분석 모형마다 적용)대로 세 분석 모형(M0·M1·M2)을 모두 20,000회로 보고한다.
#       D-048 확대는 사전 고정 설계 뒤에 추가한 사후 규칙이라(보고서 3.7절, 부록 D) 10,000회 값을 함께 적는다.
# 규칙 C는 통계분석계획 제안의 세트 (i)로 인용하고(명목), 세트 (ii)~(iv)는 세트별로 표시한다. 점추정 칸 수와 Wilson 분류 칸 수는 이름을 붙여 따로 적는다.
s19_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?%?$", p)) p else paste0("+", p)
# 분류 코드(영문)를 추적 행으로 남기고 한국어 표기로 돌려준다
s19_cls <- function(rel, where, item, col = "class") {
  r <- row1(rel, where); x <- as.character(r[[col]]); .record(item, rel, sprintf("%s :: %s", where, col), x, x)
  lab <- DK$txt$S19$cls[[x]]; premise(!is.null(lab), paste("class label for", x)); lab
}
# 문자열 열(영문)을 그대로 추적해 인쇄
s19_dtxt <- function(rel, where, col, item) { r <- row1(rel, where); x <- as.character(r[[col]]); .record(item, rel, sprintf("%s :: %s", where, col), x, x); x }

slide_S19 <- function() {
  T1 <- "oc_models/type1_models.csv"; DE <- "oc_models/p2_decomposition_models.csv"; EX <- "oc_models/extension_decision_models.csv"
  SS <- "oc_models/sd_se_models.csv"; CG <- "criteria/criteria_g2_type1.csv"; EC <- "oc_models/expectations_check.csv"; CB <- "criteria/criteria_bias.csv"
  W <- "pk_model=='k2020' & scenario=='V2_up_080'"
  WD <- function(am, m = "k2020") sprintf("pk_model=='%s' & mechanism=='V2' & analysis_model=='%s'", m, am)    # 보고서 6.2절·Q7과 같은 행 조건
  WT <- function(am, cf, w = W) sprintf("%s & analysis_model=='%s' & config=='%s'", w, am, cf)
  deck_slide("S19", tag = "sim")
  L <- DK$txt$S19
  nom_num <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2

  # ---- 전제: 문장이 기대는 사실 ----
  premise(nrow(rows(T1, "analysis_model=='M1' & config=='P2'")) == 16, "16 boundary cells under M1")
  ov <- rows(T1, "analysis_model=='M1' & config=='P2' & pass_pct > 5")
  premise(nrow(ov) == 1 && ov$pk_model == "k2020" && ov$scenario == "V2_up_080", "the only boundary cell with P2 above 5% under M1 is the 2020 model V2 up cell (title: one cell)")
  premise(nrow(rows(T1, "analysis_model=='M1' & config=='P2' & class!='conservative'")) == 1, "every other M1 cell is conservative")
  tw <- rows(DE, "auclast_bias_dir=='toward_1'")
  premise(nrow(tw) == 3 && all(tw$pk_model == "k2020" & tw$scenario == "V2_up_080"), "AUC0-last is biased toward 1 only in this cell, under each analysis model (card)")
  v2 <- rows(DE, "mechanism=='V2' & analysis_model=='M1'")
  premise(nrow(v2) == 2 && all(abs(v2$target - 0.8) < 1e-9), "the two V2 cells (2016 and 2020 models) target the same true AUC0-inf ratio")
  premise(row1(T1, WT("M1", "P2", "pk_model=='k2016' & scenario=='V2_up_080'"))$class == "conservative", "the 2016 model V2 cell is conservative under M1")
  ex <- rows(EX, W)
  premise(isTRUE(ex[analysis_model == "M0", triggers]) && !isTRUE(ex[analysis_model == "M1", triggers]) && all(ex$selected) && ex[analysis_model == "M1", lo] > 5 &&
            ex[analysis_model == "M0", lo] <= 5 && ex[analysis_model == "M0", hi] >= 5, "extension triggered by M0 (interval includes 5%); M1 did not trigger (lower bound above 5%); both extended")
  premise(all(rows(T1, sprintf("%s & analysis_model %%in%% c('M0','M1') & config %%in%% c('G2_Ai','G2_Aii','G2_B')", W))$class == "conservative"), "rules A and B conservative in this cell under M0 and M1")
  premise(all(rows(SS, sprintf("%s & analysis_model=='M1' & endpoint %%in%% c('AUCinf_A','AUCinf_Ai','AUCinf_B')", W))$bias_dir == "away_from_1"), "NCA AUC0-inf rules A and B biased away from 1 in this cell (M1)")
  dd <- rows(DE, W); premise(nrow(dd) == 3 && all(abs(dd$ref_minus_5_pp + dd$auclast_minus_ref_pp + dd$p2_minus_auclast_pp - dd$p2_minus_5_pp) < 1e-9), "the decomposition adds up to P2 minus 5 in each analysis model")
  premise(all(dd$n_trials == 20000) && all(rows(T1, sprintf("%s & analysis_model %%in%% c('M0','M1')", W))$n_trials == 20000), "20,000 trials in this cell")
  premise(row1(DE, WD("M1"))$auclast_minus_ref_pp > row1(DE, WD("M1"))$p2_minus_5_pp, "the AUC0-last term exceeds the excess over 5% (title: the excess comes from AUC0-last)")

  # ---- 제목 ----
  nom <- f_nominal()
  m20 <- dv(DE, WD("M0"), "multiplier", 2, "", "V2 multiplier k2020")
  premise(row1(DE, WD("M1"))$ref_minus_5_pp < 0 && row1(DE, WD("M1"))$p2_minus_5_pp > 0, "the unbiased reference is below 5% and P2 above 5% in this cell (title: the excess comes from AUC0-last)")
  f <- list(nom = nom, m20 = m20, p2 = dv(T1, WT("M1", "P2"), "pass_pct", 2, "%", "P2 boundary type I error, 2020 model V2 cell, M1, 20,000 trials"),
            bias = s19_signed(dv(DE, WD("M1"), "auclast_bias_pct", 2, "%", "2020 V2 AUC0-last bias M1")))
  deck_kicker(tx("S19.kicker", list(nom = nom))); deck_title(tx("S19.title", f))

  # ---- 왼쪽: 같은 칸의 구성별 경계 1종 오류(M1, M0), 20,000회와 사전 고정 10,000회(P2) ----
  CF <- c("P2", "P2_10k", "AUClast_only", "AUCinf_true_only", "G2_Ai", "G2_Aii", "G2_B", "G2_Ci", "G2_Cii")
  d <- rows(T1, sprintf("%s & analysis_model %%in%% c('M0','M1')", W))
  fd <- rbind(d[config %in% setdiff(CF, "P2_10k"), .(config, am = analysis_model, est = pass_pct, lo, hi)],
              d[config == "P2", .(config = "P2_10k", am = analysis_model, est = pass_pct_10k, lo = lo_10k, hi = hi_10k)])
  premise(nrow(fd) == 2 * length(CF), "figure rows: nine configurations x two analysis models")
  n20 <- unique(d$n_trials); n10 <- unique(d$n_trials_10k); premise(length(n20) == 1 && length(n10) == 1, "one trial count per column")
  lab <- vapply(CF, function(k) fill(L$fig$cfg[[k]], list(n = fint(if (k == "P2_10k") n10 else n20))), "")
  gap <- c(0, 0, 0.35, 0.35, 0.7, 0.7, 0.7, 0.7, 0.7)                     # 묶음 사이 간격: P2 | 단일 평가변수 | AUC0-inf + Cmax
  ypos <- setNames(length(CF) - seq_along(CF) + 1 - gap, CF)
  AML <- unlist(DK$txt$common$analysis_models[c("M1", "M0")])
  fd[, y := ypos[config] + ifelse(am == "M1", 0.17, -0.17)]
  fd[, amf := factor(AML[am], levels = AML)]
  # 오른쪽 수치 열: 행마다 M1(파랑), M0(주황) 점추정. 그림 안 숫자는 자료에서 바로 그린다
  xmax <- max(ceiling(max(fd$hi)), nom_num + 2); xc <- c(M1 = xmax + 0.95, M0 = xmax + 2.05)
  txt <- fd[, .(amf, x = xc[am], y = ypos[config], lab = fnum(est, 2))]
  hdr <- data.table(amf = factor(AML, levels = AML), x = unname(xc[c("M1", "M0")]), y = max(ypos) + 0.75, lab = c("M1", "M0"))
  FW <- 6.4; FH <- 2.85
  p <- ggplot(fd, aes(x = est, y = y, colour = amf, shape = amf)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = min(ypos[c("P2", "P2_10k")]) - 0.45, ymax = max(ypos[c("P2", "P2_10k")]) + 0.45, fill = PAL$tint_blue, alpha = 0.8) +
    geom_vline(xintercept = nom_num, linetype = "22", colour = PAL$ink2, linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7) +
    geom_point(size = 2.6) +
    geom_text(data = txt, aes(x = x, y = y, label = lab), hjust = 1, size = 3.8, family = FONT, show.legend = FALSE, inherit.aes = FALSE, colour = ifelse(txt$amf == AML[1], PAL$blue, PAL$orange)) +
    geom_text(data = hdr, aes(x = x, y = y, label = lab), hjust = 1, size = 3.8, family = FONT, fontface = "bold", show.legend = FALSE, inherit.aes = FALSE, colour = PAL$ink2) +
    scale_colour_manual(values = c(PAL$blue, PAL$orange)) + scale_shape_manual(values = c(16, 17)) +
    scale_y_continuous(breaks = unname(ypos), labels = unname(lab), expand = expansion(add = 0)) +
    scale_x_continuous(limits = c(0, max(xc) + 0.05), breaks = seq(0, xmax, by = 1), expand = expansion(add = c(0.1, 0.1))) +
    coord_cartesian(ylim = c(min(ypos) - 0.45, max(ypos) + 0.95), clip = "off") +
    labs(x = fill(L$fig$xlab, list(nom = nom)), y = NULL) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          axis.title.x = element_text(hjust = 1, margin = margin(4, 0, 0, 0)),
          panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), axis.text.y = element_text(colour = PAL$ink, size = 12))
  deck_figure(p, "s19_v2_cell", c(GEO$ML, GEO$BODY_TOP, FW, FH), src = T1)

  # ---- 왼쪽 아래: 같은 칸의 AUC0-inf 구성(세트별), 16칸 전체(점추정, Wilson 초과) ----
  WG <- function(am) sprintf("%s & analysis_model=='%s'", W, am)
  GA <- c("G2_A_i", "G2_A_ii", "G2_A_iii", "G2_A_iv", "G2_B"); GC <- c("G2_C_ii", "G2_C_iii", "G2_C_iv")
  cg1 <- rows(CG, WG("M1"))
  premise(all(cg1[config %in% GA]$class == "conservative"), "rules A (i) to (iv) and B + Cmax are conservative in this cell (M1)")
  premise(cg1[config == "G2_C_i"]$class == "nominal" && all(cg1[config %in% GC]$class == "exceeding"), "rule C (i) nominal and C (ii) to (iv) exceeding in this cell (M1)")
  premise(all(rows(CB, sprintf("%s & analysis_model=='M1' & config %%in%% c('A_i','A_ii','A_iii','A_iv','B')", W))$bias_dir == "away_from_1") &&
          all(rows(CB, sprintf("%s & analysis_model=='M1' & config %%in%% c('C_i','C_ii','C_iii','C_iv')", W))$bias_dir == "toward_1"),
          "NCA AUC0-inf biased away from 1 for rules A (all sets) and B, toward 1 for rule C (all sets), this cell, M1")
  premise(cg1[config == "G2_C_iv"]$pass_pct == max(cg1[config %in% GC]$pass_pct), "C (iv) is the largest of rule C (ii) to (iv) in this cell")
  cnt <- function(cf, rel = T1, col = "pass_pct") dcount(rel, sprintf("analysis_model=='M1' & config=='%s' & %s > 5", cf, col), sprintf("M1 %s boundary cells above 5%% (%s)", cf, if (col == "lo") "Wilson lower bound" else "point"))
  rng_int <- function(k, item, rel, loc) { k <- as.integer(k); dderived(item, rel, loc, range(k), if (min(k) == max(k)) as.character(k[1]) else sprintf("%d~%d", min(k), max(k))) }
  nab <- vapply(c("G2_Ai", "G2_Aii", "G2_B"), cnt, ""); nabw <- vapply(c("G2_Ai", "G2_Aii", "G2_B"), cnt, "", col = "lo")
  b <- list(nom = nom, n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"),
            ab = drange(CG, sprintf("%s & config %%in%% c(%s)", WG("M1"), paste(sprintf("'%s'", GA), collapse = ",")), "pass_pct", 2, "%", "AUC0-inf rules A (i) to (iv) and B + Cmax in the 2020 V2 cell, M1"),
            ci = dv(CG, sprintf("%s & config=='G2_C_i'", WG("M1")), "pass_pct", 2, "%", "G2_C_i, 2020 V2 cell, M1"),
            civ = dv(CG, sprintf("%s & config=='G2_C_iv'", WG("M1")), "pass_pct", 2, "%", "G2_C_iv, 2020 V2 cell, M1"),
            nab = rng_int(nab, "M1 boundary cells above 5% (point), rules A (i), A (ii), B: range of the three counts", T1, "count of rows [analysis_model=='M1' & config in G2_Ai, G2_Aii, G2_B & pass_pct > 5] per config :: range"),
            nabw = rng_int(nabw, "M1 boundary cells with Wilson lower bound above 5%, rules A (i), A (ii), B: range of the three counts", T1, "count of rows [analysis_model=='M1' & config in G2_Ai, G2_Aii, G2_B & lo > 5] per config :: range"),
            one = cnt("P2"))
  premise(b$one == "1", "one P2 cell above 5% under M1 (bullet)")
  BY <- GEO$BODY_TOP + FH + 0.1
  deck_bullets(tx("S19.bullets", b), box = c(GEO$ML, BY, FW, GEO$BODY_BOTTOM - BY), size = 16)

  # ---- 오른쪽: 확대 경위, 기전, 같은 시험 안의 분해(M1, M0) ----
  XR <- GEO$ML + FW + 0.3; WR <- GEO$W - GEO$MR - XR
  ex_v <- function(am, col, item) dv(EX, sprintf("%s & analysis_model=='%s'", W, am), col, 2, "%", item)
  ext <- list(nom = nom, n10 = dint(EX, sprintf("%s & analysis_model=='M1'", W), "n_before", "pre-specified trials"), n20 = dint(EX, sprintf("%s & analysis_model=='M1'", W), "n_after", "trials after extension"),
              m1k = ex_v("M1", "pass_pct", "P2 at the pre-specified 10,000 trials, 2020 V2 cell, M1"), m0k = ex_v("M0", "pass_pct", "P2 at the pre-specified 10,000 trials, 2020 V2 cell, M0"),
              m0ci = sprintf("%s~%s", sub("%$", "", ex_v("M0", "lo", "Wilson lower bound, P2 at 10,000 trials, M0")), ex_v("M0", "hi", "Wilson upper bound, P2 at 10,000 trials, M0")))
  premise(all(rows(T1, sprintf("%s & config=='P2'", W))$n_trials == 20000) && nrow(rows(T1, sprintf("%s & config=='P2'", W))) == 3, "all three analysis models are reported at 20,000 trials in this cell (card)")
  EH <- 2.0
  deck_text(tx("S19.ext", ext), c(XR, GEO$BODY_TOP, WR, EH), size = 16, label = "text_ext", bg = PAL$tint_blue, geom = "roundRect", gap_pt = 4)
  mech <- list(m20 = m20, m16 = dv(DE, WD("M0", "k2016"), "multiplier", 2, "", "V2 multiplier k2016"),
               tr = dv(DE, WD("M1"), "target", 2, "", "true AUC0-inf ratio targeted in the V2 cells"),
               p16 = dv(T1, WT("M1", "P2", "pk_model=='k2016' & scenario=='V2_up_080'"), "pass_pct", 2, "%", "P2, 2016 model V2 cell, M1"))
  MY <- GEO$BODY_TOP + EH + 0.12; MH <- 1.35
  deck_text(tx("S19.mech", mech), c(XR, MY, WR, MH), size = 16, label = "text_mech", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 4)
  dec <- function(col, am, item) s19_signed(dv(DE, WD(am), col, 3, "", sprintf("%s, 2020 V2 cell, %s", item, am)))
  cols <- c(ref = "ref_minus_5_pp", last = "auclast_minus_ref_pp", cmax = "p2_minus_auclast_pp", p2 = "p2_minus_5_pp")
  itm <- c(ref = "unbiased reference minus 5 (points)", last = "AUC0-last effect (points)", cmax = "Cmax effect (points)", p2 = "P2 minus 5 (points)")
  df <- data.frame(a = vapply(names(cols), function(k) fill(L$table$rows[[k]], list(nom = nom)), ""),
                   b = vapply(names(cols), function(k) dec(cols[[k]], "M1", itm[[k]]), ""),
                   c = vapply(names(cols), function(k) dec(cols[[k]], "M0", itm[[k]]), ""), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S19.table.head")
  TY <- MY + MH + 0.12
  deck_table(df, box = c(XR, TY, WR, GEO$BODY_BOTTOM - TY), widths = c(3.3, 0.95, 0.95), size = 12, highlight = 2, label = "table_decomp")

  # ---- 노트 ----
  g2 <- function(cf, am) dv(CG, sprintf("%s & analysis_model=='%s' & config=='%s'", W, am, cf), "pass_pct", 2, "%", sprintf("%s, 2020 V2 cell, %s", cf, am))
  nb <- function(ep) dv(SS, sprintf("%s & analysis_model=='M1' & endpoint=='%s'", W, ep), "bias_pct", 2, "%", sprintf("NCA %s bias, 2020 V2 cell, M1", ep))
  p2ci <- function(am) dci(DE, WD(am), "p2_pct", "lo", "hi", 2, "%", sprintf("2020 V2 cell P2 %s", am))
  p2k <- function(am) dci(T1, WT(am, "P2"), "pass_pct_10k", "lo_10k", "hi_10k", 2, "%", sprintf("P2 first 10,000 trials, 2020 V2 cell, %s", am))
  nts <- list(
    nom = nom, m20 = m20, m16 = mech$m16, p16 = mech$p16, n = b$n, npm = dcount(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020'", "boundary cells per model"),
    auc = dv(DE, WD("M1"), "auc_ratio", 3, "", "true AUC0-inf ratio, 2020 V2 cell"), cmr = dv(DE, WD("M1"), "cmax_ratio", 3, "", "true Cmax ratio, 2020 V2 cell"),
    n20 = ext$n20, n10 = ext$n10,
    p0 = p2ci("M0"), p1 = p2ci("M1"), p2m = p2ci("M2"),
    c0 = s19_cls(DE, WD("M0"), "class M0"), c1 = s19_cls(DE, WD("M1"), "class M1"), c2 = s19_cls(DE, WD("M2"), "class M2"),
    k0 = p2k("M0"), k1 = p2k("M1"), k2 = p2k("M2"),
    ck0 = s19_cls(T1, WT("M0", "P2"), "P2 class at the first 10,000 trials, M0", "class_10k"), ck1 = s19_cls(T1, WT("M1", "P2"), "P2 class at the first 10,000 trials, M1", "class_10k"),
    ck2 = s19_cls(T1, WT("M2", "P2"), "P2 class at the first 10,000 trials, M2", "class_10k"),
    b1 = dci(DE, WD("M1"), "auclast_bias_pct", "auclast_bias_lo_pct", "auclast_bias_hi_pct", 2, "%", "AUC0-last bias with 95% interval, M1"),
    b0 = dci(DE, WD("M0"), "auclast_bias_pct", "auclast_bias_lo_pct", "auclast_bias_hi_pct", 2, "%", "AUC0-last bias with 95% interval, M0"),
    r0 = dv(DE, WD("M0"), "ref_pct", 2, "%", "unbiased reference M0"), r1 = dv(DE, WD("M1"), "ref_pct", 2, "%", "2020 V2 unbiased reference M1"),
    l1 = dv(DE, WD("M1"), "auclast_pct", 2, "%", "2020 V2 AUC0-last alone M1"),
    s0 = dv(DE, WD("M0"), "ref_sd_se", 3, "", "SD/SE of the unbiased reference, M0"), s1 = dv(DE, WD("M1"), "ref_sd_se", 3, "", "SD/SE of the unbiased reference, M1"),
    ai1 = g2("G2_A_i", "M1"), aii1 = g2("G2_A_ii", "M1"), aiii1 = g2("G2_A_iii", "M1"), aiv1 = g2("G2_A_iv", "M1"), bb1 = g2("G2_B", "M1"),
    ci1 = g2("G2_C_i", "M1"), cii1 = g2("G2_C_ii", "M1"), ciii1 = g2("G2_C_iii", "M1"), civ1 = g2("G2_C_iv", "M1"),
    ai0 = g2("G2_A_i", "M0"), aii0 = g2("G2_A_ii", "M0"), aiii0 = g2("G2_A_iii", "M0"), aiv0 = g2("G2_A_iv", "M0"), bb0 = g2("G2_B", "M0"),
    ci0 = g2("G2_C_i", "M0"), cii0 = g2("G2_C_ii", "M0"), ciii0 = g2("G2_C_iii", "M0"), civ0 = g2("G2_C_iv", "M0"),
    nba = nb("AUCinf_A"), nbai = nb("AUCinf_Ai"), nbb = nb("AUCinf_B"), nbc = nb("AUCinf_C"),
    nai = nab[["G2_Ai"]], naii = nab[["G2_Aii"]], nbB = nab[["G2_B"]], wai = nabw[["G2_Ai"]], waii = nabw[["G2_Aii"]], wbB = nabw[["G2_B"]],
    nci = cnt("G2_Ci"), wci = cnt("G2_Ci", col = "lo"), ncii = dcount(CG, "analysis_model=='M1' & config=='G2_C_ii' & pass_pct > 5", "M1 rule C (ii) cells above 5% (point)"),
    ka16 = dv(T1, "pk_model=='k2016' & scenario=='ka_down_080' & analysis_model=='M1' & config=='AUCinf_true_only'", "pass_pct", 2, "%", "unbiased reference, 2016 ka down, M1"),
    ka20 = dv(T1, "pk_model=='k2020' & scenario=='ka_down_080' & analysis_model=='M1' & config=='AUCinf_true_only'", "pass_pct", 2, "%", "unbiased reference, 2020 ka down, M1"),
    kap2 = drange(T1, "scenario=='ka_down_080' & analysis_model=='M1' & config=='P2'", "pass_pct", 2, "%", "P2 in the ka down cells, M1"),
    exp_c = s19_dtxt(EC, "grepl('V2', expectation)", "criterion", "registered expectation for the V2 cell (criterion)"),
    exp_ok = s19_dtxt(EC, "grepl('V2', expectation)", "consistent", "registered expectation consistent"))
  premise(nrow(rows(T1, "scenario=='ka_down_080' & analysis_model=='M1' & config=='P2' & pass_pct > 0")) == 0, "P2 passes in no trial in the ka down cells (notes)")
  premise(nrow(rows(T1, "analysis_model=='M1' & config=='G2_Ci' & pass_pct > 5 & pk_model=='k2020' & scenario=='V2_up_080'")) == 1 && row1(T1, WT("M1", "G2_Ci"))$class == "nominal", "the rule C (i) point-estimate cell above 5% is this cell and is nominal (notes)")
  premise(row1(T1, WT("M0", "P2"))$class_10k == "nominal" && row1(T1, WT("M1", "P2"))$class_10k == "exceeding", "10,000-trial classes: M0 nominal, M1 exceeding (notes)")
  deck_notes(tx("S19.notes", nts))
  deck_end()
}
