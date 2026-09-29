# S15 논거 ② 경계 1종 오류, 기전별: AUC0-last + Cmax(P2, 제안) 대 AUC0-inf + Cmax(규칙 A (ii) 사전 규정 지침 기본값, 규칙 B).
# 그림: results/oc_models/type1_models.csv, M1(체중 층 포함), 경계 16칸(모델당 8칸, 기전 x 방향), Wilson 95% 구간, 색 + 표식 모양.
# 오른쪽: P2 헤드라인(M1에서 명목 이하 칸 수; 보고서 요약 헤드라인 p2_atnom과 같은 조건), AUC0-inf 구성의 점추정 > 5% 칸 수와 최대(M1, M0 병기).
# 시험 모집단 재생성 시험(건강인, 체중 층화)만. 연장 칸(2020 모델 V2)은 20,000회, 10,000회 값은 노트.

# 보고서 p2_atnom(am)과 같은 조건·locator: 명목 이하(보수적 + 명목) 칸 수
s15_atnom <- function(am) {
  T1f <- "oc_models/type1_models.csv"; r <- rows(T1f, sprintf("analysis_model=='%s' & config=='P2'", am)); premise(nrow(r) == 16, "16 P2 cells")
  k <- sum(r$class %in% c("conservative", "nominal"))
  dderived(sprintf("P2 cells at or below nominal, %s", am), T1f, sprintf("analysis_model=='%s' & config=='P2' :: count class in (conservative, nominal)", am), k, as.character(k))
}
# 여러 구성 각각의 점추정 > 5% 칸 수를 범위로(보고서 g_rng_txt, g_rng_m1_txt와 같은 계산)
s15_cnt_rng <- function(am, cfs, item) {
  T1f <- "oc_models/type1_models.csv"
  k <- vapply(cfs, function(cf) nrow(rows(T1f, sprintf("analysis_model=='%s' & config=='%s' & pass_pct > 5", am, cf))), 1L)
  dderived(item, T1f, sprintf("analysis_model=='%s' & config in (%s) & pass_pct > 5 :: range over configs of row counts", am, paste(cfs, collapse = ", ")), k, rng_fmt(min(k), max(k), 0))
}

slide_S15 <- function() {
  T1 <- "oc_models/type1_models.csv"; CB <- "criteria/criteria_bias.csv"; SS <- "oc_models/sd_se_models.csv"; PC <- "oc_models/type1_paired_change.csv"
  CF <- c("P2", "G2_Aii", "G2_B"); USUAL <- c("G2_Ai", "G2_Aii", "G2_B")
  SC <- c("F_down_080", "F_up_125", "ke_up_080", "ke_down_125", "Vmax_up_080", "Vmax_down_125", "V2_up_080", "ka_down_080")   # 그림 위에서 아래 순서
  deck_slide("S15", tag = "sim")

  # ---- 전제 ----
  d <- copy(rows(T1, sprintf("analysis_model=='M1' & config %%in%% c(%s)", paste(sprintf("'%s'", CF), collapse = ","))))
  for (pk in c("k2016", "k2020")) for (cf in CF) premise(setequal(d[pk_model == pk & config == cf, scenario], SC), sprintf("8 boundary cells, %s %s", pk, cf))
  p2 <- d[config == "P2"]; ex <- p2[class == "exceeding"]
  premise(sum(p2$class == "conservative") == 15 && sum(p2$class == "nominal") == 0 && nrow(ex) == 1 && ex$pk_model == "k2020" && ex$scenario == "V2_up_080",
          "M1 P2: 15 conservative, the 16th cell (2020 V2 up) exceeding (text)")
  premise(all(rows(T1, "grepl('^ka_', scenario) & config!='AUCinf_true_only' & config!='AUClast_only'")$pass_pct == 0), "ka-down cells: every Cmax configuration 0% (text)")
  premise(all(rows(T1, "grepl('^ka_', scenario)")$cmax_ratio < min(unlist(.read("config/trial_design.yaml")$be$limits))), "ka-down cells: true Cmax ratio below the lower equivalence limit (notes)")
  for (pk in c("k2016", "k2020")) for (cf in c("G2_Aii", "G2_B")) {
    r <- d[pk_model == pk & config == cf]; premise(r[which.max(pass_pct)]$mechanism == "Vmax", sprintf("%s %s highest in a Vmax cell (text)", pk, cf))
    premise(r[mechanism == "V2"]$class == "conservative", sprintf("%s %s conservative in the V2 cell (text)", pk, cf))
    premise(rows(CB, sprintf("analysis_model=='M1' & config=='%s' & pk_model=='%s' & scenario=='V2_up_080'", sub("G2_", "", sub("Aii", "A_ii", cf)), pk))$bias_dir == "away_from_1", sprintf("%s %s: AUC0-inf GMR bias away from 1 in the V2 cell (notes)", pk, cf))
  }

  # ---- 제목 ----
  f <- list(n = dcount(T1, "analysis_model=='M1' & config=='P2'", "M1 P2 boundary cells"),
            k = headline(s15_atnom("M1")), nom = f_nominal(),
            r = s15_cnt_rng("M1", USUAL, "M1 AUC0-inf + Cmax, rules A (i), A (ii) and B: cells above 5% (point), range over configurations"))
  deck_kicker(tx("S15.kicker")); deck_title(tx("S15.title", f))

  # ---- 그림: 경계 16칸 x 구성 3개(M1), Wilson 95% 구간 ----
  L <- DK$txt$S15$fig
  nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  d[, row := length(SC) + 1 - match(scenario, SC)]
  off <- c(P2 = 0.24, G2_Aii = 0, G2_B = -0.24); d[, y := row + off[config]]
  d[, cfg := factor(config, levels = CF, labels = unlist(L$cfg[CF]))]
  lab_y <- unique(d[, .(row, mechanism, direction, target)])[order(row)]
  premise(nrow(lab_y) == length(SC), "one label per boundary row (same in both models)")
  lab_y[, lab := vapply(seq_len(.N), function(i) fill(L$cell, list(mech = L$mech[[mechanism[i]]], dir = L$dir[[direction[i]]], tgt = fnum(target[i], 2))), "")]
  COL <- setNames(c(PAL$blue, PAL$orange, PAL$ink2), levels(d$cfg)); SHP <- setNames(c(16, 17, 15), levels(d$cfg))
  kar <- unique(d[mechanism == "ka", .(model, row)])
  exl <- d[config == "P2" & class == "exceeding"]
  xmax <- ceiling(max(d$hi) / 5) * 5
  p <- ggplot(d, aes(y = y, colour = cfg, shape = cfg)) +
    geom_hline(yintercept = seq(1.5, length(SC) - 0.5, 1), colour = PAL$grid, linewidth = 0.35) +
    geom_vline(xintercept = nomv, linetype = "22", colour = PAL$ink2, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.8) +
    geom_point(aes(x = pass_pct), size = 2.4) +
    geom_text(data = kar, aes(x = nomv + 0.7, y = row + 0.02, label = L$ka), inherit.aes = FALSE, hjust = 0, size = 3.7, family = FONT, colour = PAL$ink2) +
    geom_text(data = exl, aes(x = hi + 0.5, y = y + 0.05, label = paste0(fnum(pass_pct, 2), "%")), hjust = 0, size = 3.9, family = FONT, colour = PAL$blue, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ model, nrow = 1) +
    scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
    scale_y_continuous(breaks = lab_y$row, labels = lab_y$lab, limits = c(0.55, length(SC) + 0.45), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, 5), expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = fill(L$xlab, list(nom = paste0(fnum(nomv, 0), "%"))), y = NULL) + theme_deck(13) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "top", legend.justification = "left", legend.text = element_text(size = 12),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"), strip.text = element_text(hjust = 0, size = 13),
          axis.text.y = element_text(size = 12, colour = PAL$ink), panel.spacing.x = grid::unit(14, "pt"), axis.title.x = element_text(size = 12, margin = margin(4, 0, 0, 0))) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE))
  FW <- 7.1; FH <- 4.4
  deck_figure(p, "s15_boundary_by_mechanism", c(GEO$ML, GEO$BODY_TOP, FW, FH), src = T1)
  deck_text(tx("S15.gloss"), c(GEO$ML, GEO$BODY_TOP + FH + 0.05, FW, GEO$BODY_BOTTOM - GEO$BODY_TOP - FH - 0.05), size = 16, color = PAL$ink2, label = "gloss", gap_pt = 0)

  # ---- 오른쪽: 헤드라인 카드, 구성별 표, 요점 ----
  xr <- GEO$ML + FW + 0.3; wr <- GEO$W - GEO$MR - xr
  p2ci <- dci(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M1 P2, 2020 model V2 up cell")
  sh <- 1.62
  deck_stat(tx("S15.stat.value", f), tx("S15.stat.label", list(exc = dcount(T1, "analysis_model=='M1' & config=='P2' & class=='exceeding'", "M1 P2 cells classified exceeding"),
                                                              p2 = dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2, 2020 model V2 up cell"),
                                                              m0 = s15_atnom("M0"))),
            c(xr, GEO$BODY_TOP, wr, sh), value_size = 36)
  TR <- DK$txt$S15$table
  cell <- function(am, cf) sprintf("%s · %s", dcount(T1, sprintf("analysis_model=='%s' & config=='%s' & pass_pct > 5", am, cf), sprintf("%s %s cells above 5%% (point)", am, cf)),
                                   dext(T1, sprintf("analysis_model=='%s' & config=='%s'", am, cf), "pass_pct", max, 2, "%", sprintf("%s %s largest boundary pass rate", am, cf)))
  rc <- c("G2_Aii", "G2_Ai", "G2_B", "G2_Cii", "P2")
  df <- data.frame(a = unlist(TR$rows[rc]), b = vapply(rc, function(cf) cell("M1", cf), ""), c = vapply(rc, function(cf) cell("M0", cf), ""), stringsAsFactors = FALSE)
  names(df) <- tx("S15.table.head", f)
  ty <- GEO$BODY_TOP + sh + 0.12; th <- 2.3
  deck_table(df, box = c(xr, ty, wr, th), widths = c(2.4, 1.21, 1.21), size = 13, highlight = 5, highlight_fill = PAL$tint_blue)
  by <- ty + th + 0.08
  deck_bullets(tx("S15.bullets", list(vmax = drange(T1, "analysis_model=='M1' & config=='G2_Aii' & mechanism=='Vmax'", "pass_pct", 2, "%", "M1 G2_Aii, Vmax cells, range"))),
               box = c(xr, by, wr, GEO$BODY_BOTTOM - by), size = 16, gap_pt = 4)

  # ---- 노트 ----
  W20 <- "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'"
  mx <- function(am, cf) { r <- rows(T1, sprintf("analysis_model=='%s' & config=='%s'", am, cf))[which.max(pass_pct)]
    dci(T1, sprintf("analysis_model=='%s' & config=='%s' & pk_model=='%s' & scenario=='%s'", am, cf, r$pk_model, r$scenario), "pass_pct", "lo", "hi", 2, "%", sprintf("%s %s largest", am, cf)) }
  premise(rows(T1, "analysis_model=='M1' & config=='G2_Aii'")[which.max(pass_pct)][, pk_model == "k2020" & scenario == "Vmax_down_125"] &&
          rows(T1, "analysis_model=='M1' & config=='G2_B'")[which.max(pass_pct)][, pk_model == "k2016" & scenario == "Vmax_up_080"], "cells of the M1 maxima named in the notes")
  lo5 <- function(am, cf) dcount(T1, sprintf("analysis_model=='%s' & config=='%s' & lo > 5", am, cf), sprintf("%s %s cells with lower bound above 5%%", am, cf))
  cls <- function(am, k) dcount(T1, sprintf("analysis_model=='%s' & config=='P2' & class=='%s'", am, k), sprintf("%s P2 cells %s", am, k))
  deck_notes(tx("S15.notes", c(f, list(
    tgt = dcfg("oc_design.yaml", "boundary_targets", "boundary true AUC0-inf ratios", function(x) paste(fnum(x, 2), collapse = ", ")),
    vm = dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016' & scenario=='Vmax_up_080'", "multiplier", 2, "", "Vmax multiplier, 2016 model, ratio 0.80 cell"),
    vmt = dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016' & scenario=='Vmax_up_080'", "target", 2, "", "target ratio of the Vmax up cell"),
    npm = dcount(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016'", "boundary cells per PK model"),
    c1 = cls("M1", "conservative"), n1 = cls("M1", "nominal"), e1 = cls("M1", "exceeding"),
    rng = drange(T1, "analysis_model=='M1' & config=='P2'", "pass_pct", 2, "%", "M1 P2 boundary range"),
    v2m = dv(T1, W20, "multiplier", 2, "", "V2 multiplier, 2020 model"), v2r = dv(T1, W20, "auc_ratio", 4, "", "true AUC0-inf ratio, 2020 model V2 cell"),
    p2 = p2ci, n20 = dint(T1, W20, "n_trials", "trials, 2020 model V2 cell"), reps = f_reps("boundary"),
    p2k = dci(T1, W20, "pass_pct_10k", "lo_10k", "hi_10k", 2, "%", "M1 P2, 2020 model V2 cell, first 10,000 trials"),
    v2m16 = dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016' & scenario=='V2_up_080'", "multiplier", 2, "", "V2 multiplier, 2016 model"),
    p216 = dv(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2, 2016 model V2 cell"),
    m0c = cls("M0", "conservative"), m0n = cls("M0", "nominal"),
    p2m0 = dci(T1, "analysis_model=='M0' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M0 P2, 2020 model V2 cell"),
    cmr = drange(T1, "analysis_model=='M1' & config=='P2' & grepl('^ka_', scenario)", "cmax_ratio", 3, "", "true Cmax ratio, ka-down cells"),
    zero = dext(T1, "analysis_model=='M1' & config=='P2' & grepl('^ka_', scenario)", "pass_pct", max, 0, "%", "P2 in the ka-down cells"),
    ninf = dcount(T1, "analysis_model=='M1' & config=='P2' & !grepl('^ka_', scenario)", "informative boundary cells"),
    aii = mx("M1", "G2_Aii"), b = mx("M1", "G2_B"), aii0 = mx("M0", "G2_Aii"), b0 = mx("M0", "G2_B"),
    lo_aii = lo5("M1", "G2_Aii"), lo_ai = lo5("M1", "G2_Ai"), lo_b = lo5("M1", "G2_B"), lo_aii0 = lo5("M0", "G2_Aii"), lo_ai0 = lo5("M0", "G2_Ai"), lo_b0 = lo5("M0", "G2_B"),
    r0 = s15_cnt_rng("M0", USUAL, "M0 AUC0-inf + Cmax, rules A (i), A (ii) and B: cells above 5% (point), range over configurations"),
    ref = drange(T1, "analysis_model=='M1' & config=='AUCinf_true_only'", "pass_pct", 2, "%", "M1 unbiased reference (true AUC0-inf) boundary range"),
    sd0 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M0' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE AUCinf_true M0"),
    sd1 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M1' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE AUCinf_true M1"),
    pc = drange(PC, "config=='P2' & comparison=='M1 minus M0'", "diff_pp", 2, "", "P2 paired change M1 minus M0 (points)"),
    nv2 = dcount(CB, "analysis_model=='M1' & config %in% c('A_ii','B') & scenario=='V2_up_080' & bias_dir=='away_from_1'", "V2 cells with AUC0-inf GMR bias away from 1, rules A (ii) and B, M1")))))
  deck_end()
}
