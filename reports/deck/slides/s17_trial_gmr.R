# S17 논거 ③-2 시험 수준: 시나리오별 시험 GMR의 기하평균(exp(mean_log_gmr))을 참 AUC0-inf 비와 비교. 평가변수: 참 AUC0-inf(모델 적분, 잔차 없음),
# AUC0-last, NCA AUC0-inf 규칙 A (i)·A (ii)·B·C (i). 분석 모형 M1(주), M0 병기(노트). 경계 16칸 + 제품 시나리오(S00 동일 제품, F097 생체이용률 감소).
# 자료: results/oc_models/sd_se_models.csv(시험 모집단 재생성 시험, 칸당 10,000회; 2020 모델 V2 증가는 20,000회).
# 끝점 코드: AUCinf_A = 규칙 A 세트 (ii), AUCinf_Ai = 규칙 A 세트 (i), AUCinf_Ci = 규칙 C 세트 (i)(R/oc.R, R/oc_models.R).
# 미산출: AUC0-last 판정과 이상적 AUC0-inf(모델 적분값) 판정의 시험별 일치율은 프로젝트에서 계산하지 않았다(있는 일치율은 AUC0-last 대 NCA 규칙 A (ii)뿐).
s17_maxabs <- function(rel, where, col, d, unit, item) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- max(abs(r[[col]]))
  dderived(item, rel, sprintf("max(abs(%s)) over %d rows [%s]", col, nrow(r), where), x, paste0(fnum(x, d), unit))
}
s17_signed <- function(p) if (!startsWith(p, "-")) paste0("+", p) else p
slide_S17 <- function() {
  SSf <- "oc_models/sd_se_models.csv"
  EP <- c("AUCinf_true", "AUClast", "AUCinf_Ai", "AUCinf_A", "AUCinf_B", "AUCinf_Ci")
  BND <- "!scenario %in% c('S00','F097')"
  wb <- function(ep, am = "M1") sprintf("analysis_model=='%s' & endpoint=='%s' & %s", am, ep, BND)
  deck_slide("S17", tag = "sim")
  L <- DK$txt$S17
  premise(nrow(rows(SSf, wb("AUClast"))) == 16 && all(vapply(EP, function(e) nrow(rows(SSf, wb(e))) == 16, TRUE)), "16 boundary cells per endpoint under M1")

  # ---- 제목 ----
  nb <- dcount(SSf, wb("AUClast"), "boundary cells (both models), M1")
  f <- list(a = s17_maxabs(SSf, wb("AUClast"), "bias_pct", 2, "%", "AUC0-last: largest absolute bias of the geometric mean trial GMR vs the true AUC0-inf ratio, 16 boundary cells, M1"),
            b = s17_maxabs(SSf, sprintf("analysis_model=='M1' & endpoint %%in%% c('AUCinf_Ai','AUCinf_A','AUCinf_B','AUCinf_Ci') & %s", BND), "bias_pct", 2, "%",
                           "NCA AUC0-inf rules A (i), A (ii), B, C (i): largest absolute bias, 16 boundary cells, M1"))
  deck_kicker(tx("S17.kicker")); deck_title(tx("S17.title", f))

  # ---- 왼쪽: 경계 16칸의 편향(1 쪽 = 동등 판정이 쉬워지는 방향을 양수로), 평가변수별, 두 모델 ----
  d <- copy(rows(SSf, sprintf("analysis_model=='M1' & endpoint %%in%% c(%s) & %s", paste(sprintf("'%s'", EP), collapse = ","), BND)))
  premise(all(d$true_ratio != 1), "boundary cells have a true ratio different from 1")
  d[, tw := ifelse(true_ratio < 1, bias_pct, -bias_pct)]
  EL <- unlist(L$ep[EP]); d[, ep := factor(EL[endpoint], levels = rev(EL))]
  d[, y := as.numeric(ep) + ifelse(pk_model == "k2016", 0.13, -0.13)]
  d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  xr <- range(d$tw); xl <- c(floor(xr[1]) - 0.5, ceiling(xr[2]) + 0.5)
  FW <- 6.55; FH <- 3.15
  p <- ggplot(d, aes(x = tw, y = y, colour = model, shape = model)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = PAL$tint_orange, alpha = 0.7) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_point(size = 2.7, alpha = 0.9, stroke = 0.8) +
    annotate("text", x = xl[2], y = length(EL) + 0.62, label = L$fig$right, hjust = 1, vjust = 1, size = 3.9, family = FONT, colour = PAL$orange, fontface = "bold") +
    annotate("text", x = xl[1], y = length(EL) + 0.62, label = L$fig$left, hjust = 0, vjust = 1, size = 3.9, family = FONT, colour = PAL$ink2) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = c(16, 17)) +
    scale_y_continuous(breaks = seq_along(levels(d$ep)), labels = levels(d$ep), limits = c(0.5, length(EL) + 0.65), expand = expansion(0)) +
    scale_x_continuous(limits = xl, breaks = seq(ceiling(xl[1] / 2) * 2, floor(xl[2] / 2) * 2, by = 2)) +
    labs(x = L$fig$xlab, y = NULL) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          panel.grid.major.y = element_line(colour = PAL$grid, linewidth = 0.35), axis.text.y = element_text(colour = PAL$ink, size = 12))
  deck_figure(p, "s17_bias_strip", c(GEO$ML, GEO$BODY_TOP, FW, FH), src = SSf)

  # ---- 왼쪽 아래: 요점 ----
  v2 <- dv(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2020' & scenario=='V2_up_080'", "bias_pct", 2, "%", "AUC0-last bias, 2020 model V2 up (true ratio 0.80), M1")
  aw <- dcount(SSf, sprintf("%s & bias_dir=='away_from_1'", wb("AUClast")), "AUC0-last cells biased away from 1, M1")
  premise(as.integer(aw) == as.integer(nb) - 1L && nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1' & pk_model=='k2020' & scenario=='V2_up_080'", wb("AUClast")))) == 1,
          "AUC0-last is biased away from 1 in all boundary cells but the 2020 model V2 cell (M1)")
  premise(all(vapply(c("AUCinf_Ai", "AUCinf_A", "AUCinf_B"), function(e) nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1'", wb(e)))) > as.integer(nb) / 2, TRUE)),
          "NCA rules A (i), A (ii) and B are biased toward 1 in most boundary cells (M1)")
  deck_bullets(tx("S17.bullets", list(aw = aw, n = nb, v2 = s17_signed(v2))), box = c(GEO$ML, GEO$BODY_TOP + FH + 0.08, FW, GEO$BODY_BOTTOM - GEO$BODY_TOP - FH - 0.08), size = 16)

  # ---- 오른쪽: 시나리오별 기하평균비 표(M1) ----
  XR <- GEO$ML + FW + 0.3; WR <- GEO$W - GEO$MR - XR
  gm <- function(ep, m, sc) {
    r <- row1(SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & pk_model=='%s' & scenario=='%s'", ep, m, sc)); x <- exp(r$mean_log_gmr)
    dderived(sprintf("geometric mean of trial GMRs, %s, %s, %s, M1", ep, m, sc), SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & pk_model=='%s' & scenario=='%s' :: exp(mean_log_gmr)", ep, m, sc), x, fnum(x, 4))
  }
  tr <- function(m, sc) dv(SSf, sprintf("analysis_model=='M1' & endpoint=='AUClast' & pk_model=='%s' & scenario=='%s'", m, sc), "true_ratio", 4, "", sprintf("true AUC0-inf ratio, %s, %s", m, sc))
  SC <- list(c("k2016", "Vmax_up_080"), c("k2020", "ka_down_080"), c("k2016", "F097"))
  tw_n <- function(ep) dcount(SSf, sprintf("%s & bias_dir=='toward_1'", wb(ep)), sprintf("cells biased toward 1 (Monte Carlo interval excludes 0), %s, M1", ep))
  df <- data.frame(a = EL, b = vapply(EP, tw_n, ""),
                   c = vapply(EP, gm, "", m = SC[[1]][1], sc = SC[[1]][2]), d = vapply(EP, gm, "", m = SC[[2]][1], sc = SC[[2]][2]), e = vapply(EP, gm, "", m = SC[[3]][1], sc = SC[[3]][2]),
                   stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S17.table.head", list(n = nb, t1 = tr(SC[[1]][1], SC[[1]][2]), t2 = tr(SC[[2]][1], SC[[2]][2]), t3 = tr(SC[[3]][1], SC[[3]][2])))
  TH <- 3.3
  deck_table(df, box = c(XR, GEO$BODY_TOP, WR, TH), widths = c(1.95, 0.95, 1.05, 1.05, 1.05), size = 12, highlight = 2, label = "table_gm")

  # ---- 오른쪽 아래: 미산출(판정 일치율) ----
  deck_text(tx("S17.na"), c(XR, GEO$BODY_TOP + TH + 0.12, WR, GEO$BODY_BOTTOM - GEO$BODY_TOP - TH - 0.12), size = 16, label = "text_na", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 4)

  # ---- 노트 ----
  rng <- function(ep, am = "M1") { p <- drange(SSf, wb(ep, am), "bias_pct", 2, "%", sprintf("bias range, 16 boundary cells, %s, %s", ep, am)); p }
  gs00 <- function(ep) { r <- rows(SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & scenario=='S00'", ep)); x <- range(exp(r$mean_log_gmr))
    dderived(sprintf("identical products: geometric mean of trial GMRs, two models, %s, M1", ep), SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & scenario=='S00' :: range(exp(mean_log_gmr))", ep), x, rng_fmt(x[1], x[2], 4)) }
  deck_notes(tx("S17.notes", list(
    ntr = dint(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2016' & scenario=='F_down_080'", "n_trials", "trials per boundary cell"),
    ntr_v2 = dint(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials, 2020 model V2 cell"),
    n_arm = f_n_arm(), wt = f_wt_range(),
    r_true = rng("AUCinf_true"), r_last = rng("AUClast"), r_last0 = rng("AUClast", "M0"), r_ai = rng("AUCinf_Ai"), r_aii = rng("AUCinf_A"), r_b = rng("AUCinf_B"), r_ci = rng("AUCinf_Ci"),
    s_true = gs00("AUCinf_true"), s_last = gs00("AUClast"), s_ai = gs00("AUCinf_Ai"), s_b = gs00("AUCinf_B"),
    v2 = s17_signed(v2), a = f$a,
    g_last = gm("AUClast", "k2016", "Vmax_up_080"), g_ai = gm("AUCinf_Ai", "k2016", "Vmax_up_080"), g_b = gm("AUCinf_B", "k2016", "Vmax_up_080"),
    t_v = tr("k2016", "Vmax_up_080"), n = nb,
    agr = drange("rationale/pillar2_products_B0.csv", "TRUE", "agree", 1, "%", "decision agreement AUC0-last vs NCA AUC0-inf rule A (ii), product scenarios, two models"))))
  deck_end()
}
