# 별첨 A8 Km 비가시성: 시험약 arm의 Km(표적 매개 소실의 Michaelis-Menten 반포화 농도)만 ×0.01~×100으로 바꾼 역산 결과.
#  왼쪽 그림: Km 배율 대 참 AUCinf 비(두 모델: 선 모양·표식; 선별 20,000명 곡선 results/oc/inversion_scan_<모델>_Km.csv),
#             도달 가능한 유일한 목표(1.05)의 빈 표식(200,000명, results/oc/inversion_all.csv), 동등성 범위 음영.
#  오른쪽 표: 기전별로 참 AUCinf 비가 동등 경계(0.80, 1.25)에 닿는 배율(두 모델; inversion_all.csv).
#  제목의 범위(0.999~1.059)는 두 모델을 합친 범위 끝 값과 도달 행(보고서 km_rng와 같은 계산). 결과보고 덱 S18(slides/s18_binding.R)의 자료 처리를 따른다.
a8_mult <- function(x) if (x < 1) fnum(x, 2) else fnum(x, 0)
a8_rng <- function(IA, where, item, cols = c("end_auc_ratio", "auc_ratio")) {
  r <- rows(IA, sprintf("%s & is.finite(%s)", where, cols[1])); r2 <- rows(IA, sprintf("%s & reachable==TRUE", where)); x <- range(c(r[[cols[1]]], r2[[cols[2]]]))
  dderived(item, IA, sprintf("%s :: %s, %s", where, cols[1], cols[2]), x, rng_fmt(x[1], x[2], 3))
}
slide_A8 <- function() {
  IA <- "oc/inversion_all.csv"; SC <- function(m) sprintf("oc/inversion_scan_%s_Km.csv", m)
  M <- c("k2016", "k2020"); MECH <- c("F", "ka", "ke", "Vmax", "V2", "Km")
  deck_slide("A8", tag = "sim")
  L <- DK$txt$A8; ML <- DK$txt$common$models_short

  # ---- 제목: Km 배율 범위(config), 참 AUCinf 비 범위(두 모델) ----
  kmr <- dcfg("oc_design.yaml", c("mechanisms", "Km", "range"), "Km multiplier search range, test arm", function(x) sprintf("%s~%s", a8_mult(x[1]), a8_mult(x[2])))
  premise(all(range(rows(IA, "mechanism=='Km' & is.finite(end_multiplier)")$end_multiplier) == unlist(.read("config/oc_design.yaml")$mechanisms$Km$range)),
          "Km range ends in the result file equal the configured search range")
  km_rng <- a8_rng(IA, "mechanism=='Km'", "true AUCinf ratio range over the Km search range (range ends and reached rows), two models")
  y0 <- core_title(tx("A8.title", list(m = kmr, r = km_rng)), tx("A8.kicker"))

  # ---- 캡션·본문(아래) ----
  lim <- unlist(.read("config/trial_design.yaml")$be$limits)
  cm <- a8_rng(IA, "mechanism=='Km'", "true Cmax ratio range over the Km search range (range ends and reached rows), two models", c("end_cmax_ratio", "cmax_ratio"))
  rch <- rows(IA, "mechanism=='Km' & reachable==TRUE")
  premise(nrow(rch) == 2 && all(abs(rch$target - 1.05) < 1e-9) && setequal(rch$model, M), "the only reachable Km target is 1.05 in both models")
  premise(nrow(rows(IA, sprintf("mechanism=='Km' & reachable==TRUE & (abs(target-%s)<1e-9 | abs(target-%s)<1e-9)", lim[1], lim[2]))) == 0, "Km reaches neither equivalence limit")
  vk <- c(.read("config/params_variability.yaml")$iiv_omega2$Km$omega2, .read("config/params_k2020_model1.yaml")$iiv$sd$Km$omega2)
  ntruth <- dcfg("oc_design.yaml", c("estimand", "population", "n_subjects"), "common virtual subjects for the true ratio", function(x) fnum(as.numeric(x), 0, big = TRUE))
  capy <- core_caption(tx("A8.caption", list(n = ntruth, lo = dcfg("trial_design.yaml", c("be", "limits"), "equivalence limits", function(x) sprintf("%s~%s", fnum(x[1], 2), fnum(x[2], 2))))),
                       GEO$BODY_BOTTOM, size = 14)
  body <- tx("A8.body", list(cm = cm, t = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "target", 2, "", "reachable Km target"),
                             m16 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "multiplier", 0, "", "Km multiplier reaching 1.05, k2016"),
                             m20 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2020'", "multiplier", 0, "", "Km multiplier reaching 1.05, k2020")))

  # ---- 오른쪽 표: 기전별 동등 경계 도달 배율 ----
  WT <- c(2.5, 1.55, 1.55); WR <- sum(WT); XR <- GEO$ML + GEO$CW - WR; FW <- XR - 0.3 - GEO$ML
  lo <- dcfg("trial_design.yaml", c("be", "limits"), "lower equivalence limit", function(x) fnum(x[1], 2))
  hi <- dcfg("trial_design.yaml", c("be", "limits"), "upper equivalence limit", function(x) fnum(x[2], 2))
  cell <- function(k, tg) {
    w <- sprintf("mechanism=='%s' & abs(target-%s)<1e-9 & reachable==TRUE", k, tg); r <- rows(IA, w)
    if (!nrow(r)) { dcount(IA, w, sprintf("%s: rows reaching target %s (none)", k, tg)); return(L$table$none) }
    premise(nrow(r) == 2 && setequal(r$model, M), sprintf("%s reaches %s in both models", k, tg))
    sprintf("×%s / ×%s", dv(IA, sprintf("%s & model=='k2016'", w), "multiplier", 2, "", sprintf("%s multiplier reaching %s, k2016", k, tg)),
            dv(IA, sprintf("%s & model=='k2020'", w), "multiplier", 2, "", sprintf("%s multiplier reaching %s, k2020", k, tg)))
  }
  df <- data.frame(a = unlist(L$table$mech[MECH]), b = vapply(MECH, cell, "", tg = fnum(lim[1], 2)), c = vapply(MECH, cell, "", tg = fnum(lim[2], 2)), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("A8.table.head", list(lo = lo, hi = hi))
  TH <- 2.75
  deck_table(df, box = c(XR, y0 + 0.05, WR, TH), widths = WT, size = 14, highlight = length(MECH), label = "table_mech")
  bh <- core_body_h(body, WR); by <- y0 + 0.05 + TH + 0.12; premise(by + bh < capy - 0.04, "A8: body text fits between the table and the caption")
  deck_text(body, c(XR, by, WR, bh), size = SZ$body, label = "body", gap_pt = 4)

  # ---- 왼쪽 그림: Km 배율 대 참 AUCinf 비 ----
  km <- rbindlist(lapply(M, function(m) rows(SC(m))))
  premise(all(km$auc_ratio_screen > lim[1] & km$auc_ratio_screen < lim[2]), "the Km curves stay inside the equivalence limits (figure)")
  km[, mod := factor(unlist(ML[model]), levels = unlist(ML))]; setorder(km, model, direction, multiplier)
  rp <- rch[, .(model, x = multiplier, y = auc_ratio)][, mod := factor(unlist(ML[model]), levels = unlist(ML))]
  F <- L$fig
  rlab <- fill(F$reach, list(t = fnum(rch$target[1], 2), m16 = a8_mult(rch[model == "k2016", multiplier]), m20 = a8_mult(rch[model == "k2020", multiplier])))
  XL <- c(0.0065, 155)
  p <- ggplot() +
    annotate("rect", xmin = XL[1], xmax = XL[2], ymin = lim[1], ymax = lim[2], fill = PAL$tint_blue, alpha = 0.7) +
    geom_hline(yintercept = lim, linetype = "22", colour = PAL$ink2, linewidth = 0.5) +
    geom_hline(yintercept = 1, colour = PAL$muted, linewidth = 0.4) +
    geom_line(data = km, aes(multiplier, auc_ratio_screen, linetype = mod, group = interaction(model, direction)), colour = PAL$orange, linewidth = 1.2) +
    geom_point(data = rp, aes(x, y, shape = mod), colour = PAL$ink, size = 3.6, stroke = 1.1) +
    annotate("segment", x = 12, y = 1.135, xend = min(rp$x) * 0.9, yend = max(rp$y) + 0.012, colour = PAL$ink2, linewidth = 0.4) +
    annotate("text", x = 11, y = 1.135, label = rlab, hjust = 1, vjust = 0.5, size = PT(14), family = FONT, colour = PAL$ink) +
    annotate("text", x = 0.009, y = lim[2] - 0.015, label = F$band, hjust = 0, vjust = 1, size = PT(14), family = FONT, colour = PAL$ink2) +
    scale_linetype_manual(values = setNames(unname(CORE_MODEL_LT[M]), unlist(ML[M])), name = NULL) +
    scale_shape_manual(values = setNames(c(1, 2), unlist(ML[M])), name = NULL) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = function(x) paste0("×", formatC(x, format = "fg"))) +
    scale_y_continuous(breaks = c(lim[1], 1, lim[2]), labels = function(x) fnum(x, 2)) +
    coord_cartesian(xlim = XL, ylim = c(lim[1] - 0.05, lim[2] + 0.04), expand = FALSE) +
    labs(x = F$xlab, y = NULL, subtitle = F$ylab) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.minor = element_blank(),
          plot.subtitle = element_text(colour = PAL$ink2, size = 14, margin = margin(0, 0, 2, 0)), plot.margin = margin(4, 14, 4, 4))
  deck_figure(p, "a8_km_inversion", c(GEO$ML, y0, FW, capy - 0.1 - y0), src = c(IA, SC("k2016"), SC("k2020")))

  # ---- 노트 ----
  premise(length(vk) == 2 && all(vk == 0), "Km has no between-subject variability in either model (notes)")
  e <- function(m, dir, col, item) dv(IA, sprintf("mechanism=='Km' & model=='%s' & direction=='%s' & target==min(target)", m, dir), col, 3, "", item)
  mx1 <- function(v) if (v < 1) sub("(\\.[0-9]*[1-9])0+$", "\\1", fnum(v, 2)) else fnum(v, 0)
  rng_x <- function(x) sprintf("×%s~×%s", mx1(x[1]), mx1(x[2]))
  premise(all(vapply(setdiff(MECH, "Km"), function(k) identical(unlist(.read("config/oc_design.yaml")$mechanisms[[k]]$range), unlist(.read("config/oc_design.yaml")$mechanisms$F$range)), TRUE)),
          "the other mechanisms share one search range (notes)")
  deck_notes(tx("A8.notes", list(
    ntruth = ntruth, nscr = dcfg("oc_design.yaml", c("inversion", "screening_subjects"), "screening subjects", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    wt = f_wt_range(), m = kmr,
    r16 = a8_rng(IA, "mechanism=='Km' & model=='k2016'", "true AUCinf ratio range over the Km search range, k2016"),
    r20 = a8_rng(IA, "mechanism=='Km' & model=='k2020'", "true AUCinf ratio range over the Km search range, k2020"),
    a16d = e("k2016", "down", "end_auc_ratio", "true AUCinf ratio at Km x0.01, k2016"), a16u = e("k2016", "up", "end_auc_ratio", "true AUCinf ratio at Km x100, k2016"),
    a20d = e("k2020", "down", "end_auc_ratio", "true AUCinf ratio at Km x0.01, k2020"), a20u = e("k2020", "up", "end_auc_ratio", "true AUCinf ratio at Km x100, k2020"),
    c16u = e("k2016", "up", "end_cmax_ratio", "true Cmax ratio at Km x100, k2016"), c20u = e("k2020", "up", "end_cmax_ratio", "true Cmax ratio at Km x100, k2020"),
    nt = dcount(IA, "mechanism=='Km' & model=='k2016' & direction=='up'", "pre-specified target ratios"),
    tr = drange(IA, "mechanism=='Km'", "target", 2, "", "pre-specified target ratios, range"),
    q16 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "multiplier", 2, "", "Km multiplier reaching 1.05, k2016 (2 decimals)"),
    q20 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2020'", "multiplier", 2, "", "Km multiplier reaching 1.05, k2020 (2 decimals)"),
    kmx = dcfg("oc_design.yaml", c("mechanisms", "Km", "range"), "Km multiplier search range (notes)", rng_x),
    othx = dcfg("oc_design.yaml", c("mechanisms", "F", "range"), "search range of the other mechanisms (F, ka, ke, Vmax, V2)", rng_x),
    olo = dext(IA, "mechanism!='Km' & reachable==TRUE & (abs(target-0.8)<1e-9 | abs(target-1.25)<1e-9)", "multiplier", min, 2, "", "smallest multiplier reaching 0.80 or 1.25, other mechanisms, two models"),
    ohi = dext(IA, "mechanism!='Km' & reachable==TRUE & (abs(target-0.8)<1e-9 | abs(target-1.25)<1e-9)", "multiplier", max, 2, "", "largest multiplier reaching 0.80 or 1.25, other mechanisms, two models"),
    kv = dcfg("params_typical.yaml", c("theta", "Km", "value"), "Km typical value, 2016 model (mg/L)", function(x) fnum(as.numeric(x), 2)),
    kv20 = dcfg("params_k2020_model1.yaml", c("theta", "Km", "value"), "Km typical value, 2020 model (mg/L)", function(x) fnum(as.numeric(x), 2)),
    lq = f_lloq())))
  deck_end()
}
