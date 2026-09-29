# S18 논거 ③ 보충: 표적 매개 소실 Km(보고서 5.6절의 표적 결합 상수; Michaelis-Menten 근사 상수, 6.5절) 차이는 PK로 보이지 않는다. 시험약 arm의 Km만 ×0.01~×100으로 바꾼 역산 결과(참 AUC0-inf 비, 참 Cmax 비)와
# 다른 기전이 경계(0.80, 1.25)에 닿는 배율의 비교. 결합 유사성은 분석적 유사성 평가가 담당한다(보고서 5.6절, 예상 질의 Q10; 의뢰자 입장).
# 자료: results/oc/inversion_all.csv(범위 끝과 도달 행: 공통 가상 대상자 200,000명), results/oc/inversion_scan_<모델>_<기전>.csv(선별 20,000명 곡선).
# 주의: 0.999~1.059는 두 모델을 합친 범위다(2016 모델 0.999~1.059, 2020 모델 0.999~1.054). 도달 가능한 목표는 1.05뿐.
s18_mult <- function(x) if (x < 1) fnum(x, 2) else fnum(x, 0)
slide_S18 <- function() {
  IA <- "oc/inversion_all.csv"; SC <- function(m, k) sprintf("oc/inversion_scan_%s_%s.csv", m, k)
  M <- c("k2016", "k2020"); MECH <- c("F", "ka", "ke", "Vmax", "V2", "Km")
  deck_slide("S18", tag = "sim")
  L <- DK$txt$S18

  # ---- 제목: Km 배율 범위(config), 참 AUC0-inf 비 범위(범위 끝 + 도달 행, 두 모델; 보고서 km_rng와 같은 계산) ----
  kmr <- dcfg("oc_design.yaml", c("mechanisms", "Km", "range"), "Km multiplier search range, test arm", function(x) sprintf("%s~%s", s18_mult(x[1]), s18_mult(x[2])))
  km_rng <- local({ r <- rows(IA, "mechanism=='Km' & is.finite(end_auc_ratio)"); r2 <- rows(IA, "mechanism=='Km' & reachable==TRUE"); x <- range(c(r$end_auc_ratio, r2$auc_ratio))
    dderived("true AUC0-inf ratio range over Km x0.01 to x100 (range ends and reached rows)", IA, "mechanism=='Km' :: end_auc_ratio, auc_ratio", x, rng_fmt(x[1], x[2], 3)) })
  premise(all(range(rows(IA, "mechanism=='Km' & is.finite(end_multiplier)")$end_multiplier) == c(0.01, 100)), "Km range ends in the result file equal the configured search range")
  deck_kicker(tx("S18.kicker")); deck_title(tx("S18.title", list(m = kmr, r = km_rng)))

  # ---- 그림: Km 배율 대 참 AUC0-inf 비(두 모델, 선별 곡선), 범위 끝(200,000명), 도달 목표(빈 표식). 다른 기전은 표에 둔다 ----
  lim <- unlist(.read("config/trial_design.yaml")$be$limits)
  km <- rbindlist(lapply(M, function(m) rows(SC(m, "Km"))))
  ends <- rows(IA, "mechanism=='Km' & is.finite(end_auc_ratio)")[, .(model, x = end_multiplier, y = end_auc_ratio)] |> unique()
  rch <- rows(IA, "mechanism=='Km' & reachable==TRUE")[, .(model, x = multiplier, y = auc_ratio, target)]
  premise(nrow(rch) == 2 && all(abs(rch$target - 1.05) < 1e-9), "the only reachable Km target is 1.05 in both models")
  premise(all(km$auc_ratio_screen > lim[1] & km$auc_ratio_screen < lim[2]), "the Km curves stay inside the equivalence limits (figure)")
  rlab <- fill(L$fig$reach, list(t = fnum(rch$target[1], 2), m16 = s18_mult(rch[model == "k2016", x]), m20 = s18_mult(rch[model == "k2020", x])))
  km[, mod := factor(model_lab()[model], levels = model_lab())]; ends[, mod := factor(model_lab()[model], levels = model_lab())]; rch[, mod := factor(model_lab()[model], levels = model_lab())]
  setorder(km, model, direction, multiplier)
  FW <- 6.45; FH <- 2.9; XL <- c(0.0065, 155)
  p <- ggplot() +
    annotate("rect", xmin = XL[1], xmax = XL[2], ymin = lim[1], ymax = lim[2], fill = PAL$tint_blue, alpha = 0.6) +
    geom_hline(yintercept = lim, linetype = "22", colour = PAL$ink2, linewidth = 0.5) +
    geom_hline(yintercept = 1, colour = PAL$muted, linewidth = 0.4) +
    geom_line(data = km, aes(multiplier, auc_ratio_screen, colour = mod, linetype = mod, group = interaction(model, direction)), linewidth = 1.1) +
    geom_point(data = ends, aes(x, y, colour = mod, shape = mod), size = 2.8) +
    geom_point(data = rch, aes(x, y, colour = mod), shape = ifelse(rch$model == "k2016", 1, 2), size = 3.4, stroke = 1.1, show.legend = FALSE) +
    annotate("text", x = max(rch$x) * 1.08, y = 1 - 0.018, label = rlab, hjust = 1, vjust = 1, size = 3.7, family = FONT, colour = PAL$ink) +
    annotate("text", x = 0.009, y = lim[2] - 0.012, label = L$fig$band, hjust = 0, vjust = 1, size = 3.7, family = FONT, colour = PAL$ink2) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_linetype_manual(values = unname(MODEL_LT)) + scale_shape_manual(values = unname(MODEL_SHAPE)) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = function(x) paste0("×", formatC(x, format = "fg"))) +
    scale_y_continuous(breaks = c(lim[1], 1, lim[2]), labels = function(x) fnum(x, 2)) +
    coord_cartesian(xlim = XL, ylim = c(lim[1] - 0.05, lim[2] + 0.04), expand = FALSE) +
    labs(x = L$fig$xlab, y = NULL, subtitle = L$fig$ylab) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.minor = element_blank(),
          plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)), plot.margin = margin(6, 16, 6, 6))
  deck_figure(p, "s18_km_inversion", c(GEO$ML, GEO$BODY_TOP, FW, FH), src = c(IA, SC("k2016", "Km"), SC("k2020", "Km")))

  # ---- 오른쪽: 표 제목, 기전별로 참 비 0.80, 1.25에 닿는 배율(두 모델) ----
  XR <- GEO$ML + FW + 0.3; WR <- GEO$W - GEO$MR - XR
  lo <- dcfg("trial_design.yaml", c("be", "limits"), "lower equivalence limit", function(x) fnum(x[1], 2))
  hi <- dcfg("trial_design.yaml", c("be", "limits"), "upper equivalence limit", function(x) fnum(x[2], 2))
  cell <- function(k, tg) {
    w <- sprintf("mechanism=='%s' & abs(target-%s)<1e-9 & reachable==TRUE", k, tg); r <- rows(IA, w)
    if (!nrow(r)) { dcount(IA, w, sprintf("%s: rows reaching target %s (none)", k, tg)); return(L$table$none) }
    premise(nrow(r) == 2 && setequal(r$model, M), sprintf("%s reaches %s in both models", k, tg))
    sprintf("×%s / ×%s", dv(IA, sprintf("%s & model=='k2016'", w), "multiplier", 2, "", sprintf("%s multiplier reaching %s, k2016", k, tg)),
            dv(IA, sprintf("%s & model=='k2020'", w), "multiplier", 2, "", sprintf("%s multiplier reaching %s, k2020", k, tg)))
  }
  df <- data.frame(a = unlist(L$table$mech[MECH]), b = vapply(MECH, cell, "", tg = "0.80"), c = vapply(MECH, cell, "", tg = "1.25"), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S18.table.head", list(lo = lo, hi = hi))
  ntruth <- dcfg("oc_design.yaml", c("estimand", "population", "n_subjects"), "common virtual subjects for the true ratio", function(x) fnum(as.numeric(x), 0, big = TRUE))
  mx1 <- function(v) if (v < 1) sub("(\\.[0-9]*[1-9])0+$", "\\1", fnum(v, 2)) else fnum(v, 0)
  rng_x <- function(x) sprintf("×%s~×%s", mx1(x[1]), mx1(x[2]))
  premise(all(vapply(setdiff(MECH, "Km"), function(k) identical(unlist(.read("config/oc_design.yaml")$mechanisms[[k]]$range), unlist(.read("config/oc_design.yaml")$mechanisms$F$range)), TRUE)),
          "the other mechanisms share one search range (caption)")
  kmx <- dcfg("oc_design.yaml", c("mechanisms", "Km", "range"), "Km multiplier search range (caption)", rng_x)
  othx <- dcfg("oc_design.yaml", c("mechanisms", "F", "range"), "search range of the other mechanisms (F, ka, ke, Vmax, V2)", rng_x)
  CH <- 0.98
  deck_text(tx("S18.table.caption", list(n = ntruth, km = kmx, oth = othx)), c(XR, GEO$BODY_TOP, WR, CH), size = 16, label = "text_caption", color = PAL$ink2, gap_pt = 0)
  TY <- GEO$BODY_TOP + CH + 0.04; TH <- 2.4
  deck_table(df, box = c(XR, TY, WR, TH), widths = c(2.3, 1.35, 1.35), size = 12, highlight = length(MECH), label = "table_mech")

  # ---- 왼쪽 아래: 요점 ----
  nt <- dcount(IA, "mechanism=='Km' & model=='k2016' & direction=='up'", "pre-specified target ratios")
  cm <- local({ r <- rows(IA, "mechanism=='Km' & is.finite(end_cmax_ratio)"); r2 <- rows(IA, "mechanism=='Km' & reachable==TRUE"); x <- range(c(r$end_cmax_ratio, r2$cmax_ratio))
    dderived("true Cmax ratio range over Km x0.01 to x100 (range ends and reached rows)", IA, "mechanism=='Km' :: end_cmax_ratio, cmax_ratio", x, rng_fmt(x[1], x[2], 3)) })
  kmm <- function(m) { r <- rows(IA, sprintf("mechanism=='Km' & model=='%s' & is.finite(end_auc_ratio)", m)); r2 <- rows(IA, sprintf("mechanism=='Km' & model=='%s' & reachable==TRUE", m))
    x <- range(c(r$end_auc_ratio, r2$auc_ratio)); dderived(sprintf("true AUC0-inf ratio range over Km x0.01 to x100, %s", m), IA, sprintf("mechanism=='Km' & model=='%s' :: end_auc_ratio, auc_ratio", m), x, rng_fmt(x[1], x[2], 3)) }
  b <- list(nt = nt, r16 = kmm("k2016"), r20 = kmm("k2020"), tr = drange(IA, "mechanism=='Km'", "target", 2, "", "pre-specified target ratios, range"),
            t = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "target", 2, "", "reachable Km target"),
            m16 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "multiplier", 0, "", "Km multiplier reaching 1.05, k2016"),
            m20 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2020'", "multiplier", 0, "", "Km multiplier reaching 1.05, k2020"),
            cm = cm)
  BY <- GEO$BODY_TOP + FH + 0.1
  deck_bullets(tx("S18.bullets", b), box = c(GEO$ML, BY, FW, GEO$BODY_BOTTOM - BY), size = 16)

  # ---- 오른쪽 아래: 결합 유사성은 분석적 유사성 평가가 담당(의뢰자 입장) ----
  AY <- TY + TH + 0.14
  deck_text(tx("S18.analytic"), c(XR, AY, WR, GEO$BODY_BOTTOM - AY), size = 16, label = "text_analytic", bg = PAL$tint_blue, geom = "roundRect", gap_pt = 6)

  # ---- 노트 ----
  e <- function(m, dir, col, item) dv(IA, sprintf("mechanism=='Km' & model=='%s' & direction=='%s' & target==min(target)", m, dir), col, 3, "", item)
  deck_notes(tx("S18.notes", list(
    ntruth = ntruth,
    nscr = dcfg("oc_design.yaml", c("inversion", "screening_subjects"), "screening subjects", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    wt = f_wt_range(), m = kmr,
    a16d = e("k2016", "down", "end_auc_ratio", "true AUC0-inf ratio at Km x0.01, k2016"), a16u = e("k2016", "up", "end_auc_ratio", "true AUC0-inf ratio at Km x100, k2016"),
    a20d = e("k2020", "down", "end_auc_ratio", "true AUC0-inf ratio at Km x0.01, k2020"), a20u = e("k2020", "up", "end_auc_ratio", "true AUC0-inf ratio at Km x100, k2020"),
    c16u = e("k2016", "up", "end_cmax_ratio", "true Cmax ratio at Km x100, k2016"), c20u = e("k2020", "up", "end_cmax_ratio", "true Cmax ratio at Km x100, k2020"),
    c16d = e("k2016", "down", "end_cmax_ratio", "true Cmax ratio at Km x0.01, k2016"), c20d = e("k2020", "down", "end_cmax_ratio", "true Cmax ratio at Km x0.01, k2020"),
    cm = cm, kmr2 = kmx, oth = othx,
    r16 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "multiplier", 2, "", "Km multiplier reaching 1.05, k2016 (2 decimals)"),
    r20 = dv(IA, "mechanism=='Km' & reachable==TRUE & model=='k2020'", "multiplier", 2, "", "Km multiplier reaching 1.05, k2020 (2 decimals)"),
    t = b$t, olo = dext(IA, "mechanism!='Km' & reachable==TRUE & (abs(target-0.8)<1e-9 | abs(target-1.25)<1e-9)", "multiplier", min, 2, "", "smallest multiplier reaching 0.80 or 1.25, other mechanisms, two models"),
    ohi = dext(IA, "mechanism!='Km' & reachable==TRUE & (abs(target-0.8)<1e-9 | abs(target-1.25)<1e-9)", "multiplier", max, 2, "", "largest multiplier reaching 0.80 or 1.25, other mechanisms, two models"))))
  deck_end()
}
