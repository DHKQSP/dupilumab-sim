# S5 ① 신뢰할 수 있는 AUCinf(지시 §2): 막대 = 기준 미달 비율(adjusted R² ≥ 0.80 세트 (i), ≥ 0.90 세트 (iii); 두 모델; 0.90 강조),
# 오른쪽 = arm당 평가 가능 인원 중 세트 (iii) 충족 인원(중앙값, 5~95백분위; 시험 10,000회). 본문 3줄: λz 산출률, Phoenix·공개 SAP, arm당 인원.
slide_S5 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TRA <- "trialpop/tp_retained_per_arm.csv"
  deck_slide("S5", tag = "litsim")
  f <- list(r2 = f_set("iii", "r2"), fail = drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models"))
  y0 <- core_title(tx("S5.title", f), tx("S5.kicker"))
  L <- DK$txt$S5$fig; ML <- DK$txt$common$models_short

  # ---- 왼쪽: 미달 비율 막대 ----
  d <- copy(rows(TPF, "set %in% c('i','iii')"))[, .(pk_model, set, fail_pct, fail_lo, fail_hi)]
  premise(nrow(d) == 4, "sets (i) and (iii), two models")
  premise(all(rows(TPF, "set=='iii'")[, est_fail_pct > lz_pct]), "set (iii): most failures have an estimable lambda-z (text: mostly below the fit criterion)")
  d[, x := c(i = 0, iii = 2.6)[set] + c(k2016 = 1, k2020 = 2)[pk_model]]
  d[, lab := paste0(fnum(fail_pct, 1), "%")]
  r2i <- fnum(.read("config/prereg_20260926.yaml")$section4$criteria_sets$i$adj_r2_min, 2); r2iii <- fnum(.read("config/prereg_20260926.yaml")$section4$criteria_sets$iii$adj_r2_min, 2)
  grp <- data.table(x = c(1.5, 4.1), lab = c(fill(L$g80, list(v = r2i)), fill(L$g90, list(v = r2iii))), face = c("plain", "bold"))
  ymax <- max(d$fail_hi) * 1.45
  p1 <- ggplot(d) +
    geom_col(aes(x, fail_pct, fill = set), width = 0.78) +
    geom_errorbar(aes(x, ymin = fail_lo, ymax = fail_hi), width = 0.18, colour = PAL$ink2, linewidth = 0.5) +
    geom_text(aes(x, fail_hi, label = lab, fontface = ifelse(set == "iii", "bold", "plain")), vjust = -0.45, size = PT(18), family = FONT, colour = PAL$ink) +
    geom_text(data = grp, aes(x, ymax * 0.97, label = lab, fontface = face), vjust = 1, size = PT(16), family = FONT, colour = PAL$ink, lineheight = 0.95) +
    scale_fill_manual(values = c(i = ORANGE_LIGHT, iii = PAL$orange), guide = "none") +
    scale_x_continuous(breaks = d$x, labels = unlist(ML[d$pk_model]), expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, 40, 10), labels = function(v) paste0(v, "%"), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") + labs(x = NULL, y = NULL, subtitle = L$ylab) + theme_core(16) +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 14, colour = PAL$ink2, margin = margin(2, 0, 0, 0)))
  # ---- 오른쪽: arm당 신뢰할 수 있는 AUCinf 인원 ----
  na <- as.numeric(f_n_arm())
  r <- copy(rows(TRA, "set=='iii'"))[, .(pk_model, med = retained_median, lo = retained_p05, hi = retained_p95)]
  premise(nrow(r) == 2 && all(r$hi <= na), "set (iii) retained per arm: two models, at most the evaluable count")
  r[, y := c(k2016 = 2, k2020 = 1)[pk_model]][, lab := vapply(seq_len(.N), function(i) fill(L$ret, list(m = fnum(med[i], 0), lo = fnum(lo[i], 0), hi = fnum(hi[i], 0))), "")]
  p2 <- ggplot(r) +
    geom_rect(aes(xmin = 0, xmax = na, ymin = y - 0.32, ymax = y + 0.32), fill = PAL$tint_grey) +
    geom_rect(aes(xmin = 0, xmax = med, ymin = y - 0.32, ymax = y + 0.32), fill = PAL$orange) +
    geom_errorbarh(aes(y = y, xmin = lo, xmax = hi), height = 0.2, colour = PAL$ink, linewidth = 0.6) +
    geom_text(aes(x = 2, y = y + 0.52, label = unlist(ML[pk_model])), hjust = 0, size = PT(15), family = FONT, colour = PAL$ink2) +
    geom_text(aes(x = med / 2, y = y, label = lab), size = PT(16), family = FONT, colour = "white", fontface = "bold") +
    annotate("text", x = na, y = 2.55, label = fill(L$total, list(n = fnum(na, 0))), hjust = 1, size = PT(15), family = FONT, colour = PAL$ink2) +
    scale_x_continuous(limits = c(0, na), breaks = c(0, 50, 100), expand = expansion(add = c(0, 2))) +
    scale_y_continuous(limits = c(0.55, 2.75), breaks = NULL, expand = expansion(mult = 0)) +
    labs(x = L$ret_x, y = NULL, subtitle = fill(L$ret_title, list(n = fnum(na, 0)))) + theme_core(16) + theme(panel.grid.major.y = element_blank())
  p <- patchwork::wrap_plots(p1, p2, widths = c(1.25, 1))

  lzr <- { x <- rows(TPF, "set=='iii'"); v <- 100 - x$lz_pct
    dderived("share with an estimable lambda-z (100 - lz_pct), two models", TPF, "set=='iii' :: 100 - lz_pct, range over models", range(v), rng_fmt(min(v), max(v), 1, "%")) }
  body <- tx("S5.body", list(lz = lzr, sap = core_sap("NCT04117607", "at least", "public SAP NCT04117607: adjusted R-squared at least"),
                             n = f_n_arm(), ret = drange(TRA, "set=='iii'", "retained_median", 0, "", "median evaluable subjects per arm with a reliable AUCinf, set (iii), two models")))
  premise(core_sap("NCT04441905", "at least", "public SAP NCT04441905: adjusted R-squared at least") == f$r2, "second public SAP uses the same threshold")
  cap <- tx("S5.caption", list(r2 = f$r2, ex = f_set("iii", "extrap")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "s5_reliable", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(TPF, TRA))

  deck_notes(tx("S5.notes", list(
    f16i = dv(TPF, "pk_model=='k2016' & set=='i'", "fail_pct", 1, "%", "set (i) failing share, 2016"), f20i = dv(TPF, "pk_model=='k2020' & set=='i'", "fail_pct", 1, "%", "set (i) failing share, 2020"),
    f16 = dv(TPF, "pk_model=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing share, 2016"), f20 = dv(TPF, "pk_model=='k2020' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing share, 2020"),
    e16 = dv(TPF, "pk_model=='k2016' & set=='iii'", "est_fail_pct", 1, "%", "set (iii) failing with an estimable lambda-z, 2016"), e20 = dv(TPF, "pk_model=='k2020' & set=='iii'", "est_fail_pct", 1, "%", "set (iii) failing with an estimable lambda-z, 2020"),
    r2i = f_set("i", "r2"), r2 = f$r2, ex = f_set("iii", "extrap"),
    s3 = core_sap("NCT04700163", "above", "public SAP NCT04700163: adjusted R-squared above"),
    ri = drange(TRA, "set=='i'", "retained_median", 0, "", "median evaluable subjects per arm meeting set (i), two models"),
    lo = dspan(TRA, "set=='iii'", "retained_p05", "retained_p95", 0, "", "5th to 95th percentile subjects per arm meeting set (iii), two models"),
    wt = f_wt_range(), nsub = dint(TPF, "pk_model=='k2016' & set=='iii'", "n", "virtual subjects per model"), ntr = dint(TRA, "pk_model=='k2016' & set=='iii'", "n_trials", "simulated trials per model (retained per arm)"))))
  deck_end()
}
