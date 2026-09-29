# S4 전제: 말단 절벽(지시 §2). 대표 대상자 1명(2016 모델, 60~90 kg, 정량한계 도달일 중앙값 대상자; results/deck_inputs, scripts/62에서 커밋된 요약과 대조)의
# 참 농도 반로그 곡선, 현행 채혈일(B0) 점(정량값은 채운 점, 정량한계 미만은 빈 점), 정량한계 선, 절벽 음영(순간 반감기 1일 미만 ~ 정량한계 도달).
# 제목 수치는 두 모델의 중앙값(cliff_summary.csv): 절벽 시작 농도(약 0.9 mg/L), 절벽 길이(1.4일).
# 인용: FDA BLA 761055 임상약리 리뷰(config/literature_core_deck.yaml, 검색 발췌·쪽수 미대조)와 Kovalenko 2020(literature_qualitative.csv).
slide_S4 <- function() {
  CS <- "cliff/cliff_summary.csv"; RP <- "deck_inputs/rep_profiles.csv"; RS <- "deck_inputs/rep_subjects.csv"; CP <- "cliff/cliff_points.csv"; LC <- "config/literature_core_deck.yaml"
  deck_slide("S4", tag = "litsim")
  f <- list(c = drange(CS, "weight=='base'", "c_start1_median", 1, "", "median concentration at cliff start (mg/L), two models"),
            len = drange(CS, "weight=='base'", "len1_median", 1, "", "median cliff length (days), two models"))
  y0 <- core_title(tx("S4.title", f), tx("S4.kicker"))

  # ---- 그림 ----
  L <- DK$txt$S4$fig
  subj <- rows(RS, "percentile==50"); premise(nrow(subj) == 1, "one median representative subject")
  pc <- dderived("representative subject: percentile of the LLOQ day", RS, "percentile==50 :: percentile", subj$percentile, fnum(subj$percentile, 0))
  prof <- rows(RP, sprintf("id==%d", subj$id))[, day := time + 1]
  lloq_v <- .read("config/assay.yaml")$lloq_mg_L$value; lloq <- f_lloq()
  b0 <- unlist(.read("config/trial_design.yaml")$schedules$B0$days)
  full <- as.data.table(expand.grid(time = b0))[, C := approx(prof$time, prof$C, time, rule = 2)$y][, day := time + 1]
  premise(all(round(b0, 2) <= max(prof$time) | TRUE), "B0 days inside the profile grid or after it")
  full[time > max(prof$time), C := NA_real_]                        # 격자 밖(0.01 mg/L 미만으로 잘린 뒤)은 정량한계 미만
  full[, quant := !is.na(C) & C >= lloq_v]
  premise(sum(full$quant) >= 8 && sum(!full$quant) >= 1, "the median subject has quantifiable and BLQ planned samples")
  ymin <- 0.01; full[quant == FALSE, Cp := ymin * 1.25][quant == TRUE, Cp := C]
  x0 <- subj$start1 + 1; x1 <- subj$t_lloq + 1; c0 <- approx(prof$time, prof$C, subj$start1)$y
  brk <- c(1, sort(unique(round(b0[b0 >= 7] + 1))))
  lab_ <- function(x, y, s, hjust, col = PAL$ink2, face = "plain", vjust = 0.5)
    annotate("label", x = x, y = y, label = s, hjust = hjust, vjust = vjust, size = PT(15), family = FONT, colour = col, fontface = face, lineheight = 0.95,
             fill = "white", label.size = 0, label.padding = grid::unit(0.10, "lines"), label.r = grid::unit(0, "lines"))
  p <- ggplot() +
    annotate("rect", xmin = x0, xmax = x1, ymin = ymin, ymax = 100, fill = PAL$orange, alpha = 0.22) +
    geom_hline(yintercept = lloq_v, linetype = "22", colour = PAL$ink2, linewidth = 0.6) +
    geom_line(data = prof, aes(day, C), colour = PAL$ink, linewidth = 1.0) +
    geom_point(data = full[quant == TRUE], aes(day, Cp), shape = 21, fill = PAL$blue, colour = "white", size = 3.4, stroke = 0.7) +
    geom_point(data = full[quant == FALSE], aes(day, Cp), shape = 21, fill = "white", colour = PAL$blue, size = 3.2, stroke = 1.1) +
    lab_(x1 + 0.8, 20, fill(L$cliff, list(len = fnum(subj$len1, 2))), hjust = 0, col = PAL$orange, face = "bold") +
    lab_(x1 + 0.8, 7, fill(L$start, list(c = fnum(c0, 2), d1 = fnum(.read("config/oc_design.yaml")$cliff$definition_days[1], 0))), hjust = 0, col = PAL$orange) +
    lab_(max(b0) + 1 + 2.5, lloq_v, fill(L$lloq, list(v = lloq)), hjust = 1, vjust = -0.15) +
    lab_(max(b0) + 1 + 2.5, ymin * 1.25, L$blq, hjust = 1, vjust = -0.5, col = PAL$blue) +
    lab_(3, 0.03, L$legend, hjust = 0, col = PAL$ink2) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = function(x) formatC(x, format = "fg"), limits = c(ymin, 100), expand = expansion(mult = 0)) +
    scale_x_continuous(breaks = brk, limits = c(0, max(b0) + 1 + 2.5), expand = expansion(mult = 0)) +
    labs(x = L$xlab, y = NULL, subtitle = L$ylab) + theme_core(16) +
    theme(plot.subtitle = element_text(margin = margin(0, 0, 4, 0)), panel.grid.minor = element_blank())
  body <- tx("S4.body", list(d1 = dcfg("oc_design.yaml", c("cliff", "definition_days"), "cliff start: instantaneous half-life threshold (days)", function(x) fnum(x[1], 0)),
                             int = dint(CP, "model=='k2016' & weight=='base' & timing=='nominal' & definition_day==1 & schedule=='current'", "min_interval_day", "sampling interval of the current schedule in the cliff window (days)"),
                             ge1 = dv(CP, "model=='k2016' & weight=='base' & timing=='nominal' & definition_day==1 & schedule=='current'", "pct_ge1", 0, "%", "share with one or more current samples on the cliff, 2016"),
                             ge2 = dv(CP, "model=='k2016' & weight=='base' & timing=='nominal' & definition_day==1 & schedule=='current'", "pct_ge2", 0, "%", "share with two or more current samples on the cliff, 2016"),
                             lloq = lloq))
  cap <- tx("S4.caption", list(pg = dcfg("literature_core_deck.yaml", c("fda_bla761055_clinpharm", "page"), "FDA BLA 761055 clinical pharmacology review page (as given, unverified)", num_fmt(0))))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  fh <- by - 0.08 - y0
  deck_figure(p, "s4_cliff", c(GEO$ML, y0, GEO$CW, fh), src = c(RP, RS))
  dsrc("literature statement Kovalenko 2020 Results", "literature/literature_qualitative.csv", "(table)")
  premise(grepl("steep target-mediated", .read(LC)$fda_bla761055_clinpharm$phrase), "FDA phrase recorded in config/literature_core_deck.yaml")

  deck_notes(tx("S4.notes", c(f, list(pc = pc, id = dderived("representative subject id", RS, "percentile==50 :: id", subj$id, fnum(subj$id, 0)),
    lday = dv(RS, "percentile==50", "study_day_lloq", 1, "", "representative subject: study day at the LLOQ"),
    c16 = dv(CS, "model=='k2016' & weight=='base'", "c_start1_median", 2, "", "median concentration at cliff start, 2016 (mg/L)"),
    c20 = dv(CS, "model=='k2020' & weight=='base'", "c_start1_median", 2, "", "median concentration at cliff start, 2020 (mg/L)"),
    l16 = dv(CS, "model=='k2016' & weight=='base'", "len1_median", 2, "", "median cliff length, 2016 (days)"),
    l20 = dv(CS, "model=='k2020' & weight=='base'", "len1_median", 2, "", "median cliff length, 2020 (days)"),
    p95 = dspan(CS, "model=='k2016' & weight=='base'", "len1_p05", "len1_p95", 2, "", "5th to 95th percentile cliff length, 2016"),
    d16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at the LLOQ, 2016"),
    d20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at the LLOQ, 2020"),
    nsub = dint(CS, "model=='k2016' & weight=='base'", "n", "virtual subjects per model (cliff analysis)"), wt = f_wt_range()))))
  deck_end()
}
