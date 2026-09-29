# S05 전제: 약의 특성, 말단 절벽. 대표 대상자 3명(2016 모델, 60~90 kg, LLOQ 도달일 25·50·75백분위)의 참 농도 곡선(반로그),
# 절벽 구간 음영(start1 ~ t_lloq), 현행 채혈일(B0), 정량한계 선. 절벽 길이·LLOQ 도달일 카드, 원개발사 문헌 서술 표.
# 곡선 자료: results/deck_inputs/rep_profiles.csv, rep_subjects.csv(scripts/62, 커밋된 요약과 대조 확인, checks.csv).
# 문헌 수치는 results/literature/literature_qualitative.csv의 서술에서 정규식으로 읽고 출처를 남긴다(s05_litnum).
s05_litnum <- function(src_rx, k, n_expected, item) {
  LQ <- "literature/literature_qualitative.csv"; w <- sprintf("grepl('%s', source)", src_rx); rx <- "[0-9]+(\\.[0-9]+)?"
  r <- row1(LQ, w); m <- regmatches(r$statement, gregexpr(rx, r$statement))[[1]]
  premise(length(m) == n_expected, sprintf("literature statement [%s] holds %d numbers", w, n_expected))
  dderived(item, LQ, sprintf("%s :: statement, number %d of %d (regex %s)", w, k, n_expected, rx), as.numeric(m[k]), m[k])
}
slide_S05 <- function() {
  CS <- "cliff/cliff_summary.csv"; RP <- "deck_inputs/rep_profiles.csv"; RS <- "deck_inputs/rep_subjects.csv"; CP <- "cliff/cliff_points.csv"
  deck_slide("S05", tag = "litsim")
  d1 <- dcfg("oc_design.yaml", c("cliff", "definition_days"), "cliff start: instantaneous half-life threshold, primary definition (days)", function(x) fnum(x[1], 0))
  len <- drange(CS, "weight=='base'", "len1_median", 2, "", "median cliff length, 1-day definition, two models (60-90 kg)")
  lloq <- f_lloq()
  deck_kicker(tx("S05.kicker")); deck_title(tx("S05.title", list(d1 = d1, len = len)))

  # ---- 그림: 대표 대상자 3명의 참 농도(반로그), 절벽 음영, B0 채혈점, LLOQ 선 ----
  L <- DK$txt$S05$fig
  prof <- copy(rows(RP)); subj <- copy(rows(RS))
  premise(nrow(subj) == 3 && setequal(unique(prof$id), subj$id), "three representative subjects with profiles")
  lloq_v <- .read("config/assay.yaml")$lloq_mg_L$value
  b0 <- unlist(.read("config/trial_design.yaml")$schedules$B0$days)
  subj <- subj[order(percentile)]
  subj[, lab := vapply(seq_len(.N), function(i) fill(L$subj, list(p = fnum(percentile[i], 0), day = fnum(study_day_lloq[i], 1))), "")]
  subj[, lab := factor(lab, levels = lab)]
  prof <- merge(prof, subj[, .(id, lab)], by = "id"); prof[, day := time + 1]
  pts <- prof[round(time, 2) %in% round(b0, 2) & C >= lloq_v]
  premise(nrow(pts) >= 3 * 8, "planned B0 samples above the LLOQ found on the profile grid")
  shade <- subj[, .(lab, x0 = start1 + 1, x1 = t_lloq + 1)]
  COL <- setNames(c("#8fb8ea", PAL$blue, "#15406f"), levels(subj$lab)); SHP <- setNames(c(21, 22, 24), levels(subj$lab)); LT <- setNames(c("solid", "42", "13"), levels(subj$lab))
  sd_ <- b0 + 1; brk <- sd_[abs(sd_ - round(sd_)) < 1e-9 & (round(sd_) - 1) %% 7 == 0]; brk <- c(1, brk)
  ymax <- max(prof$C) * 1.6
  FH <- 3.3
  # 글자 표지는 흰 바탕 상자(격자선이 글자를 지나지 않게). 범례는 곡선이 없는 왼쪽 아래(Day 20 이전, 1 mg/L 아래)에 둔다
  premise(all(prof[day >= 1.5 & day <= 29, C] > 2), "profiles stay above 2 mg/L between study days 1.5 and 29 (legend placed there between the LLOQ line and 2 mg/L)")
  premise(all(prof[day >= max(shade$x1) + 1, C] < 1), "no profile above 1 mg/L after the last cliff (labels placed in the empty right area)")
  lab_ <- function(x, y, s, hjust, col = PAL$ink2, face = "plain", vjust = 0.5)
    annotate("label", x = x, y = y, label = s, hjust = hjust, vjust = vjust, size = 4.2, family = FONT, colour = col, fontface = face, lineheight = 0.95,
             fill = "white", label.size = 0, label.padding = grid::unit(0.12, "lines"), label.r = grid::unit(0, "lines"))
  p <- ggplot() +
    geom_rect(data = shade, aes(xmin = x0, xmax = x1, ymin = 0.01, ymax = ymax), fill = PAL$orange, alpha = 0.28) +
    geom_hline(yintercept = lloq_v, linetype = "22", colour = PAL$ink2, linewidth = 0.6) +
    geom_line(data = prof, aes(day, C, colour = lab, linetype = lab), linewidth = 0.95) +
    geom_point(data = pts, aes(day, C, colour = lab, fill = lab, shape = lab), size = 2.3, stroke = 0.6, colour = "white") +
    lab_(max(sd_) + 2, lloq_v, fill(L$lloq, list(v = lloq)), hjust = 1, vjust = -0.12) +
    lab_(max(shade$x1) + 0.9, 1.6, L$cliff, hjust = 0, col = PAL$orange, face = "bold") +
    lab_(max(sd_) + 2, 12, L$b0, hjust = 1) +
    scale_colour_manual(values = COL, name = L$legend) + scale_fill_manual(values = COL, name = L$legend) + scale_shape_manual(values = SHP, name = L$legend) + scale_linetype_manual(values = LT, name = L$legend) +
    guides(colour = guide_legend(title.position = "top", ncol = 1), fill = guide_legend(title.position = "top", ncol = 1), shape = guide_legend(title.position = "top", ncol = 1, override.aes = list(size = 2.6)),
           linetype = guide_legend(title.position = "top", ncol = 1)) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10), labels = function(x) formatC(x, format = "fg"), limits = c(0.01, ymax), expand = expansion(mult = 0)) +
    scale_x_continuous(breaks = brk, limits = c(0, max(sd_) + 2), expand = expansion(mult = c(0, 0))) +
    labs(x = L$xlab, y = NULL, subtitle = L$ylab) + theme_deck(13) +
    theme(legend.position = c(0.035, 0.285), legend.justification = c(0, 0), legend.key.width = grid::unit(1.8, "lines"), legend.key.height = grid::unit(1.0, "lines"),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, colour = PAL$ink2), legend.title.align = 0,
          legend.background = element_rect(fill = "white", colour = NA), legend.margin = margin(3, 6, 3, 4), legend.spacing.y = grid::unit(1, "pt"),
          plot.title.position = "plot", plot.subtitle = element_text(colour = PAL$ink2, size = 13, margin = margin(0, 0, 4, 0)))
  deck_figure(p, "s05_rep_profiles", c(GEO$ML, GEO$BODY_TOP, 7.25, FH), src = c(RP, RS))
  pc <- dderived("representative subjects: LLOQ-day percentiles", RS, "TRUE :: percentile", subj$percentile, paste(fnum(subj$percentile, 0), collapse = "·"))

  # ---- 카드: 절벽 길이, LLOQ 도달 연구일 ----
  xr <- GEO$ML + 7.25 + 0.3; wr <- GEO$W - GEO$MR - xr
  day16 <- dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, 2016 model")
  day20 <- dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, 2020 model")
  ch <- (FH - 0.12) / 2
  deck_stat(tx("S05.stat1.value", list(len = len)), tx("S05.stat1.label", list(d1 = d1)), c(xr, GEO$BODY_TOP, wr, ch), color = PAL$orange, bg = PAL$tint_orange, value_size = 32)
  deck_stat(tx("S05.stat2.value", list(day = day16)), tx("S05.stat2.label", list(lloq = lloq, day20 = day20)), c(xr, GEO$BODY_TOP + ch + 0.12, wr, ch), value_size = 32)

  # ---- 표: 원개발사 문헌의 서술(프로젝트 기록) ----
  k21 <- "^Kovalenko 2021"
  lit <- list(km = s05_litnum(k21, 1, 5, "Kovalenko 2021: Km (mg/L)"),
              c90 = s05_litnum(k21, 2, 5, "Kovalenko 2021: concentration removing the stated share of circulating target (mg/L)"),
              p90 = s05_litnum(k21, 3, 5, "Kovalenko 2021: share of circulating target removed (%)"),
              hl = s05_litnum(k21, 4, 5, "Kovalenko 2021: instantaneous half-life in the beta phase (days)"),
              hl0 = s05_litnum(k21, 5, 5, "Kovalenko 2021: instantaneous half-life in the target-mediated phase (days)"))
  # 머리글 문구의 전제: 보고서 1.1절에는 Kovalenko 2020, Li 2020, Cohen 2022가 있고 Kovalenko 2021은 없다
  rmd <- readLines(proj_path("regulatory", "src", "MS_report.Rmd"), encoding = "UTF-8"); i11 <- grep("^## 1\\.1 ", rmd); i12 <- grep("^## 1\\.2 ", rmd)
  premise(length(i11) == 1 && length(i12) == 1, "report section 1.1 found"); s11 <- paste(rmd[i11:i12], collapse = "\n")
  premise(all(vapply(c("Kovalenko 2020", "Li 2020", "Cohen 2022"), grepl, TRUE, x = s11, fixed = TRUE)) && !grepl("Kovalenko 2021", s11, fixed = TRUE),
          "report section 1.1 cites Kovalenko 2020, Li 2020 and Cohen 2022 but not Kovalenko 2021 (table header)")
  for (s_ in c("^Kovalenko 2020 Results", "^Kovalenko 2020 Discussion", "^Li 2020", "^Cohen 2022")) dsrc(sprintf("literature statement %s", s_), "literature/literature_qualitative.csv", "(table)")
  T <- DK$txt$S05$table
  df <- data.frame(a = unlist(T$src), b = c(T$k20, fill(T$k21, lit), T$li, T$cohen), stringsAsFactors = FALSE)
  names(df) <- unlist(T$head)
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP + FH + 0.15, GEO$CW, GEO$BODY_BOTTOM - GEO$BODY_TOP - FH - 0.15), widths = c(1.75, 10.48), size = 13, align_num = FALSE)

  # ---- 노트 ----
  CPW <- "model=='k2016' & weight=='base' & timing=='nominal' & definition_day==1 & schedule=='current'"
  cs <- rows(CS, "weight=='base'"); premise(nrow(cs) == 2 && length(unique(fnum(cs$len1_median, 2))) == 1 && length(unique(fnum(cs$len2_median, 2))) == 1, "cliff lengths (1- and 2-day definitions) equal in both models at 2 decimals (notes: the same in both models)")
  premise(row1(CP, CPW)$pct_ge2 == 0, "no subject has 2 or more B0 samples on the cliff (notes: lambda-z estimated before the cliff)")
  deck_notes(tx("S05.notes", list(
    d1 = d1, d2 = dcfg("oc_design.yaml", c("cliff", "definition_days"), "cliff start: sensitivity definition (days)", function(x) fnum(x[2], 0)),
    nsub = dint(CS, "model=='k2016' & weight=='base'", "n", "virtual subjects per model (cliff analysis)"),
    len16 = dv(CS, "model=='k2016' & weight=='base'", "len1_median", 2, "", "median cliff length, 2016 model"),
    p95 = dspan(CS, "model=='k2016' & weight=='base'", "len1_p05", "len1_p95", 2, "", "5th to 95th percentile cliff length, 2016 model"),
    len2 = drange(CS, "weight=='base'", "len2_median", 2, "", "median cliff length, 2-day definition, two models"),
    cst = dv(CS, "model=='k2016' & weight=='base'", "c_start1_median", 2, "", "median concentration at cliff start (mg/L), 2016 model"),
    day16 = day16, day20 = day20,
    r16 = sprintf("%s~%s", dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_p05", 1, "", "5th percentile study day at LLOQ, 2016"), dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_p95", 1, "", "95th percentile study day at LLOQ, 2016")),
    r20 = sprintf("%s~%s", dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_p05", 1, "", "5th percentile study day at LLOQ, 2020"), dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_p95", 1, "", "95th percentile study day at LLOQ, 2020")),
    d58 = dderived("study-day threshold in column name lloq_after_day58_pct", CS, "column name lloq_after_day58_pct (100 x mean(t_lloq + 1 > 58))", 58, "58"),
    a16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_after_day58_pct", 1, "%", "share above LLOQ after study day 58, 2016"),
    a20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_after_day58_pct", 1, "%", "share above LLOQ after study day 58, 2020"),
    int = dint(CP, CPW, "min_interval_day", "minimum sampling interval of the current schedule in the cliff window (days)"),
    ge1 = dv(CP, CPW, "pct_ge1", 1, "%", "share with 1 or more B0 samples on the cliff, 2016, nominal, 1-day definition"),
    ge2 = dv(CP, CPW, "pct_ge2", 1, "%", "share with 2 or more B0 samples on the cliff, 2016, nominal, 1-day definition"),
    pc = pc, lloq = lloq, km = lit$km, wt = f_wt_range(),
    minpts = dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0)))))
  deck_end()
}
