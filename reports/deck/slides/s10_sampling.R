# S10 논거 ①-2 채혈로 해결되지 않음. 시험 모집단(건강인 60~90 kg, B0 기준)만.
#  그림 1: 참 농도가 LLOQ에 닿는 연구일(절벽 위치)의 1일 간격 분포, 두 모델, 현행 B0 채혈일 표시, 절벽 길이 대 채혈 간격 비교 막대.
#          자료 results/deck_inputs/cliff_lloq_day_hist.csv(로컬 대상자 수준 결과에서 만들고 커밋된 요약과 대조: deck_inputs/checks.csv).
#  그림 2: 일정별 절벽 안 채혈점 수(정확히 1점, 2점, 3점 이상), 두 모델, 명목일, 1일 정의. results/cliff/cliff_points.csv의
#          pct_ge1/ge2/ge3는 누적이므로 차이로 나눈다. weight=='base'만(비만 체중 행은 부록 A5 자료).
#  표: B0 대비 AUC0-inf 신뢰 충족률 변화(대응 비교, %p), 세트 (i)·(ii)별, results/reliability/reliability_paired_vs_B0.csv의
#      variant base·struct2020(시험 모집단의 두 모델)만. 요청의 "-3.5~+2.8%p"는 세트 (ii) 최솟값과 세트 (i) 최댓값을 합친 범위라 세트별로 쓴다.
#      보고서의 "세트 (ii)가 떨어지는 칸" 수는 여섯 민감도 변이(체중 50~90 kg 분포 포함)까지 센 값이라 쓰지 않는다.
S10_SCHED <- c(current = "B0", plus_39_46 = "D1", plus_39_46_53 = "D2", plus_32_39_46_53 = "D3", plus_40_47 = "D4", daily_29_57 = "daily")

# config 채혈일(투여 후 일)을 연구일로(연구일 = 투여 후 일 + 1). yaml이 정수·실수 혼합 목록으로 읽으므로 unlist(S08과 같은 사정)
s10_days <- function(schedule) { d <- unlist(.read("config/trial_design.yaml")$schedules[[schedule]]$days); premise(length(d) > 0, paste("schedule", schedule)); d + 1 }
# 후보 일정이 B0에 더하는 연구일(예 "39, 46"). 표 칸 폭보다 긴 한 단어는 deck_lib est_lines가 실수형을 돌려 표 높이 추정이 멈추므로 띄어 쓴다
s10_added <- function(schedule, sep = ", ") {
  a <- sort(setdiff(round(s10_days(schedule), 6), round(s10_days("B0"), 6))); premise(length(a) >= 1 && all(a == round(a)), paste("added whole study days,", schedule))
  dderived(sprintf("schedule %s: study days added to B0", schedule), "config/trial_design.yaml", sprintf("setdiff(schedules.%s.days, schedules.B0.days) + 1", schedule), a, paste(fnum(a, 0), collapse = sep))
}
# 부호 표시(+). 인쇄값(추적 행)은 d* 함수가 남긴 그대로이고 "+"만 붙인다
s10_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?$", p)) p else paste0("+", p)
s10_signed_rng <- function(p) { p <- s10_signed(p); sub("~([0-9])", "~+\\1", p) }
# 매일 채혈(참고) 일정의 연구일 범위: 일정 코드 daily_<첫날>_<끝날>에서 읽고, cliff_schedules.csv의 추가일이 그 사이에 있는지 확인
s10_daily_range <- function(CP) {
  code <- unique(rows(CP, "schedule=='daily_29_57'")$schedule); m <- as.numeric(regmatches(code, gregexpr("[0-9]+", code))[[1]])
  add <- as.numeric(strsplit(rows("cliff/cliff_schedules.csv", "schedule=='daily_29_57'")$added_study_days, " ")[[1]])
  premise(length(m) == 2 && min(add) == m[1] + 1 && max(add) == m[2] - 1, "daily schedule code gives its first and last study day (added days lie between)")
  dderived("daily sampling schedule: first and last study day (from schedule code)", "cliff/cliff_schedules.csv", "schedule=='daily_29_57' :: digits of the schedule code; added_study_days lie between", m, rng_fmt(m[1], m[2], 0))
}

slide_S10 <- function() {
  CH <- "deck_inputs/cliff_lloq_day_hist.csv"; CS <- "cliff/cliff_summary.csv"; CP <- "cliff/cliff_points.csv"; CSC <- "cliff/cliff_schedules.csv"
  RP <- "reliability/reliability_paired_vs_B0.csv"; LW <- "reliability/reliability_lz_window_by_schedule.csv"
  NOM <- "weight=='base' & timing=='nominal' & definition_day==1"; WIN <- "weight=='base' & timing=='windowed' & definition_day==1"
  FIX <- "schedule!='daily_29_57'"; DS <- c("D1", "D2", "D3", "D4"); RW <- "variant %in% c('base','struct2020') & schedule %in% c('D1','D2','D3','D4')"
  VAR <- c(k2016 = "base", k2020 = "struct2020")
  deck_slide("S10", tag = "sim")
  L <- DK$txt$S10

  # ---- 전제: 문장이 기대는 사실 ----
  sc <- rows(CSC)
  for (k in names(S10_SCHED)[2:5]) {                      # cliff 파일의 일정 코드와 config D1~D4의 추가일이 같다
    a <- sort(round(setdiff(s10_days(S10_SCHED[[k]]), s10_days("B0")), 6)); b <- sort(as.numeric(strsplit(sc[schedule == k, added_study_days], " ")[[1]]))
    premise(isTRUE(all.equal(a, b)), sprintf("cliff schedule %s adds the study days of config schedule %s", k, S10_SCHED[[k]]))
  }
  cp <- copy(rows(CP, NOM)); premise(nrow(cp) == 12 && setequal(cp$schedule, names(S10_SCHED)), "points on the cliff: six schedules x two models (base, nominal, 1-day definition)")
  premise(all(cp[schedule != "daily_29_57", pct_ge2] == 0), "no fixed schedule puts two or more nominal samples on the cliff (both models)")
  premise(all(rows(CP, paste(WIN, "&", FIX))$pct_ge2 < 0.1), "windowed fixed schedules: two or more samples on the cliff below 0.1% (title bound)")
  premise(all(cp$pct_ge3 == 0), "no schedule, even daily sampling, puts three nominal samples on the cliff")
  rp <- rows(RP, RW); premise(nrow(rp) == 8, "reliability change: two models x D1 to D4")
  premise(!any(rp$crit_c_i %in% TRUE) && !any(rp$crit_c_ii %in% TRUE), "criterion (c) met in no cell under either criteria set")
  premise(max(rp$gain_i_pp) > max(rp$gain_ii_pp), "the largest gain is under set (i) (title: largest gain)")
  premise(rows(CS, "weight=='base'")[, all(lloq_never_pct == 0)], "every subject reaches the LLOQ (histogram sums to 100%)")

  # ---- 제목 ----
  f <- list(ge2 = dext(CP, paste(WIN, "&", FIX), "pct_ge2", max, 2, "%", "two or more samples on the cliff, fixed schedules, windowed days, max over schedules and models"),
            gmax = dext(RP, RW, "gain_i_pp", max, 2, "", "largest reliability gain versus B0 (set (i); set (ii) lower), trial population, D1 to D4"))
  deck_kicker(tx("S10.kicker")); deck_title(tx("S10.title", f))

  # ---- 그림 1: 절벽 위치 분포(두 모델), B0 채혈일 ----
  XL <- GEO$ML; WL <- 6.75; FH1 <- 2.62
  h <- copy(rows(CH)); premise(all(abs(h[, sum(pct), by = pk_model]$V1 - 100) < 1e-6), "LLOQ-day histogram sums to 100% per model")
  h[, model := factor(model_lab()[pk_model], levels = model_lab())]; h[, x := day_bin + 0.5]
  b0 <- s10_days("B0"); XMIN <- 14; XMAX <- 72; vis <- b0[b0 >= XMIN & abs(b0 - round(b0)) < 1e-9]
  len <- drange(CS, "weight=='base'", "len1_median", 2, "", "median cliff length, 1-day definition, two models (60-90 kg)")
  int0 <- row1(CP, paste(NOM, "& model=='k2016' & schedule=='current'"))$min_interval_day
  lenv <- rows(CS, "model=='k2016' & weight=='base'")$len1_median
  ymax <- max(h$pct) * 1.28; x0 <- 58.5
  cmp <- data.table(y = ymax * c(0.9, 0.74), x1 = x0 + c(int0, lenv), col = c("b0", "cliff"),
                    lab = c(fill(L$fig$cmp_b0, list(v = fnum(int0, 0))), fill(L$fig$cmp_cliff, list(v = fnum(lenv, 2)))))
  p1 <- ggplot(h, aes(x = x, y = pct)) +
    geom_vline(xintercept = vis, colour = PAL$muted, linewidth = 0.45, linetype = "22") +
    geom_line(aes(colour = model, linetype = model), linewidth = 0.95) +
    geom_segment(data = cmp, aes(x = x0, xend = x1, y = y, yend = y, colour = NULL, linetype = NULL), inherit.aes = FALSE,
                 colour = c(PAL$ink2, PAL$orange), linewidth = c(1.2, 3.2), lineend = "butt") +
    geom_text(data = cmp, aes(x = x0, y = y, label = lab), inherit.aes = FALSE, vjust = -0.7, hjust = 0, size = 3.9, family = FONT, colour = c(PAL$ink2, PAL$orange)) +
    annotate("text", x = min(vis) + 0.4, y = ymax * 0.99, label = L$fig$b0, hjust = 0, vjust = 1, size = 3.9, family = FONT, colour = PAL$ink2) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_linetype_manual(values = unname(MODEL_LT)) +
    scale_x_continuous(breaks = vis, expand = expansion(add = 0)) + coord_cartesian(xlim = c(XMIN, XMAX), ylim = c(0, ymax), expand = FALSE) +
    labs(x = L$fig$xlab1, y = NULL, subtitle = L$fig$ylab1) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.major.x = element_blank(), plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))
  deck_figure(p1, "s10_lloq_day_hist", c(XL, GEO$BODY_TOP, WL, FH1), src = c(CH, CS, CP))

  # ---- 그림 2: 일정별 절벽 안 채혈점 수(정확히 k점), 두 모델 ----
  kk <- as.integer(sub("pct_ge", "", c("pct_ge1", "pct_ge2", "pct_ge3")))           # 열 이름의 점 수
  premise(kk[3] == .read("config/nca_rules.yaml")$standard$lambda_z$min_points, "the top category (three or more points) is the lambda-z minimum (legend)")
  CAT <- c(fill(L$fig$k1, list(k = kk[1])), fill(L$fig$k2, list(k = kk[2])), fill(L$fig$k3, list(k = kk[3])))
  cp[, sch := factor(ifelse(schedule == "daily_29_57", L$fig$daily, S10_SCHED[schedule]), levels = c(unname(S10_SCHED[1:5]), L$fig$daily))]
  cp[, model := factor(model_lab()[model], levels = model_lab())]
  lg <- rbind(cp[, .(model, sch, cat = CAT[1], v = pct_ge1 - pct_ge2)], cp[, .(model, sch, cat = CAT[2], v = pct_ge2 - pct_ge3)], cp[, .(model, sch, cat = CAT[3], v = pct_ge3)])
  lg[, cat := factor(cat, levels = rev(CAT))]
  p2 <- ggplot(lg, aes(x = sch, y = v, fill = cat)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    geom_text(data = cp, aes(x = sch, y = pct_ge1, label = fnum(pct_ge1, 1)), inherit.aes = FALSE, vjust = -0.4, size = 3.7, family = FONT, colour = PAL$ink) +
    facet_wrap(~model, nrow = 1) +
    scale_fill_manual(values = setNames(c(PAL$orange, PAL$blue, "#9cc3ef"), rev(CAT)), breaks = CAT) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = NULL, subtitle = L$fig$ylab2) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          panel.grid.major.x = element_blank(), panel.spacing = grid::unit(16, "pt"), plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))
  FY2 <- GEO$BODY_TOP + FH1 + 0.08
  deck_figure(p2, "s10_points_in_cliff", c(XL, FY2, WL, GEO$BODY_BOTTOM - FY2), src = CP)

  # ---- 오른쪽 위: 요점 ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  minpts <- dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0))
  b <- list(p16 = dspan(CS, "model=='k2016' & weight=='base'", "lloq_studyday_p05", "lloq_studyday_p95", 1, "", "study day at LLOQ, 5th to 95th percentile, 2016 model"),
            len = len, int = dint(CP, paste(NOM, "& model=='k2016' & schedule=='current'"), "min_interval_day", "minimum sampling interval of B0 in the cliff window (days)"),
            mi = dint(CSC, "schedule=='plus_39_46'", "min_interval_day_from_day22", "minimum interval after Day 22 with added samples (days)"),
            b0 = drange(CP, paste(NOM, "& schedule=='current'"), "pct_ge1", 1, "%", "one or more samples on the cliff, B0, two models"),
            d3 = drange(CP, paste(NOM, "& schedule=='plus_32_39_46_53'"), "pct_ge1", 1, "%", "one or more samples on the cliff, D3, two models"),
            minpts = minpts, dly = drange(CP, paste(NOM, "& schedule=='daily_29_57'"), "pct_ge3", 1, "%", "three or more samples on the cliff, daily sampling, two models"),
            ri = s10_signed_rng(drange(RP, RW, "gain_i_pp", 2, "", "reliability change versus B0, set (i), two models, D1 to D4")),
            rii = s10_signed_rng(drange(RP, RW, "gain_ii_pp", 2, "", "reliability change versus B0, set (ii), two models, D1 to D4")),
            ndown = dcount(RP, paste(RW, "& gain_ii_pp < 0"), "cells where reliability under set (ii) falls, trial population"),
            nall = dcount(RP, RW, "added-sampling cells, trial population (two models x D1 to D4)"),
            thr = dcfg("trial_design.yaml", c("decision_rule", "reliability_gain_pp_min"), "criterion (c) threshold (points)", num_fmt(0)))
  BH <- 3.02
  deck_bullets(tx("S10.bullets", b), box = c(XR, GEO$BODY_TOP, WR, BH), size = 16, gap_pt = 7)

  # ---- 오른쪽 아래: 신뢰 충족률 변화 표(세트별, 모델별, D1~D4) ----
  T <- L$table
  cell <- function(m, s_, set) s10_signed(dv(RP, sprintf("variant=='%s' & schedule=='%s'", VAR[[m]], s_), sprintf("gain_%s_pp", set), 2, "",
                                              sprintf("reliability change versus B0, set (%s), %s, %s", set, m, s_)))
  rowv <- function(m, set) vapply(DS, function(s_) cell(m, s_, set), "")
  df <- as.data.frame(rbind(c(T$days, vapply(DS, s10_added, "")),
                            c(T$i16, rowv("k2016", "i")), c(T$i20, rowv("k2020", "i")),
                            c(T$ii16, rowv("k2016", "ii")), c(T$ii20, rowv("k2020", "ii"))), stringsAsFactors = FALSE)
  names(df) <- c(T$head, DS)
  TY <- GEO$BODY_TOP + BH + 0.1
  deck_table(df, box = c(XR, TY, WR, GEO$BODY_BOTTOM - TY), widths = c(1.72, 0.8, 0.8, 0.92, 0.8), size = 12, highlight = 1, highlight_fill = PAL$tint_grey)

  # ---- 노트 ----
  d3ci <- function(m, set) dci(RP, sprintf("variant=='%s' & schedule=='D3'", VAR[[m]]), sprintf("gain_%s_pp", set), sprintf("gain_%s_lo", set), sprintf("gain_%s_hi", set), 2, "%p",
                               sprintf("reliability change versus B0, D3, set (%s), %s, with paired 95%% CI", set, m))
  mn <- rp[which.min(gain_ii_pp)]; premise(mn$variant == "struct2020" && mn$schedule == "D4", "lowest set (ii) change is the 2020 model D4 (notes)")
  lw <- rows(LW, "variant %in% c('base','struct2020') & schedule %in% c('B0','D1','D2','D3','D4')")
  premise(all(lw[schedule == "B0", min_3pt_nominal_window_days] > lw[schedule != "B0", max(min_3pt_nominal_window_days)]), "added samples shorten the shortest three-point window (notes)")
  premise(all(lw[schedule == "B0", window_median] > lw[schedule == "D3", window_median]), "D3 gives a shorter median lambda-z window than B0 in both models (notes)")
  deck_notes(tx("S10.notes", list(
    wt = f_wt_range(), n = dint(CS, "model=='k2016' & weight=='base'", "n", "virtual subjects per model (cliff analysis)"), lloq = f_lloq(),
    d1def = dcfg("oc_design.yaml", c("cliff", "definition_days"), "cliff start: instantaneous half-life threshold, primary definition (days)", function(x) fnum(x[1], 0)),
    d2def = dcfg("oc_design.yaml", c("cliff", "definition_days"), "cliff start: sensitivity definition (days)", function(x) fnum(x[2], 0)),
    dyr = s10_daily_range(CP),
    m16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, 2016 model"),
    m20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, 2020 model"),
    r16 = b$p16, r20 = dspan(CS, "model=='k2020' & weight=='base'", "lloq_studyday_p05", "lloq_studyday_p95", 1, "", "study day at LLOQ, 5th to 95th percentile, 2020 model"),
    len = len, int = b$int, mi = b$mi, d22 = dderived("study day in column name min_interval_day_from_day22", CSC, "column name min_interval_day_from_day22", 22, "22"),
    dmi = dint(CSC, "schedule=='daily_29_57'", "min_interval_day_from_day22", "minimum interval with daily sampling (days)"),
    b0 = b$b0, d1 = drange(CP, paste(NOM, "& schedule=='plus_39_46'"), "pct_ge1", 1, "%", "one or more samples on the cliff, D1, two models"),
    d2 = drange(CP, paste(NOM, "& schedule=='plus_39_46_53'"), "pct_ge1", 1, "%", "one or more samples on the cliff, D2, two models"),
    d3 = b$d3, d4 = drange(CP, paste(NOM, "& schedule=='plus_40_47'"), "pct_ge1", 1, "%", "one or more samples on the cliff, D4, two models"),
    g2n = dext(CP, paste(NOM, "&", FIX), "pct_ge2", max, 1, "%", "two or more samples on the cliff, fixed schedules, nominal days, max"),
    g2w = f$ge2,
    dy2 = drange(CP, paste(NOM, "& schedule=='daily_29_57'"), "pct_ge2", 1, "%", "two or more samples on the cliff, daily sampling, nominal, two models"),
    dy3 = b$dly, dy3w = drange(CP, paste(WIN, "& schedule=='daily_29_57'"), "pct_ge3", 2, "%", "three or more samples on the cliff, daily sampling, windowed, two models"),
    g2d2 = dext(CP, "weight=='base' & timing=='nominal' & definition_day==2 & schedule!='daily_29_57'", "pct_ge2", max, 1, "%", "two or more samples on the cliff, fixed schedules, nominal, 2-day definition, max"),
    g2d2w = dext(CP, "weight=='base' & timing=='windowed' & definition_day==2 & schedule!='daily_29_57'", "pct_ge2", max, 2, "%", "two or more samples on the cliff, fixed schedules, windowed, 2-day definition, max"),
    minpts = minpts,
    nrel = dint(RP, "variant=='base' & schedule=='D3'", "n", "subjects per model, paired reliability comparison"),
    ci16 = d3ci("k2016", "i"), ci20 = d3ci("k2020", "i"), cii16 = d3ci("k2016", "ii"), cii20 = d3ci("k2020", "ii"),
    mn = dci(RP, "variant=='struct2020' & schedule=='D4'", "gain_ii_pp", "gain_ii_lo", "gain_ii_hi", 2, "%p", "lowest set (ii) change, 2020 model D4, with paired 95% CI"),
    thr = b$thr, ndown = b$ndown, nall = b$nall, sp2 = f_set("ii", "span"),
    w0 = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "min_3pt_nominal_window_days", 0, "", "shortest three-point nominal window after Day 22, B0 (days)"),
    w1 = drange(LW, "variant %in% c('base','struct2020') & schedule %in% c('D1','D2','D3','D4')", "min_3pt_nominal_window_days", 0, "", "shortest three-point nominal window, D1 to D4 (days)"),
    wm0 = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "window_median", 1, "", "median lambda-z window, B0, two models (days)"),
    wm3 = drange(LW, "variant %in% c('base','struct2020') & schedule=='D3'", "window_median", 1, "", "median lambda-z window, D3, two models (days)"))))
  deck_end()
}
