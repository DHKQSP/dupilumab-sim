# S16 논거 ③-1 개인 수준: AUC0-last는 총노출을 담는다. 창 포착률(채혈 구간 커버리지) 분포와 요약(중앙값, 5번째 백분위수, 최소, 80% 미만),
# 관측 대 참 비(AUC0-last 대 NCA AUC0-inf 규칙 A 세트 (i)·규칙 B)를 같은 분석군끼리 나란히. 시험 모집단(건강인 60~90 kg, B0), 모델당 20,000명.
# 수치: results/trialpop/tp_coverage_individual.csv, results/rationale/pillar1_coverage_B0.csv(group 'all'; reg_helpers가 한글 라벨을 영문 코드로 바꿈).
# 분포 그림: results/deck_inputs/coverage_hist.csv(로컬 대상자 수준 결과에서 만든 0.02 간격 분포, 커밋된 요약과 대조: deck_inputs/checks.csv).
# 주의: coverage의 p95 열은 위쪽 꼬리(~1.000)라 정보가 적다(요청 항목이라 표에 함께 싣는다). 의미 있는 값은 아래쪽 5번째 백분위수(p05)다.
# 색: 오른쪽 그림의 파랑·주황은 평가변수(AUC0-last, NCA AUC0-inf)라서 왼쪽 분포 그림의 모델은 회색 두 톤으로 구분한다.
#       'AUC0-inf 개인 최대 5배 이상'은 2016 모델 규칙 B(λz 산출 가능군)에만 해당한다. 2020 모델 규칙 B는 1.87, 규칙 A 세트 (i)는 1.40/1.29.
# 창 포착률 80% 미만 대상자: 수(최소가 80% 이상이면 0명, 개인 수준 요약의 최솟값으로 확인)와 Wilson 95% CI(보고서 열 'est [lo, hi]')
s16_lt80 <- function(m) {
  PC <- "rationale/pillar1_coverage_B0.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; w <- sprintf("model=='%s' & group=='all'", m)
  x <- row1(PC, w)$coverage_lt80_pct_ci; mm <- regmatches(x, regexec("^([0-9.]+) \\[([0-9.]+), ([0-9.]+)\\]$", x))[[1]]
  premise(length(mm) == 4, "coverage_lt80_pct_ci format 'est [lo, hi]'")
  mn <- row1(TCV, sprintf("metric=='window coverage (true AUC0-tlast / true AUC0-inf)' & pk_model=='%s'", m))$min
  premise(mn >= 0.80 && as.numeric(mm[2]) == 0, sprintf("no subject below 80%% window coverage (%s): smallest subject %.4f, estimate %s", m, mn, mm[2]))
  n0 <- dderived(sprintf("subjects with window coverage below 80 percent, %s", m), TCV, sprintf("metric=='window coverage ...' & pk_model=='%s' :: min = %.4f >= 0.80, so the count below 0.80 is 0", m, mn), 0, "0")
  ci <- dderived(sprintf("share of subjects with window coverage below 80 percent, Wilson 95%% CI, %s", m), PC, sprintf("%s :: coverage_lt80_pct_ci (interval)", w), x, sprintf("%s~%s%%", mm[3], mm[4]))
  fill(DK$txt$S16$cov_table$lt, list(n = n0, ci = ci))
}
slide_S16 <- function() {
  TCV <- "trialpop/tp_coverage_individual.csv"; PC <- "rationale/pillar1_coverage_B0.csv"; CH <- "deck_inputs/coverage_hist.csv"
  MW <- "window coverage (true AUC0-tlast / true AUC0-inf)"
  MET <- c(lz_last = "observed-to-true, AUClast (lambda-z estimable)", lz_B = "observed-to-true, AUCinf rule B (lambda-z estimable)",
           i_last = "observed-to-true, AUClast (meeting set (i))", i_A = "observed-to-true, AUCinf rule A (meeting set (i))",
           all_last = "observed-to-true, AUClast (all subjects)")
  M <- c("k2016", "k2020")
  premise(all(unlist(MET) %in% rows(TCV)$metric) && MW %in% rows(TCV)$metric, "observed-to-true metric labels present in tp_coverage_individual.csv")
  wm <- function(met, m) sprintf("metric=='%s' & pk_model=='%s'", met, m)
  deck_slide("S16", tag = "sim")
  L <- DK$txt$S16

  # ---- 제목: 창 포착률 중앙값(두 모델), 모든 대상자의 최소 ----
  f <- list(med = drange(TCV, sprintf("metric=='%s'", MW), "median", 1, "%", "window coverage, median, two models", scale = 100),
            min = dext(TCV, sprintf("metric=='%s'", MW), "min", min, 1, "%", "window coverage, smallest subject over both models", scale = 100))
  deck_kicker(tx("S16.kicker")); deck_title(tx("S16.title", f))

  # ---- 왼쪽 위: 정의 ----
  XL <- GEO$ML; WL <- 5.6
  DH <- 2.14
  deck_text(tx("S16.defs"), c(XL, GEO$BODY_TOP, WL, DH), size = 16, label = "text_defs", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 6)

  # ---- 왼쪽 가운데: 창 포착률 분포(두 모델, 0.02 간격 막대, 로그 세로축) ----
  # 색: 이 슬라이드에서 파랑·주황은 평가변수(오른쪽 그림: AUC0-last, NCA AUC0-inf)이므로 모델은 회색 두 톤으로 구분한다.
  thr <- dderived("window coverage threshold in column name coverage_lt80_pct_ci (percent)", PC, "column name coverage_lt80_pct_ci", 80, "80")
  hw <- copy(rows(CH, "metric=='window'"))
  premise(nrow(hw) > 0 && all(abs(hw[, sum(pct), by = pk_model]$V1 - 100) < 1e-6), "window coverage histogram sums to 100% per model")
  bw <- unique(round(diff(sort(unique(hw$bin_lo))), 6)); premise(length(bw) == 1, "window coverage bins have one width")
  premise(all(hw$pct > 0.005), "every window coverage bin is above the log-axis floor")
  GREY <- setNames(c(PAL$ink2, "#b9b8b2"), model_lab())
  hw[, model := factor(model_lab()[pk_model], levels = model_lab())]
  hw[, `:=`(x0 = 100 * bin_lo + ifelse(pk_model == "k2016", 0.12, 1.0), x1 = 100 * bin_lo + ifelse(pk_model == "k2016", 0.98, 1.86))]
  YF <- 0.004
  lab2 <- dcast(hw, bin_lo ~ pk_model, value.var = "pct")[order(-bin_lo)][1:2]      # 위 두 구간만 값 표시(2016 / 2020 모델, %)
  lab2[, lab := vapply(seq_len(.N), function(i) fill(L$fig$cov_lab, list(a = fnum(k2016[i], 1), b = fnum(k2020[i], 1))), "")]
  # 가장 높은 구간 값은 막대 바로 위(패널 오른쪽 끝에 맞춤). 그다음 구간 값은 막대 왼쪽 위 빈 곳에 두고 짧은 지시선으로 그 막대를 가리킨다
  # (두 값이 가로로 겹치지 않게). 흰 바탕 글상자로 눈금선이 글자를 지나가지 않게 한다
  lab2[, `:=`(x = 100 * bin_lo + c(2.25, -0.75), y = pmax(k2016, k2020) * c(1.8, 3), vj = c(0, 0.5))]
  lead2 <- lab2[2][, .(x = x + 0.1, y = y * 0.85, xend = 100 * bin_lo + 0.35, yend = pmax(k2016, k2020) * 1.25)]
  mn <- vapply(M, function(m) row1(TCV, wm(MW, m))$min, 1)
  premise(all(mn >= min(hw$bin_lo)) && all(mn < min(hw$bin_lo) + bw), "the smallest subject of each model lies in the lowest histogram bin")
  low <- dcast(hw[bin_lo == min(bin_lo)], bin_lo ~ pk_model, value.var = "n")                     # 가장 낮은 구간의 대상자 수(아래 꼬리)
  mlab <- data.table(x = 100 * low$bin_lo + 0.1, y = max(hw[bin_lo == min(bin_lo), pct]) * 2.2, lab = fill(L$fig$cov_min, list(a = fint(low$k2016), b = fint(low$k2020))))
  pw <- ggplot(hw) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = YF, ymax = pct, fill = model), colour = NA) +
    geom_segment(data = lead2, aes(x = x, y = y, xend = xend, yend = yend), colour = PAL$ink2, linewidth = 0.35) +
    geom_label(data = lab2, aes(x = x, y = y, label = lab, vjust = vj), hjust = 1, size = 3.8, family = FONT, colour = PAL$ink, fill = "white",
               label.size = 0, label.padding = grid::unit(0.06, "lines"), label.r = grid::unit(0, "lines")) +
    geom_text(data = mlab, aes(x = x, y = y, label = lab), vjust = 0, hjust = 0, size = 3.8, family = FONT, colour = PAL$ink) +
    scale_fill_manual(values = GREY, name = L$fig$cov_leg) +
    scale_x_continuous(breaks = seq(84, 100, by = 4), labels = function(x) paste0(fnum(x, 0), "%"), limits = c(100 * min(hw$bin_lo) - 0.3, 100.3), expand = expansion(0)) +
    scale_y_log10(limits = c(YF, 1000), breaks = c(0.01, 0.1, 1, 10, 100), labels = function(x) ifelse(x %in% c(0.01, 1, 100), formatC(x, format = "fg"), ""), expand = expansion(0)) +
    labs(x = L$fig$cov_xlab, y = NULL) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          legend.title = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 10, 0, 0)), panel.grid.major.x = element_blank())
  FY <- GEO$BODY_TOP + DH + 0.08; TH <- 1.18; FH <- GEO$BODY_BOTTOM - FY - TH - 0.08      # TH = 머리글 두 줄 + 두 행(행 높이가 같게)
  deck_figure(pw, "s16_window_coverage_hist", c(XL, FY, WL, FH), src = CH)

  # ---- 왼쪽 아래: 창 포착률 요약 표(두 모델) ----
  cv <- function(m, col) dv(TCV, wm(MW, m), col, 1, "%", sprintf("window coverage, %s, %s", col, m), scale = 100)
  dfc <- data.frame(a = unlist(L$cov_table$rows[M]), b = vapply(M, cv, "", col = "median"), c = vapply(M, cv, "", col = "p05"), c2 = vapply(M, cv, "", col = "p95"),
                    d = vapply(M, cv, "", col = "min"), e = vapply(M, s16_lt80, ""), stringsAsFactors = FALSE, check.names = FALSE)
  names(dfc) <- tx("S16.cov_table.head", list(thr = thr))
  TY <- GEO$BODY_BOTTOM - TH
  deck_table(dfc, box = c(XL, TY, WL, TH), widths = c(0.96, 0.73, 0.76, 0.8, 0.73, 1.62), size = 12, label = "table_cov")

  # ---- 오른쪽 위: 관측 대 참 비 분포(AUC0-last 전체, NCA AUC0-inf 규칙 A 세트 (i), 규칙 B), 모델별 ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  ho <- copy(rows(CH, "metric %in% c('auclast','aucinf_A','aucinf_B')"))
  premise(all(abs(ho[, sum(pct), by = .(pk_model, metric)]$V1 - 100) < 1e-6), "observed-to-true histograms sum to 100% per model and metric")
  premise(all(ho[metric == "auclast", sum(n), by = pk_model]$V1 == vapply(M, function(m) row1(TCV, wm(MET[["all_last"]], m))$n, 1)), "AUClast histogram covers all subjects (legend)")
  top <- max(ho$bin_lo)                                       # 마지막 구간은 top 이상 전부(scripts/62의 구간: 0.5~1.5, 1.5 이상은 한 구간)
  premise(abs(top - 1.5) < 1e-9, "last histogram bin is the open bin from 1.5 (scripts/62 breaks)")
  SER <- c(auclast = L$fig$s_last, aucinf_A = L$fig$s_A, aucinf_B = L$fig$s_B)
  ho[, ser := factor(SER[metric], levels = SER)]; ho[, model := factor(model_lab()[pk_model], levels = model_lab())]
  body_ <- ho[bin_lo < top]; body_[, x := bin_lo + bw / 2]
  ovf <- ho[bin_lo >= top]; ovf[, x := top + bw * 2]
  premise(nrow(ovf) > 0 && all(ovf$metric == "aucinf_B"), "only rule B has subjects in the open bin (label marker is the rule B square)")
  ovf[, lab := vapply(seq_len(.N), function(i) fill(L$fig$ovf, list(v = fnum(top, 1), n = fint(n[i]), p = fnum(pct[i], if (pct[i] < 0.01) 3 else 2))), "")]
  LY <- 30                                                     # 열린 구간 값: 패널 위쪽 오른쪽(점 바로 위에서 끝나게), 점까지 세로 지시선
  COLS <- setNames(c(PAL$blue, "#a8441c", PAL$orange), SER); LTS <- setNames(c("solid", "22", "solid"), SER); SHP <- setNames(c(16, 17, 15), SER)
  po <- ggplot(body_, aes(x = x, y = pct, colour = ser, linetype = ser)) +
    geom_vline(xintercept = 1, colour = PAL$ink2, linewidth = 0.5) +
    geom_line(linewidth = 0.85) +
    geom_segment(data = ovf, aes(x = x, xend = x, y = LY * 0.62, yend = pct * 1.7), linetype = "solid", linewidth = 0.35, show.legend = FALSE) +
    geom_point(data = ovf, aes(shape = ser), size = 2.8, show.legend = FALSE) +
    geom_text(data = ovf, aes(x = x + bw * 0.6, y = LY, label = lab), hjust = 1, size = 3.7, family = FONT, show.legend = FALSE) +
    facet_wrap(~model, nrow = 1) +
    scale_colour_manual(values = COLS) + scale_linetype_manual(values = LTS) + scale_shape_manual(values = SHP) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10), labels = function(x) formatC(x, format = "fg"), limits = c(0.004, 100)) +
    scale_x_continuous(breaks = c(0.6, 0.8, 1, 1.2, top + bw * 2), labels = c(fnum(c(0.6, 0.8, 1, 1.2), 1), fill(L$fig$ovf_tick, list(v = fnum(top, 1)))), limits = c(0.55, top + bw * 3.4)) +
    labs(x = L$fig$ot_xlab, y = NULL, subtitle = L$fig$ot_ylab) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.0, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.spacing = grid::unit(18, "pt"), plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))

  # ---- 오른쪽 아래: 같은 분석군끼리 AUC0-last 대 NCA AUC0-inf (중앙값, 최대; 2016 / 2020 모델) ----
  q <- function(k, m, col, d) dv(TCV, wm(MET[[k]], m), col, d, "", sprintf("%s, %s, %s", MET[[k]], col, m))
  two <- function(k, col, d) sprintf("%s / %s", q(k, "k2016", col, d), q(k, "k2020", col, d))
  T <- L$ot_table; KS <- c("lz_last", "lz_B", "i_last", "i_A")
  dfo <- data.frame(a = unlist(T[KS]), b = vapply(KS, two, "", col = "median", d = 3), c = vapply(KS, two, "", col = "max", d = 2), stringsAsFactors = FALSE, check.names = FALSE)
  names(dfo) <- unlist(T$head)
  # 전제: 개인 최대가 5배를 넘는 것은 2016 모델 규칙 B뿐(나머지 NCA AUC0-inf 최대는 2 미만)
  r_ <- rows(TCV, "grepl('AUCinf', metric) & grepl('observed-to-true', metric)")
  premise(r_[pk_model == "k2016" & metric == MET[["lz_B"]], max] > 5 && all(r_[!(pk_model == "k2016" & metric == MET[["lz_B"]]), max] < 2),
          "only the 2016 model rule B maximum exceeds 5; the other NCA AUC0-inf maxima are below 2")
  KH <- 0.98; OTH <- 1.78
  OTY <- GEO$BODY_BOTTOM - KH - 0.06 - OTH; OH <- OTY - 0.08 - GEO$BODY_TOP
  deck_figure(po, "s16_observed_to_true_hist", c(XR, GEO$BODY_TOP, WR, OH), src = CH)
  deck_table(dfo, box = c(XR, OTY, WR, OTH), widths = c(2.25, 2.04, 2.04), size = 12, highlight = 2, label = "table_ot")
  deck_text(tx("S16.takeaway", list(b16 = q("lz_B", "k2016", "max", 2), b20 = q("lz_B", "k2020", "max", 2), a16 = q("i_A", "k2016", "max", 2), a20 = q("i_A", "k2020", "max", 2))),
            c(XR, GEO$BODY_BOTTOM - KH, WR, KH), size = 16, label = "text_takeaway")

  # ---- 노트 ----
  n16 <- dint(TCV, wm(MW, "k2016"), "n", "subjects per model, trial population B0")
  deck_notes(tx("S16.notes", list(
    n = n16, wt = f_wt_range(), dose = f_dose(),
    p95a = dv(TCV, wm(MW, "k2016"), "p95", 3, "", "window coverage, p95 (upper tail), k2016"),
    p95b = dv(TCV, wm(MW, "k2020"), "p95", 3, "", "window coverage, p95 (upper tail), k2020"),
    ex16 = dv(PC, "model=='k2016' & group=='all'", "extrap_true_p95", 2, "%", "true extrapolation, 95th percentile, k2016"),
    ex20 = dv(PC, "model=='k2020' & group=='all'", "extrap_true_p95", 2, "%", "true extrapolation, 95th percentile, k2020"),
    exm16 = dv(PC, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "true extrapolation, median, k2016"),
    exm20 = dv(PC, "model=='k2020' & group=='all'", "extrap_true_median", 2, "%", "true extrapolation, median, k2020"),
    thr = thr,
    nB16 = dint(TCV, wm(MET[["lz_B"]], "k2016"), "n", "lambda-z estimable subjects, k2016"),
    nB20 = dint(TCV, wm(MET[["lz_B"]], "k2020"), "n", "lambda-z estimable subjects, k2020"),
    nA16 = dint(TCV, wm(MET[["i_A"]], "k2016"), "n", "subjects meeting set (i), k2016"),
    nA20 = dint(TCV, wm(MET[["i_A"]], "k2020"), "n", "subjects meeting set (i), k2020"),
    sdl16 = q("lz_last", "k2016", "sd_log", 3), sdB16 = q("lz_B", "k2016", "sd_log", 3),
    sdl20 = q("lz_last", "k2020", "sd_log", 3), sdB20 = q("lz_B", "k2020", "sd_log", 3),
    sdi16 = q("i_last", "k2016", "sd_log", 3), sdA16 = q("i_A", "k2016", "sd_log", 3),
    lmin16 = q("all_last", "k2016", "min", 3), lmin20 = q("all_last", "k2020", "min", 3),
    fx16 = dv(TCV, wm("extrapolation factor AUCinf / AUClast (lambda-z estimable)", "k2016"), "max", 2, "", "extrapolation factor AUCinf/AUClast, max, k2016"),
    fx20 = dv(TCV, wm("extrapolation factor AUCinf / AUClast (lambda-z estimable)", "k2020"), "max", 2, "", "extrapolation factor AUCinf/AUClast, max, k2020"),
    r2i = f_set("i", "r2"), exi = f_set("i", "extrap"), b16 = q("lz_B", "k2016", "max", 2),
    top = dderived("lower edge of the open histogram bin of the observed-to-true ratio", CH, "max(bin_lo) of metric aucinf_B", top, fnum(top, 1)),
    o16 = dint(CH, "pk_model=='k2016' & metric=='aucinf_B' & bin_lo >= 1.5", "n", "subjects in the open bin, rule B, k2016"),
    o20 = dint(CH, "pk_model=='k2020' & metric=='aucinf_B' & bin_lo >= 1.5", "n", "subjects in the open bin, rule B, k2020"))))
  deck_end()
}
