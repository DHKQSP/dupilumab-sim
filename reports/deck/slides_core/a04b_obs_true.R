# 별첨 A4b ② 세부: 관측 대 참 비. 시험 모집단(건강인, 체중 층화 범위, 현행 채혈 B0), 모델당 20,000명.
#  왼쪽 그림: 대상자 수준 분포(results/deck_inputs/coverage_hist.csv, 0.02 간격, 1.5 이상은 한 구간; 2016 모델, 2020 모델 나란히, 로그 세로축):
#        AUClast/AUCinf(참값)(참 AUClast ÷ 참 AUCinf) 옅은 파랑 막대, 관측 AUClast ÷ 참 AUCinf 파랑 계단선(전체),
#        NCA AUCinf ÷ 참 AUCinf 규칙 B(λz 산출 가능 전원) 주황 실선, 규칙 A 세트 (i)(기준 충족자만) 주황 점선. 빈 구간에서 선을 끊는다.
#  오른쪽 표 두 개(모델마다, 위아래): results/trialpop/tp_coverage_individual.csv(대상자 수, 중앙값, 5~95백분위, 최대). 최솟값은 노트.
#  세트 (iii)(신뢰할 수 있는 AUCinf)의 관측 대 참 비는 결과 파일에 없어 세트 (i)로 적는다(캡션, 노트). 외삽 비율 중앙값은 pillar1_coverage_B0.csv(S7과 같은 열).
# 전제와 분포 자료 처리는 결과보고 덱 S16(slides/s16_coverage.R)을 따른다.
A4B_TCV <- "trialpop/tp_coverage_individual.csv"
# 두 모델 값(2016, 2020 순서의 문자 벡터)
a4b_two <- function(met, col, d, what, big = FALSE)
  vapply(c("k2016", "k2020"), function(m) { w <- sprintf("metric=='%s' & pk_model=='%s'", met, m); it <- sprintf("%s, %s, %s", what, col, m)
    if (big) dint(A4B_TCV, w, col, it) else dv(A4B_TCV, w, col, d, "", it) }, "", USE.NAMES = FALSE)
slide_A4b <- function() {
  TCV <- A4B_TCV; CH <- "deck_inputs/coverage_hist.csv"; P1 <- "rationale/pillar1_coverage_B0.csv"
  MW <- "window coverage (true AUC0-tlast / true AUC0-inf)"
  MET <- c(win = MW, last = "observed-to-true, AUClast (all subjects)", B = "observed-to-true, AUCinf rule B (lambda-z estimable)",
           A = "observed-to-true, AUCinf rule A (meeting set (i))")
  K <- c("win", "last", "B", "A"); M <- c("k2016", "k2020")
  premise(all(MET %in% rows(TCV)$metric), "metric labels present in tp_coverage_individual.csv")
  premise(nrow(rows(TCV, "grepl('set \\\\(iii\\\\)', metric)")) == 0, "no set (iii) observed-to-true metric in the file (set (i) is shown; caption, notes)")
  wm <- function(k, m) sprintf("metric=='%s' & pk_model=='%s'", MET[[k]], m)
  deck_slide("A4b", tag = "sim")
  L <- DK$txt$A4b

  # ---- 전제: 제목과 본문의 방향 ----
  x <- rows(TCV)
  for (m in M) {
    r0 <- x[metric == MET[["last"]] & pk_model == m]
    for (k in c("B", "A")) { r1 <- x[metric == MET[[k]] & pk_model == m]
      premise(r1$p95 > r0$p95 && r1$max > r0$max, sprintf("NCA AUCinf %s has a longer upper tail than observed AUClast (%s): title", k, m)) }
    premise(x[metric == MET[["B"]] & pk_model == m, n] < x[metric == MET[["last"]] & pk_model == m, n], sprintf("rule B covers fewer subjects than AUClast (%s): body", m))
    premise(r0$median < 1, sprintf("observed AUClast median below the true AUCinf (%s): body", m))
  }
  premise(all(rows(P1, "group=='all'")[, extrap_nca_median > extrap_true_median]), "NCA extrapolated share above the true share at the median, both models (body)")
  premise(all(as.numeric(sub(" .*", "", rows(P1, "group=='all'")$lambda_ok_pct_ci)) < 100), "NCA extrapolation median over lambda-z-estimable subjects, true over all (caption denominators)")
  bmax <- x[metric == MET[["B"]], .(pk_model, max)]
  premise(bmax[which.max(max), pk_model] == "k2016" && all(x[metric == MET[["A"]], max] < max(bmax$max)), "largest NCA-to-true ratio: rule B, 2016 model (title; rule A lower)")
  f <- list(med = drange(TCV, sprintf("metric=='%s'", MW), "median", 0, "%", "window coverage, median, two models", scale = 100),
            mx = dv(TCV, wm("B", "k2016"), "max", 1, "", "largest NCA AUCinf (rule B) / true AUCinf, 2016 model (largest of both models, premise)"))
  y0 <- core_title(tx("A4b.title", f), tx("A4b.kicker"))

  # ---- 캡션(전체 폭, 아래) ----
  h <- copy(rows(CH, "metric %in% c('window','auclast','aucinf_B','aucinf_A')"))
  bw <- unique(round(diff(sort(unique(h$bin_lo))), 6)); premise(length(bw) == 1, "one bin width")
  top <- max(h$bin_lo); premise(abs(top - 1.5) < 1e-9 && all(h[bin_lo >= top, metric] == "aucinf_B"), "open bin from 1.5 holds rule B only")
  bwp <- dderived("histogram bin width", CH, "unique(diff(sort(unique(bin_lo))))", bw, fnum(bw, 2))
  cap <- tx("A4b.caption", list(r2 = f_set("i", "r2"), ex = f_set("i", "extrap"), bw = bwp,
                                n = drange(P1, "group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, two models"),
                                t = drange(P1, "group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, two models")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)

  # ---- 오른쪽: 표 두 개(모델마다, 위아래) ----
  TB <- L$table; WT <- c(2.22, 0.84, 0.71, 1.36, 0.6); XT <- GEO$W - GEO$MR - sum(WT); GAPV <- 0.16
  TH <- (capy - 0.12 - y0 - GAPV) / 2
  for (i in seq_along(M)) {
    m <- M[i]; w1 <- function(k) sprintf("metric=='%s' & pk_model=='%s'", MET[[k]], m); it <- function(k, col) sprintf("%s, %s, %s", k, col, m)
    df <- data.frame(a = unlist(TB$rows[K]),
                     b = vapply(K, function(k) dint(TCV, w1(k), "n", it(k, "n")), ""),
                     c = vapply(K, function(k) dv(TCV, w1(k), "median", 3, "", it(k, "median")), ""),
                     d = vapply(K, function(k) dspan(TCV, w1(k), "p05", "p95", 3, "", it(k, "p05 to p95")), ""),
                     e = vapply(K, function(k) dv(TCV, w1(k), "max", 2, "", it(k, "max")), ""), stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- c(DK$txt$common$models_short[[m]], unlist(TB$head))
    deck_table(df, box = c(XT, y0 + (i - 1) * (TH + GAPV), sum(WT), TH), widths = WT, size = 14, label = sprintf("table_ratio_%s", m))
  }

  # ---- 왼쪽 그림: 분포(두 모델 나란히, 로그 세로축). 계단선은 비어 있는 구간에서 끊는다(0으로 이어 그리지 않음) ----
  premise(all(abs(h[, sum(pct), by = .(pk_model, metric)]$V1 - 100) < 1e-6), "histograms sum to 100% per model and metric")
  premise(all(h[metric == "auclast", sum(n), by = pk_model]$V1 == vapply(M, function(m) row1(TCV, wm("last", m))$n, 1)), "AUClast histogram covers all subjects")
  premise(all(h[metric == "aucinf_B", sum(n), by = pk_model]$V1 == vapply(M, function(m) row1(TCV, wm("B", m))$n, 1)) &&
          all(h[metric == "aucinf_A", sum(n), by = pk_model]$V1 == vapply(M, function(m) row1(TCV, wm("A", m))$n, 1)), "NCA histograms cover the rule B and rule A (set (i)) subjects")
  premise(all(h$n > 0), "histogram file lists only non-empty bins (empty bins are gaps)")
  ML <- DK$txt$common$models_short; F <- L$fig
  h[, mod := factor(unlist(ML[pk_model]), levels = unlist(ML))]
  SER <- c(auclast = F$s_last, aucinf_B = F$s_B, aucinf_A = F$s_A)
  win <- h[metric == "window"]
  # 계단 경로: 연속한 구간끼리만 잇는다(구간 번호 차이가 1이 아니면 새 경로)
  st <- h[metric != "window" & bin_lo < top][order(pk_model, metric, bin_lo)]
  st[, k := round(bin_lo / bw)][, run := cumsum(c(1, diff(k) != 1)), by = .(pk_model, metric)]
  stp <- st[, .(x = c(rbind(bin_lo, bin_lo + bw)), y = rep(pct, each = 2)), by = .(pk_model, mod, metric, run)][, grp := paste(pk_model, metric, run)]
  stp[, ser := factor(SER[metric], levels = SER)]
  iso <- st[, if (.N == 1) .(x = bin_lo + bw / 2, y = pct), by = .(pk_model, mod, metric, run)][, ser := factor(SER[metric], levels = SER)]
  ox <- top + bw * 3
  ov <- h[bin_lo >= top][, x := ox][, ser := factor(SER[metric], levels = SER)]
  ov[, lab := vapply(seq_len(.N), function(i) fill(F$ovf, list(n = fint(n[i]))), "")]
  YF <- 0.003
  p <- ggplot() +
    geom_rect(data = win, aes(xmin = bin_lo, xmax = bin_lo + bw, ymin = YF, ymax = pct, fill = F$s_win), colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 1, colour = PAL$ink2, linewidth = 0.5) +
    geom_path(data = stp, aes(x, y, colour = ser, linetype = ser, group = grp), linewidth = 0.8) +
    geom_point(data = iso, aes(x, y, colour = ser), size = 1.6, show.legend = FALSE) +
    geom_point(data = ov, aes(x, pct), colour = PAL$orange, shape = 15, size = 2.6) +
    geom_text(data = ov, aes(x, pct * 4, label = lab), hjust = 0.5, vjust = 0, size = PT(14), family = FONT, colour = PAL$ink) +
    facet_wrap(~ mod, nrow = 1) + coord_cartesian(clip = "off") +
    scale_fill_manual(values = setNames(BLUE_LIGHT, F$s_win), name = NULL) +
    scale_colour_manual(values = setNames(c(PAL$blue, PAL$orange, PAL$orange), SER), name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "solid", "22"), SER), name = NULL) +
    scale_y_log10(limits = c(YF, 150), breaks = c(0.01, 1, 100), labels = function(v) paste0(formatC(v, format = "fg"), "%"), expand = expansion(0)) +
    scale_x_continuous(limits = c(0.56, ox + bw * 2.5), breaks = c(0.6, 0.8, 1, 1.2, ox),
                       labels = c(fnum(c(0.6, 0.8, 1, 1.2), 1), fill(F$ovf_tick, list(v = fnum(top, 1)))), expand = expansion(0)) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, nrow = 1), linetype = guide_legend(order = 2, nrow = 1)) +
    labs(x = F$xlab, y = NULL) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(1.6, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box = "vertical", legend.box.just = "left", legend.box.spacing = grid::unit(2, "pt"), legend.spacing.x = grid::unit(3, "pt"), legend.spacing.y = grid::unit(0, "pt"),
          panel.grid.minor = element_blank(), panel.spacing.x = grid::unit(18, "pt"), plot.margin = margin(2, 8, 2, 2))
  deck_figure(p, "a4b_obs_true_hist", c(GEO$ML, y0, XT - 0.25 - GEO$ML, capy - 0.12 - y0), src = c(CH, TCV))

  # ---- 노트 ----
  q <- function(k, m, col, d, what = col) dv(TCV, wm(k, m), col, d, "", sprintf("%s, %s, %s", MET[[k]], what, m))
  TF <- "trialpop/tp_failure_by_set.csv"
  deck_notes(tx("A4b.notes", list(
    wt = f_wt_range(), dose = f_dose(), n = dint(TCV, wm("win", "k2016"), "n", "subjects per model (notes)"),
    b16 = q("B", "k2016", "median", 3), b20 = q("B", "k2020", "median", 3),
    sl16 = q("last", "k2016", "sd_log", 3, "sd of log"), sb16 = q("B", "k2016", "sd_log", 3, "sd of log"), sa16 = q("A", "k2016", "sd_log", 3, "sd of log"),
    sl20 = q("last", "k2020", "sd_log", 3, "sd of log"), sb20 = q("B", "k2020", "sd_log", 3, "sd of log"), sa20 = q("A", "k2020", "sd_log", 3, "sd of log"),
    mn = local({ v <- vapply(K, function(k) paste(a4b_two(MET[[k]], "min", 3, k), collapse = " / "), ""); paste(sprintf("%s %s", unlist(L$table$rows[K]), v), collapse = "; ") }),
    fx16 = dv(TCV, "metric=='extrapolation factor AUCinf / AUClast (lambda-z estimable)' & pk_model=='k2016'", "max", 2, "", "extrapolation factor AUCinf/AUClast, max, k2016"),
    fx20 = dv(TCV, "metric=='extrapolation factor AUCinf / AUClast (lambda-z estimable)' & pk_model=='k2020'", "max", 2, "", "extrapolation factor AUCinf/AUClast, max, k2020"),
    top = dderived("lower edge of the open histogram bin", CH, "max(bin_lo) of metric aucinf_B", top, fnum(top, 1)), bw = bwp,
    o16 = dint(CH, "pk_model=='k2016' & metric=='aucinf_B' & bin_lo >= 1.5", "n", "subjects in the open bin, rule B, k2016"),
    o20 = dint(CH, "pk_model=='k2020' & metric=='aucinf_B' & bin_lo >= 1.5", "n", "subjects in the open bin, rule B, k2020"),
    r2i = f_set("i", "r2"), exi = f_set("i", "extrap"), r2 = f_set("iii", "r2"),
    f3 = drange(TF, "set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), two models"))))
  deck_end()
}
