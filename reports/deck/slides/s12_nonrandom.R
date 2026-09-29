# S12 논거 ① 비무작위 탈락: 탈락자의 참 AUC0-inf가 낮고(세트별 탈락자 대 유지자 기하평균비), 무거운 층이 더 탈락하며,
# 규칙 A 분석군은 동일 제품 시험에서도 arm 간 층 구성이 달라진다(AUC0-last 분석군은 무작위배정 균형 유지). 시험 모집단만(건강인 60~90 kg, B0).
# 수치: results/trialpop/tp_characteristics.csv(개인 수준 탈락자 대 유지자), tp_strata_individual.csv(층별 탈락),
#       tp_strata_composition.csv(동일 제품 S00 시험 2,000회의 arm 간 무거운 층 비율 차이 요약).
# 그림: 커밋된 시험별 인원 results/trialpop/trialpop_counts_<model>.csv.gz에서 scripts/55와 같은 식으로 시험마다 차이를 다시 계산해
#       정수 %p 중심 구간(round)으로 그린다. 요약 파일의 백분위수·5%p 초과 비율과 results/deck_inputs/strata_armdiff_hist.csv(floor 구간)를 코드에서 대조한다.
# 주의: 인원 수 기반 결과라 분석 모형(M0/M1)과 무관하다(M1로 표시하지 않는다). "19~21% 낮다"는 세트 (i)만의 정수 반올림이므로 GMR과 정확한 %를 쓴다.
s12_armdiff <- function() {
  SX <- c(auclast = "n_auclast", i = "n_i", iii = "n_iii")
  comp <- function(d, col) { w <- dcast(d, trial ~ stratum, value.var = col); setnames(w, c("1", "2"), c("a1", "a2")); w[, .(trial, h = 100 * a2 / (a1 + a2))] }
  rbindlist(lapply(c("k2016", "k2020"), function(m) {
    d <- rows(sprintf("trialpop/trialpop_counts_%s.csv.gz", m), "scenario %in% c('REF', 'S00')")
    premise(all(d[scenario == "REF", arm] == "R") && all(d[scenario == "S00", arm] == "T"), "counts file: REF is the reference arm, S00 the identical-product test arm")
    rbindlist(lapply(names(SX), function(k) {
      z <- comp(d[scenario == "S00"], SX[[k]])[comp(d[scenario == "REF"], SX[[k]]), on = "trial", .(trial, armdiff = h - i.h)]
      z[, `:=`(pk_model = m, analysis_set = k)]
    }))
  }))
}
slide_S12 <- function() {
  TCH <- "trialpop/tp_characteristics.csv"; TSI <- "trialpop/tp_strata_individual.csv"; TSC <- "trialpop/tp_strata_composition.csv"
  TPF <- "trialpop/tp_failure_by_set.csv"; HIST <- "deck_inputs/strata_armdiff_hist.csv"
  GZ <- c("trialpop/trialpop_counts_k2016.csv.gz", "trialpop/trialpop_counts_k2020.csv.gz")
  deck_slide("S12", tag = "sim")
  L <- DK$txt$S12

  # ---- 전제(문장이 기대는 자료의 부호·크기) ----
  ch <- rows(TCH); si <- rows(TSI); sc <- rows(TSC, "scenario=='S00'")
  premise(all(ch[set %in% c("i", "iii"), true_aucinf_gmr_hi] < 1) && all(ch[set %in% c("i", "iii"), auclast_gmr_hi] < 1), "failing subjects have lower true AUC0-inf and AUC0-last (upper 95% limit below 1), sets i and iii")
  premise(all(ch[set == "i", wt_diff_lo] > 0), "failing subjects are heavier, set i (lower 95% limit above 0)")
  premise(all(si[set %in% c("i", "iii"), diff_lo] > 0), "heavier stratum fails more often, sets i and iii (Newcombe lower limit above 0)")
  premise(all(ch[set == "iii", true_aucinf_gmr] > ch[set == "i", true_aucinf_gmr]) && all(ch[set == "iv", true_aucinf_gmr] > ch[set == "i", true_aucinf_gmr]) &&
            all(ch[set %in% c("iii", "iv"), true_aucinf_gmr] > max(ch[set %in% c("i", "ii"), true_aucinf_gmr])), "sets iii and iv ratios are closer to 1 than sets i and ii (smaller difference), both models")
  premise(all(ch[set == "iii", wt_diff_kg] < ch[set == "i", wt_diff_kg]), "set iii weight difference is smaller than set i, both models")
  n_arm_num <- as.numeric(.read("config/trial_design.yaml")$n_per_arm)
  premise(all(abs(sc$rand_armdiff_abs_max - 100 / n_arm_num) < 1e-9), "randomization bound equals one subject of n_per_arm")

  # ---- 시험별 arm 간 차이(커밋된 인원 파일에서 다시 계산) + 요약·deck_inputs와 대조 ----
  AD <- s12_armdiff()
  gtc <- grep("^armdiff_abs_gt[0-9]+_pct$", names(sc), value = TRUE); premise(length(gtc) == 1, "one threshold share column in tp_strata_composition.csv")
  thr_v <- as.numeric(sub("^armdiff_abs_gt([0-9]+)_pct$", "\\1", gtc))                    # 기준(%p)은 열 이름에서 읽는다
  for (m in c("k2016", "k2020")) for (k in c("auclast", "i", "iii")) {
    v <- AD[pk_model == m & analysis_set == k, armdiff]; r <- sc[pk_model == m & analysis_set == k]
    premise(nrow(r) == 1 && length(v) == r$n_trials && abs(median(v) - r$armdiff_median) < 1e-9 && abs(quantile(v, 0.05, names = FALSE) - r$armdiff_p05) < 1e-9 &&
              abs(quantile(v, 0.95, names = FALSE) - r$armdiff_p95) < 1e-9 && abs(100 * mean(abs(v) > thr_v) - r[[gtc]]) < 1e-9,
            sprintf("recomputed arm differences reproduce tp_strata_composition (%s, %s)", m, k))
  }
  hf <- AD[, .(n = .N), by = .(pk_model, analysis_set, bin = floor(armdiff))][, pct := 100 * n / sum(n), by = .(pk_model, analysis_set)]
  hc <- rows(HIST)[hf, on = c("pk_model", "analysis_set", "bin")]
  premise(nrow(hc) == nrow(rows(HIST)) && all(hc$n == hc$i.n), "recomputed floor-binned histogram equals results/deck_inputs/strata_armdiff_hist.csv")
  la <- AD[analysis_set == "auclast", max(abs(armdiff))]
  premise(la <= max(sc$rand_armdiff_abs_max) + 1e-9, "AUC0-last analysis set stays within the randomization bound in every trial")
  premise(all(sc[analysis_set %in% c("i", "iii"), armdiff_p95] > max(sc$rand_armdiff_abs_max)), "rule A analysis sets exceed the randomization bound")

  # ---- 제목 ----
  gmr_i <- drange(TCH, "set=='i'", "true_aucinf_gmr", 3, "", "failing versus retained, set i, true_aucinf_gmr")
  deck_kicker(tx("S12.kicker")); deck_title(tx("S12.title", list(gmr = gmr_i)))

  # ---- 왼쪽: 수치 카드 두 개 ----
  XL <- GEO$ML; WL <- 5.6
  lo_pct <- 100 * (1 - ch[set == "i", true_aucinf_gmr]); pl <- rng_fmt(min(lo_pct), max(lo_pct), 1, "%")
  pl <- dderived("true AUC0-inf of failing subjects lower than retained, set i, two models (percent)", TCH, "range over 2 rows [set=='i'] :: 100 x (1 - true_aucinf_gmr)", range(lo_pct), pl)
  g3 <- drange(TCH, "set=='iii'", "true_aucinf_gmr", 3, "", "failing versus retained, set iii, true_aucinf_gmr")
  C1H <- 1.77; GAP <- 0.08
  VS <- 36                                                                                  # 큰 수치 36pt: 세 줄 설명 아래 여백 확보(검토: 카드 아래 가장자리에 붙음)
  deck_stat(gmr_i, tx("S12.card_gmr", list(pl = pl, g3 = g3)), c(XL, GEO$BODY_TOP, WL, C1H), value_size = VS)
  sd_i <- drange(TSI, "set=='i'", "diff_pp", 1, "", "stratum difference, set i, diff_pp")
  sd_iii <- drange(TSI, "set=='iii'", "diff_pp", 1, "", "stratum difference, set iii, diff_pp")
  wd <- drange(TCH, "set=='i'", "wt_diff_kg", 2, "", "failing versus retained, set i, wt_diff_kg")
  C2Y <- GEO$BODY_TOP + C1H + GAP; C2H <- 1.77
  deck_stat(paste0(sd_i, "%p"), tx("S12.card_strata", list(split = f_split(), sd3 = sd_iii, wd = wd)), c(XL, C2Y, WL, C2H), color = PAL$ink, bg = PAL$tint_grey, value_size = VS)

  # ---- 왼쪽 아래: 분석군별 arm 간 차이 요약(두 모델) ----
  thr <- dderived(sprintf("threshold in column name %s (points)", gtc), TSC, sprintf("column name %s", gtc), thr_v, fnum(thr_v, 0))
  SETS <- c("auclast", "i", "iii")
  w_ <- function(k) sprintf("scenario=='S00' & analysis_set=='%s'", k)
  df <- data.frame(a = unlist(L$table$rows[SETS]),
                   b = vapply(SETS, function(k) paste0(dspan(TSC, w_(k), "armdiff_p05", "armdiff_p95", 1, "", sprintf("armdiff, S00, %s: 5th to 95th percentile, two models", k)), "%p"), ""),
                   c = vapply(SETS, function(k) drange(TSC, w_(k), gtc, 1, "%", sprintf("share of trials with a between-arm heavier-stratum difference above 5 points, S00, %s", k)), ""),
                   stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S12.table.head", list(thr = thr))
  TY <- C2Y + C2H + GAP
  deck_table(df, box = c(XL, TY, WL, GEO$BODY_BOTTOM - TY), widths = c(1.75, 2.25, 1.3), size = 12, highlight = 3, label = "table_balance")

  # ---- 오른쪽: 동일 제품 시험의 arm 간 무거운 층 비율 차이 분포(분석군 × 두 모델) ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  ntr <- dint(TSC, "pk_model=='k2016' & scenario=='S00' & analysis_set=='auclast'", "n_trials", "regenerated trials per scenario")
  premise(all(sc$n_trials == sc$n_trials[1]), "same number of trials for both models and all analysis sets")
  h <- AD[, .(n = .N), by = .(pk_model, analysis_set, bin = floor(armdiff + 0.5))][, pct := 100 * n / sum(n), by = .(pk_model, analysis_set)]
  SL <- unlist(L$fig$sets[SETS])
  h[, set := factor(SL[analysis_set], levels = SL)]; h[, model := factor(model_lab()[pk_model], levels = model_lab())]
  xr <- range(h$bin) + c(-1, 1); thr_n <- thr_v
  p <- ggplot(h, aes(x = bin, y = pct, fill = model)) +
    annotate("rect", xmin = -Inf, xmax = -thr_n, ymin = -Inf, ymax = Inf, fill = PAL$tint_grey) +
    annotate("rect", xmin = thr_n, xmax = Inf, ymin = -Inf, ymax = Inf, fill = PAL$tint_grey) +
    geom_vline(xintercept = c(-thr_n, thr_n), colour = PAL$ink2, linewidth = 0.5, linetype = "22") +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.9), width = 0.86, colour = NA) +
    facet_wrap(~set, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = unname(MODEL_COL)) +
    scale_x_continuous(breaks = seq(-20, 20, by = 5), limits = xr, expand = expansion(add = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)), breaks = function(l) pretty(c(0, l[2]), n = 3)) +
    labs(x = L$fig$xlab, y = NULL, subtitle = fill(L$fig$ylab, list(n = fint(sc$n_trials[1]), thr = thr))) + theme_deck(13) +          # 회색 영역 설명은 부제에(패널 구석 표지 대신)
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          panel.grid.major.x = element_blank(), panel.spacing = grid::unit(8, "pt"), strip.text = element_text(hjust = 0, size = 13, face = "bold", colour = PAL$ink),
          plot.subtitle = element_text(colour = PAL$ink2, size = 13, margin = margin(0, 0, 2, 0)))
  FH <- 4.1
  deck_figure(p, "s12_strata_armdiff", c(XR, GEO$BODY_TOP, WR, FH), src = c(GZ, HIST, TSC))
  rb <- dext(TSC, "analysis_set=='auclast'", "rand_armdiff_abs_max", max, 1, "", "largest between-arm difference in the heavier-stratum share at randomization (points)")
  deck_text(tx("S12.takeaway", list(rb = rb)), c(XR, GEO$BODY_TOP + FH + 0.05, WR, GEO$BODY_BOTTOM - GEO$BODY_TOP - FH - 0.05), size = 16, label = "text_takeaway")

  # ---- 노트 ----
  q <- function(where, col, d, unit, item) dv(TSI, where, col, d, unit, item)
  wci <- function(m) dci(TCH, sprintf("pk_model=='%s' & set=='i'", m), "wt_diff_kg", "wt_diff_lo", "wt_diff_hi", 2, " kg", sprintf("failing minus retained body weight with Welch 95%% CI, set i, %s", m))
  ci <- function(m) dci(TSI, sprintf("pk_model=='%s' & set=='i'", m), "diff_pp", "diff_lo", "diff_hi", 1, "%p", sprintf("stratum difference with Newcombe 95%% CI, set i, %s", m))
  deck_notes(tx("S12.notes", list(
    wt = f_wt_range(), dose = f_dose(), n = dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"), ntr = ntr,
    g16 = dv(TCH, "pk_model=='k2016' & set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf GMR failing to retained, set i, k2016"),
    g20 = dv(TCH, "pk_model=='k2020' & set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf GMR failing to retained, set i, k2020"),
    pl = pl, g3 = g3,
    g2 = drange(TCH, "set=='ii'", "true_aucinf_gmr", 3, "", "failing versus retained, set ii, true_aucinf_gmr, two models"),
    g4 = drange(TCH, "set=='iv'", "true_aucinf_gmr", 3, "", "failing versus retained, set iv, true_aucinf_gmr, two models"),
    w16 = wci("k2016"), w20 = wci("k2020"), w3 = drange(TCH, "set=='iii'", "wt_diff_kg", 2, " kg", "failing minus retained body weight, set iii, two models"),
    al = drange(TCH, "set=='i'", "auclast_gmr", 3, "", "AUC0-last GMR failing to retained, set i, two models"),
    l16 = q("pk_model=='k2016' & set=='i'", "fail_light_pct", 1, "%", "failing, lighter stratum, set i, k2016"),
    h16 = q("pk_model=='k2016' & set=='i'", "fail_heavy_pct", 1, "%", "failing, heavier stratum, set i, k2016"),
    l20 = q("pk_model=='k2020' & set=='i'", "fail_light_pct", 1, "%", "failing, lighter stratum, set i, k2020"),
    h20 = q("pk_model=='k2020' & set=='i'", "fail_heavy_pct", 1, "%", "failing, heavier stratum, set i, k2020"),
    c16 = ci("k2016"), c20 = ci("k2020"), rb = rb, n_arm = f_n_arm(),
    med = drange(TSC, "scenario=='S00' & analysis_set %in% c('i','iii')", "armdiff_median", 2, "", "armdiff median, S00, rule A sets i and iii, two models"))))
  deck_end()
}
