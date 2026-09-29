# S7 ② 총노출 포착(지시 §2, 그림 4-1): 2016 모델 대표 대상자 3명(창 포착률 중앙값·5백분위·최솟값; 사전 등록 7c)의 음영 그림.
# 제목: 중앙값 대상자의 AUClast/AUCinf(반올림)와 최솟값 대상자의 하한(내림). 캡션: 모집단 외삽 비율 중앙값(비구획 대 참, pillar1_coverage_B0.csv).
slide_S7 <- function() {
  RSf <- "core_deck/rep_subjects.csv"; P1 <- "rationale/pillar1_coverage_B0.csv"
  deck_slide("S7", tag = "sim")
  med <- dv(RSf, "model=='k2016' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2016 median-coverage subject", scale = 100)
  mn <- { x <- 100 * row1(RSf, "model=='k2016' & role=='min'")$coverage_true
    dderived("window coverage (true), 2016 minimum-coverage subject, rounded down", RSf, "model=='k2016' & role=='min' :: floor(coverage_true x 100)", x, paste0(fnum(fl(x, 0), 0), "%")) }
  y0 <- core_title(tx("S7.title", list(med = med, min = mn)), tx("S7.kicker"))
  L <- DK$txt$S7$fig
  p <- core_shaded_panels("k2016", L)
  cap <- tx("S7.caption", list(nca = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, 2016"),
                               true = dv(P1, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, 2016")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(tx("S7.body", list(nsub = dint(P1, "model=='k2016' & group=='all'", "n_subjects", "virtual subjects, 2016"))), capy - 0.06)
  deck_figure(p, "s7_shaded_k2016", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(RSf, "core_deck/rep_profiles.csv", "core_deck/rep_obs.csv"))
  r <- function(role, col, d, item, unit = "") dv(RSf, sprintf("model=='k2016' & role=='%s'", role), col, d, unit, item)
  deck_notes(tx("S7.notes", list(
    i1 = r("median", "id", 0, "subject id, median"), i2 = r("p05", "id", 0, "subject id, 5th percentile"), i3 = r("min", "id", 0, "subject id, minimum"),
    t1 = r("median", "study_day_tlast", 1, "study day of the observed tlast, median subject"), t2 = r("p05", "study_day_tlast", 1, "study day of tlast, 5th percentile subject"), t3 = r("min", "study_day_tlast", 1, "study day of tlast, minimum subject"),
    a1 = r("median", "extrap_area_nca", 1, "NCA extrapolated area, median subject"), b1 = r("median", "extrap_area_true", 1, "true extrapolated area, median subject"),
    a3 = r("min", "extrap_area_nca", 1, "NCA extrapolated area, minimum subject"), b3 = r("min", "extrap_area_true", 1, "true extrapolated area, minimum subject"),
    q3 = r("min", "extrap_ratio_nca_to_true", 2, "NCA-to-true extrapolated area ratio, minimum subject"),
    r3 = r("min", "adj_r2", 2, "adjusted R-squared, minimum subject"), e3 = r("min", "pct_extrap_nca", 1, "NCA extrapolated share, minimum subject", "%"),
    wt = f_wt_range(), p95n = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_p95", 1, "%", "95th percentile NCA extrapolated share, 2016"),
    p95t = dv(P1, "model=='k2016' & group=='all'", "extrap_true_p95", 1, "%", "95th percentile true extrapolated share, 2016"))))
  deck_end()
}
