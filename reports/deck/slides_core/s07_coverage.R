# S7 ② 총노출 포착(지시 §2, 그림 4-1): 2016 모델 대표 대상자 3명(창 포착률 중앙값·5백분위·최솟값; 사전 등록 7c)의 음영 그림.
# 제목: 중앙값 대상자와 최솟값 대상자(패널과 같은 대상자)의 AUClast/AUCinf. 캡션: 모집단 외삽 비율 중앙값(비구획 대 참, pillar1_coverage_B0.csv).
slide_S7 <- function() {
  RSf <- "core_deck/rep_subjects.csv"; P1 <- "rationale/pillar1_coverage_B0.csv"
  deck_slide("S7", tag = "sim")
  med <- dv(RSf, "model=='k2016' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2016 median-coverage subject", scale = 100)
  mn <- dv(RSf, "model=='k2016' & role=='min'", "coverage_true", 1, "%", "window coverage (true), 2016 minimum-coverage subject", scale = 100)
  y0 <- core_title(tx("S7.title", list(med = med, min = mn)), tx("S7.kicker"))
  L <- DK$txt$S7$fig
  p <- core_shaded_panels("k2016", L)
  rr <- drange("core_deck/extrap_ratio_by_group.csv", "group=='reliable_iii'", "median", 1, "", "median NCA-to-true extrapolated area ratio, subjects with a reliable AUCinf, two models")
  premise(all(!rows(RSf, "model=='k2016'")$reliable_iii), "all three 2016 representative subjects fail set (iii) (notes)")
  cap <- tx("S7.caption", list(nsub = dint(P1, "model=='k2016' & group=='all'", "n_subjects", "virtual subjects, 2016"), rr = rr, nca = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, 2016"),
                               true = dv(P1, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, 2016"),
                               q3 = dv(RSf, "model=='k2016' & role=='min'", "extrap_ratio_nca_to_true", 2, "", "NCA-to-true extrapolated area ratio, minimum subject")))
  premise(row1(RSf, "model=='k2016' & role=='min'")$extrap_ratio_nca_to_true < 1 && row1(RSf, "model=='k2016' & role=='median'")$extrap_ratio_nca_to_true > 1, "median subject NCA over, minimum subject NCA under (caption)")
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(tx("S7.body"), capy - 0.06)
  deck_figure(p, "s7_shaded_k2016", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(RSf, "core_deck/rep_profiles.csv", "core_deck/rep_obs.csv"))
  r <- function(role, col, d, item, unit = "") dv(RSf, sprintf("model=='k2016' & role=='%s'", role), col, d, unit, item)
  deck_notes(tx("S7.notes", list(
    i1 = r("median", "id", 0, "subject id, median"), i2 = r("p05", "id", 0, "subject id, 5th percentile"), i3 = r("min", "id", 0, "subject id, minimum"),
    t1 = r("median", "study_day_tlast", 1, "study day of the observed tlast, median subject"), t2 = r("p05", "study_day_tlast", 1, "study day of tlast, 5th percentile subject"), t3 = r("min", "study_day_tlast", 1, "study day of tlast, minimum subject"),
    a1 = r("median", "extrap_area_nca", 1, "NCA extrapolated area, median subject"), b1 = r("median", "extrap_area_true", 1, "true extrapolated area, median subject"),
    a3 = r("min", "extrap_area_nca", 1, "NCA extrapolated area, minimum subject"), b3 = r("min", "extrap_area_true", 1, "true extrapolated area, minimum subject"),
    q3 = r("min", "extrap_ratio_nca_to_true", 2, "NCA-to-true extrapolated area ratio, minimum subject"), nsub = dint(P1, "model=='k2016' & group=='all'", "n_subjects", "virtual subjects, 2016"),
    r3 = r("min", "adj_r2", 2, "adjusted R-squared, minimum subject"), e3 = r("min", "pct_extrap_nca", 1, "NCA extrapolated share, minimum subject", "%"),
    cobs = r("min", "Clast", 2, "observed Clast, minimum subject (mg/L)"),
    ctrue = { s_ <- row1(RSf, "model=='k2016' & role=='min'"); o <- rows("core_deck/rep_obs.csv", sprintf("model=='k2016' & id==%d", s_$id))
      x <- o[abs(time_after_dose - s_$tlast) < 1e-9, conc_true]; premise(length(x) == 1, "true concentration at the observed tlast")
      dderived("true concentration at the observed tlast, minimum subject (mg/L)", "core_deck/rep_obs.csv", "model=='k2016' & id==<minimum subject> & time_after_dose==tlast :: conc_true", x, fnum(x, 2)) },
    cr = { s_ <- row1(RSf, "model=='k2016' & role=='min'"); o <- rows("core_deck/rep_obs.csv", sprintf("model=='k2016' & id==%d", s_$id)); x <- s_$Clast / o[abs(time_after_dose - s_$tlast) < 1e-9, conc_true]
      dderived("observed Clast / true concentration at tlast, minimum subject", "core_deck/rep_obs.csv", "Clast (rep_subjects.csv) / conc_true at tlast", x, fnum(x, 2)) },
    xr = { s_ <- row1(RSf, "model=='k2016' & role=='min'"); o <- rows("core_deck/rep_obs.csv", sprintf("model=='k2016' & id==%d", s_$id)); x <- o[abs(time_after_dose - s_$tlast) < 1e-9, conc_true] / s_$lambda_z / s_$extrap_area_true
      dderived("extrapolation from the true concentration with the NCA lambda-z / true extrapolated area, minimum subject", "core_deck/rep_obs.csv", "conc_true at tlast / lambda_z / extrap_area_true (rep_subjects.csv)", x, fnum(x, 2)) },
    r1 = r("median", "adj_r2", 2, "adjusted R-squared, median subject"), r2v = r("p05", "adj_r2", 2, "adjusted R-squared, 5th percentile subject"),
    ea = drange("core_deck/extrap_ratio_by_group.csv", "group=='estimable'", "median", 1, "", "median NCA-to-true extrapolated area ratio, lambda-z estimable, two models"),
    rr = rr, b1p = dv("core_deck/extrap_ratio_by_group.csv", "model=='k2016' & group=='reliable_iii'", "pct_below_1", 1, "%", "share with NCA below true, reliable subjects, 2016"),
    b2p = dv("core_deck/extrap_ratio_by_group.csv", "model=='k2020' & group=='reliable_iii'", "pct_below_1", 2, "%", "share with NCA below true, reliable subjects, 2020"),
    nca = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, 2016"),
    true = dv(P1, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, 2016"),
    wt = f_wt_range(), p95n = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_p95", 1, "%", "95th percentile NCA extrapolated share, 2016"),
    p95t = dv(P1, "model=='k2016' & group=='all'", "extrap_true_p95", 1, "%", "95th percentile true extrapolated share, 2016"))))
  deck_end()
}
