# 별첨 A4a ② 세부(S6의 2016 모델판, v1.2): 2016 모델(민감도 모델) 대표 대상자 3명(AUClast/AUCinf(참값) 중앙값·5백분위·최솟값; 사전 등록 7c)의 음영 그림.
# 그림은 S6과 같은 core_shaded_panels("k2016", L). v1.1 본문 S7(2016 모델)의 설명을 옮겼다: 중앙값·5백분위 대상자의 비구획 외삽 면적 배수,
# 최솟값 대상자의 비구획 외삽이 오히려 작은 이유(관측 Clast가 측정 잔차로 참 농도보다 낮음), 모집단 배수(results/core_deck/extrap_ratio_by_group.csv, 7f).
# 제목: 중앙값 대상자의 AUClast/AUCinf(반올림)와 최솟값 대상자(소수 한 자리, 그림과 같은 대상자). 캡션: 외삽 비율 중앙값의 모수(비구획 = λz 산출 대상자, 참 = 전원).
slide_A4a <- function() {
  RSf <- "core_deck/rep_subjects.csv"; ROf <- "core_deck/rep_obs.csv"; P1 <- "rationale/pillar1_coverage_B0.csv"; ER <- "core_deck/extrap_ratio_by_group.csv"
  deck_slide("A4a", tag = "sim")
  S <- rows(RSf, "model=='k2016'")
  premise(nrow(S) == 3 && setequal(S$role, c("median", "p05", "min")), "three representative subjects of the 2016 model")
  premise(all(S$lambda_ok) && all(!S$reliable_iii), "2016 model: lambda-z estimable in all three subjects, none meets set (iii) (body)")
  premise(all(S[role != "min", extrap_ratio_nca_to_true] > 1) && S[role == "min", extrap_ratio_nca_to_true] < 1,
          "2016 model: NCA extrapolated area above the true one for the median and 5th percentile subjects, below it for the minimum subject (body)")
  premise(row1(P1, "model=='k2016' & group=='all'")[, extrap_nca_median > extrap_true_median], "2016 model: NCA extrapolated share above the true share at the median (caption)")
  premise(as.numeric(sub(" .*", "", row1(P1, "model=='k2016' & group=='all'")$lambda_ok_pct_ci)) < 100,
          "not every subject has lambda-z: the NCA median is over lambda-z-estimable subjects, the true median over all subjects (caption denominators; scripts/15 median(pct_extrap, na.rm = TRUE))")
  med <- dv(RSf, "model=='k2016' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2016 median-coverage subject", scale = 100)
  mn <- dv(RSf, "model=='k2016' & role=='min'", "coverage_true", 1, "%", "window coverage (true), 2016 minimum subject, one decimal", scale = 100)   # 그림의 최솟값 대상자 값 그대로
  y0 <- core_title(tx("A4a.title", list(med = med, min = mn)), tx("A4a.kicker"))
  L <- DK$txt$A4a$fig
  p <- core_shaded_panels("k2016", L)
  r <- function(role, col, d, item, unit = "") dv(RSf, sprintf("model=='k2016' & role=='%s'", role), col, d, unit, item)

  # 최솟값 대상자: 관측 Clast 대 같은 시점 참 농도(측정 잔차), 참 농도에서 같은 λz로 외삽한 면적 대 참 외삽 면적
  s3 <- row1(RSf, "model=='k2016' & role=='min'"); o3 <- rows(ROf, sprintf("model=='k2016' & id==%d", s3$id))
  ct <- o3[abs(time_after_dose - s3$tlast) < 1e-9, conc_true]; premise(length(ct) == 1, "true concentration at the observed tlast of the minimum subject")
  premise(s3$Clast / ct < 1, "minimum subject: observed Clast below the true concentration at tlast (body)")
  cr <- dderived("observed Clast / true concentration at tlast, 2016 minimum subject", ROf, "Clast (rep_subjects.csv) / conc_true at tlast", s3$Clast / ct, fnum(s3$Clast / ct, 2))

  cap <- tx("A4a.caption", list(nca = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, 2016"),
                                true = dv(P1, "model=='k2016' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, 2016"),
                                nsub = dint(P1, "model=='k2016' & group=='all'", "n_subjects", "virtual subjects, 2016"),
                                ea = drange(ER, "group=='estimable'", "median", 1, "", "median NCA-to-true extrapolated area ratio, lambda-z estimable, two models"),
                                rr = drange(ER, "group=='reliable_iii'", "median", 1, "", "median NCA-to-true extrapolated area ratio, subjects with a reliable AUCinf, two models")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  body <- tx("A4a.body", list(q1 = r("median", "extrap_ratio_nca_to_true", 1, "NCA-to-true extrapolated area ratio, median subject"),
                              q2 = r("p05", "extrap_ratio_nca_to_true", 1, "NCA-to-true extrapolated area ratio, 5th percentile subject"),
                              q3 = r("min", "extrap_ratio_nca_to_true", 2, "NCA-to-true extrapolated area ratio, minimum subject"), cr = cr))
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "a4a_shaded_k2016", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(RSf, "core_deck/rep_profiles.csv", ROf))

  deck_notes(tx("A4a.notes", list(
    i1 = r("median", "id", 0, "subject id, median"), i2 = r("p05", "id", 0, "subject id, 5th percentile"), i3 = r("min", "id", 0, "subject id, minimum"),
    t1 = r("median", "study_day_tlast", 1, "study day of the observed tlast, median subject"), t2 = r("p05", "study_day_tlast", 1, "study day of tlast, 5th percentile subject"),
    t3 = r("min", "study_day_tlast", 1, "study day of tlast, minimum subject"),
    a1 = r("median", "extrap_area_nca", 1, "NCA extrapolated area, median subject"), b1 = r("median", "extrap_area_true", 1, "true extrapolated area, median subject"),
    a2 = r("p05", "extrap_area_nca", 1, "NCA extrapolated area, 5th percentile subject"), b2 = r("p05", "extrap_area_true", 1, "true extrapolated area, 5th percentile subject"),
    a3 = r("min", "extrap_area_nca", 1, "NCA extrapolated area, minimum subject"), b3 = r("min", "extrap_area_true", 1, "true extrapolated area, minimum subject"),
    r1 = r("median", "adj_r2", 2, "adjusted R-squared, median subject"), r2v = r("p05", "adj_r2", 2, "adjusted R-squared, 5th percentile subject"),
    r3 = r("min", "adj_r2", 2, "adjusted R-squared, minimum subject"), e3 = r("min", "pct_extrap_nca", 1, "NCA extrapolated share, minimum subject", "%"),
    r2iii = f_set("iii", "r2"), exiii = f_set("iii", "extrap"),
    cobs = r("min", "Clast", 2, "observed Clast, minimum subject (mg/L)"),
    ctrue = dderived("true concentration at the observed tlast, 2016 minimum subject (mg/L)", ROf, "model=='k2016' & id==<minimum subject> & time_after_dose==tlast :: conc_true", ct, fnum(ct, 2)),
    cr = cr,
    xr = { x <- ct / s3$lambda_z / s3$extrap_area_true
      dderived("extrapolation from the true concentration with the NCA lambda-z / true extrapolated area, 2016 minimum subject", ROf, "conc_true at tlast / lambda_z / extrap_area_true (rep_subjects.csv)", x, fnum(x, 2)) },
    rr16 = dv(ER, "model=='k2016' & group=='reliable_iii'", "median", 1, "", "median NCA-to-true extrapolated area ratio, reliable subjects, 2016"),
    rr20 = dv(ER, "model=='k2020' & group=='reliable_iii'", "median", 1, "", "median NCA-to-true extrapolated area ratio, reliable subjects, 2020"),
    fa16 = dv(ER, "model=='k2016' & group=='failing_iii_estimable'", "median", 1, "", "median NCA-to-true extrapolated area ratio, lambda-z estimable subjects failing set (iii), 2016"),
    fa20 = dv(ER, "model=='k2020' & group=='failing_iii_estimable'", "median", 1, "", "median NCA-to-true extrapolated area ratio, lambda-z estimable subjects failing set (iii), 2020"),
    b1p = dv(ER, "model=='k2016' & group=='reliable_iii'", "pct_below_1", 1, "%", "share with NCA below true, reliable subjects, 2016"),
    b2p = dv(ER, "model=='k2020' & group=='reliable_iii'", "pct_below_1", 2, "%", "share with NCA below true, reliable subjects, 2020"),
    nsub = dint(P1, "model=='k2016' & group=='all'", "n_subjects", "virtual subjects, 2016"),
    wt = f_wt_range(), p95n = dv(P1, "model=='k2016' & group=='all'", "extrap_nca_p95", 1, "%", "95th percentile NCA extrapolated share, 2016"),
    p95t = dv(P1, "model=='k2016' & group=='all'", "extrap_true_p95", 1, "%", "95th percentile true extrapolated share, 2016"))))
  deck_end()
}
