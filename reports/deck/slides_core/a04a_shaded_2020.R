# 별첨 A4a ② 세부(S7의 2020 모델판): 2020 모델 대표 대상자 3명(AUClast/AUCinf(참값) 중앙값·5백분위·최솟값; 사전 등록 7c)의 음영 그림.
# 그림은 S7과 같은 core_shaded_panels("k2020", L). 최솟값 대상자는 λz를 산출하지 못해 비구획 외삽(AUCinf)이 없다(삽입도 제목 L$no_lz).
# 제목: 중앙값 대상자의 AUClast/AUCinf(반올림)와 최솟값 대상자의 하한(내림). 캡션: 모집단 외삽 비율 중앙값(비구획 대 참, pillar1_coverage_B0.csv).
slide_A4a <- function() {
  RSf <- "core_deck/rep_subjects.csv"; ROf <- "core_deck/rep_obs.csv"; P1 <- "rationale/pillar1_coverage_B0.csv"; NR <- "config/nca_rules.yaml"
  deck_slide("A4a", tag = "sim")
  S <- rows(RSf, "model=='k2020'")
  premise(nrow(S) == 3 && setequal(S$role, c("median", "p05", "min")), "three representative subjects of the 2020 model")
  premise(!isTRUE(S[role == "min", lambda_ok]) && all(S[role != "min", lambda_ok]), "2020 model: only the minimum-coverage subject has no estimable lambda-z (text)")
  premise(all(S[role != "min", reliable_iii]) && all(S[role != "min", extrap_ratio_nca_to_true] > 1),
          "2020 model: median and 5th percentile subjects meet set (iii) and their NCA extrapolated area exceeds the true one (text)")
  premise(row1(P1, "model=='k2020' & group=='all'")[, extrap_nca_median > extrap_true_median], "2020 model: NCA extrapolated share above the true share at the median (caption: over-estimates)")
  med <- dv(RSf, "model=='k2020' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2020 median-coverage subject", scale = 100)
  mn <- dv(RSf, "model=='k2020' & role=='min'", "coverage_true", 1, "%", "window coverage (true), minimum subject, one decimal", scale = 100)   # S7과 같이 그림의 최솟값 대상자 값 그대로
  y0 <- core_title(tx("A4a.title", list(med = med, min = mn)), tx("A4a.kicker"))
  L <- DK$txt$A4a$fig
  p <- core_shaded_panels("k2020", L)
  r <- function(role, col, d, item, unit = "") dv(RSf, sprintf("model=='k2020' & role=='%s'", role), col, d, unit, item)
  premise(as.numeric(sub(" .*", "", row1(P1, "model=='k2020' & group=='all'")$lambda_ok_pct_ci)) < 100,
          "not every subject has lambda-z: the NCA median is over lambda-z-estimable subjects, the true median over all subjects (caption denominators; scripts/15 median(pct_extrap, na.rm = TRUE))")
  cap <- tx("A4a.caption", list(nca = dv(P1, "model=='k2020' & group=='all'", "extrap_nca_median", 1, "%", "median NCA extrapolated share, 2020"),
                                true = dv(P1, "model=='k2020' & group=='all'", "extrap_true_median", 2, "%", "median true extrapolated share, 2020"),
                                nsub = dint(P1, "model=='k2020' & group=='all'", "n_subjects", "virtual subjects, 2020")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  body <- tx("A4a.body", list(q1 = r("median", "extrap_ratio_nca_to_true", 1, "NCA-to-true extrapolated area ratio, median subject"),
                              q2 = r("p05", "extrap_ratio_nca_to_true", 1, "NCA-to-true extrapolated area ratio, 5th percentile subject")))
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "a4a_shaded_k2020", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(RSf, "core_deck/rep_profiles.csv", ROf))

  # 노트: 최솟값 대상자가 λz를 얻지 못한 이유(최고 농도 뒤 정량 시료 수 < λz 최소 점 수; nca_rules.yaml)
  s3 <- row1(RSf, "model=='k2020' & role=='min'"); ob <- rows(ROf, sprintf("model=='k2020' & id==%d", s3$id))[blq == FALSE][order(time_after_dose)]
  npost <- nrow(ob) - which.max(ob$conc_obs); mp <- as.integer(.read(NR)$standard$lambda_z$min_points)
  premise(npost < mp && isTRUE(.read(NR)$standard$lambda_z$after_tmax_only), "2020 minimum subject: fewer quantified samples after the highest concentration than the lambda-z minimum (notes)")
  deck_notes(tx("A4a.notes", list(
    i1 = r("median", "id", 0, "subject id, median"), i2 = r("p05", "id", 0, "subject id, 5th percentile"), i3 = r("min", "id", 0, "subject id, minimum"),
    t1 = r("median", "study_day_tlast", 1, "study day of the observed tlast, median subject"), t2 = r("p05", "study_day_tlast", 1, "study day of tlast, 5th percentile subject"),
    t3 = r("min", "study_day_tlast", 1, "study day of tlast, minimum subject"),
    a1 = r("median", "extrap_area_nca", 1, "NCA extrapolated area, median subject"), b1 = r("median", "extrap_area_true", 1, "true extrapolated area, median subject"),
    a2 = r("p05", "extrap_area_nca", 1, "NCA extrapolated area, 5th percentile subject"), b2 = r("p05", "extrap_area_true", 1, "true extrapolated area, 5th percentile subject"),
    r1 = r("median", "adj_r2", 3, "adjusted R-squared, median subject"), r2 = r("p05", "adj_r2", 3, "adjusted R-squared, 5th percentile subject"),
    e1 = r("median", "pct_extrap_nca", 1, "NCA extrapolated share, median subject", "%"), e2 = r("p05", "pct_extrap_nca", 1, "NCA extrapolated share, 5th percentile subject", "%"),
    r2iii = f_set("iii", "r2"), exiii = f_set("iii", "extrap"),
    b3 = r("min", "extrap_area_true", 1, "true extrapolated area, minimum subject"),
    c3 = dv(RSf, "model=='k2020' & role=='min'", "coverage_true", 1, "%", "window coverage (true), minimum subject, one decimal", scale = 100),
    np = dderived("quantified samples after the highest observed concentration, 2020 minimum subject", ROf,
                  sprintf("model=='k2020' & id==%d & blq==FALSE :: count after which.max(conc_obs)", s3$id), npost, fnum(npost, 0)),
    mp = dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum number of points", num_fmt(0)),
    wt = f_wt_range(), p95n = dv(P1, "model=='k2020' & group=='all'", "extrap_nca_p95", 1, "%", "95th percentile NCA extrapolated share, 2020"),
    p95t = dv(P1, "model=='k2020' & group=='all'", "extrap_true_p95", 1, "%", "95th percentile true extrapolated share, 2020"))))
  deck_end()
}
