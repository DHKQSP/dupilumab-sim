# S2 결론(지시 §2): 카드 3개와 결론 한 줄. 카드 1 = 세트 (iii) 미달 범위(두 모델), 카드 2 = 기본 조건 창 포착률 최솟값 하한(두 모델)과 중앙값,
# 카드 3 = 경계 1종 오류 최대(M1): AUCinf(세트 (iii) 미달 제외) + Cmax 대 AUClast + Cmax. 표시 규칙: 사전 등록 section7 7b·7e.
slide_S2 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; CG <- "criteria/criteria_g2_type1.csv"; T1 <- "oc_models/type1_models.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  deck_slide("S2", tag = "sim")
  premise(nrow(rows(CG, "analysis_model=='M1' & config=='G2_A_iii'")) == 16, "S9 configuration G2_A_iii under M1 exists (no fallback to set (i))")
  y0 <- core_title(tx("S2.title"), tx("S2.kicker"))

  f <- list(fail = headline(drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models")),
            r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"))
  cmin <- { r <- rows(TCV, MW); premise(nrow(r) == 2, "window coverage: two models"); x <- 100 * min(r$min)
    headline(dderived("window coverage (true), base case, smallest subject over both models, rounded down", TCV, sprintf("%s :: floor(min(min) x 100)", MW), x, paste0(fnum(fl(x, 0), 0), "%"))) }
  cmed <- drange(TCV, MW, "median", 0, "%", "window coverage (true), median, two models", scale = 100)
  t_inf <- headline(dext(CG, "analysis_model=='M1' & config=='G2_A_iii'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii, rule A) + Cmax, M1"))
  t_last <- headline(dext(T1, "analysis_model=='M1' & config=='P2'", "pass_pct", max, 1, "%", "largest boundary type I error, AUClast + Cmax, M1"))
  nom <- f_nominal()

  # ---- 카드 3개(주 시각 요소) ----
  gap <- 0.25; cw <- (GEO$CW - 2 * gap) / 3; ch <- 3.45; cy <- y0 + 0.10
  C <- DK$txt$S2$cards
  core_card(C$c1$head, list(list(f$fail, PAL$orange)), fill(C$c1$label, f), c(GEO$ML, cy, cw, ch), bg = PAL$tint_orange)
  core_card(C$c2$head, list(list(fill(C$c2$value, list(v = cmin)), PAL$blue)), fill(C$c2$label, list(med = cmed)), c(GEO$ML + cw + gap, cy, cw, ch), bg = PAL$tint_blue)
  core_card(fill(C$c3$head, list(nom = nom)), list(list(list(t_inf, PAL$orange, 40), list(C$c3$inf, PAL$ink2, 18)), list(list(t_last, PAL$blue, 40), list(C$c3$last, PAL$ink2, 18))),
            fill(C$c3$label, list(nom = nom, n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"),
                                  z = dcount(T1, "analysis_model=='M1' & config=='P2' & pass_pct <= 5", "AUClast + Cmax cells at or below 5%, M1"))), c(GEO$ML + 2 * (cw + gap), cy, cw, ch), bg = PAL$tint_grey)
  deck_visual(c(GEO$ML, cy, GEO$CW, ch))

  # ---- 결론 한 줄 ----
  by <- cy + ch + 0.35
  deck_text(tx("S2.conclusion"), c(GEO$ML, by, GEO$CW, 0.66), size = 24, bold = TRUE, color = PAL$ink, label = "cmid_body_conclusion", bg = PAL$tint_grey, align = "center")
  DK$cur$body_lines <- DK$cur$body_lines + 1L                         # 결론 한 줄(이름이 cmid_로 시작해 자동으로 세지 않음)

  deck_notes(tx("S2.notes", c(f, list(cmin = cmin, cmed = cmed, t_inf = t_inf, t_last = t_last, nom = nom,
    fail16 = dv(TPF, "pk_model=='k2016' & set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), 2016 model"),
    fail20 = dv(TPF, "pk_model=='k2020' & set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), 2020 model"),
    min16 = dv(TCV, paste(MW, "& pk_model=='k2016'"), "min", 1, "%", "window coverage minimum, 2016 model", scale = 100),
    min20 = dv(TCV, paste(MW, "& pk_model=='k2020'"), "min", 1, "%", "window coverage minimum, 2020 model", scale = 100),
    n_inf = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii' & pass_pct > 5", "M1 G2_A_iii cells above 5%"),
    n_last = dcount(T1, "analysis_model=='M1' & config=='P2' & pass_pct > 5", "M1 P2 cells above 5%"),
    n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"), wt = f_wt_range(),
    nsub = dint(TPF, "pk_model=='k2016' & set=='iii'", "n", "virtual subjects per model (trial population)")))))
  deck_end()
}
