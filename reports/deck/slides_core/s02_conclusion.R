# S2 결론(지시 §2): 카드 3개와 결론 한 줄. 카드 1 = 세트 (iii) 미달 범위(두 모델), 카드 2 = 본문 S6의 주 모델(2020) 대표 대상자 AUClast/AUCinf,
# 카드 3(v1.2, 지시 2026-09-29 "S9 재설계" §5) = 본문 S7의 동일 제품 2종 오류 AUCinf(A) + Cmax 대 AUClast + Cmax(두 모델 범위), 부제에 경계 1종 오류 최대.
#   사전 등록 규칙 4a(config/prereg_20260929_oc.yaml section8 wording)가 성립하지 않으면 큰 숫자를 1종 오류 최대로 바꾼다. 표시 규칙: 사전 등록 section7 7b·7e.
slide_S2 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; CG <- "criteria/criteria_g2_type1.csv"; T1 <- "oc_models/type1_models.csv"
  PF <- "oc_curves/oc_curves_pass.csv"; OT1 <- "oc_curves/oc_type1_summary.csv"; WR <- "oc_curves/oc_wording_rules.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  deck_slide("S2", tag = "sim")
  premise(nrow(rows(CG, "analysis_model=='M1' & config=='G2_A_iii'")) == 16, "S9 configuration G2_A_iii under M1 exists (no fallback to set (i))")
  premise(all(rows(TPF, "set=='iii'")$fail_pct > 30) && all(rows(TPF, "set=='iii'")$fail_pct < 40), "set (iii) failing share between 30% and 40% in both models (title: about one in three)")
  y0 <- core_title(tx("S2.title"), tx("S2.kicker"))

  f <- list(fail = headline(drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models")),
            r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"))
  RSf <- "core_deck/rep_subjects.csv"                                # 카드 2 = 본문 S6의 수치(주 모델 2020 대표 대상자; v1.2 지시: 본문 수치만)
  cmin <- dv(RSf, "model=='k2020' & role=='min'", "coverage_true", 1, "%", "window coverage (true), 2020 minimum-coverage subject", scale = 100)
  cmed <- dv(RSf, "model=='k2020' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2020 median-coverage subject", scale = 100)
  t_inf <- dext(CG, "analysis_model=='M1' & config=='G2_A_iii'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii, rule A) + Cmax, M1")
  t_last <- dext(T1, "analysis_model=='M1' & config=='P2'", "pass_pct", max, 1, "%", "largest boundary type I error, AUClast + Cmax, M1")
  nom <- f_nominal()
  r4a <- all(rows(WR, "TRUE")$rule_4a_keep_type2)                     # 사전 등록 규칙 4a: 두 모델 모두 G2 동일 제품 2종 오류가 P2보다 2.0%p 이상 크면 2종 오류 카드
  t2id <- function(cf) { r <- rows(PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s'", cf)); x <- range(100 - r$pass_pct)
    dderived(sprintf("type II error, identical product, %s, M1, range over two models", cf), PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s' :: range(100 - pass_pct)", cf), x, rng_fmt(x[1], x[2], 1, "%")) }
  t2id20 <- function(cf) { w <- sprintf("code=='S00' & analysis_model=='M1' & config=='%s' & pk_model=='k2020'", cf); x <- 100 - row1(PF, w)$pass_pct
    dderived(sprintf("type II error, identical product, %s, M1, 2020 model", cf), PF, paste(w, ":: 100 - pass_pct"), x, paste0(fnum(x, 1), "%")) }
  g1 <- dv(OT1, "config=='G2A_iii' & scope=='k2020'", "max_pct", 1, "%", "largest boundary type I error, G2A_iii, 2020 model, M1")      # 부제도 큰 숫자와 같은 주 모델(2020) 값
  p1 <- dv(OT1, "config=='P2' & scope=='k2020'", "max_pct", 1, "%", "largest boundary type I error, P2, 2020 model, M1")
  g1r <- drange(OT1, "config=='G2A_iii' & scope!='both'", "max_pct", 1, "%", "largest boundary type I error per model, AUCinf (set iii) + Cmax, M1, range over two models")

  # ---- 카드 3개(주 시각 요소) ----
  gap <- 0.25; cw <- (GEO$CW - 2 * gap) / 3; ch <- 3.45; cy <- y0 + 0.10
  C <- DK$txt$S2$cards
  core_card(C$c1$head, list(list(f$fail, PAL$orange)), fill(C$c1$label, f), c(GEO$ML, cy, cw, ch), bg = PAL$tint_orange)
  nsub <- dint(TPF, "pk_model=='k2020' & set=='iii'", "n", "virtual subjects per model (trial population)")
  core_card(C$c2$head, list(list(fill(C$c2$value, list(v = cmin)), PAL$blue)), fill(C$c2$label, list(med = cmed, nsub = nsub)), c(GEO$ML + cw + gap, cy, cw, ch), bg = PAL$tint_blue)
  K3 <- if (r4a) C$c3$type2 else C$c3$type1
  if (r4a) {                                                           # 큰 숫자 = 주 모델(2020) 동일 제품 2종 오류(본문 S7 표의 2020 모델 값), 두 모델 범위는 노트
    v_inf <- headline(t2id20("G2A_iii")); v_last <- headline(t2id20("P2"))
    core_card(K3$head, list(list(list(v_inf, PAL$orange, 40), list(C$c3$inf, PAL$ink2, 18)), list(list(v_last, PAL$blue, 40), list(C$c3$last, PAL$ink2, 18))),
              fill(K3$label, list(g1 = g1, p1 = p1)), c(GEO$ML + 2 * (cw + gap), cy, cw, ch), bg = PAL$tint_grey)
  } else {
    core_card(K3$head, list(list(list(headline(t_inf), PAL$orange, 40), list(C$c3$inf, PAL$ink2, 18)), list(list(headline(t_last), PAL$blue, 40), list(C$c3$last, PAL$ink2, 18))),
              fill(K3$label, list(nom = nom, n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"))), c(GEO$ML + 2 * (cw + gap), cy, cw, ch), bg = PAL$tint_grey)
  }
  deck_visual(c(GEO$ML, cy, GEO$CW, ch))

  # ---- 결론 한 줄 ----
  by <- cy + ch + 0.35
  deck_text(tx("S2.conclusion"), c(GEO$ML, by, GEO$CW, 0.66), size = 24, bold = TRUE, color = PAL$ink, label = "cmid_body_conclusion", bg = PAL$tint_grey, align = "center")
  DK$cur$body_lines <- DK$cur$body_lines + 1L                         # 결론 한 줄(이름이 cmid_로 시작해 자동으로 세지 않음)

  deck_notes(tx("S2.notes", c(f, list(cmin = cmin, cmed = cmed, t_inf = t_inf, t_last = t_last, nom = nom, g1 = g1, p1 = p1,
    c3 = DK$txt$S2$c3note[[if (r4a) "type2" else "type1"]],
    g2id = t2id("G2A_iii"), p2id = t2id("P2"),
    g2id20 = t2id20("G2A_iii"), p2id20 = t2id20("P2"),
    g1r = g1r,
    p1_16 = dv(OT1, "config=='P2' & scope=='k2016'", "max_pct", 1, "%", "largest boundary type I error, P2, 2016 model, M1"),
    fail16 = dv(TPF, "pk_model=='k2016' & set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), 2016 model"),
    fail20 = dv(TPF, "pk_model=='k2020' & set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), 2020 model"),
    min16 = dv(TCV, paste(MW, "& pk_model=='k2016'"), "min", 1, "%", "window coverage minimum, 2016 model", scale = 100),
    min20 = dv(TCV, paste(MW, "& pk_model=='k2020'"), "min", 1, "%", "window coverage minimum, 2020 model", scale = 100),
    n_inf = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii' & pass_pct > 5", "M1 G2_A_iii cells above 5%"),
    n_last = dcount(T1, "analysis_model=='M1' & config=='P2' & pass_pct > 5", "M1 P2 cells above 5%"),
    n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"), wt = f_wt_range(),
    nsub = nsub, redraw = dv("core_deck/coverage_by_case.csv", "case=='curve_base'", "min", 1, "%", "window coverage minimum, curve-shape study base draw", scale = 100)))))
  deck_end()
}
