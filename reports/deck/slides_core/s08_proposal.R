# S10 결론·제안(지시 §2): 상자 4개(주 시각 요소) = 1차 AUClast + Cmax; 2차 AUCinf 계산·보고(λz 산출 전원과 기준 충족자, 두 분석군);
# 채혈 현행 유지(후보 일정 D1~D4, B- 모두 사전 기준 미충족: schedule_decision_*.csv, config/schedule_decision.yaml); AUCinf 공동 1차 요구 시 AUClast 유지하고 추가(별첨 A10).
# 본문 한 줄: 세 근거의 대표 수치(S2와 같은 파일·조건).
slide_S8 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; CG <- "criteria/criteria_g2_type1.csv"; T1 <- "oc_models/type1_models.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  deck_slide("S8", tag = "sim")
  for (v in c("base", "struct2020")) { r <- rows(sprintf("trials/schedule_decision_%s.csv", v))
    premise(!any(r$crit_a %in% TRUE | r$crit_b %in% TRUE | r$crit_c %in% TRUE | r$crit_d %in% TRUE | r$recommend %in% TRUE), paste("no candidate schedule meets a criterion,", v)) }
  premise(.read("config/schedule_decision.yaml")$final_schedule == "B0", "final schedule decision is B0")
  f <- list(nB0 = f_study_days("B0", "n"), last = f_study_days("B0", "last"))
  y0 <- core_title(tx("S8.title", f), tx("S8.kicker"))
  B <- DK$txt$S8$boxes
  gap <- 0.25; bw <- (GEO$CW - gap) / 2; body <- tx("S8.body", list(
    fail = drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models"),
    cmed = dv("core_deck/rep_subjects.csv", "model=='k2020' & role=='median'", "coverage_true", 0, "%", "window coverage (true), 2020 median-coverage subject", scale = 100),
    cmin = dv("core_deck/rep_subjects.csv", "model=='k2020' & role=='min'", "coverage_true", 1, "%", "window coverage (true), 2020 minimum-coverage subject", scale = 100),
    t_inf = dext(CG, "analysis_model=='M1' & config=='G2_A_iii'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii) + Cmax, M1"),
    t_last = dext(T1, "analysis_model=='M1' & config=='P2'", "pass_pct", max, 1, "%", "largest boundary type I error, AUClast + Cmax, M1")))
  by <- GEO$BODY_BOTTOM - core_body_h(body)
  bh <- (by - 0.20 - y0 - gap) / 2
  fills <- c(PAL$tint_blue, PAL$tint_orange, PAL$tint_grey, PAL$tint_grey); hc <- c(PAL$blue, PAL$orange, PAL$ink, PAL$ink)
  for (k in 1:4) {
    x <- GEO$ML + ((k - 1) %% 2) * (bw + gap); y <- y0 + ((k - 1) %/% 2) * (bh + gap); b <- B[[k]]
    fit_check("proposal", c(b$head, b$main, fill(b$sub, f)), c(x, y, bw, bh), SZ$body, gap_pt = 6, card = TRUE)
    subs <- strsplit(fill(b$sub, f), "\n", fixed = TRUE)[[1]]
    ps <- c(list(para(b$head, SZ$body, PAL$ink2, TRUE, "left", gap_pt = 4), para(b$main, 28, hc[k], TRUE, "left", gap_pt = 6, line = 1.0)),
            lapply(subs, function(z) para(z, SZ$body, PAL$ink, FALSE, "left", gap_pt = 2)))
    DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(c(x, y, bw, bh), "proposal", bg = fills[k], geom = "roundRect", ln = no_line()))   # 위 정렬(상자끼리 줄 맞춤)
  }
  deck_visual(c(GEO$ML, y0, GEO$CW, 2 * bh + gap))
  core_body(body, GEO$BODY_BOTTOM)
  deck_notes(tx("S8.notes", c(f, list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap")))))
  deck_end()
}
