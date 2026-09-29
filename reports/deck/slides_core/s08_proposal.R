# S10 결론·제안(지시 §2): 상자 4개(주 시각 요소) = 1차 AUClast + Cmax; 2차 AUCinf 계산·보고(λz 산출 전원과 기준 충족자, 두 분석군);
# 채혈 현행 유지(후보 일정 D1~D4, B- 모두 사전 기준 미충족: schedule_decision_*.csv, config/schedule_decision.yaml); AUCinf 공동 1차 요구 시 AUClast 유지하고 추가(별첨 A10).
# 본문 한 줄: 세 근거(본문 수치만, v1.2). ③은 사전 등록 규칙 T(config/prereg_20260929_oc.yaml section8)가 성립하면 "AUClast + Cmax만 두 오류를 모두 통제"(S7 제목과 같은 주장),
# 성립하지 않으면 경계 1종 오류 최대(S7 표)로 적는다.
slide_S8 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"
  deck_slide("S8", tag = "sim")
  for (v in c("base", "struct2020")) { r <- rows(sprintf("trials/schedule_decision_%s.csv", v))
    premise(!any(r$crit_a %in% TRUE | r$crit_b %in% TRUE | r$crit_c %in% TRUE | r$crit_d %in% TRUE | r$recommend %in% TRUE), paste("no candidate schedule meets a criterion,", v)) }
  premise(.read("config/schedule_decision.yaml")$final_schedule == "B0", "final schedule decision is B0")
  f <- list(nB0 = f_study_days("B0", "n"), last = f_study_days("B0", "last"))
  y0 <- core_title(tx("S8.title", f), tx("S8.kicker"))
  B <- DK$txt$S8$boxes
  OT1 <- "oc_curves/oc_type1_summary.csv"; ruleT <- all(rows("oc_curves/oc_wording_rules.csv", "TRUE")$rule_T)
  e3 <- if (ruleT) DK$txt$S8$ev3$ruleT else fill(DK$txt$S8$ev3$type1, list(
    t_inf = dv(OT1, "config=='G2A_iii' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, G2A_iii, M1"),
    t_last = dv(OT1, "config=='P2' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, P2, M1")))
  gap <- 0.25; bw <- (GEO$CW - gap) / 2; body <- tx("S8.body", list(
    fail = drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models"),
    cmin = dv("core_deck/rep_subjects.csv", "model=='k2020' & role=='min'", "coverage_true", 1, "%", "window coverage (true), 2020 minimum-coverage subject", scale = 100),
    e3 = e3))
  by <- GEO$BODY_BOTTOM - core_body_h(body)
  bh <- (by - 0.20 - y0 - gap) / 2
  fills <- c(PAL$tint_blue, PAL$tint_orange, PAL$tint_grey, PAL$tint_grey); hc <- c(PAL$blue, PAL$orange, PAL$ink, PAL$ink)
  for (k in 1:4) {
    x <- GEO$ML + ((k - 1) %% 2) * (bw + gap); y <- y0 + ((k - 1) %/% 2) * (bh + gap); b <- B[[k]]
    fit_check("proposal", c(b$head, b$main, fill(b$sub, f)), c(x, y, bw, bh), SZ$body, gap_pt = 6, card = TRUE)
    subs <- strsplit(fill(b$sub, f), "\n", fixed = TRUE)[[1]]
    ps <- c(list(para(b$head, SZ$body, PAL$ink2, TRUE, "left", gap_pt = 4), para(b$main, 26, hc[k], TRUE, "left", gap_pt = 4, line = 1.0)),
            lapply(subs, function(z) para(z, SZ$body, PAL$ink, FALSE, "left", gap_pt = 2)))
    DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(c(x, y, bw, bh), "proposal", bg = fills[k], geom = "roundRect", ln = no_line()))   # 위 정렬(상자끼리 줄 맞춤)
  }
  deck_visual(c(GEO$ML, y0, GEO$CW, 2 * bh + gap))
  core_body(body, GEO$BODY_BOTTOM)
  PF <- "oc_curves/oc_curves_pass.csv"
  t2i <- function(cf) { r <- rows(PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s'", cf)); x <- range(100 - r$pass_pct)
    dderived(sprintf("type II error, identical product, %s, M1, range over two models", cf), PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s' :: range(100 - pass_pct)", cf), x, rng_fmt(x[1], x[2], 1, "%")) }
  deck_notes(tx("S8.notes", c(f, list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"),
    f3t1 = dv(OT1, "config=='F3A_iii' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, F3A_iii, M1"),
    p2t1 = dv(OT1, "config=='P2' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, P2, M1"),
    g2t1 = dv(OT1, "config=='G2A_iii' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, G2A_iii, M1"),
    f3id = t2i("F3A_iii"), p2id = t2i("P2"), g2id = t2i("G2A_iii")))))
  deck_end()
}
