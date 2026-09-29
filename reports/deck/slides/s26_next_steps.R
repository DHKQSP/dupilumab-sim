# S26 남은 결정과 다음 단계(어두운 마무리): 의뢰자 결정 4개(분석법 LLOQ, 잠정 설계값, 주분석 모형(이 덱의 제안 M1)과 SAP 문안, 목표 검정력과 n)
# + 후속 작업 2개(결과 버전 고정 = v1.0.1 태그, BPD 미팅·과학자문 질의 문안) + 제안 요약. [모의]
# 자리표시자 여부는 config의 status 'assumption'으로 검사한다. 태그 상태는 빌드 시점의 로컬 git 태그로 판단한다(D-061: 분석 환경에서 태그 push 거부).
# BPD·과학자문 질의 문안은 프로젝트에 초안이 없다(다음 단계로 표시).

s26_col_thr <- function(rel, col, item) {
  premise(col %in% names(.read(rel)), sprintf("%s has column %s", rel, col))
  x <- as.numeric(sub("^[a-z_]*lt([0-9]+)_pct$", "\\1", col)); premise(is.finite(x), sprintf("threshold in column name %s", col))
  dderived(item, rel, sprintf("column name %s :: threshold", col), x, fnum(x, 0))
}
# 방문 허용 범위 규칙의 한 행(max_day 기준): what = "window_h"(시간), "window_d"(일), "max_day"
s26_window <- function(k, what, item) {
  r <- .read("config/trial_design.yaml")$sampling_window_days$rule; premise(length(r) >= k, "visit-window rule rows")
  x <- switch(what, window_h = 24 * r[[k]]$window, window_d = r[[k]]$window, max_day = r[[k]]$max_day)
  dderived(item, "config/trial_design.yaml", sprintf("sampling_window_days.rule[%d].%s", k, switch(what, window_h = "window x 24", window_d = "window", max_day = "max_day")), x, fnum(x, 0))
}
# 카드 추정 높이(표시층 nobreak 문구와 원문 중 큰 값)
s26_card_h <- function(p, w) max(est_height(p, w, 16, 5, card = TRUE), est_height(vapply(p, nobreak, ""), w, 16, 5, card = TRUE))
s26_tag_exists <- function(tag) length(suppressWarnings(system2("git", c("-C", PROJ_ROOT, "tag", "-l", tag), stdout = TRUE, stderr = FALSE))) > 0

slide_S26 <- function() {
  LIf <- "lloq/lloq_individual_table.csv"; LTf <- "lloq/lloq_trial_type1.csv"; TPf <- "sample_size/ss_table_power.csv"; NNf <- "sample_size/ss_table_n_needed.csv"
  deck_slide("S26", tag = "sim", dark = TRUE)
  deck_kicker(tx("S26.kicker")); deck_title(tx("S26.title"))

  # ---- 전제: 자리표시자 상태, 창 규칙 순서 -----------------------------------------------------------------------------------------
  TD <- .read("config/trial_design.yaml"); AS <- .read("config/assay.yaml")
  premise(AS$lloq_mg_L$status == "assumption", "assay LLOQ is a placeholder (status assumption)")
  premise(TD$day1_postdose_time$status == "assumption" && TD$sampling_window_days$status == "assumption" && TD$stratification$split_kg$status == "assumption" &&
          TD$weight$base$status == "assumption" && TD$weight$sex_ratio_male$status == "assumption",
          "Day 1 post-dose time, visit windows, weight distribution, share male and stratum split are placeholders (status assumption)")
  b0 <- unlist(TD$schedules$B0$days); premise(min(b0[b0 > 0]) == TD$day1_postdose_time$value, "the Day 1 post-dose time is the first post-dose sample of B0 (card: first sample)")
  # LLOQ 시험(경계 1종 오류)은 2016 모델의 경계 3칸, 칸당 5,000회뿐이다(파일에 PK 모델 열이 없고, 배율이 2016 모델 경계 칸과 같다)
  lt <- rows(LTf); t1k <- rows("oc_models/type1_models.csv", "pk_model=='k2016' & analysis_model=='M1' & config=='P2'")
  premise(!"pk_model" %in% names(lt) && all(lt$scenario %in% t1k$scenario) &&
          all(abs(lt$multiplier - t1k$multiplier[match(lt$scenario, t1k$scenario)]) < 1e-9), "LLOQ trials are 2016-model boundary cells (same multipliers as the 2016 model cells)")
  wr <- TD$sampling_window_days$rule; premise(length(wr) == 4 && wr[[1]]$window == 0 && wr[[4]]$max_day > 100, "visit-window rule: pre-dose, two early windows, then one window thereafter")
  cv_v <- as.numeric(.read("config/prereg_20260926.yaml")$section3$base_cv_pct); n_v <- as.numeric(TD$n_per_arm)
  wp <- function(am) sprintf("input_model=='k2016' & cv==%s & gmr==0.95 & n==%s & analysis_model=='%s'", cv_v, n_v, am)
  wn <- function(am, tg) sprintf("cv==%s & gmr==0.95 & analysis_model=='%s' & target_pct==%s", cv_v, am, tg)
  tag <- DK$version; has_tag <- s26_tag_exists(tag)
  # BPD·과학자문 질의 초안이 규제 문서 원본(regulatory/src)에 없다(예상 심사 질의응답, SAP 문안, 보고서만 있음)
  src_txt <- unlist(lapply(list.files(proj_path("regulatory", "src"), pattern = "\\.Rmd$", full.names = TRUE), readLines, warn = FALSE, encoding = "UTF-8"))
  premise(length(src_txt) > 0 && !any(grepl("BPD|Biosimilar Product Development|scientific advice|Type 2 meeting|\uacfc\ud559\uc790\ubb38", src_txt, ignore.case = TRUE)),
          "no BPD meeting or scientific-advice question draft in regulatory/src")

  f <- list(
    lloq = f_lloq(), grid = f_lloq_grid(),
    thr = s26_col_thr(LIf, "coverage_lt80_pct", "window coverage threshold (%) of the column coverage_lt80_pct"),
    cov = drange(LIf, "model=='k2016' & resid=='fixed'", "coverage_lt80_pct", 2, "%", "LLOQ grid, window coverage below 80%, 2016 model"),
    p2 = drange(LTf, "model=='M1' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid M1 P2 range"),
    t1 = dcfg("trial_design.yaml", c("day1_postdose_time", "value"), "Day 1 post-dose sampling time (day after dose)", function(x) format(x)),
    h0 = dcfg("trial_design.yaml", c("day1_postdose_time", "value"), "Day 1 post-dose sampling time (hours after dose = value x 24)", function(x) fnum(24 * x, 0)),
    wm = dcfg("trial_design.yaml", c("weight", "base", "mean"), "body weight mean (kg)", num_fmt(0)), wsd = dcfg("trial_design.yaml", c("weight", "base", "sd"), "body weight SD (kg)", num_fmt(0)),
    nc = dderived("boundary cells in the LLOQ trials", LTf, "length(unique(scenario))", length(unique(lt$scenario)), as.character(length(unique(lt$scenario)))),
    nt = dint(LTf, "model=='M1' & resid=='fixed' & config=='P2' & scenario=='F_down_080' & abs(lloq - 0.078) < 1e-9", "n_trials", "LLOQ trials per cell"),
    n_rand = f_n_rand(),
    h1 = s26_window(2, "window_h", "visit window, first post-dose interval (hours)"), d3 = s26_window(4, "window_d", "visit window thereafter (days)"),
    split = f_split(), male = dcfg("trial_design.yaml", c("weight", "sex_ratio_male", "value"), "share male (%)", function(x) fnum(100 * as.numeric(x), 0)),
    n_arm = f_n_arm(), cv = dcfg("prereg_20260926.yaml", c("section3", "base_cv_pct"), "protocol CV (%)", num_fmt(0)),
    g = dv(TPf, wp("M1"), "gmr", 2, "", "true GMR of the power row"),
    pw1 = dv(TPf, wp("M1"), "analytic_pct", 1, "%", "P2 power at the protocol n, M1"),
    tg90 = dint(NNf, wn("M1", 90), "target_pct", "target power (%)"), n90 = dint(NNf, wn("M1", 90), "n_evaluable_per_arm", "n per arm for the higher target, M1"),
    tg85 = dint(NNf, wn("M1", 85), "target_pct", "target power (%)"), n85 = dint(NNf, wn("M1", 85), "n_evaluable_per_arm", "n per arm for the lower target, M1"),
    tag = tag)

  # ---- 1행: 의뢰자 결정 4개 ------------------------------------------------------------------------------------------------------
  CARD <- "#223246"; INK <- PAL$dark_ink; LAB <- "#8fb8ea"
  # 두 행 모두 카드 간격 gap. 카드 폭은 글 양과 렌더링된 줄바꿈을 보고 나눈다(1행: LLOQ 카드가 가장 길고, 설계값·검정력 카드는 한 줄 항목이
  # 줄바꿈되지 않는 폭; 2행: 버전 고정 카드가 가장 길고, 제안 카드는 2차 평가변수 줄이 한 줄에 들어가는 폭).
  # 2행 높이 = 2행 카드의 추정 높이 중 최댓값 + 아래 여백, 1행이 나머지
  y1 <- GEO$BODY_TOP; gap <- 0.2; lh <- 0.42
  keys <- c("lloq", "design", "model", "power"); row1 <- lapply(keys, function(k) tx(sprintf("S26.cards.%s", k), f))
  w4 <- c(3.53, 2.85, 2.45, 2.80); w4 <- w4 / sum(w4) * (GEO$CW - 3 * gap); x4 <- GEO$ML + c(0, cumsum(w4 + gap))[1:4]
  row2 <- list(tx(if (has_tag) "S26.cards.tag_done" else "S26.cards.tag_none", f), tx("S26.cards.bpd"), tx("S26.cards.proposal"))
  w3 <- c(4.32, 3.71, 3.80); w3 <- w3 / sum(w3) * (GEO$CW - 2 * gap); x3 <- GEO$ML + c(0, cumsum(w3 + gap))[1:3]
  h2 <- max(mapply(s26_card_h, row2, w3)) + 0.1; h1 <- GEO$BODY_BOTTOM - y1 - lh - gap - h2
  deck_text(tx("S26.row1"), c(GEO$ML, y1 - 0.06, GEO$CW, lh), size = 16, bold = TRUE, color = LAB, label = "label_row1")
  for (k in seq_along(keys))
    deck_text(row1[[k]], c(x4[k], y1 + lh, w4[k], h1), size = 16, color = INK, bg = CARD, geom = "roundRect",
              label = sprintf("card_%s", keys[k]), gap_pt = 5)
  # ---- 2행: 후속 작업 2개 + 제안 요약 ------------------------------------------------------------------------------------------------
  y2 <- y1 + lh + h1 + gap
  deck_text(row2[[1]], c(x3[1], y2, w3[1], h2), size = 16, color = INK, bg = CARD, geom = "roundRect", label = "card_tag", gap_pt = 5)
  deck_text(row2[[2]], c(x3[2], y2, w3[2], h2), size = 16, color = INK, bg = CARD, geom = "roundRect", label = "card_bpd", gap_pt = 5)
  deck_text(row2[[3]], c(x3[3], y2, w3[3], h2), size = 16, color = "#ffffff", bg = PAL$blue, geom = "roundRect", label = "card_proposal", gap_pt = 5)

  deck_notes(tx(if (has_tag) "S26.notes_tag_done" else "S26.notes", c(f, list(
    lloq_st = dcfg("assay.yaml", c("lloq_mg_L", "status"), "assay LLOQ status", function(x) as.character(x)),
    r02 = dv(LIf, "model=='k2016' & resid=='fixed' & abs(lloq - 0.02) < 1e-9", "reliable_no_span_pct", 1, "%", "LLOQ 0.02 fixed k2016 reliable_no_span_pct"),
    r078 = dv(LIf, "model=='k2016' & resid=='fixed' & abs(lloq - 0.078) < 1e-9", "reliable_no_span_pct", 1, "%", "LLOQ 0.078 fixed k2016 reliable_no_span_pct"),
    r05 = dv(LIf, "model=='k2016' & resid=='fixed' & abs(lloq - 0.5) < 1e-9", "reliable_no_span_pct", 1, "%", "LLOQ 0.5 fixed k2016 reliable_no_span_pct"),
    p2m0 = drange(LTf, "model=='M0' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid M0 P2 range"),
    a1 = s26_window(2, "max_day", "visit window, end of first interval (day after dose)"), a2 = s26_window(3, "max_day", "visit window, end of second interval (day after dose)"),
    h2 = s26_window(3, "window_h", "visit window, second interval (hours)"),
    wt = f_wt_range(),
    pw0 = dv(TPf, wp("M0"), "analytic_pct", 1, "%", "P2 power at the protocol n, M0"),
    n90_0 = dint(NNf, wn("M0", 90), "n_evaluable_per_arm", "n per arm for the higher target, M0"), n85_0 = dint(NNf, wn("M0", 85), "n_evaluable_per_arm", "n per arm for the lower target, M0"),
    cv50 = dcfg("prereg_20260926.yaml", c("section3", "sensitivity_cv_pct"), "sensitivity CV (%)", num_fmt(0)),
    n90_50 = dint(NNf, sprintf("cv==%s & gmr==0.95 & analysis_model=='M1' & target_pct==90", .read("config/prereg_20260926.yaml")$section3$sensitivity_cv_pct), "n_evaluable_per_arm", "n per arm, higher target, sensitivity CV, M1")))))
  deck_end()
}
