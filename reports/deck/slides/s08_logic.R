# S08 논리 구조 도식: 전제(말단 절벽) → 논거 ①, ②, ③(병렬, 세 줄) → 결론. 왼쪽 전제 카드, 가운데 논거 세 줄, 오른쪽 결론 카드, 화살표는 전제에서 각 논거로, 각 논거에서 결론으로.
# 단계마다 대표 수치 하나(headline).
# 시험 모집단(건강인, 체중 층화, B0)만. 주분석 M1, M0은 노트에 병기.
# 글자 없는 도형: deck_box는 8pt 자리 글자를 넣어 검사 6(글자 크기)에 걸리므로 16pt 빈 글상자로 그린다
s08_panel <- function(box, fill, geom = "roundRect", label = "panel") deck_text(" ", box, size = 16, bg = fill, geom = geom, label = label)
# f_study_days()와 같은 계산(config 채혈일이 정수·실수 혼합이라 yaml이 목록으로 읽으므로 unlist; 공용 함수 수정 요청)
s08_study_days <- function(schedule = "B0", which = c("all", "last", "n")) {
  which <- match.arg(which); d <- unlist(.read("config/trial_design.yaml")$schedules[[schedule]]$days); premise(length(d) > 0, paste("schedule", schedule))
  sd <- d + 1; p <- switch(which, all = paste(fnum(sd[sd == round(sd)], 0), collapse = ", "), last = fnum(max(sd), 0), n = as.character(length(d)))
  dderived(sprintf("schedule %s, %s (study day = days after dose + 1)", schedule, which), "config/trial_design.yaml", sprintf("schedules.%s.days :: %s", schedule, which), sd, p)
}
s08_inst <- function(am, what = c("median", "max")) {
  what <- match.arg(what); CSf <- "criteria/criteria_instability.csv"
  r <- rows(CSf, sprintf("analysis_model=='%s' & !scenario %%in%% c('S00','F097')", am)); premise(nrow(r) == 16, "16 boundary cells in the instability file")
  x <- if (what == "median") median(r$inst_all) else max(r$inst_all)
  dderived(sprintf("decision instability inst_all, %s over 16 boundary cells, %s", what, am), CSf,
           sprintf("analysis_model=='%s' & !scenario in (S00, F097) :: %s of inst_all", am, what), x, paste0(fnum(x, 1), "%"))
}

# 창 포착률 최솟값: "이상"이라고 쓰므로 0.1 단위로 내림(보고서 tpcv_min과 같은 파일·조건·계산)
s08_cov_min <- function() {
  TCV <- "trialpop/tp_coverage_individual.csv"; w <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"
  r <- rows(TCV, w); premise(nrow(r) == 2, "window coverage: two models"); x <- min(r$min) * 100
  dderived("window coverage (true AUC0-tlast / true AUC0-inf), smallest subject-level value over both models", TCV, sprintf("%s :: min(min)", w), x, sprintf("%s%%", fnum(floor(x * 10) / 10, 1)))
}

# 현행 B0에서 대상자당 절벽 안 채혈점 수의 최댓값(두 모델, 명목일과 허용창 반영, 1일 정의): pct_ge{k} > 0인 가장 큰 k
s08_cliff_max <- function() {
  CP <- "cliff/cliff_points.csv"; w <- "model %in% c('k2016','k2020') & weight=='base' & definition_day==1 & schedule=='current'"
  r <- rows(CP, w); premise(nrow(r) == 4 && setequal(r$timing, c("nominal", "windowed")), "cliff points: two models x nominal and windowed timing (B0, trial weight range)")
  k <- max(c(0L, which(c(any(r$pct_ge1 > 0), any(r$pct_ge2 > 0), any(r$pct_ge3 > 0)))))
  premise(k == 1L && all(r$pct_ge2 == 0), "B0: no subject has two or more samples on the cliff (premise card: at most one)")
  dderived("largest number of B0 samples on the cliff per subject, both models, nominal and windowed timing (1-day definition)", CP,
           sprintf("%s :: largest k with pct_ge{k} > 0 in any row (pct_ge1 > 0, pct_ge2 == 0, pct_ge3 == 0)", w), k, fnum(k, 0))
}

slide_S08 <- function() {
  CS <- "cliff/cliff_summary.csv"; TPF <- "trialpop/tp_failure_by_set.csv"; TCV <- "trialpop/tp_coverage_individual.csv"; TCH <- "trialpop/tp_characteristics.csv"
  CGf <- "criteria/criteria_g2_type1.csv"; T1f <- "oc_models/type1_models.csv"
  MW <- "metric=='window coverage (true AUC0-tlast / true AUC0-inf)'"; CB <- "model %in% c('k2016','k2020') & weight=='base'"
  deck_slide("S08", tag = "litsim")
  deck_kicker(tx("S08.kicker")); deck_title(tx("S08.title"))

  # 전제 검사
  premise(nrow(rows(CS, CB)) == 2, "cliff summary: two models, trial weight range only (base rows)")
  premise(all(rows(TCH, "set=='i'")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUC0-inf (text: not random)")
  premise(sum(rows(CGf, "analysis_model=='M1' & config=='G2_A_i'")$pass_pct > 5) > sum(rows(CGf, "analysis_model=='M1' & config=='G2_C_i'")$pass_pct > 5),
          "rule A (i) exceeds 5% in more boundary cells than rule C (i) under M1 (text: decision depends on the rule)")
  cpt <- rows(CGf, "analysis_model=='M1' & config=='G2_C_i' & pass_pct > 5"); premise(nrow(cpt) >= 1 && all(cpt$class == "nominal"), "rule C (i) cells above 5% (point) are Wilson-nominal under M1 (caption)")
  premise(all(rows(CGf, "analysis_model=='M1' & config=='G2_A_i' & class=='exceeding'")$pass_pct > 5), "Wilson-exceeding cells are a subset of point-exceeding cells (caption: of which)")
  p2 <- rows(T1f, "analysis_model=='M1' & config=='P2'"); premise(nrow(p2) == 16 && sum(p2$class == "exceeding") == 1, "one exceeding P2 cell under M1 (caption)")
  pm <- p2[which.max(pass_pct)]; premise(pm$pk_model == "k2020" && pm$scenario == "V2_up_080", "the exceeding P2 cell under M1 is the 2020 model V2 up (notes)")
  sets_ <- .read("config/prereg_20260926.yaml")$section4$criteria_sets
  r2s <- vapply(c("i", "ii", "iii", "iv"), function(s_) as.numeric(sets_[[s_]]$adj_r2_min), 0)
  premise(r2s[["i"]] == min(r2s) && is.null(sets_$i$span_ratio_min), "set (i) has the lowest adjusted R-squared and no span condition (caption: most lenient)")

  f <- list(
    nom = f_nominal(), ns = s08_study_days("B0", "n"), last = s08_study_days("B0", "last"), cmax = s08_cliff_max(), r2iii = f_set("iii", "r2"),
    len = headline(drange(CS, CB, "len1_median", 2, "", "cliff length (1-day definition), median, two models")),
    iii = headline(drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models")),
    i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
    a = headline(dcount(CGf, "analysis_model=='M1' & config=='G2_A_i' & pass_pct > 5", "M1 G2_A_i cells above 5% (point)")),
    aw = dcount(CGf, "analysis_model=='M1' & config=='G2_A_i' & class=='exceeding'", "M1 G2_A_i cells exceeding (Wilson lower bound above 5%)"),
    gmr = drange(TCH, "set=='i'", "true_aucinf_gmr", 3, "", "true AUC0-inf ratio failing to retained, set (i)"),
    ncell = dcount(CGf, "analysis_model=='M1' & config=='G2_A_i'", "boundary cells per configuration (M1)"),
    c = dcount(CGf, "analysis_model=='M1' & config=='G2_C_i' & pass_pct > 5", "M1 G2_C_i cells above 5% (point)"),
    cov = headline(drange(TCV, MW, "median", 1, "%", "window coverage, median, two models", scale = 100)),
    cmin = s08_cov_min(),
    n_cons = headline(dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='conservative'", "M1 P2 cells classified conservative")),
    n_all = dcount(T1f, "analysis_model=='M1' & config=='P2'", "M1 P2 boundary cells"),
    n_exc = dcount(T1f, "analysis_model=='M1' & config=='P2' & class=='exceeding'", "M1 P2 cells classified exceeding"),
    p2max = dv(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", 2, "%", "M1 P2 maximum"))

  # 도식: 왼쪽 전제 카드 -> 가운데 논거 세 줄(서로 병렬) -> 오른쪽 결론 카드. 화살표는 전제와 각 논거, 각 논거와 결론 사이
  y <- GEO$BODY_TOP + 0.05; bh <- GEO$BODY_BOTTOM - y - 0.04
  swl <- 2.1; swr <- 2.3; aw <- 0.22; ag <- 0.04; slot <- aw + 2 * ag; ins <- 0.06          # 겹친 글상자는 카드 안쪽으로 ins만큼 들여 좌우 여백을 맞춘다
  mx <- GEO$ML + swl + slot; mw <- GEO$CW - swl - swr - 2 * slot; xc <- mx + mw + slot
  T_ <- function(s_, k) tx(sprintf("S08.steps.%s.%s", s_, k), f)
  # 양쪽 카드: 라벨, 주장, 수치, 설명, 근거 슬라이드(아래)
  xi_ <- function(x) x + ins
  hh <- function(s_, k, wi) est_height(vapply(T_(s_, k), nobreak, ""), wi, 16)
  y1 <- y + 0.08; y2 <- y1 + 0.46
  side <- function(s_, x, sw, fillc, col) {
    s08_panel(c(x, y, sw, bh), fillc, "roundRect", sprintf("box_%s", s_)); wi <- sw - 2 * ins
    hcl <- hh(s_, "claim", wi); hcap <- hh(s_, "caption", wi); hr <- hh(s_, "ref", wi); yr <- y + bh - hr - 0.04
    # 수치와 설명은 주장 끝과 근거 슬라이드 사이 남는 높이의 가운데에 둔다(카드 안 빈 곳이 한쪽에 몰리지 않게)
    slack <- max(0, (yr - 0.02) - (y2 + hcl + 0.02) - 0.5 - hcap); y3 <- y2 + hcl + 0.02 + slack / 2; y4 <- y3 + 0.5
    premise(y4 + hcap <= yr + 0.02, sprintf("S08 %s card: caption above the slide reference", s_))
    deck_text(T_(s_, "label"), c(xi_(x), y1, wi, 0.42), size = 18, bold = TRUE, color = col, label = sprintf("step_%s", s_))
    deck_text(T_(s_, "claim"), c(xi_(x), y2, wi, hcl), size = 16, label = sprintf("claim_%s", s_))
    deck_text(T_(s_, "value"), c(xi_(x), y3, wi, 0.5), size = 22, bold = TRUE, color = col, label = sprintf("value_%s", s_))
    deck_text(T_(s_, "caption"), c(xi_(x), y4, wi, hcap), size = 16, color = PAL$ink2, label = sprintf("caption_%s", s_))
    deck_text(T_(s_, "ref"), c(xi_(x), yr, wi, hr), size = 16, color = PAL$muted, label = sprintf("ref_%s", s_))
  }
  side("premise", GEO$ML, swl, PAL$tint_grey, PAL$ink2)
  side("concl", xc, swr, PAL$tint_orange, PAL$orange)
  # 논거 세 줄: 왼쪽 글(라벨 + 주장, 그 아래 수치 설명), 오른쪽 수치와 근거 슬라이드. 줄 높이는 글 양에 맞추고 남는 높이는 고르게 나눈다
  args_ <- c("a1", "a2", "a3"); vw <- 2.05; lw <- mw - vw - 3 * ins; gp <- 0.1
  cl <- vapply(args_, function(s_) sprintf("__%s__ %s", T_(s_, "label"), T_(s_, "claim")), "")
  hcl <- vapply(cl, function(z) est_height(nobreak(z), lw, 16), 0); hcap <- vapply(args_, function(s_) est_height(nobreak(T_(s_, "caption")), lw, 16), 0)
  need <- pmax(hcl + hcap, 0.52 + 0.42) + 0.12; free <- bh - 2 * gp - sum(need)
  premise(free >= 0, "S08 argument rows fit the body height")
  rh <- need + free / 3; ry <- y + c(0, cumsum(rh + gp))[1:3]
  for (k in 1:3) {
    s_ <- args_[k]; yy <- ry[k]; off <- (rh[k] - hcl[k] - hcap[k]) / 2
    s08_panel(c(mx, yy, mw, rh[k]), PAL$tint_blue, "roundRect", sprintf("box_%s", s_))
    deck_text(cl[k], c(mx + ins, yy + off, lw, hcl[k]), size = 16, label = sprintf("claim_%s", s_))
    deck_text(T_(s_, "caption"), c(mx + ins, yy + off + hcl[k], lw, hcap[k]), size = 16, color = PAL$ink2, label = sprintf("caption_%s", s_))
    vo <- (rh[k] - 0.94) / 2; xv <- mx + mw - vw - ins
    deck_text(T_(s_, "value"), c(xv, yy + vo, vw, 0.5), size = 22, bold = TRUE, color = PAL$blue, label = sprintf("value_%s", s_))
    deck_text(T_(s_, "ref"), c(xv, yy + vo + 0.52, vw, 0.42), size = 16, color = PAL$muted, label = sprintf("ref_%s", s_))
    yc <- yy + rh[k] / 2 - 0.21
    s08_panel(c(GEO$ML + swl + ag, yc, aw, 0.42), PAL$muted, "rightArrow", sprintf("arrow_in_%d", k))
    s08_panel(c(mx + mw + ag, yc, aw, 0.42), PAL$muted, "rightArrow", sprintf("arrow_out_%d", k))
  }

  deck_notes(tx("S08.notes", c(f, list(
    cst = drange(CS, CB, "c_start1_median", 2, "", "cliff start concentration (mg/L), median, two models"),
    len595 = dspan(CS, CB, "len1_p05", "len1_p95", 2, "", "cliff length, 5th to 95th percentile, two models"),
    minpts = dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0)),
    cw = dcount(CGf, "analysis_model=='M1' & config=='G2_C_i' & class=='exceeding'", "M1 G2_C_i cells exceeding (Wilson lower bound above 5%)"),
    aii = dcount(CGf, "analysis_model=='M1' & config=='G2_A_ii' & pass_pct > 5", "M1 G2_A_ii cells above 5% (point)"),
    b = dcount(CGf, "analysis_model=='M1' & config=='G2_B' & pass_pct > 5", "M1 G2_B cells above 5% (point)"),
    inst = s08_inst("M1", "median"), inst0 = s08_inst("M0", "median"),
    p2ci = dci(T1f, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "M1 P2 maximum with 95% CI"),
    m0c = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='conservative'", "M0 P2 cells classified conservative"),
    m0n = dcount(T1f, "analysis_model=='M0' & config=='P2' & class=='nominal'", "M0 P2 cells classified nominal")))))
  deck_end()
}
