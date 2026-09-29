# A9 표본 수(별첨): arm당 무작위배정 130명(평가 가능 117명) 유지. 왼쪽 표 = AUClast + Cmax 동시 통과 검정력(분석식, 참 기하평균비 0.95, 2016 모델 입력, M1)
# 변동계수 격자 x arm당 평가 가능 인원(results/sample_size/ss_table_power.csv), 오른쪽 카드 = 프로토콜 변동계수(43%)와 민감도(50%)의 현행 인원 검정력과
# 목표 검정력에 필요한 인원(ss_table_n_needed.csv). 변동계수 근거: ss_inputs.csv(PK 모델), config/design_clot2021.yaml(Li 2020), fallback/sample_size_logsd.csv(Cohen 2022).
# 자료 논리: 결과보고 덱 s22_sample_size.R(보고서 5.11절과 같은 행 조건).
A9_TP <- "sample_size/ss_table_power.csv"; A9_NN <- "sample_size/ss_table_n_needed.csv"
a9_w <- function(cv, n, am = "M1") sprintf("input_model=='k2016' & cv==%s & gmr==0.95 & n==%s & analysis_model=='%s'", cv, n, am)
a9_wn <- function(cv, am, tg) sprintf("cv==%s & gmr==0.95 & analysis_model=='%s' & target_pct==%s", cv, am, tg)

slide_A9 <- function() {
  TP <- A9_TP; NN <- A9_NN; SI <- "sample_size/ss_inputs.csv"; MC <- "sample_size/ss_power_mc.csv"; LS <- "fallback/sample_size_logsd.csv"; PK <- "sample_size/ss_pk_check.csv"
  deck_slide("A9", tag = "litsim")
  L <- DK$txt$A9

  # ---- 격자와 전제 ----
  g <- rows(TP, "input_model=='k2016' & gmr==0.95 & analysis_model=='M1'")
  cvs <- sort(unique(g$cv)); ns <- sort(unique(g$n)); premise(length(cvs) == 6 && length(ns) == 4 && nrow(g) == 24, "power grid: 6 CVs x 4 evaluable n (M1, GMR 0.95, 2016 model inputs)")
  n_arm <- as.numeric(.read("config/trial_design.yaml")$n_per_arm); premise(n_arm %in% ns, "the protocol evaluable n per arm is a grid column")
  pr <- .read("config/prereg_20260926.yaml")$section3; cv_b <- pr$base_cv_pct; cv_s <- pr$sensitivity_cv_pct
  premise(cv_b %in% cvs && cv_s %in% cvs, "protocol and sensitivity CVs are grid rows")
  tg <- 90; premise(nrow(rows(NN, sprintf("gmr==0.95 & analysis_model=='M1' & target_pct==%s", tg))) == length(cvs), "n needed for the target power at every grid CV")
  at <- g[n == n_arm][order(cv)]; cv_ok <- max(at$cv[at$analytic_pct >= tg])
  premise(all(at[cv <= cv_ok]$analytic_pct >= tg) && all(at[cv > cv_ok]$analytic_pct < tg) && cv_b <= cv_ok && cv_s > cv_ok,
          "power at the protocol n reaches the target up to cv_ok, including the protocol CV, and not at the sensitivity CV (body, cards)")
  premise(all(rows(PK, "analysis_model != 'M2' & auc_ratio > 0.85 & auc_ratio < 0.92")$diff_pp < 0), "at a true GMR of 0.90 the analytic power is above the PK-model trials in every comparison (caption: optimistic)")

  f <- list(n_rand = f_n_rand(), n_arm = f_n_arm(),
            p_b = dv(TP, a9_w(cv_b, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_b)))
  y0 <- core_title(tx("A9.title", f), tx("A9.kicker"))

  # ---- 아래: 캡션(전체 폭) ----
  cap <- tx("A9.caption", list(
    g = drange(TP, "input_model=='k2016' & gmr==0.95 & analysis_model=='M1'", "gmr", 2, "", "true GMR of the power table"),
    mc_in = dcount(MC, "analytic_in_ci == TRUE", "grid cells with the analytic power inside the Wilson interval of the statistical simulation"),
    mc_n = dcount(MC, "TRUE", "grid cells simulated"),
    cv16 = dv(SI, "input_model=='k2016'", "cv_auc_pct", 1, "%", "AUClast CV, 2016 model (%)"), cv20 = dv(SI, "input_model=='k2020'", "cv_auc_pct", 1, "%", "AUClast CV, 2020 model (%)"),
    li = dcfg("design_clot2021.yaml", c("variability_checks", "li2020", "cv_pct_range"), "Li 2020 Table 3, 300 mg arms, SD/mean range (%)", function(x) rng_fmt(x[1], x[2], 0)),
    coh = dv(LS, "startsWith(variant, 'Cohen')", "AUClast_logcv", 0, "%", "Cohen 2022 implied AUClast CV (%)"),
    cvc = dv(SI, "input_model=='k2016'", "cv_cmax_pct", 1, "%", "Cmax CV, 2016 model (%)")))
  premise(any(grepl("^held at the simulated value", unlist(pr))), "Cmax CV held at the simulated value (caption)")
  premise(any(grepl("^equal to the AUClast GMR \\(primary, conservative\\)", unlist(pr))), "true Cmax ratio equal to the AUClast ratio, primary conservative assumption (caption)")
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)

  # ---- 왼쪽: 검정력 표와 본문 ----
  TW <- 7.6; hh <- 0.40
  body <- tx("A9.body", list(g = dv(TP, a9_w(cv_b, n_arm), "gmr", 2, "", "true GMR of the power table (body)"), n_arm = f$n_arm, n_rand = f$n_rand, cv_ok = dint(TP, a9_w(cv_ok, n_arm), "cv", "largest grid CV with power at the target, protocol n (%)"),
                             tgt = dint(NN, a9_wn(cv_b, "M1", tg), "target_pct", "target power (%)")))
  by <- core_body(body, capy - 0.06, width = TW)
  deck_text(tx("A9.head", list(n_rand = f$n_rand, n_arm = f$n_arm, g = dv(TP, a9_w(cv_b, n_arm), "gmr", 2, "", "true GMR of the power table (table label)"))), c(GEO$ML, y0, TW, hh), size = 16, bold = TRUE, label = "label_table", gap_pt = 0)
  tag <- function(cv) if (cv == cv_b) L$table$base else if (cv == cv_s) L$table$sens else ""
  row_lab <- function(cv) { s <- fill(L$table$row, list(cv = dint(TP, a9_w(cv, n_arm), "cv", sprintf("grid CV %s (%%)", cv)))); t_ <- tag(cv); if (nzchar(t_)) paste(s, t_) else s }
  cell <- function(cv, n) dv(TP, a9_w(cv, n), "analytic_pct", 1, "", sprintf("P2 power n %s CV %s GMR 0.95 M1", n, cv))
  df <- data.frame(a = vapply(cvs, row_lab, ""), check.names = FALSE, stringsAsFactors = FALSE)
  for (n in ns) df[[as.character(n)]] <- vapply(cvs, function(cv) cell(cv, n), "")
  nlab <- function(n) dint(TP, a9_w(cv_b, n), "n", "evaluable subjects per arm (grid column)")
  names(df) <- c(L$table$h_cv, vapply(ns, function(n) fill(if (n == n_arm) L$table$h_now else L$table$h_n, list(n = nlab(n))), ""))
  ty <- y0 + hh; th <- by - 0.06 - ty
  deck_table(df, box = c(GEO$ML, ty, TW, th), widths = c(2.3, 1.25, 1.55, 1.25, 1.25), size = 14, highlight = which(cvs %in% c(cv_b, cv_s)), highlight_fill = PAL$tint_blue, label = "table_power", pad = 2)   # 셀 위아래 여백 2pt: 표 아래 본문·캡션 세 줄 자리
  dsrc("power table", TP, "(table)")

  # ---- 오른쪽: 카드 두 개(프로토콜, 민감도), 캡션 위까지 ----
  xr <- GEO$ML + TW + 0.3; wr <- GEO$W - GEO$MR - xr; gap <- 0.2; ch <- (capy - 0.10 - y0 - gap) / 2
  nn_ <- function(cv, col = "n_evaluable_per_arm", it = "n evaluable") dint(NN, a9_wn(cv, "M1", tg), col, sprintf("%s M1 CV %s GMR 0.95 target %s", it, cv, tg))
  for (k in 1:2) {
    cv <- c(cv_b, cv_s)[k]; C <- L$cards[[c("base", "sens")[k]]]
    v <- list(cv = dint(TP, a9_w(cv, n_arm), "cv", sprintf("card CV (%%), %s", c("protocol", "sensitivity")[k])), n_arm = f$n_arm, tgt = dint(NN, a9_wn(cv, "M1", tg), "target_pct", "target power (%)"),
              p = dv(TP, a9_w(cv, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv)), ne = nn_(cv), nr = nn_(cv, "n_randomized_per_arm", "n randomized"))
    core_card(fill(C$head, v), list(list(v$p, if (k == 1) PAL$blue else PAL$ink)), fill(C$label, v), c(xr, y0 + (k - 1) * (ch + gap), wr, ch), bg = if (k == 1) PAL$tint_blue else PAL$tint_grey, value_size = 40)
  }

  # ---- 노트 ----
  s_unique <- function(rel, where, col, item) { x <- unique(rows(rel, where)[[col]]); premise(length(x) == 1, sprintf("%s [%s] :: %s has one value", rel, where, col))
    dderived(item, rel, sprintf("%s :: unique(%s)", where, col), x, fnum(x, 0, big = TRUE)) }
  pk_rng <- function(lo_, hi_) { w <- sprintf("analysis_model != 'M2' & auc_ratio > %s & auc_ratio < %s", lo_, hi_); r <- rows(PK, w); premise(nrow(r) >= 4, "PK-model check rows")
    x <- range(-r$diff_pp); dderived(sprintf("analytic minus PK-model power, true GMR %s to %s (points)", lo_, hi_), PK, sprintf("%s :: range(-diff_pp)", w), x, rng_fmt(x[1], x[2], 2)) }
  pk_abs <- function(lo_, hi_) { w <- sprintf("analysis_model != 'M2' & auc_ratio > %s & auc_ratio < %s", lo_, hi_); r <- rows(PK, w); x <- max(abs(r$diff_pp))
    dderived(sprintf("largest absolute analytic minus PK-model power, true GMR %s to %s (points)", lo_, hi_), PK, sprintf("%s :: max(abs(diff_pp))", w), x, fnum(x, 2)) }
  w90 <- sprintf("input_model=='k2016' & cv==%s & abs(gmr - 0.9) < 1e-9 & n==%s & analysis_model=='M1'", cv_b, n_arm)
  premise(nrow(rows(TP, w90)) == 1 && row1(TP, w90)$analytic_pct < tg, "power at a true GMR of 0.90 (protocol CV and n, M1) is below the target (notes)")
  deck_notes(tx("A9.notes", c(f, list(
    g90r = dv(TP, w90, "gmr", 2, "", "true GMR of the lower power row"), p90 = dv(TP, w90, "analytic_pct", 1, "%", "P2 power n 117 CV 43 GMR 0.90 M1 (analytic)"),
    p90mc = dv(TP, w90, "mc_pct", 1, "%", "P2 power n 117 CV 43 GMR 0.90 M1 (statistical simulation)"),
    g = dv(TP, a9_w(cv_b, n_arm), "gmr", 2, "", "true GMR of the power table (notes)"),
    ne_b = nn_(cv_b), nr_b = nn_(cv_b, "n_randomized_per_arm", "n randomized"), ne_s = nn_(cv_s), nr_s = nn_(cv_s, "n_randomized_per_arm", "n randomized"),
    cv_b = dint(TP, a9_w(cv_b, n_arm), "cv", "protocol CV (%)"), cv_s = dint(TP, a9_w(cv_s, n_arm), "cv", "sensitivity CV (%)"),
    p_b0 = dv(TP, a9_w(cv_b, n_arm, "M0"), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M0", cv_b)),
    p_s = dv(TP, a9_w(cv_s, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_s)),
    p_s0 = dv(TP, a9_w(cv_s, n_arm, "M0"), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M0", cv_s)),
    n130 = dint(TP, a9_w(cv_s, ns[ns > n_arm][1]), "n", "next grid n above the protocol n"),
    p130 = dv(TP, a9_w(cv_s, ns[ns > n_arm][1]), "analytic_pct", 1, "%", "P2 power, next grid n above the protocol n, sensitivity CV, M1"),
    n_b0 = dint(NN, a9_wn(cv_b, "M0", tg), "n_evaluable_per_arm", "n evaluable M0 protocol CV target"), n_s0 = dint(NN, a9_wn(cv_s, "M0", tg), "n_evaluable_per_arm", "n evaluable M0 sensitivity CV target"),
    cv16 = dv(SI, "input_model=='k2016'", "cv_auc_pct", 1, "%", "AUClast CV, 2016 model (%)"), cvc = dv(SI, "input_model=='k2016'", "cv_cmax_pct", 1, "%", "Cmax CV, 2016 model (%)"),
    rho = dv(SI, "input_model=='k2016'", "rho_total", 2, "", "correlation of log AUClast and log Cmax, 2016 model"),
    mc_tr = s_unique(MC, "TRUE", "n_trials", "trials per grid cell (statistical simulation)"),
    mc_max = local({ x <- max(abs(rows(MC)$diff_pp)); dderived("largest absolute difference, analytic versus statistical simulation (points)", MC, "max(abs(diff_pp))", x, fnum(x, 2)) }),
    g90 = drange(PK, "analysis_model != 'M2' & auc_ratio > 0.85 & auc_ratio < 0.92", "auc_ratio", 2, "", "true GMR of the PK-model check band 0.85 to 0.92"), rng90 = pk_rng(0.85, 0.92),
    g95 = drange(PK, "analysis_model != 'M2' & auc_ratio > 0.93 & auc_ratio < 1.01", "auc_ratio", 2, "", "true GMRs of the PK-model check band 0.93 to 1.01"), amax = pk_abs(0.93, 1.01)))))
  deck_end()
}
