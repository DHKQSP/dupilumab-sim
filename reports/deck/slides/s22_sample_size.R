# S22 표본 수: 참 GMR 0.95에서 CV x arm당 평가 가능 인원 검정력 표(AUC0-last + Cmax, M1, 2016 모델 입력; n = 117 열은 M0 병기),
# 분석식과 PK 모델 시험의 차이 그림(results/sample_size/ss_pk_check.csv, 분석식 - PK 모델 시험 = -diff_pp). [문헌+모의] 보고서 5.11절.
# 표의 n은 평가 가능 인원이다(130 열 = 평가 가능 130명). 격자 CV는 35/40/43/46/50/52뿐이므로 "CV 46%까지"는 격자점 기준이다.
S22_TP <- "sample_size/ss_table_power.csv"; S22_PK <- "sample_size/ss_pk_check.csv"; S22_NN <- "sample_size/ss_table_n_needed.csv"
# 보고서(5.11절 p117, nn)와 같은 행 조건 문자열(대조 검사 8이 같은 locator를 찾는다)
s22_w <- function(cv, n, am = "M1") sprintf("input_model=='k2016' & cv==%s & gmr==0.95 & n==%s & analysis_model=='%s'", cv, n, am)
s22_wn <- function(cv, am, tg) sprintf("cv==%s & gmr==0.95 & analysis_model=='%s' & target_pct==%s", cv, am, tg)
# 분석식 - PK 모델 시험 검정력(보고서 ss_pk와 같은 조건과 계산): 참 GMR 구간의 범위 또는 절댓값 최대
s22_pk <- function(lo_, hi_, what = c("range", "absmax")) {
  what <- match.arg(what); w <- sprintf("analysis_model != 'M2' & auc_ratio > %s & auc_ratio < %s", lo_, hi_); r <- rows(S22_PK, w)
  premise(nrow(r) >= 4, "PK-model check rows in the GMR band"); d <- -r$diff_pp
  if (what == "absmax") { x <- max(abs(d)); return(dderived(sprintf("largest absolute analytic minus PK-model power, true GMR %s to %s (points)", lo_, hi_), S22_PK, sprintf("%s :: max(abs(diff_pp))", w), x, fnum(x, 2))) }
  x <- range(d); dderived(sprintf("analytic minus PK-model power, true GMR %s to %s (points)", lo_, hi_), S22_PK, sprintf("%s :: range(-diff_pp)", w), x, rng_fmt(x[1], x[2], 2))
}

# 여러 행에서 같은 값(시행 횟수 등)을 천 단위 쉼표로
s22_unique <- function(rel, where, col, item) {
  x <- unique(rows(rel, where)[[col]]); premise(length(x) == 1, sprintf("%s [%s] :: %s has one value", rel, where, col))
  dderived(item, rel, sprintf("%s :: unique(%s)", where, col), x, fnum(x, 0, big = TRUE))
}

slide_S22 <- function() {
  TP <- S22_TP; PK <- S22_PK; NN <- S22_NN; SI <- "sample_size/ss_inputs.csv"; MC <- "sample_size/ss_power_mc.csv"; LS <- "fallback/sample_size_logsd.csv"
  deck_slide("S22", tag = "litsim")
  L <- DK$txt$S22

  # ---- 격자와 전제 ----
  g <- rows(TP, "input_model=='k2016' & gmr==0.95 & analysis_model=='M1'")
  cvs <- sort(unique(g$cv)); ns <- sort(unique(g$n)); premise(length(cvs) == 6 && length(ns) == 4 && nrow(g) == 24, "power grid: 6 CVs x 4 evaluable n (M1, GMR 0.95, 2016 model inputs)")
  n_arm <- as.numeric(.read("config/trial_design.yaml")$n_per_arm); premise(n_arm %in% ns, "the protocol evaluable n per arm is a grid column")
  tg <- 90; tgt_rows <- rows(NN, sprintf("gmr==0.95 & analysis_model=='M1' & target_pct==%s", tg)); premise(nrow(tgt_rows) == 6, "n needed for the target power at every grid CV")
  at <- g[n == n_arm][order(cv)]
  cv_ok <- max(at$cv[at$analytic_pct >= tg])
  premise(all(at[cv <= cv_ok]$analytic_pct >= tg) && all(at[cv > cv_ok]$analytic_pct < tg), "power at the protocol n is at least the target at every grid CV up to cv_ok and below it above (title)")
  pr <- .read("config/prereg_20260926.yaml")$section3; cv_b <- pr$base_cv_pct; cv_s <- pr$sensitivity_cv_pct
  premise(cv_b %in% cvs && cv_s %in% cvs && cv_s > cv_ok && cv_b <= cv_ok, "base and sensitivity CVs are grid rows; the sensitivity CV lies above cv_ok (title)")
  w90 <- "analysis_model != 'M2' & auc_ratio > 0.85 & auc_ratio < 0.92"
  premise(all(rows(PK, w90)$diff_pp < 0), "at a true GMR of 0.90 the analytic power is above the PK-model trials in every comparison (text: optimistic)")

  f <- list(n_rand = f_n_rand(), n_arm = f_n_arm(),
            g = drange(TP, "input_model=='k2016' & gmr==0.95 & analysis_model=='M1'", "gmr", 2, "", "true GMR of the power table"),
            tgt = dint(NN, s22_wn(cv_b, "M1", tg), "target_pct", "target power (%)"),
            cv_ok = dint(TP, s22_w(cv_ok, n_arm), "cv", "largest grid CV with power at the target, protocol n (%)"),
            p_ok = dv(TP, s22_w(cv_ok, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_ok)),
            cv_s = dint(TP, s22_w(cv_s, n_arm), "cv", "sensitivity CV (%)"),
            p_s = dv(TP, s22_w(cv_s, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_s)),
            cv_b = dint(TP, s22_w(cv_b, n_arm), "cv", "protocol CV (%)"),
            p_b = dv(TP, s22_w(cv_b, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_b)))
  premise(as.numeric(sub("%", "", f$p_b)) >= tg && as.numeric(sub("%", "", f$p_s)) < tg, "protocol CV power at or above the target, sensitivity CV power below it (bullet 1, notes)")
  deck_kicker(tx("S22.kicker")); deck_title(tx("S22.title", f))

  # ---- 왼쪽: 검정력 표 ----
  tw <- 7.6; y0 <- GEO$BODY_TOP
  deck_text(tx("S22.caption", list(g = f$g)), c(GEO$ML, y0, tw, 0.4), size = 16, color = PAL$ink2, label = "caption", gap_pt = 0)
  tag <- function(cv) if (cv == cv_b) L$table$base else if (cv == cv_s) L$table$sens else ""
  row_lab <- function(cv) { s <- fill(L$table$row, list(cv = dint(TP, s22_w(cv, n_arm), "cv", sprintf("grid CV %s (%%)", cv)))); t_ <- tag(cv); if (nzchar(t_)) paste(s, t_) else s }
  cell <- function(cv, n, am = "M1") dv(TP, s22_w(cv, n, am), "analytic_pct", 1, "", sprintf("P2 power n %s CV %s GMR 0.95 %s", n, cv, am))
  df <- data.frame(a = vapply(cvs, row_lab, ""), check.names = FALSE, stringsAsFactors = FALSE)
  for (n in ns) df[[as.character(n)]] <- vapply(cvs, function(cv) cell(cv, n), "")
  df[["m0"]] <- vapply(cvs, function(cv) cell(cv, n_arm, "M0"), "")
  nlab <- function(n) dint(TP, s22_w(cv_b, n), "n", "evaluable subjects per arm (grid column)")
  names(df) <- c(L$table$h_cv, vapply(ns, function(n) fill(if (n == n_arm) L$table$h_now else L$table$h_n, list(n = nlab(n))), ""),
                 fill(L$table$h_m0, list(n = nlab(n_arm))))
  yt <- y0 + 0.42; TH <- 2.5
  deck_table(df, box = c(GEO$ML, yt, tw, TH), widths = c(1.85, 1.05, 1.45, 1.05, 1.05, 1.15), size = 13, highlight = which(cvs %in% c(cv_b, cv_s)))

  # ---- 왼쪽 아래: 요점 ----
  nn_ <- function(cv, am, col = "n_evaluable_per_arm", it = "n") dint(NN, s22_wn(cv, am, tg), col, sprintf("%s %s CV %s GMR 0.95 target %s", it, am, cv, tg))
  b <- list(tgt = f$tgt, cv_b = f$cv_b, cv_s = f$cv_s, cv_ok = f$cv_ok,
            n_b = nn_(cv_b, "M1"), n_s = nn_(cv_s, "M1"), nr_s = nn_(cv_s, "M1", "n_randomized_per_arm", "n randomized"), nr_b = nn_(cv_b, "M1", "n_randomized_per_arm", "n randomized"),
            cv16 = dv(SI, "input_model=='k2016'", "cv_auc_pct", 1, "%", "AUC0-last CV, 2016 model (%)"),
            cv20 = dv(SI, "input_model=='k2020'", "cv_auc_pct", 1, "%", "AUC0-last CV, 2020 model (%)"),
            li = dcfg("design_clot2021.yaml", c("variability_checks", "li2020", "cv_pct_range"), "Li 2020 Table 3, 300 mg arms, SD/mean range (%)", function(x) rng_fmt(x[1], x[2], 0)),
            coh = dv(LS, "startsWith(variant, 'Cohen')", "AUClast_logcv", 0, "%", "Cohen 2022 implied AUC0-last CV (%)"))
  yb <- yt + TH + 0.1
  deck_bullets(tx("S22.bullets", b), c(GEO$ML, yb, tw, GEO$BODY_BOTTOM - yb), size = 16, gap_pt = 5)

  # ---- 오른쪽: 분석식 - PK 모델 시험 그림(참 GMR별, 두 모델, M1과 M0) ----
  d <- rows(PK, "analysis_model != 'M2' & scenario != 'F097'")[, .(pk_model, scenario, auc_ratio, analysis_model, empirical_pct, lo, hi, analytic_pct)]
  premise(nrow(d) == 12 && !anyNA(d$analytic_pct), "PK-model check: 3 true GMRs x 2 models x M0/M1")
  d[, `:=`(dd = analytic_pct - empirical_pct, dlo = analytic_pct - hi, dhi = analytic_pct - lo)]
  xl <- d[, .(r = mean(auc_ratio)), by = scenario][order(r)]; xl[, lab := fnum(r, 2)]; premise(!anyDuplicated(xl$lab), "distinct true GMR labels")
  d[, x := factor(xl$lab[match(scenario, xl$scenario)], levels = xl$lab)]
  d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  d[, am := factor(unlist(L$fig$am)[analysis_model], levels = unlist(L$fig$am)[c("M1", "M0")])]
  pd <- position_dodge(width = 0.5)
  p <- ggplot(d, aes(x = x, y = dd, colour = model, shape = model)) + geom_hline(yintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_errorbar(aes(ymin = dlo, ymax = dhi), position = pd, width = 0.25, linewidth = 0.55) + geom_point(position = pd, size = 3) +
    facet_wrap(~am, nrow = 1) + scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = unname(MODEL_SHAPE)) +
    labs(x = L$fig$xlab, y = NULL, subtitle = fill(L$fig$ylab, list(n = f$n_arm))) + theme_deck(12) +
    theme(panel.grid.major.x = element_blank(), legend.position = "top", legend.margin = margin(0, 0, 0, 0), plot.subtitle = element_text(colour = PAL$ink2, size = 12, lineheight = 1.1, margin = margin(0, 0, 3, 0)))
  xr <- GEO$ML + tw + 0.3; wr <- GEO$W - GEO$MR - xr
  hf <- 3.6
  deck_figure(p, "s22_analytic_vs_pk", c(xr, GEO$BODY_TOP, wr, hf), src = PK)
  f2 <- list(g90 = drange(PK, w90, "auc_ratio", 2, "", "true GMR of the PK-model check band 0.85 to 0.92"), rng = s22_pk(0.85, 0.92, "range"),
             g95 = drange(PK, "analysis_model != 'M2' & auc_ratio > 0.93 & auc_ratio < 1.01", "auc_ratio", 2, "", "true GMRs of the PK-model check band 0.93 to 1.01"),
             amax = s22_pk(0.93, 1.01, "absmax"))
  yc <- GEO$BODY_TOP + hf + 0.12
  deck_text(tx("S22.caution", f2), c(xr, yc, wr, GEO$BODY_BOTTOM - yc), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "caution", gap_pt = 0)

  # ---- 노트 ----
  deck_notes(tx("S22.notes", c(f, b, f2, list(
    cvc = dv(SI, "input_model=='k2016'", "cv_cmax_pct", 1, "%", "Cmax CV, 2016 model (%)"),
    rho = dv(SI, "input_model=='k2016'", "rho_total", 2, "", "correlation of log AUC0-last and log Cmax, 2016 model"),
    r2s = dv(SI, "input_model=='k2016'", "r2_stratum_auc", 1, "%", "share of AUC0-last variance explained by the weight stratum (%)", scale = 100),
    mc_in = dcount(MC, "analytic_in_ci == TRUE", "grid cells with the analytic power inside the Wilson interval of the statistical simulation"),
    mc_n = dcount(MC, "TRUE", "grid cells simulated"),
    mc_tr = s22_unique(MC, "TRUE", "n_trials", "trials per grid cell (statistical simulation)"),
    mc_max = local({ x <- max(abs(rows(MC)$diff_pp)); dderived("largest absolute difference, analytic versus statistical simulation (points)", MC, "max(abs(diff_pp))", x, fnum(x, 2)) }),
    n_b0 = nn_(cv_b, "M0"), nr_b0 = nn_(cv_b, "M0", "n_randomized_per_arm", "n randomized"), nr_b = nn_(cv_b, "M1", "n_randomized_per_arm", "n randomized"),
    n_s0 = nn_(cv_s, "M0"), nr_s0 = nn_(cv_s, "M0", "n_randomized_per_arm", "n randomized"),
    p_b = dv(TP, s22_w(cv_b, n_arm), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M1", cv_b)),
    p_b0 = dv(TP, s22_w(cv_b, n_arm, "M0"), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M0", cv_b)),
    p_s0 = dv(TP, s22_w(cv_s, n_arm, "M0"), "analytic_pct", 1, "%", sprintf("P2 power n 117 CV %s GMR 0.95 M0", cv_s)),
    p130 = dv(TP, s22_w(cv_s, ns[ns > n_arm][1]), "analytic_pct", 1, "%", "P2 power, next grid n above the protocol n, sensitivity CV, M1"),
    n130 = dint(TP, s22_w(cv_s, ns[ns > n_arm][1]), "n", "next grid n above the protocol n"),
    ntr = s22_unique(PK, "analysis_model != 'M2' & scenario %in% c('F_down_090','F_down_095')", "n_trials", "PK-model trials per check cell, F-decrease cells (true GMR 0.90 and 0.95)"),
    ntr2 = s22_unique(PK, "analysis_model != 'M2' & scenario %in% c('S00','F097')", "n_trials", "PK-model trials per check cell, S00 and F097")))))
  deck_end()
}
