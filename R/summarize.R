# summarize.R — 산출 지표 (지시서 §5–§6)

q05 <- function(x) quantile(x, 0.05, na.rm = TRUE, names = FALSE)
q95 <- function(x) quantile(x, 0.95, na.rm = TRUE, names = FALSE)

# 개인 수준 지표. nca: run_nca + attach_truth 결과(+ WT 등). by: 그룹 열(예: "schedule")
summarize_individual <- function(nca, by = NULL, last_planned = NULL) {
  nca[, .(
    n = .N,
    tlast_median = as.numeric(median(tlast, na.rm = TRUE)), tlast_p05 = q05(tlast), tlast_p95 = q95(tlast),
    quant_at_last_pct = if (!is.null(last_planned)) 100 * mean(tlast >= last_planned - 1.5, na.rm = TRUE) else NA_real_,
    lambda_ok_pct = 100 * mean(lambda_ok), adj_r2_ge080_pct = 100 * mean(!is.na(adj_r2) & adj_r2 >= 0.80), reliable_pct = 100 * mean(reliable),
    n_lambda_median = as.numeric(median(n_lambda, na.rm = TRUE)), n_lambda_eq3_pct = 100 * mean(n_lambda == 3, na.rm = TRUE),
    extrap_median = median(pct_extrap, na.rm = TRUE), extrap_p95 = q95(pct_extrap), extrap_gt20_pct = 100 * mean(pct_extrap > 20, na.rm = TRUE),
    extrap_true_median = median(pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(pct_extrap_true), coverage_lt80_pct = 100 * mean(coverage_true < 0.80, na.rm = TRUE),
    err_tlast_median = median(err_AUClast_vs_true_tlast, na.rm = TRUE), err_tlast_sd = sd(err_AUClast_vs_true_tlast, na.rm = TRUE),
    err_inf_sd = sd(err_AUClast_vs_true_inf, na.rm = TRUE), err_inf_median = median(err_AUClast_vs_true_inf, na.rm = TRUE),
    AUClast_geo = geo_mean(AUClast), AUClast_logcv = log_cv_pct(AUClast), AUClast_mean = mean(AUClast, na.rm = TRUE),
    Cmax_geo = geo_mean(Cmax), Cmax_logcv = log_cv_pct(Cmax), Cmax_mean = mean(Cmax, na.rm = TRUE),
    AUCinf_geo = geo_mean(AUCinf[reliable]), AUCinf_true_geo = geo_mean(AUCinf_true)
  ), by = by]
}

# λz 점 수 분포
summarize_n_lambda <- function(nca, by = "schedule") nca[, .N, by = c(by, "n_lambda")][order(get(by), n_lambda)][, pct := 100 * N / sum(N), by = by][]

# 시험 수준 지표. be: run_trials()$be (trial, scenario, schedule, endpoint, method, GMR, CI, pass, width, n_R, n_T)
summarize_trials <- function(be, method = "pooled_t") {
  m_ <- method                                  # data.table 스코프: 인자명이 열이름과 같으면 열로 해석되므로 지역 변수로 옮긴다(D-022)
  b <- be[be[["method"]] == m_]
  per_ep <- b[, .(n_trials = .N, pass_rate = 100 * mean(pass, na.rm = TRUE), GMR_mean = mean(GMR, na.rm = TRUE), GMR_sd = sd(GMR, na.rm = TRUE),
                  logGMR_sd = sd(log(GMR), na.rm = TRUE), width_mean_pp = 100 * mean(width, na.rm = TRUE), width_p95_pp = 100 * q95(width),
                  n_R_mean = mean(n_R), n_T_mean = mean(n_T)), by = .(scenario, schedule, endpoint)]
  w <- dcast(b, scenario + schedule + trial ~ endpoint, value.var = c("pass", "GMR"))
  for (ep in c("Cmax", "AUClast", "AUCinf_all", "AUCinf_reliable", "AUCinf_true")) {      # 없는 평가변수는 NA로 채운다
    if (!paste0("pass_", ep) %in% names(w)) w[, paste0("pass_", ep) := NA]
    if (!paste0("GMR_", ep) %in% names(w)) w[, paste0("GMR_", ep) := NA_real_]
  }
  safe_cor <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) NA_real_ else cor(x[ok], y[ok]) }
  conc <- w[, .(
    n_trials = .N,
    pass_both_primary = 100 * mean(pass_AUClast & pass_Cmax, na.rm = TRUE),
    agree_last_infrel = 100 * mean(pass_AUClast == pass_AUCinf_reliable, na.rm = TRUE),
    agree_last_infall = 100 * mean(pass_AUClast == pass_AUCinf_all, na.rm = TRUE),
    last_pass_infrel_fail = 100 * mean(pass_AUClast & !pass_AUCinf_reliable, na.rm = TRUE),
    last_fail_infrel_pass = 100 * mean(!pass_AUClast & pass_AUCinf_reliable, na.rm = TRUE),
    cor_logGMR_last_infrel = safe_cor(log(GMR_AUClast), log(GMR_AUCinf_reliable)),
    cor_logGMR_last_true = safe_cor(log(GMR_AUClast), log(GMR_AUCinf_true)),
    dlogGMR_last_minus_infrel_mean = mean(log(GMR_AUClast) - log(GMR_AUCinf_reliable), na.rm = TRUE),
    dlogGMR_last_minus_infrel_sd = sd(log(GMR_AUClast) - log(GMR_AUCinf_reliable), na.rm = TRUE)
  ), by = .(scenario, schedule)]
  list(per_endpoint = per_ep, concordance = conc, wide = w)
}

# 개인 수준 대응 비교(같은 대상자, 공통 난수): 일정 s vs 기준 일정. 이진 지표는 정확 McNemar, 연속 지표는 대상자별 차이.
paired_individual_vs_ref <- function(nca, ref = "B0") {
  ref_ <- ref
  base <- nca[schedule == ref_, .(id, r0 = reliable, x0 = !is.na(pct_extrap) & pct_extrap > 20, l0 = lambda_ok, e0 = pct_extrap_true, a0 = err_AUClast_vs_true_inf)]
  rbindlist(lapply(setdiff(unique(nca$schedule), ref_), function(sh) {
    m <- merge(base, nca[schedule == sh, .(id, r1 = reliable, x1 = !is.na(pct_extrap) & pct_extrap > 20, l1 = lambda_ok, e1 = pct_extrap_true, a1 = err_AUClast_vs_true_inf)], by = "id")
    mcn <- function(b, c) if (b + c == 0) 1 else binom.test(min(b, c), b + c, 0.5)$p.value
    b_x <- sum(m$x0 & !m$x1); c_x <- sum(!m$x0 & m$x1)          # b: B0에서만 >20%, c: s에서만 >20%
    b_r <- sum(!m$r0 & m$r1); c_r <- sum(m$r0 & !m$r1)          # b: s에서만 신뢰 충족
    dr <- m$r1 - m$r0
    data.table(schedule = sh, n = nrow(m),
               reliable_gain_pp = 100 * mean(dr), reliable_gain_lo = 100 * (mean(dr) - 1.96 * sd(dr) / sqrt(nrow(m))), reliable_gain_hi = 100 * (mean(dr) + 1.96 * sd(dr) / sqrt(nrow(m))),
               reliable_mcnemar_p = mcn(b_r, c_r),
               extrap_gt20_change_pp = 100 * (mean(m$x1) - mean(m$x0)), extrap_gt20_only_ref = b_x, extrap_gt20_only_sched = c_x,
               extrap_gt20_mcnemar_p = mcn(b_x, c_x),
               lambda_ok_change_pp = 100 * (mean(m$l1) - mean(m$l0)),
               extrap_true_median_change = median(m$e1 - m$e0, na.rm = TRUE),
               auclast_err_inf_sd_ratio = sd(m$a1, na.rm = TRUE) / sd(m$a0, na.rm = TRUE))
  }))
}

# 시험 수준 대응 비교: 같은 시험의 일정 s vs 기준 일정 AUClast CI 폭 상대 변화(평균과 95% CI)
paired_trial_width_vs_ref <- function(be, ref = "B0", scenario = "S00", endpoint = "AUClast", method = "pooled_t") {
  sc_ <- scenario; ep_ <- endpoint; m_ <- method; ref_ <- ref   # data.table 스코프 회피(D-022)
  b <- be[be[["scenario"]] == sc_ & be[["endpoint"]] == ep_ & be[["method"]] == m_, .(trial, schedule, width, pass)]
  if (anyDuplicated(b[, .(trial, schedule)])) stop("paired_trial_width_vs_ref: (trial, schedule) 중복 — 필터 확인")
  w0 <- b[schedule == ref_, .(trial, w0 = width, p0 = pass)]
  rbindlist(lapply(setdiff(unique(b$schedule), ref_), function(sh) {
    m <- merge(w0, b[schedule == sh, .(trial, w1 = width, p1 = pass)], by = "trial")
    rel <- 1 - m$w1 / m$w0
    data.table(schedule = sh, scenario = sc_, n_trials = nrow(m), width_rel_decrease = mean(rel),
               width_rel_decrease_lo = mean(rel) - 1.96 * sd(rel) / sqrt(nrow(m)), width_rel_decrease_hi = mean(rel) + 1.96 * sd(rel) / sqrt(nrow(m)),
               pass_change_pp = 100 * (mean(m$p1) - mean(m$p0)), pass_discordant = sum(m$p1 != m$p0))
  }))
}

# 채혈 일정 판정 규칙(지시서 §5, SPEC §7.3): B0 대비. 입력은 위 두 대응 비교 표.
schedule_decision <- function(paired_ind, paired_trial, design, n_points) {
  dr <- design$decision_rule
  out <- merge(paired_ind, paired_trial[scenario == "S00", .(schedule, width_rel_decrease, width_rel_decrease_lo, width_rel_decrease_hi)], by = "schedule", all.x = TRUE)
  out <- merge(out, n_points, by = "schedule", all.x = TRUE)
  out[, `:=`(
    crit_a_width = !is.na(width_rel_decrease_lo) & width_rel_decrease >= dr$ci_width_relative_decrease_min & width_rel_decrease_lo > 0,
    crit_b_reliable = reliable_gain_pp >= dr$reliability_gain_pp_min,
    crit_c_extrap = extrap_gt20_change_pp < 0 & extrap_gt20_mcnemar_p < dr$extrap_gt20_test_alpha)]
  out[, recommend := crit_a_width | crit_b_reliable | crit_c_extrap]
  out[, added_visits_total := added_points * dr$cost_per_added_point$subjects * dr$cost_per_added_point$visits]
  out[]
}
