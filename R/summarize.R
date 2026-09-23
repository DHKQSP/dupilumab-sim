# summarize.R — 산출 지표 (지시서 §5–§6)

q05 <- function(x) quantile(x, 0.05, na.rm = TRUE, names = FALSE)
q95 <- function(x) quantile(x, 0.95, na.rm = TRUE, names = FALSE)

# 개인 수준 지표. nca: run_nca + attach_truth 결과(+ WT 등). by: 그룹 열(예: "schedule")
summarize_individual <- function(nca, by = NULL, last_planned = NULL) {
  nca[, .(
    n = .N,
    tlast_median = as.numeric(median(tlast, na.rm = TRUE)), tlast_p05 = q05(tlast), tlast_p95 = q95(tlast),
    quant_at_last_pct = if (!is.null(last_planned)) 100 * mean(tlast >= last_planned - 1.5, na.rm = TRUE) else NA_real_,
    lambda_ok_pct = 100 * mean(lambda_ok), adj_r2_ge080_pct = 100 * mean(adj_r2 >= 0.80, na.rm = FALSE), reliable_pct = 100 * mean(reliable),
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
  b <- be[be$method == method]
  per_ep <- b[, .(n_trials = .N, pass_rate = 100 * mean(pass, na.rm = TRUE), GMR_mean = mean(GMR, na.rm = TRUE), GMR_sd = sd(GMR, na.rm = TRUE),
                  logGMR_sd = sd(log(GMR), na.rm = TRUE), width_mean_pp = 100 * mean(width, na.rm = TRUE), width_p95_pp = 100 * q95(width),
                  n_R_mean = mean(n_R), n_T_mean = mean(n_T)), by = .(scenario, schedule, endpoint)]
  w <- dcast(b, scenario + schedule + trial ~ endpoint, value.var = c("pass", "GMR"))
  conc <- w[, .(
    n_trials = .N,
    pass_both_primary = 100 * mean(pass_AUClast & pass_Cmax),
    agree_last_infrel = 100 * mean(pass_AUClast == pass_AUCinf_reliable, na.rm = TRUE),
    agree_last_infall = 100 * mean(pass_AUClast == pass_AUCinf_all, na.rm = TRUE),
    last_pass_infrel_fail = 100 * mean(pass_AUClast & !pass_AUCinf_reliable, na.rm = TRUE),
    last_fail_infrel_pass = 100 * mean(!pass_AUClast & pass_AUCinf_reliable, na.rm = TRUE),
    cor_logGMR_last_infrel = cor(log(GMR_AUClast), log(GMR_AUCinf_reliable), use = "complete.obs"),
    cor_logGMR_last_true = cor(log(GMR_AUClast), log(GMR_AUCinf_true), use = "complete.obs"),
    dlogGMR_last_minus_infrel_mean = mean(log(GMR_AUClast) - log(GMR_AUCinf_reliable), na.rm = TRUE),
    dlogGMR_last_minus_infrel_sd = sd(log(GMR_AUClast) - log(GMR_AUCinf_reliable), na.rm = TRUE)
  ), by = .(scenario, schedule)]
  list(per_endpoint = per_ep, concordance = conc, wide = w)
}

# 채혈 일정 판정 규칙(지시서 §5): B0 대비
schedule_decision <- function(ind_by_sched, trial_per_ep, design, ref = "B0", scenario = "S00") {
  dr <- design$decision_rule
  base_ind <- ind_by_sched[schedule == ref]
  base_tr <- trial_per_ep[trial_per_ep$schedule == ref & trial_per_ep$scenario == scenario & trial_per_ep$endpoint == "AUClast"]
  n_ind <- ind_by_sched[schedule == ref, n]
  out <- ind_by_sched[, {
    tr <- trial_per_ep[trial_per_ep$schedule == schedule & trial_per_ep$scenario == scenario & trial_per_ep$endpoint == "AUClast"]
    width_rel <- if (nrow(tr)) 1 - tr$width_mean_pp / base_tr$width_mean_pp else NA_real_
    rel_gain <- reliable_pct - base_ind$reliable_pct
    # 외삽 20% 초과 비율: 두 비율 z 검정(같은 모집단이라 보수적)
    p1 <- extrap_gt20_pct / 100; p0 <- base_ind$extrap_gt20_pct / 100
    z <- if (p1 < p0) { pp <- (p1 + p0) / 2; (p0 - p1) / sqrt(pp * (1 - pp) * (2 / n)) } else 0
    pval <- 1 - pnorm(z)
    .(n_points = n_points, added_points = n_points - base_ind$n_points,
      ci_width_pp = if (nrow(tr)) tr$width_mean_pp else NA_real_, ci_width_rel_decrease = width_rel,
      reliable_pct = reliable_pct, reliable_gain_pp = rel_gain,
      extrap_gt20_pct = extrap_gt20_pct, extrap_gt20_pval = pval,
      crit_a = isTRUE(width_rel >= dr$ci_width_relative_decrease_min),
      crit_b = isTRUE(rel_gain >= dr$reliability_gain_pp_min),
      crit_c = isTRUE(pval < dr$extrap_gt20_test_alpha & p1 < p0))
  }, by = schedule]
  out[, recommend := crit_a | crit_b | crit_c]
  out[, marginal_visits := added_points * dr$cost_per_added_point$subjects * dr$cost_per_added_point$visits]
  out[]
}
