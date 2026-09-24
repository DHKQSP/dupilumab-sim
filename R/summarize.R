# summarize.R — 산출 지표 (지시서 §5–§6)

q05 <- function(x) quantile(x, 0.05, na.rm = TRUE, names = FALSE)

# Monte Carlo 정밀도 (검토 의견 통합본 §4-1)
# 비율: Wilson 95% 구간(%). x = 성공 수, n = 반복 수
wilson_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n
  den <- 1 + z^2 / n; ctr <- (p + z^2 / (2 * n)) / den; hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.table(est = 100 * p, lo = 100 * pmax(0, ctr - hw), hi = 100 * pmin(1, ctr + hw))
}
# 평균(예: GMR 평균): Monte Carlo 표준오차 기반 95% 구간
mc_mean_ci <- function(x, conf = 0.95) { x <- x[is.finite(x)]; z <- qnorm(1 - (1 - conf) / 2); m <- mean(x); se <- sd(x) / sqrt(length(x)); data.table(est = m, lo = m - z * se, hi = m + z * se, se = se, n = length(x)) }
# 쌍대 비율 차이(같은 시험의 두 이진 판정): 평균 차이와 95% 구간(%p)
paired_prop_diff_ci <- function(a, b, conf = 0.95) { d <- as.numeric(a) - as.numeric(b); d <- d[is.finite(d)]; z <- qnorm(1 - (1 - conf) / 2); m <- mean(d); se <- sd(d) / sqrt(length(d)); data.table(est = 100 * m, lo = 100 * (m - z * se), hi = 100 * (m + z * se), n = length(d)) }

# 탈락자(신뢰 기준 미충족)와 잔류자의 체중 비교(§6-2)
summarize_dropout <- function(nca, by = NULL) {
  nca[, .(n_reliable = sum(reliable %in% TRUE), n_dropout = sum(!(reliable %in% TRUE)),
          wt_retained_mean = mean(WT[reliable %in% TRUE]), wt_dropout_mean = if (any(!(reliable %in% TRUE))) mean(WT[!(reliable %in% TRUE)]) else NA_real_,
          wt_dropout_median = if (any(!(reliable %in% TRUE))) median(WT[!(reliable %in% TRUE)]) else NA_real_), by = by]
}
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
               extrap_gt20_ratio = if (sum(m$x0) > 0) sum(m$x1) / sum(m$x0) else NA_real_,
               lambda_ok_change_pp = 100 * (mean(m$l1) - mean(m$l0)),
               extrap_true_median_change = median(m$e1 - m$e0, na.rm = TRUE),
               auclast_err_inf_sd_ratio = sd(m$a1, na.rm = TRUE) / sd(m$a0, na.rm = TRUE))
  }))
}

# 판정 (c)(d)의 몬테카를로 불확실성: 대응 2×2 칸(대상자 수)을 다항 재표집해 신뢰 충족률 차이와 외삽 20% 초과 비율의 비에 대한 95% 구간(D-028)
paired_bootstrap_cd <- function(nca, ref = "B0", B = 2000, seed = 20260924L) {
  ref_ <- ref
  base <- nca[schedule == ref_, .(id, r0 = reliable, x0 = !is.na(pct_extrap) & pct_extrap > 20)]
  rbindlist(lapply(setdiff(unique(nca$schedule), ref_), function(sh) {
    m <- merge(base, nca[schedule == sh, .(id, r1 = reliable, x1 = !is.na(pct_extrap) & pct_extrap > 20)], by = "id")
    n <- nrow(m)
    cr <- c(sum(!m$r0 & !m$r1), sum(!m$r0 & m$r1), sum(m$r0 & !m$r1), sum(m$r0 & m$r1))   # 00, 01, 10, 11
    cx <- c(sum(!m$x0 & !m$x1), sum(!m$x0 & m$x1), sum(m$x0 & !m$x1), sum(m$x0 & m$x1))
    bs <- with_seed(derive_seed(seed, sh, "cd"), {
      R <- rmultinom(B, n, cr / n); X <- rmultinom(B, n, cx / n)
      list(gain = 100 * (R[2, ] - R[3, ]) / n, ratio = (X[2, ] + X[4, ]) / pmax(X[3, ] + X[4, ], 1))
    })
    data.table(schedule = sh, c_gain_boot_lo = quantile(bs$gain, 0.025, names = FALSE), c_gain_boot_hi = quantile(bs$gain, 0.975, names = FALSE),
               d_ratio_boot_lo = quantile(bs$ratio, 0.025, names = FALSE), d_ratio_boot_hi = quantile(bs$ratio, 0.975, names = FALSE),
               d_n_ref = cx[3] + cx[4], d_n_sched = cx[2] + cx[4], d_n_subjects = n)
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
               pass_change_pp = 100 * (mean(m$p1) - mean(m$p0)), pass_discordant = sum(m$p1 != m$p0),
               pass_change_lo = paired_prop_diff_ci(m$p1, m$p0)$lo, pass_change_hi = paired_prop_diff_ci(m$p1, m$p0)$hi)
  }))
}

# 채혈 일정 판정 규칙(검토 의견 3차 §5, SPEC §7.3, D-026). B0 대비, 하나 이상 충족이면 추가 채혈 권고.
#  (a) AUClast 90% CI 평균 폭 상대 감소 ≥ 2% (S00): 1 − mean(폭_s)/mean(폭_B0). 대응 비교 95% CI 병기.
#  (b) ke ×1.10(KE110)에서 AUClast 통과율 +2%p 이상.
#  (c) AUCinf 신뢰 기준 충족률 +5%p 이상(같은 20,000명).
#  (d) 비구획 외삽 20% 초과 비율이 B0의 절반 이하.
# 입력: ind(summarize_individual 일정별 표), paired_ind(paired_individual_vs_ref), per_ep(summarize_trials()$per_endpoint, pooled_t),
#       paired_tr(paired_trial_width_vs_ref, 시나리오별), n_points(schedule, n_points, added_points)
schedule_decision <- function(ind, paired_ind, per_ep, paired_tr, design, n_points, ref = "B0") {
  dr <- design$decision_rule; ref_ <- ref
  ind_ref <- ind[schedule == ref_]
  w_ref <- per_ep[scenario == "S00" & endpoint == "AUClast" & schedule == ref_, width_mean_pp]
  ke_ref <- per_ep[scenario == "KE110" & endpoint == "AUClast" & schedule == ref_, pass_rate]
  scheds <- setdiff(ind$schedule, ref_)
  out <- rbindlist(lapply(scheds, function(sh) {
    w_s <- per_ep[scenario == "S00" & endpoint == "AUClast" & schedule == sh, width_mean_pp]
    ke_s <- per_ep[scenario == "KE110" & endpoint == "AUClast" & schedule == sh, pass_rate]
    pt_s <- paired_tr[scenario == "S00" & schedule == sh]
    pk_s <- paired_tr[scenario == "KE110" & schedule == sh]
    pi_s <- paired_ind[schedule == sh]
    i_s <- ind[schedule == sh]
    a <- if (length(w_s) && length(w_ref)) 1 - w_s / w_ref else NA_real_
    b <- if (length(ke_s) && length(ke_ref)) ke_s - ke_ref else NA_real_
    d_ratio <- if (ind_ref$extrap_gt20_pct > 0) i_s$extrap_gt20_pct / ind_ref$extrap_gt20_pct else NA_real_
    data.table(schedule = sh,
      a_width_B0_pp = if (length(w_ref)) w_ref else NA_real_, a_width_pp = if (length(w_s)) w_s else NA_real_, a_mean_width_rel_decrease = a,
      a_paired_lo = if (nrow(pt_s)) pt_s$width_rel_decrease_lo else NA_real_, a_paired_hi = if (nrow(pt_s)) pt_s$width_rel_decrease_hi else NA_real_,
      b_ke110_pass_B0 = if (length(ke_ref)) ke_ref else NA_real_, b_ke110_pass = if (length(ke_s)) ke_s else NA_real_, b_ke110_gain_pp = b,
      b_discordant_trials = if (nrow(pk_s)) pk_s$pass_discordant else NA_integer_,
      b_paired_lo = if (nrow(pk_s) && "pass_change_lo" %in% names(pk_s)) pk_s$pass_change_lo else NA_real_, b_paired_hi = if (nrow(pk_s) && "pass_change_hi" %in% names(pk_s)) pk_s$pass_change_hi else NA_real_,
      c_reliable_B0 = ind_ref$reliable_pct, c_reliable = i_s$reliable_pct, c_reliable_gain_pp = pi_s$reliable_gain_pp,
      d_extrap20_B0 = ind_ref$extrap_gt20_pct, d_extrap20 = i_s$extrap_gt20_pct, d_extrap20_ratio = d_ratio, d_mcnemar_p = pi_s$extrap_gt20_mcnemar_p,
      c_gain_boot_lo = if ("c_gain_boot_lo" %in% names(pi_s)) pi_s$c_gain_boot_lo else NA_real_, c_gain_boot_hi = if ("c_gain_boot_hi" %in% names(pi_s)) pi_s$c_gain_boot_hi else NA_real_,
      d_ratio_boot_lo = if ("d_ratio_boot_lo" %in% names(pi_s)) pi_s$d_ratio_boot_lo else NA_real_, d_ratio_boot_hi = if ("d_ratio_boot_hi" %in% names(pi_s)) pi_s$d_ratio_boot_hi else NA_real_,
      d_n_ref = if ("d_n_ref" %in% names(pi_s)) pi_s$d_n_ref else NA_integer_, d_n_sched = if ("d_n_sched" %in% names(pi_s)) pi_s$d_n_sched else NA_integer_,
      d_abs_change_per_arm = if ("d_n_subjects" %in% names(pi_s)) (pi_s$d_n_sched - pi_s$d_n_ref) / pi_s$d_n_subjects * design$n_per_arm else NA_real_)
  }))
  out[, `:=`(crit_a = a_mean_width_rel_decrease >= dr$ci_width_mean_rel_decrease_min,
             crit_b = b_ke110_gain_pp >= dr$ke110_pass_gain_pp_min,
             crit_c = c_reliable_gain_pp >= dr$reliability_gain_pp_min,
             crit_d = d_extrap20_ratio <= dr$extrap_gt20_ratio_max)]
  # 권고: 하나라도 TRUE면 권고. 판정 불가(NA) 기준은 권고 근거로 쓰지 않되 표에 NA로 남긴다.
  out[, recommend := (crit_a %in% TRUE) | (crit_b %in% TRUE) | (crit_c %in% TRUE) | (crit_d %in% TRUE)]
  out[, a_label := fifelse(is.na(a_paired_lo), NA_character_, fifelse(a_paired_lo <= 0 & a_paired_hi >= 0, "변화 없음", fifelse(a_mean_width_rel_decrease > 0, "감소", "확대")))]   # §2
  out[, n_criteria_evaluable := (!is.na(crit_a)) + (!is.na(crit_b)) + (!is.na(crit_c)) + (!is.na(crit_d))]
  out <- merge(n_points, out, by = "schedule", all.y = TRUE)
  out[, added_visits_total := added_points * dr$cost_per_added_point$subjects * dr$cost_per_added_point$visits]
  setcolorder(out, c("schedule", "n_points", "added_points", "added_visits_total"))
  out[]
}

# 시험 수준 요약 + Monte Carlo 95% 구간(§4-1): 통과율·일치율은 Wilson, GMR 평균은 MC SE 구간. 반복 수 명시.
summarize_trials_ci <- function(be, method = "pooled_t") {
  m_ <- method
  b <- be[be[["method"]] == m_]
  pe <- b[, {
    w <- wilson_ci(sum(pass %in% TRUE), sum(!is.na(pass))); g <- mc_mean_ci(GMR)
    .(n_trials = .N, pass_rate = w$est, pass_lo = w$lo, pass_hi = w$hi, GMR_mean = g$est, GMR_lo = g$lo, GMR_hi = g$hi,
      width_mean_pp = 100 * mean(width, na.rm = TRUE), n_R_mean = mean(n_R), n_T_mean = mean(n_T))
  }, by = .(scenario, schedule, endpoint)]
  w <- dcast(b, scenario + schedule + trial ~ endpoint, value.var = c("pass", "GMR"))
  for (ep in c("Cmax", "AUClast", "AUCinf_all", "AUCinf_reliable", "AUCinf_true", "AUCinf_subC")) {
    if (!paste0("pass_", ep) %in% names(w)) w[, paste0("pass_", ep) := NA]; if (!paste0("GMR_", ep) %in% names(w)) w[, paste0("GMR_", ep) := NA_real_]
  }
  prop <- function(v) { v <- v[!is.na(v)]; wilson_ci(sum(v), length(v)) }
  cc <- w[, {
    a <- prop(pass_AUClast == pass_AUCinf_reliable); d1 <- prop(pass_AUClast & !pass_AUCinf_reliable); d2 <- prop(!pass_AUClast & pass_AUCinf_reliable)
    j2 <- prop(pass_AUClast & pass_Cmax); j3 <- prop(pass_AUClast & pass_Cmax & pass_AUCinf_reliable); j3a <- prop(pass_AUClast & pass_Cmax & pass_AUCinf_all)
    ok <- is.finite(GMR_AUClast) & is.finite(GMR_AUCinf_reliable)
    .(n_trials = .N, agree = a$est, agree_lo = a$lo, agree_hi = a$hi, last_pass_inf_fail = d1$est, last_pass_inf_fail_lo = d1$lo, last_pass_inf_fail_hi = d1$hi,
      last_fail_inf_pass = d2$est, joint_last_cmax = j2$est, joint_last_cmax_lo = j2$lo, joint_last_cmax_hi = j2$hi,
      joint_3rel = j3$est, joint_3rel_lo = j3$lo, joint_3rel_hi = j3$hi, joint_3all = j3a$est, joint_3all_lo = j3a$lo, joint_3all_hi = j3a$hi,
      cor_logGMR_last_infrel = if (sum(ok) > 2) cor(log(GMR_AUClast[ok]), log(GMR_AUCinf_reliable[ok])) else NA_real_)
  }, by = .(scenario, schedule)]
  list(per_endpoint = pe, concordance = cc, wide = w)
}
