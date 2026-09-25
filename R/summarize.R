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
    flag_rsq_pct = if ("flag_rsq" %in% names(.SD)) 100 * mean(flag_rsq %in% TRUE) else NA_real_, flag_extrap_pct = if ("flag_extrap" %in% names(.SD)) 100 * mean(flag_extrap %in% TRUE) else NA_real_,
    flag_span_pct = if ("flag_span" %in% names(.SD)) 100 * mean(flag_span %in% TRUE) else NA_real_,
    reliable_no_span_pct = if ("flag_span" %in% names(.SD)) 100 * mean(lambda_ok & !(flag_rsq %in% TRUE) & !(flag_extrap %in% TRUE)) else NA_real_,
    n_lambda_median = as.numeric(median(n_lambda, na.rm = TRUE)), n_lambda_eq3_pct = 100 * mean(n_lambda == 3, na.rm = TRUE),
    extrap_median = median(pct_extrap, na.rm = TRUE), extrap_p95 = q95(pct_extrap), extrap_gt20_pct = 100 * mean(pct_extrap > 20, na.rm = TRUE),
    extrap_true_median = median(pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(pct_extrap_true), coverage_lt80_pct = 100 * mean(coverage_true < 0.80, na.rm = TRUE),
    err_tlast_median = median(err_AUClast_vs_true_tlast, na.rm = TRUE), err_tlast_sd = sd(err_AUClast_vs_true_tlast, na.rm = TRUE),
    err_inf_sd = sd(err_AUClast_vs_true_inf, na.rm = TRUE), err_inf_median = median(err_AUClast_vs_true_inf, na.rm = TRUE),
    AUClast_geo = geo_mean(AUClast), AUClast_logcv = log_cv_pct(AUClast), AUClast_mean = mean(AUClast, na.rm = TRUE),
    Cmax_geo = geo_mean(Cmax), Cmax_logcv = log_cv_pct(Cmax), Cmax_mean = mean(Cmax, na.rm = TRUE),
    AUCinf_geo = geo_mean(AUCinf[reliable]), AUCinf_true_geo = geo_mean(AUCinf_true)
  ), by = by, .SDcols = intersect(c("flag_rsq", "flag_extrap", "flag_span"), names(nca))]
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

# ---------------------------------------------------------------------------------------------
# 보고 문구용 AUCinf 신뢰 충족률·탈락률과 그 증감(검토 의견 W2 §1). 보고서(report.Rmd)·summary_en(18)·key_numbers_en(37)이 같은 값과 같은 전제 검사를 쓴다.
# 플래그 세트 (i) = λz 산출 & adj R² ≥ 0.80 & 외삽 ≤ 20%, (ii) = (i) & span ratio ≥ 2(일부 통계분석계획만 쓰는 관행, Phoenix 기능 아님, D-039).
# 문구는 (i)을 주값, (ii)를 병기한다. 입력: results/reliability/(scripts/39), rationale/pillar1_two_model_range.csv,
# trials/schedule_decision_<변형>.csv, nca_engine/engine_difference_individual_B0.csv. 문구의 정성 주장(아래 premise)이 결과와 어긋나면 stop().
# 39 산출물이 없으면 NULL(호출 쪽은 신뢰 충족률 문구를 쓰지 않는다).
reliability_facts <- function(root = proj_path("results")) {
  rd_ <- function(...) { f <- file.path(root, ...); if (file.exists(f)) fread(f) else NULL }
  two <- rd_("reliability", "reliability_two_flag_sets_summary.csv"); prd <- rd_("reliability", "reliability_paired_vs_B0.csv")
  drs <- rd_("reliability", "dropout_reasons_by_schedule.csv"); lzw <- rd_("reliability", "reliability_lz_window_by_schedule.csv")
  if (is.null(two) || is.null(prd) || is.null(drs)) return(NULL)
  premise <- function(ok, msg) if (!isTRUE(all(ok))) stop("신뢰 충족률 문구의 전제가 결과와 다르다: ", msg, call. = FALSE)
  M2 <- c(base = "k2016", struct2020 = "k2020"); DENSE <- c("D1", "D2", "D3", "D4")
  # 1. B0 두 모델: (i)·(ii) 충족률(Wilson 95% 구간), 미충족률, 탈락 사유(중복 포함, (ii))
  t0 <- two[match(names(M2), two$variant)]; d0 <- drs[drs$schedule == "B0"]; d0 <- d0[match(names(M2), d0$variant)]
  premise(!anyNA(t0$variant) && !anyNA(d0$variant), "두 모델(base, struct2020)의 B0 행이 없다")
  b0 <- data.table(variant = names(M2), model = unname(M2), rel_i = t0$reliable_i_pct, rel_i_lo = t0$reliable_i_lo, rel_i_hi = t0$reliable_i_hi,
                   rel_ii = t0$reliable_ii_pct, rel_ii_lo = t0$reliable_ii_lo, rel_ii_hi = t0$reliable_ii_hi, span_loss = t0$span_only_loss_pct,
                   fail_i = 100 - t0$reliable_i_pct, fail_ii = 100 - t0$reliable_ii_pct, lambda_fail = d0$lambda_fail_pct, flag_rsq = d0$flag_rsq_pct,
                   flag_extrap = d0$flag_extrap_pct, flag_span = d0$flag_span_pct, any_ii = d0$any_flag_or_fail_pct, excl_span_only = d0$excl_span_pct)
  premise(abs(b0$fail_ii - b0$any_ii) < 0.01, "(ii) 미충족률 = 탈락 사유 합집합")
  premise(abs(b0$span_loss - b0$excl_span_only) < 0.01, "span 단독 손실 = (i) − (ii) = span 플래그만 있는 대상자")
  premise(b0$rel_i > b0$rel_ii, "B0에서 (i) 충족률 > (ii) 충족률")
  p1r <- rd_("rationale", "pillar1_two_model_range.csv")
  if (!is.null(p1r)) premise(abs(min(b0$rel_ii) - p1r$reliable_min) < 0.01 && abs(max(b0$rel_ii) - p1r$reliable_max) < 0.01,
                             "(ii) 두 모델 범위 = rationale/pillar1_two_model_range.csv")
  rng <- list(rel_i = range(b0$rel_i), rel_ii = range(b0$rel_ii), fail_i = range(b0$fail_i), fail_ii = range(b0$fail_ii), span_loss = range(b0$span_loss),
              gt20 = if (!is.null(p1r)) c(p1r$gt20_min, p1r$gt20_max) else NULL)
  # 2. B0 대비 쌍대 증감(같은 20,000명, %p). 두 모델 D1–D4에서 span 플래그가 증감을 낮춘다((ii) < (i))는 문구의 전제
  premise(!(prd$crit_c_i %in% TRUE) & !(prd$crit_c_ii %in% TRUE), "기준 (c)는 어느 변형·일정에서도 두 세트 모두 미충족")
  premise(prd[!is.na(prd$recommend_ii), recommend_i == recommend_ii], "두 플래그 세트의 권고 판정이 같다")
  for (v_ in unique(prd$variant)) { f <- file.path(root, "trials", sprintf("schedule_decision_%s.csv", v_)); if (!file.exists(f)) next
    m <- merge(prd[prd$variant == v_, .(schedule, gain_ii_pp)], fread(f)[, .(schedule, c_reliable_gain_pp)], by = "schedule")
    premise(abs(m$gain_ii_pp - m$c_reliable_gain_pp) < 0.01, sprintf("(ii) 증감 = trials/schedule_decision_%s.csv의 c_reliable_gain_pp", v_)) }
  pr2 <- prd[prd$variant %in% names(M2)][, model := M2[variant]][]
  dn <- pr2[schedule %in% DENSE]
  premise(nrow(dn) == 2 * length(DENSE), "두 모델 × D1–D4 쌍대 증감 행")
  premise(dn$gain_ii_pp < dn$gain_i_pp, "두 모델 D1–D4 모두 (ii) 증감 < (i) 증감(span 단독 손실 증가)")
  flips <- dn[sign(round(gain_i_pp, 2)) != sign(round(gain_ii_pp, 2)), .(variant, model, schedule, gain_i_pp, gain_ii_pp)]
  o3 <- dn[schedule != "D3"]
  dense <- list(rise = range(dn$span_only_loss_pct - dn$span_only_loss_B0_pct), other_i = range(o3$gain_i_pp), other_ii = range(o3$gain_ii_pp),
                d3 = dn[schedule == "D3"], bminus = pr2[schedule == "Bminus"])
  # 3. 이전 엔진(D-010) 대비 하락이 거의 전부 span 플래그라는 문구의 전제(같은 관측치, 엔진 비교용 B0 단독 표본)
  eng <- rd_("nca_engine", "engine_difference_individual_B0.csv")
  if (!is.null(eng)) {
    eo <- eng[grepl("D-010", eng$engine)]; en <- eng[grepl("D-039", eng$engine)]
    eng <- merge(eo[, .(variant = model, old = reliable_pct)], en[, .(variant = model, new_i = reliable_rsq_extrap_only_pct, new_ii = reliable_pct)], by = "variant")
    eng[, span_share_pct := 100 * (new_i - new_ii) / (old - new_ii)]
    eng <- eng[match(names(M2), eng$variant)][, model := M2[variant]][]
    premise(eng$old > eng$new_ii & eng$span_share_pct >= 90, "이전 엔진 대비 (ii) 하락의 90% 이상이 span 플래그 단독")
  }
  # 4. 잔차 가정에 민감하다는 문구의 전제: 잔차 비례 12% 변형과 2016 모델의 Wilson 구간이 두 세트 모두 겹치지 않고, 외삽 20% 초과 비율이 0.8배 미만 또는 1.25배 초과
  resid <- NULL; r12 <- two[two$variant == "resid12"]; d12 <- rd_("trials", "schedule_decision_resid12.csv"); db <- rd_("trials", "schedule_decision_base.csv")
  if (nrow(r12) && !is.null(d12) && !is.null(db)) {
    bb <- two[two$variant == "base"]
    resid <- data.table(rel_i = r12$reliable_i_pct, rel_ii = r12$reliable_ii_pct, gt20 = d12$d_extrap20_B0[1], gt20_base = db$d_extrap20_B0[1])
    premise(c(r12$reliable_i_lo > bb$reliable_i_hi || r12$reliable_i_hi < bb$reliable_i_lo, r12$reliable_ii_lo > bb$reliable_ii_hi || r12$reliable_ii_hi < bb$reliable_ii_lo),
            "잔차 비례 12% 변형의 신뢰 충족률 구간이 2016 모델과 겹치지 않는다((i)·(ii))")
    premise(resid$gt20 / resid$gt20_base < 0.8 || resid$gt20 / resid$gt20_base > 1.25, "잔차 비례 12% 변형의 외삽 20% 초과 비율이 2016 모델의 0.8–1.25배 밖")
  }
  # 5. λz 창 구조: Day 22 이후 명목일로 만들 수 있는 가장 짧은 3점 창(일)
  lz_min <- if (!is.null(lzw)) c(B0 = unique(lzw[lzw$schedule == "B0", min_3pt_nominal_window_days]), dense = unique(lzw[lzw$schedule %in% DENSE, min_3pt_nominal_window_days])) else NULL
  if (!is.null(lz_min)) premise(length(lz_min) == 2 && lz_min[["dense"]] < lz_min[["B0"]], "추가 채혈 일정의 최단 3점 창이 B0보다 짧다")
  list(b0 = b0, rng = rng, paired = pr2, dense = dense, flips = flips, eng = eng, resid = resid, lz_min = lz_min, all = prd)
}

# 반올림: 20,000명 비율은 0.005의 배수라 CSV에서 읽은 값(예: -1.015 → -1.01499…)이 아래로 잘린다. 0에서 먼 쪽으로 반올림하고 음의 0을 없앤다
round_half_away <- function(x, d = 1) round(x + sign(x) * 1e-9, d) + 0
fmt_num <- function(x, d = 1, signed = FALSE) { x <- round_half_away(x, d); if (signed) sprintf(paste0("%+.", d, "f"), x) else formatC(x, format = "f", digits = d) }
# "(i) x [(ii) y]" 형식. each = FALSE: i·ii 각각을 값 하나 또는 범위(최솟값–최댓값)로, each = TRUE: 원소마다 한 쌍(벡터 반환). sep: 범위 구분, signed: 부호 표시(증감)
fmt_flag_pair <- function(i, ii, d = 1, unit = "%", sep = "–", signed = FALSE, each = FALSE) {
  f <- function(x) fmt_num(x, d, signed)
  if (each) return(sprintf("(i) %s%s [(ii) %s%s]", f(i), unit, f(ii), unit))
  one <- function(x) if (length(x) == 1 || isTRUE(all.equal(min(x), max(x)))) paste0(f(x[1]), unit) else paste0(f(min(x)), sep, f(max(x)), unit)
  sprintf("(i) %s [(ii) %s]", one(i), one(ii))
}

# ---------------------------------------------------------------------------------------------
# 보고 우선순위와 핵심 문장(지시 2026-09-25 §4)의 수치. 보고서(요약표, Pillar 3 첫 문장, 역산 경계표)·summary_en(18)·key_numbers_en(37)이
# 같은 값과 같은 전제 검사를 쓴다. 문구의 정성 주장(아래 premise)이 결과와 어긋나면 stop(). 입력 파일이 없으면 그 원소는 NULL.
#   km      : 결합 상수(Km) 역산. oc/inversion_all.csv(없으면 oc/inversion_k20*_*.csv)의 탐색 범위 끝·도달 행(200,000명)과
#             oc/inversion_scan_<모델>_Km.csv(스크리닝 20,000명)에서 참 AUC0-inf 비의 범위, 도달 가능한 목표와 그 배율
#   inv_bnd : 모델 × 기전 × 경계 목표(0.80, 1.25): 도달 배율과 그 배율의 참 Cmax 비, 도달 불가면 목표 쪽 탐색 범위 끝 배율과 참값 비
#   vm150, ke120: 예비 제품 시나리오(trials5000/products5000_props_base.csv, 주 모델). AUCinf 신뢰군은 플래그 세트 (ii).
#             참값 비 = 같은 시험 대상자의 개인 모델 AUC0-inf 시험 GMR의 평균(fallback/consumer_risk.csv·discordance_classification.csv의
#             true_ratio = rationale/pillar2_products_B0.csv GMR_mean_AUCinf_true). 200,000명 적분 참값(역산)과는 다른 정의
#   p2max, g2max: 경계 1종 오류(oc/boundary_type1.csv) 최대 칸과 적응적 연장 판정(oc/extension_decision.csv가 있으면)
key_facts <- function(root = proj_path("results"), oc = read_cfg("oc_design.yaml")) {
  rd_ <- function(...) { f <- file.path(root, ...); if (file.exists(f)) fread(f) else NULL }
  premise <- function(ok, msg) if (!isTRUE(all(ok))) stop("핵심 수치 문구의 전제가 결과와 다르다: ", msg, call. = FALSE)
  bnd <- as.numeric(unlist(oc$boundary_targets)); MODELS <- c("k2016", "k2020")
  out <- list(km = NULL, inv_bnd = NULL, vm150 = NULL, ke120 = NULL, p2max = NULL, g2max = NULL, bnd = bnd)
  # 1. 역산(200,000명 공통 난수)
  iv <- rd_("oc", "inversion_all.csv")
  if (is.null(iv)) { fs <- list.files(file.path(root, "oc"), pattern = "^inversion_k20(16|20)_.*\\.csv$", full.names = TRUE); if (length(fs)) iv <- rbindlist(lapply(fs, fread), fill = TRUE) }
  if (!is.null(iv) && nrow(iv)) {
    ib <- iv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9]
    out$inv_bnd <- rbindlist(lapply(split(ib, by = c("model", "mechanism", "target")), function(g) {
      lab <- sprintf("%s %s %.2f", g$model[1], g$mechanism[1], g$target[1]); tg <- g$target[1]
      r <- g[reachable %in% TRUE]
      premise(nrow(r) <= 1, paste0(lab, ": 경계 목표에 닿는 방향은 많아야 하나"))
      if (nrow(r)) return(r[, .(model, mechanism, target, reachable = TRUE, direction, multiplier, auc_ratio, cmax_ratio, end_multiplier = NA_real_, end_auc_ratio = NA_real_, end_cmax_ratio = NA_real_)])
      # 도달 불가: 목표 쪽(참값 비가 1에서 목표 방향으로 움직이는) 탐색 범위 끝 중 목표에 가장 가까운 것. 그 비가 목표에 못 미쳐야 한다
      u <- g[is.finite(end_auc_ratio) & sign(log(end_auc_ratio)) == sign(log(tg))]
      premise(nrow(u) >= 1, paste0(lab, ": 도달 불가인데 목표 쪽 탐색 범위 끝이 없다"))
      u <- u[which.max(abs(log(end_auc_ratio)))]
      premise(abs(log(u$end_auc_ratio)) < abs(log(tg)), paste0(lab, ": 도달 불가 방향의 범위 끝 참값 비가 목표에 못 미친다"))
      u[, .(model, mechanism, target, reachable = FALSE, direction, multiplier = NA_real_, auc_ratio = NA_real_, cmax_ratio = NA_real_, end_multiplier, end_auc_ratio, end_cmax_ratio)] }))
    kr <- iv[mechanism == "Km"]
    if (nrow(kr)) {
      rng_cfg <- as.numeric(unlist(oc$mechanisms$Km$range))
      premise(setequal(unique(kr$model), MODELS), "Km 역산이 두 모델 모두에 있다")
      ends <- unique(kr[is.finite(end_auc_ratio), .(model, direction, end_multiplier, end_auc_ratio, end_cmax_ratio)])
      premise(nrow(ends) == 2 * length(MODELS) && !anyDuplicated(ends, by = c("model", "direction")), "Km 탐색 범위 끝(모델 × 방향)의 참값 비가 하나씩")
      premise(abs(ends[direction == "down", end_multiplier] - rng_cfg[1]) < 1e-12 & abs(ends[direction == "up", end_multiplier] - rng_cfg[2]) < 1e-12,
              sprintf("Km 탐색 범위 끝 = config 범위(×%s–×%s)", format(rng_cfg[1]), format(rng_cfg[2])))
      reach <- kr[reachable %in% TRUE]
      tg <- unique(reach$target)
      premise(length(tg) == 1 && nrow(reach) == length(MODELS) && setequal(reach$model, MODELS), "Km으로 도달 가능한 사전 고정 목표는 하나이고 두 모델 모두 도달")
      premise(reach$within_tol %in% TRUE, "Km 도달 행이 허용 오차(±0.1%) 안")
      sc <- rbindlist(lapply(MODELS, function(m_) rd_("oc", sprintf("inversion_scan_%s_Km.csv", m_))))
      if (nrow(sc)) premise(sc[order(multiplier), .(ok = all(diff(auc_ratio_screen) >= 0) || all(diff(auc_ratio_screen) <= 0)), by = model]$ok,
                            "스크리닝 스캔에서 참 AUC0-inf 비가 Km 배율에 단조(범위 = 양 끝)")
      vals <- c(ends$end_auc_ratio, reach$auc_ratio, if (nrow(sc)) sc$auc_ratio_screen)
      rr <- range(vals)
      premise(rr[1] > bnd[1] && rr[2] < bnd[2], "Km 탐색 범위 전체에서 참 AUC0-inf 비가 동등성 한계 안")
      tl <- as.numeric(unlist(oc$targets)); tl <- tl[abs(tl - 1) > 1e-9]
      premise(all(abs(unique(kr$target) - 1) > 1e-9) && setequal(round(unique(kr$target), 6), round(tl, 6)), "Km 역산 목표 = config targets(1.00 제외)")
      premise(uniqueN(kr$n_subjects) == 1, "Km 역산(범위 끝·도달 행)의 참값 대상자 수가 하나(문구 'N명')")
      out$km <- list(range_mult = rng_cfg, ratio = rr, target = tg, direction = unique(reach$direction),
                     mult = setNames(reach$multiplier, reach$model)[MODELS], cmax = setNames(reach$cmax_ratio, reach$model)[MODELS],
                     ends = ends, targets = range(tl), n_targets = length(tl), n_subjects = unique(kr$n_subjects),
                     n_screen = if (nrow(sc)) as.integer(oc$inversion$screening_subjects) else NA_integer_)
    }
  }
  # 2. 예비 제품 시나리오(주 모델, 5,000회 이상)
  pr5 <- rd_("trials5000", "products5000_props_base.csv"); cr <- rd_("fallback", "consumer_risk.csv"); dc <- rd_("fallback", "discordance_classification.csv")
  p2p <- rd_("rationale", "pillar2_products_B0.csv")
  g5 <- function(s_, m_) { x <- pr5[pr5$scenario == s_ & pr5$metric == m_]; premise(nrow(x) == 1, sprintf("products5000_props_base.csv에 %s %s 행이 하나", s_, m_)); x }
  if (!is.null(pr5) && !is.null(cr)) {
    ri <- g5("VM150", "AUCinf_reliable"); rl <- g5("VM150", "AUClast"); tr <- cr[cr$scenario == "VM150"]
    premise(nrow(tr) == 1 && tr$n_trials == ri$n_trials, "VM150 참값 비(consumer_risk.csv)가 같은 시험 수에서")
    if (!is.null(p2p)) { pp <- p2p[p2p$model == "k2016" & p2p$scenario == "VM150" & p2p$n_trials == ri$n_trials]
      premise(nrow(pp) == 1 && abs(pp$GMR_mean_AUCinf_true - tr$true_ratio) < 1e-9 && abs(pp$pass_rate_AUCinf_reliable - ri$est) < 1e-9 && abs(pp$pass_rate_AUClast - rl$est) < 1e-9,
              "VM150: consumer_risk.csv true_ratio = pillar2_products_B0.csv GMR_mean_AUCinf_true, 통과율 = products5000_props_base.csv") }
    premise(ri$est > 5, "VM150 AUCinf(신뢰군) 단독 통과율 > 명목 5%")
    premise(tr$true_ratio < bnd[1], "VM150 참 AUC0-inf 비 < 0.80(범위 밖)")
    premise(rl$hi < 5, "VM150 AUClast 단독 통과율의 Wilson 상한 < 5%")
    premise(ri$n_trials >= 5000L && rl$n_trials == ri$n_trials, "VM150 시험 수 ≥ 5,000(문구 '5,000회 이상'), AUCinf·AUClast 같은 시험 수")
    out$vm150 <- list(n_trials = ri$n_trials, rel = ri[, .(est, lo, hi)], last = rl[, .(est, lo, hi)], true_ratio = tr$true_ratio, rel_lo_gt5 = ri$lo > 5, flag_set = "(ii)")
  }
  if (!is.null(pr5) && !is.null(dc)) {
    d_ <- g5("KE120", "discord_last_pass_inf_fail"); c_ <- dc[dc$scenario == "KE120"]
    premise(nrow(c_) == 1 && c_$n_trials == d_$n_trials && abs(c_$last_pass_inf_fail - d_$est) < 1e-9, "KE120 불일치율: discordance_classification.csv = products5000_props_base.csv")
    premise(c_$true_ratio >= bnd[1] && c_$true_ratio <= bnd[2] && isTRUE(c_$true_inside), "KE120 참 AUC0-inf 비가 80–125% 안")
    premise(startsWith(c_$classification, "AUCinf 위음성"), "KE120 분류 = AUCinf 위음성(discordance_classification.csv)")
    premise(d_$n_trials >= 5000L, "KE120 시험 수 ≥ 5,000(문구 '5,000회 이상')")
    out$ke120 <- list(n_trials = d_$n_trials, est = d_$est, lo = d_$lo, hi = d_$hi, true_ratio = c_$true_ratio, flag_set = "(ii)")
  }
  # 3. 경계 1종 오류 최대 칸(1차 지표)과 적응적 연장
  bt <- rd_("oc", "boundary_type1.csv"); dec <- rd_("oc", "extension_decision.csv")
  if (!is.null(bt) && nrow(bt)) {
    cls <- function(r) classify_type1(r$lo, r$hi, 5)                  # Wilson 95% 구간 대 5% (R/oc_interpret.R, scripts/43과 같은 분류)
    one <- function(cf) { x <- bt[bt$config == cf]; if (!nrow(x)) return(NULL); r <- x[which.max(x$pass_pct)]
      e <- if (!is.null(dec)) dec[dec$model == r$model & dec$scenario == r$scenario & dec$selected %in% TRUE] else NULL
      if (!is.null(e) && !nrow(e)) e <- NULL
      if (!is.null(e)) {
        premise(r$n_trials >= e$n_before, sprintf("%s %s: 시험 수 ≥ 연장 전 시험 수", r$model, r$scenario))
        # 문구("미실행", "진행 중", 완료)는 boundary_type1.csv의 시험 수에서 고른다: 연장 파일의 시험이 모두 반영되어 있어야 한다(아니면 scripts/33을 다시 실행)
        ef <- file.path(root, "oc", sprintf("oc_trials_ext_be_%s.csv.gz", r$model))
        n_ext <- if (file.exists(ef)) { z <- fread(ef, select = c("trial", "scenario")); uniqueN(z$trial[z$scenario == r$scenario]) } else 0L
        premise(r$n_trials == e$n_before + n_ext, sprintf("%s %s: boundary_type1.csv 시험 수 %d = 사전 고정 %d + 연장 파일 %d회(scripts/33이 연장 파일을 반영)",
                                                          r$model, r$scenario, as.integer(r$n_trials), as.integer(e$n_before), as.integer(n_ext)))
      }
      list(row = r, class = cls(r), ext = e, n_gt5 = sum(x$pass_pct > 5), n_cells = nrow(x)) }
    out$p2max <- one("P2"); out$g2max <- one("G2")
  }
  out
}

# 무작위 제품 공간의 가정 분포(config/oc_design.yaml random_space.ranges, 기전별 배율 로그 균등) 문구. 그림 3-C 부제·결론 문단(scripts/33),
# 보고서·summary_en의 보조 지표 절이 같은 문구를 쓴다. 예: "F 0.80–1.25, ka 0.67–1.50, …, Km 0.20–5"
random_space_ranges_text <- function(oc = read_cfg("oc_design.yaml"), sep = "–") {
  r <- oc$random_space$ranges; stopifnot(length(r) > 0, all(lengths(r) == 2))
  fr_ <- function(v) sub("\\.00$", "", formatC(as.numeric(v), format = "f", digits = 2))
  paste(vapply(names(r), function(k) sprintf("%s %s%s%s", k, fr_(r[[k]][[1]]), sep, fr_(r[[k]][[2]])), ""), collapse = ", ")
}

# 시험 수 문구: 적응적 연장(extension_decision.csv) 대상이면 사전 고정 수·그때 값과 연장 상태를 붙인다
key_trials_text <- function(r, e, lang = c("ko", "en")) {
  lang <- match.arg(lang); fi <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
  if (is.null(e)) return(if (lang == "ko") sprintf("시험 %s회", fi(r$n_trials)) else sprintf("%s trials", fi(r$n_trials)))
  pre_ko <- sprintf("사전 고정 %s회 %s%% [%s, %s]", fi(e$n_before), fmt_num(e$pass_pct, 2), fmt_num(e$lo, 2), fmt_num(e$hi, 2))
  pre_en <- sprintf("pre-registered %s trials: %s%% (95%% CI %s to %s)", fi(e$n_before), fmt_num(e$pass_pct, 2), fmt_num(e$lo, 2), fmt_num(e$hi, 2))
  if (r$n_trials <= e$n_before) return(if (lang == "ko") sprintf("시험 %s회(적응적 연장 대상, 목표 %s회, 미실행)", fi(r$n_trials), fi(e$n_after))
                                       else sprintf("%s trials (selected for adaptive extension to %s trials, not yet run)", fi(r$n_trials), fi(e$n_after)))
  prog <- r$n_trials < e$n_after
  if (lang == "ko") sprintf("시험 %s회(적응적 연장%s; %s)", fi(r$n_trials), if (prog) sprintf(" 진행 중, 목표 %s회", fi(e$n_after)) else "", pre_ko)
  else sprintf("%s trials (adaptive extension%s; %s)", fi(r$n_trials), if (prog) sprintf(" in progress, target %s", fi(e$n_after)) else "", pre_en)
}

# Pillar 3 첫 문장(결합 상수): key_facts()$km에서
key_km_sentence <- function(km, lang = c("ko", "en")) {
  lang <- match.arg(lang); fg <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  scan_ko <- if (is.na(km$n_screen)) "" else sprintf(", 그 사이 스크리닝 스캔은 %s명", format(km$n_screen, big.mark = ","))
  scan_en <- if (is.na(km$n_screen)) "" else sprintf(", screening scan in between with %s", format(km$n_screen, big.mark = ","))
  if (lang == "ko") sprintf("결합 상수(Km)를 %s–%s배로 바꿔도 참 AUC0-inf 비는 %s–%s에 머문다(두 모델; 탐색 범위 끝과 도달 행은 %s명%s, 공통 난수). 사전 고정 목표 %d개(%s–%s) 중 도달 가능한 것은 %s 하나다(Km 약 ×%s(2016 모델), ×%s(Model 1)).",
                            fg(km$range_mult[1]), fg(km$range_mult[2]), fmt_num(km$ratio[1], 3), fmt_num(km$ratio[2], 3), format(km$n_subjects, big.mark = ","), scan_ko, km$n_targets,
                            fmt_num(km$targets[1], 2), fmt_num(km$targets[2], 2), fmt_num(km$target, 2), fmt_num(km$mult[["k2016"]], 0), fmt_num(km$mult[["k2020"]], 0))
  else sprintf("Changing the binding constant (Km) from %s to %s times keeps the true AUC0-inf ratio within %s to %s (both models; range ends and reachable rows with %s common-random-number subjects%s). Of the %d pre-specified targets (%s to %s) the only reachable one is %s (Km about x%s in the 2016 model and x%s in Model 1).",
               fg(km$range_mult[1]), fg(km$range_mult[2]), fmt_num(km$ratio[1], 3), fmt_num(km$ratio[2], 3), format(km$n_subjects, big.mark = ","), scan_en, km$n_targets,
               fmt_num(km$targets[1], 2), fmt_num(km$targets[2], 2), fmt_num(km$target, 2), fmt_num(km$mult[["k2016"]], 0), fmt_num(km$mult[["k2020"]], 0))
}

# 역산 경계표(모델 × 기전): 참 AUC0-inf 비 0.80·1.25를 주는 배율과 그 배율의 참 Cmax 비. 도달 불가는 목표 쪽 범위 끝 배율과 그 참 AUC0-inf 비
key_inv_bnd_table <- function(ib, lang = c("ko", "en"), mech = c("F", "ka", "ke", "Vmax", "Km", "V2")) {
  lang <- match.arg(lang); fg <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  bnd <- sort(unique(ib$target)); stopifnot(length(bnd) == 2)
  ML <- if (lang == "ko") c(k2016 = "2016 모델", k2020 = "Model 1") else c(k2016 = "2016 model", k2020 = "Model 1")
  mult_cell <- function(r) if (!nrow(r)) "" else if (isTRUE(r$reachable)) sprintf(if (lang == "ko") "×%s" else "x%s", fmt_num(r$multiplier, 3))
    else sprintf(if (lang == "ko") "도달 불가(범위 끝 ×%s: %s)" else "unreachable (range end x%s: %s)", fg(r$end_multiplier), fmt_num(r$end_auc_ratio, 3))
  cmax_cell <- function(r) if (nrow(r) && isTRUE(r$reachable)) fmt_num(r$cmax_ratio, 3) else ""
  keys <- unique(ib[, .(model, mechanism)])[order(match(model, names(ML)), match(mechanism, mech))]
  rows <- lapply(seq_len(nrow(keys)), function(i) { k <- keys[i]; a <- ib[model == k$model & mechanism == k$mechanism & abs(target - bnd[1]) < 1e-9]
    b <- ib[model == k$model & mechanism == k$mechanism & abs(target - bnd[2]) < 1e-9]
    data.table(v1 = ML[[k$model]], v2 = k$mechanism, v3 = mult_cell(a), v4 = cmax_cell(a), v5 = mult_cell(b), v6 = cmax_cell(b)) })
  x <- rbindlist(rows)
  setnames(x, if (lang == "ko") c("모델", "기전", sprintf("참값 %s: 배율", fmt_num(bnd[1], 2)), sprintf("참값 %s: 참 Cmax 비", fmt_num(bnd[1], 2)), sprintf("참값 %s: 배율", fmt_num(bnd[2], 2)), sprintf("참값 %s: 참 Cmax 비", fmt_num(bnd[2], 2)))
              else c("Model", "Mechanism", sprintf("True %s: multiplier", fmt_num(bnd[1], 2)), sprintf("True %s: true Cmax ratio", fmt_num(bnd[1], 2)), sprintf("True %s: multiplier", fmt_num(bnd[2], 2)), sprintf("True %s: true Cmax ratio", fmt_num(bnd[2], 2))))
  x[]
}

# 요약표(보고서 요약·summary_en 핵심 결론): 1차 지표(경계 1종 오류 최대), 예비 제품 시나리오 VM150·KE120, Pillar 3(Km). 모든 수치는 key_facts()에서
key_summary_table <- function(kf, lang = c("ko", "en")) {
  lang <- match.arg(lang); ko <- lang == "ko"
  pc <- function(e, lo, hi) { d <- if (e < 1) 3 else 2; sprintf(if (ko) "%s%% [%s, %s]" else "%s%% (95%% CI %s to %s)", fmt_num(e, d), fmt_num(lo, d), fmt_num(hi, d)) }
  ML <- if (ko) c(k2016 = "2016 모델", k2020 = "Model 1") else c(k2016 = "2016 model", k2020 = "Model 1")
  CL <- if (ko) c(conservative = "보수적(Wilson 상한 < 5%)", nominal = "명목(Wilson 구간이 5% 포함)", exceeding = "초과(Wilson 하한 > 5%)")
        else c(conservative = "conservative (Wilson upper bound below 5%)", nominal = "nominal (Wilson interval includes 5%)", exceeding = "exceeding (Wilson lower bound above 5%)")
  DL <- if (ko) c(down = "하향", up = "상향") else c(down = "down", up = "up")
  rows <- list()
  bt_row <- function(m, what) { r <- m$row
    data.table(a = what, b = sprintf("%s; %s", pc(r$pass_pct, r$lo, r$hi), key_trials_text(r, m$ext, lang)),
               c = sprintf(if (ko) "%s(%s %s %s, 목표 %s)" else "%s (%s, %s %s, target %s)", fmt_num(r$auc_ratio, 3), ML[[r$model]], r$mechanism, DL[[r$direction]], fmt_num(r$target, 2)),
               d = if (ko) sprintf("%s; 경계 칸 %d개 중 점추정 5%% 초과 %d개", CL[[m$class]], m$n_cells, m$n_gt5)
                   else sprintf("%s; point estimate above 5%% in %d of %d boundary cells", CL[[m$class]], m$n_gt5, m$n_cells),
               e = "oc/boundary_type1.csv") }
  if (!is.null(kf$p2max)) rows$p2 <- bt_row(kf$p2max, if (ko) "1차 지표: P2(AUClast + Cmax) 경계 1종 오류, 최대 칸" else "Primary metric: P2 (AUClast + Cmax) boundary type I error, largest cell")
  if (!is.null(kf$g2max)) rows$g2 <- bt_row(kf$g2max, if (ko) "1차 지표: G2(AUCinf 규칙 A + Cmax) 경계 1종 오류, 최대 칸" else "Primary metric: G2 (AUCinf rule A + Cmax) boundary type I error, largest cell")
  v <- kf$vm150; k <- kf$ke120; fi <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
  tr_src <- if (ko) "같은 시험 대상자의 개인 모델 AUC0-inf 시험 GMR 평균(fallback/consumer_risk.csv true_ratio = rationale/pillar2_products_B0.csv GMR_mean_AUCinf_true)"
            else "mean trial GMR of the individual model AUC0-inf of the same trial subjects (fallback/consumer_risk.csv true_ratio = rationale/pillar2_products_B0.csv GMR_mean_AUCinf_true)"
  if (!is.null(v)) {
    rows$vm_rel <- data.table(a = if (ko) sprintf("VM150(Vmax ×1.50) AUCinf 신뢰군 단독 통과율(플래그 세트 %s, 규칙 A)", v$flag_set) else sprintf("VM150 (Vmax x1.50): pass rate of AUCinf alone, reliable subjects (flag set %s, rule A)", v$flag_set),
                              b = sprintf(if (ko) "%s; 시험 %s회" else "%s; %s trials", pc(v$rel$est, v$rel$lo, v$rel$hi), fi(v$n_trials)),
                              c = sprintf(if (ko) "%s(범위 밖); %s" else "%s (outside the limits); %s", fmt_num(v$true_ratio, 3), tr_src),
                              d = if (ko) sprintf("참값이 한계 밖인데 명목 5%% 초과%s", if (v$rel_lo_gt5) "(Wilson 하한도 5% 초과)" else "(점추정)")
                                  else sprintf("exceeds the nominal 5%% although the truth is outside the limits%s", if (v$rel_lo_gt5) " (Wilson lower bound also above 5%)" else " (point estimate)"),
                              e = "trials5000/products5000_props_base.csv")
    rows$vm_last <- data.table(a = if (ko) "VM150 AUClast 단독 통과율" else "VM150: pass rate of AUClast alone",
                               b = sprintf(if (ko) "%s; 시험 %s회" else "%s; %s trials", pc(v$last$est, v$last$lo, v$last$hi), fi(v$n_trials)),
                               c = sprintf(if (ko) "%s(범위 밖)" else "%s (outside the limits)", fmt_num(v$true_ratio, 3)),
                               d = if (ko) "같은 범위 밖 제품을 AUClast는 5% 미만으로 통과(Wilson 상한 < 5%)" else "the same out-of-limits product passes on AUClast below 5% (Wilson upper bound below 5%)",
                               e = "trials5000/products5000_props_base.csv") }
  if (!is.null(k)) rows$ke <- data.table(a = if (ko) sprintf("KE120(ke ×1.20) AUClast 통과·AUCinf 신뢰군 불통과(플래그 세트 %s)", k$flag_set) else sprintf("KE120 (ke x1.20): AUClast pass, AUCinf (reliable subjects, flag set %s) fail", k$flag_set),
                                         b = sprintf(if (ko) "%s; 시험 %s회" else "%s; %s trials", pc(k$est, k$lo, k$hi), fi(k$n_trials)),
                                         c = sprintf(if (ko) "%s(범위 안); 같은 정의(discordance_classification.csv true_ratio)" else "%s (inside the limits); same definition (discordance_classification.csv true_ratio)", fmt_num(k$true_ratio, 3)),
                                         d = if (ko) "참값이 한계 안이므로 AUCinf의 위음성으로 분류" else "truth inside the limits, so classified as false negatives of AUCinf",
                                         e = "trials5000/products5000_props_base.csv, fallback/discordance_classification.csv")
  if (!is.null(kf$km)) { km <- kf$km; fg <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
    rows$km <- data.table(a = if (ko) "Pillar 3: 결합 상수 Km 배율의 참 AUC0-inf 비" else "Pillar 3: true AUC0-inf ratio across binding-constant (Km) multipliers",
                          b = if (ko) sprintf("Km ×%s–×%s에서 %s–%s", fg(km$range_mult[1]), fg(km$range_mult[2]), fmt_num(km$ratio[1], 3), fmt_num(km$ratio[2], 3))
                              else sprintf("%s to %s for Km x%s to x%s", fmt_num(km$ratio[1], 3), fmt_num(km$ratio[2], 3), fg(km$range_mult[1]), fg(km$range_mult[2])),
                          c = sprintf(if (ko) "모델 적분 참값(%s명 공통 난수)" else "model-integrated truth (%s common-random-number subjects)", format(km$n_subjects, big.mark = ",")),
                          d = sprintf(if (ko) "동등성 한계 안; 도달 가능한 사전 고정 목표는 %s 하나(Km 약 ×%s(2016 모델), ×%s(Model 1))" else "inside the limits; the only reachable pre-specified target is %s (Km about x%s in the 2016 model, x%s in Model 1)",
                                      fmt_num(km$target, 2), fmt_num(km$mult[["k2016"]], 0), fmt_num(km$mult[["k2020"]], 0)),
                          e = "oc/inversion_all.csv, oc/inversion_scan_k20*_Km.csv") }
  x <- rbindlist(rows)
  if (!nrow(x)) return(x)
  setnames(x, if (ko) c("항목", "추정 (95% 구간)", "참 AUC0-inf 비", "해석", "출처 (results/)") else c("Item", "Estimate (95% CI)", "True AUC0-inf ratio", "Interpretation", "Source (results/)"))
  x[]
}
