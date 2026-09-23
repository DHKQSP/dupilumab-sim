# summarize.R — 단계 4 산출 지표 (SPEC §7)

# 개인 수준: 외삽 비율 분포, λz 산출 불가 비율 (NCA 표 필요)
summarize_individual <- function(nca) {
  if (is.null(nca) || nrow(nca) == 0) return(NULL)
  nca[, .(
    n = .N,
    lambda_fail_pct = 100 * mean(!lambda_ok, na.rm = TRUE),
    extrap_median = median(pct_extrap, na.rm = TRUE),
    extrap_p90 = quantile(pct_extrap, 0.9, na.rm = TRUE, names = FALSE),
    extrap_max = max(pct_extrap, na.rm = TRUE),
    extrap_gt20_pct = 100 * mean(pct_extrap > 20, na.rm = TRUE),
    bias_AUClast_vs_true_pct = 100 * (exp(mean(log(AUClast / AUClast_true), na.rm = TRUE)) - 1),
    bias_AUCinf_vs_true_pct  = 100 * (exp(mean(log(AUCinf / AUCinf_true), na.rm = TRUE)) - 1)
  ), by = .(arm)]
}

# 시험 수준: 판정 일치율, log GMR 차이 분포, 공동 1차 조합별 검정력
summarize_trials <- function(be, endpoint_sets) {
  if (is.null(be) || nrow(be) == 0) return(NULL)
  w <- dcast(be, outer + inner ~ endpoint, value.var = c("GMR", "pass"))
  agree <- w[, .(
    n_trials = .N,
    pass_AUClast = mean(pass_AUClast), pass_AUCinf = mean(pass_AUCinf), pass_Cmax = mean(pass_Cmax),
    agreement_AUClast_AUCinf = mean(pass_AUClast == pass_AUCinf),
    discord_last_pass_inf_fail = mean(pass_AUClast & !pass_AUCinf),
    discord_last_fail_inf_pass = mean(!pass_AUClast & pass_AUCinf),
    dlogGMR_mean = mean(log(GMR_AUClast) - log(GMR_AUCinf)),
    dlogGMR_sd = sd(log(GMR_AUClast) - log(GMR_AUCinf)),
    dlogGMR_p05 = quantile(log(GMR_AUClast) - log(GMR_AUCinf), 0.05, names = FALSE),
    dlogGMR_p95 = quantile(log(GMR_AUClast) - log(GMR_AUCinf), 0.95, names = FALSE)
  )]
  power <- rbindlist(lapply(endpoint_sets, function(s) {
    cols <- paste0("pass_", s$endpoints)
    data.table(set = s$name, power = mean(Reduce(`&`, lapply(cols, function(cc) w[[cc]]))))
  }))
  list(agreement = agree, power = power, wide = w)
}

# 참값 기반(진단): AUC(0-tlast)/AUC(0-inf) 참값 분포, tlast 분포
summarize_truth <- function(subjects) {
  subjects[, .(
    n = .N,
    ratio_true_median = median(ratio_true, na.rm = TRUE),
    ratio_true_p05 = quantile(ratio_true, 0.05, na.rm = TRUE, names = FALSE),
    ratio_true_p95 = quantile(ratio_true, 0.95, na.rm = TRUE, names = FALSE),
    tlast_median = median(tlast_planned, na.rm = TRUE),
    tlast_min = min(tlast_planned, na.rm = TRUE), tlast_max = max(tlast_planned, na.rm = TRUE)
  ), by = .(arm)]
}
