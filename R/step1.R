# step1.R — 단계 1 검증 도우미 (지시서 §1). 데이터셋별 체중 매칭 모집단을 만들어 Clot 일정으로 관측치·NCA를 산출한다.
# 채혈 편차 없음(문헌 재현), 잔차·BLQ·IIV 포함.

# wt: list(mean, sd, bounds_male = c(50,100), bounds_female = c(40,90)) — 지시서 §1 절단 규칙
simulate_dataset <- function(p, dose_mg, wt_mean, wt_sd, sched_days, n, seed_tag, sex_ratio_male = 0.5,
                             bounds_male = c(50, 100), bounds_female = c(40, 90), design_for_windows = NULL) {
  wt_spec <- list(mean = wt_mean, sd = wt_sd, trunc = NULL, male_bounds = bounds_male, female_bounds = bounds_female)
  subj <- with_seed(derive_seed(seed_tag, "subj"), make_subjects(n, p, wt_spec, sex_ratio_male, 0))
  planned <- sort(unique(c(0, sched_days)))
  obs <- CJ(id = subj$id, planned = planned)[, time := planned]
  eps <- with_seed(derive_seed(seed_tag, "eps"), draw_eps(subj$id, planned, p$sigma))
  ip <- individual_params(p, subj)
  sim <- simulate_observations(ip, obs, dose_mg, p$lloq, eps, model_id = p$model_id)
  nca <- attach_truth(run_nca(sim$obs), sim$obs, sim$truth)
  nca <- merge(nca, subj[, .(id, WT, sex)], by = "id")
  tl <- observed_tlast(sim$obs)
  list(subj = subj, obs = sim$obs, nca = merge(nca, tl[, .(id, tlast_planned)], by = "id"))
}

# 코호트(8명) 중앙값 tlast 분포
cohort_median_tlast <- function(nca, n_per_cohort = 8) {
  x <- copy(nca)[, cohort := (id - 1L) %/% n_per_cohort]
  x[, .(tlast_median = median(tlast_planned, na.rm = TRUE), tlast_min = min(tlast_planned, na.rm = TRUE), tlast_max = max(tlast_planned, na.rm = TRUE)), by = cohort]
}

# 데이터셋 gate 표
gate_row <- function(ds, nca, tol_pct, cv_range) {
  auc_mean <- mean(nca$AUClast, na.rm = TRUE); cv <- log_cv_pct(nca$AUClast)
  data.table(id = ds$id, dose_mg = ds$dose_mg, n_obs = if (is.null(ds$n_subj)) NA_integer_ else as.integer(ds$n_subj), wt_mean = ds$weight_mean,
             AUClast_obs_mean = ds$auclast_mean, AUClast_sim_mean = auc_mean, AUClast_ratio = auc_mean / ds$auclast_mean,
             AUClast_sim_geo = geo_mean(nca$AUClast), AUClast_obs_geo = if (is.null(ds$auclast_geo)) NA_real_ else ds$auclast_geo,
             AUClast_sim_logcv = cv, AUClast_obs_cv_arith = 100 * ds$auclast_sd / ds$auclast_mean,
             Cmax_obs_mean = ds$cmax_mean, Cmax_sim_mean = mean(nca$Cmax, na.rm = TRUE), Cmax_ratio = mean(nca$Cmax, na.rm = TRUE) / ds$cmax_mean,
             Cmax_sim_logcv = log_cv_pct(nca$Cmax),
             AUCinf_obs_mean = if (is.null(ds$aucinf_mean)) NA_real_ else ds$aucinf_mean, AUCinf_sim_mean_reliable = mean(nca[reliable == TRUE, AUCinf], na.rm = TRUE),
             tmax_sim_median = median(nca$tmax, na.rm = TRUE), tlast_sim_median = median(nca$tlast_planned, na.rm = TRUE),
             pass_mean = abs(auc_mean / ds$auclast_mean - 1) <= tol_pct / 100,
             pass_cv = cv >= cv_range[1] & cv <= cv_range[2])
}

ds_weight_sd <- function(ds) if (is.list(ds$weight_sd)) as.numeric(ds$weight_sd$value) else as.numeric(ds$weight_sd)
ds_weight_bounds <- function(ds) if (!is.null(ds$weight_bounds)) list(male = as.numeric(ds$weight_bounds), female = as.numeric(ds$weight_bounds)) else list(male = c(50, 100), female = c(40, 90))
