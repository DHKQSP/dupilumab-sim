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

# 데이터셋 표: gate_role = gate 이면 pass/fail, external 이면 관측 대비 비·z만(판정 NA) (D-023)
gate_row <- function(ds, nca, tol_pct, cv_range) {
  auc_mean <- mean(nca$AUClast, na.rm = TRUE); cv <- log_cv_pct(nca$AUClast)
  role <- if (is.null(ds$gate_role)) "gate" else ds$gate_role
  n_obs <- if (is.null(ds$n_subj)) NA_real_ else as.numeric(ds$n_subj)
  se_obs <- ds$auclast_sd / sqrt(n_obs)
  data.table(id = ds$id, gate_role = role, dose_mg = ds$dose_mg, presentation = if (is.null(ds$presentation)) NA_character_ else ds$presentation,
             n_obs = n_obs, wt_mean = ds$weight_mean,
             AUClast_obs_mean = ds$auclast_mean, AUClast_sim_mean = auc_mean, AUClast_ratio = auc_mean / ds$auclast_mean,
             AUClast_z = (auc_mean - ds$auclast_mean) / se_obs, AUClast_se_obs_pct = 100 * se_obs / ds$auclast_mean,
             AUClast_sim_geo = geo_mean(nca$AUClast), AUClast_obs_geo = if (is.null(ds$auclast_geo)) NA_real_ else ds$auclast_geo,
             AUClast_sim_logcv = cv, AUClast_obs_cv_arith = 100 * ds$auclast_sd / ds$auclast_mean,
             Cmax_obs_mean = ds$cmax_mean, Cmax_sim_mean = mean(nca$Cmax, na.rm = TRUE), Cmax_ratio = mean(nca$Cmax, na.rm = TRUE) / ds$cmax_mean,
             Cmax_sim_logcv = log_cv_pct(nca$Cmax),
             tmax_obs_median = if (is.null(ds$tmax_median)) NA_real_ else ds$tmax_median, tmax_sim_median = median(nca$tmax, na.rm = TRUE),
             AUCinf_obs_mean = if (is.null(ds$aucinf_mean)) NA_real_ else ds$aucinf_mean, AUCinf_sim_mean_reliable = mean(nca[reliable == TRUE, AUCinf], na.rm = TRUE),
             tlast_sim_median = median(nca$tlast_planned, na.rm = TRUE),
             pass_mean = if (role == "gate") abs(auc_mean / ds$auclast_mean - 1) <= tol_pct / 100 else NA,
             pass_cv = if (role == "gate") (cv >= cv_range[1] & cv <= cv_range[2]) else NA)
}

ds_weight_sd <- function(ds) if (is.list(ds$weight_sd)) as.numeric(ds$weight_sd$value) else as.numeric(ds$weight_sd)
ds_weight_bounds <- function(ds) if (!is.null(ds$weight_bounds)) list(male = as.numeric(ds$weight_bounds), female = as.numeric(ds$weight_bounds)) else list(male = c(50, 100), female = c(40, 90))

# 연구별 일정 조회(검토 의견 통합본 §0-3). key → 투여 후 경과일 벡터
schedule_days <- function(dz, key) { s <- dz$dataset_schedules[[key]]; if (is.null(s)) stop("알 수 없는 채혈 일정: ", key); as.numeric(unlist(if (is.list(s) && !is.null(s$days)) s$days else s)) }

# 데이터셋 모의: components가 있으면 구성 연구별 일정으로 나눠 모의하고 합친다(동일 가중, 확정 전 가정)
simulate_ds <- function(p, ds, dz, nsim, tag, sr) {
  b <- ds_weight_bounds(ds)
  if (!is.null(ds$components)) {
    k <- length(ds$components); n_k <- ceiling(nsim / k)
    parts <- lapply(seq_len(k), function(i) {
      cp <- ds$components[[i]]
      sim <- simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), schedule_days(dz, cp$schedule), n_k, paste0(tag, "_c", i), sr, b$male, b$female)
      sim$nca[, `:=`(id = id + (i - 1L) * n_k, component = cp$study)]
    })
    return(rbindlist(parts, fill = TRUE))
  }
  simulate_dataset(p, ds$dose_mg, ds$weight_mean, ds_weight_sd(ds), schedule_days(dz, ds$schedule), nsim, tag, sr, b$male, b$female)$nca
}

# 문헌 정합성용 커버리지 비(§8): 평균비·기하평균비, 비구획(산출 가능 전체·신뢰군)과 참값
coverage_ratios <- function(nca) {
  a <- nca[lambda_ok == TRUE]; r <- nca[reliable == TRUE]
  data.table(n = nrow(nca), n_aucinf = nrow(a), n_reliable = nrow(r),
             mean_ratio_nca_all = mean(a$AUClast) / mean(a$AUCinf), geo_ratio_nca_all = geo_mean(a$AUClast) / geo_mean(a$AUCinf),
             mean_ratio_nca_reliable = mean(r$AUClast) / mean(r$AUCinf), geo_ratio_nca_reliable = geo_mean(r$AUClast) / geo_mean(r$AUCinf),
             mean_ratio_true = mean(nca$AUClast_true, na.rm = TRUE) / mean(nca$AUCinf_true), geo_ratio_true = geo_mean(nca$AUClast_true) / geo_mean(nca$AUCinf_true),
             mean_ratio_obs_last_vs_true_inf = mean(nca$AUClast, na.rm = TRUE) / mean(nca$AUCinf_true))
}
