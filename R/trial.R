# trial.R — 가상 시험 엔진 (지시서 §4–§6)
# 시험 j마다: (1) 2n명 피험자 추출(체중·성별·IIV·ADA) (2) 층화 1:1 배정 (3) 합집합 격자 채혈시각·잔차 추출
# (4) 대조군 1회 풀이, 시험군은 시나리오별 배율로 풀이(같은 eta·같은 잔차 = 공통 난수, D-013)
# (5) 일정 × 시나리오별 부분집합 → NCA → 참값 부착 → BE.
# 시드 파생: (master, j, "subj"), (master, j, "alloc"), (master, j, arm, "jitter"), (master, j, arm, "eps") — 시나리오·일정과 무관.

simulate_trial_arms <- function(j, p, design, scenarios, schedules, master_seed, wt_spec, jitter = TRUE,
                                n_per_arm = NULL, sex_ratio_male = NULL, model_id = NULL) {
  if (is.null(n_per_arm)) n_per_arm <- design$n_per_arm
  if (is.null(sex_ratio_male)) sex_ratio_male <- design$weight$sex_ratio_male$value
  if (is.null(model_id)) model_id <- p$model_id
  dose <- design$dose_mg; lloq <- p$lloq
  grid <- union_grid(design, schedules)
  # (1) 피험자 (2) 배정
  subj <- with_seed(derive_seed(master_seed, j, "subj"), make_subjects(2 * n_per_arm, p, wt_spec, sex_ratio_male, p$ada$fraction))
  st <- strata_from_weight_spec(wt_spec, as.numeric(design$stratification$split_kg$value))
  subj <- with_seed(derive_seed(master_seed, j, "alloc"), assign_arms_stratified(subj, st$breaks, st$labels))
  arms <- list()
  for (a in c("R", "T")) {
    e <- subj[arm == a]
    ids <- e$id
    obs <- with_seed(derive_seed(master_seed, j, a, "jitter"), make_obs_times(ids, grid, design, jitter = jitter))
    eps <- with_seed(derive_seed(master_seed, j, a, "eps"), draw_eps(ids, sort(unique(c(0, grid))), p$sigma))
    if (a == "R") {
      ip <- individual_params(p, e)
      sim <- simulate_observations(ip, obs, dose, lloq, eps, model_id = model_id)
      arms[[a]] <- list(S00 = list(obs = sim$obs, truth = sim$truth, subj = e, ipar = ip))
    } else {
      # 시나리오 스택: id 오프셋으로 한 번에 풀이
      res <- list(); ip_all <- list(); obs_all <- list(); eps_all <- list(); off <- 0L; offs <- list()
      for (sc in names(scenarios)) {
        ps <- apply_multipliers(p, scenarios[[sc]]$T_multipliers)
        ip <- individual_params(ps, e); ip[, id := id + off]
        o2 <- copy(obs)[, id := id + off]; e2 <- copy(eps)[, id := id + off]
        ip_all[[sc]] <- ip; obs_all[[sc]] <- o2; eps_all[[sc]] <- e2; offs[[sc]] <- off
        off <- off + max(ids)
      }
      sim <- simulate_observations(rbindlist(ip_all), rbindlist(obs_all), dose, lloq, rbindlist(eps_all), model_id = model_id)
      for (sc in names(scenarios)) {
        o <- offs[[sc]]
        ob <- sim$obs[id > o & id <= o + max(ids)][, id := id - o]
        tr <- sim$truth[id > o & id <= o + max(ids)][, id := id - o]
        arms[[a]][[sc]] <- list(obs = ob, truth = tr, subj = e, ipar = copy(ip_all[[sc]])[, id := id - o])
      }
    }
  }
  list(j = j, subj = subj, arms = arms, grid = grid)
}

# 한 (시나리오, 일정)의 NCA + BE. 반환 list(be, nca)
analyze_arms <- function(sim_arms, scenario, sched_days, design, nca_mode = "standard", methods = "pooled_t") {
  parts <- list()
  for (a in c("R", "T")) {
    x <- sim_arms$arms[[a]][[if (a == "R") "S00" else scenario]]
    ob <- subset_schedule(x$obs, sched_days)
    nca <- run_nca(ob, mode = nca_mode)
    nca <- attach_truth(nca, ob, x$truth)
    nca <- merge(nca, x$subj[, .(id, WT, sex, stratum)], by = "id")
    nca[, AUCinf_subC := fifelse(reliable %in% TRUE, AUCinf, AUClast)]
    nca[, arm := a]
    parts[[a]] <- nca
  }
  nca_all <- rbindlist(parts)
  be <- be_analyze(nca_all, BE_ENDPOINTS, ci_level = design$be$ci_level, limits = as.numeric(design$be$limits), methods = methods)
  list(be = be, nca = nca_all)
}

# 시험 j 전체: 시나리오 × 일정 조합의 BE 표(+ 개체 지표 요약). combos: data.table(scenario, schedule)
run_trial <- function(j, p, design, scenarios, combos, master_seed, wt_spec, jitter = TRUE, methods = "pooled_t", model_id = NULL,
                      keep_nca = FALSE, n_per_arm = NULL) {
  scen_needed <- scenarios[unique(combos$scenario)]
  sched_needed <- unique(combos$schedule)
  sa <- simulate_trial_arms(j, p, design, scen_needed, sched_needed, master_seed, wt_spec, jitter = jitter, model_id = model_id, n_per_arm = n_per_arm)
  be_list <- list(); ind_list <- list(); nca_list <- list()
  for (r in seq_len(nrow(combos))) {
    sc <- combos$scenario[r]; sh <- combos$schedule[r]
    an <- analyze_arms(sa, sc, get_schedule(design, sh), design, methods = methods)
    be_list[[r]] <- an$be[, `:=`(trial = j, scenario = sc, schedule = sh)]
    ind_list[[r]] <- merge(summarize_individual(an$nca, by = "arm"), summarize_dropout(an$nca, by = "arm"), by = "arm")[, `:=`(trial = j, scenario = sc, schedule = sh)]
    if (keep_nca) nca_list[[r]] <- an$nca[, `:=`(trial = j, scenario = sc, schedule = sh)]
  }
  list(be = rbindlist(be_list), ind = rbindlist(ind_list), nca = if (keep_nca) rbindlist(nca_list) else NULL)
}
