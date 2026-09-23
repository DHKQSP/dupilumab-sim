# mc.R — 시험 반복 루프 (단일층 MC, 시험당 n_trials; D-012). 병렬은 fork(parallel::mclapply).
run_trials <- function(n_trials, p, design, scenarios, combos, master_seed, wt_spec, jitter = TRUE, methods = "pooled_t",
                       model_id = NULL, cores = 1L, progress_every = 50) {
  invisible(get_model(if (is.null(model_id)) p$model_id else model_id))   # fork 전에 컴파일(자식 프로세스 간 경합 방지)
  one <- function(j) {
    r <- run_trial(j, p, design, scenarios, combos, master_seed, wt_spec, jitter = jitter, methods = methods, model_id = model_id)
    if (progress_every > 0 && j %% progress_every == 0) cat(sprintf("  trial %d/%d %s\n", j, n_trials, format(Sys.time(), "%H:%M:%S")))
    r
  }
  res <- if (cores > 1) parallel::mclapply(seq_len(n_trials), one, mc.cores = cores, mc.preschedule = TRUE) else lapply(seq_len(n_trials), one)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
  if (any(bad)) stop("실패한 시험 반복: ", paste(which(bad), collapse = ","))
  list(be = rbindlist(lapply(res, `[[`, "be")), ind = rbindlist(lapply(res, `[[`, "ind")))
}

# 개인 수준 대규모 모집단(지시서 §5): n명, 합집합 격자 1회 풀이 → 일정별 NCA
run_individual_population <- function(n, p, design, schedules, master_seed, wt_spec, jitter = TRUE, sex_ratio_male = NULL,
                                      model_id = NULL, dose_mg = NULL, multipliers = list(), tag = "pop") {
  if (is.null(sex_ratio_male)) sex_ratio_male <- design$weight$sex_ratio_male$value
  if (is.null(model_id)) model_id <- p$model_id
  if (is.null(dose_mg)) dose_mg <- design$dose_mg
  grid <- union_grid(design, schedules)
  subj <- with_seed(derive_seed(master_seed, tag, "subj"), make_subjects(n, p, wt_spec, sex_ratio_male, p$ada$fraction))
  obs <- with_seed(derive_seed(master_seed, tag, "jitter"), make_obs_times(subj$id, grid, design, jitter = jitter))
  eps <- with_seed(derive_seed(master_seed, tag, "eps"), draw_eps(subj$id, sort(unique(c(0, grid))), p$sigma))
  ip <- individual_params(apply_multipliers(p, multipliers), subj)
  sim <- simulate_observations(ip, obs, dose_mg, p$lloq, eps, model_id = model_id)
  nca_by_sched <- lapply(schedules, function(sh) {
    ob <- subset_schedule(sim$obs, get_schedule(design, sh))
    nca <- attach_truth(run_nca(ob), ob, sim$truth)
    nca <- merge(nca, subj[, .(id, WT, sex, ada)], by = "id")
    nca[, schedule := sh][]
  })
  names(nca_by_sched) <- schedules
  list(subj = subj, obs = sim$obs, truth = sim$truth, nca = rbindlist(nca_by_sched))
}
