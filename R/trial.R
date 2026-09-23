# trial.R — 가상 시험 1회 (SPEC §5.1–5.4)
# 1) 두 arm 피험자 생성(같은 IIV, T arm만 theta 차이) 2) 채혈 편차 3) 관측치(오차·BLQ)
# 4) 참값 AUC(진단) 5) NCA(BEmaster) 6) BE 판정(BEmaster)
# nca_fn/be_fn 은 테스트에서 skip 용도로만 교체하며, 결과 보고에는 BEmaster 어댑터만 쓴다.

run_trial <- function(pR, pT, design, seed, schedule_name = NULL, n_per_arm = NULL,
                      lloq = NULL, dose_mg = NULL, nca_fn = run_nca_bemaster, be_fn = run_be_bemaster,
                      endpoints = c("AUClast", "AUCinf", "Cmax"), keep_profiles = FALSE, do_be = TRUE) {
  if (is.null(n_per_arm)) n_per_arm <- design$n_per_arm$value
  if (is.null(lloq)) lloq <- design$lloq_mg_L$value
  if (is.null(dose_mg)) dose_mg <- design$dose_mg
  planned <- get_schedule(design, schedule_name)

  sims <- with_seed(seed, {
    sR <- make_subjects(n_per_arm, pR, design, id_offset = 0L)
    sT <- make_subjects(n_per_arm, pT, design, id_offset = n_per_arm)
    oR <- make_obs_times(sR$ipar$id, planned, design, jitter = TRUE)
    oT <- make_obs_times(sT$ipar$id, planned, design, jitter = TRUE)
    simR <- simulate_observations(sR$ipar, oR, dose_mg, pR$sigma, lloq)
    simT <- simulate_observations(sT$ipar, oT, dose_mg, pT$sigma, lloq)
    list(sR = sR, sT = sT, simR = simR, simT = simT)
  })
  ipar <- rbind(sims$sR$ipar[, arm := "R"], sims$sT$ipar[, arm := "T"])
  sim  <- rbind(sims$simR[, arm := "R"], sims$simT[, arm := "T"])
  ipar[, WT_stratum := weight_stratum(WT, design)]

  # 참값(진단): AUCinf_true, AUC(0-tlast_obs)_true
  tr <- true_auc(ipar, dose_mg, t_grid = c(0, max(planned)))
  tl <- observed_tlast(sim)
  # tlast_actual 시각의 참값 AUC = 그 시각까지의 auc 상태 (관측 시각에 정확히 존재)
  auc_at_tlast <- sim[!is.na(conc) & time > 0, .SD[which.max(time), .(tlast = time, AUClast_true = auc)], by = id]
  subj <- merge(ipar[, .(id, arm, WT, sex, WT_stratum, F_i = F, Vc_i, ke, Vmax, Km, ka)], tr$inf, by = "id", all.x = TRUE)
  subj <- merge(subj, auc_at_tlast, by = "id", all.x = TRUE)
  subj <- merge(subj, tl, by = "id", all.x = TRUE)
  subj[, ratio_true := AUClast_true / AUCinf_true]

  out <- list(seed = seed, n_per_arm = n_per_arm, schedule = if (is.null(schedule_name)) design$default_schedule else schedule_name,
              subjects = subj, sim = if (keep_profiles) sim else NULL, nca = NULL, be = NULL)
  if (!is.null(nca_fn)) {
    nca <- nca_fn(sim, dose_mg, lloq)
    nca <- merge(nca, subj[, .(id, arm, WT, WT_stratum, AUCinf_true, AUClast_true, ratio_true)], by = "id")
    out$nca <- nca
    if (do_be && !is.null(be_fn)) out$be <- be_fn(nca, endpoints = endpoints,
                                                  ci_level = design$ci_level,
                                                  limits = c(design$be_limits$lower, design$be_limits$upper))
  }
  out
}
