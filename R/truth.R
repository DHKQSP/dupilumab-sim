# truth.R — 모델 진적분 참값(진단 전용). AUC(0-t)는 auc 상태, AUC(0-inf)는 t_inf_day까지 적분 + 잔여량 검사.
TRUTH_T_INF_DAY <- 400

true_auc <- function(ipar, dose_mg, t_grid = NULL, t_inf_day = TRUTH_T_INF_DAY, model_id = "k2016_2cmt_linMM_ka") {
  ids <- ipar$id
  if (is.null(t_grid)) t_grid <- c(0, t_inf_day)
  obs <- CJ(id = ids, time = sort(unique(c(t_grid, t_inf_day))))
  sol <- solve_model(ipar, obs, dose_mg, model_id = model_id)
  inf <- sol[time == t_inf_day, .(id, AUCinf_true = auc, resid_amt = central + periph + depot)]
  if (length(dose_mg) == 1) dose_mg <- rep(dose_mg, length(ids))
  chk <- inf$resid_amt / dose_mg
  if (any(chk > 1e-6)) warning("AUCinf_true: t_inf_day=", t_inf_day, " 에서 잔여량 비율 최대 ", signif(max(chk), 3))
  list(profile = sol[time != t_inf_day], inf = inf[, .(id, AUCinf_true)])
}

true_auc_to_tlast <- function(sol, lloq) {
  q <- sol[time > 0 & C >= lloq]
  if (nrow(q) == 0) return(data.table(id = unique(sol$id), tlast_true = NA_real_, AUClast_true = NA_real_))
  q[, .SD[which.max(time), .(tlast_true = time, AUClast_true = auc)], by = id]
}
