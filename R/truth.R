# truth.R — 모델 기반 참값(진단 전용). NCA 대체가 아니다 (SPEC §6, DECISIONS D-004).
# AUC(0-t)는 ODE의 auc 상태, AUC(0-inf)는 긴 시간 지평(t_inf_day)까지의 적분 + 잔여량 검사.

TRUTH_T_INF_DAY <- 400

true_auc <- function(ipar, dose_mg, t_grid = NULL, t_inf_day = TRUTH_T_INF_DAY) {
  ids <- ipar$id
  if (is.null(t_grid)) t_grid <- c(0, t_inf_day)
  obs <- CJ(id = ids, time = sort(unique(c(t_grid, t_inf_day))))
  sol <- solve_model(ipar, obs, dose_mg)
  inf <- sol[time == t_inf_day, .(id, AUCinf_true = auc, resid_amt = central + periph + depot)]
  # 잔여량이 투여량의 1e-6 이상이면 지평이 짧다
  if (length(dose_mg) == 1) dose_mg <- rep(dose_mg, length(ids))
  chk <- inf$resid_amt / dose_mg
  if (any(chk > 1e-6)) warning("AUCinf_true: t_inf_day=", t_inf_day, " 에서 잔여량 비율 최대 ", signif(max(chk), 3))
  list(profile = sol[time != t_inf_day], inf = inf[, .(id, AUCinf_true)])
}

# 관측 프로파일(진농도 C, 시각)에서 정량 가능한 마지막 시각과 그 시각까지의 참값 AUC
true_auc_to_tlast <- function(sol, lloq) {
  q <- sol[time > 0 & C >= lloq]
  if (nrow(q) == 0) return(data.table(id = unique(sol$id), tlast_true = NA_real_, AUClast_true = NA_real_))
  q[, .SD[which.max(time), .(tlast_true = time, AUClast_true = auc)], by = id]
}
