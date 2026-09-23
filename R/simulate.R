# simulate.R — 관측치 생성 (지시서 §2–§4)
# 합집합 격자에서 한 번 풀고, 일정별로 부분집합을 취한다(공통 난수: 같은 피험자·같은 잔차가 모든 일정에 쓰임).
# 잔차 y = C·(1+eps_p) + eps_a, LLOQ 미만은 BLQ(conc NA). 투여 전(0)은 항상 BLQ.

# eps 표 미리 추출: data.table(id, planned, eps_p, eps_a)
draw_eps <- function(ids, planned, sigma) {
  g <- CJ(id = ids, planned = planned)
  n <- nrow(g)
  g[, `:=`(eps_p = rnorm(n, 0, sigma[["prop"]]), eps_a = rnorm(n, 0, sigma[["add"]]))][]
}

# ipar: individual_params(); obs: data.table(id, planned, time); eps: draw_eps(); 반환 obs + C, auc, y_raw, blq, conc
simulate_observations <- function(ipar, obs, dose_mg, lloq, eps, model_id = "k2016_2cmt_linMM_ka", t_inf_day = TRUTH_T_INF_DAY) {
  stopifnot(all(c("id", "planned", "time") %in% names(obs)))
  obs2 <- rbind(obs[, .(id, time)], data.table(id = ipar$id, time = t_inf_day))
  obs2 <- unique(obs2)
  sol <- solve_model(ipar, obs2, dose_mg, model_id = model_id)
  inf <- sol[time == t_inf_day, .(id, AUCinf_true = auc)]
  sol <- sol[time != t_inf_day, .(id, time, C, auc)]
  out <- merge(obs, sol, by = c("id", "time"), all.x = TRUE, sort = TRUE)
  out <- merge(out, eps, by = c("id", "planned"), all.x = TRUE, sort = TRUE)
  out[, y_raw := C * (1 + eps_p) + eps_a]
  out[, blq := planned == 0 | y_raw < lloq]
  out[, conc := fifelse(blq, NA_real_, y_raw)]
  out[, c("eps_p", "eps_a") := NULL]
  setorder(out, id, time)
  list(obs = out[], truth = inf)
}

# 일정별 부분집합(투여 전 0 포함)
subset_schedule <- function(obs, sched_days) obs[planned %in% c(0, sched_days)]

# 관측 tlast(계획·실제) — NCA 없이
observed_tlast <- function(obs) {
  q <- obs[!blq & planned > 0]
  out <- q[, .(tlast_actual = max(time), tlast_planned = max(planned), n_quant = .N), by = id]
  miss <- setdiff(unique(obs$id), out$id)
  if (length(miss)) out <- rbind(out, data.table(id = miss, tlast_actual = NA_real_, tlast_planned = NA_real_, n_quant = 0L))
  setorder(out, id)[]
}

# NCA 결과에 참값 부착: AUClast_true = 관측 tlast(실제 시각) 시점의 auc, AUCinf_true
attach_truth <- function(nca, obs, truth) {
  at <- merge(nca[, .(id, tlast)], obs[, .(id, time, auc)], by.x = c("id", "tlast"), by.y = c("id", "time"), all.x = TRUE)
  setnames(at, "auc", "AUClast_true")
  subj <- merge(at[, .(id, AUClast_true)], truth, by = "id", all.x = TRUE)
  add_truth_metrics(nca, subj)
}
