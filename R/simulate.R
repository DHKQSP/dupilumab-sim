# simulate.R — 관측치 생성: 진농도 → 잔차 오차(가산+비례) → LLOQ 검열 (SPEC §3.5, §5.1)
# 반환 data.table(id, planned, time, C, auc, y_raw, blq, conc)
#   C     : 진농도(모델), auc: 참값 AUC(0-t) (진단 전용)
#   y_raw : 오차 포함 관측값(검열 전), blq: LLOQ 미만 플래그, conc: 검열 후 값(BLQ는 NA)
# BLQ의 NCA 취급(0 대체, 제외 등)은 BEmaster 규칙을 따른다 — 여기서는 값만 NA로 두고 플래그를 보존.

simulate_observations <- function(ipar, obs, dose_mg, sigma, lloq) {
  stopifnot(all(c("id", "time") %in% names(obs)))
  if (!"planned" %in% names(obs)) obs <- copy(obs)[, planned := time]
  sol <- solve_model(ipar, obs[, .(id, time)], dose_mg)
  # 동일 id 내 시각 중복 방지 가정(sampling.R의 min_gap). 계획시각 부착.
  sol <- merge(sol, obs[, .(id, time, planned)], by = c("id", "time"), all.x = TRUE, sort = TRUE)
  n <- nrow(sol)
  eps_p <- rnorm(n, 0, sigma[["prop"]])
  eps_a <- rnorm(n, 0, sigma[["add"]])
  sol[, y_raw := C * (1 + eps_p) + eps_a]
  sol[, blq := y_raw < lloq]
  sol[, conc := fifelse(blq, NA_real_, y_raw)]
  setcolorder(sol, c("id", "planned", "time", "C", "auc", "y_raw", "blq", "conc"))
  sol[]
}

# 관측 데이터에서 마지막 정량 시점(투여 후 일) — 채혈 계획시각 기준과 실제시각 기준 모두 반환
observed_tlast <- function(sim) {
  q <- sim[!is.na(conc) & time > 0]
  out <- q[, .(tlast_actual = max(time), tlast_planned = max(planned), n_quant = .N), by = id]
  # 정량치가 하나도 없는 개체
  miss <- setdiff(unique(sim$id), out$id)
  if (length(miss)) out <- rbind(out, data.table(id = miss, tlast_actual = NA_real_, tlast_planned = NA_real_, n_quant = 0L))
  setorder(out, id)[]
}
