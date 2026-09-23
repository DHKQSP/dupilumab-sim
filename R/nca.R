# nca.R — 비구획분석 (지시서 2026-09-23 §3, config/nca_rules.yaml). BEmaster 바인딩 없음(D-010).
# 입력 sim: data.table(id, time, conc) — conc는 BLQ이면 NA. time 0(투여 전)은 있어도 되고 없어도 됨.
# 벡터화: data.table 그룹 연산으로 수만 명을 한 번에 처리한다.
#
# 표준 모드 규칙
#  - Cmax/tmax: 정량 가능 농도의 최대와 첫 도달 시각. tlast/Clast: 마지막 정량 시점과 농도.
#  - AUClast: linear-up/log-down. tmax 이전 BLQ = 0 (투여 전 (0,0) 포함), tmax–tlast 사이 BLQ = 결측(인접점 연결), tlast 이후 제외.
#  - λz: tmax 이후 정량 가능 점(Cmax 제외) 중 마지막 k점(k ≥ 3) 후보, adjusted R² 최대(동률 ±1e-4는 점 수 많은 쪽), 기울기 < 0.
#  - AUCinf = AUClast + Clast/λz, 외삽률 = (AUCinf − AUClast)/AUCinf·100, 신뢰 = adj R² ≥ 0.80 & 외삽 ≤ 20%.
# BEmaster 호환 모드(교차검증 전용): 선형 사다리꼴, BLQ = 0, 마지막 양수점까지 AUCt, λz는 마지막 3점 고정 OLS.

NCA_RULES <- NULL
get_nca_rules <- function() { if (is.null(NCA_RULES)) NCA_RULES <<- read_cfg("nca_rules.yaml"); NCA_RULES }

NCA_OUT_COLS <- c("id", "n_quant", "Cmax", "tmax", "tlast", "Clast", "AUClast", "lambda_z", "adj_r2", "n_lambda",
                  "lambda_span", "t_half", "AUCinf", "pct_extrap", "lambda_ok", "reliable")

# --- 사다리꼴 세그먼트 (벡터) ---
.seg_auc <- function(t1, c1, t2, c2, method = c("linear_up_log_down", "linear")) {
  method <- match.arg(method)
  dt <- t2 - t1
  lin <- (c1 + c2) / 2 * dt
  if (method == "linear") return(lin)
  logdown <- c2 < c1 & c1 > 0 & c2 > 0
  out <- lin
  out[logdown] <- (c1[logdown] - c2[logdown]) / log(c1[logdown] / c2[logdown]) * dt[logdown]
  out
}

# --- λz 후보 창 탐색 (벡터: id별 후보점을 시간 내림차순으로 두고 누적합) ---
.lambda_z_search <- function(cand, min_points = 3, tie_tol = 1e-4) {
  # cand: data.table(id, time, conc) — tmax 이후, Cmax 제외, 정량 가능. 반환: id별 lambda_z, adj_r2, n_lambda, lambda_span
  if (nrow(cand) == 0) return(data.table(id = integer(0), lambda_z = numeric(0), adj_r2 = numeric(0), n_lambda = integer(0), lambda_span = numeric(0)))
  setorder(cand, id, -time)
  cand[, `:=`(x = time, y = log(conc))]
  cand[, `:=`(k = seq_len(.N), Sx = cumsum(x), Sy = cumsum(y), Sxx = cumsum(x * x), Syy = cumsum(y * y), Sxy = cumsum(x * y),
              tmin = cummin(time), tmax_w = time[1]), by = id]
  w <- cand[k >= min_points]
  if (nrow(w) == 0) return(data.table(id = integer(0), lambda_z = numeric(0), adj_r2 = numeric(0), n_lambda = integer(0), lambda_span = numeric(0)))
  w[, `:=`(sxx = k * Sxx - Sx^2, syy = k * Syy - Sy^2, sxy = k * Sxy - Sx * Sy)]
  w[, slope := fifelse(sxx > 0, sxy / sxx, NA_real_)]
  w[, r2 := fifelse(sxx > 0 & syy > 0, sxy^2 / (sxx * syy), NA_real_)]
  w[, adj_r2 := 1 - (1 - r2) * (k - 1) / (k - 2)]
  w <- w[is.finite(slope) & slope < 0 & is.finite(adj_r2)]
  if (nrow(w) == 0) return(data.table(id = integer(0), lambda_z = numeric(0), adj_r2 = numeric(0), n_lambda = integer(0), lambda_span = numeric(0)))
  # 선택: adj R² 최대, 동률(±tie_tol)이면 k 큰 쪽
  w[, best := max(adj_r2), by = id]
  w <- w[adj_r2 >= best - tie_tol]
  w <- w[order(id, -k)][, .SD[1], by = id]
  w[, .(id, lambda_z = -slope, adj_r2, n_lambda = as.integer(k), lambda_span = tmax_w - tmin)]
}

run_nca <- function(sim, mode = c("standard", "bemaster_compat"), rules = get_nca_rules()) {
  mode <- match.arg(mode)
  stopifnot(all(c("id", "time", "conc") %in% names(sim)))
  d <- sim[, .(id, time, conc)]
  d <- d[time >= 0]
  ids <- sort(unique(d$id))
  if (mode == "bemaster_compat") d[is.na(conc), conc := 0]
  # 투여 전 (0, NA/0) 앵커 보장
  has0 <- d[time == 0, unique(id)]
  if (length(setdiff(ids, has0))) d <- rbind(d, data.table(id = setdiff(ids, has0), time = 0, conc = if (mode == "bemaster_compat") 0 else NA_real_))
  setorder(d, id, time)
  d[, quant := !is.na(conc) & conc > 0]

  # 요약: Cmax/tmax/tlast/Clast
  s <- d[quant == TRUE, {
    i <- which.max(conc); j <- which.max(time)
    .(n_quant = .N, Cmax = conc[i], tmax = time[i], tlast = time[j], Clast = conc[j])
  }, by = id]
  none <- setdiff(ids, s$id)
  # --- AUClast ---
  dd <- merge(d, s[, .(id, tmax, tlast)], by = "id")
  if (mode == "standard") {
    dd <- dd[time <= tlast]
    dd[time < tmax & !quant, conc := 0]                 # tmax 이전 BLQ = 0 (투여 전 앵커 포함)
    dd <- dd[quant | time < tmax]                         # 중간 BLQ = 결측(행 제거 → 인접점 연결)
    method <- rules$standard$auc_method
  } else {
    dd <- dd[time <= tlast]                               # 마지막 양수점까지, BLQ = 0
    method <- rules$bemaster_compat$auc_method
  }
  setorder(dd, id, time)
  dd[, `:=`(t_prev = shift(time), c_prev = shift(conc)), by = id]
  seg <- dd[!is.na(t_prev), .(AUClast = sum(.seg_auc(t_prev, c_prev, time, conc, method))), by = id]
  s <- merge(s, seg, by = "id", all.x = TRUE)
  s[is.na(AUClast), AUClast := 0]

  # --- λz ---
  if (mode == "standard") {
    lr <- rules$standard$lambda_z
    cand <- merge(d[quant == TRUE, .(id, time, conc)], s[, .(id, tmax)], by = "id")[time > tmax]
    lz <- .lambda_z_search(cand[, .(id, time, conc)], min_points = lr$min_points, tie_tol = lr$tie_tolerance)
  } else {
    cand <- merge(d[quant == TRUE, .(id, time, conc)], s[, .(id, tmax)], by = "id")
    setorder(cand, id, -time)
    cand <- cand[, head(.SD, 3), by = id][, if (.N == 3) .SD, by = id]
    lz <- cand[, {
      x <- time; y <- log(conc); sl <- if (var(x) > 0) cov(x, y) / var(x) else NA_real_
      .(lambda_z = -sl, adj_r2 = NA_real_, n_lambda = 3L, lambda_span = max(x) - min(x))
    }, by = id][is.finite(lambda_z) & lambda_z > 0]
  }
  s <- merge(s, lz, by = "id", all.x = TRUE)
  s[, lambda_ok := !is.na(lambda_z)]
  s[, t_half := log(2) / lambda_z]
  s[, AUCinf := AUClast + Clast / lambda_z]
  s[, pct_extrap := (AUCinf - AUClast) / AUCinf * 100]
  rel <- rules$standard$reliability
  s[, reliable := lambda_ok & !is.na(adj_r2) & adj_r2 >= rel$adj_r2_min & pct_extrap <= rel$extrap_max_pct]
  if (mode == "bemaster_compat") s[, reliable := NA]
  if (length(none)) s <- rbind(s, data.table(id = none, n_quant = 0L, Cmax = NA_real_, tmax = NA_real_, tlast = NA_real_, Clast = NA_real_,
                                             AUClast = NA_real_, lambda_z = NA_real_, adj_r2 = NA_real_, n_lambda = NA_integer_, lambda_span = NA_real_,
                                             t_half = NA_real_, AUCinf = NA_real_, pct_extrap = NA_real_, lambda_ok = FALSE, reliable = FALSE), fill = TRUE)
  setorder(s, id)
  setcolorder(s, NCA_OUT_COLS)
  s[]
}

# 참값 기반 외삽 지표: subj에 AUClast_true(관측 tlast까지 진적분), AUCinf_true 가 있어야 한다.
add_truth_metrics <- function(nca, subj) {
  m <- merge(nca, subj[, .(id, AUClast_true, AUCinf_true)], by = "id", all.x = TRUE)
  m[, coverage_true := AUClast_true / AUCinf_true]
  m[, pct_extrap_true := (1 - coverage_true) * 100]
  m[, err_AUClast_vs_true_tlast := log(AUClast / AUClast_true)]
  m[, err_AUClast_vs_true_inf := log(AUClast / AUCinf_true)]
  m[]
}
