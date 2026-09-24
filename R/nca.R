# nca.R — 비구획분석. 표준 엔진은 Phoenix WinNonlin 8.x 호환(검토 의견 2026-09-24 §1, D-039). 참조 구현 NonCompart 0.8.4, 보조 PKNCA 0.12.1과
# 검증(scripts/28_nca_engine_validation.R). 입력 sim: data.table(id, time, conc) — conc는 BLQ이면 NA. time 0(투여 전)은 있어도 되고 없어도 됨.
# 벡터화: data.table 그룹 연산으로 수만 명을 한 번에 처리한다.
#
# BLQ 전처리(Phoenix는 자동 처리하지 않으므로 명시, 통계분석계획 관행): preprocess_blq()
#  - 첫 정량값 이전 BLQ = 0 (투여 전 (0, 0) 앵커 포함)
#  - 첫 정량값 이후 BLQ가 2회 연속이면 그 뒤의 정량값은 결측(= 첫 번째 연속 BLQ 시점 이후 전부 제외)
#  - 정량값 사이 BLQ = 결측(행 제거 → 인접 정량점 연결), 마지막 정량값 이후 BLQ = 제외
# NCA(run_nca, Phoenix 표기 출력):
#  - AUClast: Linear Up Log Down(농도가 감소하고 두 값이 양수인 구간만 로그 사다리꼴), 실제 채혈 시각.
#  - λz Best Fit(NonCompart::BestSlope와 같은 순서): Cmax 이후 양수 농도의 마지막 3, 4, 5, ... 점으로 비가중 OLS(ln C),
#    기울기 > 0인 창 제외, adjusted R² = 1 − (1 − R²)(n − 1)/(n − 2) 최대 창, |최대 − 창| < 1e-4이면 점 수가 많은 창.
#    Cmax 이후 양수 농도 3개 미만이거나 선택된 기울기가 음이 아니면 산출 불가.
#  - Clast_pred = exp(b0 − λz·Tlast), AUCINF_obs/pred = AUClast + Clast(_pred)/λz, AUC_%Extrap = (1 − AUClast/AUCINF)·100.
#  - 신뢰 플래그(통계분석계획 관행, Phoenix 기능 아님): Rsq_adjusted < 0.80, AUC_%Extrap_obs > 20%, span ratio < 2.
#    reliable = λz 산출 가능 & 플래그 없음. 규칙 A 플래그 대상 제외 / B 플래그 후 포함 / C 플래그(또는 산출 불가)면 AUClast 대입.
#  - 하위 호환 별칭(내부 요약 코드용): tmax, tlast, lambda_z, adj_r2, n_lambda, lambda_span, t_half, AUCinf(= AUCINF_obs), pct_extrap(= AUC_%Extrap_obs).
# 이전 엔진(run_nca_legacy): tmax 이전 BLQ = 0, 2회 연속 BLQ 규칙 없음, 동률 판정 ≤ 1e-4, span 플래그 없음 — 차이표(scripts/29) 전용.
# BEmaster 호환 모드(교차검증 전용): 선형 사다리꼴, BLQ = 0, 마지막 양수점까지 AUCt, λz는 마지막 3점 고정 OLS.

NCA_RULES <- NULL
get_nca_rules <- function() { if (is.null(NCA_RULES)) NCA_RULES <<- read_cfg("nca_rules.yaml"); NCA_RULES }

NCA_OUT_COLS_LEGACY <- c("id", "n_quant", "Cmax", "tmax", "tlast", "Clast", "AUClast", "lambda_z", "adj_r2", "n_lambda",
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
.lambda_z_search_legacy <- function(cand, min_points = 3, tie_tol = 1e-4) {
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

run_nca_legacy <- function(sim, mode = c("standard", "bemaster_compat"), rules = get_nca_rules()) {
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
    lz <- .lambda_z_search_legacy(cand[, .(id, time, conc)], min_points = lr$min_points, tie_tol = lr$tie_tolerance)
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
  setcolorder(s, NCA_OUT_COLS_LEGACY)
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


# ---------------------------------------------------------------------------------------------
# Phoenix WinNonlin 호환 표준 엔진 (D-039)
PHOENIX_COLS <- c("Cmax", "Tmax", "Tlast", "Clast", "AUClast", "Lambda_z", "HL_Lambda_z", "Rsq", "Rsq_adjusted", "No_points_lambda_z",
                  "Lambda_z_lower", "Lambda_z_upper", "Clast_pred", "AUCINF_obs", "AUCINF_pred", "AUC_%Extrap_obs", "AUC_%Extrap_pred", "Span_ratio")
NCA_OUT_COLS <- c("id", "n_quant", PHOENIX_COLS, "flag_rsq", "flag_extrap", "flag_span", "lambda_ok", "reliable",
                  "tmax", "tlast", "lambda_z", "adj_r2", "n_lambda", "lambda_span", "t_half", "AUCinf", "AUCinf_pred", "pct_extrap")

# BLQ 전처리. 입력 data.table(id, time, conc[NA = BLQ]). 반환: 분석에 쓰는 행(id, time, conc)만 — 첫 정량 이전 0, 이후 양수 정량값(≤ Tlast).
preprocess_blq <- function(sim) {
  d <- sim[time >= 0, .(id, time, conc)]
  ids <- unique(d$id)
  has0 <- d[time == 0, unique(id)]
  if (length(setdiff(ids, has0))) d <- rbind(d, data.table(id = setdiff(ids, has0), time = 0, conc = NA_real_))   # 투여 전 앵커(BLQ → 0)
  setorder(d, id, time)
  d[, q := !is.na(conc) & conc > 0]
  d[, `:=`(k = seq_len(.N), first_q = if (any(q)) which(q)[1] else NA_integer_), by = id]
  # 첫 정량 이후 처음 나타나는 연속 BLQ 2회의 시작 위치: 그 위치부터 전부 제외
  d[, cut2 := {
    b <- !q; pos <- NA_integer_
    if (!is.na(first_q[1]) && .N >= 2) { cand <- which(b[-.N] & b[-1] & seq_len(.N - 1) > first_q[1]); if (length(cand)) pos <- cand[1] }
    pos }, by = id]
  d <- d[is.na(cut2) | k < cut2]
  d[, last_q := if (any(q)) max(which(q)) else NA_integer_, by = id]
  out <- d[!is.na(first_q) & ((k < first_q) | (q & k <= last_q))]
  out[k < first_q, conc := 0]
  out[, .(id, time, conc)]
}

# 창별 OLS를 위한 누적합(중심 이동으로 상쇄 오차 억제). cand: id, time, conc(양수, 시간 내림차순으로 정렬됨)
.lambda_z_best_fit <- function(cand, min_points = 3L, tie_tol = 1e-4) {
  empty <- data.table(id = integer(0), Lambda_z = numeric(0), Rsq = numeric(0), Rsq_adjusted = numeric(0), No_points_lambda_z = integer(0),
                      Lambda_z_lower = numeric(0), Lambda_z_upper = numeric(0), b0 = numeric(0))
  if (!nrow(cand)) return(empty)
  setorder(cand, id, -time)
  cand[, `:=`(x0 = time[1], y0 = log(conc[1])), by = id]
  cand[, `:=`(x = time - x0, y = log(conc) - y0)]
  cand[, `:=`(k = seq_len(.N), Sx = cumsum(x), Sy = cumsum(y), Sxx = cumsum(x * x), Syy = cumsum(y * y), Sxy = cumsum(x * y), lo = time, up = time[1]), by = id]
  w <- cand[k >= min_points]
  if (!nrow(w)) return(empty)
  w[, `:=`(sxx = Sxx - Sx^2 / k, syy = Syy - Sy^2 / k, sxy = Sxy - Sx * Sy / k)]
  w[, b1 := sxy / sxx]
  w[, r2 := b1 * sxy / syy]
  w[, r2adj := 1 - (1 - r2) * (k - 1) / (k - 2)]
  w <- w[is.finite(b1) & b1 <= 0 & is.finite(r2adj)]            # NonCompart::Slope: b1 > 0 창 제외, R2ADJ 유한
  if (!nrow(w)) return(empty)
  w[, best := max(r2adj), by = id]
  w <- w[abs(best - r2adj) < tie_tol]                            # NonCompart: abs(max − R2ADJ) < TOL
  w <- w[w[, .I[which.max(k)], by = id]$V1]                      # 점 수 최대(같은 k는 한 창뿐)
  w[, b0 := (y0 + Sy / k - b1 * Sx / k) - b1 * x0]               # 원래 척도의 절편: ln C = b0 + b1·t
  w[b1 < 0, .(id, Lambda_z = -b1, Rsq = r2, Rsq_adjusted = r2adj, No_points_lambda_z = as.integer(k), Lambda_z_lower = lo, Lambda_z_upper = up, b0)]
}

run_nca <- function(sim, mode = c("standard", "bemaster_compat"), rules = get_nca_rules()) {
  mode <- match.arg(mode)
  if (mode == "bemaster_compat") return(run_nca_legacy(sim, mode = "bemaster_compat", rules = rules))
  stopifnot(all(c("id", "time", "conc") %in% names(sim)))
  ids <- sort(unique(sim$id))
  rel <- rules$standard$reliability; lr <- rules$standard$lambda_z
  d <- preprocess_blq(sim)
  s <- d[conc > 0, { i <- which.max(conc); .(n_quant = .N, Cmax = conc[i], Tmax = time[i], Tlast = time[.N], Clast = conc[.N]) }, by = id]
  # AUClast: Linear Up Log Down
  setorder(d, id, time)
  d[, `:=`(t_prev = shift(time), c_prev = shift(conc)), by = id]
  auc <- d[!is.na(t_prev), .(AUClast = sum(.seg_auc(t_prev, c_prev, time, conc, "linear_up_log_down"))), by = id]
  s <- merge(s, auc, by = "id", all.x = TRUE); s[is.na(AUClast), AUClast := 0]
  # λz: Cmax(첫 최대) 이후 양수 농도
  cand <- merge(d[conc > 0, .(id, time, conc)], s[, .(id, Tmax)], by = "id")[time > Tmax, .(id, time, conc)]
  lz <- .lambda_z_best_fit(cand, min_points = as.integer(lr$min_points), tie_tol = lr$tie_tolerance)
  s <- merge(s, lz, by = "id", all.x = TRUE)
  s[, HL_Lambda_z := log(2) / Lambda_z]
  s[, Clast_pred := exp(b0 - Lambda_z * Tlast)]
  s[, `:=`(AUCINF_obs = AUClast + Clast / Lambda_z, AUCINF_pred = AUClast + Clast_pred / Lambda_z)]
  s[, `:=`(`AUC_%Extrap_obs` = (1 - AUClast / AUCINF_obs) * 100, `AUC_%Extrap_pred` = (1 - AUClast / AUCINF_pred) * 100)]
  s[, Span_ratio := (Lambda_z_upper - Lambda_z_lower) / HL_Lambda_z]
  s[, b0 := NULL]
  s[, lambda_ok := !is.na(Lambda_z)]
  s[, `:=`(flag_rsq = lambda_ok & Rsq_adjusted < rel$adj_r2_min, flag_extrap = lambda_ok & `AUC_%Extrap_obs` > rel$extrap_max_pct,
           flag_span = lambda_ok & Span_ratio < rel$span_ratio_min)]
  s[, reliable := lambda_ok & !flag_rsq & !flag_extrap & !flag_span]
  none <- setdiff(ids, s$id)
  if (length(none)) s <- rbind(s, data.table(id = none, n_quant = 0L, lambda_ok = FALSE, reliable = FALSE, flag_rsq = FALSE, flag_extrap = FALSE, flag_span = FALSE), fill = TRUE)
  s[, `:=`(tmax = Tmax, tlast = Tlast, lambda_z = Lambda_z, adj_r2 = Rsq_adjusted, n_lambda = No_points_lambda_z, lambda_span = Lambda_z_upper - Lambda_z_lower,
           t_half = HL_Lambda_z, AUCinf = AUCINF_obs, AUCinf_pred = AUCINF_pred, pct_extrap = `AUC_%Extrap_obs`)]
  setorder(s, id)
  setcolorder(s, NCA_OUT_COLS)
  s[, NCA_OUT_COLS, with = FALSE]
}

# 탈락 사유(중복 포함) 요약: λz 산출 불가, 각 플래그, 조합
dropout_reasons <- function(nca, by = NULL) {
  x <- copy(nca)[n_quant > 0]
  x[, reason := fifelse(!lambda_ok, "λz 산출 불가",
                 fifelse(reliable, "신뢰(플래그 없음)",
                         paste0(fifelse(flag_rsq, "Rsq_adj<0.80 ", ""), fifelse(flag_extrap, "Extrap>20% ", ""), fifelse(flag_span, "Span<2", ""))))]
  x[, reason := trimws(reason)]
  x[, .(n = .N, lambda_fail_pct = 100 * mean(!lambda_ok), flag_rsq_pct = 100 * mean(flag_rsq), flag_extrap_pct = 100 * mean(flag_extrap),
        flag_span_pct = 100 * mean(flag_span), any_flag_or_fail_pct = 100 * mean(!reliable)), by = by]
}
