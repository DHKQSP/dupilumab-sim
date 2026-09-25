# criteria.R — λz(말단 소실 속도 상수) 신뢰 기준 세트 (추가 지시 2026-09-26 두 번째 지시 §1, config/prereg_20260926.yaml section4, D-057).
# Phoenix WinNonlin의 Lambda Z Acceptance Criteria(Rules 탭)는 선택 항목이며 사용자가 값을 넣고 미달 프로필은 표시만 된다. 아래 값은 SAP 관행값이다.
#   (i)   adjusted R² ≥ 0.80, AUC_%Extrap_obs ≤ 20%              (이 저장소의 1차 세트; rel_i)
#   (ii)  (i) + span ratio ≥ 2                                   (기존 사전 고정 세트; NCA 출력 reliable)
#   (iii) adjusted R² ≥ 0.90, AUC_%Extrap_obs ≤ 20%              (공개 SAP에서 흔한 0.90)
#   (iv)  (iii) + span ratio ≥ 3
# 모든 세트는 λz 산출 가능을 전제로 한다(산출 불가·NA는 미충족).
CRIT_SETS <- list(
  i   = list(adj_r2_min = 0.80, extrap_max_pct = 20, span_ratio_min = NA_real_),
  ii  = list(adj_r2_min = 0.80, extrap_max_pct = 20, span_ratio_min = 2),
  iii = list(adj_r2_min = 0.90, extrap_max_pct = 20, span_ratio_min = NA_real_),
  iv  = list(adj_r2_min = 0.90, extrap_max_pct = 20, span_ratio_min = 3))
CRIT_LABEL_EN <- c(i = "(i) adjusted R-squared at least 0.80, extrapolation at most 20%", ii = "(ii) set (i) and span ratio at least 2",
                   iii = "(iii) adjusted R-squared at least 0.90, extrapolation at most 20%", iv = "(iv) set (iii) and span ratio at least 3")

# 대상자별 충족 여부(논리 벡터). nca: run_nca 출력(lambda_ok, Rsq_adjusted, AUC_%Extrap_obs, Span_ratio)
crit_ok <- function(nca, set) {
  s <- CRIT_SETS[[set]]; if (is.null(s)) stop("crit_ok: 알 수 없는 기준 세트 ", set)
  ok <- nca$lambda_ok %in% TRUE & !is.na(nca$Rsq_adjusted) & nca$Rsq_adjusted >= s$adj_r2_min &
    !is.na(nca[["AUC_%Extrap_obs"]]) & nca[["AUC_%Extrap_obs"]] <= s$extrap_max_pct
  if (!is.na(s$span_ratio_min)) ok <- ok & !is.na(nca$Span_ratio) & nca$Span_ratio >= s$span_ratio_min
  ok %in% TRUE
}
# 시험 평가변수(세트 (iii)·(iv)): 규칙 A(미충족 제외), 규칙 C(미충족이면 AUClast 대입). 규칙 B는 세트와 무관(AUCinf_B)
OC_ENDPOINTS_CRIT <- c("AUCinf_Aiii", "AUCinf_Aiv", "AUCinf_Ciii", "AUCinf_Civ")

# 기준 세트별 미달 요약(개인 수준). nca: run_nca + attach_truth 결과(+ WT). by: 그룹 열(예: 체중 밴드). 반환: 그룹 × 세트 행
crit_summary <- function(nca, by = NULL, per_arm = 117L) {
  geo <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
  rbindlist(lapply(names(CRIT_SETS), function(s_) {
    x <- copy(nca); x[, ok_ := crit_ok(x, s_)]
    x[, {
      f <- !ok_; lz <- !(lambda_ok %in% TRUE)
      tt <- if (sum(f) >= 2 && sum(!f) >= 2) t.test(WT[f], WT[!f])$conf.int else c(NA_real_, NA_real_)
      r <- AUClast / AUCinf_true
      .(set = s_, n = .N, fail_pct = 100 * mean(f), lz_fail_pct = 100 * mean(lz), est_fail_pct = 100 * mean(f & !lz), fail_per_arm = per_arm * mean(f),
        wt_fail_mean = if (any(f)) mean(WT[f]) else NA_real_, wt_retained_mean = if (any(!f)) mean(WT[!f]) else NA_real_,
        wt_diff = if (any(f) && any(!f)) mean(WT[f]) - mean(WT[!f]) else NA_real_, wt_diff_lo = tt[1], wt_diff_hi = tt[2],
        auclast_gm_ratio_fail_to_retained = if (any(f) && any(!f)) geo(AUClast[f]) / geo(AUClast[!f]) else NA_real_,
        extrap_true_median_fail = if (any(f)) median(pct_extrap_true[f], na.rm = TRUE) else NA_real_, extrap_true_median_retained = if (any(!f)) median(pct_extrap_true[!f], na.rm = TRUE) else NA_real_,
        auclast_over_true_fail_median = if (any(f)) median(r[f], na.rm = TRUE) else NA_real_, auclast_over_true_fail_p05 = if (any(f)) quantile(r[f], 0.05, na.rm = TRUE, names = FALSE) else NA_real_) }, by = by]
  }))
}
# Pillar 1(총노출 포착) 요약
pillar1_summary <- function(nca, by = NULL) nca[, .(n = .N, extrap_true_median = median(pct_extrap_true, na.rm = TRUE), extrap_true_p95 = quantile(pct_extrap_true, 0.95, na.rm = TRUE, names = FALSE),
                                                    extrap_true_max = max(pct_extrap_true, na.rm = TRUE), coverage_lt80_pct = 100 * mean(coverage_true < 0.80, na.rm = TRUE),
                                                    coverage_min = min(coverage_true, na.rm = TRUE), wt_mean = mean(WT), wt_median = median(WT)), by = by]
weight_band <- function(WT) factor(ifelse(WT < 60, "below 60", ifelse(WT <= 90, "60-90", ifelse(WT <= 100, "above 90-100", "above 100"))), c("below 60", "60-90", "above 90-100", "above 100"))
