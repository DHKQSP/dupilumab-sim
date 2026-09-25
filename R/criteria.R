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
