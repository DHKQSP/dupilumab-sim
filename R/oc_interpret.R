# oc_interpret.R — P2 경계 1종 오류 해석 도우미 (지시 2026-09-25 §3, scripts/43_p2_interpretation.R).
# (1) 불편 추정량의 이론 기준: 두 단측 검정(90% CI가 80–125% 안)의 경계 통과 확률(t 기반, se 고정)
# (2) 정규 근사 예측 통과율(시험 간 분포 평균·SD와 시험 안 se가 다를 때)
# (3) Wilson 95% 구간 대 5% 분류, (4) 편향 방향(1 쪽/1에서 멀어짐/0 포함)

# (1) 이론 기준. log GMR 추정량 D가 참값 mu = log(true_ratio)에 대해 불편이고 (D − mu)/se ~ t_df(se 고정)이면
#   통과 = {log L + q·se < D < log U − q·se},  q = qt(1 − (1 − ci_level)/2, df)
#   P(통과) = P(T > (log L − mu)/se + q) − P(T > (log U − mu)/se − q)   (0 미만이면 0: 구간이 너무 넓어 통과 불가)
# 가까운 한계(near): 참값에 가까운 쪽 한계를 넘어 CI가 안쪽에 들어올 확률. 참값이 정확히 한계(0.80 또는 1.25)면 (1 − ci_level)/2 = 5%.
# 먼 한계(far): 반대쪽 한계에서 불통과할 확률(참값 0.80이면 상한 > 1.25, 1.25면 하한 < 0.80). 경계에서 P(통과) = 5% − P(far)(정확한 항등식).
# 반환: data.table(q, p_pass, p_near, p_far) — 확률(0–1), 입력 길이만큼(재활용 규칙)
be_boundary_type1_theory <- function(se, df, true_ratio, ci_level = 0.90, limits = c(0.80, 1.25)) {
  if (!length(se) || any(!is.finite(se) | se <= 0)) stop("be_boundary_type1_theory: se는 양의 유한값이어야 합니다")
  if (any(!is.finite(df) | df <= 0)) stop("be_boundary_type1_theory: df는 양수여야 합니다")
  if (any(!is.finite(true_ratio) | true_ratio <= 0)) stop("be_boundary_type1_theory: true_ratio는 양수여야 합니다")
  if (length(limits) != 2 || !(limits[1] < limits[2]) || limits[1] <= 0) stop("be_boundary_type1_theory: limits = c(하한, 상한), 0 < 하한 < 상한")
  n <- max(length(se), length(df), length(true_ratio))
  se <- rep_len(se, n); df <- rep_len(df, n); true_ratio <- rep_len(true_ratio, n)
  q <- qt(1 - (1 - ci_level) / 2, df)
  mu <- log(true_ratio); a <- log(limits[1]); b <- log(limits[2])
  zl <- (a - mu) / se + q                                      # 하한 > L  <=>  T > zl
  zu <- (b - mu) / se - q                                      # 상한 < U  <=>  T < zu
  p_lo_ok <- pt(zl, df, lower.tail = FALSE); p_up_fail <- pt(zu, df, lower.tail = FALSE)
  p_up_ok <- pt(zu, df); p_lo_fail <- pt(zl, df)
  near_low <- abs(mu - a) <= abs(mu - b)
  data.table(q = q, p_pass = pmax(0, p_lo_ok - p_up_fail),   # 두 꼬리 확률의 차(1 근처 뺄셈 없음: 작은 p_far도 정밀)
             p_near = ifelse(near_low, p_lo_ok, p_up_ok), p_far = ifelse(near_low, p_up_fail, p_lo_fail))
}

# (2) 정규 근사 예측 통과율: D ~ N(mean_log, sd_log²), CI 반폭 q·se(se 고정). 시험 간 SD(sd_log)가 시험 안 se와 다르거나 D가 치우칠 때의 설명용.
# sd_log = se, df → ∞이면 (1)과 같다(테스트). 반환: 확률(0–1)
be_pass_prob_normal <- function(mean_log, sd_log, se, df, ci_level = 0.90, limits = c(0.80, 1.25)) {
  if (any(!is.finite(sd_log) | sd_log <= 0) || any(!is.finite(se) | se <= 0)) stop("be_pass_prob_normal: sd_log, se는 양의 유한값이어야 합니다")
  q <- qt(1 - (1 - ci_level) / 2, df)
  lo <- (log(limits[1]) + q * se - mean_log) / sd_log; hi <- (log(limits[2]) - q * se - mean_log) / sd_log
  pmax(0, pnorm(lo, lower.tail = FALSE) - pnorm(hi, lower.tail = FALSE))
}

# (3) 경계 1종 오류 분류(Wilson 95% 구간, %): 상한 < alpha → conservative, 하한 > alpha → exceeding, 그 밖(구간이 alpha 포함, 끝점 같음 포함) → nominal.
# NA 구간은 NA. 하한 > 상한이면 중단.
OC_TYPE1_CLASSES <- c("conservative", "nominal", "exceeding")
classify_type1 <- function(lo, hi, alpha_pct = 5) {
  if (length(lo) != length(hi)) stop("classify_type1: lo, hi 길이가 다릅니다")
  if (!identical(is.na(lo), is.na(hi))) stop("classify_type1: lo, hi의 NA 위치가 다릅니다")
  if (any(lo > hi, na.rm = TRUE)) stop("classify_type1: 하한 > 상한")
  fifelse(hi < alpha_pct, "conservative", fifelse(lo > alpha_pct, "exceeding", "nominal"))
}

# (4) 편향 방향: 참값 비 < 1(예: 0.80)이면 음의 편향이 1에서 멀어짐(범위 바깥쪽), 양이 1 쪽. 참값 비 > 1이면 반대.
# 95% 구간 [lo, hi]가 0을 포함하면 "none". 반환: "toward_1" | "away_from_1" | "none"
bias_direction <- function(bias_log, true_ratio, lo = bias_log, hi = bias_log) {
  if (any(!is.finite(true_ratio) | true_ratio <= 0 | true_ratio == 1)) stop("bias_direction: 참값 비는 1이 아닌 양수여야 합니다")
  s <- sign(bias_log) * sign(-log(true_ratio))                 # +1 = 1 쪽
  fifelse(lo <= 0 & hi >= 0, "none", fifelse(s > 0, "toward_1", "away_from_1"))
}
