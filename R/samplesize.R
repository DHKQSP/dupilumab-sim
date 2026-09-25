# samplesize.R — P2 동시 검정력 분석식 (추가 지시 2026-09-26 §3, config/prereg_20260926.yaml section3, D-055).
# scripts/47_sample_size.R의 같은 이름 계산과 같은 식이며, scripts/49가 저장 격자(ss_power_analytic.csv)와 1e-10 이내 일치를 검사한다.
# 자료 생성 모형: log Y = tau·T + beta·(log WT − 평균) + e (AUClast, Cmax). 층 = WT > 분할점, 층 안 1:1 배정.

# log WT의 전체·층 간·층 안 분산(절단 정규 체중, 수치 적분)
ss_weight_moments <- function(wt_mean, wt_sd, lo, hi, split) {
  f <- function(w) dnorm(w, wt_mean, wt_sd)
  I <- function(a, b, g) integrate(function(w) g(w) * f(w), a, b, rel.tol = 1e-12)$value
  Z <- I(lo, hi, function(w) 1); p1 <- I(lo, split, function(w) 1) / Z; p2 <- 1 - p1
  m1 <- I(lo, split, log) / (p1 * Z); m2 <- I(split, hi, log) / (p2 * Z)
  v1 <- I(lo, split, function(w) log(w)^2) / (p1 * Z) - m1^2; v2 <- I(split, hi, function(w) log(w)^2) / (p2 * Z) - m2^2
  ma <- p1 * m1 + p2 * m2
  data.table(p_stratum2 = p2, mean_lwt = ma, var_between = p1 * (m1 - ma)^2 + p2 * (m2 - ma)^2, var_within = p1 * v1 + p2 * v2)[, var_total := var_between + var_within][]
}
# 칸의 분산 구조. inp: sd_cmax, ratio_sd_cmax_auc, rho_total, beta_auc, beta_cmax. wm: ss_weight_moments
ss_cell_cov <- function(inp, cv, cmax_rule, wm) {
  sA <- sqrt(log(1 + (cv / 100)^2)); sC <- if (cmax_rule == "fixed") inp$sd_cmax else inp$ratio_sd_cmax_auc * sA
  cov_t <- inp$rho_total * sA * sC; vb <- wm$var_between; vt <- wm$var_total
  list(sA = sA, sC = sC, cov_t = cov_t, vwA = sA^2 - inp$beta_auc^2 * vb, vwC = sC^2 - inp$beta_cmax^2 * vb, cw = cov_t - inp$beta_auc * inp$beta_cmax * vb,
       veA = sA^2 - inp$beta_auc^2 * vt, veC = sC^2 - inp$beta_cmax^2 * vt, ce = cov_t - inp$beta_auc * inp$beta_cmax * vt)
}
# 이변량 정규 직사각형(SE는 기댓값 고정). M1: 추정량·SE 층 안 분산, df 2n − 3. M0: 추정량 층 안, SE 총분산, df 2n − 2. 확률(0–1)
ss_p2_power <- function(n, cv, gmr_a, gmr_c, inp, am, cmax_rule, wm, ci_level = 0.90, limits = c(0.80, 1.25)) {
  cc <- ss_cell_cov(inp, cv, cmax_rule, wm); k <- 2 / n
  sdA <- sqrt(cc$vwA * k); sdC <- sqrt(cc$vwC * k); rho <- cc$cw / sqrt(cc$vwA * cc$vwC)
  if (am == "M1") { seA <- sdA; seC <- sdC; df <- 2 * n - 3 } else { seA <- sqrt(cc$sA^2 * k); seC <- sqrt(cc$sC^2 * k); df <- 2 * n - 2 }
  q <- qt(1 - (1 - ci_level) / 2, df)
  lA <- (log(limits[1]) + q * seA - log(gmr_a)) / sdA; uA <- (log(limits[2]) - q * seA - log(gmr_a)) / sdA
  lC <- (log(limits[1]) + q * seC - log(gmr_c)) / sdC; uC <- (log(limits[2]) - q * seC - log(gmr_c)) / sdC
  if (lA >= uA || lC >= uC) return(0)
  r1 <- sqrt(1 - rho^2)
  integrate(function(z) dnorm(z) * (pnorm((uC - rho * z) / r1) - pnorm((lC - rho * z) / r1)), lA, uA, rel.tol = 1e-10)$value
}
