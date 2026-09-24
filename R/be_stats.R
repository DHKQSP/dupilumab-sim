# be_stats.R — 평행군 생물학적 동등성 통계 (지시서 §3). log 변환 후 두 표본 pooled t 90% CI, 80.00–125.00%.
# 민감도: 체중 공변량 ANCOVA (lm(log y ~ arm + WT)).

be_pooled_t <- function(y, arm, ci_level = 0.90, limits = c(0.80, 1.25)) {
  ok <- is.finite(y) & y > 0 & !is.na(arm)
  ly <- log(y[ok]); grp <- arm[ok]                     # 인자명(arm)을 부분집합 조건에 직접 쓰지 않는다(D-025)
  yR <- ly[grp == "R"]; yT <- ly[grp == "T"]
  nR <- length(yR); nT <- length(yT)
  if (nR < 2 || nT < 2) return(list(GMR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, pass = NA, n_R = nR, n_T = nT, width = NA_real_, se = NA_real_, df = NA_real_))
  df <- nR + nT - 2
  sp2 <- ((nR - 1) * var(yR) + (nT - 1) * var(yT)) / df
  se <- sqrt(sp2 * (1 / nR + 1 / nT))
  diff <- mean(yT) - mean(yR)
  q <- qt(1 - (1 - ci_level) / 2, df)
  lo <- exp(diff - q * se); hi <- exp(diff + q * se)
  list(GMR = exp(diff), CI_lower = lo, CI_upper = hi, pass = (lo >= limits[1] & hi <= limits[2]), n_R = nR, n_T = nT,
       width = hi - lo, se = se, df = df)
}

be_ancova_weight <- function(y, arm, WT, ci_level = 0.90, limits = c(0.80, 1.25)) {
  ok <- is.finite(y) & y > 0 & !is.na(arm) & is.finite(WT)
  d <- data.frame(ly = log(y[ok]), arm = factor(arm[ok], levels = c("R", "T")), WT = WT[ok])
  if (sum(d$arm == "R") < 2 || sum(d$arm == "T") < 2) return(list(GMR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, pass = NA, n_R = sum(d$arm == "R"), n_T = sum(d$arm == "T"), width = NA_real_, se = NA_real_, df = NA_real_))
  fit <- lm(ly ~ arm + WT, data = d)
  est <- coef(fit)[["armT"]]; se <- sqrt(vcov(fit)["armT", "armT"]); df <- fit$df.residual
  q <- qt(1 - (1 - ci_level) / 2, df)
  lo <- exp(est - q * se); hi <- exp(est + q * se)
  list(GMR = exp(est), CI_lower = lo, CI_upper = hi, pass = (lo >= limits[1] & hi <= limits[2]), n_R = sum(d$arm == "R"), n_T = sum(d$arm == "T"),
       width = hi - lo, se = se, df = df)
}

# nca_arm: data.table with arm, WT, and endpoint columns. endpoints: named list endpoint -> list(col, subset_col)
# 반환: data.table(endpoint, method, GMR, CI_lower, CI_upper, pass, n_R, n_T, width)
be_analyze <- function(nca_arm, endpoints, ci_level = 0.90, limits = c(0.80, 1.25), methods = c("pooled_t")) {
  rbindlist(lapply(names(endpoints), function(ep) {
    spec <- endpoints[[ep]]
    d <- nca_arm
    if (!is.null(spec$subset_col)) d <- d[get(spec$subset_col) %in% TRUE]
    y <- d[[spec$col]]
    rbindlist(lapply(methods, function(m) {
      r <- if (m == "pooled_t") be_pooled_t(y, d$arm, ci_level, limits) else be_ancova_weight(y, d$arm, d$WT, ci_level, limits)
      data.table(endpoint = ep, method = m, GMR = r$GMR, CI_lower = r$CI_lower, CI_upper = r$CI_upper, pass = r$pass,
                 n_R = r$n_R, n_T = r$n_T, width = r$width)
    }))
  }))
}

BE_ENDPOINTS <- list(
  Cmax            = list(col = "Cmax"),
  AUClast         = list(col = "AUClast"),
  AUCinf_all      = list(col = "AUCinf", subset_col = "lambda_ok"),
  AUCinf_reliable = list(col = "AUCinf", subset_col = "reliable"),
  AUCinf_true     = list(col = "AUCinf_true")
)
