# population.R — 가상 피험자 생성 (SPEC §5.1, §3.5)
# 체중: 성별 혼합 + 포함기준 절단 정규분포. IIV: 로그정규(eta ~ N(0, omega^2)), F는 logit-normal.
# 모든 난수는 호출자가 set.seed/with_seed로 제어한다.

rtrunc_norm <- function(n, mean, sd, lo, hi) {
  if (requireNamespace("truncnorm", quietly = TRUE)) return(truncnorm::rtruncnorm(n, a = lo, b = hi, mean = mean, sd = sd))
  # 역변환 샘플링 (truncnorm 없을 때)
  u <- runif(n, pnorm(lo, mean, sd), pnorm(hi, mean, sd)); qnorm(u, mean, sd)
}

draw_weights <- function(n, design) {
  pop <- design$population
  male <- rbinom(n, 1, pop$sex_ratio_male$value) == 1
  lo <- pop$weight_kg$inclusion$min; hi <- pop$weight_kg$inclusion$max
  wt <- numeric(n)
  wt[male]  <- rtrunc_norm(sum(male),  pop$weight_kg$male$mean,   pop$weight_kg$male$sd,   lo, hi)
  wt[!male] <- rtrunc_norm(sum(!male), pop$weight_kg$female$mean, pop$weight_kg$female$sd, lo, hi)
  data.table(sex = ifelse(male, "M", "F"), WT = wt)
}

# eta 추출: omega(명명 SD 벡터), 상관행렬(선택). 반환 data.table(id, eta_*)
draw_etas <- function(n, omega, corr = NULL) {
  nm <- names(omega)
  sd <- as.numeric(omega); sd[is.na(sd)] <- 0
  if (is.null(corr)) {
    m <- sapply(seq_along(nm), function(j) if (sd[j] > 0) rnorm(n, 0, sd[j]) else rep(0, n))
  } else {
    stopifnot(all(dim(corr) == length(nm)))
    Sig <- diag(sd) %*% corr %*% diag(sd)
    m <- MASS::mvrnorm(n, mu = rep(0, length(nm)), Sigma = Sig)
  }
  m <- matrix(m, nrow = n)
  colnames(m) <- paste0("eta_", nm)
  cbind(data.table(id = seq_len(n)), as.data.table(m))
}

# 한 arm의 가상 피험자: 체중 + eta → individual_params
make_subjects <- function(n, p, design, id_offset = 0L) {
  w  <- draw_weights(n, design)
  et <- draw_etas(n, p$omega, corr = NULL)
  et[, id := id + id_offset]
  et[, `:=`(WT = w$WT, sex = w$sex)]
  ip <- individual_params(p, et)
  ip[, sex := w$sex]
  list(ipar = ip, eta = et)
}

# 체중 층
weight_stratum <- function(WT, design) {
  s <- design$population$weight_strata_kg
  cut(WT, breaks = s$breaks, labels = s$labels, right = FALSE)
}
