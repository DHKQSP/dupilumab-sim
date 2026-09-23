# population.R — 가상 피험자 (지시서 §4). 체중: 절단 정규분포(성별 하한 옵션), IIV: 로그정규, 층화 무작위배정.
# 모든 난수는 호출자가 with_seed로 제어한다.

rtrunc_norm <- function(n, mean, sd, lo, hi) {
  if (n == 0) return(numeric(0))
  if (requireNamespace("truncnorm", quietly = TRUE)) return(truncnorm::rtruncnorm(n, a = lo, b = hi, mean = mean, sd = sd))
  u <- runif(n, pnorm(lo, mean, sd), pnorm(hi, mean, sd)); qnorm(u, mean, sd)
}

# wt_spec: list(mean, sd, trunc = c(lo, hi), female_min = NULL, male_bounds = NULL, female_bounds = NULL)
draw_weights <- function(n, wt_spec, sex_ratio_male = 0.5) {
  male <- rbinom(n, 1, sex_ratio_male) == 1
  bounds <- function(is_male) {
    b <- if (!is.null(wt_spec$trunc)) as.numeric(wt_spec$trunc) else c(-Inf, Inf)
    if (is_male && !is.null(wt_spec$male_bounds)) b <- as.numeric(wt_spec$male_bounds)
    if (!is_male && !is.null(wt_spec$female_bounds)) b <- as.numeric(wt_spec$female_bounds)
    if (!is_male && !is.null(wt_spec$female_min)) b[1] <- max(b[1], as.numeric(wt_spec$female_min))
    b
  }
  wt <- numeric(n)
  bm <- bounds(TRUE); bf <- bounds(FALSE)
  wt[male]  <- rtrunc_norm(sum(male),  wt_spec$mean, wt_spec$sd, bm[1], bm[2])
  wt[!male] <- rtrunc_norm(sum(!male), wt_spec$mean, wt_spec$sd, bf[1], bf[2])
  data.table(sex = ifelse(male, "M", "F"), WT = wt)
}

draw_etas <- function(n, omega) {
  nm <- names(omega); sd <- as.numeric(omega); sd[is.na(sd)] <- 0
  m <- sapply(seq_along(nm), function(j) if (sd[j] > 0) rnorm(n, 0, sd[j]) else rep(0, n))
  m <- matrix(m, nrow = n); colnames(m) <- paste0("eta_", nm)
  cbind(data.table(id = seq_len(n)), as.data.table(m))
}

# n명의 가상 피험자: id, sex, WT, eta_*, ada
make_subjects <- function(n, p, wt_spec, sex_ratio_male = 0.5, ada_fraction = 0) {
  w <- draw_weights(n, wt_spec, sex_ratio_male)
  e <- draw_etas(n, p$omega)
  e[, `:=`(sex = w$sex, WT = w$WT)]
  e[, ada := if (ada_fraction > 0) as.integer(runif(n) < ada_fraction) else 0L]
  e[]
}

# 층 경계의 바깥 끝은 열어 둔다(-Inf, Inf): 선정 기준 밖 체중(민감도·교차검증 모집단)도 반드시 층에 들어가야 배정에서 빠지지 않는다(D-020)
weight_stratum <- function(WT, breaks, labels) {
  b <- as.numeric(breaks); b[1] <- -Inf; b[length(b)] <- Inf
  cut(WT, breaks = b, labels = labels, include.lowest = TRUE, right = TRUE)
}

# 층화 무작위배정 1:1. subj에 WT 필요. 반환: arm 열 추가(R/T), 각 arm 정확히 n/2 (홀수 층은 층 간 교대 배정으로 보정)
assign_arms_stratified <- function(subj, breaks, labels) {
  subj <- copy(subj)
  subj[, stratum := weight_stratum(WT, breaks, labels)]
  subj[, arm := NA_character_]
  leftover <- integer(0)
  for (s in levels(subj$stratum)) {
    idx <- subj[stratum == s, which = TRUE]
    idx <- idx[sample.int(length(idx))]
    npair <- floor(length(idx) / 2)
    if (npair > 0) {
      subj[idx[seq_len(npair)], arm := "R"]
      subj[idx[npair + seq_len(npair)], arm := "T"]
    }
    if (length(idx) %% 2 == 1) leftover <- c(leftover, idx[length(idx)])
  }
  if (length(leftover)) {
    leftover <- leftover[sample.int(length(leftover))]
    subj[leftover, arm := rep(c("R", "T"), length.out = length(leftover))]
  }
  if (anyNA(subj$arm)) stop("층화 배정에서 배정되지 않은 대상자가 있습니다: ", sum(is.na(subj$arm)))
  subj[]
}

# config의 weight 항목 → wt_spec
weight_spec_from_design <- function(design, which = c("base", "sensitivity")) {
  which <- match.arg(which)
  w <- design$weight[[which]]
  list(mean = w$mean, sd = w$sd, trunc = as.numeric(w$trunc), female_min = w$female_min)
}
