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
  if (identical(wt_spec$dist, "lognormal")) {                # 아토피 성인 분포(config/population_atopic.yaml): 산술 평균·SD → meanlog·sdlog, 절단 구간 역누적분포
    b <- as.numeric(wt_spec$trunc); sl <- sqrt(log(1 + (wt_spec$sd / wt_spec$mean)^2)); ml <- log(wt_spec$mean) - sl^2 / 2
    u <- runif(n, plnorm(b[1], ml, sl), plnorm(b[2], ml, sl))
    return(data.table(sex = ifelse(male, "M", "F"), WT = qlnorm(u, ml, sl)))
  }
  if (identical(wt_spec$dist, "uniform")) {                  # §6 체중 밴드 균등분포
    b <- as.numeric(wt_spec$trunc); return(data.table(sex = ifelse(male, "M", "F"), WT = runif(n, b[1], b[2])))
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
  ada_vec <- if (ada_fraction > 0) as.integer(runif(nrow(e)) < ada_fraction) else 0L   # 난수 순서 동일(D-025)
  e[, ada := ada_vec]
  if (!is.null(wt_spec$height)) {                            # §6: 키·BMI (지정된 분포에서만 추출 → 기존 실행의 난수 순서 불변)
    h <- wt_spec$height; ht <- rtrunc_norm(nrow(e), h$mean, h$sd, h$trunc[1], h$trunc[2])
    e[, `:=`(HT = ht, BMI = WT / (ht / 100)^2)]
  }
  e[]
}

# 층 경계: 활성 체중 분포의 절단 범위에서 자동 생성 [lo, split], (split, hi] (검토 의견 3차 §4, D-024).
weight_range_of_spec <- function(wt_spec) {
  b <- c(if (!is.null(wt_spec$trunc)) as.numeric(wt_spec$trunc), if (!is.null(wt_spec$male_bounds)) as.numeric(wt_spec$male_bounds),
         if (!is.null(wt_spec$female_bounds)) as.numeric(wt_spec$female_bounds))
  if (!length(b) || any(!is.finite(b))) stop("층화: 체중 분포에 유한한 절단 범위가 없습니다")
  lo <- min(b); if (!is.null(wt_spec$female_min) && is.null(wt_spec$trunc)) lo <- min(lo, as.numeric(wt_spec$female_min))
  c(lo, max(b))
}
strata_from_weight_spec <- function(wt_spec, split_kg) {
  r <- weight_range_of_spec(wt_spec)
  if (!(split_kg > r[1] && split_kg < r[2])) stop(sprintf("층화: 분할점 %g kg이 체중 범위 [%g, %g] 안에 있지 않습니다", split_kg, r[1], r[2]))
  list(breaks = c(r[1], split_kg, r[2]), labels = c(sprintf("%g-%g", r[1], split_kg), sprintf(">%g-%g", split_kg, r[2])))
}
weight_stratum <- function(WT, breaks, labels) cut(WT, breaks = as.numeric(breaks), labels = labels, include.lowest = TRUE, right = TRUE)

# 층화 무작위배정 1:1. subj에 WT 필요. 반환: arm 열 추가(R/T), 각 arm 정확히 n/2 (홀수 층은 층 간 교대 배정으로 보정)
assign_arms_stratified <- function(subj, breaks, labels) {
  subj <- copy(subj)
  subj[, stratum := weight_stratum(WT, breaks, labels)]
  if (anyNA(subj$stratum)) stop(sprintf("층화 배정: 층 밖 대상자 %d명 (체중 범위 [%g, %g] 밖)", sum(is.na(subj$stratum)), min(breaks), max(breaks)))
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
