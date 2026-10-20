# m15.R — ICH M15 분석(plans/M15_analysis_plan_20261020.md, config/prereg_20261020_m15.yaml, D-068) 전용 함수.
# source_project()에는 넣지 않는다(보류 중인 nlmixr2 작업의 R/00_setup.R 수정과 충돌 방지). scripts/70~75와 tests/testthat/test-m15.R가 직접 source한다.
# 모델 구조·추정값은 바꾸지 않는다. 개체 파라미터는 기존 individual_params(), 풀이는 기존 solve_model()/get_model(), NCA 규칙은 기존 R/nca.R를 쓴다.

M15_ETA_NAMES <- c("Vc", "ke", "ka", "Vmax", "MTT", "k12", "k21", "Km", "F")   # 표준 정규 점수의 고정 순서(두 모델이 같은 대상자를 공유)

m15_cfg <- function() read_cfg("prereg_20261020_m15.yaml")

# ---------------------------------------------------------------- 체중·키 ----------------------------------------------------------------
# w: list(dist = range_truncnorm | truncnorm | uniform, lo, hi[, mean, sd]) → list(dist = truncnorm | uniform, mean, sd, lo, hi)
m15_weight_spec <- function(w) {
  d <- w$dist
  if (identical(d, "range_truncnorm")) return(list(dist = "truncnorm", lo = as.numeric(w$lo), hi = as.numeric(w$hi),
                                                   mean = (as.numeric(w$lo) + as.numeric(w$hi)) / 2, sd = (as.numeric(w$hi) - as.numeric(w$lo)) / 4))
  if (identical(d, "uniform")) return(list(dist = "uniform", lo = as.numeric(w$lo), hi = as.numeric(w$hi)))
  if (identical(d, "truncnorm")) return(list(dist = "truncnorm", mean = as.numeric(w$mean), sd = as.numeric(w$sd),
                                             lo = if (is.null(w$lo)) -Inf else as.numeric(w$lo), hi = if (is.null(w$hi)) Inf else as.numeric(w$hi)))
  stop("m15_weight_spec: 알 수 없는 분포 ", d)
}

# 역누적분포 체중: u(0,1) → 성별 범위 [lo_s, hi_s]와 분포 범위의 교집합에서 절단한 분포의 분위수
m15_q_weight <- function(u, spec, lo_s, hi_s) {
  L <- pmax(spec$lo, lo_s); U <- pmin(spec$hi, hi_s)
  if (any(L >= U)) stop("m15_q_weight: 빈 체중 범위")
  if (spec$dist == "uniform") return(L + u * (U - L))
  a <- pnorm((L - spec$mean) / spec$sd); b <- pnorm((U - spec$mean) / spec$sd)
  spec$mean + spec$sd * qnorm(a + u * (b - a))
}

# n명(등록 순서): 성별, 체중(선정 기준 안), 키(체중 정규 점수와 상관 rho, Gaussian copula), BMI, 표준 정규 점수 z_<파라미터>.
# 난수 순서 고정: u_sex(n), u_w(n), z_h(n), z_eta(n x 9). 집단과 무관 → 같은 시드면 집단 간 공통 난수.
# pop: prereg populations$S1 등(weight, height). incl: list(male_kg, female_kg). wkey: "weight" 또는 "weight_sensitivity"
m15_subjects <- function(n, pop, incl, rho, sex_ratio_male, seed, wkey = "weight") {
  r <- with_seed(seed, list(u_sex = runif(n), u_w = runif(n), z_h = rnorm(n), z = matrix(rnorm(n * length(M15_ETA_NAMES)), n)))
  spec <- m15_weight_spec(pop[[wkey]])
  male <- r$u_sex < sex_ratio_male
  wt <- numeric(n)
  im <- as.numeric(unlist(incl$male_kg)); iff <- as.numeric(unlist(incl$female_kg))
  wt[male] <- m15_q_weight(r$u_w[male], spec, im[1], im[2])
  wt[!male] <- m15_q_weight(r$u_w[!male], spec, iff[1], iff[2])
  zw <- qnorm(r$u_w); zh <- rho * zw + sqrt(1 - rho^2) * r$z_h
  h <- pop$height
  if (isTRUE(h$by_sex)) { mu <- ifelse(male, h$male$mean_cm, h$female$mean_cm); sdv <- ifelse(male, h$male$sd_cm, h$female$sd_cm) } else { mu <- h$mean_cm; sdv <- h$sd_cm }
  ht <- mu + sdv * zh
  out <- data.table(id = seq_len(n), sex = ifelse(male, "M", "F"), WT = wt, HT = ht, BMI = wt / (ht / 100)^2)
  zt <- as.data.table(r$z); setnames(zt, paste0("z_", M15_ETA_NAMES))
  cbind(out, zt)
}

# 표준 정규 점수 → 모델별 η(로그 척도 SD = p$omega). IIV가 0인 파라미터는 0. individual_params() 입력 형태(id, WT, eta_*, ada).
m15_eta_table <- function(subj, p) {
  e <- data.table(id = subj$id, WT = subj$WT)
  for (nm in names(p$omega)) {
    s <- as.numeric(p$omega[[nm]]); if (is.na(s)) s <- 0
    zc <- paste0("z_", nm)
    e[, (paste0("eta_", nm)) := if (s > 0) { if (!zc %in% names(subj)) stop("표준 정규 점수 없음: ", nm); s * subj[[zc]] } else 0]
  }
  e[, ada := 0]
  e[]
}

# 등록 집단 체중 분위수(층 경계 중앙값과 [8] 백분위수): 200,000명
m15_population_quantiles <- function(pop, incl, rho, sex_ratio_male, seed, probs, n = 200000L, wkey = "weight") {
  s <- m15_subjects(n, pop, incl, rho, sex_ratio_male, seed, wkey)
  data.table(prob = probs, WT = as.numeric(quantile(s$WT, probs, names = FALSE, type = 7)))
}

# ---------------------------------------------------------------- 배정·탈락 ----------------------------------------------------------------
# 층 안 크기 2의 순열 블록(등록 순서). u_alloc: 대상자별 균등 난수(블록 첫 대상자의 값이 블록 순서를 정한다).
# 반환: stratum(0 = 경계 이하, 1 = 초과), arm, arm_order(군 안 등록 순서)
m15_allocate <- function(WT, boundary, u_alloc) {
  s <- as.integer(WT > boundary)
  arm <- character(length(WT))
  for (k in 0:1) {
    idx <- which(s == k); if (!length(idx)) next
    pos <- seq_along(idx); blk <- (pos + 1L) %/% 2L; first <- pos %% 2L == 1L
    first_arm <- ifelse(u_alloc[idx[first]] < 0.5, "R", "T")
    fa <- first_arm[blk]
    arm[idx] <- ifelse(first, fa, ifelse(fa == "R", "T", "R"))
  }
  ord <- integer(length(WT))
  for (a in c("R", "T")) { i <- which(arm == a); ord[i] <- seq_along(i) }
  data.table(stratum = s, arm = arm, arm_order = ord)
}

# ---------------------------------------------------------------- 채혈 시각·잔차 ----------------------------------------------------------------
# 벡터화 허용창 편차(R/sampling.R jitter_times와 같은 규칙: 시점 0 고정, 0 이상, 계획 순서 유지·최소 간격 1시간). U: n x K 균등 난수.
m15_obs_times <- function(ids, planned, design, U, jitter = TRUE, min_gap = 1 / 24) {
  pl <- sort(unique(c(0, planned))); K <- length(pl); n <- length(ids)
  if (!jitter) return(data.table(id = rep(ids, each = K), planned = rep(pl, n), time = rep(pl, n)))
  stopifnot(nrow(U) == n, ncol(U) == K)
  w <- window_for_time(pl, design)
  A <- matrix(pl, n, K, byrow = TRUE) + (2 * U - 1) * matrix(w, n, K, byrow = TRUE)
  A[, pl == 0] <- 0
  A <- pmax(A, 0)
  if (K > 1) for (k in 2:K) A[, k] <- pmax(A[, k], A[, k - 1] + min_gap)
  data.table(id = rep(ids, each = K), planned = rep(pl, n), time = as.vector(t(A)))
}

# 잔차 표: E_p, E_a (n x K 표준 정규) → data.table(id, planned, eps_p, eps_a)
m15_eps_table <- function(ids, planned, sigma, Zp, Za) {
  pl <- sort(unique(c(0, planned))); K <- length(pl); n <- length(ids)
  data.table(id = rep(ids, each = K), planned = rep(pl, n), eps_p = as.vector(t(Zp)) * sigma[["prop"]], eps_a = as.vector(t(Za)) * sigma[["add"]])
}

# 관측 농도: y = g * C * (1 + eps_p) + eps_a, BLQ = 투여 전 또는 y < LLOQ (R/simulate.R와 같은 식, g = 1이면 동일)
m15_observe <- function(C, planned, eps_p, eps_a, lloq, g = 1) {
  y <- g * C * (1 + eps_p) + eps_a
  blq <- planned == 0 | y < lloq
  fifelse(blq, NA_real_, y)
}

# ---------------------------------------------------------------- NCA 보조 ----------------------------------------------------------------
# AUC0-last, Cmax, Tlast만(run_nca와 같은 BLQ 전처리·Linear Up Log Down). sim: data.table(id, time, conc).
m15_nca_basic <- function(sim) {
  ids <- sort(unique(sim$id))
  d <- preprocess_blq(sim)
  s <- d[conc > 0, { i <- which.max(conc); .(Cmax = conc[i], Tlast = time[.N]) }, by = id]
  setorder(d, id, time)
  d[, `:=`(t_prev = shift(time), c_prev = shift(conc)), by = id]
  auc <- d[!is.na(t_prev), .(AUClast = sum(.seg_auc(t_prev, c_prev, time, conc, "linear_up_log_down"))), by = id]
  s <- merge(s, auc, by = "id", all.x = TRUE); s[is.na(AUClast), AUClast := 0]
  none <- setdiff(ids, s$id)
  if (length(none)) s <- rbind(s, data.table(id = none, Cmax = NA_real_, Tlast = NA_real_, AUClast = NA_real_))
  setorder(s, id)[]
}

# 부분 AUC: BLQ 전처리 뒤 Linear Up Log Down, 구간 경계는 그 구간 선분 안에서 보간(상승 선형, 하강 로그 선형), Tlast 뒤 기여 없음.
# intervals: 이름 있는 list(c(a, b)). 정량값이 하나도 없는 대상자는 NA.
m15_pauc <- function(sim, intervals) {
  ids <- sort(unique(sim$id))
  d <- preprocess_blq(sim); setorder(d, id, time)
  d[, `:=`(t1 = shift(time), c1 = shift(conc)), by = id]
  sg <- d[!is.na(t1), .(id, t1, c1, t2 = time, c2 = conc)]
  out <- data.table(id = ids)
  has <- unique(d[conc > 0, id])
  for (nm in names(intervals)) {
    a <- intervals[[nm]][1]; b <- intervals[[nm]][2]
    u <- pmax(sg$t1, a); v <- pmin(sg$t2, b); ok <- v > u & sg$t2 > sg$t1
    lg <- sg$c2 < sg$c1 & sg$c1 > 0 & sg$c2 > 0
    dt <- sg$t2 - sg$t1
    k <- ifelse(lg, log(sg$c1 / sg$c2) / dt, 0)
    cu <- ifelse(lg, sg$c1 * exp(-k * (u - sg$t1)), sg$c1 + (sg$c2 - sg$c1) * (u - sg$t1) / dt)
    cv <- ifelse(lg, sg$c1 * exp(-k * (v - sg$t1)), sg$c1 + (sg$c2 - sg$c1) * (v - sg$t1) / dt)
    area <- ifelse(lg, (cu - cv) / k, (cu + cv) / 2 * (v - u)); area[!ok] <- 0
    r <- data.table(id = sg$id, A = area)[, .(A = sum(A)), by = id]
    out[, (nm) := 0]
    out[r, on = "id", (nm) := i.A]
    out[!id %in% has, (nm) := NA_real_]
  }
  out[]
}

# ---------------------------------------------------------------- 분석 모형(중첩 표본 수) ----------------------------------------------------------------
# d: data.table(<keys>, arm R/T, arm_order 1..n_max(군마다 연속), ok(평가 가능·유효 y), y(양수 원척도), s(0/1 층)).
# 군 안 등록 순서 누적합으로 표본 수 n(군당 무작위배정)마다 M0(pooled t)와 M1(처리 + 층)을 닫힌 꼴로 푼다(R/oc_models.R와 같은 정규방정식).
# M1에서 층이 비면(분석 집합의 층 하나가 0명) 층 항을 뺀다(= M0, m1_dropped = TRUE).
m15_be_nested <- function(d, keys, n_grid, ci_level = 0.90, limits = c(0.80, 1.25)) {
  x <- d[, c(keys, "arm", "arm_order", "ok", "y", "s"), with = FALSE]
  x[, ok := ok & is.finite(y) & y > 0]
  x[, ly := fifelse(ok, log(pmax(y, 1e-300)), 0)]
  x[, cen := { v <- ly[ok]; if (length(v)) mean(v) else 0 }, by = keys]
  x[, `:=`(w = as.numeric(ok), lyc = fifelse(ok, ly - cen, 0))]
  setorderv(x, c(keys, "arm", "arm_order"))
  x[, `:=`(cn = cumsum(w), cs = cumsum(w * s), cy = cumsum(lyc), csy = cumsum(s * lyc), cyy = cumsum(lyc^2)), by = c(keys, "arm")]
  sel <- x[arm_order %in% n_grid, c(keys, "arm", "arm_order", "cn", "cs", "cy", "csy", "cyy"), with = FALSE]
  fml <- as.formula(paste(paste(c(keys, "arm_order"), collapse = " + "), "~ arm"))
  w <- dcast(sel, fml, value.var = c("cn", "cs", "cy", "csy", "cyy"))
  setnames(w, "arm_order", "n_rand")
  q_ <- function(df) qt(1 - (1 - ci_level) / 2, pmax(df, 1))
  w[, `:=`(n_R = cn_R, n_T = cn_T)]
  # M0
  w[, `:=`(mR = cy_R / n_R, mT = cy_T / n_T)]
  w[, df0 := n_R + n_T - 2]
  w[, ss0 := (cyy_R - n_R * mR^2) + (cyy_T - n_T * mT^2)]
  w[, se0 := sqrt(pmax(ss0, 0) / df0 * (1 / n_R + 1 / n_T))]
  w[, est0 := mT - mR]
  # M1
  N <- w$n_R + w$n_T; nT <- w$n_T; ns <- w$cs_R + w$cs_T; nts <- w$cs_T
  sy <- w$cy_R + w$cy_T; sty <- w$cy_T; ssy <- w$csy_R + w$csy_T; syy <- w$cyy_R + w$cyy_T
  sol <- .solve_sym3(N, nT, ns, nT, nts, ns, sy, sty, ssy)
  df1 <- N - 3
  rss <- syy - (sol$beta1 * sy + sol$beta2 * sty + sol$beta3 * ssy); rss[rss < 0 & rss > -1e-9] <- 0
  se1 <- sqrt(rss / df1 * sol$inv22)
  drop1 <- ns == 0 | ns == N | !is.finite(se1)
  w[, `:=`(est1 = fifelse(drop1, est0, sol$beta2), se1 = fifelse(drop1, se0, se1), df1 = fifelse(drop1, df0, df1), m1_dropped = drop1)]
  ok0 <- w$n_R >= 2 & w$n_T >= 2 & w$df0 >= 1
  mk <- function(est, se, df, model) {
    est[!ok0] <- NA_real_; se[!ok0] <- NA_real_
    lo <- exp(est - q_(df) * se); hi <- exp(est + q_(df) * se)
    o <- w[, c(keys, "n_rand", "n_R", "n_T"), with = FALSE]
    o[, `:=`(model = model, est = est, se = se, df = df, CI_lower = lo, CI_upper = hi, pass = !is.na(lo) & lo >= limits[1] & hi <= limits[2])]
    o
  }
  out <- rbind(mk(w$est0, w$se0, w$df0, "M0"), mk(w$est1, w$se1, w$df1, "M1")[, m1_dropped := w$m1_dropped], fill = TRUE)
  out[model == "M0", m1_dropped := NA]
  out[]
}

# ---------------------------------------------------------------- AUROC ----------------------------------------------------------------
# Mann-Whitney 추정(동점 1/2)과 DeLong 분산. pos: 양성 통계량, neg: 음성 통계량
m15_auroc <- function(pos, neg, conf = 0.95) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  m <- length(pos); n <- length(neg)
  if (m < 2 || n < 2) return(list(auroc = NA_real_, se = NA_real_, lo = NA_real_, hi = NA_real_, n_pos = m, n_neg = n))
  sn <- sort(neg); sp <- sort(pos)
  v10 <- (findInterval(pos, sn, left.open = TRUE) + (findInterval(pos, sn) - findInterval(pos, sn, left.open = TRUE)) / 2) / n
  v01 <- ((m - findInterval(neg, sp)) + (findInterval(neg, sp) - findInterval(neg, sp, left.open = TRUE)) / 2) / m
  a <- mean(v10)
  se <- sqrt(var(v10) / m + var(v01) / n)
  z <- qnorm(1 - (1 - conf) / 2)
  list(auroc = a, se = se, lo = max(0, a - z * se), hi = min(1, a + z * se), n_pos = m, n_neg = n)
}

# ---------------------------------------------------------------- 참값 ----------------------------------------------------------------
# 격자 풀이(id, time, C, auc; time 오름차순)에서 연속 시간 참 Tlast: Cmax 뒤 참 농도가 LLOQ 아래로 처음 내려가는 격자 구간을 로그 선형 보간.
# 격자 끝까지 LLOQ 이상이면 t_true_last = Inf. 반환 data.table(id, t_true_last, auc_true_last)
.m15_tlast_grid <- function(g, lloq) {
  L_ <- lloq
  g[, {
    i <- which.max(C)
    if (C[i] < L_) .(t_true_last = 0, auc_true_last = 0) else {
      post <- if (i < .N) (i + 1L):.N else integer(0); k2 <- post[which(C[post] < L_)[1]]
      if (is.na(k2)) .(t_true_last = Inf, auc_true_last = NA_real_) else {
        k1 <- k2 - 1L; kk <- log(C[k1] / C[k2]) / (time[k2] - time[k1])
        .(t_true_last = time[k1] + log(C[k1] / L_) / kk, auc_true_last = auc[k1] + (C[k1] - L_) / kk) }
    }
  }, by = id]
}

# 개체 파라미터 표(ip)의 참값(잔차 없음, R/oc.R와 같은 허용오차): AUCinf(400일, 잔여량이 투여량의 1e-7 초과면 4,000일), Cmax(격자 + 포물선),
# 일정별 조밀 채혈 AUC0-last 참값 AUClast_dense_<일정>(0부터 min(연속 시간 참 Tlast, 일정의 마지막 계획 시각)까지의 모델 AUC; 계획서 8.2(a)·8.4),
# 연속 시간 참 Tlast(t_true_last, 0.25일 격자(t_dense_max까지)에서 로그 선형 보간; t_dense_max 이후면 Inf),
# 부분 AUC 참값(고정 구간 경계에서 모델 AUC 차, Tlast와 무관).
m15_truth <- function(ip, dose_mg, model_id, schedules, lloq, intervals, chunk = 20000L, t_dense_max = 70) {
  mod <- get_model(model_id)
  bnd <- unique(unlist(intervals)); last_t <- vapply(schedules, function(x) max(as.numeric(x)), 1)
  stopifnot(all(last_t <= t_dense_max))
  times <- sort(unique(c(seq(0.25, t_dense_max, by = 0.25), unlist(schedules), bnd[bnd > 0], 400)))
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  res <- vector("list", ceiling(nrow(ip) / chunk))
  for (b in seq_along(res)) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, nrow(ip))
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    s <- .solve_param_table(mod, prm, times, dose_mg)
    inf <- s[time == 400, .(id, AUCinf = auc, resid = central + periph + depot + (if ("tr1" %in% names(s)) tr1 + tr2 + tr3 + absc else 0))]
    d_ <- dose_mg; long <- inf[resid / d_ > 1e-7, id]
    if (length(long)) {
      s2 <- .solve_param_table(mod, prm[prm$id %in% long, , drop = FALSE], 4000, dose_mg)
      inf[s2[, .(id, a2 = auc)], on = "id", AUCinf := i.a2]
    }
    g <- s[time <= 60, .(id, time, C)]; setorder(g, id, time)
    cm <- g[, { i <- which.max(C); v <- C[i]; if (i > 1 && i < .N) v <- .parabola_max(time[(i - 1):(i + 1)], C[(i - 1):(i + 1)]); .(Cmax = v) }, by = id]
    gd <- s[time <= t_dense_max, .(id, time, C, auc)]; setorder(gd, id, time)
    tl <- .m15_tlast_grid(gd, lloq)
    out <- merge(merge(inf[, .(id, AUCinf)], cm, by = "id"), tl, by = "id")
    for (sn in names(schedules)) {
      lt <- last_t[[sn]]
      a_end <- s[abs(time - lt) < 1e-9, .(id, a_end = auc)]
      out[a_end, on = "id", (paste0("AUClast_dense_", sn)) := fifelse(t_true_last >= lt, i.a_end, auc_true_last)]
    }
    for (nm in names(intervals)) {
      lo_ <- intervals[[nm]][1]; hi_ <- intervals[[nm]][2]
      a_lo <- if (lo_ == 0) s[time == min(time), .(id, a0 = 0)] else s[time == lo_, .(id, a0 = auc)]
      a_hi <- s[time == hi_, .(id, a1 = auc)]
      pa <- merge(a_lo, a_hi, by = "id")[, .(id, v = a1 - a0)]
      out[pa, on = "id", (paste0(nm, "_true")) := i.v]
    }
    out[, auc_true_last := NULL]
    res[[b]] <- out
  }
  rbindlist(res)[]
}

# 연속 시간 참 Tlast와 그때까지의 AUC: 0.1일 격자(최대 100일)에서 Cmax 뒤 참 농도가 LLOQ 아래로 처음 내려가는 구간을 로그 선형 보간.
# 반환 data.table(id, t_true_last, AUClast_true_cont, censored = 100일 안에 LLOQ 미만이 되지 않음)
m15_true_tlast_auc <- function(ip, dose_mg, model_id, lloq, step = 0.1, t_end = 100, chunk = 2000L) {
  mod <- get_model(model_id); d_ <- dose_mg; L_ <- lloq
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  times <- seq(step, t_end, by = step)
  res <- vector("list", ceiling(nrow(ip) / chunk))
  for (b in seq_along(res)) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, nrow(ip))
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    s <- .solve_param_table(mod, prm, times, d_)[, .(id, time, C, auc)]
    setorder(s, id, time)
    res[[b]] <- s[, {
      i <- which.max(C); post <- if (i < .N) (i + 1):.N else integer(0)
      k2 <- post[which(C[post] < L_)[1]]
      if (is.na(k2)) .(t_true_last = time[.N], AUClast_true_cont = auc[.N], censored = TRUE) else {
        k1 <- k2 - 1L; kk <- log(C[k1] / C[k2]) / (time[k2] - time[k1]); ts <- time[k1] + log(C[k1] / L_) / kk
        .(t_true_last = ts, AUClast_true_cont = auc[k1] + (C[k1] - L_) / kk, censored = FALSE) }
    }, by = id]
  }
  rbindlist(res)[]
}

# 시험 단위 짝 부트스트랩 AUROC 신뢰구간: pos, neg는 같은 시험 순서(길이 같음; 양성이 여러 판이면 행렬, 열 = 판). B회 시험 재표집.
m15_auroc_boot <- function(pos, neg, B = 2000L, seed = 1L, conf = 0.95) {
  pos <- as.matrix(pos); nt <- length(neg); stopifnot(nrow(pos) == nt)
  est <- m15_auroc(as.vector(pos), neg)$auroc
  bs <- with_seed(seed, vapply(seq_len(B), function(b) { i <- sample.int(nt, nt, replace = TRUE); m15_auroc(as.vector(pos[i, , drop = FALSE]), neg[i])$auroc }, 1))
  a <- (1 - conf) / 2
  list(auroc = est, lo = quantile(bs, a, names = FALSE, na.rm = TRUE), hi = quantile(bs, 1 - a, names = FALSE, na.rm = TRUE), se_boot = sd(bs, na.rm = TRUE))
}

# 부트스트랩 재표집 색인(시험 단위, 모든 조건이 같은 색인을 공유 → 조건 간 차이의 부트스트랩도 같은 재표집)
m15_boot_index <- function(nt, B = 2000L, seed = 1L) with_seed(seed, lapply(seq_len(B), function(b) sample.int(nt, nt, replace = TRUE)))
# AUROC 부트스트랩 값: pos 행렬(시험 x 판), neg 벡터(시험), idx = m15_boot_index(). 시험 j의 양성·음성은 함께 재표집된다.
m15_auroc_draws <- function(pos, neg, idx) {
  pos <- as.matrix(pos); stopifnot(nrow(pos) == length(neg))
  vapply(idx, function(i) .m15_mw(as.vector(pos[i, , drop = FALSE]), neg[i]), 1)
}
# 같은 값의 빠른 계산: 재표집(복원 추출)은 시험별 가중치(뽑힌 횟수)와 같으므로, 정렬한 음성 값의 누적 가중치로 가중 Mann-Whitney를 구한다.
# pos 행렬(시험 x 판), neg 벡터(시험). 결측이 있는 시험은 미리 빼서 넣는다(호출자). m15_auroc_draws와 같은 값(단위 시험).
m15_auroc_draws_w <- function(pos, neg, idx) {
  pos <- as.matrix(pos); nt <- length(neg); stopifnot(nrow(pos) == nt, all(is.finite(pos)), all(is.finite(neg)))
  o <- order(neg); sn <- neg[o]; k <- ncol(pos)
  pv <- as.vector(pos); pt <- rep(seq_len(nt), k)                     # 양성 값과 그 시험
  lt <- findInterval(pv, sn, left.open = TRUE); le <- findInterval(pv, sn)   # 음성 중 미만·이하 개수(정렬 순서 위치)
  vapply(idx, function(i) {
    w <- tabulate(i, nt); wn <- w[o]; C <- c(0, cumsum(wn))
    below <- C[lt + 1L]; eq <- C[le + 1L] - C[lt + 1L]
    sum(w[pt] * (below + 0.5 * eq)) / (sum(w[pt]) * sum(wn))
  }, 1)
}

# Mann-Whitney AUROC(순위 기반, 동점 1/2)
.m15_mw <- function(pos, neg) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]; m <- length(pos); n <- length(neg)
  if (m < 1 || n < 1) return(NA_real_)
  r <- rank(c(pos, neg)); (sum(r[seq_len(m)]) - m * (m + 1) / 2) / (m * n)
}

# 짝지은 이항 비율 차이(같은 시험에서 두 설계의 판정)의 Tango 점수 신뢰구간(Tango 1998, Stat Med 17:891).
# b = 설계 A만 통과한 시험 수, c = 설계 B만 통과한 시험 수, n = 시험 수. 차이 = (b - c) / n = P(A) - P(B).
m15_tango_ci <- function(b, c, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2); d_hat <- (b - c) / n
  score <- function(D) {
    W <- -b - c + (2 * n - b + c) * D
    q21 <- (sqrt(pmax(W^2 + 8 * n * c * D * (1 - D), 0)) - W) / (4 * n)
    v <- n * (2 * q21 + D * (1 - D))
    (b - c - n * D) / sqrt(pmax(v, 1e-300))
  }
  eps <- 1e-12
  lo <- if (score(-1 + eps) <= z) -1 else uniroot(function(D) score(D) - z, c(-1 + eps, d_hat), tol = 1e-12)$root
  hi <- if (score(1 - eps) >= -z) 1 else uniroot(function(D) score(D) + z, c(d_hat, 1 - eps), tol = 1e-12)$root
  c(diff = d_hat, lo = lo, hi = hi, b = b, c = c, n = n)
}

# λz 외삽 부분 AUC(민감도): Tlast < 구간 끝 b이고 λz가 있으면 Clast·exp(-λz (t - Tlast))를 max(a, Tlast)부터 b까지 더한다. 없으면 원래 값.
m15_pauc_ext <- function(pauc, Tlast, Clast, lambda_z, a, b) {
  ok <- is.finite(pauc) & is.finite(Tlast) & is.finite(Clast) & is.finite(lambda_z) & lambda_z > 0 & Tlast < b
  st <- pmax(a, Tlast)
  add <- ifelse(ok, Clast / lambda_z * (exp(-lambda_z * (st - Tlast)) - exp(-lambda_z * (b - Tlast))), 0)
  pauc + add
}

# 짝지은 모집단 기하평균비(같은 대상자): exp(mean(log(X_T / X_R)))와 log 척도 MC 표준오차
m15_paired_gmr <- function(xT, xR) {
  ok <- is.finite(xT) & is.finite(xR) & xT > 0 & xR > 0
  la <- log(xT[ok] / xR[ok])
  c(gmr = exp(mean(la)), se_log = sd(la) / sqrt(length(la)), n = length(la))
}

# 전형 대상자(η = 0)의 AUC0-inf, Cmax(촘촘한 격자 0.01일, 포물선 보정), 조밀 채혈 AUC0-last 참값(0부터 min(참 Tlast, last_day)까지) — [8] 백분위수 비와 정규화 민감도 계수
m15_typical_metrics <- function(p, WT, dose_mg, mult = list(), step = 0.01, lloq = study_lloq(), last_day = 56) {
  wts_ <- WT; d_ <- dose_mg                                          # data.table [ 안에서는 지역 변수로(D-025 lint)
  pp <- apply_multipliers(p, mult)
  ip <- individual_params(pp, data.table(id = seq_along(wts_), WT = wts_, ada = 0))
  mod <- get_model(pp$model_id)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  prm <- as.data.frame(ip[, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
  s <- .solve_param_table(mod, prm, c(seq(step, 60, by = step), 400), d_)
  inf <- s[time == 400, .(id, AUCinf = auc, resid = central + periph + depot + (if ("tr1" %in% names(s)) tr1 + tr2 + tr3 + absc else 0))]
  long <- inf[resid / d_ > 1e-7, id]
  if (length(long)) { s2 <- .solve_param_table(mod, prm[prm$id %in% long, , drop = FALSE], 4000, d_); inf[s2[, .(id, a2 = auc)], on = "id", AUCinf := i.a2] }
  g <- s[time <= 60, .(id, time, C)]; setorder(g, id, time)
  cm <- g[, { i <- which.max(C); v <- C[i]; if (i > 1 && i < .N) v <- .parabola_max(time[(i - 1):(i + 1)], C[(i - 1):(i + 1)]); .(Cmax = v) }, by = id]
  ld_ <- last_day
  tl <- .m15_tlast_grid(s[time <= 60, .(id, time, C, auc)][order(id, time)], lloq)
  ae <- s[abs(time - ld_) < 1e-6, .(id, a_end = auc)]
  tl <- merge(tl, ae, by = "id")[, .(id, AUClast_dense = fifelse(t_true_last >= ld_, a_end, auc_true_last))]
  out <- merge(merge(inf[, .(id, AUCinf)], cm, by = "id"), tl, by = "id")
  set(out, j = "WT", value = wts_[out$id]); out[]
}

# 정규화 민감도 계수 (ln X(1.01 θ) - ln X(0.99 θ)) / (ln 1.01 - ln 0.99)
m15_sensitivity <- function(p, WT, dose_mg, params = c("F", "ke", "Vmax"), h = 0.01) {
  rbindlist(lapply(params, function(nm) {
    up <- m15_typical_metrics(p, WT, dose_mg, setNames(list(1 + h), nm)); dn <- m15_typical_metrics(p, WT, dose_mg, setNames(list(1 - h), nm))
    den <- log(1 + h) - log(1 - h)
    data.table(param = nm, WT = up$WT, S_AUCinf = (log(up$AUCinf) - log(dn$AUCinf)) / den, S_Cmax = (log(up$Cmax) - log(dn$Cmax)) / den,
               S_AUClast_dense = (log(up$AUClast_dense) - log(dn$AUClast_dense)) / den)
  }))
}

# ---------------------------------------------------------------- 시드·격자·체중 모멘트 ----------------------------------------------------------------
# 시험 시드 표(계획서 7.6절): derive_seed(ms, "trial", j, 흐름). 28비트 시드 공간의 충돌은 (j, 흐름) 순서상 뒤의 것에 "r1"(다음 "r2" ...)을 붙여 다시 파생.
# 표 크기는 등록된 최대 시험 수로 고정(모드와 무관하게 같은 시드). 반환 data.table(j, stream, seed, resolved)
m15_trial_seeds <- function(ms, n_trials, streams = c("subj", "alloc", "drop", "jitter", "eps")) {
  tab <- CJ(j = seq_len(n_trials), stream = streams, sorted = FALSE)
  tab[, `:=`(seed = mapply(function(a, b) derive_seed(ms, "trial", a, b), j, stream), resolved = 0L)]
  repeat {
    dup <- duplicated(tab$seed); if (!any(dup)) break
    i <- which(dup); set(tab, i, "resolved", tab$resolved[i] + 1L)
    set(tab, i, "seed", mapply(function(a, b, r) derive_seed(ms, "trial", a, b, paste0("r", r)), tab$j[i], tab$stream[i], tab$resolved[i]))
  }
  tab[]
}

# 시드 표 → 행렬(행 = 시험 j, 열 = 흐름) 조회. data.table 안에서 인자 이름이 열 이름(j, stream)으로 해석되는 일을 피하려고 행렬로 바꾼다
# (scripts/73의 SEEDS[.(j, s), seed]는 j를 열로 읽어 모든 시험에 같은 시드를 주었다; D-068 이탈 기록).
m15_seed_matrix <- function(tab) {
  w <- dcast(tab, j ~ stream, value.var = "seed"); stopifnot(identical(w$j, seq_len(nrow(w))))
  m <- as.matrix(w[, setdiff(names(w), "j"), with = FALSE]); storage.mode(m) <- "integer"; m
}
m15_seed_of <- function(sm, trial_j, stream_s) { x <- unname(sm[trial_j, stream_s]); stopifnot(length(x) == 1L, is.finite(x)); x }

# 참 GMR 격자: log 대칭 n점(lo, hi 포함; n 홀수면 1 포함) + 추가 점
m15_gmr_grid <- function(spec) {
  x <- exp(seq(log(as.numeric(spec$lo)), log(as.numeric(spec$hi)), length.out = as.integer(spec$n_points)))
  if (!is.null(spec$extra)) x <- c(x, as.numeric(unlist(spec$extra)))
  x <- sort(unique(signif(x, 12))); x[abs(x - 1) < 1e-12] <- 1; x
}

# 경계별 log 체중 모멘트(R/samplesize.R ss_p2_power 입력): wt = 등록 체중 표본, bk = 경계(kg)
m15_weight_moments <- function(wt, bk) {
  lw <- log(wt); st <- wt > bk; p2 <- mean(st); m1 <- mean(lw[!st]); m2 <- mean(lw[st]); ma <- mean(lw)
  vb <- (1 - p2) * (m1 - ma)^2 + p2 * (m2 - ma)^2
  vw <- (1 - p2) * var(lw[!st]) * (sum(!st) - 1) / sum(!st) + p2 * var(lw[st]) * (sum(st) - 1) / sum(st)
  data.table(p_stratum2 = p2, mean_lwt = ma, var_between = vb, var_within = vw, var_total = vb + vw)
}

# ---------------------------------------------------------------- 요약 보조 ----------------------------------------------------------------
m15_wilson <- function(x, n, conf = 0.95) { w <- wilson_ci(x, n, conf); list(pct = w$est, lo = w$lo, hi = w$hi) }

# 조건부 표시(계획서 6.2절 Action): 주 모델이 고정 체중 가정에서 재현 실패 규칙을 충족했으면 이후 결과 파일마다 같은 문구 열을 붙인다
m15_conditional_label <- function(rd, primary) {
  f <- file.path(rd, "validation", "validation_rule.csv"); if (!file.exists(f)) return("")
  r <- data.table::fread(f); r <- r[r$pk_model == primary & r$wmode == "primary"]
  if (isTRUE(r$failure[1])) "Conditional on a primary model that did not pass the pre-specified external check (plan section 6.2)." else ""
}
m15_fwrite_cond <- function(tab, path, cond) {                       # 호출자 표는 바꾸지 않는다(복사본에 열 추가)
  tab <- data.table::copy(data.table::as.data.table(tab)); if (nzchar(cond)) tab[, conditionality := cond]; data.table::fwrite(tab, path)
}
