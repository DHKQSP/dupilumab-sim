# oc.R — 운용특성(operating characteristics) 분석 (검토 의견 2026-09-24 §3–§4, config/oc_design.yaml, D-040).
# (1) 참값: 같은 가상 대상자(공통 난수)에 대조·시험 파라미터를 적용해 모집단 기하평균비(AUC0-inf, Cmax; 잔차 없음)
# (2) 역산: 기전·방향별 로그 배율 이분법으로 목표 참값 비를 주는 배율
# (3) 시험: 대조군 NCA 1회 + 시험군 시나리오 전부를 한 번에 NCA(id 오프셋), 평가변수별 pooled t, 구성별 통과

OC_ENDPOINTS <- c("Cmax", "AUClast", "AUCinf_A", "AUCinf_B", "AUCinf_C", "AUCinf_true")

# --- (1) 참값 -------------------------------------------------------------------------------
# 공통 난수 대상자(체중·성별·IIV). 200,000명의 앞 n명이 선별용과 같다(같은 시드, 같은 순서).
truth_subjects <- function(n, p, design, seed, weight = "base") {
  wt <- weight_spec_from_design(design, weight)
  with_seed(seed, make_subjects(n, p, wt, design$weight$sex_ratio_male$value, 0))
}

# 세 점을 지나는 포물선의 최댓값(간격이 달라도 됨). 꼭짓점이 세 점 구간 밖이거나 위로 볼록하지 않으면 격자 최댓값
.parabola_max <- function(x, y) {
  d1 <- (y[2] - y[1]) / (x[2] - x[1]); d2 <- (y[3] - y[2]) / (x[3] - x[2]); a <- (d2 - d1) / (x[3] - x[1])
  if (!is.finite(a) || a >= 0) return(max(y))
  b <- d1 - a * (x[1] + x[2]); xs <- -b / (2 * a)
  if (xs < x[1] || xs > x[3]) return(max(y))
  y[1] + d1 * (xs - x[1]) + a * (xs - x[1]) * (xs - x[2])
}
TRUTH_ATOL <- 1e-8; TRUTH_RTOL <- 1e-6   # 참값 전용. 1e-10/1e-8 대비 개인별 최대 상대 차이 5e-6, 평균 log 차이 5e-7(두 모델, 극단 배율 포함) — D-042
.truth_grid <- function(cmax = TRUE) if (cmax) c(seq(0.25, 30, by = 0.25), seq(30.5, 60, by = 0.5), 400) else 400

.solve_param_table <- function(mod, prm, times, dose_mg) {
  ev <- rxode2::add.sampling(rxode2::et(amt = dose_mg, cmt = "depot", time = 0), times)
  s <- rxode2::rxSolve(mod, params = prm, events = ev, atol = TRUTH_ATOL, rtol = TRUTH_RTOL, maxsteps = 2000000L, cores = 1L, returnType = "data.table", addDosing = FALSE)
  if ("sim.id" %in% names(s)) s[, id := prm$id[sim.id]] else if ("id" %in% names(s)) s[, id := prm$id[as.integer(id)]] else if (nrow(prm) == 1) s[, id := prm$id[1]] else stop("rxSolve 출력에 id가 없습니다")
  s
}

# 빠른 풀이: 파라미터 표 + 단일 사건표(모든 대상자 같은 용량). 반환 data.table(id, AUCinf[, Cmax, tmax_grid])
# 400일 잔여량/투여량 > 1e-7인 대상자는 4,000일까지 다시 적분(느린 흡수 배율에서 AUC0-inf 절단 방지).
truth_metrics <- function(ip, dose_mg, model_id, cmax = TRUE, chunk = 20000L) {
  mod <- get_model(model_id)
  times <- .truth_grid(cmax)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  out <- vector("list", ceiling(nrow(ip) / chunk))
  for (b in seq_along(out)) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, nrow(ip))
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    s <- .solve_param_table(mod, prm, times, dose_mg)
    inf <- s[time == 400, .(id, AUCinf = auc, resid = central + periph + depot + (if ("tr1" %in% names(s)) tr1 + tr2 + tr3 + absc else 0))]
    long <- inf[resid / dose_mg > 1e-7, id]
    if (length(long)) {
      s2 <- .solve_param_table(mod, prm[prm$id %in% long, , drop = FALSE], 4000, dose_mg)
      inf[s2[, .(id, a2 = auc)], on = "id", AUCinf := i.a2]
    }
    if (cmax) {
      g <- s[time < 400, .(id, time, C)]
      setorder(g, id, time)
      cm <- g[, {
        i <- which.max(C); v <- C[i]
        if (i > 1 && i < .N) v <- .parabola_max(time[(i - 1):(i + 1)], C[(i - 1):(i + 1)])
        .(Cmax = v, tmax_grid = time[i], at_edge = i == .N) }, by = id]
      if (any(cm$at_edge)) warning(sprintf("참값 Cmax: %d명이 격자 끝(60일)에서 최대", sum(cm$at_edge)))
      inf <- merge(inf, cm[, .(id, Cmax, tmax_grid)], by = "id")
    }
    out[[b]] <- inf[, resid := NULL]
  }
  rbindlist(out)[]
}

# 짝지은 모집단 기하평균비와 MC 표준오차
truth_ratio <- function(ref, test) {
  m <- merge(ref, test, by = "id", suffixes = c(".R", ".T"))
  la <- log(m$AUCinf.T / m$AUCinf.R)
  out <- list(auc_ratio = exp(mean(la)), auc_se_log = sd(la) / sqrt(length(la)), n = length(la))
  if ("Cmax.R" %in% names(m)) { lc <- log(m$Cmax.T / m$Cmax.R); out$cmax_ratio <- exp(mean(lc)); out$cmax_se_log <- sd(lc) / sqrt(length(lc)) }
  out
}

# 시험 파라미터로 개체 파라미터 표(대상자 순서·id 동일)
ip_for_mult <- function(p, subj, mult) individual_params(apply_multipliers(p, mult), subj)

mech_multiplier <- function(mech, m) setNames(list(m), mech)

# 기전의 방향별 탐색 구간(F는 시험군 F <= 1로 절단)
mech_range <- function(oc, mech, p) {
  r <- as.numeric(unlist(oc$mechanisms[[mech]]$range))
  if (mech == "F") r[2] <- min(r[2], (1 - 1e-9) / p$theta[["F"]])
  r
}

# --- (2) 역산 -------------------------------------------------------------------------------
# 한 기전·방향에서 목표 참값 비(AUC0-inf)를 주는 배율을 로그 배율 이분법으로 찾는다.
#   eval_fn(m) -> 참 AUC 비(선별 대상자). bracket: c(m_lo, m_hi) with ratio(m_lo), ratio(m_hi)가 목표를 사이에 둠.
bisect_log <- function(eval_fn, target, bracket, r_bracket, tol_rel, max_iter = 60) {
  lo <- log(bracket[1]); hi <- log(bracket[2]); rlo <- r_bracket[1]; rhi <- r_bracket[2]
  stopifnot((rlo - target) * (rhi - target) <= 0)
  it <- 0; m <- NA_real_; r <- NA_real_
  repeat {
    it <- it + 1
    mid <- (lo + hi) / 2; m <- exp(mid); r <- eval_fn(m)
    if (abs(r / target - 1) <= tol_rel || it >= max_iter) break
    if ((rlo - target) * (r - target) <= 0) { hi <- mid; rhi <- r } else { lo <- mid; rlo <- r }
  }
  list(m = m, ratio = r, iter = it, bracket = exp(c(lo, hi)))
}

# --- (3) 시험 -------------------------------------------------------------------------------
# 대조군 NCA 1회, 시험군 시나리오 전부를 id 오프셋으로 한 번의 NCA. 반환: be(trial, scenario, endpoint, GMR, CI, pass, n_R, n_T), drop(시나리오·군별 탈락 수)
run_trial_oc <- function(j, p, design, scenarios, master_seed, wt_spec, schedule = "B0", jitter = TRUE, model_id = NULL) {
  sa <- simulate_trial_arms(j, p, design, scenarios, schedule, master_seed, wt_spec, jitter = jitter, model_id = model_id)
  sd_days <- get_schedule(design, schedule)
  # 대조군(id 1..nmax)과 시험군 시나리오 k(id + k·nmax)를 한 번의 NCA로. 개인 AUCinf_true는 키 조인(진단용 평가변수)
  scs <- names(scenarios); nmax <- max(sa$subj$id)
  R <- sa$arms$R$S00
  obs_l <- list(subset_schedule(R$obs, sd_days)[, .(id, time, conc)]); tr_l <- list(R$truth[, .(id, AUCinf_true)])
  map <- list(data.table(id = R$subj$id, scenario = "REF"))
  for (k in seq_along(scs)) {
    x <- sa$arms$T[[scs[k]]]; o <- k * nmax
    obs_l[[k + 1]] <- subset_schedule(x$obs, sd_days)[, .(id = id + o, time, conc)]
    tr_l[[k + 1]] <- x$truth[, .(id = id + o, AUCinf_true)]
    map[[k + 1]] <- data.table(id = x$subj$id + o, scenario = scs[k])
  }
  nca <- run_nca(rbindlist(obs_l))
  nca <- rbindlist(tr_l)[nca, on = "id"]
  nca <- rbindlist(map)[nca, on = "id"]
  ncaR <- nca[scenario == "REF"]; ncaT <- nca[scenario != "REF"]
  addC <- function(n) n[, AUCinf_C := fifelse(reliable %in% TRUE, AUCinf, AUClast)][]
  ncaR <- addC(ncaR); ncaT <- addC(ncaT)
  endpoint_vals <- function(n, ep) switch(ep,
    Cmax = n$Cmax, AUClast = n$AUClast, AUCinf_A = fifelse(n$reliable %in% TRUE, n$AUCinf, NA_real_),
    AUCinf_B = fifelse(n$lambda_ok %in% TRUE, n$AUCinf, NA_real_), AUCinf_C = n$AUCinf_C, AUCinf_true = n$AUCinf_true)
  lim <- as.numeric(unlist(design$be$limits)); cl <- design$be$ci_level
  be <- rbindlist(lapply(scs, function(sc) {
    t_ <- ncaT[scenario == sc]
    rbindlist(lapply(OC_ENDPOINTS, function(ep) {
      yR <- endpoint_vals(ncaR, ep); yT <- endpoint_vals(t_, ep)
      r <- be_pooled_t(c(yR, yT), c(rep("R", length(yR)), rep("T", length(yT))), cl, lim)
      data.table(trial = j, scenario = sc, endpoint = ep, GMR = r$GMR, CI_lower = r$CI_lower, CI_upper = r$CI_upper, pass = r$pass, n_R = r$n_R, n_T = r$n_T)
    }))
  }))
  drop <- rbind(data.table(trial = j, scenario = "REF", arm = "R", n_reliable = sum(ncaR$reliable %in% TRUE), n_lambda = sum(ncaR$lambda_ok %in% TRUE),
                           flag_rsq = sum(ncaR$flag_rsq %in% TRUE), flag_extrap = sum(ncaR$flag_extrap %in% TRUE), flag_span = sum(ncaR$flag_span %in% TRUE), n = nrow(ncaR)),
                ncaT[, .(trial = j, arm = "T", n_reliable = sum(reliable %in% TRUE), n_lambda = sum(lambda_ok %in% TRUE), flag_rsq = sum(flag_rsq %in% TRUE),
                         flag_extrap = sum(flag_extrap %in% TRUE), flag_span = sum(flag_span %in% TRUE), n = .N), by = scenario], use.names = TRUE)
  list(be = be, drop = drop)
}

run_trials_oc <- function(trial_ids, p, design, scenarios, master_seed, wt_spec, cores = 1L, model_id = NULL, progress_every = 250) {
  invisible(get_model(if (is.null(model_id)) p$model_id else model_id))
  one <- function(j) {
    r <- run_trial_oc(j, p, design, scenarios, master_seed, wt_spec, model_id = model_id)
    if (progress_every > 0 && j %% progress_every == 0) cat(sprintf("  trial %d %s\n", j, format(Sys.time(), "%H:%M:%S")))
    r
  }
  res <- if (cores > 1) parallel::mclapply(trial_ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(trial_ids, one)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
  if (any(bad)) stop("실패한 시험 반복: ", paste(trial_ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
  list(be = rbindlist(lapply(res, `[[`, "be")), drop = rbindlist(lapply(res, `[[`, "drop")))
}

# 구성별 통과: be(long) -> wide(trial, scenario, pass_<endpoint>) -> 구성 pass
config_pass <- function(be, oc) {
  w <- dcast(be, trial + scenario ~ endpoint, value.var = "pass")
  for (cf in names(oc$configurations)) {
    eps <- unlist(oc$configurations[[cf]]$endpoints)
    w[, (paste0("cfg_", cf)) := Reduce(`&`, lapply(eps, function(e) w[[e]] %in% TRUE))]
  }
  w[]
}
