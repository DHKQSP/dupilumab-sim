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
    d_ <- dose_mg; long <- inf[resid / d_ > 1e-7, id]   # 인자는 지역 변수로(D-025 lint)
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

# --- (4) 재판정 확장: AUCinf 규칙 × 플래그 세트 (검토 의견 W2 §2) -------------------------------------
# 플래그 세트 (i):  λz 산출 가능 & Rsq_adjusted ≥ 0.80 & AUC_%Extrap_obs ≤ 20% (span ratio 무시; NCA 출력 flag_rsq·flag_extrap 열)
# 플래그 세트 (ii): (i) & Span_ratio ≥ 2 (= reliable 열, D-039). 원래 AUCinf_A·AUCinf_C는 (ii), AUCinf_B는 플래그와 무관.
# 추가 평가변수: AUCinf_Ai(규칙 A, 세트 (i) 충족자만), AUCinf_Ci(규칙 C, 세트 (i) 미충족자는 AUClast 대입)
OC_ENDPOINTS_EXT <- c(OC_ENDPOINTS, "AUCinf_Ai", "AUCinf_Ci")

# run_trial_oc와 같은 simulate_trial_arms 호출·시드·NCA에 평가변수 2개를 더한다. 원래 6개 평가변수는 run_trial_oc와 비트 단위로 같다
# (tests/testthat/test-rejudge.R, scripts/40_oc_rejudge.R가 저장본과 대조). drop에는 세트 (i) 충족 수(n_reliable_i)를 더한다.
run_trial_oc_ext <- function(j, p, design, scenarios, master_seed, wt_spec, schedule = "B0", jitter = TRUE, model_id = NULL) {
  sa <- simulate_trial_arms(j, p, design, scenarios, schedule, master_seed, wt_spec, jitter = jitter, model_id = model_id)
  sd_days <- get_schedule(design, schedule)
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
  addC <- function(n) {
    n[, AUCinf_C := fifelse(reliable %in% TRUE, AUCinf, AUClast)]
    n[, rel_i := (lambda_ok & !flag_rsq & !flag_extrap) %in% TRUE]      # 세트 (i). 플래그가 NA면 미충족(세트 (ii)의 reliable %in% TRUE와 같은 처리)
    n[, AUCinf_Ci := fifelse(rel_i, AUCinf, AUClast)][]
  }
  ncaR <- addC(ncaR); ncaT <- addC(ncaT)
  endpoint_vals <- function(n, ep) switch(ep,
    Cmax = n$Cmax, AUClast = n$AUClast, AUCinf_A = fifelse(n$reliable %in% TRUE, n$AUCinf, NA_real_),
    AUCinf_B = fifelse(n$lambda_ok %in% TRUE, n$AUCinf, NA_real_), AUCinf_C = n$AUCinf_C, AUCinf_true = n$AUCinf_true,
    AUCinf_Ai = fifelse(n$rel_i, n$AUCinf, NA_real_), AUCinf_Ci = n$AUCinf_Ci)
  lim <- as.numeric(unlist(design$be$limits)); cl <- design$be$ci_level
  be <- rbindlist(lapply(scs, function(sc) {
    t_ <- ncaT[scenario == sc]
    rbindlist(lapply(OC_ENDPOINTS_EXT, function(ep) {
      yR <- endpoint_vals(ncaR, ep); yT <- endpoint_vals(t_, ep)
      r <- be_pooled_t(c(yR, yT), c(rep("R", length(yR)), rep("T", length(yT))), cl, lim)
      data.table(trial = j, scenario = sc, endpoint = ep, GMR = r$GMR, CI_lower = r$CI_lower, CI_upper = r$CI_upper, pass = r$pass, n_R = r$n_R, n_T = r$n_T)
    }))
  }))
  drop <- rbind(data.table(trial = j, scenario = "REF", arm = "R", n_reliable = sum(ncaR$reliable %in% TRUE), n_reliable_i = sum(ncaR$rel_i), n_lambda = sum(ncaR$lambda_ok %in% TRUE),
                           flag_rsq = sum(ncaR$flag_rsq %in% TRUE), flag_extrap = sum(ncaR$flag_extrap %in% TRUE), flag_span = sum(ncaR$flag_span %in% TRUE), n = nrow(ncaR)),
                ncaT[, .(trial = j, arm = "T", n_reliable = sum(reliable %in% TRUE), n_reliable_i = sum(rel_i), n_lambda = sum(lambda_ok %in% TRUE), flag_rsq = sum(flag_rsq %in% TRUE),
                         flag_extrap = sum(flag_extrap %in% TRUE), flag_span = sum(flag_span %in% TRUE), n = .N), by = scenario], use.names = TRUE)
  list(be = be, drop = drop)
}

# 재판정 구성(규칙 × 플래그 세트). pass = 구성 안의 모든 평가변수 통과(평가변수 pass가 NA면 불통과, config_pass와 같음).
# 규칙 B는 플래그와 무관하므로 (i)·(ii) 구분 없이 하나. P2·G2_Aii·AUCinf_Aii는 사전 고정 구성 P2·G2·AUCinf_only와 같은 정의.
OC_REJUDGE_CONFIGS <- list(
  P2 = c("AUClast", "Cmax"),
  G2_Aii = c("AUCinf_A", "Cmax"), G2_B = c("AUCinf_B", "Cmax"), G2_Cii = c("AUCinf_C", "Cmax"),
  G2_Ai = c("AUCinf_Ai", "Cmax"), G2_Ci = c("AUCinf_Ci", "Cmax"),
  AUCinf_Aii = "AUCinf_A", AUCinf_B = "AUCinf_B", AUCinf_Cii = "AUCinf_C", AUCinf_Ai = "AUCinf_Ai", AUCinf_Ci = "AUCinf_Ci")

# be(long) -> wide(trial, scenario, 평가변수별 0/1, cfg_<구성>). 평가변수 행 자체가 없는 (시험, 시나리오)의 구성은 NA(미평가)로 두어
# 판정 불가(pass NA = 불통과)와 구분한다(예: 세트 (i) 평가변수가 재판정 파일에만 있을 때).
config_pass_ext <- function(be, configs = OC_REJUDGE_CONFIGS) {
  b <- be[, .(trial, scenario, endpoint, ok = as.integer(pass %in% TRUE))]
  if (anyDuplicated(b, by = c("trial", "scenario", "endpoint"))) stop("config_pass_ext: (trial, scenario, endpoint) 중복")
  w <- dcast(b, trial + scenario ~ endpoint, value.var = "ok", fill = NA_integer_)
  for (cf in names(configs)) {
    eps <- configs[[cf]]
    v <- rep(NA, nrow(w))
    if (all(eps %in% names(w))) {
      miss <- Reduce(`|`, lapply(eps, function(e) is.na(w[[e]])))
      v <- Reduce(`&`, lapply(eps, function(e) w[[e]] %in% 1L)); v[miss] <- NA
    }
    w[, (paste0("cfg_", cf)) := v]
  }
  w[]
}

# --- (5) 적응적 연장 (지시 2026-09-25 §3) ----------------------------------------------------------------------------
# 규칙(scripts/20_products_5000.R의 5,000회 제품 시나리오와 같은 적응적 상향을 모든 모델에 일반 적용): 사전 고정 reps_boundary(10,000)회의
# 경계 1종 오류(구성 P2) Wilson 95% 구간이 임계 5%를 포함하면(lo ≤ 5 ≤ hi, 경계 포함) 그 (모델, 경계 시나리오)만 n_to(20,000)회로 늘린다.
# 연장분은 별도 파일 oc_trials_ext_be_<model>.csv.gz(scripts/42)에 두어 사전 고정 10,000회 파일을 바꾸지 않는다.
oc_extension_rule_text <- function(n_from, n_to, cfg = "P2", threshold = 5)
  sprintf("%s boundary type I error Wilson 95%% CI at %d trials includes %g%% (lo <= %g <= hi) -> extend that boundary scenario to %d trials",
          cfg, as.integer(n_from), threshold, threshold, as.integer(n_to))

# bt: boundary_type1.csv 형식(model, scenario, config, n_trials, pass_pct, lo, hi[, mechanism, direction, target, multiplier]).
# 반환: 구성 cfg의 모델 × 경계 시나리오별 한 행(n_before, pass_pct·lo·hi = n_from회 값, threshold, rule, selected, n_after = 목표 시험 수).
# 모든 행의 n_trials가 n_from이어야 한다(규칙은 사전 고정 시험 수에서 판정한다. 연장 뒤의 표로 다시 판정하지 않는다).
oc_extension_select <- function(bt, n_from, n_to, cfg = "P2", threshold = 5) {
  need <- c("model", "scenario", "config", "n_trials", "pass_pct", "lo", "hi")
  if (!all(need %in% names(bt))) stop("oc_extension_select: 경계 표에 열이 없습니다: ", paste(setdiff(need, names(bt)), collapse = ", "))
  n_from <- as.integer(n_from); n_to <- as.integer(n_to)
  if (is.na(n_to) || n_to <= n_from) stop(sprintf("oc_extension_select: 목표 시험 수(%s)가 %d보다 커야 합니다", n_to, n_from))
  x <- as.data.table(bt)[bt$config == cfg]
  if (!nrow(x)) stop(sprintf("oc_extension_select: 구성 %s 행이 없습니다", cfg))
  if (anyDuplicated(x, by = c("model", "scenario"))) stop("oc_extension_select: (model, scenario) 중복")
  if (anyNA(x$lo) || anyNA(x$hi) || anyNA(x$pass_pct)) stop("oc_extension_select: 통과율·구간에 NA")
  if (any(x$n_trials != n_from)) stop(sprintf("oc_extension_select: 판정은 사전 고정 %d회 값으로 한다(n_trials가 다른 행 %d개)", n_from, sum(x$n_trials != n_from)))
  keep <- intersect(c("model", "scenario", "mechanism", "direction", "target", "multiplier"), names(x))
  out <- x[, c(keep, "n_trials", "pass_pct", "lo", "hi"), with = FALSE]
  setnames(out, "n_trials", "n_before")
  sel <- out$lo <= threshold & out$hi >= threshold
  out[, `:=`(config = cfg, threshold = threshold, rule = oc_extension_rule_text(n_from, n_to, cfg, threshold), selected = sel)]
  out[, n_after := ifelse(sel, n_to, n_from)]
  setcolorder(out, c(keep, "config", "n_before", "pass_pct", "lo", "hi", "threshold", "rule", "selected", "n_after"))
  out[order(model, scenario)][]
}

# 시험 수 문구(그림 캡션·자동 문구·보고서 공용): 시나리오 집합의 최빈 시험 수를 "N each"(ko "각 N회")로, 그와 다른 시나리오는 실제 수로 따로 적는다.
# x: (model, scenario, mechanism, direction, target, n_trials) 행(시나리오당 여러 행이면 같은 n_trials여야 함). dec: extension_decision.csv(scripts/42) 또는 NULL.
# 연장 선택되어 사전 고정 수보다 많은 시나리오는 "adaptive extension rule"(목표 미달이면 진행 중과 목표)을 붙인다. model_labels가 있으면 예외 항목 앞에 모델 이름.
# exceptions_only = TRUE: 예외 항목만("; "로 연결, 없으면 "").
oc_n_trials_text <- function(x, lang = c("en", "ko"), dec = NULL, model_labels = NULL, exceptions_only = FALSE) {
  lang <- match.arg(lang)
  u <- unique(as.data.table(x)[, c("model", "scenario", "mechanism", "direction", "target", "n_trials"), with = FALSE])
  if (!nrow(u)) stop("oc_n_trials_text: 행이 없습니다")
  if (anyDuplicated(u, by = c("model", "scenario"))) stop("oc_n_trials_text: 한 (모델, 시나리오)에 시험 수가 둘 이상")
  fmt <- function(n) format(as.integer(n), big.mark = ",", trim = TRUE)
  tb <- u[, .N, by = "n_trials"][order(-N, n_trials)]; n0 <- tb$n_trials[1]
  base <- if (lang == "ko") sprintf("각 %s회", fmt(n0)) else sprintf("%s each", fmt(n0))
  ex <- u[u$n_trials != n0][order(model, target, mechanism, direction)]
  if (!nrow(ex)) return(if (exceptions_only) "" else base)
  dir_l <- if (lang == "ko") c(down = "하향", up = "상향") else c(down = "down", up = "up")
  items <- vapply(seq_len(nrow(ex)), function(i) {
    r <- ex[i]
    lab <- if (r$scenario == "S00") "S00" else sprintf("%s %s %.2f", r$mechanism, if (r$direction %in% names(dir_l)) dir_l[[r$direction]] else r$direction, r$target)
    if (!is.null(model_labels)) lab <- paste(model_labels[[r$model]], lab)
    d <- if (is.null(dec)) NULL else dec[dec$model == r$model & dec$scenario == r$scenario]
    ext <- !is.null(d) && nrow(d) == 1 && isTRUE(d$selected) && r$n_trials > d$n_before
    n_ <- fmt(r$n_trials)
    if (!ext) { if (lang == "ko") sprintf("%s: %s회", lab, n_) else sprintf("%s: %s", lab, n_) }
    else if (r$n_trials >= d$n_after) { if (lang == "ko") sprintf("%s: 적응적 연장 규칙으로 %s회", lab, n_) else sprintf("%s: %s by the adaptive extension rule", lab, n_) }
    else if (lang == "ko") sprintf("%s: 적응적 연장 진행 중 %s회(목표 %s회)", lab, n_, fmt(d$n_after)) else sprintf("%s: %s so far, adaptive extension toward %s", lab, n_, fmt(d$n_after))
  }, "")
  if (exceptions_only) return(paste(items, collapse = "; "))
  sprintf(if (lang == "ko") "%s(%s)" else "%s (%s)", base, paste(items, collapse = "; "))
}
