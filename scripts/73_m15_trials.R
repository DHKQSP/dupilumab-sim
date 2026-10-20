#!/usr/bin/env Rscript
# 73_m15_trials.R — [3] OC 곡선(농도 배율 g), [4] 기전별 OC(시험군 F, ka, ke, Vmax 배율), [5] AUROC 입력(추정 log GMR).
#   계획서 plans/M15_analysis_plan_20261020.md 7.2~7.6·8.1~8.3·11절, 사전 등록 config/prereg_20261020_m15.yaml oc, mechanisms, auroc, compute (D-068).
# 모드(모델 k2020|k2016, 집단 S1|S2|S1_uniform|S2_uniform, 용량 300|200)
#   oc       <model> <pop> <dose> [from] [to] [cores]   g 21개 x 층 경계 2 x 일정 3 x 표본 수 9 x 분석 모형 2 → 시험별 통과 비트(두 평가변수, AUC0-last만, Cmax만)
#   ocext    <model> <pop> <dose> [from] [to] [cores]   g = 0.80, 1.25(와 대조 1.00)만(소비자 위험 정밀도, 시험 1,001~10,000; k2020만)
#   ocshape  <model> <pop> <dose> [from] [to] [cores]   체중 모양 민감도: 균등 집단, k2020, 300 mg, 확정 일정, g = 0.80, 0.95, 1.00, 1.25
#   mech     <model> <pop> <dose> [from] [to] [cores]   기전 판 33개(층 경계 75 kg) → 표본 수별 통과 수, 130명에서 시험별 통과·M1 추정 log GMR·부분 AUC 산술 평균
#   mechsupp <model> <pop> <dose> [from] [to] [cores]   보충 배율(ke, Vmax, ka x 0.60, 0.70, 1.30, 1.40, 1.50; k2020, 300 mg, 확정 일정)
#   pilotcheck <model> <pop> <dose> [from] [to]         11절 예비 실행 점검(oc g = 1.00 대 mech 항등의 시험별 일치, 분석식 검정력 대조, 단조성) → pilot_check_*.csv
# 시험 j의 대상자·배정·탈락·채혈 시각·잔차 시드는 m15_trial_seeds(master, 10,000)의 (j, 흐름)이며 모델·용량·집단·일정·표본 수·g·판과 무관(공통 난수).
# 묶음 파일(50회)은 임시 이름으로 쓴 뒤 바꾼다(이미 있는 묶음은 건너뜀 = 이어서 실행).
# 사용법: nice -n 5 Rscript scripts/73_m15_trials.R <oc|ocext|ocshape|mech|mechsupp|pilotcheck> <model> <pop> <dose> [from] [to] [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("oc", "ocext", "ocshape", "mech", "mechsupp", "pilotcheck")
if (length(args) < 4 || !args[1] %in% MODES) stop("사용법: Rscript scripts/73_m15_trials.R <", paste(MODES, collapse = "|"), "> <model> <pop> <dose> [from] [to] [cores]")
mode <- args[1]; model <- args[2]; popv <- args[3]; dose <- as.numeric(args[4])
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
stopifnot(model %in% c(CFG$models$primary, CFG$models$sensitivity), dose %in% c(CFG$doses_mg$primary, CFG$doses_mg$comparison))
NT_MAIN <- as.integer(CFG$oc$n_trials); NT_EXT <- as.integer(CFG$oc$extension$n_trials)
from <- if (length(args) >= 5) as.integer(args[5]) else if (mode == "ocext") NT_MAIN + 1L else 1L
to <- if (length(args) >= 6) as.integer(args[6]) else if (mode == "ocext") NT_EXT else NT_MAIN
cores <- if (length(args) >= 7) as.integer(args[7]) else 4L
POPS <- CFG$populations; INCL <- CFG$inclusion; RHO <- as.numeric(POPS$weight_height_rho); SR <- as.numeric(POPS$sex_ratio_male)
pmap <- list(S1 = c("S1", "weight"), S2 = c("S2", "weight"), S1_uniform = c("S1", "weight_sensitivity"), S2_uniform = c("S2", "weight_sensitivity"))
stopifnot(popv %in% names(pmap)); POP <- POPS[[pmap[[popv]][1]]]; WKEY <- pmap[[popv]][2]
# 모드별 등록 범위 점검
EXT <- CFG$oc$extension; SHP <- CFG$oc$weight_shape_sensitivity; SUP <- CFG$mechanisms$supplementary
if (mode == "ocext") stopifnot(model == EXT$model, to <= NT_EXT, popv %in% c("S1", "S2"))
if (mode == "ocshape") stopifnot(model == SHP$model, dose == as.numeric(SHP$dose_mg), popv %in% unlist(SHP$populations), to <= as.integer(SHP$n_trials))
if (mode == "mechsupp") stopifnot(model == SUP$model, dose == as.numeric(SUP$dose_mg), popv %in% unlist(SUP$populations), to <= as.integer(SUP$n_trials))
if (mode %in% c("oc", "ocext", "mech", "pilotcheck")) stopifnot(popv %in% c("S1", "S2"))
design <- read_cfg("trial_design.yaml"); LLOQ <- study_lloq(); stopifnot(isTRUE(all.equal(LLOQ, as.numeric(CFG$lloq_mg_L))))
SCHED <- lapply(CFG$schedules, as.numeric)
if (mode == "ocshape") SCHED <- SCHED[unlist(SHP$schedules)]
if (mode == "mechsupp") SCHED <- SCHED[unlist(SUP$schedules)]
GRID <- sort(unique(c(0, unlist(lapply(CFG$schedules, as.numeric))))); K <- length(GRID)    # 채혈 난수는 항상 세 일정의 합집합 격자(공통 난수)
NPOOL <- as.integer(CFG$enrollment_pool); NGRID <- as.integer(unlist(CFG$n_randomized_per_arm)); NMAX <- max(NGRID)
NBASE <- as.integer(CFG$mechanisms$base_n_randomized_per_arm); PDROP <- as.numeric(CFG$dropout$p)
BND <- fread(proj_path(CFG$output_dir, "populations", "strata_boundaries.csv"))[population == popv]
stopifnot(nrow(BND) == 2)
BOUNDS <- setNames(BND$boundary_kg, BND$boundary)
if (mode %in% c("mech", "mechsupp")) BOUNDS <- BOUNDS[CFG$mechanisms$boundary]
GALL <- m15_gmr_grid(CFG$oc$gmr_grid); stopifnot(length(GALL) == as.integer(CFG$oc$gmr_grid$n_total), any(GALL == 1))
G <- switch(mode, oc = GALL, ocext = c(1, as.numeric(unlist(EXT$gmr))), ocshape = as.numeric(unlist(SHP$gmr)), NULL)
if (!is.null(G)) stopifnot(sum(abs(G - 1) < 1e-12) == 1, length(G) <= 30)
INTERVALS <- list(pAUC_0_14 = as.numeric(unlist(CFG$auroc$endpoints$pAUC_0_14)), pAUC_14_28 = as.numeric(unlist(CFG$auroc$endpoints$pAUC_14_28)),
                  pAUC_28_56 = as.numeric(unlist(CFG$auroc$endpoints$pAUC_28_56)))
C_OFF <- as.numeric(CFG$auroc$pauc_offset_mg_day_L); stopifnot(isTRUE(all.equal(C_OFF, LLOQ * 1)))
SEEDS <- m15_trial_seeds(MS, NT_EXT); SM <- m15_seed_matrix(SEEDS)
sd_ <- function(trial_j, stream_s) m15_seed_of(SM, trial_j, stream_s)
p <- load_params(model); invisible(get_model(p$model_id))
out_dir <- proj_path(CFG$output_dir, "trials", if (mode == "pilotcheck") "pilot" else mode, sprintf("%s_%s_%g", model, popv, dose)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
lf <- start_run_log(sprintf("m15_trials_%s_%s_%s_%g_%d_%d", mode, model, popv, dose, from, to), master_seed = MS, run_mode = "final",
                    extra = list(mode = mode, model = model, population = popv, dose = dose, from = from, to = to, cores = cores, boundaries = BOUNDS,
                                 seed_collisions_resolved = SEEDS[resolved > 0, paste(j, stream, collapse = ";")]))

# ---------------------------------------------------------------- 시험 j 공통 자료 ----------------------------------------------------------------
trial_base <- function(j) {
  subj <- m15_subjects(NPOOL, POP, INCL, RHO, SR, sd_(j, "subj"), WKEY)
  ua <- with_seed(sd_(j, "alloc"), runif(NPOOL))
  ud <- with_seed(sd_(j, "drop"), runif(NPOOL))
  U <- with_seed(sd_(j, "jitter"), matrix(runif(NPOOL * K), NPOOL))
  Z <- with_seed(sd_(j, "eps"), list(p = matrix(rnorm(NPOOL * K), NPOOL), a = matrix(rnorm(NPOOL * K), NPOOL)))
  ob <- m15_obs_times(subj$id, GRID, design, U); ep <- m15_eps_table(subj$id, GRID, p$sigma, Z$p, Z$a)
  obs <- merge(ob, ep, by = c("id", "planned"))
  al <- lapply(BOUNDS, function(b) { a <- m15_allocate(subj$WT, b, ua); a[, `:=`(id = subj$id, ok = ud >= PDROP)]; a[arm_order <= NMAX] })
  for (b in names(al)) stopifnot(sum(al[[b]]$arm == "R") == NMAX, sum(al[[b]]$arm == "T") == NMAX)
  list(subj = subj, e = m15_eta_table(subj, p), obs = obs, al = al)
}
# ip(id 1..NPOOL)를 판 v마다 id + (v - 1) * NPOOL로 쌓아 한 번에 푼다. 반환: obs + v, C
solve_versions <- function(tb, ip_list) {
  stk <- rbindlist(lapply(seq_along(ip_list), function(v) copy(ip_list[[v]])[, id := id + (v - 1L) * NPOOL]))
  ob <- rbindlist(lapply(seq_along(ip_list), function(v) { x <- tb$obs[id %in% ip_list[[v]]$id]; x[, `:=`(v = v, sid = id, id = id + (v - 1L) * NPOOL)] }))
  sol <- solve_model(stk, unique(ob[, .(id, time)]), dose, model_id = p$model_id)
  merge(ob, sol[, .(id, time, C)], by = c("id", "time"))
}
# 분석 행: 대조군은 판 ref_v, 시험군은 판 v. nca: data.table(sid, v, schedule, <평가변수 열>)
analysis_rows <- function(al, nca, versions, ref_v, endpoints, bname) {
  rbindlist(lapply(versions, function(ver) {
    vv <- ver
    r <- merge(al[arm == "R"], nca[v == ref_v], by.x = "id", by.y = "sid", allow.cartesian = TRUE)
    t <- merge(al[arm == "T"], nca[v == vv], by.x = "id", by.y = "sid", allow.cartesian = TRUE)
    x <- rbind(r, t)
    rbindlist(lapply(endpoints, function(ep) x[, .(boundary = bname, schedule, vi = vv, endpoint = ep, arm, arm_order, ok, y = get(ep), s = stratum)]))   # vv: nca의 v 열과 구분
  }))
}
bits <- function(pass, vi) sum(bitwAnd(as.integer(pass), 1L) * 2^(vi - 1))

# ---------------------------------------------------------------- [3] OC ----------------------------------------------------------------
run_oc_trial <- function(j) {
  tb <- trial_base(j)
  ip <- individual_params(p, tb$e)
  d <- solve_versions(tb, list(ip))
  # 판 = g(농도 배율, 잔차 전): 같은 C·같은 잔차에 g를 곱한다
  x <- rbindlist(lapply(seq_along(G), function(gi) copy(d)[, `:=`(gi = gi, id = sid + (gi - 1L) * NPOOL, conc = m15_observe(C, planned, eps_p, eps_a, LLOQ, g = G[gi]))]))
  nca <- rbindlist(lapply(names(SCHED), function(sn) {
    z <- x[planned %in% c(0, SCHED[[sn]]), .(id, time, conc)]
    m15_nca_basic(z)[, .(sid = (id - 1L) %% NPOOL + 1L, v = (id - 1L) %/% NPOOL + 1L, schedule = sn, AUClast, Cmax)]
  }))
  g1 <- which(abs(G - 1) < 1e-12)
  rows <- rbindlist(lapply(names(tb$al), function(b) analysis_rows(tb$al[[b]], nca, seq_along(G), g1, c("AUClast", "Cmax"), b)))
  be <- m15_be_nested(rows, c("boundary", "schedule", "vi", "endpoint"), NGRID)
  w <- dcast(be, boundary + schedule + vi + n_rand + model ~ endpoint, value.var = "pass")
  w[, pass := AUClast & Cmax]
  ev <- be[endpoint == "AUClast" & model == "M0" & vi == g1, .(boundary, schedule, n_rand, n_R, n_T)]
  md <- be[model == "M1" & endpoint == "AUClast", .(m1_dropped = sum(m1_dropped, na.rm = TRUE)), by = .(boundary, schedule, n_rand)]
  ns2 <- rbindlist(lapply(names(tb$al), function(b) { a <- tb$al[[b]][ok == TRUE]                           # 층 2(경계 초과) 평가 가능 인원(층 경계 해석용)
    data.table(boundary = b, n_rand = NGRID, n_stratum2 = vapply(NGRID, function(n) sum(a$stratum == 1L & a$arm_order <= n), 1L)) }))
  mk <- w[, .(mask = bits(pass, vi), mask_auc = bits(AUClast, vi), mask_cmax = bits(Cmax, vi)), by = .(boundary, schedule, n_rand, model)]
  out <- merge(merge(merge(mk, ev, by = c("boundary", "schedule", "n_rand")), md, by = c("boundary", "schedule", "n_rand")), ns2, by = c("boundary", "n_rand"))
  out[, trial := j][]
}

# ---------------------------------------------------------------- [4] 기전 ----------------------------------------------------------------
VERS <- local({
  M <- CFG$mechanisms; v <- list(list(code = "identity", mult = list()))
  if (mode == "mechsupp") { for (pm in unlist(SUP$parameters)) for (m in as.numeric(unlist(SUP$multipliers))) v[[length(v) + 1]] <- list(code = sprintf("%s_%.2f", pm, m), mult = setNames(list(m), pm))
  } else for (pm in unlist(M$parameters)) for (m in as.numeric(unlist(M$multipliers))) v[[length(v) + 1]] <- list(code = sprintf("%s_%.2f", pm, m), mult = setNames(list(m), pm))
  v })
PAUC <- names(INTERVALS)
run_mech_trial <- function(j) {
  tb <- trial_base(j); al <- tb$al[[1]]
  tids <- al[arm == "T", id]
  ip_list <- c(list(individual_params(p, tb$e)), lapply(VERS[-1], function(vv) individual_params(apply_multipliers(p, vv$mult), tb$e[id %in% tids])))
  d <- solve_versions(tb, ip_list)
  d[, conc := m15_observe(C, planned, eps_p, eps_a, LLOQ)]
  nca <- rbindlist(lapply(names(SCHED), function(sn) {
    z <- d[planned %in% c(0, SCHED[[sn]]), .(id, time, conc)]
    a <- run_nca(z); pa <- m15_pauc(z, INTERVALS)
    x <- merge(a[, .(id, AUClast, Cmax, AUCinf_A = fifelse(reliable, AUCINF_obs, NA_real_), AUCinf_B = fifelse(lambda_ok, AUCINF_obs, NA_real_), Tlast, Clast, lambda_z)], pa, by = "id")
    iv <- INTERVALS$pAUC_28_56
    x[, pAUC_28_56_ext := m15_pauc_ext(pAUC_28_56, Tlast, Clast, lambda_z, iv[1], iv[2])]          # λz 외삽 민감도
    x[, `:=`(sid = (id - 1L) %% NPOOL + 1L, v = (id - 1L) %/% NPOOL + 1L, schedule = sn)][, c("id", "Tlast", "Clast", "lambda_z") := NULL]
  }))
  PA <- c(PAUC, "pAUC_28_56_ext")
  EPS <- c("AUClast", "Cmax", "AUCinf_A", "AUCinf_B", PA)
  rows <- analysis_rows(al, nca, seq_along(VERS), 1L, EPS, names(tb$al)[1])
  # 부분 AUC: 0도 유효한 값(유지). 산술 평균비 통계량과 0 비율(130명)은 원척도로, 분석 모형은 log(pAUC + c)
  am <- rows[endpoint %in% PA & arm_order <= NBASE & ok & is.finite(y),
             .(mean_R = mean(y[arm == "R"]), mean_T = mean(y[arm == "T"]), zero_R = sum(y[arm == "R"] == 0), zero_T = sum(y[arm == "T"] == 0),
               n_R = sum(arm == "R"), n_T = sum(arm == "T")), by = .(schedule, vi, endpoint)]
  rows[endpoint %in% PA, y := y + C_OFF]
  be <- m15_be_nested(rows, c("schedule", "vi", "endpoint"), NGRID)
  w <- dcast(be[endpoint %in% c("AUClast", "Cmax")], schedule + vi + n_rand + model ~ endpoint, value.var = "pass")
  w[, pass := AUClast & Cmax]
  cnt <- w[, .(pass = as.integer(pass), pass_auc = as.integer(AUClast), pass_cmax = as.integer(Cmax)), by = .(schedule, vi, n_rand, model)]
  est <- dcast(be[n_rand == NBASE & model == "M1"], schedule + vi ~ endpoint, value.var = "est")
  setnames(est, EPS, paste0("est_", EPS))
  pb <- dcast(w[n_rand == NBASE], schedule + vi ~ model, value.var = "pass"); setnames(pb, c("M0", "M1"), c("pass_M0", "pass_M1"))
  amw <- dcast(am, schedule + vi ~ endpoint, value.var = c("mean_R", "mean_T", "zero_R", "zero_T", "n_R", "n_T"))
  list(cnt = cnt[, trial := j], base = merge(merge(pb, est, by = c("schedule", "vi")), amw, by = c("schedule", "vi"), all.x = TRUE)[, trial := j])
}

# ---------------------------------------------------------------- 11절 예비 실행 점검 ----------------------------------------------------------------
if (mode == "pilotcheck") {
  PC <- CFG$compute$pilot_checks
  rd <- function(m, pat) { dd <- proj_path(CFG$output_dir, "trials", m, sprintf("%s_%s_%g", model, popv, dose)); f <- list.files(dd, pattern = pat, full.names = TRUE)
    list(x = rbindlist(lapply(f, fread)), v = fread(file.path(dd, "versions.csv"))) }
  oc <- rd("oc", "^oc_\\d+_\\d+\\.csv\\.gz$"); oc$x <- oc$x[trial >= from & trial <= to]
  me <- rd("mech", "^mech_\\d+_\\d+\\.csv\\.gz$"); me$x <- me$x[trial >= from & trial <= to]
  mc <- rd("mech", "^cnt_mech_\\d+_\\d+\\.csv\\.gz$")
  ntr <- to - from + 1L; stopifnot(uniqueN(oc$x$trial) == ntr, uniqueN(me$x$trial) == ntr, sum(mc$x[vi == 1 & n_rand == NBASE & model == "M0" & schedule == "B0", n_trials]) == ntr)
  g1 <- oc$v[abs(g - 1) < 1e-12, vi]; bit1 <- as.integer(2^(g1 - 1))
  # (1) 시험별 일치: oc g = 1.00, 75 kg 경계, 130명 대 mech 항등
  a <- oc$x[boundary == "fixed75" & n_rand == NBASE, .(trial, schedule, model, oc_pass = bitwAnd(as.integer(mask), bit1) > 0)]
  b <- melt(me$x[vi == 1], id.vars = c("trial", "schedule"), measure.vars = c("pass_M0", "pass_M1"), variable.name = "model", value.name = "mech_pass")[, model := sub("pass_", "", model)]
  ab <- merge(a, b, by = c("trial", "schedule", "model"))
  stopifnot(nrow(ab) == nrow(a))
  # 표본 수별 통과 수 일치(모든 n): oc 집계 대 mech 항등 집계(묶음 파일 합)
  ca <- oc$x[boundary == "fixed75", .(oc_n = sum(bitwAnd(as.integer(mask), bit1) > 0)), by = .(schedule, n_rand, model)]
  cb <- mc$x[vi == 1, .(mech_n = sum(n_pass)), by = .(schedule, n_rand, model)]
  cab <- merge(ca, cb, by = c("schedule", "n_rand", "model"))
  # (2) 분석식 검정력 대조(g = 0.95, 1.00; 확정 일정, 75 kg 경계, 130명, M0와 M1)
  SSI <- fread(proj_path(CFG$output_dir, "individual", "ss_inputs_m15.csv"))
  SSI <- SSI[SSI$pk_model == model & SSI$population == popv & SSI$dose_mg == dose & SSI$schedule == "B0"]
  stopifnot(nrow(SSI) == 1)
  wt <- m15_subjects(200000L, POP, INCL, RHO, SR, derive_seed(MS, "popq"), WKEY)$WT
  wm <- m15_weight_moments(wt, BOUNDS[["fixed75"]])
  inp <- list(sd_cmax = SSI$sd_cmax, ratio_sd_cmax_auc = SSI$ratio_sd_cmax_auc, rho_total = SSI$rho_total, beta_auc = SSI$beta_auc, beta_cmax = SSI$beta_cmax)
  ne <- (1 - PDROP) * NBASE
  AP <- rbindlist(lapply(c(0.95, 1.00), function(gg) { vi_ <- oc$v[abs(g - gg) < 1e-9, vi]; stopifnot(length(vi_) == 1)
    rbindlist(lapply(c("M0", "M1"), function(am) {
      k <- oc$x[boundary == "fixed75" & schedule == "B0" & n_rand == NBASE & model == am, sum(bitwAnd(as.integer(mask), as.integer(2^(vi_ - 1))) > 0)]
      sim <- 100 * k / ntr
      an <- 100 * ss_p2_power(ne, SSI$cv_auc_pct, gg, gg, inp, am, "fixed", wm)
      se <- 100 * sqrt((an / 100) * (1 - an / 100) / ntr)                 # MC SE at the analytic value (defined when the simulated share is 0 or 1)
      data.table(true_gmr = gg, model = am, sim_pct = sim, mc_se_pct = se, analytic_pct = an, tol_pct = max(as.numeric(PC$analytic_mc_se) * se, as.numeric(PC$analytic_min_pp)))
    })) }))[, pass := abs(sim_pct - analytic_pct) <= tol_pct]
  # (3) 단조성: 1.00을 향해 인접 g 사이 감소가 2 MC SE를 넘지 않는다(확정 일정, 두 경계, 130명, M0와 M1)
  MO <- oc$x[n_rand == NBASE & schedule == "B0", {
    pr <- vapply(oc$v$vi, function(v_) mean(bitwAnd(as.integer(mask), as.integer(2^(v_ - 1))) > 0), 1)
    gg <- oc$v$g; o <- order(gg); pr <- pr[o]; gg <- gg[o]; i1 <- which(gg == 1)
    se <- sqrt(pmax(pr * (1 - pr), 1 / ntr) / ntr)
    dn <- c(pr[2:i1] - pr[1:(i1 - 1)], pr[i1:(length(pr) - 1)] - pr[(i1 + 1):length(pr)])       # 1.00을 향해 증가해야 함
    sd2 <- c(sqrt(se[2:i1]^2 + se[1:(i1 - 1)]^2), sqrt(se[i1:(length(pr) - 1)]^2 + se[(i1 + 1):length(pr)]^2))
    .(max_drop_in_mc_se = max(-dn / sd2), pass = all(dn >= -as.numeric(PC$monotone_mc_se) * sd2)) }, by = .(boundary, model)]
  R1 <- data.table(check = "identity_trial_match", value = sum(ab$oc_pass != ab$mech_pass), n = nrow(ab), pass = all(ab$oc_pass == ab$mech_pass))
  R2 <- data.table(check = "identity_count_match_all_n", value = sum(cab$oc_n != cab$mech_n), n = nrow(cab), pass = all(cab$oc_n == cab$mech_n))
  fwrite(rbind(R1, R2), file.path(out_dir, sprintf("pilot_check_identity_%d_%d.csv", from, to)))
  fwrite(AP, file.path(out_dir, sprintf("pilot_check_analytic_%d_%d.csv", from, to))); fwrite(MO, file.path(out_dir, sprintf("pilot_check_monotone_%d_%d.csv", from, to)))
  print(rbind(R1, R2)); print(AP); print(MO)
  ok <- all(R1$pass, R2$pass, AP$pass, MO$pass)
  append_run_log(lf, sprintf("pilot checks: %s", if (ok) "all pass" else "FAILED (full run must not start; record a deviation)"))
  quit(save = "no", status = if (ok) 0L else 3L)
}

# 판 번호(vi)의 뜻: oc·ocext·ocshape는 g, mech·mechsupp는 기전·배율
meta <- if (mode %in% c("mech", "mechsupp")) data.table(vi = seq_along(VERS), code = vapply(VERS, `[[`, "", "code"),
                                                           mechanism = c("identity", vapply(VERS[-1], function(x) names(x$mult), "")),
                                                           multiplier = c(1, vapply(VERS[-1], function(x) as.numeric(x$mult[[1]]), 1))) else data.table(vi = seq_along(G), g = G)
fwrite(meta, file.path(out_dir, "versions.csv"))

# ---------------------------------------------------------------- 실행(50회 묶음) ----------------------------------------------------------------
BATCH <- 50L
starts <- seq((from - 1L) %/% BATCH * BATCH + 1L, to, by = BATCH)
fun <- if (mode %in% c("mech", "mechsupp")) run_mech_trial else run_oc_trial
for (b0 in starts) {
  js <- max(b0, from):min(b0 + BATCH - 1L, to)
  f <- file.path(out_dir, sprintf("%s_%05d_%05d.csv.gz", mode, min(js), max(js)))
  if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
  t0 <- Sys.time()
  res <- parallel::mclapply(js, function(j) tryCatch(fun(j), error = function(e) e), mc.cores = cores, mc.preschedule = TRUE)
  bad <- vapply(res, function(x) is.null(x) || inherits(x, "error") || inherits(x, "try-error"), TRUE)       # NULL = 강제 종료된 작업자(결과 없음)
  if (length(res) != length(js) || any(bad)) stop("시험 실패 또는 작업자 결과 없음: ", paste(js[bad], collapse = ","), " ",
                                                  if (is.null(res[bad][[1]])) "worker died" else conditionMessage(res[bad][[1]]))
  tmp <- file.path(out_dir, paste0("tmp_", basename(f)))
  if (mode %in% c("mech", "mechsupp")) {
    cnt <- rbindlist(lapply(res, `[[`, "cnt"))[, .(n_trials = .N, n_pass = sum(pass), n_pass_auc = sum(pass_auc), n_pass_cmax = sum(pass_cmax)), by = .(schedule, vi, n_rand, model)]
    fc <- file.path(out_dir, paste0("cnt_", basename(f))); fct <- file.path(out_dir, paste0("tmp_cnt_", basename(f))); fwrite(cnt, fct)
    base <- rbindlist(lapply(res, `[[`, "base"), fill = TRUE)
    for (cc in grep("^(est_|mean_)", names(base), value = TRUE)) set(base, j = cc, value = signif(base[[cc]], 7))
    stopifnot(setequal(unique(base$trial), js), all(cnt$n_trials == length(js)))
    fwrite(base, tmp); file.rename(fct, fc)
  } else { out <- rbindlist(res); stopifnot(setequal(unique(out$trial), js)); fwrite(out, tmp) }
  file.rename(tmp, f)
  msg <- sprintf("%s %s %s %g trials %d-%d written in %s", mode, model, popv, dose, min(js), max(js), format(Sys.time() - t0)); cat(msg, "\n"); append_run_log(lf, msg)
}
append_run_log(lf, "done")
