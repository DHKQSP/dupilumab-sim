#!/usr/bin/env Rscript
# 76_m15_param_assurance.R — 개정 1(계획서 16.8절, D-069): 파라미터 불확실성 보증(assurance)과 8.6절 3번 체중 효과 구간.
#   입력 config/param_uncertainty_m15.yaml(요청자 제공 값), 함수 R/pu.R. 사전 등록 config/prereg_20261020_m15.yaml amendment_20261021.param_uncertainty.
# 모드
#   check0 [cores]                    세트 0(추정값), 20,000명, 태그 'ind' → results/m15/individual/ss_inputs_m15.csv의 B0 300 mg 행 재현(상대 1e-12)
#   sets <set_id> [from] [to] [cores] 세트 from~to(50개 묶음, 이어서 실행): 모든 세트가 같은 2,000명(시드 derive_seed(20261020, 'pu', 흐름); 모델·집단과 무관, 세트 간 차이 = 파라미터만),
#                                     S1·S2에서 [7] 분석식 입력 재추정과 검정력(군당 107, 130명 무작위배정, 참 GMR 0.95, M1, 75 kg 경계),
#                                     전형 대상자 체중 효과 비(P5, P95 체중 / 75 kg; AUC0-inf, Cmax) → results/m15/amend1/param_assurance/
#                                     set_id: k2016 | k2016_thetas(구조 파라미터 8개만) | k2020_b16(Model 1, Kovalenko 2016 상대 폭 차용; 정정 2026-10-21)
#                                     검정력 보정(제어 변량): power_corr(s) = power(s) - power(0; 같은 2,000명) + power(0; 저장된 20,000명)
#   mcnoise [cores]                   세트 0에서 2,000명 표본 30개(시드 derive_seed(20261020, 'pu', 'mc', r, 흐름)): Monte Carlo 잡음
#   summary                           assurance_param.csv(보정·미보정), weight_ratio_param.csv
# 사용법: nice -n 5 Rscript scripts/76_m15_param_assurance.R <check0|sets|mcnoise|summary> ...
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R")); source(proj_path("R", "pu.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("check0", "sets", "mcnoise", "summary")
if (!length(args) || !args[1] %in% MODES) stop("사용법: Rscript scripts/76_m15_param_assurance.R <", paste(MODES, collapse = "|"), "> ...")
mode <- args[1]
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
AU <- CFG$amendment_20261021$param_uncertainty; stopifnot(!is.null(AU))
PU <- pu_cfg(AU$input); stopifnot(identical(PU$status, "registered"))
design <- read_cfg("trial_design.yaml"); LLOQ <- study_lloq()
POPS <- CFG$populations; INCL <- CFG$inclusion; RHO <- as.numeric(POPS$weight_height_rho); SR <- as.numeric(POPS$sex_ratio_male)
GRID <- sort(unique(c(0, unlist(lapply(CFG$schedules, as.numeric))))); B0 <- as.numeric(CFG$schedules$B0)   # [6]과 같은 합집합 격자, NCA는 확정 일정
DOSE <- as.numeric(AU$dose_mg); NSUB <- as.integer(AU$n_subj_per_set); NSETS <- as.integer(AU$n_sets); NRAND <- as.integer(unlist(AU$n_randomized_per_arm))
GMR <- as.numeric(AU$true_gmr); AM <- AU$analysis_model; BNDY <- AU$boundary; PV <- unlist(AU$populations)
SETS <- AU$sets                                                            # set_id -> {model, entries}
WM <- fread(proj_path(CFG$output_dir, "samplesize", "weight_moments.csv"))
WQ <- fread(proj_path(CFG$output_dir, "populations", "population_quantiles.csv"))
od <- proj_path(CFG$output_dir, "amend1", "param_assurance"); dir.create(od, showWarnings = FALSE, recursive = TRUE)

entries_of <- function(sid) { s <- SETS[[sid]]; if (is.null(s)) stop("set_id 없음: ", sid)
  E <- if (s$model == "k2016") pu_entries_k2016(PU) else pu_entries_k2020(PU, s$entries)
  if (!is.null(s$only)) E <- E[paste(group, name, sep = ".") %in% unlist(s$only)]                  # 범위 맞춤(k2016_thetas: 구조 파라미터 8개)
  E }
m15_ss_inputs <- function(x) {                                    # = scripts/72_m15_individual.R ss_inputs
  z <- x[is.finite(AUClast) & AUClast > 0 & is.finite(Cmax) & Cmax > 0]
  la <- log(z$AUClast); lc <- log(z$Cmax); lw <- log(z$WT) - mean(log(z$WT))
  ba <- unname(coef(lm(la ~ lw))[2]); bc <- unname(coef(lm(lc ~ lw))[2])
  data.table(n = nrow(z), sd_auc = sd(la), sd_cmax = sd(lc), rho_total = cor(la, lc), beta_auc = ba, beta_cmax = bc, ratio_sd_cmax_auc = sd(lc) / sd(la),
             cv_auc_pct = 100 * sqrt(exp(var(la)) - 1), cv_cmax_pct = 100 * sqrt(exp(var(lc)) - 1))
}
# [6]의 대상자·관측 생성(scripts/72 sim_individual과 같은 난수 순서)에서 확정 일정 NCA → 분석식 입력
set_inputs <- function(p, pop, seeds, n) {
  subj <- m15_subjects(n, pop, INCL, RHO, SR, seeds[["subj"]]); ip <- individual_params(p, m15_eta_table(subj, p)); K <- length(GRID)
  U <- with_seed(seeds[["jitter"]], matrix(runif(n * K), n))
  Z <- with_seed(seeds[["eps"]], list(p = matrix(rnorm(n * K), n), a = matrix(rnorm(n * K), n)))
  ob <- m15_obs_times(subj$id, GRID, design, U); ep <- m15_eps_table(subj$id, GRID, p$sigma, Z$p, Z$a)
  sol <- solve_model(ip, unique(ob[, .(id, time)]), DOSE, model_id = p$model_id)
  d <- merge(merge(ob, sol[, .(id, time, C)], by = c("id", "time")), ep, by = c("id", "planned")); d[, conc := m15_observe(C, planned, eps_p, eps_a, LLOQ)]
  m15_ss_inputs(merge(m15_nca_basic(d[planned %in% c(0, B0), .(id, time, conc)]), subj[, .(id, WT)], by = "id"))
}
pw <- function(ss, pv, nr) { wm <- WM[population == pv & boundary == BNDY, .(p_stratum2, mean_lwt, var_between, var_within, var_total)]; stopifnot(nrow(wm) == 1)
  inp <- list(sd_cmax = ss$sd_cmax, ratio_sd_cmax_auc = ss$ratio_sd_cmax_auc, rho_total = ss$rho_total, beta_auc = ss$beta_auc, beta_cmax = ss$beta_cmax)
  100 * ss_p2_power(0.9 * nr, ss$cv_auc_pct, GMR, GMR, inp, AM, "fixed", wm) }
seeds_of <- function(...) c(subj = derive_seed(MS, "pu", ..., "subj"), jitter = derive_seed(MS, "pu", ..., "jitter"), eps = derive_seed(MS, "pu", ..., "eps"))
# 8.6절 3번: 전형 대상자(η = 0) 체중 효과 비, 집단별 5·95번째 백분위수 체중 / 75 kg
WPCT <- as.numeric(unlist(AU$weight_ratio_percentiles))
wratio <- function(p) { w <- WQ[abs(prob - WPCT[1]) < 1e-9 | abs(prob - WPCT[2]) < 1e-9][population %in% PV]; stopifnot(nrow(w) == 2 * length(PV))
  tm <- m15_typical_metrics(p, c(w$WT, as.numeric(CFG$weight_effect$reference_kg)), DOSE, last_day = max(B0)); r <- tm[.N]
  data.table(population = w$population, prob = w$prob, WT = w$WT, ratio_AUCinf = tm$AUCinf[seq_len(nrow(w))] / r$AUCinf, ratio_Cmax = tm$Cmax[seq_len(nrow(w))] / r$Cmax) }

if (mode == "check0") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
  lf <- start_run_log("m15_pu_check0", master_seed = MS, run_mode = "final", extra = list(cores = cores))
  SSI <- fread(proj_path(CFG$output_dir, "individual", "ss_inputs_m15.csv"))
  jobs <- CJ(pk_model = c("k2020", "k2016"), population = PV, sorted = FALSE)
  for (m in unique(jobs$pk_model)) invisible(get_model(load_params(m)$model_id))
  sk <- c(subj = derive_seed(MS, "ind", "subj"), jitter = derive_seed(MS, "ind", "jitter"), eps = derive_seed(MS, "ind", "eps"))
  R <- rbindlist(parallel::mclapply(seq_len(nrow(jobs)), function(k) { m <- jobs$pk_model[k]; pv <- jobs$population[k]
    ss <- set_inputs(load_params(m), POPS[[pv]], sk, as.integer(CFG$individual$n_subj))
    st <- SSI[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0"]; stopifnot(nrow(st) == 1)
    cols <- setdiff(names(ss), "n"); rd <- max(abs(unlist(ss[, ..cols]) / unlist(st[, ..cols]) - 1))
    data.table(pk_model = m, population = pv, n = ss$n, n_stored = st$n, max_rel_diff = rd, power_n107 = pw(ss, pv, 107), power_n130 = pw(ss, pv, 130)) }, mc.cores = cores, mc.preschedule = FALSE))
  R[, pass := n == n_stored & max_rel_diff < 1e-10]; fwrite(R, file.path(od, "check0_set0_reproduction.csv")); print(R)
  if (!all(R$pass)) stop("세트 0 재현 실패")
  append_run_log(lf, "done: set 0 reproduces ss_inputs_m15.csv")
}

if (mode == "sets") {
  sid <- args[2]; from <- if (length(args) >= 3) as.integer(args[3]) else 0L; to <- if (length(args) >= 4) as.integer(args[4]) else NSETS
  cores <- if (length(args) >= 5) as.integer(args[5]) else 4L; stopifnot(from >= 0, to <= NSETS, from <= to)
  m <- SETS[[sid]]$model; p <- load_params(m); invisible(get_model(p$model_id)); E <- entries_of(sid); pu_check_estimates(p, E)
  st0 <- pu_param_set(p, E, 0, MS); stopifnot(isTRUE(all.equal(unclass(st0$p)[c("theta", "cov", "omega", "omega2", "sigma")], unclass(p)[c("theta", "cov", "omega", "omega2", "sigma")])))
  lf <- start_run_log(paste0("m15_pu_sets_", sid), master_seed = MS, run_mode = "final", extra = list(set_id = sid, from = from, to = to, cores = cores, n_subj = NSUB))
  BATCH <- 50L; starts <- seq(from %/% BATCH * BATCH, to, by = BATCH)
  for (b0 in starts) {
    ss_ <- max(b0, from):min(b0 + BATCH - 1L, to); f <- file.path(od, sprintf("sets_%s_%04d_%04d.csv.gz", sid, min(ss_), max(ss_)))
    if (file.exists(f)) next
    t0 <- Sys.time()
    sk <- seeds_of()                                                                           # 모든 세트가 같은 대상자(제어 변량 보정)
    R <- parallel::mclapply(ss_, function(s) { st <- pu_param_set(p, E, s, MS, max_redraw = as.integer(AU$redraw$max))
      wr <- wratio(st$p)
      rbindlist(lapply(PV, function(pv) { ss <- set_inputs(st$p, POPS[[pv]], sk, NSUB)
        w <- wr[population == pv]
        cbind(data.table(set_id = sid, pk_model = m, set = s, population = pv, redraws = st$redraws), ss,
              as.data.table(setNames(as.list(vapply(NRAND, function(nr) pw(ss, pv, nr), 1)), paste0("power_n", NRAND))),
              data.table(ratio_AUCinf_p05 = w[abs(prob - WPCT[1]) < 1e-9, ratio_AUCinf], ratio_AUCinf_p95 = w[abs(prob - WPCT[2]) < 1e-9, ratio_AUCinf],
                         ratio_Cmax_p05 = w[abs(prob - WPCT[1]) < 1e-9, ratio_Cmax], ratio_Cmax_p95 = w[abs(prob - WPCT[2]) < 1e-9, ratio_Cmax]),
              as.data.table(st$values)) })) }, mc.cores = cores, mc.preschedule = FALSE)
    bad <- vapply(R, function(x) is.null(x) || inherits(x, "try-error"), TRUE); if (length(R) != length(ss_) || any(bad)) stop("세트 실패: ", paste(ss_[bad], collapse = ","))
    tmp <- file.path(od, paste0("tmp_", basename(f))); fwrite(rbindlist(R), tmp); file.rename(tmp, f)
    append_run_log(lf, sprintf("%s sets %d-%d in %.1f s", sid, min(ss_), max(ss_), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  append_run_log(lf, "done")
}

if (mode == "mcnoise") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 4L; NR <- as.integer(AU$mc_noise_replicates)
  lf <- start_run_log("m15_pu_mcnoise", master_seed = MS, run_mode = "final", extra = list(cores = cores, replicates = NR))
  for (m in c("k2020", "k2016")) invisible(get_model(load_params(m)$model_id))
  R <- rbindlist(parallel::mclapply(seq_len(NR), function(r) rbindlist(lapply(c("k2020", "k2016"), function(m) { p <- load_params(m); sk <- seeds_of("mc", r)
    rbindlist(lapply(PV, function(pv) { ss <- set_inputs(p, POPS[[pv]], sk, NSUB)
      cbind(data.table(pk_model = m, replicate = r, population = pv), ss, as.data.table(setNames(as.list(vapply(NRAND, function(nr) pw(ss, pv, nr), 1)), paste0("power_n", NRAND)))) })) })),
    mc.cores = cores, mc.preschedule = FALSE))
  fwrite(R, file.path(od, "mc_noise_set0.csv")); append_run_log(lf, "done")
}

if (mode == "summary") {
  lf <- start_run_log("m15_pu_summary", master_seed = MS, run_mode = "final")
  COND <- m15_conditional_label(proj_path(CFG$output_dir), CFG$models$primary)
  SSI <- fread(proj_path(CFG$output_dir, "individual", "ss_inputs_m15.csv")); MC <- fread(file.path(od, "mc_noise_set0.csv"))
  WR <- fread(proj_path(CFG$output_dir, "weight", "weight_ratio.csv"))
  S <- rbindlist(lapply(names(SETS), function(sid) { f <- list.files(od, sprintf("^sets_%s_\\d+_\\d+\\.csv\\.gz$", sid), full.names = TRUE); if (!length(f)) return(NULL); rbindlist(lapply(f, fread)) }), fill = TRUE)
  stopifnot(nrow(S) > 0)
  chk <- S[, .(n = uniqueN(set), lo = min(set), hi = max(set)), by = .(set_id, population)]; if (any(chk$n != NSETS + 1 | chk$lo != 0 | chk$hi != NSETS)) { print(chk); stop("세트 누락") }
  A <- rbindlist(lapply(NRAND, function(nr) { pc <- paste0("power_n", nr)
    z0 <- S[set == 0, .(set_id, population, power_set0_2000_pct = get(pc))]
    p20 <- unique(S[, .(set_id, pk_model, population)])
    p20$power_set0_20000_pct <- vapply(seq_len(nrow(p20)), function(i) { mi <- p20$pk_model[i]; pv_ <- p20$population[i]
      st <- SSI[SSI$pk_model == mi & SSI$population == pv_ & SSI$dose_mg == DOSE & SSI$schedule == "B0"]; stopifnot(nrow(st) == 1); pw(st, pv_, nr) }, 1)
    x <- merge(merge(S[set >= 1, .(set_id, pk_model, population, set, redraws, pw = get(pc))], z0, by = c("set_id", "population")), p20, by = c("set_id", "pk_model", "population"))
    x[, pw_corr := pw - power_set0_2000_pct + power_set0_20000_pct]                                 # 제어 변량 보정(같은 대상자의 세트 0과의 차이를 20,000명 값에 더함)
    a <- x[, .(power_set0_2000_pct = power_set0_2000_pct[1], power_set0_20000_pct = power_set0_20000_pct[1],
               assurance_param_pct = mean(pw_corr), se_assurance_pp = sd(pw_corr) / sqrt(.N), median_pct = median(pw_corr),
               p025_pct = quantile(pw_corr, 0.025, names = FALSE), p975_pct = quantile(pw_corr, 0.975, names = FALSE), sd_between_sets_pp = sd(pw_corr),
               share_sets_ge90_pct = 100 * mean(pw_corr >= 90), assurance_uncorrected_pct = mean(pw), median_uncorrected_pct = median(pw),
               p025_uncorrected_pct = quantile(pw, 0.025, names = FALSE), p975_uncorrected_pct = quantile(pw, 0.975, names = FALSE),
               share_sets_ge90_uncorrected_pct = 100 * mean(pw >= 90), n_sets = .N, n_sets_redrawn = sum(redraws > 0)), by = .(set_id, pk_model, population)]
    mc <- MC[, .(mc_noise_mean_pct = mean(get(pc)), mc_noise_sd_pp = sd(get(pc))), by = .(pk_model, population)]
    a <- merge(a, mc, by = c("pk_model", "population"))
    a[, `:=`(n_rand = nr, n_eval = 0.9 * nr, true_gmr = GMR, distribution_label = "parameter-uncertainty distribution (not a confidence interval)",
             truncated_distribution = n_sets_redrawn > as.integer(AU$redraw_label_threshold))][] }))
  A[, `:=`(dose_mg = DOSE, schedule = "B0", boundary = BNDY, model = AM)]
  setcolorder(A, c("set_id", "pk_model", "population", "dose_mg", "schedule", "boundary", "model", "n_rand", "n_eval", "true_gmr"))
  m15_fwrite_cond(A, file.path(od, "assurance_param.csv"), COND)
  # 8.6절 3번: 전형 대상자 체중 효과 비의 파라미터 불확실성 구간(세트 1~N의 2.5~97.5번째 백분위수), 세트 0 = 저장 값과 대조
  W <- melt(S[population %in% PV, .(set_id, pk_model, set, population, ratio_AUCinf_p05, ratio_AUCinf_p95, ratio_Cmax_p05, ratio_Cmax_p95)], id.vars = c("set_id", "pk_model", "set", "population"))
  W[, `:=`(metric = sub("^ratio_(.*)_p\\d+$", "\\1", variable), prob = as.numeric(sub(".*_p", "", variable)) / 100)]
  WS <- W[, .(set0 = value[set == 0], p025 = quantile(value[set >= 1], 0.025, names = FALSE), median = median(value[set >= 1]), p975 = quantile(value[set >= 1], 0.975, names = FALSE)),
          by = .(set_id, pk_model, population, metric, prob)]
  WS[, stored := WR[match(paste(WS$pk_model, DOSE, WS$population, WS$prob, WS$metric), paste(WR$pk_model, WR$dose_mg, WR$population, WR$prob, WR$metric)), typical_ratio]]
  WS[, set0_matches_stored := abs(set0 / stored - 1) < 1e-6]
  m15_fwrite_cond(WS[, dose_mg := DOSE], file.path(od, "weight_ratio_param.csv"), COND)
  print(A); print(WS); append_run_log(lf, sprintf("done: %d sets", NSETS))
}
