#!/usr/bin/env Rscript
# 72_m15_individual.R — [6] 채혈 평가(개인 수준 20,000명, 세 일정 짝지음)와 6.1절 3번 기존 엔진 대조(legacy).
#   계획서 plans/M15_analysis_plan_20261020.md 8.4·6.1절, 사전 등록 config/prereg_20261020_m15.yaml individual (D-068).
# 모드
#   run [cores]     모델 x 집단 x 용량(8묶음): NCA 편향·정밀도, 외삽 분포, λz 성공률, Day 64/71 정량 몫 → results/m15/individual/
#                   [7] 분석식 입력(로그 SD, 상관, log 체중 기울기)도 여기서 낸다(ss_inputs_m15.csv)
#   legacy [cores]  기존 시험 모집단(정규 75 kg, SD 9, 60~90 kg 절단)·B0로 저장 결과(results/individual/individual_*.csv)와 3 MC SE 이내 대조
#   verify_pauc     m15_pauc 대 PKNCA 0.12.1 pk.calc.aucint(lin up/log down, AUClast; preprocess_blq 출력): 일정마다 허용창·잔차 있는 1,000개 프로필(6.1절 1번)
# AUC0-last 참값(8.4절): 주 = 조밀 채혈 참값(0부터 min(연속 시간 참 Tlast, 일정의 마지막 계획 시각)까지의 모델 AUC), 보조 = 관측 Tlast까지의 모델 AUC(적분 오차만)
# 사용법: nice -n 5 Rscript scripts/72_m15_individual.R <run|legacy|verify_pauc> [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("run", "legacy", "verify_pauc")) stop("사용법: Rscript scripts/72_m15_individual.R <run|legacy|verify_pauc> [cores]")
mode <- args[1]; cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
design <- read_cfg("trial_design.yaml"); LLOQ <- study_lloq()
MODELS <- c(CFG$models$primary, CFG$models$sensitivity); DOSES <- c(CFG$doses_mg$primary, CFG$doses_mg$comparison)
POPS <- CFG$populations; INCL <- CFG$inclusion; RHO <- as.numeric(POPS$weight_height_rho); SR <- as.numeric(POPS$sex_ratio_male)
SCHED <- lapply(CFG$schedules, as.numeric); GRID <- sort(unique(c(0, unlist(SCHED))))
NI <- as.integer(CFG$individual$n_subj); stopifnot(NI == 20000L)
od <- proj_path(CFG$output_dir, "individual"); dir.create(od, showWarnings = FALSE, recursive = TRUE)

# 한 묶음: 대상자(공통 난수) → 관측(허용창·잔차) → 세 일정 NCA + 참값
sim_individual <- function(m, pop, dose, schedules, n = NI, tag = "ind", incl = INCL) {
  p <- load_params(m)
  subj <- m15_subjects(n, pop, incl, RHO, SR, derive_seed(MS, tag, "subj"))
  e <- m15_eta_table(subj, p); ip <- individual_params(p, e)
  grid <- sort(unique(c(0, unlist(schedules)))); K <- length(grid)
  U <- with_seed(derive_seed(MS, tag, "jitter"), matrix(runif(n * K), n))
  Z <- with_seed(derive_seed(MS, tag, "eps"), list(p = matrix(rnorm(n * K), n), a = matrix(rnorm(n * K), n)))
  ob <- m15_obs_times(subj$id, grid, design, U)
  ep <- m15_eps_table(subj$id, grid, p$sigma, Z$p, Z$a)
  sol <- solve_model(ip, unique(ob[, .(id, time)]), dose, model_id = p$model_id)
  d <- merge(ob, sol[, .(id, time, C, auc)], by = c("id", "time"))
  d <- merge(d, ep, by = c("id", "planned"))
  d[, conc := m15_observe(C, planned, eps_p, eps_a, LLOQ)]
  tr <- m15_truth(ip, dose, p$model_id, schedules, LLOQ, list())      # AUC0-inf, Cmax, 일정별 조밀 채혈 AUC0-last 참값
  out <- rbindlist(lapply(names(schedules), function(sn) {
    x <- d[planned %in% c(0, schedules[[sn]])]
    nca <- run_nca(x[, .(id, time, conc)])
    at <- merge(nca[, .(id, Tlast)], x[, .(id, time, auc, C)], by.x = c("id", "Tlast"), by.y = c("id", "time"), all.x = TRUE)
    nca <- merge(nca, at[, .(id, AUClast_obs_tlast_true = auc, C_true_at_tlast = C)], by = "id")
    trs <- tr[, c("id", "AUCinf", "Cmax", "t_true_last", paste0("AUClast_dense_", sn)), with = FALSE]
    setnames(trs, c("id", "AUCinf_true", "Cmax_true", "t_true_last", "AUClast_true"))
    nca <- merge(nca, trs, by = "id")
    q <- x[planned %in% c(63, 70), .(obs_quant_63 = any(planned == 63 & !is.na(conc)), true_quant_63 = any(planned == 63 & C >= LLOQ),
                                     obs_quant_70 = any(planned == 70 & !is.na(conc)), true_quant_70 = any(planned == 70 & C >= LLOQ)), by = id]
    if (nrow(q)) nca <- merge(nca, q, by = "id", all.x = TRUE)
    nca[, schedule := sn]
  }), fill = TRUE)
  merge(out, subj[, .(id, sex, WT, HT, BMI)], by = "id")
}

pe <- function(est, tru) 100 * (est / tru - 1)
summ_schedule <- function(x) {
  ok_last <- is.finite(x$AUClast) & x$AUClast > 0; ok_inf <- x$lambda_ok & is.finite(x$AUCINF_obs); ok_A <- ok_inf & x$reliable
  e_last <- pe(x$AUClast[ok_last], x$AUClast_true[ok_last]); e_cmax <- pe(x$Cmax[ok_last], x$Cmax_true[ok_last])
  e_int <- pe(x$AUClast[ok_last], x$AUClast_obs_tlast_true[ok_last])
  e_inf <- pe(x$AUCINF_obs[ok_inf], x$AUCinf_true[ok_inf]); e_A <- pe(x$AUCINF_obs[ok_A], x$AUCinf_true[ok_A])
  bp <- function(e, lab) data.table(endpoint = lab, n = length(e), bias_median_pct = median(e), rmse_pct = sqrt(mean(e^2)), rmse_log = sqrt(mean(log(1 + e / 100)^2)),
                                    median_abs_pct = median(abs(e)), p05_pct = quantile(e, 0.05, names = FALSE), p95_pct = quantile(e, 0.95, names = FALSE))
  ex <- x$`AUC_%Extrap_obs`[x$lambda_ok]
  succ <- x$lambda_ok & x$`AUC_%Extrap_obs` < 20
  w_all <- wilson_ci(sum(succ, na.rm = TRUE), nrow(x)); w_c <- wilson_ci(sum(succ, na.rm = TRUE), sum(x$lambda_ok))
  ext <- data.table(n_subj = nrow(x), lambda_ok_pct = 100 * mean(x$lambda_ok), success_pct = w_all$est, success_lo = w_all$lo, success_hi = w_all$hi,
                    success_computable_pct = w_c$est, reliable_pct = 100 * mean(x$reliable),
                    extrap_p05 = quantile(ex, 0.05, names = FALSE), extrap_p25 = quantile(ex, 0.25, names = FALSE), extrap_p50 = quantile(ex, 0.50, names = FALSE),
                    extrap_p75 = quantile(ex, 0.75, names = FALSE), extrap_p95 = quantile(ex, 0.95, names = FALSE), extrap_p99 = quantile(ex, 0.99, names = FALSE),
                    extrap_gt20_pct = 100 * mean(ex > 20), tlast_median = median(x$Tlast, na.rm = TRUE),
                    tlast_noise_pct = 100 * mean(x$C_true_at_tlast < LLOQ, na.rm = TRUE),
                    obs_quant_63_pct = if ("obs_quant_63" %in% names(x)) 100 * mean(x$obs_quant_63, na.rm = TRUE) else NA_real_,
                    true_quant_63_pct = if ("true_quant_63" %in% names(x)) 100 * mean(x$true_quant_63, na.rm = TRUE) else NA_real_,
                    obs_quant_70_pct = if ("obs_quant_70" %in% names(x)) 100 * mean(x$obs_quant_70, na.rm = TRUE) else NA_real_,
                    true_quant_70_pct = if ("true_quant_70" %in% names(x)) 100 * mean(x$true_quant_70, na.rm = TRUE) else NA_real_)
  list(bias = rbind(bp(e_last, "AUClast"), bp(e_int, "AUClast_integration_only"), bp(e_inf, "AUCinf_computable"), bp(e_A, "AUCinf_ruleA"), bp(e_cmax, "Cmax")), ext = ext)
}
# [7] 분석식 입력: R/samplesize.R와 같은 정의(scripts/47 est_inputs와 같은 회귀)
ss_inputs <- function(x) {
  z <- x[is.finite(AUClast) & AUClast > 0 & is.finite(Cmax) & Cmax > 0]
  la <- log(z$AUClast); lc <- log(z$Cmax); lw <- log(z$WT) - mean(log(z$WT))
  ba <- unname(coef(lm(la ~ lw))[2]); bc <- unname(coef(lm(lc ~ lw))[2])
  data.table(n = nrow(z), sd_auc = sd(la), sd_cmax = sd(lc), rho_total = cor(la, lc), beta_auc = ba, beta_cmax = bc, ratio_sd_cmax_auc = sd(lc) / sd(la),
             cv_auc_pct = 100 * sqrt(exp(var(la)) - 1), cv_cmax_pct = 100 * sqrt(exp(var(lc)) - 1))
}

if (mode == "run") {
  lf <- start_run_log("m15_individual", master_seed = MS, run_mode = "final", extra = list(n = NI, cores = cores))
  jobs <- CJ(pk_model = MODELS, population = c("S1", "S2"), dose = DOSES, sorted = FALSE)
  for (m in MODELS) invisible(get_model(load_params(m)$model_id))
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
    x <- sim_individual(jobs$pk_model[k], POPS[[jobs$population[k]]], jobs$dose[k], SCHED)
    s <- lapply(names(SCHED), function(sn) { r <- summ_schedule(x[schedule == sn]); list(bias = r$bias[, schedule := sn], ext = r$ext[, schedule := sn], ss = ss_inputs(x[schedule == sn])[, schedule := sn]) })
    lab <- function(d) d[, `:=`(pk_model = jobs$pk_model[k], population = jobs$population[k], dose_mg = jobs$dose[k])]
    list(bias = lab(rbindlist(lapply(s, `[[`, "bias"))), ext = lab(rbindlist(lapply(s, `[[`, "ext"))), ss = lab(rbindlist(lapply(s, `[[`, "ss"))),
         subj = lab(x[, .(id, schedule, WT, sex, Cmax, AUClast, AUCINF_obs, lambda_ok, reliable, pct_extrap = `AUC_%Extrap_obs`, Tlast, t_true_last, AUClast_true,
                          AUClast_obs_tlast_true, AUCinf_true, Cmax_true)]))
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE); if (any(bad)) stop("개인 수준 실패: ", res[bad][[1]])
  B <- rbindlist(lapply(res, `[[`, "bias")); E <- rbindlist(lapply(res, `[[`, "ext")); S <- rbindlist(lapply(res, `[[`, "ss"))
  for (d in list(B, E, S)) setcolorder(d, c("pk_model", "population", "dose_mg", "schedule"))
  fwrite(B, file.path(od, "bias_precision.csv")); fwrite(E, file.path(od, "extrapolation_lambda_z.csv")); fwrite(S, file.path(od, "ss_inputs_m15.csv"))
  saveRDS(rbindlist(lapply(res, `[[`, "subj")), file.path(od, "individual_subjects.rds"))   # 상자 그림용(.gitignore 대상)
  # 일정 간 짝지은 비교(B0 기준): λz 성공률 차이(같은 대상자)
  print(B); print(E)
  append_run_log(lf, sprintf("done in %s", format(Sys.time() - t0)))
}

if (mode == "legacy") {
  lf <- start_run_log("m15_individual_legacy", master_seed = MS, run_mode = "final", extra = list(n = NI))
  w <- design$weight$base
  pop <- list(weight = list(dist = "truncnorm", mean = w$mean, sd = w$sd, lo = w$trunc[[1]], hi = w$trunc[[2]]), height = list(by_sex = FALSE, mean_cm = 170, sd_cm = 9))
  incl <- list(male_kg = as.numeric(unlist(w$trunc)), female_kg = as.numeric(unlist(w$trunc)))
  stored <- list(k2016 = fread(proj_path("results", "individual", "individual_base.csv"))[schedule == "B0"],
                 k2020 = fread(proj_path("results", "individual", "individual_struct2020.csv"))[schedule == "B0"])
  res <- parallel::mclapply(MODELS, function(m) {
    x <- sim_individual(m, pop, as.numeric(design$dose_mg), list(B0 = SCHED$B0), tag = "legacy", incl = incl)
    ex <- function(x_) x_[is.finite(x_)]
    st <- function(xx) c(lambda_ok_pct = 100 * mean(xx$lambda_ok), extrap_median = median(ex(xx$`AUC_%Extrap_obs`)), extrap_p95 = quantile(ex(xx$`AUC_%Extrap_obs`), 0.95, names = FALSE))
    est <- st(x)
    bs <- with_seed(derive_seed(MS, "legacy", "boot", m), t(replicate(200, st(x[sample.int(.N, .N, replace = TRUE)]))))
    se <- apply(bs, 2, sd)
    s <- stored[[m]]
    old <- c(lambda_ok_pct = s$lambda_ok_pct, extrap_median = s$extrap_median, extrap_p95 = s$extrap_p95)
    data.table(pk_model = m, metric = names(est), new = est, stored = old[names(est)], se_new = se, tol = 3 * sqrt(2) * se)[, agree := abs(new - stored) <= tol]
  }, mc.cores = 2L, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE); if (any(bad)) stop("legacy 실패: ", res[bad][[1]])
  L <- rbindlist(res); fwrite(L, file.path(od, "legacy_check.csv")); print(L)
  append_run_log(lf, sprintf("done: agree all %s", all(L$agree)))
}

if (mode == "verify_pauc") {
  lf <- start_run_log("m15_verify_pauc", master_seed = MS, run_mode = "final", extra = list(n = 1000))
  stopifnot(requireNamespace("PKNCA", quietly = TRUE)); stopifnot(as.character(packageVersion("PKNCA")) == "0.12.1")
  IV <- list(pAUC_0_14 = c(0, 14), pAUC_14_28 = c(14, 28), pAUC_28_56 = c(28, 56))
  p <- load_params(CFG$models$primary); n <- 1000L
  subj <- m15_subjects(n, POPS$S1, INCL, RHO, SR, derive_seed(MS, "verify_pauc", "subj"))
  ip <- individual_params(p, m15_eta_table(subj, p)); K <- length(GRID)
  U <- with_seed(derive_seed(MS, "verify_pauc", "jitter"), matrix(runif(n * K), n))
  Z <- with_seed(derive_seed(MS, "verify_pauc", "eps"), list(p = matrix(rnorm(n * K), n), a = matrix(rnorm(n * K), n)))
  ob <- m15_obs_times(subj$id, GRID, design, U); ep <- m15_eps_table(subj$id, GRID, p$sigma, Z$p, Z$a)
  sol <- solve_model(ip, unique(ob[, .(id, time)]), as.numeric(CFG$doses_mg$primary), model_id = p$model_id)
  d <- merge(merge(ob, sol[, .(id, time, C)], by = c("id", "time")), ep, by = c("id", "planned"))
  d[, conc := m15_observe(C, planned, eps_p, eps_a, LLOQ)]
  V <- rbindlist(lapply(names(SCHED), function(sn) {
    x <- d[planned %in% c(0, SCHED[[sn]]), .(id, time, conc)]
    a <- m15_pauc(x, IV); pre <- preprocess_blq(x)
    rbindlist(lapply(names(IV), function(nm) {
      b <- pre[, .(pk = if (any(conc > 0, na.rm = TRUE)) suppressWarnings(PKNCA::pk.calc.aucint(conc = conc, time = time, start = IV[[nm]][1], end = IV[[nm]][2],
                                                                                                 method = "lin up/log down", auc.type = "AUClast")) else NA_real_), by = id]
      z <- merge(a[, c("id", nm), with = FALSE], b, by = "id"); setnames(z, nm, "m15")
      both <- is.finite(z$pk) & is.finite(z$m15)
      data.table(schedule = sn, interval = nm, n_profiles = nrow(z), n_compared = sum(both), max_rel_diff = max(abs(z$m15[both] - z$pk[both]) / pmax(abs(z$pk[both]), 1e-12)),
                 n_pknca_na = sum(!is.finite(z$pk) & is.finite(z$m15)), n_pknca_na_m15_zero = sum(!is.finite(z$pk) & is.finite(z$m15) & z$m15 == 0),
                 n_m15_zero = sum(z$m15 == 0, na.rm = TRUE))
    }))
  }))
  V[, pass := max_rel_diff <= 1e-8 & n_pknca_na == n_pknca_na_m15_zero]
  fwrite(V, file.path(od, "pauc_pknca_check.csv")); print(V)
  append_run_log(lf, sprintf("done: all pass %s", all(V$pass)))
}
