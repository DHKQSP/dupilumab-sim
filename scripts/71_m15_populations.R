#!/usr/bin/env Rscript
# 71_m15_populations.R — [2] 집단(체중 분위수, 층 경계 중앙값, 집단 특성), [4](a)·[5] 참값(유도 참 GMR, 부분 AUC 참값), [8] 체중 효과.
#   계획서 plans/M15_analysis_plan_20261020.md 7.1·7.2·8.2·8.6절, 사전 등록 config/prereg_20261020_m15.yaml (D-068).
# 모드
#   populations                집단 S1, S2(및 균등 모양 민감도) 200,000명: 분위수, 층 경계 중앙값, 특성 → results/m15/populations/
#   truth [cores]              모델 x 집단 x 용량 x 기전 판(33개, k2020 300 mg은 보충 15개 추가) 20,000명 참값과 짝지은 모집단 GMR → results/m15/truth/
#   truthx [cores]             개정 1(계획서 16.3절 5, D-069): mechsuppx 판 가운데 truth_gmr.csv에 없는 판의 유도 참 GMR(같은 참값 대상자) → results/m15/truth/truth_gmr_ext.csv
#   weight [cores]             [8] 백분위 체중의 전형 대상자 비, 개인 짝 비 90% 구간, 개인차 예측구간, 정규화 민감도 계수 → results/m15/weight/
# 사용법: nice -n 5 Rscript scripts/71_m15_populations.R <populations|truth|weight> [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("populations", "truth", "truthx", "weight")) stop("사용법: Rscript scripts/71_m15_populations.R <populations|truth|truthx|weight> [cores]")
mode <- args[1]; cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
design <- read_cfg("trial_design.yaml"); LLOQ <- study_lloq(); stopifnot(isTRUE(all.equal(LLOQ, as.numeric(CFG$lloq_mg_L))))
MODELS <- c(CFG$models$primary, CFG$models$sensitivity); DOSES <- c(CFG$doses_mg$primary, CFG$doses_mg$comparison)
POPS <- CFG$populations; INCL <- CFG$inclusion; RHO <- as.numeric(POPS$weight_height_rho); SR <- as.numeric(POPS$sex_ratio_male)
POP_IDS <- c("S1", "S2")
pop_variants <- list(S1 = list(pop = "S1", wkey = "weight"), S2 = list(pop = "S2", wkey = "weight"),
                     S1_uniform = list(pop = "S1", wkey = "weight_sensitivity"), S2_uniform = list(pop = "S2", wkey = "weight_sensitivity"))
SCHED <- lapply(CFG$schedules, as.numeric)
INTERVALS <- list(pAUC_0_14 = c(0, 14), pAUC_14_28 = c(14, 28), pAUC_28_56 = c(28, 56))
stopifnot(isTRUE(all.equal(as.numeric(unlist(CFG$auroc$endpoints$pAUC_28_56)), c(28, 56))), isTRUE(all.equal(as.numeric(unlist(CFG$auroc$endpoints$pAUC_0_14)), c(0, 14))))

# ---------------------------------------------------------------- 집단 ----------------------------------------------------------------
if (mode == "populations") {
  od <- proj_path(CFG$output_dir, "populations"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
  lf <- start_run_log("m15_populations", master_seed = MS, run_mode = "final", extra = list(n = 200000))
  N <- 200000L; probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  seed <- derive_seed(MS, "popq")                                     # 모든 집단이 같은 균등 난수(공통 난수)
  Q <- list(); CH <- list(); BND <- list()
  for (pv in names(pop_variants)) {
    v <- pop_variants[[pv]]; s <- m15_subjects(N, POPS[[v$pop]], INCL, RHO, SR, seed, v$wkey)
    Q[[pv]] <- data.table(population = pv, prob = probs, WT = quantile(s$WT, probs, names = FALSE, type = 7))
    med <- round(median(s$WT), 1)
    BND[[pv]] <- data.table(population = pv, boundary = c("fixed75", "median"), boundary_kg = c(as.numeric(CFG$strata_boundary$fixed75), med),
                            share_upper_stratum = c(mean(s$WT > as.numeric(CFG$strata_boundary$fixed75)), mean(s$WT > med)))
    CH[[pv]] <- data.table(population = pv, n = N, share_male = mean(s$sex == "M"), WT_mean = mean(s$WT), WT_sd = sd(s$WT), WT_min = min(s$WT), WT_max = max(s$WT),
                           WT_mean_male = mean(s[sex == "M", WT]), WT_mean_female = mean(s[sex == "F", WT]),
                           HT_mean = mean(s$HT), HT_sd = sd(s$HT), cor_WT_HT = cor(s$WT, s$HT), BMI_mean = mean(s$BMI), BMI_sd = sd(s$BMI),
                           BMI_p05 = quantile(s$BMI, 0.05, names = FALSE), BMI_p95 = quantile(s$BMI, 0.95, names = FALSE), share_BMI_outside_18_32 = mean(s$BMI < 18 | s$BMI > 32))
  }
  Q <- rbindlist(Q); CH <- rbindlist(CH); BND <- rbindlist(BND)
  fwrite(Q, file.path(od, "population_quantiles.csv")); fwrite(CH, file.path(od, "population_characteristics.csv")); fwrite(BND, file.path(od, "strata_boundaries.csv"))
  print(Q); print(CH); print(BND)
  append_run_log(lf, sprintf("done: medians %s", paste(sprintf("%s %.1f", BND[boundary == "median", population], BND[boundary == "median", boundary_kg]), collapse = ", ")))
}

# ---------------------------------------------------------------- 참값 ----------------------------------------------------------------
# 판: 항등 + 4 파라미터 x 8 배율(33개). 보충 배율(ke, Vmax, ka x 0.60, 0.70, 1.30, 1.40, 1.50)은 등록된 모델·용량(k2020, 300 mg)에서만 더한다.
mech_versions <- function(supp = FALSE) {
  M <- CFG$mechanisms; v <- list(list(mechanism = "identity", multiplier = 1, mult = list(), set = "main"))
  for (pm in unlist(M$parameters)) for (m in as.numeric(unlist(M$multipliers))) v[[length(v) + 1]] <- list(mechanism = pm, multiplier = m, mult = setNames(list(m), pm), set = "main")
  if (supp) { S <- M$supplementary
    for (pm in unlist(S$parameters)) for (m in as.numeric(unlist(S$multipliers))) v[[length(v) + 1]] <- list(mechanism = pm, multiplier = m, mult = setNames(list(m), pm), set = "supplementary") }
  v
}
if (mode == "truth") {
  od <- proj_path(CFG$output_dir, "truth"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
  NT <- as.integer(CFG$mechanisms$truth$n_subj); stopifnot(NT == 20000L)
  lf <- start_run_log("m15_truth", master_seed = MS, run_mode = "final", extra = list(n = NT, cores = cores))
  stopifnot(length(mech_versions()) == 33L)
  SUPP <- CFG$mechanisms$supplementary
  jobs <- CJ(pk_model = MODELS, population = POP_IDS, dose = DOSES, sorted = FALSE)
  jobs[, supp := pk_model == SUPP$model & dose == as.numeric(SUPP$dose_mg)]
  for (m in MODELS) invisible(get_model(load_params(m)$model_id))
  one <- function(k) {
    m <- jobs$pk_model[k]; pop <- jobs$population[k]; dose <- jobs$dose[k]; p <- load_params(m); VERS <- mech_versions(jobs$supp[k])
    subj <- m15_subjects(NT, POPS[[pop]], INCL, RHO, SR, derive_seed(MS, "truth"))   # 집단·모델·용량 무관 시드(공통 난수)
    e <- m15_eta_table(subj, p)
    tt <- lapply(VERS, function(v) m15_truth(individual_params(apply_multipliers(p, v$mult), e), dose, p$model_id, SCHED, LLOQ, INTERVALS))
    ref <- tt[[1]]
    mets <- setdiff(names(ref), c("id", "t_true_last"))
    G <- rbindlist(lapply(seq_along(VERS), function(i) rbindlist(lapply(mets, function(mt) {
      r <- m15_paired_gmr(tt[[i]][[mt]], ref[[mt]])
      data.table(mechanism = VERS[[i]]$mechanism, multiplier = VERS[[i]]$multiplier, set = VERS[[i]]$set, metric = mt, gmr = r[["gmr"]], se_log = r[["se_log"]], n_used = r[["n"]],
                 share_zero_test = mean(tt[[i]][[mt]] <= 0), share_true_tlast_ge56 = mean(tt[[i]]$t_true_last >= 56)) }))))
    REF <- data.table(metric = mets, geo_mean = vapply(mets, function(mt) { x <- ref[[mt]]; exp(mean(log(x[x > 0]))) }, 1),
                      share_zero = vapply(mets, function(mt) mean(ref[[mt]] <= 0), 1))
    TL <- data.table(stat = c("median", "p05", "p95", "share_ge_28", "share_ge_56", "share_ge_70"),
                     value = c(quantile(pmin(ref$t_true_last, 1e6), c(0.5, 0.05, 0.95), names = FALSE), mean(ref$t_true_last >= 28), mean(ref$t_true_last >= 56), mean(ref$t_true_last >= 70)))
    lab <- function(d) d[, `:=`(pk_model = m, population = pop, dose_mg = dose)]
    list(gmr = lab(G), ref = lab(REF), tl = lab(TL))
  }
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) is.null(x) || inherits(x, "try-error"), TRUE); if (length(res) != nrow(jobs) || any(bad)) stop("참값 실패 또는 작업자 결과 없음: ", paste(which(bad), collapse = ","))
  G <- rbindlist(lapply(res, `[[`, "gmr")); R <- rbindlist(lapply(res, `[[`, "ref")); TL <- rbindlist(lapply(res, `[[`, "tl"))
  for (d in list(G, R, TL)) setcolorder(d, c("pk_model", "population", "dose_mg"))
  stopifnot(all(abs(G[mechanism == "identity", gmr] - 1) < 1e-12))
  fwrite(G, file.path(od, "truth_gmr.csv")); fwrite(R, file.path(od, "truth_reference.csv")); fwrite(TL, file.path(od, "truth_tlast.csv"))
  print(G[metric %in% c("AUCinf", "Cmax", "AUClast_dense_B0") & multiplier %in% c(0.8, 1.2)]); print(TL)
  append_run_log(lf, sprintf("done in %s", format(Sys.time() - t0)))
}

if (mode == "truthx") {
  od <- proj_path(CFG$output_dir, "truth"); X <- CFG$amendment_20261021$mechsuppx; stopifnot(!is.null(X))
  NT <- as.integer(CFG$mechanisms$truth$n_subj); TG <- fread(file.path(od, "truth_gmr.csv")); REF0 <- fread(file.path(od, "truth_reference.csv"))
  lf <- start_run_log("m15_truthx", master_seed = MS, run_mode = "final", extra = list(n = NT, cores = cores))
  want <- rbind(CJ(mechanism = unlist(X$parameters), multiplier = as.numeric(unlist(X$multipliers))),
                rbindlist(lapply(X$appended, function(ex) data.table(mechanism = ex$parameter, multiplier = as.numeric(unlist(ex$multipliers))))))
  jobs <- CJ(pk_model = unlist(X$models), population = unlist(X$populations), dose = as.numeric(X$dose_mg), sorted = FALSE)
  one <- function(k) {
    m <- jobs$pk_model[k]; pop <- jobs$population[k]; dose <- jobs$dose[k]; p <- load_params(m)
    have <- unique(TG[TG$pk_model == m & TG$population == pop & TG$dose_mg == dose, .(mechanism, multiplier)])
    need <- want[!have, on = c("mechanism", "multiplier")]; if (!nrow(need)) return(NULL)
    subj <- m15_subjects(NT, POPS[[pop]], INCL, RHO, SR, derive_seed(MS, "truth")); e <- m15_eta_table(subj, p)
    ref <- m15_truth(individual_params(p, e), dose, p$model_id, SCHED, LLOQ, INTERVALS)
    mets <- setdiff(names(ref), c("id", "t_true_last"))
    r0 <- REF0[REF0$pk_model == m & REF0$population == pop & REF0$dose_mg == dose]                      # 항등 참값이 저장된 것과 같아야 한다
    chk <- vapply(mets, function(mt) { x <- ref[[mt]]; abs(exp(mean(log(x[x > 0]))) / r0[metric == mt, geo_mean] - 1) }, 1); stopifnot(max(chk) < 1e-10)
    rbindlist(lapply(seq_len(nrow(need)), function(i) { v <- need[i]; tt <- m15_truth(individual_params(apply_multipliers(p, setNames(list(v$multiplier), v$mechanism)), e), dose, p$model_id, SCHED, LLOQ, INTERVALS)
      rbindlist(lapply(mets, function(mt) { r <- m15_paired_gmr(tt[[mt]], ref[[mt]])
        data.table(pk_model = m, population = pop, dose_mg = dose, mechanism = v$mechanism, multiplier = v$multiplier, set = "supplementary_ext", metric = mt, gmr = r[["gmr"]], se_log = r[["se_log"]],
                   n_used = r[["n"]], share_zero_test = mean(tt[[mt]] <= 0), share_true_tlast_ge56 = mean(tt$t_true_last >= 56)) })) }))
  }
  for (m in unique(jobs$pk_model)) invisible(get_model(load_params(m)$model_id))
  res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE); if (length(res) != nrow(jobs) || any(bad)) stop("truthx 실패: ", paste(which(bad), collapse = ","))
  GX <- rbindlist(res); fwrite(GX, file.path(od, "truth_gmr_ext.csv"))
  print(GX[metric %in% c("AUCinf", "Cmax", "AUClast_dense_B0")]); append_run_log(lf, "done")
}

# ---------------------------------------------------------------- [8] 체중 효과 ----------------------------------------------------------------
# 점추정: 전형 대상자(η = 0) X(w) / X(75 kg). 구간(주): 같은 η의 개인 짝 비 X_i(w) / X_i(75 kg)의 5~95번째 백분위수(무게 효과 자체의 개인 간 퍼짐),
# 짝 비의 기하평균도 낸다. 별도 열: 개인차 예측구간 X_i(w) / 전형 X(75 kg). 지표: AUC0-inf, Cmax, 조밀 채혈 AUC0-last 참값(확정 일정).
if (mode == "weight") {
  od <- proj_path(CFG$output_dir, "weight"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
  Q <- fread(proj_path(CFG$output_dir, "populations", "population_quantiles.csv"))
  WE <- CFG$weight_effect; REFKG <- as.numeric(WE$reference_kg); NIIV <- as.integer(WE$n_subj_iiv); stopifnot(NIIV == 20000L)
  lf <- start_run_log("m15_weight", master_seed = MS, run_mode = "final", extra = list(n_iiv = NIIV, cores = cores))
  wts <- unique(c(REFKG, Q$WT)); MET <- c("AUCinf", "Cmax", "AUClast_dense_B0")
  jobs <- CJ(pk_model = MODELS, dose = as.numeric(unlist(WE$doses_mg)), sorted = FALSE)
  for (m in MODELS) invisible(get_model(load_params(m)$model_id))
  one <- function(k) {
    m <- jobs$pk_model[k]; dose <- jobs$dose[k]; p <- load_params(m)
    typ <- m15_typical_metrics(p, wts, dose, last_day = max(SCHED$B0))   # 전형 대상자(η = 0), 0.01일 격자
    setnames(typ, "AUClast_dense", "AUClast_dense_B0")
    t75 <- typ[abs(WT - REFKG) < 1e-9]
    base <- m15_subjects(NIIV, POPS$S1, INCL, RHO, SR, derive_seed(MS, "weight_iiv"))   # η 점수만 쓴다(모든 체중에 같은 점수)
    tr_at <- function(w) { e <- m15_eta_table(copy(base)[, WT := w], p); m15_truth(individual_params(p, e), dose, p$model_id, list(B0 = SCHED$B0), LLOQ, list()) }
    r75 <- tr_at(REFKG)
    iiv <- rbindlist(lapply(wts, function(w) {
      tr <- tr_at(w); stopifnot(identical(tr$id, r75$id))
      rbindlist(lapply(MET, function(mt) { pr <- tr[[mt]] / r75[[mt]]; pi_ <- tr[[mt]] / t75[[mt]]
        data.table(WT = w, metric = mt, paired_gm = exp(mean(log(pr))), paired_q05 = quantile(pr, 0.05, names = FALSE), paired_q50 = quantile(pr, 0.50, names = FALSE),
                   paired_q95 = quantile(pr, 0.95, names = FALSE), pi_q05 = quantile(pi_, 0.05, names = FALSE), pi_q50 = quantile(pi_, 0.50, names = FALSE), pi_q95 = quantile(pi_, 0.95, names = FALSE)) }))
    }))
    pt <- rbindlist(lapply(MET, function(mt) data.table(WT = typ$WT, metric = mt, typical_ratio = typ[[mt]] / t75[[mt]], typical_value = typ[[mt]])))
    sc <- m15_sensitivity(p, wts, dose, unlist(WE$sensitivity_coefficients$parameters), as.numeric(WE$sensitivity_coefficients$step))
    list(ratio = merge(pt, iiv, by = c("WT", "metric"))[, `:=`(pk_model = m, dose_mg = dose)], sens = sc[, `:=`(pk_model = m, dose_mg = dose)])
  }
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = min(cores, nrow(jobs)), mc.preschedule = FALSE)
  bad <- vapply(res, function(x) is.null(x) || inherits(x, "try-error"), TRUE); if (length(res) != nrow(jobs) || any(bad)) stop("체중 효과 실패 또는 작업자 결과 없음: ", paste(which(bad), collapse = ","))
  RT <- rbindlist(lapply(res, `[[`, "ratio")); SC <- rbindlist(lapply(res, `[[`, "sens"))
  map <- rbind(data.table(population = "reference", prob = NA_real_, WT = REFKG), Q)
  RTm <- merge(map, RT, by = "WT", allow.cartesian = TRUE); SCm <- merge(map, SC, by = "WT", allow.cartesian = TRUE)
  setcolorder(RTm, c("pk_model", "dose_mg", "population", "prob", "WT", "metric")); setcolorder(SCm, c("pk_model", "dose_mg", "population", "prob", "WT", "param"))
  setorder(RTm, pk_model, dose_mg, population, prob, metric); setorder(SCm, pk_model, dose_mg, population, prob, param)
  COND <- m15_conditional_label(proj_path(CFG$output_dir), CFG$models$primary)   # 조건부 표시(계획서 6.2절)
  m15_fwrite_cond(RTm, file.path(od, "weight_ratio.csv"), COND); m15_fwrite_cond(SCm, file.path(od, "sensitivity_coefficients.csv"), COND)
  print(RTm[population %in% c("S1", "S2") & prob %in% c(0.05, 0.95)]); print(SCm[population %in% c("S1", "S2") & prob %in% c(0.05, 0.95)])
  append_run_log(lf, sprintf("done in %s", format(Sys.time() - t0)))
}
