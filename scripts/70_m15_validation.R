#!/usr/bin/env Rscript
# 70_m15_validation.R — [1] 모델 확인과 외부 검증(사후 예측 점검). 계획서 plans/M15_analysis_plan_20261020.md 4.2·5·6.2절,
#   사전 등록 config/prereg_20261020_m15.yaml validation (D-068). 모델 구조·추정값은 그대로 쓴다.
# 판정 단위는 독립 연구 5개(TDU12265 300, 600 mg, PKM12350, PKM14161, R668-HV-1108; 평행군은 두 군 합침), 95% 예측구간(R4는 99%).
# 체중 정보가 없는 연구: 주 판정은 등록 가정 B14(고정 평균 78 kg, PKM14271 80 kg, SD 10). 두 번째 판정은 반복마다 연구 수준 평균 체중을
#   Uniform(가정 평균 ± 8 kg)에서 뽑아 두 군이 공유(가정 체중의 불확실성). 고정 70, 86 kg은 기술 민감도.
# 산출: results/m15/validation/ — covariate_inputs.csv, validation_reps.csv.gz(반복별 군 통계), validation_arms.csv(군, 기술 통계),
#   validation_units.csv(판정 단위 예측구간 대 관측), validation_pooled.csv(비아시아 합동군, 기술 통계), validation_gmr.csv(평행군 GMR·연구 내 변동, 기술 통계),
#   validation_rule.csv(R1, R2, R4, R5 판정, 모의 오경보율, 검출 확률; 주 가정과 평균 불확실 가정), validation_rule_weight_sensitivity.csv, step1_arm_check_prior.csv, validation_summary_ko.md
# 사용법: nice -n 5 Rscript scripts/70_m15_validation.R [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE); cores <- if (length(args)) as.integer(args[1]) else 4L
CFG <- m15_cfg(); V <- CFG$validation; dz <- read_cfg("design_clot2021.yaml")
MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
REPS <- as.integer(V$replicates); REPS_S <- as.integer(V$replicates_sensitivity); stopifnot(REPS_S <= REPS)
PI <- as.numeric(V$pi_level); qlo <- (1 - PI) / 2; qhi <- 1 - qlo
SR <- as.numeric(V$sex_ratio_male); LLOQ <- as.numeric(V$lloq_mg_L); stopifnot(abs(LLOQ - as.numeric(dz$lloq_mg_L)) < 1e-12)   # 문헌 자료(원개발사 분석법)의 LLOQ
MODELS <- c(CFG$models$primary, CFG$models$sensitivity)
WSENS <- as.numeric(unlist(V$weight_sensitivity_means_kg)); HW <- as.numeric(V$weight_mean_uncertainty$half_width_kg)
MODES <- c("primary", "wunc", paste0("w", WSENS))                     # primary = 고정 가정 평균(B14), wunc = 평균 불확실, w70·w86 = 고정 민감도
FULL_MODES <- c("primary", "wunc")                                     # REPS회; 나머지는 앞 REPS_S회
out_dir <- proj_path(CFG$output_dir, "validation"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
lf <- start_run_log("m15_validation", master_seed = MS, run_mode = "final", extra = list(replicates = REPS, replicates_sensitivity = REPS_S, models = MODELS, cores = cores))

# ---------------------------------------------------------------- 공변량 입력 보고(4.2절) ----------------------------------------------------------------
cov_tab <- rbindlist(lapply(MODELS, function(m) {
  p <- load_params(m)
  data.table(pk_model = m, covariate = c("Body weight", "EASI", "Albumin", "Race", "BMI", "Sex"),
             in_model = c(sprintf("Yes: Vc x (WT / %g kg)^%g", p$cov$WT_ref, p$cov$theta_WT), "No", "No", "No",
                          if (!is.null(p$cov$ke_bmi_exp) && p$cov$ke_bmi_exp != 0) "Yes" else "No", "No"),
             healthy_input = c("simulated from each study's or scenario's weight distribution", "none (not a model input)", "none (not a model input)",
                               "none (ethnicity enters only through the weight distribution)", "reported only (no effect on PK)", "sex ratio 0.5 (assumption); weight bounds and height only"),
             flag = c("none", "none: no patient EASI value is used", "none", "none", "none", "none"))
}))
stopifnot(all(cov_tab[covariate == "BMI", in_model] == "No"))
fwrite(cov_tab, file.path(out_dir, "covariate_inputs.csv"))

# ---------------------------------------------------------------- 군 모의 ----------------------------------------------------------------
dev <- dz$model_development_data
internal_of <- function(study, model) study %in% unlist(dev[[model]])
wdef <- V$weight_default_bounds; BM <- as.numeric(unlist(wdef$male)); BF <- as.numeric(unlist(wdef$female))
is_assumed <- function(a) identical(a$weight$dist, "truncnorm")         # 체중 정보 없음(가정 평균)
# 연구 수준 평균 체중(반복별): primary는 가정 평균(고정), wunc는 Uniform(가정 평균 ± HW), w<kg>는 그 값. 연구 이름만 시드 태그에 넣는다(두 군·두 모델·모드 공통).
study_mu <- function(a, mode) {
  m0 <- as.numeric(a$weight$mean)
  if (mode == "primary") return(rep(m0, REPS))
  if (mode == "wunc") return(with_seed(derive_seed(MS, "val", a$study, "wmean"), runif(REPS, m0 - HW, m0 + HW)))
  rep(as.numeric(sub("^w", "", mode)), REPS)
}
# 한 군 REPS 반복 x n_subj명. 대상자·잔차 난수는 군 이름만 시드 태그(모델·모드와 무관 → 공통 난수); 체중만 모드에 따라 다르다.
sim_arm <- function(m, aid, mode) {
  a <- V$arms[[aid]]; p <- load_params(m); n_s <- as.integer(a$n_subj); N <- n_s * REPS
  r <- with_seed(derive_seed(MS, "val", aid, "subj"), list(u_sex = runif(N), u_w = runif(N), z = matrix(rnorm(N * length(M15_ETA_NAMES)), N)))
  male <- r$u_sex < SR; rep_ <- (seq_len(N) - 1L) %/% n_s + 1L; wt <- numeric(N)
  if (is_assumed(a)) {
    mu <- study_mu(a, mode)[rep_]; sdw <- as.numeric(a$weight$sd)
    wt[male] <- m15_q_weight(r$u_w[male], list(dist = "truncnorm", mean = mu[male], sd = sdw, lo = -Inf, hi = Inf), BM[1], BM[2])
    wt[!male] <- m15_q_weight(r$u_w[!male], list(dist = "truncnorm", mean = mu[!male], sd = sdw, lo = -Inf, hi = Inf), BF[1], BF[2])
  } else {
    stopifnot(mode == "primary"); s <- m15_weight_spec(a$weight)
    wt <- m15_q_weight(r$u_w, s, s$lo, s$hi)                            # 범위만 알려진 연구: 남녀 같은 범위(5절)
  }
  keep <- rep_ <= (if (mode %in% FULL_MODES) REPS else REPS_S)
  subj <- data.table(id = seq_len(N), WT = wt, rep = rep_)
  zt <- as.data.table(r$z); setnames(zt, paste0("z_", M15_ETA_NAMES)); subj <- cbind(subj, zt)[keep]
  ip <- individual_params(p, m15_eta_table(subj, p))
  planned <- sort(unique(c(0, schedule_days(dz, a$schedule)))); K <- length(planned)
  Z <- with_seed(derive_seed(MS, "val", aid, "eps"), list(p = matrix(rnorm(N * K), N), a = matrix(rnorm(N * K), N)))
  ob <- CJ(id = subj$id, planned = planned)[, time := planned]
  ep <- m15_eps_table(subj$id, planned, p$sigma, Z$p[keep, , drop = FALSE], Z$a[keep, , drop = FALSE])
  sim <- simulate_observations(ip, ob, as.numeric(a$dose_mg), LLOQ, ep, model_id = p$model_id)
  nca <- m15_nca_basic(sim$obs[, .(id, time, conc)])
  merge(nca, subj[, .(id, rep, WT)], by = "id")[, `:=`(pk_model = m, arm = aid, wmode = mode)][]
}
jobs <- rbindlist(lapply(MODELS, function(m) rbindlist(lapply(names(V$arms), function(aid)
  data.table(pk_model = m, arm = aid, wmode = if (is_assumed(V$arms[[aid]])) MODES else "primary")))))
for (m in MODELS) invisible(get_model(load_params(m)$model_id))     # 포크 전 컴파일
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) sim_arm(jobs$pk_model[k], jobs$arm[k], jobs$wmode[k]), mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE); if (any(bad)) stop("모의 실패: ", paste(which(bad), collapse = ","), " ", res[bad][[1]])
append_run_log(lf, sprintf("%d arm jobs simulated in %s", nrow(jobs), format(Sys.time() - t0)))
NCA <- rbindlist(res); rm(res)
cv_ <- function(x) 100 * sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
# 반복 통계: 평균, CV(합친 SD), 군 안 CV(군마다 중심화한 제곱합 / (N - k)). by는 rep를 포함하고 arm 열이 있어야 한다.
stat_of <- function(d, by) {
  bya <- unique(c(by, "arm"))
  x <- copy(d)[, `:=`(ac = AUClast - mean(AUClast, na.rm = TRUE), cc = Cmax - mean(Cmax, na.rm = TRUE)), by = bya]
  x[, .(n = .N, AUClast_mean = mean(AUClast, na.rm = TRUE), AUClast_cv = cv_(AUClast),
        AUClast_cv_within = 100 * sqrt(sum(ac^2, na.rm = TRUE) / (sum(is.finite(ac)) - uniqueN(arm[is.finite(ac)]))) / mean(AUClast, na.rm = TRUE),
        Cmax_mean = mean(Cmax, na.rm = TRUE), Cmax_cv = cv_(Cmax),
        Cmax_cv_within = 100 * sqrt(sum(cc^2, na.rm = TRUE) / (sum(is.finite(cc)) - uniqueN(arm[is.finite(cc)]))) / mean(Cmax, na.rm = TRUE),
        n_auc_missing = sum(!is.finite(AUClast)), WT_mean = mean(WT)), by = by]
}
REP <- stat_of(NCA, c("pk_model", "arm", "wmode", "rep"))
stopifnot(REP[wmode %in% FULL_MODES, .N, by = .(pk_model, arm, wmode)]$N == REPS, REP[!wmode %in% FULL_MODES, .N, by = .(pk_model, arm, wmode)]$N == REPS_S)
fwrite(REP, file.path(out_dir, "validation_reps.csv.gz"))

# ---------------------------------------------------------------- 관측값(군, 합동 단위) ----------------------------------------------------------------
ob_arm <- function(aid) { o <- V$arms[[aid]]$observed; if (is.null(o)) return(NULL)
  g <- function(nm) if (is.null(o[[nm]])) NA_real_ else as.numeric(o[[nm]])
  list(n = as.numeric(V$arms[[aid]]$n_subj), AUClast_mean = g("AUClast_mean"), AUClast_sd = g("AUClast_sd"), Cmax_mean = g("Cmax_mean"), Cmax_sd = g("Cmax_sd"), source = o$source) }
# 여러 군 합침: 평균 = n 가중 평균, CV = 합동 제곱합(군 안 + 군 간) / (N - 1)의 SD / 평균, 군 안 CV = sqrt(sum (n_i - 1) SD_i^2 / (N - k)) / 평균(R2용: 제품 차이 항 제외)
combine_obs <- function(lst, e) {
  n <- vapply(lst, `[[`, 1, "n"); mu <- vapply(lst, `[[`, 1, paste0(e, "_mean")); sdv <- vapply(lst, `[[`, 1, paste0(e, "_sd"))
  if (anyNA(mu) || anyNA(sdv)) return(c(mean = NA_real_, cv = NA_real_, cv_within = NA_real_))
  N <- sum(n); m <- sum(n * mu) / N; ssw <- sum((n - 1) * sdv^2); ss <- ssw + sum(n * (mu - m)^2)
  c(mean = m, cv = 100 * sqrt(ss / (N - 1)) / m, cv_within = 100 * sqrt(ssw / (N - length(n))) / m)
}
pi_row <- function(x, obs) {
  x <- x[is.finite(x)]; lo <- quantile(x, qlo, names = FALSE, type = 7); hi <- quantile(x, qhi, names = FALSE, type = 7)
  if (is.na(obs)) return(data.table(pred_median = median(x), pi_lo = lo, pi_hi = hi, observed = NA_real_, pct_rank = NA_real_, rank_lo = NA_real_, rank_hi = NA_real_,
                                    outside_pi = NA, borderline = NA, n_rep = length(x)))
  k <- sum(x <= obs); w <- wilson_ci(k, length(x))
  data.table(pred_median = median(x), pi_lo = lo, pi_hi = hi, observed = obs, pct_rank = w$est, rank_lo = w$lo, rank_hi = w$hi,
             outside_pi = obs < lo | obs > hi, borderline = (w$lo <= 100 * qlo & w$hi >= 100 * qlo) | (w$lo <= 100 * qhi & w$hi >= 100 * qhi), n_rep = length(x))
}
dev_status <- function(studies, m) { inter <- vapply(studies, internal_of, TRUE, model = m); if (all(inter)) "internal" else if (any(inter)) "partly internal" else "external" }
METS <- c("AUClast_mean", "AUClast_cv", "AUClast_cv_within", "Cmax_mean", "Cmax_cv", "Cmax_cv_within")
four_rows <- function(st, o) rbindlist(lapply(METS, function(mt) cbind(data.table(metric = mt), pi_row(st[[mt]], o[[mt]]))))
obs_vec <- function(lst) { a <- combine_obs(lst, "AUClast"); c_ <- combine_obs(lst, "Cmax")
  c(AUClast_mean = a[["mean"]], AUClast_cv = a[["cv"]], AUClast_cv_within = a[["cv_within"]], Cmax_mean = c_[["mean"]], Cmax_cv = c_[["cv"]], Cmax_cv_within = c_[["cv_within"]]) }
NA4 <- setNames(rep(NA_real_, length(METS)), METS)

# 군 표(기술 통계): 모든 군 x 모드
ARMS <- rbindlist(lapply(seq_len(nrow(jobs)), function(k) {
  m <- jobs$pk_model[k]; aid <- jobs$arm[k]; md <- jobs$wmode[k]; a <- V$arms[[aid]]; o <- ob_arm(aid)
  four_rows(REP[pk_model == m & arm == aid & wmode == md], if (is.null(o)) NA4 else obs_vec(list(o)))[,
    `:=`(pk_model = m, arm = aid, wmode = md, study = a$study, dose_mg = a$dose_mg, n_subj = as.numeric(a$n_subj), dev_status = dev_status(a$study, m),
         obs_source = if (is.null(o)) "not in repository" else o$source)]
}))
setcolorder(ARMS, c("pk_model", "arm", "wmode", "study", "dose_mg", "n_subj", "dev_status", "metric")); fwrite(ARMS, file.path(out_dir, "validation_arms.csv"))

# 판정 단위 표: 단위의 군들을 반복마다 합쳐 통계(평행군은 두 군 합침)
UNITS <- V$rule_units
unit_rep <- function(m, md, arms) stat_of(NCA[pk_model == m & wmode == md & arm %in% arms], "rep")
UT <- rbindlist(lapply(MODELS, function(m) rbindlist(lapply(names(UNITS), function(u) {
  arms <- unlist(UNITS[[u]]$arms); md_set <- if (any(vapply(arms, function(aid) is_assumed(V$arms[[aid]]), TRUE))) MODES else "primary"
  studies <- unique(vapply(arms, function(aid) V$arms[[aid]]$study, ""))
  rbindlist(lapply(md_set, function(md) four_rows(unit_rep(m, md, arms), obs_vec(lapply(arms, ob_arm)))[,
    `:=`(pk_model = m, unit = u, wmode = md, arms = paste(arms, collapse = ";"), n_subj = sum(vapply(arms, function(aid) as.numeric(V$arms[[aid]]$n_subj), 1)),
         dev_status = dev_status(studies, m))]))
}))))
setcolorder(UT, c("pk_model", "unit", "wmode", "arms", "n_subj", "dev_status", "metric")); fwrite(UT, file.path(out_dir, "validation_units.csv"))

# 비아시아 합동군(기술 통계): 같은 반복의 구성 군을 합침(PKM12350 시험군·PKM14161 시험군과 같은 모의 대상자)
PA <- rbindlist(lapply(names(V$pooled_arms), function(pid) { pa <- V$pooled_arms[[pid]]; arms <- unlist(pa$from); o <- pa$observed
  rbindlist(lapply(MODELS, function(m) four_rows(unit_rep(m, "primary", arms),
    c(AUClast_mean = as.numeric(o$AUClast_mean), AUClast_cv = 100 * as.numeric(o$AUClast_sd) / as.numeric(o$AUClast_mean), AUClast_cv_within = NA_real_,
      Cmax_mean = as.numeric(o$Cmax_mean), Cmax_cv = 100 * as.numeric(o$Cmax_sd) / as.numeric(o$Cmax_mean), Cmax_cv_within = NA_real_))[, `:=`(pk_model = m, pooled_arm = pid, from = paste(arms, collapse = ";"),
      n_subj = sum(vapply(arms, function(aid) as.numeric(V$arms[[aid]]$n_subj), 1)), role = "descriptive", composition_note = pa$composition_note)])) }))
fwrite(PA, file.path(out_dir, "validation_pooled.csv"))

# ---------------------------------------------------------------- 평행군 GMR과 연구 내 변동(기술 통계) ----------------------------------------------------------------
# 예측: (i) 대상자 로그 pooled t GMR, (ii) 관측과 같은 추정량: 반복의 군 산술 평균·SD를 GM = mean / sqrt(1 + CV^2)로 바꾼 GMR.
# 연구 내 합동 로그 분산 s2_w = sum((n_i - 1) log(1 + CV_i^2)) / sum(n_i - 1)도 같은 반복에서.
gm_of <- function(mu, cvp) mu / sqrt(1 + (cvp / 100)^2)
GMR <- rbindlist(lapply(V$parallel_studies, function(ps) rbindlist(lapply(MODELS, function(m) rbindlist(lapply(MODES, function(md) {
  if (!nrow(REP[pk_model == m & arm == ps$T & wmode == md])) return(NULL)
  d <- NCA[pk_model == m & wmode == md & arm %in% c(ps$T, ps$R)][, grp := fifelse(arm == ps$T, "T", "R")]
  pt <- d[is.finite(AUClast) & AUClast > 0 & is.finite(Cmax) & Cmax > 0, .(mT_a = mean(log(AUClast[grp == "T"])) - mean(log(AUClast[grp == "R"])),
                                                                          mT_c = mean(log(Cmax[grp == "T"])) - mean(log(Cmax[grp == "R"]))), by = rep]
  sT <- REP[pk_model == m & arm == ps$T & wmode == md]; sR <- REP[pk_model == m & arm == ps$R & wmode == md]; setkey(sT, rep); setkey(sR, rep)
  nT <- as.numeric(V$arms[[ps$T]]$n_subj); nR <- as.numeric(V$arms[[ps$R]]$n_subj)
  mom_a <- gm_of(sT$AUClast_mean, sT$AUClast_cv) / gm_of(sR$AUClast_mean, sR$AUClast_cv); mom_c <- gm_of(sT$Cmax_mean, sT$Cmax_cv) / gm_of(sR$Cmax_mean, sR$Cmax_cv)
  s2w <- ((nT - 1) * log(1 + (sT$AUClast_cv / 100)^2) + (nR - 1) * log(1 + (sR$AUClast_cv / 100)^2)) / (nT + nR - 2)
  oT <- ob_arm(ps$T); oR <- ob_arm(ps$R)
  o_mom_a <- gm_of(oT$AUClast_mean, 100 * oT$AUClast_sd / oT$AUClast_mean) / gm_of(oR$AUClast_mean, 100 * oR$AUClast_sd / oR$AUClast_mean)
  o_mom_c <- if (is.na(oT$Cmax_mean) || is.na(oR$Cmax_mean)) NA_real_ else gm_of(oT$Cmax_mean, 100 * oT$Cmax_sd / oT$Cmax_mean) / gm_of(oR$Cmax_mean, 100 * oR$Cmax_sd / oR$Cmax_mean)
  o_s2w <- ((nT - 1) * log(1 + (oT$AUClast_sd / oT$AUClast_mean)^2) + (nR - 1) * log(1 + (oR$AUClast_sd / oR$AUClast_mean)^2)) / (nT + nR - 2)
  rg <- if (!is.null(ps$reported_gmr)) as.numeric(ps$reported_gmr$AUClast) else NA_real_
  rbind(cbind(data.table(metric = "AUClast_gmr_moment"), pi_row(mom_a, o_mom_a)), cbind(data.table(metric = "Cmax_gmr_moment"), pi_row(mom_c, o_mom_c)),
        cbind(data.table(metric = "AUClast_gmr_pooled_t_vs_reported"), pi_row(exp(pt$mT_a), rg)), cbind(data.table(metric = "Cmax_gmr_pooled_t"), pi_row(exp(pt$mT_c), NA_real_)),
        cbind(data.table(metric = "AUClast_within_study_log_var"), pi_row(s2w, o_s2w)))[, `:=`(pk_model = m, study = ps$study, wmode = md, role = "descriptive")]
}))))))
setcolorder(GMR, c("pk_model", "study", "wmode", "role", "metric")); fwrite(GMR, file.path(out_dir, "validation_gmr.csv"))

# ---------------------------------------------------------------- 실패 규칙 R1, R2, R4, R5 ----------------------------------------------------------------
# 규칙마다 등록된 예측구간 수준(R4 99%, 나머지 95%)으로 단위별 구간을 반복에서 다시 구한다.
# 오경보율: 각 반복을 차례로 '관측값'으로 놓고 같은 구간으로 판정한 비율(모델과 가정이 맞을 때).
# 검출 확률: 각 반복의 평균(AUC0-last, Cmax)에 배율 f(0.75, 0.85, 1.15, 1.25)를 곱해 '관측값'으로 놓았을 때 규칙이 충족되는 비율(CV는 그대로).
fr <- V$failure_rule; DET <- as.numeric(unlist(fr$detection$factors))
unit_set <- function(spec, m) {
  if (identical(spec, "all")) return(names(UNITS))
  if (identical(spec, "external")) return(UT[pk_model == m & wmode == "primary" & metric == "AUClast_mean" & dev_status == "external", unit])
  unlist(spec)
}
rule_eval <- function(m, md) {
  # 단위별 반복 통계: 체중 범위가 알려진 단위(TDU12265)는 모드와 무관하게 primary
  ur <- rbindlist(lapply(names(UNITS), function(un) { arms <- unlist(UNITS[[un]]$arms); md_ <- if (md %in% UT[pk_model == m & unit == un, wmode]) md else "primary"
    unit_rep(m, md_, arms)[, unit := un] }))
  if (!md %in% FULL_MODES) ur <- ur[rep <= REPS_S]
  ob <- rbindlist(lapply(names(UNITS), function(un) { o <- obs_vec(lapply(unlist(UNITS[[un]]$arms), ob_arm)); data.table(unit = un, metric = names(o), obs = unname(o)) }))
  defs <- list(R1 = fr$R1, R2 = fr$R2, R4 = fr$R4, R5 = fr$R5)
  ev <- lapply(names(defs), function(rn) {
    d <- defs[[rn]]; us <- unit_set(d$units, m); thr <- d$outside_pi_min; thr <- as.integer(if (is.list(thr)) thr[[m]] else thr)
    of <- d$of; of <- as.integer(if (is.list(of)) of[[m]] else of); lev <- if (is.null(d$pi_level)) PI else as.numeric(d$pi_level); ql <- (1 - lev) / 2
    z <- ur[unit %in% us, c("unit", "rep", d$metric), with = FALSE]; setnames(z, d$metric, "v")
    lim <- z[, .(lo = quantile(v, ql, names = FALSE, na.rm = TRUE), hi = quantile(v, 1 - ql, names = FALSE, na.rm = TRUE)), by = unit]
    o <- merge(lim, ob[metric == d$metric & unit %in% us], by = "unit")[, out := obs < lo | obs > hi]
    if (nrow(o) != of || anyNA(o$out)) stop(sprintf("%s %s: 비교 가능한 단위 %d != 등록 %d", m, rn, sum(!is.na(o$out)), of))
    z <- merge(z, lim, by = "unit")
    is_mean <- grepl("_mean$", d$metric)
    k_of <- function(f) z[, .(k = sum((if (is_mean) f * v else v) < lo | (if (is_mean) f * v else v) > hi)), by = rep]
    fa <- k_of(1)
    det <- vapply(DET, function(f) mean(k_of(f)$k >= thr), 1)
    list(row = data.table(rule = rn, metric = d$metric, pi_level = lev, units = paste(us, collapse = ";"), units_outside = paste(o[out == TRUE, unit], collapse = ";"),
                          n_outside = sum(o$out), n_compared = nrow(o), threshold = thr, met = sum(o$out) >= thr, false_alarm_sim_pct = 100 * mean(fa$k >= thr)),
         fa_rep = fa[k >= thr, rep], det = setNames(det, paste0("detect_x", formatC(DET, format = "f", digits = 2), "_pct")),
         det_rep = lapply(DET, function(f) { kk <- k_of(f); kk[k >= thr, rep] }))
  })
  rows <- rbindlist(lapply(seq_along(ev), function(i) cbind(ev[[i]]$row, as.data.table(as.list(100 * ev[[i]]$det)))))
  nrep <- max(ur$rep)
  any_fa <- 100 * length(unique(unlist(lapply(ev, `[[`, "fa_rep")))) / nrep
  any_det <- vapply(seq_along(DET), function(j) 100 * length(unique(unlist(lapply(ev, function(e) e$det_rep[[j]])))) / nrep, 1)
  anyrow <- data.table(rule = "any", metric = "", pi_level = NA_real_, units = "", units_outside = "", n_outside = NA_integer_, n_compared = NA_integer_, threshold = NA_integer_,
                       met = any(rows$met), false_alarm_sim_pct = any_fa)
  anyrow <- cbind(anyrow, as.data.table(as.list(setNames(any_det, names(ev[[1]]$det)))))
  out <- rbind(rows, anyrow)
  out[, `:=`(pk_model = m, wmode = md, failure = any(rows$met), n_rep = nrep)][]
}
RULE <- rbindlist(lapply(MODELS, function(m) rbindlist(lapply(FULL_MODES, function(md) rule_eval(m, md)))))
# 결과 범주: 두 가정 모두 실패 / 고정 가정에서만 실패 / 평균 불확실 가정에서만 실패 / 두 가정 모두 통과
RULE[, outcome := { f1 <- failure[wmode == "primary"][1]; f2 <- failure[wmode == "wunc"][1]
  if (f1 && f2) "fails under both weight assumptions" else if (f1) "fails under the fixed weight assumption only" else if (f2) "fails under the uncertain-mean weight assumption only" else "passes under both weight assumptions" }, by = pk_model]
setcolorder(RULE, c("pk_model", "wmode", "rule")); fwrite(RULE, file.path(out_dir, "validation_rule.csv"))
RULE_WS <- rbindlist(lapply(MODELS, function(m) rbindlist(lapply(setdiff(MODES, FULL_MODES), function(md) rule_eval(m, md)))))
fwrite(RULE_WS, file.path(out_dir, "validation_rule_weight_sensitivity.csv"))
# 이전 Step 1 군별 결과(규칙 작성 전에 알려져 있던 값) 나란히
S1F <- rbindlist(lapply(c(k2016 = "step1", k2020 = "step1_k2020"), function(d) { f <- proj_path("results", d, "step1f_300mg_arm_check.csv"); if (file.exists(f)) fread(f) }), idcol = "pk_model", fill = TRUE)
if (nrow(S1F)) fwrite(S1F, file.path(out_dir, "step1_arm_check_prior.csv"))
print(cov_tab); print(UT[wmode == "primary"]); print(PA); print(GMR[wmode == "primary"]); print(RULE); print(RULE_WS[rule == "any"])

# ---------------------------------------------------------------- 한국어 요약(내부 보고용; 수치는 위 표에서) ----------------------------------------------------------------
f1 <- function(x, d = 1) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d, big.mark = ","))
ML <- c(k2020 = "Kovalenko 2020 Model 1(주 모델)", k2016 = "Kovalenko 2016(민감도 모델)")
L <- c("# [1] 모델 확인과 외부 검증 (D-068, config/prereg_20261020_m15.yaml validation)", "",
       sprintf("반복 %s회(체중 고정 민감도 %s회), %d%% 예측구간. 관측값은 저장소에 있는 값만 썼다(문헌 사이트 접근 불가). 수치는 results/m15/validation/*.csv에서 자동 생성.",
               format(REPS, big.mark = ","), format(REPS_S, big.mark = ","), round(100 * PI)), "",
       "## 공변량 입력", "", "- 두 모델 모두 체중(→Vc)만 공변량이다. EASI·알부민·인종·BMI는 모델 입력이 아니어서 환자 EASI 값이 들어가는 곳이 없다(플래그 없음).", "",
       "## 판정(단위 5개: TDU12265 300·600 mg, PKM12350, PKM14161, R668-HV-1108 두 군 합침)", "")
WL <- c(primary = "고정 체중 가정(B14)", wunc = "평균 불확실 체중 가정")
for (m in MODELS) for (md in FULL_MODES) { r <- RULE[pk_model == m & wmode == md]; g <- function(rn, col) r[rule == rn][[col]]
  L <- c(L, sprintf("- %s, %s: R1(평균, 95%%) %d/%d, R2(군 안 CV, 95%%) %d/%d, R4(외부 단위, 99%%) %d/%d, R5(TDU Cmax, 95%%) %d/%d 구간 밖 → %s (전체 오경보율 모의 %.1f%%, 25%% 편향 검출 확률 %.0f%%/%.0f%%)",
                    ML[[m]], WL[[md]], g("R1", "n_outside"), g("R1", "n_compared"), g("R2", "n_outside"), g("R2", "n_compared"), g("R4", "n_outside"), g("R4", "n_compared"),
                    g("R5", "n_outside"), g("R5", "n_compared"), if (r$failure[1]) "재현 실패(규칙 충족)" else "재현 실패 아님", r[rule == "any", false_alarm_sim_pct],
                    r[rule == "any", detect_x0.75_pct], r[rule == "any", detect_x1.25_pct])) }
for (m in MODELS) L <- c(L, sprintf("- %s 결과 범주: %s", ML[[m]], RULE[pk_model == m, outcome][1]))
L <- c(L, "", "## 판정 단위 AUC0-last 평균(관측 / 예측 중앙값 [95% 예측구간], 백분위 순위)", "")
for (m in MODELS) for (un in names(UNITS)) { a <- UT[pk_model == m & unit == un & wmode == "primary" & metric == "AUClast_mean"]
  L <- c(L, sprintf("- %s %s (n %s, %s): 관측 %s / 예측 %s [%s~%s], 순위 %s%%%s%s", ML[[m]], un, a$n_subj, a$dev_status, f1(a$observed, 0), f1(a$pred_median, 0),
                    f1(a$pi_lo, 0), f1(a$pi_hi, 0), f1(a$pct_rank, 1), if (isTRUE(a$outside_pi)) " 밖" else "", if (isTRUE(a$borderline)) " (경계)" else "")) }
L <- c(L, "", "## 체중 가정 민감도(고정 평균 체중; 판정 충족 여부)", "")
for (m in MODELS) for (md in setdiff(MODES, FULL_MODES)) { r <- RULE_WS[pk_model == m & wmode == md]; L <- c(L, sprintf("- %s, 평균 %s kg 고정: %s", ML[[m]], sub("^w", "", md), if (r$failure[1]) "규칙 충족" else "규칙 미충족")) }
writeLines(L, file.path(out_dir, "validation_summary_ko.md"))
append_run_log(lf, sprintf("done: k2020 %s; k2016 %s", RULE[pk_model == "k2020", outcome][1], RULE[pk_model == "k2016", outcome][1]))
