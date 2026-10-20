#!/usr/bin/env Rscript
# 74_m15_samplesize.R — [7] 표본 수: 모의 검정력(8.1절 OC 결과), 분석식 검정력(R/samplesize.R ss_p2_power, 변경 없음), CV 불확실성 assurance,
#   프로토콜 가정(CV 43%, 평가 가능 117명/무작위배정 130명) 대조. 계획서 8.5절, 사전 등록 config/prereg_20261020_m15.yaml samplesize (D-068).
# 입력: results/m15/individual/ss_inputs_m15.csv(72), results/m15/trials/oc/*(73 oc), results/m15/validation/validation_reps.csv.gz(70)
# 산출: results/m15/samplesize/ — power_sim.csv(모의 검정력, Wilson CI), power_analytic.csv(분석식 검정력·assurance), assurance_cv_inputs.csv, syneos_comparison.csv(필요 표본 수 포함), weight_moments.csv
# 사용법: nice -n 5 Rscript scripts/74_m15_samplesize.R [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE); cores <- if (length(args)) as.integer(args[1]) else 4L
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); SSC <- CFG$samplesize
POPS <- CFG$populations; INCL <- CFG$inclusion; RHO <- as.numeric(POPS$weight_height_rho); SR <- as.numeric(POPS$sex_ratio_male)
NGRID <- as.integer(unlist(CFG$n_randomized_per_arm)); EVAL <- 1 - as.numeric(CFG$dropout$p)
od <- proj_path(CFG$output_dir, "samplesize"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
lf <- start_run_log("m15_samplesize", master_seed = MS, run_mode = "final", extra = list(cores = cores))
GA <- c(0.95, 1.00); AM <- c("M0", "M1")

# ---------------------------------------------------------------- 1. 모의 검정력(8.1절) ----------------------------------------------------------------
read_oc <- function(dir) {
  f <- list.files(dir, pattern = "^oc_\\d+_\\d+\\.csv\\.gz$", full.names = TRUE); if (!length(f)) return(NULL)
  x <- rbindlist(lapply(f, fread)); v <- fread(file.path(dir, "versions.csv")); list(x = x, v = v)
}
ocd <- list.dirs(proj_path(CFG$output_dir, "trials", "oc"), recursive = FALSE)
PS <- rbindlist(lapply(ocd, function(dd) {
  r <- read_oc(dd); if (is.null(r)) return(NULL); lab <- strsplit(basename(dd), "_")[[1]]
  vi <- r$v[g %in% GA]; if (nrow(vi) != length(GA)) return(NULL)
  rbindlist(lapply(seq_len(nrow(vi)), function(i) r$x[, .(n_trials = .N, n_pass = sum(bitwAnd(as.integer(mask), as.integer(2^(vi$vi[i] - 1))) > 0)),
                                                       by = .(boundary, schedule, n_rand, model)][, true_gmr := vi$g[i]]))[,
    `:=`(pk_model = lab[1], population = paste(lab[2:(length(lab) - 1)], collapse = "_"), dose_mg = as.numeric(lab[length(lab)]))]
}))
if (nrow(PS)) {
  PS[, c("power_pct", "lo", "hi") := { w <- wilson_ci(n_pass, n_trials); list(w$est, w$lo, w$hi) }]
  setcolorder(PS, c("pk_model", "population", "dose_mg", "boundary", "schedule", "model", "n_rand", "true_gmr")); fwrite(PS, file.path(od, "power_sim.csv"))
}

# ---------------------------------------------------------------- 2. 체중 모멘트(경계별, 등록 체중 200,000명) ----------------------------------------------------------------
BND <- fread(proj_path(CFG$output_dir, "populations", "strata_boundaries.csv"))
WM <- rbindlist(lapply(c("S1", "S2"), function(pv) {
  s <- m15_subjects(200000L, POPS[[pv]], INCL, RHO, SR, derive_seed(MS, "popq"))
  rbindlist(lapply(BND[population == pv, boundary], function(b) { bk <- BND[population == pv & boundary == b, boundary_kg]
    cbind(data.table(population = pv, boundary = b, boundary_kg = bk), m15_weight_moments(s$WT, bk)) }))
}))
fwrite(WM, file.path(od, "weight_moments.csv"))

# ---------------------------------------------------------------- 3. 분석식 검정력과 CV 불확실성 assurance ----------------------------------------------------------------
SSI <- fread(proj_path(CFG$output_dir, "individual", "ss_inputs_m15.csv"))
cells <- merge(SSI, WM, by = "population", allow.cartesian = TRUE)
pw <- function(n_eval, cv, g, inp, am, rule, wm) ss_p2_power(n_eval, cv, g, g, inp, am, rule, wm)
CVG <- seq(20, 80, by = 0.25)                                         # CV 격자(assurance 보간용)
# Li 2020 300 mg 군 AUC0-last CV(저장소 값, 계획서 8.5절 3)
LI <- data.table(arm = c("pkm12350_T", "pkm12350_R", "pkm14161_T", "pkm14161_R", "hv1108_T", "hv1108_R"))
VA <- CFG$validation$arms
LI[, `:=`(n = vapply(arm, function(a) as.numeric(VA[[a]]$n_subj), 1), cv = vapply(arm, function(a) 100 * as.numeric(VA[[a]]$observed$AUClast_sd) / as.numeric(VA[[a]]$observed$AUClast_mean), 1))]
DF <- sum(LI$n - 1); S2 <- sum((LI$n - 1) * log(1 + (LI$cv / 100)^2)) / DF; stopifnot(DF == 96)
# 보정 계수 c: 같은 n의 로그 정규 군 20,000세트에서 산술 CV 기반 합동 추정량의 평균 대 참 로그 분산(관측 합동값에서)
CAL_N <- 20000L
est_cal <- with_seed(derive_seed(MS, "assurance_calibration"), vapply(seq_len(CAL_N), function(r) {
  v <- vapply(LI$n, function(n_) { x <- exp(rnorm(n_, 0, sqrt(S2))); log(1 + (sd(x) / mean(x))^2) }, 1)
  sum((LI$n - 1) * v) / DF }, 1))
CAL <- S2 / mean(est_cal)
draws_raw <- with_seed(derive_seed(MS, "assurance_cv"), DF * S2 / rchisq(as.integer(SSC$assurance_cv$draws), DF))
draws <- CAL * draws_raw
# 원천(Li 2020 군) 체중의 군 안 log 분산: 연구 평균 ~ Uniform(70, 86), SD 10, 기본 범위(남 50~100, 여 40~90), 성비 0.5; 2,000 연구 x 100명
VAL <- CFG$validation; HW <- as.numeric(VAL$weight_mean_uncertainty$half_width_kg); BMv <- as.numeric(unlist(VAL$weight_default_bounds$male)); BFv <- as.numeric(unlist(VAL$weight_default_bounds$female))
VS <- with_seed(derive_seed(MS, "assurance_source_weight"), {
  mu <- runif(2000, 78 - HW, 78 + HW)
  mean(vapply(mu, function(m0) { male <- runif(100) < SR; u <- runif(100); w <- numeric(100)
    w[male] <- m15_q_weight(u[male], list(dist = "truncnorm", mean = m0, sd = 10, lo = -Inf, hi = Inf), BMv[1], BMv[2])
    w[!male] <- m15_q_weight(u[!male], list(dist = "truncnorm", mean = m0, sd = 10, lo = -Inf, hi = Inf), BFv[1], BFv[2])
    var(log(w)) }, 1)) })
n_needed <- function(cv, g, inp, am, rule, wm, target) { for (ne in 20:600) if (pw(ne, cv, g, inp, am, rule, wm) >= target / 100) return(ne); NA_integer_ }
ANL <- parallel::mclapply(seq_len(nrow(cells)), function(i) {
  c_ <- cells[i]; inp <- list(sd_cmax = c_$sd_cmax, ratio_sd_cmax_auc = c_$ratio_sd_cmax_auc, rho_total = c_$rho_total, beta_auc = c_$beta_auc, beta_cmax = c_$beta_cmax)
  wm <- c_[, .(p_stratum2, mean_lwt, var_between, var_within, var_total)]
  s2t <- draws - c_$beta_auc^2 * (VS - c_$var_total)                  # 등록 집단으로 옮긴 로그 분산
  cv_draw <- 100 * sqrt(exp(pmax(s2t, 1e-6)) - 1); clamp <- mean(cv_draw < min(CVG) | cv_draw > max(CVG))
  rbindlist(lapply(NGRID, function(nr) rbindlist(lapply(GA, function(g) rbindlist(lapply(AM, function(am) {
    ne <- EVAL * nr
    pm <- pw(ne, c_$cv_auc_pct, g, inp, am, "fixed", wm)
    curve <- vapply(CVG, function(cv) pw(ne, cv, g, inp, am, "ratio", wm), 1)
    asr <- mean(approx(CVG, curve, xout = pmin(pmax(cv_draw, min(CVG)), max(CVG)))$y)
    data.table(n_rand = nr, n_eval = ne, true_gmr = g, model = am, power_model_cv_pct = 100 * pm, assurance_cv_pct = 100 * asr,
               n_eval_85 = n_needed(c_$cv_auc_pct, g, inp, am, "fixed", wm, 85), n_eval_90 = n_needed(c_$cv_auc_pct, g, inp, am, "fixed", wm, 90),
               power_cv43_pct = 100 * pw(ne, as.numeric(SSC$syneos$auclast_cv_pct), g, inp, am, "ratio", wm), share_draws_clamped = clamp,
               cv_draw_trial_p50 = median(cv_draw))
  }))))))[, `:=`(pk_model = c_$pk_model, population = c_$population, dose_mg = c_$dose_mg, schedule = c_$schedule, boundary = c_$boundary, cv_model_pct = c_$cv_auc_pct)]
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(ANL, function(x) is.null(x) || inherits(x, "try-error"), TRUE); if (length(ANL) != nrow(cells) || any(bad)) stop("분석식 검정력 실패 또는 작업자 결과 없음: ", paste(which(bad), collapse = ","))
AN <- rbindlist(ANL); AN[, `:=`(n_rand_85 = ceiling(n_eval_85 / EVAL), n_rand_90 = ceiling(n_eval_90 / EVAL))]   # 필요 표본 수는 모델 CV(칸마다, n과 무관)
stopifnot(uniqueN(AN[, .(pk_model, population, dose_mg, schedule, boundary)]) == nrow(cells))
setcolorder(AN, c("pk_model", "population", "dose_mg", "schedule", "boundary", "model", "n_rand", "n_eval", "true_gmr"))
fwrite(AN, file.path(od, "power_analytic.csv"))
fwrite(data.table(df = DF, pooled_log_var = S2, pooled_cv_pct = 100 * sqrt(exp(S2) - 1), calibration = CAL, calibration_n_sets = CAL_N, var_logwt_source = VS,
                  cv_draw_raw_p05 = 100 * sqrt(exp(quantile(draws_raw, 0.05, names = FALSE)) - 1), cv_draw_calibrated_p50 = 100 * sqrt(exp(median(draws)) - 1),
                  cv_draw_raw_p95 = 100 * sqrt(exp(quantile(draws_raw, 0.95, names = FALSE)) - 1), n_draws = length(draws)), file.path(od, "assurance_cv_inputs.csv"))

# ---------------------------------------------------------------- 4. 프로토콜 가정 대조와 필요 표본 수 ----------------------------------------------------------------
# [1]에서 재현된 CV: 관측 합동 CV(Li 2020 6군)와 같은 군의 모델 예측 CV(70 반복 중앙값을 같은 자유도로 합동)
VR <- fread(proj_path(CFG$output_dir, "validation", "validation_reps.csv.gz"))[wmode == "primary" & arm %in% LI$arm]
REPCV <- VR[, .(cv_med = median(AUClast_cv)), by = .(pk_model, arm)][LI, on = "arm"][, .(cv_model_li_pct = 100 * sqrt(exp(sum((n - 1) * log(1 + (cv_med / 100)^2)) / sum(n - 1)) - 1)), by = pk_model]
prim <- cells[population == "S1" & dose_mg == as.numeric(CFG$doses_mg$primary) & schedule == "B0" & boundary == "fixed75"]
SY <- rbindlist(lapply(seq_len(nrow(prim)), function(i) {
  c_ <- prim[i]; inp <- list(sd_cmax = c_$sd_cmax, ratio_sd_cmax_auc = c_$ratio_sd_cmax_auc, rho_total = c_$rho_total, beta_auc = c_$beta_auc, beta_cmax = c_$beta_cmax)
  wm <- c_[, .(p_stratum2, mean_lwt, var_between, var_within, var_total)]
  s2_tr <- CAL * S2 - c_$beta_auc^2 * (VS - c_$var_total)
  srcs <- data.table(cv_source = c("protocol assumption (Syneos)", "Li 2020 observed pooled", "Li 2020 observed pooled, calibrated and transported to the enrolled population",
                                   "model prediction for the Li 2020 arms", "model, enrolled S1 population"),
                     cv_pct = c(as.numeric(SSC$syneos$auclast_cv_pct), 100 * sqrt(exp(S2) - 1), 100 * sqrt(exp(s2_tr) - 1), REPCV[pk_model == c_$pk_model, cv_model_li_pct], c_$cv_auc_pct))
  rbindlist(lapply(seq_len(nrow(srcs)), function(k) rbindlist(lapply(GA, function(g) rbindlist(lapply(AM, function(am) rbindlist(lapply(c("ratio", "fixed"), function(rl) {
    cv <- srcs$cv_pct[k]                                                 # ratio: Cmax SD = 비 x AUC0-last SD; fixed: 모델 Cmax SD(저장된 프로토콜 수치의 규칙, D-055)
    data.table(cv_source = srcs$cv_source[k], cv_pct = cv, true_gmr = g, model = am, cmax_rule = rl,
               power_n117_eval_pct = 100 * pw(as.numeric(SSC$syneos$n_evaluable_per_arm), cv, g, inp, am, rl, wm),
               n_eval_85 = n_needed(cv, g, inp, am, rl, wm, 85), n_eval_90 = n_needed(cv, g, inp, am, rl, wm, 90))
  }))))))))[, `:=`(pk_model = c_$pk_model, n_rand_85 = ceiling(n_eval_85 / EVAL), n_rand_90 = ceiling(n_eval_90 / EVAL))]
}))
setcolorder(SY, c("pk_model", "cv_source", "cv_pct", "true_gmr", "model", "cmax_rule")); fwrite(SY, file.path(od, "syneos_comparison.csv"))
print(SY[model == "M1"]); print(AN[pk_model == CFG$models$primary & population == "S1" & dose_mg == 300 & schedule == "B0" & boundary == "fixed75" & model == "M1"])
append_run_log(lf, "done")
