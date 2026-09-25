#!/usr/bin/env Rscript
# 53_atopic_summary.R — 아토피 피부염 성인 체중 분포 요약 (두 번째 추가 지시 2026-09-26 §2, prereg section5, config/population_atopic.yaml, D-058).
# 입력: scripts/52 산출(results/atopic/: 대상자 수준 NCA rds, atopic_individual_*.csv, atopic_truth.csv, atopic_trials_be/drop_<variant>.csv.gz),
#       저장된 60–90 kg 20,000명 B0 NCA(results/individual/nca_base_20000.rds, nca_struct2020_20000.rds), 비만 스트레스 시험(results/weight_generalization/obese_trials_*).
# 새 모의는 5)의 탐색 점검 1개뿐(사전 등록 밖, 라벨 exploratory): 키와 체중의 상관을 NCT03389893 요약값에서 구해 BMI 공변량 변형만 다시 모의.
# 산출(results/atopic/):
#   atopic_weight_table.csv            분포 3개 × (해석적 절단 전·후, 모의 20,000명, 지시문 기준값): 60 kg 미만, 60–90, 90 초과, 100 초과, 60–90 밖, 평균, SD
#   atopic_population_characteristics.csv  모의 체중·키·BMI 평균·SD 대 출처 요약값
#   atopic_pillar1_by_band.csv         모델 변형 × 분포 × 체중 밴드: 참 외삽 비율, AUClast/참 AUC0-inf(중앙값, 5백분위, 최소), 80% 미만 비율 (+ 60–90 kg 연구 모집단)
#   atopic_criteria_by_band.csv        모델 변형 × 분포 × 밴드 × 세트: 미달 %(Wilson 95%), λz 산출 불가 %
#   atopic_criteria_weight_trend.csv   모델 변형 × 분포 × 세트: 미달의 체중 10 kg당 오즈비(로지스틱, 기술 통계)
#   atopic_coverage_failing.csv        모델 변형 × 분포 × 세트: 미달자에서 AUClast/참 AUC0-inf(중앙값, 5·1백분위, 최소)
#   atopic_exploratory_height_corr.csv 탐색: 키-체중 상관 반영 시 BMI 분포와 세트별 미달(변형 k2016_bmi_vc0817, 주 분포)
#   atopic_trials_bias.csv, atopic_trials_pass.csv, atopic_trials_concordance.csv, atopic_trials_arm_failure.csv   시험 수준(주 분포, 시험 2,000회)
#   fig_atopic_weight.png, fig_atopic_failure_by_band.png, fig_atopic_trials.png, atopic_conclusion_en.md, atopic_conclusion_ko.md
# 사용법: Rscript scripts/53_atopic_summary.R   (환경변수 DUPI_ATOPIC_SUMMARY_OUT, DUPI_ATOPIC_ALLOW_PARTIAL=1은 시험용)
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
pr <- read_cfg("prereg_20260926.yaml")$section5; pa <- read_cfg("population_atopic.yaml"); design <- read_cfg("trial_design.yaml")
in_dir <- proj_path("results", "atopic"); out_dir <- Sys.getenv("DUPI_ATOPIC_SUMMARY_OUT", in_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
allow_partial <- nzchar(Sys.getenv("DUPI_ATOPIC_ALLOW_PARTIAL", ""))
MASTER_SEED <- 20260923L; VARIANTS <- unlist(pr$model_variants); DISTS <- unlist(pr$distributions); SETS <- names(CRIT_SETS)
MODEL_EN <- c(base = "2016 model", struct2020 = "Model 1", k2016_bmi_vc0817 = "2016 model + BMI and weight covariates")
DIST_EN <- c(primary = "primary (78/19 kg)", sens_nct03389893 = "sensitivity (80/19 kg)", sens_75_18 = "sensitivity (75/18 kg)")
BAND_EN <- c(all = "all", "below 60" = "below 60 kg", "60-90" = "60 to 90 kg", "above 90-100" = "above 90 to 100 kg", "above 100" = "above 100 kg", "above 90" = "above 90 kg")
BANDS <- names(BAND_EN)
f0 <- function(x) formatC(round(x + sign(x) * 1e-9) + 0, format = "f", digits = 0); f1 <- function(x) formatC(round(x + sign(x) * 1e-9, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x + sign(x) * 1e-9, 2) + 0, format = "f", digits = 2); f3 <- function(x) formatC(round(x + sign(x) * 1e-9, 3) + 0, format = "f", digits = 3)
rg <- function(x, f = f1, u = "%", sep = " to ") if (identical(f(min(x)), f(max(x)))) sprintf("%s%s", f(min(x)), u) else sprintf("%s%s%s%s%s", f(min(x)), u, sep, f(max(x)), u)
lnp <- function(m, s) { sl <- sqrt(log(1 + (s / m)^2)); c(ml = log(m) - sl^2 / 2, sl = sl) }
rds_of <- function(v, d) readRDS(file.path(in_dir, sprintf("atopic_nca_%s_%s.rds", v, d)))
add_above90 <- function(x) rbind(x, copy(x[band %in% c("above 90-100", "above 100")])[, band := "above 90"])

# ---- 1) 체중 표 ---------------------------------------------------------------------------------------------------------
wt_row <- function(m, s, trunc = NULL) {
  q <- lnp(m, s); P <- function(x) plnorm(x, q[["ml"]], q[["sl"]])
  lo <- if (is.null(trunc)) 0 else P(trunc[1]); hi <- if (is.null(trunc)) 1 else P(trunc[2]); Fc <- function(x) (P(x) - lo) / (hi - lo)
  # 절단 로그정규의 1·2차 적률(해석식)
  a <- if (is.null(trunc)) -Inf else (log(trunc[1]) - q[["ml"]]) / q[["sl"]]; b <- if (is.null(trunc)) Inf else (log(trunc[2]) - q[["ml"]]) / q[["sl"]]
  mom <- function(k) exp(k * q[["ml"]] + k^2 * q[["sl"]]^2 / 2) * (pnorm(b - k * q[["sl"]]) - pnorm(a - k * q[["sl"]])) / (pnorm(b) - pnorm(a))
  data.table(below_60 = 100 * Fc(60), in_60_90 = 100 * (Fc(90) - Fc(60)), above_90 = 100 * (1 - Fc(90)), above_100 = 100 * (1 - Fc(100)),
             outside_60_90 = 100 * (Fc(60) + 1 - Fc(90)), wt_mean = mom(1), wt_sd = sqrt(mom(2) - mom(1)^2))
}
WS <- fread(file.path(in_dir, "atopic_weight_simulated.csv"))
WT_TAB <- rbindlist(lapply(DISTS, function(d) { x <- pa$distributions[[d]]; m <- as.numeric(x$mean); s <- as.numeric(x$sd); tr <- as.numeric(unlist(x$trunc)); sm <- WS[distribution == d]
  rbind(cbind(data.table(distribution = d, source = "analytic, untruncated", n = NA_integer_), wt_row(m, s)),
        cbind(data.table(distribution = d, source = sprintf("analytic, truncated %g to %g kg", tr[1], tr[2]), n = NA_integer_), wt_row(m, s, tr)),
        data.table(distribution = d, source = "simulated (20,000 subjects)", n = sm$n, below_60 = sm$below_60, in_60_90 = sm$in_60_90, above_90 = sm$above_90, above_100 = sm$above_100,
                   outside_60_90 = sm$below_60 + sm$above_90, wt_mean = sm$wt_mean, wt_sd = sm$wt_sd)) }))
ref <- pa$reference_shares_primary_untruncated
WT_TAB <- rbind(WT_TAB, data.table(distribution = "primary", source = "reference (directive, independent calculation, untruncated)", n = NA_integer_, below_60 = as.numeric(ref$below_60),
                                   in_60_90 = 100 - as.numeric(ref$outside_60_90), above_90 = as.numeric(ref$above_90), above_100 = as.numeric(ref$above_100), outside_60_90 = as.numeric(ref$outside_60_90)), fill = TRUE)
chk <- WT_TAB[distribution == "primary" & source == "analytic, untruncated"]
ref_dev <- max(abs(c(chk$below_60 - ref$below_60, chk$above_90 - ref$above_90, chk$above_100 - ref$above_100, chk$outside_60_90 - ref$outside_60_90)))
if (ref_dev > 0.06) stop(sprintf("로그정규 모수화가 지시문 기준값과 다릅니다(최대 차이 %.3f %%p)", ref_dev))   # 기준값은 소수 첫째 자리 반올림
# 모의 대 해석(절단) 비율: 이항 표준오차 z
sim_z <- rbindlist(lapply(DISTS, function(d) { a <- WT_TAB[distribution == d & grepl("^analytic, truncated", source)]; s <- WT_TAB[distribution == d & grepl("^simulated", source)]
  rbindlist(lapply(c("below_60", "in_60_90", "above_90", "above_100"), function(k) { p0 <- a[[k]] / 100; data.table(distribution = d, share = k, z = (s[[k]] / 100 - p0) / sqrt(p0 * (1 - p0) / s$n)) })) }))
if (any(abs(sim_z$z) > 4)) stop("모의 체중 비율이 해석값과 맞지 않습니다: ", paste(sim_z[abs(z) > 4, sprintf("%s %s z=%.1f", distribution, share, z)], collapse = "; "))
WT_TAB[, sim_vs_analytic_max_abs_z := fifelse(grepl("^simulated", source), sim_z[, .(m = max(abs(z))), by = distribution][match(WT_TAB$distribution, distribution), m], NA_real_)]
fwrite(WT_TAB, file.path(out_dir, "atopic_weight_table.csv"))

# 모집단 특성(모의 대 출처)
PC <- rbindlist(lapply(DISTS, function(d) { x <- rds_of(VARIANTS[1], d)
  data.table(distribution = d, source = "simulated (20,000 subjects)", wt_mean = mean(x$WT), wt_sd = sd(x$WT), ht_mean = mean(x$HT), ht_sd = sd(x$HT), bmi_mean = mean(x$BMI), bmi_sd = sd(x$BMI), cor_wt_ht = cor(x$WT, x$HT)) }))
PC <- rbind(PC, data.table(distribution = "source", source = "NCT03389893 (71 US adults with atopic dermatitis)", wt_mean = 80.2, wt_sd = 19.0, ht_mean = 172.5, ht_sd = 10.4, bmi_mean = 26.9, bmi_sd = 5.45, cor_wt_ht = NA_real_),
            data.table(distribution = "source", source = "Kovalenko 2021 (adult phase 3 population PK data)", bmi_sd = 5.47), fill = TRUE)
fwrite(PC, file.path(out_dir, "atopic_population_characteristics.csv"))

# ---- 2) Pillar 1과 기준 미달: 체중 밴드별 -----------------------------------------------------------------------------------
study <- list(base = readRDS(proj_path("results", "individual", "nca_base_20000.rds"))[schedule == "B0"], struct2020 = readRDS(proj_path("results", "individual", "nca_struct2020_20000.rds"))[schedule == "B0"])
p1_of <- function(x) x[, .(n = .N, extrap_true_median = median(pct_extrap_true, na.rm = TRUE), extrap_true_p95 = quantile(pct_extrap_true, 0.95, na.rm = TRUE, names = FALSE),
                           extrap_true_max = max(pct_extrap_true, na.rm = TRUE), coverage_median = median(coverage_true, na.rm = TRUE), coverage_p05 = quantile(coverage_true, 0.05, na.rm = TRUE, names = FALSE),
                           coverage_min = min(coverage_true, na.rm = TRUE), coverage_lt80_pct = 100 * mean(coverage_true < 0.80, na.rm = TRUE), wt_mean = mean(WT)), by = band]
P1 <- list(); CB <- list(); TRD <- list(); CF <- list()
pops <- c(lapply(setNames(nm = VARIANTS), function(v) setNames(lapply(DISTS, function(d) rds_of(v, d)), DISTS)))
for (v in VARIANTS) for (d in DISTS) {
  x <- copy(pops[[v]][[d]]); x[, band := as.character(band)]; stopifnot(nrow(x) == as.integer(pr$individual$n_subjects))
  xa <- add_above90(rbind(copy(x)[, band := "all"], x))
  tag <- list(population = "adult atopic dermatitis", variant = v, distribution = d)
  P1[[paste(v, d)]] <- p1_of(xa)[, names(tag) := tag]
  for (s_ in SETS) {
    ok <- paste0("ok_", s_); stopifnot(identical(x[[ok]], crit_ok(x, s_)))                                           # 52에서 저장한 판정과 재계산 일치
    CB[[paste(v, d, s_)]] <- xa[, { k <- sum(!get(ok)); w <- wilson_ci(k, .N); .(set = s_, n = .N, n_fail = k, fail_pct = w$est, fail_lo = w$lo, fail_hi = w$hi, lz_fail_pct = 100 * mean(!(lambda_ok %in% TRUE))) }, by = band][, names(tag) := tag]
    fit <- glm(!x[[ok]] ~ I(x$WT / 10), family = binomial()); ci <- suppressMessages(confint.default(fit))[2, ]
    TRD[[paste(v, d, s_)]] <- data.table(variant = v, distribution = d, set = s_, or_per_10kg = exp(coef(fit)[[2]]), or_lo = exp(ci[[1]]), or_hi = exp(ci[[2]]))
    # 창 포착(Pillar 1 정의): 참 AUC(0–tlast)/참 AUC0-inf = coverage_true. 관측 AUClast/참 AUC0-inf는 여기에 분석 오차와 사다리꼴 근사가 더해진다(두 가지를 따로 보고)
    cw <- x[!get(ok), coverage_true]; r <- x[!get(ok), AUClast / AUCinf_true]; rr <- x[get(ok), AUClast / AUCinf_true]
    CF[[paste(v, d, s_)]] <- data.table(variant = v, distribution = d, set = s_, n_fail = length(r),
                                        window_fail_median = median(cw), window_fail_p05 = quantile(cw, 0.05, names = FALSE), window_fail_min = min(cw), window_all_min = min(x$coverage_true),
                                        obs_fail_median = median(r), obs_fail_p05 = quantile(r, 0.05, names = FALSE), obs_fail_p01 = quantile(r, 0.01, names = FALSE), obs_fail_min = min(r),
                                        obs_retained_median = median(rr), obs_retained_p05 = quantile(rr, 0.05, names = FALSE), obs_retained_min = min(rr))
  }
}
for (v in names(study)) { x <- copy(study[[v]]); x[, band := "all"]
  P1[[paste("study", v)]] <- p1_of(x)[, `:=`(population = "60-90 kg (study inclusion range)", variant = v, distribution = "study_60_90")] }
P1 <- rbindlist(P1, use.names = TRUE); CB <- rbindlist(CB, use.names = TRUE); TRD <- rbindlist(TRD); CF <- rbindlist(CF)
# 52의 요약과 대조(같은 대상자, 같은 판정)
C52 <- fread(file.path(in_dir, "atopic_individual_criteria.csv"))
cmp <- merge(CB[band != "above 90"], C52[, .(variant, distribution, band, set, fail_pct52 = fail_pct, n52 = n)], by = c("variant", "distribution", "band", "set"))
if (nrow(cmp) != nrow(CB[band != "above 90"]) || any(abs(cmp$fail_pct - cmp$fail_pct52) > 1e-9) || any(cmp$n != cmp$n52)) stop("밴드별 미달 비율이 52 요약과 다릅니다")
P1[, band := factor(band, BANDS)]; CB[, band := factor(band, BANDS)]
setorder(P1, population, variant, distribution, band); setorder(CB, variant, distribution, set, band)
fwrite(P1, file.path(out_dir, "atopic_pillar1_by_band.csv")); fwrite(CB, file.path(out_dir, "atopic_criteria_by_band.csv"))
fwrite(TRD, file.path(out_dir, "atopic_criteria_weight_trend.csv")); fwrite(CF, file.path(out_dir, "atopic_coverage_failing.csv"))

# ---- 3) 탐색 점검(사전 등록 밖): 키-체중 상관 -------------------------------------------------------------------------------
# 지시문의 키 분포는 체중과 독립(정규 172, SD 10, 150–200 cm) → 모의 BMI SD가 출처보다 크다. 상관 ρ는 NCT03389893의 평균·SD에서 로그 척도로 구한다:
# var(log BMI) = var(log W) + 4 var(log H) − 4 ρ sd(log W) sd(log H). BMI를 쓰는 변형(k2016_bmi_vc0817)만 같은 대상자·같은 관측 난수로 다시 모의한다.
sdl <- function(m, s) sqrt(log(1 + (s / m)^2))
rho <- (sdl(80.2, 19.0)^2 + 4 * sdl(172.5, 10.4)^2 - sdl(26.9, 5.45)^2) / (4 * sdl(80.2, 19.0) * sdl(172.5, 10.4))
stopifnot(rho > 0, rho < 1)
v_bmi <- "k2016_bmi_vc0817"; tag_p <- "atopic_primary"; rv <- resolve_variant(v_bmi, design); p <- rv$p
spec_p <- { x <- pa$distributions$primary; list(dist = "lognormal", mean = as.numeric(x$mean), sd = as.numeric(x$sd), trunc = as.numeric(unlist(x$trunc)),
                                                height = list(mean = as.numeric(pa$height$mean), sd = as.numeric(pa$height$sd), trunc = as.numeric(unlist(pa$height$trunc)))) }
subj0 <- with_seed(derive_seed(MASTER_SEED, tag_p, "subj"), make_subjects(as.integer(pr$individual$n_subjects), p, spec_p, as.numeric(pa$sex_ratio_male), p$ada$fraction))
run_b0 <- function(subj) {
  grid <- union_grid(design, "B0")
  obs <- with_seed(derive_seed(MASTER_SEED, tag_p, "jitter"), make_obs_times(subj$id, grid, design, jitter = TRUE))
  eps <- with_seed(derive_seed(MASTER_SEED, tag_p, "eps"), draw_eps(subj$id, sort(unique(c(0, grid))), p$sigma))
  sim <- simulate_observations(individual_params(p, subj), obs, design$dose_mg, study_lloq(), eps, model_id = p$model_id)
  ob <- subset_schedule(sim$obs, get_schedule(design, "B0"))
  merge(attach_truth(run_nca(ob), ob, sim$truth), subj[, .(id, WT, HT, BMI)], by = "id")
}
x_ind <- run_b0(subj0); st <- pops[[v_bmi]]$primary                                                                       # 재현: 저장된 52 결과와 같아야 한다
if (!isTRUE(all.equal(x_ind[order(id), .(WT, AUClast, AUCinf_true, Rsq_adjusted)], st[order(id), .(WT, AUClast, AUCinf_true, Rsq_adjusted)], check.attributes = FALSE))) stop("탐색 점검: 52 결과 재현 실패")
q <- lnp(spec_p$mean, spec_p$sd); h <- spec_p$height
subj1 <- copy(subj0)
with_seed(derive_seed(MASTER_SEED, tag_p, "height_corr"), {
  zw <- (log(subj1$WT) - q[["ml"]]) / q[["sl"]]; mu <- h$mean + h$sd * rho * zw; sdc <- h$sd * sqrt(1 - rho^2)
  u <- runif(nrow(subj1), pnorm(h$trunc[1], mu, sdc), pnorm(h$trunc[2], mu, sdc)); subj1[, HT := qnorm(u, mu, sdc)]
})
subj1[, BMI := WT / (HT / 100)^2]
x_cor <- run_b0(subj1)
EX <- rbindlist(lapply(list(independent = x_ind, correlated = x_cor), function(x) {
  r <- data.table(cor_wt_ht = cor(x$WT, x$HT), ht_sd = sd(x$HT), bmi_mean = mean(x$BMI), bmi_sd = sd(x$BMI), coverage_median = median(x$coverage_true), coverage_min = min(x$coverage_true))
  for (s_ in SETS) r[, (paste0("fail_pct_", s_)) := 100 * mean(!crit_ok(x, s_))]
  r }), idcol = "height_model")
EX[, `:=`(rho_target = rho, variant = v_bmi, distribution = "primary", status = "exploratory, not pre-registered")]
fwrite(EX, file.path(out_dir, "atopic_exploratory_height_corr.csv"))

# ---- 4) 시험 수준 --------------------------------------------------------------------------------------------------------
TRUTH <- if (file.exists(file.path(in_dir, "atopic_truth.csv"))) fread(file.path(in_dir, "atopic_truth.csv")) else NULL
tv <- VARIANTS[file.exists(file.path(in_dir, sprintf("atopic_trials_be_%s.csv.gz", VARIANTS)))]
if (length(tv) < length(VARIANTS) && !allow_partial) stop("시험 산출이 없는 변형: ", paste(setdiff(VARIANTS, tv), collapse = ", "))
NT <- as.integer(pr$trial$trials); SCN <- unlist(pr$trial$scenarios); AMS <- unlist(pr$trial$analysis_models)
G2V <- c(A_i = "AUCinf_Ai", A_ii = "AUCinf_A", A_iii = "AUCinf_Aiii", A_iv = "AUCinf_Aiv", B = "AUCinf_B", C_i = "AUCinf_Ci", C_ii = "AUCinf_C", C_iii = "AUCinf_Ciii", C_iv = "AUCinf_Civ")
BI <- PS <- CC <- AF <- NULL
if (length(tv) && !is.null(TRUTH)) {
  TRL <- rbindlist(lapply(tv, function(v) fread(file.path(in_dir, sprintf("atopic_trials_be_%s.csv.gz", v)))[, variant := v]))
  cnt <- TRL[, .(nt = uniqueN(trial)), by = variant]
  if (any(cnt$nt != NT) && !allow_partial) stop("시험 수 부족: ", paste(cnt[nt != NT, sprintf("%s %d", variant, nt)], collapse = ", "))
  tr <- melt(TRUTH[, .(variant, scenario, AUC = auc_ratio, Cmax = cmax_ratio)], id.vars = c("variant", "scenario"), variable.name = "kind", value.name = "true_ratio")
  TRL[, kind := fifelse(endpoint == "Cmax", "Cmax", "AUC")]
  BI <- TRL[is.finite(est), .(n = .N, gmr = exp(mean(est)), mean_log = mean(est), sd_log = sd(est)), by = .(variant, scenario, analysis_model = model, endpoint, kind)]
  BI <- tr[BI, on = c("variant", "scenario", "kind")]
  BI[, `:=`(bias_pct = 100 * (exp(mean_log - log(true_ratio)) - 1), bias_lo = 100 * (exp(mean_log - log(true_ratio) - 1.96 * sd_log / sqrt(n)) - 1), bias_hi = 100 * (exp(mean_log - log(true_ratio) + 1.96 * sd_log / sqrt(n)) - 1))]
  BI[, not_estimable_pct := 100 * (1 - n / TRL[, uniqueN(trial), by = variant][match(BI$variant, variant), V1])]
  fwrite(BI, file.path(out_dir, "atopic_trials_bias.csv"))
  W <- dcast(TRL[, .(variant, trial, scenario, model, endpoint, ok = pass %in% TRUE)], variant + trial + scenario + model ~ endpoint, value.var = "ok")
  W[, P2 := Cmax & AUClast]; for (g in names(G2V)) W[, (paste0("G2_", g)) := Cmax & get(G2V[[g]])]
  dcols <- c("Cmax", "AUClast", unname(G2V), "AUCinf_true", "P2", paste0("G2_", names(G2V)))
  L <- melt(W, id.vars = c("variant", "trial", "scenario", "model"), measure.vars = dcols, variable.name = "decision", value.name = "ok")
  PS <- L[, { w <- wilson_ci(sum(ok), .N); .(n_trials = .N, pass_pct = w$est, lo = w$lo, hi = w$hi) }, by = .(variant, scenario, analysis_model = model, decision)]
  fwrite(PS, file.path(out_dir, "atopic_trials_pass.csv"))
  CC <- rbindlist(lapply(names(G2V), function(g) W[, { a <- P2; b <- get(paste0("G2_", g)); .(config = g, n_trials = .N, agree_pct = 100 * mean(a == b), p2_only_pct = 100 * mean(a & !b), g2_only_pct = 100 * mean(!a & b)) },
                                                  by = .(variant, scenario, analysis_model = model)]))
  fwrite(CC, file.path(out_dir, "atopic_trials_concordance.csv"))
  DR <- rbindlist(lapply(tv, function(v) fread(file.path(in_dir, sprintf("atopic_trials_drop_%s.csv.gz", v)))[, variant := v]))
  RCOL <- c(i = "n_reliable_i", ii = "n_reliable", iii = "n_reliable_iii", iv = "n_reliable_iv")
  for (s_ in SETS) DR[, (paste0("fail_", s_)) := 100 * (1 - get(RCOL[[s_]]) / n)]
  fc <- paste0("fail_", SETS)
  R_ <- DR[arm == "R", c("variant", "trial", fc), with = FALSE]; setnames(R_, fc, paste0(fc, "_R"))
  T_ <- DR[arm == "T", c("variant", "trial", "scenario", fc), with = FALSE]
  AT <- R_[T_, on = c("variant", "trial")]
  AF <- rbindlist(lapply(SETS, function(s_) AT[, { r <- get(paste0("fail_", s_, "_R")); t <- get(paste0("fail_", s_)); d <- t - r
    .(set = s_, n_trials = .N, fail_R_mean = mean(r), fail_T_mean = mean(t), diff_mean = mean(d), diff_lo = mean(d) - 1.96 * sd(d) / sqrt(.N), diff_hi = mean(d) + 1.96 * sd(d) / sqrt(.N),
      diff_p025 = quantile(d, 0.025, names = FALSE), diff_p975 = quantile(d, 0.975, names = FALSE)) }, by = .(variant, scenario)]))
  fwrite(AF, file.path(out_dir, "atopic_trials_arm_failure.csv"))
}

# ---- 5) 그림 -------------------------------------------------------------------------------------------------------------
xs <- pops[[VARIANTS[1]]]$primary
dens <- rbindlist(lapply(DISTS, function(d) { x <- pa$distributions[[d]]; q <- lnp(as.numeric(x$mean), as.numeric(x$sd)); tr <- as.numeric(unlist(x$trunc)); w <- seq(tr[1], tr[2], by = 0.5)
  data.table(distribution = DIST_EN[[d]], WT = w, dens = dlnorm(w, q[["ml"]], q[["sl"]]) / (plnorm(tr[2], q[["ml"]], q[["sl"]]) - plnorm(tr[1], q[["ml"]], q[["sl"]]))) }))
g1 <- ggplot() + annotate("rect", xmin = 60, xmax = 90, ymin = -Inf, ymax = Inf, fill = "#e8e7e3", alpha = 0.7) +
  geom_histogram(data = xs, aes(WT, after_stat(density)), binwidth = 2.5, boundary = 40, fill = "#9bbfe9", colour = "white", linewidth = 0.2) +
  geom_line(data = dens, aes(WT, dens, colour = distribution), linewidth = 0.6) +
  scale_colour_manual(values = setNames(c("#1c5fb0", "#c4502a", "#3b8a4f"), DIST_EN), name = NULL) +
  annotate("text", x = 75, y = Inf, vjust = 1.5, label = "study inclusion range 60 to 90 kg", size = 3, colour = "#52514e") +
  labs(x = "Body weight (kg)", y = "Density", title = "Adult atopic dermatitis body weight distributions (lognormal, truncated 40 to 180 kg)",
       subtitle = sprintf("Histogram: 20,000 simulated subjects, primary distribution; %s below 60 kg, %s above 90 kg, %s above 100 kg",
                          paste0(f1(WS[distribution == "primary", below_60]), "%"), paste0(f1(WS[distribution == "primary", above_90]), "%"), paste0(f1(WS[distribution == "primary", above_100]), "%"))) +
  coord_cartesian(xlim = c(40, 160)) + theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_atopic_weight.png"), g1, width = 8, height = 4.5, dpi = 150, bg = "white")
fb <- CB[distribution == "primary" & band %in% c("below 60", "60-90", "above 90-100", "above 100")]
fb[, `:=`(vl = factor(MODEL_EN[variant], MODEL_EN), bl = factor(BAND_EN[as.character(band)], BAND_EN[c("below 60", "60-90", "above 90-100", "above 100")]),
          sl = factor(CRIT_LABEL_EN[set], CRIT_LABEL_EN))]
g2 <- ggplot(fb, aes(bl, fail_pct, colour = sl, group = sl)) + geom_line(linewidth = 0.4) + geom_pointrange(aes(ymin = fail_lo, ymax = fail_hi), size = 0.25) + facet_wrap(~vl, nrow = 1) +
  scale_colour_manual(values = setNames(c("#2a78d6", "#8fb8ea", "#c4502a", "#eb9a74"), CRIT_LABEL_EN), name = NULL) + guides(colour = guide_legend(ncol = 2)) +
  labs(x = "Body weight band", y = "Subjects without a reliable AUC0-inf (%)", title = "Share without a reliable AUC0-inf by body weight band, adult atopic dermatitis population",
       subtitle = "Primary weight distribution, planned schedule, 300 mg, 20,000 subjects per model; bars: Wilson 95% interval") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out_dir, "fig_atopic_failure_by_band.png"), g2, width = 10, height = 4.8, dpi = 150, bg = "white")
if (!is.null(BI)) {
  fbi <- BI[analysis_model == "M0" & endpoint %in% c("AUClast", "AUCinf_Ai", "AUCinf_Aiii", "AUCinf_B", "AUCinf_Ci", "AUCinf_Ciii")]
  EPL <- c(AUClast = "AUClast", AUCinf_Ai = "AUC0-inf rule A, set (i)", AUCinf_Aiii = "AUC0-inf rule A, set (iii)", AUCinf_B = "AUC0-inf rule B", AUCinf_Ci = "AUC0-inf rule C, set (i)", AUCinf_Ciii = "AUC0-inf rule C, set (iii)")
  fbi[, `:=`(el = factor(EPL[endpoint], EPL), vl = factor(MODEL_EN[variant], MODEL_EN), sc = factor(scenario, SCN))]
  g3 <- ggplot(fbi, aes(sc, bias_pct, colour = el)) + geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) +
    geom_pointrange(aes(ymin = bias_lo, ymax = bias_hi), position = position_dodge(width = 0.6), size = 0.2) + facet_wrap(~vl, nrow = 1) +
    scale_colour_manual(values = setNames(c("#1c1c1a", "#2a78d6", "#8fb8ea", "#7a5bb5", "#c4502a", "#eb9a74"), EPL), name = NULL) +
    labs(x = "Scenario (test product)", y = "Bias of the GMR against the true ratio (%)", title = "Bias of AUClast and AUC0-inf GMRs, adult atopic dermatitis population",
         subtitle = sprintf("Primary weight distribution, 117 subjects per arm, %s trials per scenario, pooled t-test (M0); bars: Monte Carlo 95%% interval", format(NT, big.mark = ","))) +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
  ggsave(file.path(out_dir, "fig_atopic_trials.png"), g3, width = 10, height = 4.8, dpi = 150, bg = "white")
}

# ---- 6) 결론 문안 --------------------------------------------------------------------------------------------------------
cp <- CB[distribution == "primary"]; cfp <- CF[distribution == "primary"]
PKP <- "base"                                                                                            # 주 PK 모델(2016, D-023)
X <- cp[band == "all" & set == "i", fail_pct]; Y <- cp[band == "all" & set == "iii", fail_pct]
X0 <- cp[band == "all" & set == "i" & variant == PKP, fail_pct]; Y0 <- cp[band == "all" & set == "iii" & variant == PKP, fail_pct]
Zmin <- min(cfp[set == "iii", window_fail_min]); Zmed <- cfp[set == "iii", window_fail_median]; Z <- floor(100 * Zmin)
stopifnot(all(cfp[set == "i", window_fail_min] >= Zmin))                                                # 세트 (i) 미달자 ⊂ 세트 (iii) 미달자
target_en <- c(sprintf("In the adult atopic dermatitis population, a reliable AUC0-inf cannot be obtained in %s%% (flag set (i)) to %s%% (flag set (iii)) of patients, concentrated in heavier patients, whereas AUClast covers at least %d%% of total exposure in the same patients.",
                       f1(X0), f1(Y0), Z), "",
               sprintf("Basis: 2016 model (primary PK model), primary weight distribution, planned schedule, 300 mg, 20,000 simulated patients. Across the three model variants the shares are %s (set (i)) and %s (set (iii)); %d%% is the smallest coverage (below) in any patient failing set (iii) in any of the three variants, rounded down (median %s).",
                       rg(X), rg(Y), Z, rg(100 * Zmed)), "",
               sprintf("Total exposure is the true AUC0-inf; coverage is the true AUC from dose to the last quantifiable sample divided by the true AUC0-inf (the Pillar 1 definition). The observed AUClast also carries assay error and the trapezoidal approximation: in the same failing patients it is a median %s of the true AUC0-inf (5th percentile %s), against a median %s (5th percentile %s) in patients meeting set (iii), so this error is not specific to the patients without a reliable AUC0-inf and largely cancels in the test-to-reference ratio.",
                       rg(100 * cfp[set == "iii", obs_fail_median]), rg(100 * cfp[set == "iii", obs_fail_p05]), rg(100 * cfp[set == "iii", obs_retained_median]), rg(100 * cfp[set == "iii", obs_retained_p05])))
target_ko <- sprintf("아토피 피부염 성인 모집단에서 신뢰할 수 있는 AUC0-inf를 얻지 못하는 환자는 %s%%(기준 세트 (i))에서 %s%%(세트 (iii))이며 무거운 환자에 몰려 있다. 같은 환자에서 채혈 구간(AUClast 구간)은 참 총노출의 %d%% 이상을 포착한다(2016 모델; 모델 변형 3개 범위 (i) %s, (iii) %s).",
                     f1(X0), f1(Y0), Z, rg(X, f1, "%", " ~ "), rg(Y, f1, "%", " ~ "))
band_line <- function(s_) paste(vapply(VARIANTS, function(v) { z <- cp[variant == v & set == s_]
  sprintf("%s: %s below 60 kg, %s at 60 to 90 kg, %s above 90 kg, %s above 100 kg", MODEL_EN[[v]], paste0(f1(z[band == "below 60", fail_pct]), "%"), paste0(f1(z[band == "60-90", fail_pct]), "%"),
          paste0(f1(z[band == "above 90", fail_pct]), "%"), paste0(f1(z[band == "above 100", fail_pct]), "%")) }, ""), collapse = "; ")
wp <- WT_TAB[distribution == "primary"]; wsim <- wp[grepl("^simulated", source)]; wana <- wp[source == "analytic, untruncated"]
p1a <- P1[population == "adult atopic dermatitis" & distribution == "primary" & band == "all"]; p1s <- P1[population != "adult atopic dermatitis"]
ci_all <- fread(file.path(in_dir, "atopic_individual_criteria.csv"))[band == "all"]
en <- c("# Adult atopic dermatitis body weight distribution", "",
  target_en, "",
  "## Weight distribution (placeholder until the phase 3 data are checked)", "",
  sprintf("- Primary distribution: lognormal with mean 78 kg and SD 19 kg, truncated to 40 to 180 kg (sources: Kamal 2022; Kovalenko 2021; NCT03389893). Untruncated: %s below 60 kg, %s above 90 kg, %s above 100 kg, %s outside 60 to 90 kg (the directive's independent reference: 16.5%%, 23.7%%, 12.4%%, 40.2%%; largest difference %s percentage points).",
          paste0(f1(wana$below_60), "%"), paste0(f1(wana$above_90), "%"), paste0(f1(wana$above_100), "%"), paste0(f1(wana$outside_60_90), "%"), f2(ref_dev)),
  sprintf("- Simulated (20,000 subjects, truncated): %s below 60 kg, %s within 60 to 90 kg, %s above 90 kg, %s above 100 kg; mean %s kg, SD %s kg. The two sensitivity distributions give %s outside 60 to 90 kg.",
          paste0(f1(wsim$below_60), "%"), paste0(f1(wsim$in_60_90), "%"), paste0(f1(wsim$above_90), "%"), paste0(f1(wsim$above_100), "%"), f1(wsim$wt_mean), f1(wsim$wt_sd),
          rg(WT_TAB[distribution != "primary" & grepl("^simulated", source), outside_60_90])),
  sprintf("- Height is drawn independently of weight (normal, mean 172 cm, SD 10 cm, 150 to 200 cm), which gives a BMI SD of %s kg/m2 against 5.45 in NCT03389893. With the weight-height correlation implied by NCT03389893 (%s), the BMI SD is %s kg/m2 and the share without a reliable AUC0-inf in the BMI-covariate model changes from %s to %s under set (i) and from %s to %s under set (iii) (exploratory check, not pre-registered).",
          f2(PC[distribution == "primary", bmi_sd]), f2(rho), f2(EX[height_model == "correlated", bmi_sd]), paste0(f2(EX[height_model == "independent", fail_pct_i]), "%"), paste0(f2(EX[height_model == "correlated", fail_pct_i]), "%"),
          paste0(f2(EX[height_model == "independent", fail_pct_iii]), "%"), paste0(f2(EX[height_model == "correlated", fail_pct_iii]), "%")),
  "- The phase 3 weight data (FDA BLA 761055 clinical pharmacology review, EMA Dupixent EPAR 2017) could not be accessed from the analysis environment; the distribution is replaced when the sponsor supplies them.", "",
  "## Pillar 1 (total exposure captured by AUClast), planned schedule, 300 mg, 20,000 subjects per model", "",
  sprintf("- True extrapolated share beyond the last sample: median %s, 95th percentile %s, maximum %s (atopic population, three models); study population 60 to 90 kg: median %s, maximum %s.",
          rg(p1a$extrap_true_median, f2), rg(p1a$extrap_true_p95, f2), rg(p1a$extrap_true_max, f1), rg(p1s$extrap_true_median, f2), rg(p1s$extrap_true_max, f1)),
  sprintf("- Coverage of the true AUC0-inf by the sampling window (true AUC to the last quantifiable sample / true AUC0-inf): median %s, minimum %s; below 80%% in %s of subjects. Above 100 kg: median true extrapolation %s, maximum %s.",
          rg(100 * p1a$coverage_median, f1), rg(100 * p1a$coverage_min, f1), rg(p1a$coverage_lt80_pct, f3), rg(P1[population == "adult atopic dermatitis" & distribution == "primary" & band == "above 100", extrap_true_median], f2),
          rg(P1[population == "adult atopic dermatitis" & distribution == "primary" & band == "above 100", extrap_true_max], f1)), "",
  "## Subjects without a reliable AUC0-inf (primary distribution)", "",
  sprintf("- Set (i) %s; set (ii) %s; set (iii) %s; set (iv) %s. Per arm of 117: %s under set (i), %s under set (iii).", rg(X), rg(cp[band == "all" & set == "ii", fail_pct]), rg(Y), rg(cp[band == "all" & set == "iv", fail_pct]),
          rg(ci_all[distribution == "primary" & set == "i", fail_per_arm], f1, ""), rg(ci_all[distribution == "primary" & set == "iii", fail_per_arm], f1, "")),
  sprintf("- By weight band under set (i): %s.", band_line("i")), sprintf("- By weight band under set (iii): %s.", band_line("iii")),
  sprintf("- Odds of failing per 10 kg of body weight (logistic, descriptive): %s under set (i), %s under set (iii). Failing subjects are %s kg heavier than retained subjects under set (i).",
          rg(TRD[distribution == "primary" & set == "i", or_per_10kg], f2, ""), rg(TRD[distribution == "primary" & set == "iii", or_per_10kg], f2, ""), rg(ci_all[distribution == "primary" & set == "i", wt_diff], f1, "")),
  sprintf("- In failing subjects the sampling window captures a median of %s of the true AUC0-inf, minimum %s, under set (iii) (set (i): median %s, minimum %s). The observed AUClast / true AUC0-inf in failing subjects: median %s, 5th percentile %s under set (i), %s and %s under set (iii) (retained subjects under set (i): median %s, 5th percentile %s).",
          rg(100 * Zmed, f1), rg(100 * cfp[set == "iii", window_fail_min], f1), rg(100 * cfp[set == "i", window_fail_median], f1), rg(100 * cfp[set == "i", window_fail_min], f1),
          rg(100 * cfp[set == "i", obs_fail_median], f1), rg(100 * cfp[set == "i", obs_fail_p05], f1), rg(100 * cfp[set == "iii", obs_fail_median], f1), rg(100 * cfp[set == "iii", obs_fail_p05], f1),
          rg(100 * cfp[set == "i", obs_retained_median], f1), rg(100 * cfp[set == "i", obs_retained_p05], f1)),
  sprintf("- Sensitivity distributions (80/19 and 75/18 kg): set (i) %s, set (iii) %s.", rg(CB[distribution != "primary" & band == "all" & set == "i", fail_pct]), rg(CB[distribution != "primary" & band == "all" & set == "iii", fail_pct])))
if (!is.null(BI)) {
  bl <- BI[analysis_model == "M0"]; cc0 <- CC[analysis_model == "M0"]; af <- AF
  en <- c(en, "", sprintf("## Trial level (primary distribution, 117 subjects per arm, %s trials per scenario; models: %s)", format(NT, big.mark = ","), paste(MODEL_EN[tv], collapse = ", ")), "",
    sprintf("- Bias of the AUClast GMR against the true AUC0-inf ratio: %s; AUC0-inf rule A: %s (set (i)), %s (set (iii)); rule B: %s; rule C: %s (set (i)), %s (set (iii)) (all scenarios, M0).",
            rg(bl[endpoint == "AUClast", bias_pct], f2), rg(bl[endpoint == "AUCinf_Ai", bias_pct], f2), rg(bl[endpoint == "AUCinf_Aiii", bias_pct], f2), rg(bl[endpoint == "AUCinf_B", bias_pct], f2),
            rg(bl[endpoint == "AUCinf_Ci", bias_pct], f2), rg(bl[endpoint == "AUCinf_Ciii", bias_pct], f2)),
    sprintf("- Agreement of the AUClast + Cmax decision with AUC0-inf + Cmax: %s (rule A, set (i)), %s (rule A, set (iii)), %s (rule C, set (i)), %s (rule C, set (iii)) of trials.",
            rg(cc0[config == "A_i", agree_pct]), rg(cc0[config == "A_iii", agree_pct]), rg(cc0[config == "C_i", agree_pct]), rg(cc0[config == "C_iii", agree_pct])),
    sprintf("- Difference between arms (test minus reference) in the share failing set (i): %s at S00, %s at F090, %s at VM125; set (iii): %s at S00, %s at F090, %s at VM125 (mean over trials, percentage points).",
            rg(af[scenario == "S00" & set == "i", diff_mean], f2, ""), rg(af[scenario == "F090" & set == "i", diff_mean], f2, ""), rg(af[scenario == "VM125" & set == "i", diff_mean], f2, ""),
            rg(af[scenario == "S00" & set == "iii", diff_mean], f2, ""), rg(af[scenario == "F090" & set == "iii", diff_mean], f2, ""), rg(af[scenario == "VM125" & set == "iii", diff_mean], f2, "")))
}
ko <- c("# 아토피 피부염 성인 체중 분포", "", target_ko, "",
  sprintf("- 주 분포(로그정규 78/19 kg, 40–180 kg 절단) 절단 전 60 kg 미만 %s, 90 kg 초과 %s, 100 kg 초과 %s, 60–90 밖 %s(지시문 기준값 16.5/23.7/12.4/40.2%%, 최대 차이 %s%%p). 모의 20,000명: 60–90 kg %s, 90 초과 %s, 100 초과 %s.",
          paste0(f1(wana$below_60), "%"), paste0(f1(wana$above_90), "%"), paste0(f1(wana$above_100), "%"), paste0(f1(wana$outside_60_90), "%"), f2(ref_dev), paste0(f1(wsim$in_60_90), "%"), paste0(f1(wsim$above_90), "%"), paste0(f1(wsim$above_100), "%")),
  sprintf("- 키를 체중과 독립으로 뽑아 BMI SD %s(출처 5.45). NCT03389893에서 구한 상관 %s을 넣으면 BMI SD %s, BMI 공변량 모델의 미달 (i) %s → %s, (iii) %s → %s(탐색, 사전 등록 밖).",
          f2(PC[distribution == "primary", bmi_sd]), f2(rho), f2(EX[height_model == "correlated", bmi_sd]), paste0(f2(EX[height_model == "independent", fail_pct_i]), "%"), paste0(f2(EX[height_model == "correlated", fail_pct_i]), "%"),
          paste0(f2(EX[height_model == "independent", fail_pct_iii]), "%"), paste0(f2(EX[height_model == "correlated", fail_pct_iii]), "%")),
  sprintf("- Pillar 1: 참 외삽 중앙값 %s, 최대 %s; AUClast/참 AUC0-inf 최소 %s, 80%% 미만 %s.", rg(p1a$extrap_true_median, f2, "%", " ~ "), rg(p1a$extrap_true_max, f1, "%", " ~ "), rg(100 * p1a$coverage_min, f1, "%", " ~ "), rg(p1a$coverage_lt80_pct, f2, "%", " ~ ")),
  sprintf("- 미달(주 분포): (i) %s, (ii) %s, (iii) %s, (iv) %s. 10 kg당 오즈비 (i) %s, (iii) %s.", rg(X, f1, "%", " ~ "), rg(cp[band == "all" & set == "ii", fail_pct], f1, "%", " ~ "), rg(Y, f1, "%", " ~ "),
          rg(cp[band == "all" & set == "iv", fail_pct], f1, "%", " ~ "), rg(TRD[distribution == "primary" & set == "i", or_per_10kg], f2, "", " ~ "), rg(TRD[distribution == "primary" & set == "iii", or_per_10kg], f2, "", " ~ ")))
if (!is.null(BI)) ko <- c(ko, sprintf("- 시험 수준(M0): AUClast 편향 %s, 규칙 A (i) %s, (iii) %s; P2와 G2 일치 A(i) %s, A(iii) %s; arm 간 미달 차이 VM125 (i) %s%%p, (iii) %s%%p.",
  rg(BI[analysis_model == "M0" & endpoint == "AUClast", bias_pct], f2, "%", " ~ "), rg(BI[analysis_model == "M0" & endpoint == "AUCinf_Ai", bias_pct], f2, "%", " ~ "), rg(BI[analysis_model == "M0" & endpoint == "AUCinf_Aiii", bias_pct], f2, "%", " ~ "),
  rg(CC[analysis_model == "M0" & config == "A_i", agree_pct], f1, "%", " ~ "), rg(CC[analysis_model == "M0" & config == "A_iii", agree_pct], f1, "%", " ~ "),
  rg(AF[scenario == "VM125" & set == "i", diff_mean], f2, "", " ~ "), rg(AF[scenario == "VM125" & set == "iii", diff_mean], f2, "", " ~ ")))
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", paste(en, collapse = "\n"))) stop("영문 결론 규칙 위반")
writeLines(en, file.path(out_dir, "atopic_conclusion_en.md")); writeLines(ko, file.path(out_dir, "atopic_conclusion_ko.md"))
cat(en, sep = "\n")
