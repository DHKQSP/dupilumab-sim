#!/usr/bin/env Rscript
# 교차검증: 검토자 독립 구현 기준값(지시서 §8, config/crossval_reference.yaml)과 비교. 채혈 편차 없음. 결과: results/crossval/
# 사용법: Rscript scripts/13_crossval.R [n_trials] [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); ref <- read_cfg("crossval_reference.yaml")
n_trials <- if (length(args) >= 1) as.integer(args[1]) else design$mc$n_trials
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores())
MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "crossval"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("crossval", master_seed = MASTER_SEED, run_mode = "final", extra = list(n_trials = n_trials))
p <- load_params("k2016")

# 1) 개인 수준 B0, 체중 60–90 (평균 75, SD 9), 20,000명, 편차 없음
wt <- weight_spec_from_design(design, "base")
pop <- run_individual_population(20000, p, design, "B0", MASTER_SEED, wt, jitter = FALSE, tag = "xv_ind")
s <- summarize_individual(pop$nca, last_planned = 56)
r <- ref$individual_B0_wt60_90_300mg_n20000
ind_cmp <- data.table(metric = c("AUClast_geo", "AUClast_logcv_pct", "Cmax_logcv_pct", "tlast_median", "tlast_p05", "tlast_p95", "quant_at_day56_pct", "lambda_ok_pct", "reliable_pct", "n_lambda_median", "extrap_median_pct", "extrap_p95_pct", "extrap_gt20_pct", "extrap_true_median_pct", "extrap_true_p95_pct", "coverage_lt80_pct"),
                      sim = c(s$AUClast_geo, s$AUClast_logcv, s$Cmax_logcv, s$tlast_median, s$tlast_p05, s$tlast_p95, s$quant_at_last_pct, s$lambda_ok_pct, s$reliable_pct, s$n_lambda_median, s$extrap_median, s$extrap_p95, s$extrap_gt20_pct, s$extrap_true_median, s$extrap_true_p95, s$coverage_lt80_pct),
                      ref = c(r$AUClast_geo, r$AUClast_logcv_pct, r$Cmax_logcv_pct, r$tlast_median_postdose_day, r$tlast_p05_postdose_day, r$tlast_p95_postdose_day, r$quant_at_day56_pct, r$lambda_ok_pct, r$reliable_pct, r$n_lambda_median, r$extrap_median_pct, r$extrap_p95_pct, r$extrap_gt20_pct, r$extrap_true_median_pct, r$extrap_true_p95_pct, r$coverage_lt80_pct))
ind_cmp[, rel_diff_pct := 100 * (sim / ref - 1)]
fwrite(ind_cmp, file.path(out_dir, "crossval_individual_B0.csv")); cat("\n개인 수준 교차검증:\n"); print(ind_cmp)
# 체중 50–115 vs 60–90 CV. 검토자 50–115 분포 형태는 미전달 → 형태별로 산출(D-019)
wt_wide <- list(mean = 82.5, sd = 1e3, trunc = c(50, 115))     # 거의 균등
wt_shapes <- list(uniform_50_115 = wt_wide, normal78_sd13_50_115 = list(mean = 78, sd = 13, trunc = c(50, 115)),
                  normal78_sd10_50_115 = list(mean = 78, sd = 10, trunc = c(50, 115)), normal75_sd9_60_90 = weight_spec_from_design(design, "base"))
cv_cmp <- rbindlist(lapply(names(wt_shapes), function(nm) {
  pp <- run_individual_population(20000, p, design, "B0", MASTER_SEED, wt_shapes[[nm]], jitter = FALSE, tag = paste0("xv_wt_", nm))
  data.table(weight_shape = nm, sim_logcv = log_cv_pct(pp$nca$AUClast), wt_mean = mean(pp$subj$WT), wt_sd = sd(pp$subj$WT))
}))
cv_cmp[, ref_logcv := fifelse(grepl("50_115", weight_shape), ref$weight_range_effect$logcv_wt50_115_pct, ref$weight_range_effect$logcv_wt60_90_pct)]
fwrite(cv_cmp, file.path(out_dir, "crossval_weight_range_cv.csv")); cat("\n체중 분포 형태별 AUClast log-CV:\n"); print(cv_cmp)
# 밀도 예상 방향: Day 46·53(경과일 45·52) 추가. 기준 일정이 B0인지 이전 일정인지 미전달 → 둘 다, 체중 균등·정규 두 형태
dens <- rbindlist(lapply(c("uniform_50_115", "normal78_sd13_50_115"), function(nm) {
  pd <- run_individual_population(20000, p, design, c("B0", "XV_B0_dens", "XV_prev", "XV_prev_dens"), MASTER_SEED, wt_shapes[[nm]], jitter = FALSE, tag = paste0("xv_dens_", nm))
  rbindlist(lapply(c("B0", "XV_B0_dens", "XV_prev", "XV_prev_dens"), function(sh) summarize_individual(pd$nca[schedule == sh], last_planned = 56)[, `:=`(schedule = sh, weight_shape = nm)]))
}))
fwrite(dens[, .(weight_shape, schedule, extrap_median, extrap_true_median, reliable_pct, lambda_ok_pct, extrap_gt20_pct)], file.path(out_dir, "crossval_density_direction.csv"))
cat("\n밀도 방향(검토자: 외삽 2.9→2.3, 참값 0.69→0.46, 신뢰 87.8→87.9):\n"); print(dens[, .(weight_shape, schedule, extrap_median, extrap_true_median, reliable_pct)])

# 2) 시험 수준: XV_prev 일정, 체중 50–115, 500회, 117명/arm, S00/F090/VM125/KM10. 체중 형태 두 가지(D-019)
scens <- sc$scenarios[c("S00", "F090", "VM125", "KM10")]
combos <- CJ(scenario = names(scens), schedule = "XV_prev")
rt <- ref$trial_XVprev_wt50_115_n117_500trials
all_cmp <- list(); all_tr <- list()
for (shape in c("normal78_sd13_50_115", "uniform_50_115")) {
  t0 <- Sys.time()
  res <- run_trials(n_trials, p, design, scens, combos, MASTER_SEED, wt_shapes[[shape]], jitter = FALSE, methods = "pooled_t", cores = cores)
  cat(shape, "trial-level elapsed:", format(Sys.time() - t0), "\n")
  st <- summarize_trials(res$be)
  tr <- st$per_endpoint[endpoint %in% c("AUClast", "AUCinf_reliable", "AUCinf_true", "Cmax"), .(scenario, endpoint, pass_rate, GMR_mean, width_mean_pp)][, weight_shape := shape]
  all_tr[[shape]] <- tr
  fwrite(st$concordance[, weight_shape := shape], file.path(out_dir, sprintf("crossval_trial_XVprev_concordance_%s.csv", shape)))
  g <- function(sc_, ep, col) tr[scenario == sc_ & endpoint == ep][[col]]
  all_cmp[[shape]] <- data.table(weight_shape = shape,
    scenario = c("S00", "S00", "S00", "F090", "F090", "VM125", "VM125", "VM125", "VM125", "VM125"),
    metric = c("pass_AUClast", "pass_AUCinf_reliable", "width_AUClast_pp", "pass_AUClast", "pass_AUCinf_reliable", "pass_AUClast", "pass_AUCinf_reliable", "gmr_AUCinf_true", "gmr_AUClast", "gmr_AUCinf_reliable"),
    sim = c(g("S00", "AUClast", "pass_rate"), g("S00", "AUCinf_reliable", "pass_rate"), g("S00", "AUClast", "width_mean_pp"), g("F090", "AUClast", "pass_rate"), g("F090", "AUCinf_reliable", "pass_rate"),
            g("VM125", "AUClast", "pass_rate"), g("VM125", "AUCinf_reliable", "pass_rate"), g("VM125", "AUCinf_true", "GMR_mean"), g("VM125", "AUClast", "GMR_mean"), g("VM125", "AUCinf_reliable", "GMR_mean")),
    ref = c(rt$S00$pass_AUClast_pct, rt$S00$pass_AUCinf_reliable_pct, rt$S00$ci_width_AUClast_pp, rt$F090$pass_AUClast_pct, rt$F090$pass_AUCinf_reliable_pct,
            rt$VM125$pass_AUClast_pct, rt$VM125$pass_AUCinf_reliable_pct, rt$VM125$gmr_AUCinf_true, rt$VM125$gmr_AUClast, rt$VM125$gmr_AUCinf_reliable))
  # 통과율의 몬테카를로 표준오차(500회): sqrt(p(1-p)/500)
  all_cmp[[shape]][grepl("^pass", metric), mc_se_pp := 100 * sqrt((ref / 100) * (1 - ref / 100) / n_trials)]
  all_cmp[[shape]][, diff := sim - ref]
}
tr_all <- rbindlist(all_tr); cmp <- rbindlist(all_cmp)
fwrite(tr_all, file.path(out_dir, "crossval_trial_XVprev.csv")); fwrite(cmp, file.path(out_dir, "crossval_trial_compare.csv"))
cat("\n시험 수준 교차검증:\n"); print(cmp)
km <- dcast(tr_all[scenario %in% c("S00", "KM10")], weight_shape + endpoint ~ scenario, value.var = "GMR_mean")[, rel_pct := 100 * (KM10 / S00 - 1)]
fwrite(km, file.path(out_dir, "crossval_km10_vs_s00.csv")); cat("\nKM10 vs S00 GMR (기준: 전 지표 ±1.1% 이내):\n"); print(km)
append_run_log(logfile, "done")
