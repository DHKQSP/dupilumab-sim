#!/usr/bin/env Rscript
# §6-1 체중 일반화(Pillar 1): 300 mg, B0, 체중 균등 밴드 × 20,000명, 키 N(170, 9) 절단 150–195 cm → BMI.
# 모델: (a) 2016 주 모델, (b) 2020 Model 1(자체 IIV·잔차), (c) 2016 + ke~BMI 0.368 (기준 26, 자리표시자), (d) (c) + Vc~체중 0.817.
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
JITTER <- !"nojitter" %in% args; args <- setdiff(args, "nojitter")   # nojitter: 채혈 허용창 편차 없이(검토자 독립 구현 조건) 같은 대상자로 재계산
BANDS_ONLY <- grep("^[0-9]+-[0-9]+$", args, value = TRUE); args <- setdiff(args, BANDS_ONLY)
only <- if (length(args)) args else NULL
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
MASTER_SEED <- 20260923L; N <- 20000L
BANDS <- list(c(40, 60), c(60, 75), c(75, 90), c(90, 110), c(110, 130), c(130, 150))
MODELS <- c(a = "base", b = "struct2020", c = "k2016_bmi", d = "k2016_bmi_vc0817"); if (!is.null(only)) MODELS <- MODELS[names(MODELS) %in% only]
if (length(BANDS_ONLY)) BANDS <- Filter(function(b) sprintf("%g-%g", b[1], b[2]) %in% BANDS_ONLY, BANDS)
SFX <- if (JITTER) "" else "_nojitter"
HEIGHT <- list(mean = 170, sd = 9, trunc = c(150, 195))
out_dir <- proj_path("results", "weight_generalization"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log(paste0("weight_generalization_", paste(names(MODELS), collapse = ""), SFX), master_seed = MASTER_SEED, run_mode = "final", extra = list(n_per_band = N, jitter = JITTER))
res <- rbindlist(lapply(names(MODELS), function(mk) rbindlist(lapply(BANDS, function(b) {
  rv <- resolve_variant(MODELS[[mk]], design, sc); p <- rv$p
  wt <- list(dist = "uniform", trunc = b, height = HEIGHT)
  pop <- run_individual_population(N, p, design, "B0", MASTER_SEED, wt, jitter = JITTER, model_id = p$model_id, tag = sprintf("wtgen_%s_%g_%g", mk, b[1], b[2]))
  x <- pop$nca
  cat(sprintf("  model %s band %g–%g %s\n", mk, b[1], b[2], format(Sys.time(), "%H:%M:%S")))
  data.table(model = mk, model_variant = MODELS[[mk]], band = sprintf("%g-%g", b[1], b[2]), n = nrow(x), BMI_median = median(x$BMI), WT_median = median(x$WT),
             extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(x$pct_extrap_true), extrap_true_max = max(x$pct_extrap_true, na.rm = TRUE),
             coverage_lt80_pct = 100 * mean(x$coverage_true < 0.80, na.rm = TRUE), extrap_nca_median = median(x$pct_extrap, na.rm = TRUE),
             extrap_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE), reliable_pct = 100 * mean(x$reliable), reliable_rsq_extrap_pct = 100 * mean(x$lambda_ok & x$adj_r2 >= 0.80 & x$pct_extrap <= 20, na.rm = TRUE), lambda_fail_pct = 100 * mean(!x$lambda_ok),
             tlast_median = median(x$tlast, na.rm = TRUE), AUClast_geo = geo_mean(x$AUClast), AUClast_logcv = log_cv_pct(x$AUClast),
             dev_range_note = if (b[1] >= 130) "개발 자료 범위 확인 불가, 외삽 가능성" else "")
}))))
fwrite(res[, jitter := JITTER], file.path(out_dir, sprintf("weight_bands_B0_%s%s.csv", paste(names(MODELS), collapse = ""), SFX)))
print(res[, .(model, band, BMI_median, extrap_true_median, extrap_true_p95, coverage_lt80_pct, reliable_pct, lambda_fail_pct, AUClast_geo)], digits = 3)
append_run_log(logfile, "done")
