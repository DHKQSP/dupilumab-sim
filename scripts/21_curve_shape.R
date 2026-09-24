#!/usr/bin/env Rscript
# §5-1 곡선 모양·곡률 민감도(Pillar 1): 개체 수준 20,000명, B0, 양 군 공통 배율. 주 모델.
# 변형: base, Km ×0.5·2·5·10, Vmax ×0.8·×1.25, 스트레스 테스트 Vmax ×0.5. 300 mg 약 78 kg 관측 pooled 544 대비 AUClast 비도 산출(gate 표기용).
source("R/00_setup.R"); source_project()
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); dz <- read_cfg("design_clot2021.yaml")
MASTER_SEED <- 20260923L; N <- 20000L
V <- c("base", "km05_both", "km2_both", "km5_both", "km10_both", "vmax080_both", "vmax125_both", "vmax050_both")
out_dir <- proj_path("results", "curve_shape"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("curve_shape", master_seed = MASTER_SEED, run_mode = "final", extra = list(variants = V, n = N))
pooled <- Filter(function(d) d$id == "clot300_nonasian_pooled", dz$datasets)[[1]]
res <- rbindlist(lapply(V, function(v) {
  rv <- resolve_variant(v, design, sc); p <- rv$p
  pop <- run_individual_population(N, p, design, "B0", MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = paste0("curve_", v))
  x <- pop$nca
  g <- simulate_ds(p, pooled, dz, 10000, paste0("curve_gate_", v), dz$cohort_sim$sex_ratio_male$value)
  cat(sprintf("  %s done %s\n", v, format(Sys.time(), "%H:%M:%S")))
  data.table(variant = v, label = sc$variants[[v]]$label, n = nrow(x), stress_test = !is.null(sc$variants[[v]]$stress_test),
             extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(x$pct_extrap_true), extrap_true_max = max(x$pct_extrap_true, na.rm = TRUE),
             coverage_lt80_pct = 100 * mean(x$coverage_true < 0.80, na.rm = TRUE),
             extrap_nca_median = median(x$pct_extrap, na.rm = TRUE), extrap_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE),
             reliable_pct = 100 * mean(x$reliable), lambda_fail_pct = 100 * mean(!x$lambda_ok), tlast_median = median(x$tlast, na.rm = TRUE),
             AUClast_geo = geo_mean(x$AUClast), AUClast_logcv = log_cv_pct(x$AUClast),
             AUClast_mean_78kg_pooled = mean(g$AUClast, na.rm = TRUE), AUClast_ratio_vs_obs544 = mean(g$AUClast, na.rm = TRUE) / pooled$auclast_mean)
}))
res[, gate_note := fifelse(abs(AUClast_ratio_vs_obs544 - 1) > 0.15, "gate 불통과(±15% 밖)", "gate 범위 안")]
res[stress_test == TRUE, gate_note := paste0(gate_note, ", 스트레스 테스트")]
# 검토자 독립 구현 기준값(8,000명): 실제 외삽 중앙값/95백분위
ref <- data.table(variant = c("km05_both", "km2_both", "km5_both", "km10_both", "vmax080_both", "vmax050_both"),
                  ref_note = c("Km 0.5–10배: 중앙값 0.38–0.67%, 95백분위 2.6–3.5%, 80% 미만 0", "", "", "", "0.53% / 2.63%, 80% 미만 0", "0.62% / 6.52%, 최대 31.7%, 80% 미만 0.11%"))
res <- merge(res, ref, by = "variant", all.x = TRUE)
fwrite(res, file.path(out_dir, "curve_shape_B0.csv")); print(res[, .(variant, extrap_true_median, extrap_true_p95, extrap_true_max, coverage_lt80_pct, extrap_nca_median, extrap_gt20_pct, reliable_pct, lambda_fail_pct, tlast_median, AUClast_geo, AUClast_ratio_vs_obs544, gate_note)], digits = 3)
append_run_log(logfile, "done")
