#!/usr/bin/env Rscript
# rationale 메모용 요약표 (검토 의견 3차 §5): B0 일정, 60–90 kg. 주 모델(2016)과 구조 민감도(2020)를 나란히.
# Pillar 1 커버리지, Pillar 2 판정 일치(제품 차이 14개), Pillar 3 결합력(Km) 비가시성. 결과: results/rationale/
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "rationale"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
design <- read_cfg("trial_design.yaml"); split <- design$stratification$split_kg$value
models <- c(k2016 = "base", k2020 = "struct2020")

# Pillar 1 — 개인 수준(20,000명, B0)
p1 <- rbindlist(lapply(names(models), function(m) {
  nca <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", models[[m]])))[schedule == "B0"]
  st <- strata_from_weight_spec(weight_spec_from_design(design, "base"), split)
  nca[, stratum := weight_stratum(WT, st$breaks, st$labels)]
  one <- function(x, grp) data.table(model = m, group = grp, n = nrow(x),
    extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(x$pct_extrap_true), coverage_lt80_pct = 100 * mean(x$coverage_true < 0.80, na.rm = TRUE),
    coverage_true_median = 100 * median(x$coverage_true, na.rm = TRUE), coverage_true_p05 = 100 * q05(x$coverage_true),
    extrap_nca_median = median(x$pct_extrap, na.rm = TRUE), extrap_nca_p95 = q95(x$pct_extrap), extrap_nca_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE),
    lambda_ok_pct = 100 * mean(x$lambda_ok), reliable_pct = 100 * mean(x$reliable),
    err_AUClast_vs_inf_true_median_pct = 100 * (exp(median(x$err_AUClast_vs_true_inf, na.rm = TRUE)) - 1))
  rbind(one(nca, "전체"), rbindlist(lapply(levels(nca$stratum), function(s) one(nca[stratum == s], paste0("체중 ", s)))))
}))
fwrite(p1, file.path(out_dir, "pillar1_coverage_B0.csv")); cat("\nPillar 1 (B0, 60–90 kg):\n"); print(p1, digits = 3)

# Pillar 2 — 제품 차이 14개(B0, 500회, 117명/arm)
p2 <- rbindlist(lapply(names(models), function(m) {
  sfx <- if (models[[m]] == "base") "" else "_struct2020"
  pe <- fread(proj_path("results", "trials", paste0("products_per_endpoint", sfx, ".csv")))
  co <- fread(proj_path("results", "trials", paste0("products_concordance", sfx, ".csv")))
  w <- dcast(pe[endpoint %in% c("Cmax", "AUClast", "AUCinf_reliable", "AUCinf_true")], scenario ~ endpoint, value.var = c("GMR_mean", "pass_rate"))
  w <- merge(w, co[, .(scenario, agree_last_infrel, last_pass_infrel_fail, last_fail_infrel_pass, cor_logGMR_last_infrel, cor_logGMR_last_true, pass_both_primary)], by = "scenario")
  w[, model := m][]
}))
setcolorder(p2, c("model", "scenario"))
fwrite(p2, file.path(out_dir, "pillar2_products_B0.csv")); cat("\nPillar 2:\n"); print(p2[, .(model, scenario, GMR_AUClast = round(GMR_mean_AUClast, 3), GMR_AUCinf_rel = round(GMR_mean_AUCinf_reliable, 3), GMR_true = round(GMR_mean_AUCinf_true, 3), pass_AUClast = pass_rate_AUClast, pass_AUCinf_rel = pass_rate_AUCinf_reliable, agree = agree_last_infrel, last_pass_inf_fail = last_pass_infrel_fail, cor = round(cor_logGMR_last_infrel, 3))])

# Pillar 3 — Km 0.5–10배 vs 동일 제품: 전 지표 변화
p3 <- rbindlist(lapply(names(models), function(m) {
  sfx <- if (models[[m]] == "base") "" else "_struct2020"
  pe <- fread(proj_path("results", "trials", paste0("products_per_endpoint", sfx, ".csv")))
  ref <- pe[scenario == "S00", .(endpoint, GMR_ref = GMR_mean, pass_ref = pass_rate, width_ref = width_mean_pp)]
  x <- merge(pe[scenario %in% c("KM05", "KM2", "KM5", "KM10"), .(scenario, endpoint, GMR_mean, pass_rate, width_mean_pp)], ref, by = "endpoint")
  x[, `:=`(model = m, GMR_change_pct = 100 * (GMR_mean / GMR_ref - 1), pass_change_pp = pass_rate - pass_ref, width_change_pp = width_mean_pp - width_ref)][]
}))
p3sum <- p3[, .(max_abs_GMR_change_pct = max(abs(GMR_change_pct)), max_abs_pass_change_pp = max(abs(pass_change_pp)), worst = .SD[which.max(abs(GMR_change_pct)), paste(scenario, endpoint)]), by = model]
fwrite(p3, file.path(out_dir, "pillar3_km_B0.csv")); fwrite(p3sum, file.path(out_dir, "pillar3_km_summary.csv"))
cat("\nPillar 3 (Km ×0.5–×10 vs 동일 제품):\n"); print(p3[, .(model, scenario, endpoint, GMR_change_pct = round(GMR_change_pct, 2), pass_change_pp = round(pass_change_pp, 1))]); print(p3sum)
