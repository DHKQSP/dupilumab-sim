#!/usr/bin/env Rscript
# rationale 메모용 요약표 (검토 의견 3차 §5, 통합본 §0-1·§4): B0, 60–90 kg. 모든 표에 반복 수·대상자 수·Monte Carlo 95% 구간.
# 주 모델(2016)과 2020 Model 1(자체 IIV·잔차)을 나란히. Pillar 2 주 모델은 5,000회 이상 결과가 있는 시나리오는 그것을, 나머지는 500회.
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "rationale"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
design <- read_cfg("trial_design.yaml"); split <- design$stratification$split_kg$value
models <- c(k2016 = "base", k2020 = "struct2020")
wprop <- function(v) { v <- v[!is.na(v)]; w <- wilson_ci(sum(v), length(v)); sprintf("%.2f [%.2f, %.2f]", w$est, w$lo, w$hi) }

# Pillar 1 — 개체 수준(20,000명, B0)
p1 <- rbindlist(lapply(names(models), function(m) {
  nca <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", models[[m]])))[schedule == "B0"]
  st <- strata_from_weight_spec(weight_spec_from_design(design, "base"), split); nca[, stratum := weight_stratum(WT, st$breaks, st$labels)]
  one <- function(x, grp) data.table(model = m, group = grp, n_subjects = nrow(x),
    extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(x$pct_extrap_true), extrap_true_max = max(x$pct_extrap_true, na.rm = TRUE),
    coverage_lt80_pct_ci = wprop(x$coverage_true < 0.80), extrap_nca_median = median(x$pct_extrap, na.rm = TRUE), extrap_nca_p95 = q95(x$pct_extrap),
    extrap_nca_gt20_pct_ci = wprop(x$pct_extrap > 20), lambda_ok_pct_ci = wprop(x$lambda_ok), reliable_pct_ci = wprop(x$reliable),
    reliable_pct = 100 * mean(x$reliable), extrap_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE))
  rbind(one(nca, "전체"), rbindlist(lapply(levels(nca$stratum), function(s) one(nca[stratum == s], paste0("체중 ", s)))))
}))
fwrite(p1, file.path(out_dir, "pillar1_coverage_B0.csv"))
rng <- p1[group == "전체", .(reliable_min = min(reliable_pct), reliable_max = max(reliable_pct), gt20_min = min(extrap_gt20_pct), gt20_max = max(extrap_gt20_pct))]
fwrite(rng, file.path(out_dir, "pillar1_two_model_range.csv"))

# Pillar 2 — 제품 차이 시나리오(B0)
p2_from <- function(raw, m, src) { st <- summarize_trials_ci(read_raw(raw)); pe <- st$per_endpoint; cc <- st$concordance
  x <- dcast(pe[endpoint %in% c("Cmax", "AUClast", "AUCinf_reliable", "AUCinf_true")], scenario + n_trials ~ endpoint, value.var = c("GMR_mean", "pass_rate", "pass_lo", "pass_hi"))
  merge(x, cc[, .(scenario, agree, agree_lo, agree_hi, last_pass_inf_fail, last_pass_inf_fail_lo, last_pass_inf_fail_hi, cor_logGMR_last_infrel, joint_last_cmax, joint_last_cmax_lo, joint_last_cmax_hi)], by = "scenario")[, `:=`(model = m, source = src)][] }
f5000 <- proj_path("results", "trials5000", "products5000_be_raw_base.csv")
p2_5000 <- p2_from(f5000, "k2016", "5,000회 이상")
p2_500 <- p2_from(proj_path("results", "trials", "products_be_raw.csv"), "k2016", "500회")[!scenario %in% p2_5000$scenario]
p2_20 <- p2_from(proj_path("results", "trials", "products_be_raw_struct2020.csv"), "k2020", "500회")
p2 <- rbind(p2_5000, p2_500, p2_20, fill = TRUE); setcolorder(p2, c("model", "scenario", "source", "n_trials"))
fwrite(p2, file.path(out_dir, "pillar2_products_B0.csv"))

# Pillar 3 — Km ×0.5–×10 vs 동일 제품: 같은 시험의 쌍대 변화(공통 난수)
p3 <- rbindlist(lapply(names(models), function(m) {
  raw <- read_raw(proj_path("results", "trials", if (m == "k2016") "products_be_raw.csv" else "products_be_raw_struct2020.csv"))[method == "pooled_t"]
  ref <- raw[scenario == "S00", .(trial, endpoint, GMR0 = GMR, pass0 = pass)]
  x <- merge(raw[scenario %in% c("KM05", "KM2", "KM5", "KM10"), .(trial, scenario, endpoint, GMR, pass)], ref, by = c("trial", "endpoint"))
  x[, {
    r <- mc_mean_ci(100 * (GMR / GMR0 - 1)); d <- paired_prop_diff_ci(pass, pass0)
    .(model = m, n_trials = .N, GMR_change_pct = r$est, GMR_change_lo = r$lo, GMR_change_hi = r$hi, pass_change_pp = d$est, pass_change_lo = d$lo, pass_change_hi = d$hi)
  }, by = .(scenario, endpoint)]
}))
p3sum <- p3[, .(max_abs_GMR_change_pct = max(abs(GMR_change_pct)), max_abs_pass_change_pp = max(abs(pass_change_pp)), worst = .SD[which.max(abs(GMR_change_pct)), paste(scenario, endpoint)]), by = model]
fwrite(p3, file.path(out_dir, "pillar3_km_B0.csv")); fwrite(p3sum, file.path(out_dir, "pillar3_km_summary.csv"))
print(p1); print(rng); print(p2[, .(model, scenario, source, n_trials, GMR_mean_AUClast, GMR_mean_AUCinf_true, pass_rate_AUClast, pass_rate_AUCinf_reliable, agree, last_pass_inf_fail)]); print(p3sum)
