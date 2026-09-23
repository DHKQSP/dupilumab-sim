#!/usr/bin/env Rscript
# 정량 gate의 관측 평균 표본오차 맥락 (진단, gate 규칙 변경 아님). 관측 평균의 SE = SD/sqrt(n).
# z = (모의 평균 − 관측 평균)/SE_obs. 관측 용량 형태 비(200/300 mg)의 근사 95% 구간(델타법).
source("R/00_setup.R"); source_project()
dz <- read_cfg("design_clot2021.yaml")
out <- rbindlist(lapply(c("step1", "step1_k2020"), function(d) {
  g <- fread(proj_path("results", d, "step1b_quant_gate.csv"))
  n_obs <- vapply(dz$datasets, function(x) if (is.null(x$n_subj)) NA_real_ else as.numeric(x$n_subj), numeric(1))
  sd_obs <- vapply(dz$datasets, function(x) as.numeric(x$auclast_sd), numeric(1))
  g <- merge(g, data.table(id = vapply(dz$datasets, `[[`, "", "id"), n_obs2 = n_obs, sd_obs = sd_obs), by = "id")
  g[, `:=`(se_obs_pct = 100 * sd_obs / sqrt(n_obs2) / AUClast_obs_mean, z = (AUClast_sim_mean - AUClast_obs_mean) / (sd_obs / sqrt(n_obs2)), model = d)]
  g[, .(model, id, dose_mg, n_obs = n_obs2, AUClast_ratio, se_obs_pct, z, pass_mean)]
}))
fwrite(out, proj_path("results", "step1", "diag_gate_sampling_error.csv"))
print(out, digits = 3)
# 관측 200/300 mg 형태 비의 근사 95% 구간
r <- (339 / 200) / (544 / 300); se_log <- sqrt((128 / sqrt(19) / 339)^2 + (218 / sqrt(40) / 544)^2)
cat(sprintf("\n관측 200/300 mg AUClast/mg 비 %.3f, 근사 95%% 구간 [%.3f, %.3f] (두 연구 평균의 표본오차만 반영)\n", r, r * exp(-1.96 * se_log), r * exp(1.96 * se_log)))
fwrite(data.table(ratio = r, lo = r * exp(-1.96 * se_log), hi = r * exp(1.96 * se_log)), proj_path("results", "step1", "diag_observed_shape_ci.csv"))
