#!/usr/bin/env Rscript
# 단계 1 정량 gate 미달 진단 (파라미터 조정 없음, SPEC §1.3). 결과: results/step1/diag_*.csv
# 1) 용량 정규화 노출: 관측 vs 모의 (200→300→600 mg 비선형성 비교)
# 2) 확정 전 가정(체중 SD·절단·LLOQ·Cohen 채혈)에 대한 실패 데이터셋의 민감도
# 3) 사전 명시된 대안 모델(Kovalenko 2020 Model 1, ω²×1.5)로 gate 재평가
# 4) 200 mg·80 kg 대표 개체 프로파일: Cmax는 맞고 AUC가 모자라는 구간
source("R/00_setup.R"); source_project()
dz <- read_cfg("design_clot2021.yaml"); design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", "step1"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("step1_diagnostics", master_seed = "derive_seed(tag)", run_mode = "final")
sr <- dz$cohort_sim$sex_ratio_male$value; tol <- dz$gate$auclast_mean_tol_pct; cvr <- as.numeric(dz$gate$log_cv_range_pct)
NSIM <- 10000
sched_of <- function(ds) { s <- dz$dataset_schedules[[ds$schedule]]; as.numeric(if (is.list(s)) s$days else s) }
run_gate <- function(p, tag, datasets = dz$datasets, nsim = NSIM, override = list()) {
  rbindlist(lapply(datasets, function(ds) {
    b <- ds_weight_bounds(ds); wsd <- ds_weight_sd(ds); sched <- sched_of(ds); pp <- p
    if (!is.null(override$weight_sd)) wsd <- override$weight_sd
    if (!is.null(override$bounds)) b <- list(male = override$bounds, female = override$bounds)
    if (!is.null(override$lloq)) pp$lloq <- override$lloq
    if (!is.null(override$sched)) sched <- override$sched
    sim <- simulate_dataset(pp, ds$dose_mg, ds$weight_mean, wsd, sched, nsim, paste0("diag_", tag, "_", ds$id), sr, b$male, b$female)
    gate_row(ds, sim$nca, tol, cvr)[, variant := tag]
  }))
}

# 1) 기본 모델 gate + 용량 정규화
p <- load_params("k2016")
cache <- function(f, expr) { fp <- file.path(out_dir, f); if (file.exists(fp) && Sys.getenv("DIAG_RECOMPUTE") != "1") fread(fp) else { x <- expr; fwrite(x, fp); x } }
g0 <- cache("diag_gate_base_raw.csv", run_gate(p, "base"))
g0[, `:=`(obs_per_mg = AUClast_obs_mean / dose_mg, sim_per_mg = AUClast_sim_mean / dose_mg)]
fwrite(g0, file.path(out_dir, "diag_gate_base.csv"))
nonasian <- g0[id %in% c("pkm14271_200_nonasian", "clot300_nonasian_pooled")]
dn <- data.table(comparison = c("비아시아인 ≈80 kg: 200 mg vs 300 mg (AUClast/mg 비)", "중국인 ≈59 kg: 600 mg vs 300 mg (AUClast/mg 비)", "Cohen 2022 두 arm 평균 vs 모의 (200 mg, 80.6 kg)"),
                 observed = c(g0[id == "pkm14271_200_nonasian", obs_per_mg] / g0[id == "clot300_nonasian_pooled", obs_per_mg],
                              g0[id == "clot600_chinese", obs_per_mg] / g0[id == "clot300_chinese", obs_per_mg],
                              mean(g0[grepl("cohen", id), AUClast_obs_mean])),
                 simulated = c(g0[id == "pkm14271_200_nonasian", sim_per_mg] / g0[id == "clot300_nonasian_pooled", sim_per_mg],
                               g0[id == "clot600_chinese", sim_per_mg] / g0[id == "clot300_chinese", sim_per_mg],
                               mean(g0[grepl("cohen", id), AUClast_sim_mean])))
dn[, sim_over_obs := simulated / observed]
fwrite(dn, file.path(out_dir, "diag_dose_normalized.csv")); cat("\n1) 용량 정규화 비교:\n"); print(dn)

# 2) 실패 데이터셋의 확정 전 가정 민감도
fail_ids <- c("pkm14271_200_nonasian", "cohen2022_200_armA", "cohen2022_200_armB")
fds <- Filter(function(d) d$id %in% fail_ids, dz$datasets)
sens <- cache("diag_assumption_sensitivity.csv", rbindlist(list(
  rbindlist(lapply(c(6, 8, 10, 12, 14), function(s) run_gate(p, paste0("wtSD_", s), fds, override = list(weight_sd = s)))),
  run_gate(p, "bounds_40_150", fds, override = list(bounds = c(40, 150))),
  rbindlist(lapply(c(0.05, 0.02), function(l) run_gate(p, paste0("lloq_", l), fds, override = list(lloq = l)))),
  run_gate(p, "sched_B0_to42", fds, override = list(sched = c(0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42))),
  run_gate(p, "sched_clot_to56", fds, override = list(sched = as.numeric(dz$dataset_schedules$clot))))))
cat("\n2) 실패 데이터셋의 확정 전 가정 민감도 (AUClast 모의/관측):\n"); print(dcast(sens[, .(variant, id, r = round(AUClast_ratio, 3))], variant ~ id, value.var = "r"))

# 3) 사전 명시 대안 모델
sc <- load_scenarios()
alt <- cache("diag_gate_alt_raw.csv", rbindlist(list(
  run_gate(resolve_variant("struct2020", design, sc)$p, "struct2020"),
  run_gate(resolve_variant("iiv150", design, sc)$p, "iiv150"))))
alt <- rbind(g0[, names(alt), with = FALSE], alt)
fwrite(alt, file.path(out_dir, "diag_gate_alt_models.csv"))
cat("\n3) 대안 모델 gate (AUClast 모의/관측):\n"); print(dcast(alt[, .(id, variant, r = round(AUClast_ratio, 3))], id ~ variant, value.var = "r"))
cat("   log-CV:\n"); print(dcast(alt[, .(id, variant, cv = round(AUClast_sim_logcv, 1))], id ~ variant, value.var = "cv"))
cat("   통과 여부(평균·CV):\n"); print(alt[, .(all_mean = all(pass_mean), all_cv = all(pass_cv), n_fail_mean = sum(!pass_mean)), by = variant])

# 4) 대표 개체 프로파일: 200 mg 80 kg vs 300 mg 78 kg
prof <- rbindlist(lapply(list(c(200, 80), c(300, 78), c(300, 60), c(600, 58.3)), function(x) {
  ip <- typical_subject(p, WT = x[2]); tr <- true_auc(ip, x[1], t_grid = seq(0, 70, by = 0.25))
  pr <- tr$profile; cross <- max(pr[C >= p$lloq, time])
  data.table(dose = x[1], WT = x[2], Cmax = max(pr$C), tmax = pr[which.max(C), time], t_cross_lloq = cross,
             AUC_0_14 = pr[time == 14, auc], AUC_14_inf = tr$inf$AUCinf_true - pr[time == 14, auc], AUCinf_true = tr$inf$AUCinf_true,
             frac_after_d14 = 1 - pr[time == 14, auc] / tr$inf$AUCinf_true)
}))
fwrite(prof, file.path(out_dir, "diag_typical_profiles.csv")); cat("\n4) 대표 개체 프로파일:\n"); print(prof)
# 5) 용량 형태(200/300 mg, 600/300 mg AUClast/mg 비)를 어떤 가정이 결정하는가 — 채택 아님, 원인 가정 식별용(SPEC §1.3)
#    변형 범위는 검토자가 사전에 그럴듯하다고 명시한 범위(지시서 §6 Km 0.005–0.1, Vmax ×0.8–×1.5)와 구조 대안(MM 소실의 체중 비례 제거)
shape_ds <- Filter(function(d) d$id %in% c("pkm14271_200_nonasian", "clot300_nonasian_pooled", "clot300_chinese", "clot600_chinese"), dz$datasets)
shape_variant <- function(tag, mult = list(), mm_not_weight_scaled = FALSE, model = "k2016") {
  pp <- apply_multipliers(load_params(model), mult)
  rbindlist(lapply(shape_ds, function(ds) {
    b <- ds_weight_bounds(ds)
    wt_spec <- list(mean = ds$weight_mean, sd = ds_weight_sd(ds), male_bounds = b$male, female_bounds = b$female)
    subj <- with_seed(derive_seed("shape", tag, ds$id, "subj"), make_subjects(4000, pp, wt_spec, sr, 0))
    planned <- sort(unique(c(0, sched_of(ds)))); obs <- CJ(id = subj$id, planned = planned)[, time := planned]
    eps <- with_seed(derive_seed("shape", tag, ds$id, "eps"), draw_eps(subj$id, planned, pp$sigma))
    ip <- individual_params(pp, subj)
    if (mm_not_weight_scaled) ip[, Vmax := Vmax * (pp$cov$WT_ref / WT)^pp$cov$theta_WT]   # MM 양 기준 속도 = Vmax·Vc_i 에서 체중 비례 성분 제거
    sim <- simulate_observations(ip, obs, ds$dose_mg, pp$lloq, eps, model_id = pp$model_id)
    nca <- run_nca(sim$obs)
    data.table(variant = tag, id = ds$id, dose = ds$dose_mg, sim_mean = mean(nca$AUClast, na.rm = TRUE), obs_mean = ds$auclast_mean)
  }))
}
shp <- cache("diag_dose_shape.csv", rbindlist(list(
  shape_variant("base"),
  shape_variant("struct2020", model = "k2020"),
  shape_variant("Km_x0.5", list(Km = 0.5)), shape_variant("Km_x2", list(Km = 2)), shape_variant("Km_x5", list(Km = 5)), shape_variant("Km_x10", list(Km = 10)),
  shape_variant("Vmax_x0.8", list(Vmax = 0.8)), shape_variant("Vmax_x1.25", list(Vmax = 1.25)),
  shape_variant("MM_not_weight_scaled", mm_not_weight_scaled = TRUE))))
w <- dcast(shp, variant ~ id, value.var = "sim_mean")
w[, `:=`(ratio_200v300_nonasian = (pkm14271_200_nonasian / 200) / (clot300_nonasian_pooled / 300),
         ratio_600v300_chinese = (clot600_chinese / 600) / (clot300_chinese / 300),
         pkm200_sim_obs = pkm14271_200_nonasian / 339, nonasian300_sim_obs = clot300_nonasian_pooled / 544,
         chinese300_sim_obs = clot300_chinese / 792, chinese600_sim_obs = clot600_chinese / 2110)]
obs_row <- data.table(variant = "OBSERVED", ratio_200v300_nonasian = (339 / 200) / (544 / 300), ratio_600v300_chinese = (2110 / 600) / (792 / 300))
shape_tab <- rbind(obs_row, w[, .(variant, ratio_200v300_nonasian, ratio_600v300_chinese, pkm200_sim_obs, nonasian300_sim_obs, chinese300_sim_obs, chinese600_sim_obs)], fill = TRUE)
fwrite(shape_tab, file.path(out_dir, "diag_dose_shape_summary.csv"))
cat("\n5) 용량 형태 진단 (AUClast/mg 비; 채택 아님):\n"); print(shape_tab, digits = 3)
append_run_log(logfile, "done")
