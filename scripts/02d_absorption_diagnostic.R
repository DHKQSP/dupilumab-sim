#!/usr/bin/env Rscript
# 부록: 제형 흡수 진단 (검토 의견 3차 §1.3, §2). 채택 아님 — 어떤 산출물에도 이 값을 쓰지 않는다.
# 200 mg 1.14 mL(175 mg/mL) 조건에서만 ka에 배율(1, 1.5, 2, 3)을 적용. 300 mg 파라미터는 변경하지 않는다(참조 행은 ka ×1).
# 또한 관측 용량 형태 비를 200 mg 4개 arm의 n 가중 평균으로 재계산하고 PKM14271 Test 단독 값과 병기한다.
source("R/00_setup.R"); source_project()
dz <- read_cfg("design_clot2021.yaml"); ad <- dz$absorption_diagnostic
p0 <- load_params("k2016")
out_dir <- proj_path("results", "step1"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("absorption_diagnostic", master_seed = "derive_seed(tag)", run_mode = "diagnostic")
sr <- dz$cohort_sim$sex_ratio_male$value; NSIM <- 20000
sched_of <- function(key) { s <- dz$dataset_schedules[[key]]; as.numeric(unlist(if (is.list(s) && !is.null(s$days)) s$days else s)) }
ds_by_id <- function(i) Filter(function(d) d$id == i, dz$datasets)[[1]]

run_cond <- function(cond, kam) {
  p <- apply_multipliers(p0, list(ka = kam))                     # 이 진단 안에서만 쓰는 사본
  b <- if (!is.null(cond$weight_bounds)) as.numeric(unlist(cond$weight_bounds)) else NULL
  bmale <- if (is.null(b)) c(50, 100) else b; bfem <- if (is.null(b)) c(40, 90) else b
  sim <- simulate_dataset(p, cond$dose_mg, cond$weight_mean, cond$weight_sd, sched_of(cond$schedule), NSIM, paste0("absdiag_", cond$id), sr, bmale, bfem)
  obs <- rbindlist(lapply(unlist(cond$obs_ids), function(i) { d <- ds_by_id(i); data.table(obs_id = i, n = d$n_subj, auc = d$auclast_mean, cmax = d$cmax_mean, tmax = if (is.null(d$tmax_median)) NA_real_ else d$tmax_median) }))
  auc_w <- sum(obs$n * obs$auc) / sum(obs$n); cmax_w <- sum(obs$n * obs$cmax) / sum(obs$n)
  data.table(condition = cond$id, dose_mg = cond$dose_mg, ka_multiplier = kam, ka = p$theta[["ka"]],
             Cmax_sim_mean = mean(sim$nca$Cmax, na.rm = TRUE), tmax_sim_median = median(sim$nca$tmax, na.rm = TRUE),
             AUClast_sim_mean = mean(sim$nca$AUClast, na.rm = TRUE),
             obs_arms = paste(sprintf("%s %g (Cmax %g, tmax %s)", obs$obs_id, obs$auc, obs$cmax, ifelse(is.na(obs$tmax), "NA", format(obs$tmax))), collapse = "; "),
             obs_tmax_median = paste(unique(obs$tmax), collapse = "/"),
             AUClast_obs_nweighted = auc_w, AUClast_ratio_nweighted = mean(sim$nca$AUClast, na.rm = TRUE) / auc_w,
             Cmax_ratio_nweighted = mean(sim$nca$Cmax, na.rm = TRUE) / cmax_w,
             AUClast_ratio_min_arm = mean(sim$nca$AUClast, na.rm = TRUE) / max(obs$auc), AUClast_ratio_max_arm = mean(sim$nca$AUClast, na.rm = TRUE) / min(obs$auc),
             n_obs_total = sum(obs$n))
}
kams <- as.numeric(unlist(ad$ka_multipliers))
tab <- rbindlist(c(lapply(ad$conditions, function(cond) rbindlist(lapply(kams, function(k) run_cond(cond, k)))),
                   list(run_cond(ad$reference_300mg, ad$reference_300mg$ka_multiplier))))
fwrite(tab, file.path(out_dir, "appendix_absorption_diagnostic.csv"))
cat("\n부록: 제형 흡수 진단 (채택 아님)\n"); print(tab[, .(condition, ka_multiplier, Cmax_sim_mean = round(Cmax_sim_mean, 1), tmax_sim_median, AUClast_sim_mean = round(AUClast_sim_mean), AUClast_obs_nweighted = round(AUClast_obs_nweighted), AUClast_ratio_nweighted = round(AUClast_ratio_nweighted, 3), Cmax_ratio_nweighted = round(Cmax_ratio_nweighted, 3))])

# 관측 용량 형태 비: 200 mg 4개 arm n 가중 vs PKM14271 Test 단독 (300 mg pooled 544 기준)
ids200 <- c("PKM14271_200mg_test", "PKM14271_200mg_reference", "Cohen2022_200mg_AI", "Cohen2022_200mg_PFSS")
a200 <- rbindlist(lapply(ids200, function(i) { d <- ds_by_id(i); data.table(id = i, n = d$n_subj, mean = d$auclast_mean, sd = d$auclast_sd) }))
m300 <- ds_by_id("clot300_nonasian_pooled")
w200 <- sum(a200$n * a200$mean) / sum(a200$n)
shape_obs <- data.table(basis = c("200 mg 4개 arm n 가중", "PKM14271 Test 단독"), AUClast_200 = c(w200, a200[id == "PKM14271_200mg_test", mean]),
                        n_200 = c(sum(a200$n), 19), AUClast_300 = m300$auclast_mean)
shape_obs[, ratio_per_mg := (AUClast_200 / 200) / (AUClast_300 / 300)]
# 모델: 200 mg 조건(PKM 유사 n=38, Cohen 유사 n=125) n 가중 모의 평균 / 300 mg 약 78 kg 모의(ka ×1)
m300_sim <- tab[condition == ad$reference_300mg$id, AUClast_sim_mean]
shape_model <- tab[dose_mg == 200, .(AUClast_200_sim = sum(AUClast_sim_mean * n_obs_total) / sum(n_obs_total)), by = ka_multiplier]
shape_model[, `:=`(AUClast_300_sim = m300_sim, ratio_per_mg = (AUClast_200_sim / 200) / (m300_sim / 300))]
# PKM14271 유사 조건 단독 기준(검토자 표기 0.73, 0.78–0.81과 같은 기준)
shape_model <- merge(shape_model, tab[condition == "pkm14271_like", .(ka_multiplier, AUClast_200_sim_pkm_only = AUClast_sim_mean)], by = "ka_multiplier")
shape_model[, ratio_per_mg_pkm_only := (AUClast_200_sim_pkm_only / 200) / (m300_sim / 300)]
fwrite(shape_obs, file.path(out_dir, "appendix_dose_shape_observed.csv")); fwrite(shape_model, file.path(out_dir, "appendix_dose_shape_model.csv"))
cat("\n관측 용량 형태 비(AUClast/mg, 200 vs 300 mg):\n"); print(shape_obs)
cat("\n모델 용량 형태 비(200 mg 조건만 ka 배율, 채택 아님):\n"); print(shape_model)
append_run_log(logfile, "done")
