#!/usr/bin/env Rscript
# ka 보정 (SPEC §3.4). 목표값(관측 Cmax·tmax, 300 mg)이 config/calibration_targets.yaml 에 있으면 최적화,
# 없으면 ka 격자 민감도만 산출한다. 결과: results/calibration/
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "dev"
p   <- load_params(run_mode)
ct  <- read_cfg("calibration_targets.yaml")
dz  <- read_cfg("design_clot2021.yaml")
out_dir <- proj_path("results", "calibration"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("calibrate_ka", master_seed = NA, run_mode = run_mode)

dose <- ct$dose_mg
WT   <- if (!is.null(ct$weight_kg_for_calibration)) ct$weight_kg_for_calibration else p$cov$WT_ref
lloq <- dz$lloq_mg_L
days <- as.numeric(dz$sample_days_post_dose)
dense <- seq(0, 84, by = 1 / 24)

sim_metrics <- function(ka_val) {
  ip <- typical_subject(p, WT = WT); ip[, ka := ka_val]   # data.table 스코프: 인자명을 열이름과 다르게
  tr <- true_auc(ip, dose, t_grid = sort(unique(c(dense, days))))
  prof <- tr$profile
  i <- which.max(prof$C)
  sched <- prof[time %in% days & time > 0]
  tl <- true_auc_to_tlast(sched, lloq)
  data.table(ka = ka_val, Cmax = prof$C[i], tmax = prof$time[i],
             tlast_sched = tl$tlast_true, AUClast_true = tl$AUClast_true, AUCinf_true = tr$inf$AUCinf_true,
             ratio_true_pct = 100 * tl$AUClast_true / tr$inf$AUCinf_true)
}

grid <- rbindlist(lapply(as.numeric(ct$ka_grid_1_day), sim_metrics))
fwrite(grid, file.path(out_dir, "ka_grid_sensitivity.csv"))
cat("\nka 격자 민감도 (대표 개체,", dose, "mg, WT =", WT, "kg):\n"); print(grid)

if (!is.null(ct$cmax_mg_L$value) && !is.null(ct$tmax_day$value)) {
  w <- if (!is.null(ct$tmax_weight)) ct$tmax_weight else 1
  loss <- function(ka) { m <- sim_metrics(ka); log(m$Cmax / ct$cmax_mg_L$value)^2 + w * log(m$tmax / ct$tmax_day$value)^2 }
  fit <- optimize(loss, interval = c(0.02, 5), tol = 1e-5)
  best <- sim_metrics(fit$minimum)
  res <- data.table(ka_hat = fit$minimum, loss = fit$objective, Cmax_sim = best$Cmax, Cmax_obs = ct$cmax_mg_L$value,
                    tmax_sim = best$tmax, tmax_obs = ct$tmax_day$value, WT = WT, dose = dose,
                    source_cmax = ct$cmax_mg_L$source, source_tmax = ct$tmax_day$source)
  fwrite(res, file.path(out_dir, "ka_fit.csv"))
  cat("\nka 보정 결과:\n"); print(res)
  cat("\n→ config/params_typical.yaml 의 theta.ka 를 다음으로 갱신하고 status를 confirmed로, source에 보정 근거를 적으세요:\n",
      sprintf("   ka: {value: %.4f, unit: 1/day, status: confirmed, source: \"scripts/01_calibrate_ka.R 보정; 목표 Cmax=%s (%s), tmax=%s (%s)\"}\n",
              fit$minimum, ct$cmax_mg_L$value, ct$cmax_mg_L$source, ct$tmax_day$value, ct$tmax_day$source))
  append_run_log(logfile, "ka_hat=", fit$minimum)
} else {
  cat("\n[PENDING] 보정 목표(Cmax, tmax)가 config/calibration_targets.yaml 에 없습니다. 격자 민감도만 산출했습니다 (SPEC Q3).\n")
  append_run_log(logfile, "calibration targets pending; grid only")
}
