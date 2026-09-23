#!/usr/bin/env Rscript
# 단계 1 검증 요약 산출 (테스트와 별도로 표를 results/dev/ 에 남긴다). BEmaster·FDA 표 미입수 항목은 상태만 기록.
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "dev"
p  <- load_params(run_mode)
dz <- read_cfg("design_clot2021.yaml"); ct <- read_cfg("calibration_targets.yaml")
days <- as.numeric(dz$sample_days_post_dose); lloq <- dz$lloq_mg_L
out_dir <- proj_path("results", "dev"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("validate_step1", master_seed = NA, run_mode = run_mode)

targets <- rbindlist(lapply(dz$targets, function(x) data.table(dose = x$dose_mg, target_ratio_pct = x$auc_ratio_pct, target_tlast = x$tlast_median_day)))
grid <- c(p$theta[["ka"]], as.numeric(ct$ka_grid_1_day))
tab <- rbindlist(lapply(targets$dose, function(d) rbindlist(lapply(grid, function(ka) {
  ip <- typical_subject(p); ip[, ka := ka]
  tr <- true_auc(ip, d, t_grid = c(days, seq(0, 84, by = 0.25)))
  prof <- tr$profile
  sched <- prof[time %in% days & time > 0]
  tl <- true_auc_to_tlast(sched, lloq)
  data.table(dose = d, ka = ka, is_config_ka = ka == p$theta[["ka"]],
             Cmax_true = max(prof$C), tmax_true = prof$time[which.max(prof$C)],
             tlast_typical = tl$tlast_true, Clast_true = sched[time == tl$tlast_true, C],
             C_next_sched = { nx <- sched[time > tl$tlast_true][1]; if (nrow(nx)) nx$C else NA_real_ },
             AUClast_true = tl$AUClast_true, AUCinf_true = tr$inf$AUCinf_true,
             ratio_true_pct = 100 * tl$AUClast_true / tr$inf$AUCinf_true)
}))))
tab <- merge(tab, targets, by = "dose")
tab[, `:=`(tlast_match = tlast_typical == target_tlast, ratio_diff_pp = ratio_true_pct - target_ratio_pct)]
setorder(tab, dose, ka)
fwrite(tab, file.path(out_dir, "step1_typical_grid.csv"))
cat("\n단계 1 (a)(b): 대표 개체, ka 격자별 tlast·참값 비율\n"); print(tab[, .(dose, ka, Cmax_true = round(Cmax_true, 2), tmax_true, tlast_typical, target_tlast, tlast_match, Clast_true = signif(Clast_true, 3), C_next_sched = signif(C_next_sched, 3), ratio_true_pct = round(ratio_true_pct, 2), target_ratio_pct, ratio_diff_pp = round(ratio_diff_pp, 2))])

status <- data.table(
  check = c("(a) 대표 개체 tlast = 목표 중앙값", "(a') ka 격자 전체에서 tlast 유지", "(b) 참값 비율 진단(±3%p)",
            "(c) NCA 기반 AUClast/AUCinf (BEmaster)", "(d) IIV 반영 tlast 분포 (FDA 표)", "Li 2020 CV 35–51% (FDA 표+BEmaster)", "Cohen 2022 log SD 0.49 (FDA 표+BEmaster)"),
  status = c(if (all(tab[is_config_ka == TRUE, tlast_match])) "PASS" else "FAIL",
             if (all(tab$tlast_match)) "PASS" else "FAIL",
             if (all(abs(tab[is_config_ka == TRUE, ratio_diff_pp]) < 3)) "PASS(진단)" else "CHECK(진단)",
             if (bemaster_available()) "RUN" else "BLOCKED: BEmaster 미입수 (Q1)",
             if (all(p$status[grepl("^omega", item), status] == "confirmed")) "RUN" else "BLOCKED: FDA 표 미입수 (Q2)",
             "BLOCKED: Q1, Q2", "BLOCKED: Q1, Q2"))
fwrite(status, file.path(out_dir, "step1_status.csv"))
cat("\n단계 1 상태:\n"); print(status)
append_run_log(logfile, "step1 status written")
