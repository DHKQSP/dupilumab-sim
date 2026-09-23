#!/usr/bin/env Rscript
# 개발 규모 가상 시험 루프 (SPEC §5.2: 바깥 20 × 안쪽 50). 사용법: Rscript scripts/03_run_dev.R [dev|final] [scenario] [schedule]
# BEmaster가 없으면 NCA 경계에서 멈추고 참값 기반 진단(§5.3 일부)만 남긴다. 결과 보고에는 쓰지 않는다.
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "dev"
scen_code <- if (length(args) >= 2) args[2] else "S0"
schedule  <- if (length(args) >= 3) args[3] else NULL
MASTER_SEED <- 20260923L

p <- load_params(run_mode); if (run_mode == "final") assert_final_ready(p)
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
scen <- c(sc$main, sc$sensitivity)[[scen_code]]; if (is.null(scen)) stop("알 수 없는 시나리오: ", scen_code)
mc <- design$mc[[run_mode]]
if (is.null(mc$n_outer)) stop("최종 실행 규모(n_outer, n_inner)가 지정되지 않았습니다 (SPEC §5.2)")
logfile <- start_run_log(paste0("mc_", scen_code), master_seed = MASTER_SEED, run_mode = run_mode,
                         extra = list(scenario = scen_code, schedule = if (is.null(schedule)) design$default_schedule else schedule,
                                      n_outer = mc$n_outer, n_inner = mc$n_inner, bemaster = bemaster_version()))
has_be <- bemaster_available()
if (!has_be) message("BEmaster 미입수: NCA/BE 생략, 참값 진단만 산출 (SPEC Q1)")
res <- run_mc(scen, p, design, n_outer = mc$n_outer, n_inner = mc$n_inner, master_seed = MASTER_SEED,
              schedule_name = schedule, nca_fn = if (has_be) run_nca_bemaster else NULL,
              be_fn = if (has_be) run_be_bemaster else NULL, do_be = has_be)
out_dir <- proj_path("results", run_mode); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- paste0(scen_code, "_", res$scenario, "_", if (is.null(schedule)) design$default_schedule else schedule)
saveRDS(res, file.path(out_dir, paste0("mc_", tag, ".rds")))
fwrite(summarize_truth(res$subjects), file.path(out_dir, paste0("truth_summary_", tag, ".csv")))
if (has_be) {
  fwrite(summarize_individual(res$nca), file.path(out_dir, paste0("individual_", tag, ".csv")))
  st <- summarize_trials(res$be, design$primary_endpoint_sets)
  fwrite(st$agreement, file.path(out_dir, paste0("agreement_", tag, ".csv")))
  fwrite(st$power, file.path(out_dir, paste0("power_", tag, ".csv")))
}
append_run_log(logfile, "done: ", tag)
print(summarize_truth(res$subjects))
