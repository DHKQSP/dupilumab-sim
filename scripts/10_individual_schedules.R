#!/usr/bin/env Rscript
# 개인 수준 지표: 가상 대상자 20,000명, 일정별 (지시서 §5). 변형(구조·변동·잔차·체중·ADA)도 실행. 결과: results/individual/
# 사용법: Rscript scripts/10_individual_schedules.R [variants...] (기본: 전체)
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
variants <- if (length(args)) args else names(sc$variants)
MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "individual"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
scheds <- design$schedule_analysis
n <- design$mc$n_individual
logfile <- start_run_log("individual_schedules", master_seed = MASTER_SEED, run_mode = "final", extra = list(variants = variants, n = n))

all_summ <- list()
for (v in variants) {
  rv <- resolve_variant(v, design, sc); p <- rv$p
  cat(sprintf("\n=== variant %s (%s) ===\n", v, rv$label))
  t0 <- Sys.time()
  pop <- run_individual_population(n, p, design, scheds, MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = paste0("ind_", v))
  nca <- pop$nca
  summ <- rbindlist(lapply(scheds, function(sh) summarize_individual(nca[schedule == sh], last_planned = max(get_schedule(design, sh)))[, schedule := sh]))
  summ[, `:=`(variant = v, n_points = vapply(schedule, function(sh) length(get_schedule(design, sh)), numeric(1)))]
  setcolorder(summ, c("variant", "schedule", "n_points"))
  fwrite(summ, file.path(out_dir, sprintf("individual_%s.csv", v)))
  fwrite(summarize_n_lambda(nca), file.path(out_dir, sprintf("n_lambda_%s.csv", v)))
  if (v == "ada10") {
    fwrite(rbindlist(lapply(scheds, function(sh) summarize_individual(nca[schedule == sh], by = "ada", last_planned = 56)[, schedule := sh])),
           file.path(out_dir, "individual_ada10_by_subgroup.csv"))
  }
  if (v == "base") saveRDS(nca, file.path(out_dir, "nca_base_20000.rds"))
  all_summ[[v]] <- summ
  print(summ[, .(schedule, n_points, tlast_median, quant_at_last_pct = round(quant_at_last_pct, 1), lambda_ok_pct = round(lambda_ok_pct, 1), reliable_pct = round(reliable_pct, 1),
                 extrap_median = round(extrap_median, 2), extrap_p95 = round(extrap_p95, 1), extrap_gt20_pct = round(extrap_gt20_pct, 2),
                 extrap_true_median = round(extrap_true_median, 2), err_tlast_sd = round(err_tlast_sd, 4), err_inf_sd = round(err_inf_sd, 4), AUClast_logcv = round(AUClast_logcv, 1))])
  append_run_log(logfile, sprintf("variant %s done in %s", v, format(Sys.time() - t0)))
}
fwrite(rbindlist(all_summ), file.path(out_dir, "individual_all_variants.csv"))
