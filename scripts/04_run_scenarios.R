#!/usr/bin/env Rscript
# 모든 시나리오 × 채혈 일정 실행 (SPEC §5.5, §7). 사용법: Rscript scripts/04_run_scenarios.R [dev|final]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "dev"
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
schedules <- names(design$schedules)
for (code in names(sc$main)) for (sch in schedules) {
  cat("\n=== ", code, " / ", sch, " ===\n")
  system2("Rscript", c("scripts/03_run_dev.R", run_mode, code, sch))
}
for (code in names(sc$sensitivity)) {
  st <- sc$sensitivity[[code]]$status
  if (!is.null(st) && st == "pending") { cat("\n[SKIP] ", code, ": 값 PENDING (SPEC Q8)\n"); next }
  system2("Rscript", c("scripts/03_run_dev.R", run_mode, code, design$default_schedule))
}
