#!/usr/bin/env Rscript
# 저장된 시험 수준 원자료에서 요약 표를 다시 만든다(시뮬레이션 재실행 없음). 사용법: Rscript scripts/14_postprocess.R
source("R/00_setup.R"); source_project()
design <- read_cfg("trial_design.yaml"); out_dir <- proj_path("results", "trials")
for (f in list.files(proj_path("results", "individual"), pattern = "^nca_.*_20000\\.rds$")) postprocess_individual_paired(sub("^nca_(.*)_20000\\.rds$", "\\1", f))
for (f in list.files(out_dir, pattern = "^schedules_be_raw_.*\\.csv(\\.gz)?$")) postprocess_schedules(sub("^schedules_be_raw_(.*)\\.csv(\\.gz)?$", "\\1", f), design)
for (f in list.files(out_dir, pattern = "^products_be_raw.*\\.csv(\\.gz)?$")) {
  v <- sub("^products_be_raw_?(.*)\\.csv(\\.gz)?$", "\\1", f); postprocess_products(if (v == "") "base" else v)
}
