#!/usr/bin/env Rscript
# §5-2 Pillar 2 곡률 견고성: 모집단 Vmax ×0.8(양 군 공통)에서 제품 시나리오 5개 × 2,000회, B0, 60–90 kg, arm당 117명.
# 시험군 Vmax ×1.25는 모집단 ×0.8 위에 곱한다(= 0.8 × 1.25).
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); cores <- if (length(args)) as.integer(args[1]) else 4L
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); rv <- resolve_variant("vmax080_both", design, sc); p <- rv$p
MASTER_SEED <- 20260923L; N <- 2000L; SCEN <- c("S00", "F090", "KE120", "VM125", "KM10")
out_dir <- proj_path("results", "pillar2_curvature"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("pillar2_curvature", master_seed = MASTER_SEED, run_mode = "final", extra = list(n = N, scenarios = SCEN, cores = cores))
res <- run_trials(N, p, design, sc$scenarios[SCEN], CJ(scenario = SCEN, schedule = "B0"), MASTER_SEED, rv$wt_spec, methods = "pooled_t", model_id = p$model_id, cores = cores, progress_every = 250)
fwrite(res$be, file.path(out_dir, "pillar2_vmax080_be_raw.csv.gz"))
st <- summarize_trials_ci(res$be)
fwrite(st$per_endpoint, file.path(out_dir, "pillar2_vmax080_per_endpoint.csv")); fwrite(st$concordance, file.path(out_dir, "pillar2_vmax080_concordance.csv"))
print(st$per_endpoint[endpoint %in% c("AUClast", "AUCinf_reliable", "AUCinf_true", "Cmax")]); print(st$concordance)
append_run_log(logfile, "done")
