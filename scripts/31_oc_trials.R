#!/usr/bin/env Rscript
# §3-5 역산 시나리오의 시험 반복 (config/oc_design.yaml, D-040). B0, arm당 117명, pooled t, 대조군 공통 난수.
# 시험 1–reps_other: 도달한 모든 시나리오 + 동일 제품(S00). 시험 reps_other+1 – reps_boundary: 경계(0.80, 1.25) 시나리오 + S00.
# 500회 묶음 체크포인트(중단 후 재시작 시 완료 묶음 건너뜀). 사용법: Rscript scripts/31_oc_trials.R <k2016|k2020> [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); model <- args[1]; cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p
out_dir <- proj_path("results", "oc")
inv <- rbindlist(lapply(list.files(out_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", model), full.names = TRUE), fread))
stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
inv <- inv[reachable == TRUE]
inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
scen_all <- c(list(S00 = list(code = "S00", T_multipliers = list())),
              setNames(lapply(seq_len(nrow(inv)), function(i) list(code = inv$code[i], T_multipliers = setNames(list(inv$multiplier[i]), inv$mechanism[i]))), inv$code))
bnd <- as.numeric(unlist(oc$boundary_targets))
scen_bnd <- scen_all[c("S00", inv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9, code])]
n_other <- as.integer(oc$trials$reps_other); n_bnd <- as.integer(oc$trials$reps_boundary); batch <- as.integer(oc$trials$batch)
MASTER_SEED <- as.integer(oc$trials$master_seed)
raw_f <- file.path(out_dir, sprintf("oc_trials_be_%s.csv.gz", model)); drop_f <- file.path(out_dir, sprintf("oc_trials_drop_%s.csv.gz", model))
logfile <- start_run_log(paste0("oc_trials_", model), master_seed = MASTER_SEED, run_mode = "final",
                         extra = list(model = model, n_scen_all = length(scen_all), n_scen_bnd = length(scen_bnd), reps_other = n_other, reps_boundary = n_bnd, cores = cores))
done <- if (file.exists(raw_f)) unique(fread(raw_f, select = c("trial", "scenario"))) else data.table(trial = integer(0), scenario = character(0))
run_batch <- function(ids, scen) {
  ids <- setdiff(ids, done[scenario == names(scen)[length(scen)], trial])   # 마지막 시나리오까지 저장된 시험은 완료
  if (!length(ids)) return(invisible())
  t0 <- Sys.time()
  r <- run_trials_oc(ids, p, design, scen, MASTER_SEED, rv$wt_spec, cores = cores, model_id = p$model_id, progress_every = 0)
  r$be[, `:=`(GMR = signif(GMR, 7), CI_lower = signif(CI_lower, 7), CI_upper = signif(CI_upper, 7))]
  fwrite(r$be, raw_f, append = file.exists(raw_f)); fwrite(r$drop, drop_f, append = file.exists(drop_f))
  msg <- sprintf("trials %d–%d (%d scenarios) done in %s", min(ids), max(ids), length(scen), format(Sys.time() - t0))
  cat(msg, "\n"); append_run_log(logfile, msg)
}
for (b0 in seq(1L, n_other, by = batch)) run_batch(b0:min(b0 + batch - 1L, n_other), scen_all)
for (b0 in seq(n_other + 1L, n_bnd, by = batch)) run_batch(b0:min(b0 + batch - 1L, n_bnd), scen_bnd)
fwrite(inv, file.path(out_dir, sprintf("oc_scenarios_%s.csv", model)))
append_run_log(logfile, "done")
