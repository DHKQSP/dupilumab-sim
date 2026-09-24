#!/usr/bin/env Rscript
# 주 모델 제품 시나리오 5,000회 (검토 의견 통합본 §4-2). B0, 60–90 kg, arm당 117명, pooled t.
# 500회 단위로 원자료를 저장(체크포인트)하고, 끝나면 적응적 상향 규칙을 적용한다:
#   인용 대상 비율의 Wilson 95% 구간이 임계(검정력 90%, 소비자 위험 5%)를 포함하면 해당 시나리오만 10,000 → 20,000회로 늘린다.
# 사용법: Rscript scripts/20_products_5000.R [cores] [variant]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
variant <- if (length(args) >= 2) args[2] else "base"
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); rv <- resolve_variant(variant, design, sc); p <- rv$p
MASTER_SEED <- 20260923L
SCEN <- c("S00", "F085", "VM150", "VM125", "KE110", "KE120", "F110", "F097", "F090")   # KE120 = 500회에서 "AUClast 통과·AUCinf 불통과" 최대
N_BASE <- 5000L; BATCH <- 500L; STEPS <- c(10000L, 20000L); THRESH <- c(90, 5)
out_dir <- proj_path("results", "trials5000"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
raw_f <- file.path(out_dir, sprintf("products5000_be_raw_%s.csv.gz", variant))
ind_f <- file.path(out_dir, sprintf("products5000_ind_%s.csv.gz", variant))
logfile <- start_run_log(paste0("products5000_", variant), master_seed = MASTER_SEED, run_mode = "final",
                         extra = list(variant = variant, scenarios = SCEN, n_base = N_BASE, batch = BATCH, cores = cores))
done_ids <- if (file.exists(raw_f)) unique(fread(raw_f)[, trial]) else integer(0)

run_batch <- function(ids, scen_names) {
  combos <- CJ(scenario = scen_names, schedule = "B0")
  r <- run_trials(length(ids), p, design, sc$scenarios[scen_names], combos, MASTER_SEED, rv$wt_spec, jitter = TRUE, methods = "pooled_t",
                  model_id = p$model_id, cores = cores, progress_every = 0, trial_ids = ids)
  fwrite(r$be, raw_f, append = file.exists(raw_f)); fwrite(r$ind, ind_f, append = file.exists(ind_f))
  append_run_log(logfile, sprintf("batch %d–%d (%s) done", min(ids), max(ids), paste(scen_names, collapse = ",")))
}
# 1) 기본 5,000회 (중단 후 재시작 시 완료된 시험은 건너뜀)
for (b0 in seq(1L, N_BASE, by = BATCH)) {
  ids <- setdiff(b0:(b0 + BATCH - 1L), done_ids)
  if (length(ids)) { t0 <- Sys.time(); run_batch(ids, SCEN); cat(sprintf("trials %d–%d done in %s\n", min(ids), max(ids), format(Sys.time() - t0))) }
}
# 2) 적응적 상향
cited_props <- function(be) {
  w <- dcast(be, scenario + trial ~ endpoint, value.var = "pass")
  w[, `:=`(joint_last_cmax = AUClast & Cmax, joint_3rel = AUClast & Cmax & AUCinf_reliable, joint_3all = AUClast & Cmax & AUCinf_all,
           discord_last_pass_inf_fail = AUClast & !AUCinf_reliable)]
  long <- melt(w, id.vars = c("scenario", "trial"), variable.name = "metric", value.name = "ok", variable.factor = FALSE)
  long[, .(x = sum(ok %in% TRUE), n = sum(!is.na(ok))), by = .(scenario, metric)][, c("est", "lo", "hi") := wilson_ci(x, n)[, .(est, lo, hi)]][]
}
for (step in STEPS) {
  be <- fread(raw_f); cp <- cited_props(be)
  cover <- cp[(lo <= THRESH[1] & hi >= THRESH[1]) | (lo <= THRESH[2] & hi >= THRESH[2])]
  need <- unique(cover$scenario)
  cur_n <- be[, .(n = uniqueN(trial)), by = scenario]
  need <- need[cur_n[match(need, scenario), n] < step]
  if (!length(need)) break
  cat(sprintf("\n적응적 상향 → %d회: %s (임계 포함 지표: %s)\n", step, paste(need, collapse = ", "), paste(unique(cover[scenario %in% need, paste(scenario, metric)]), collapse = "; ")))
  append_run_log(logfile, sprintf("escalate to %d: %s", step, paste(need, collapse = ",")))
  have <- max(be$trial)
  for (b0 in seq(have + 1L, step, by = BATCH)) run_batch(b0:min(b0 + BATCH - 1L, step), c("S00", need)[!duplicated(c("S00", need))])
}
be <- fread(raw_f); cp <- cited_props(be)
reps <- be[, .(n_trials = uniqueN(trial)), by = scenario]
fwrite(merge(cp, reps, by = "scenario"), file.path(out_dir, sprintf("products5000_props_%s.csv", variant)))
cat("\n최종 반복 수:\n"); print(reps)
append_run_log(logfile, "done")
