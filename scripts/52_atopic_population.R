#!/usr/bin/env Rscript
# 52_atopic_population.R — 아토피 피부염 성인 체중 분포 (두 번째 추가 지시 2026-09-26 §2, prereg section5, config/population_atopic.yaml, D-058).
# 모드
#   individual            모델 변형 3개(2016, Model 1, 2016 + ke~BMI + Vc~체중 0.817) × 분포 3개(주, 민감도 2개) × 20,000명, 300 mg, B0.
#                         대상자 수준 NCA(results/atopic/atopic_nca_<variant>_<dist>.rds, 저장소 제외)와 요약(atopic_individual_*.csv)
#   truth                 주 분포 20,000명에서 시나리오별 참 AUC0-inf·Cmax 비(공통 난수), 모델 변형별
#   trial <variant> [cores] [from] [to]  주 분포, arm당 117명, 시험 2,000회, S00·F090·KE120·VM125·KM10, M0·M1, 평가변수 12개(기준 세트 i–iv)
# 사용법: nice -n 15 Rscript scripts/52_atopic_population.R <individual|truth|trial <variant> [cores]>
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("individual", "truth", "trial")) stop("사용법: Rscript scripts/52_atopic_population.R <individual|truth|trial <variant> [cores] [from] [to]>")
mode <- args[1]
pr <- read_cfg("prereg_20260926.yaml")$section5; pa <- read_cfg("population_atopic.yaml"); design <- read_cfg("trial_design.yaml"); oc <- read_cfg("oc_design.yaml")
out_dir <- Sys.getenv("DUPI_ATOPIC_OUT", proj_path("results", "atopic")); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)   # 환경변수는 시험용
production <- normalizePath(out_dir, mustWork = FALSE) == normalizePath(proj_path("results", "atopic"), mustWork = FALSE)
MASTER_SEED <- 20260923L; VARIANTS <- unlist(pr$model_variants); DISTS <- unlist(pr$distributions)
stopifnot(identical(DISTS, names(pa$distributions)))
spec_of <- function(d) { x <- pa$distributions[[d]]; stopifnot(identical(x$dist, "lognormal"))
  list(dist = "lognormal", mean = as.numeric(x$mean), sd = as.numeric(x$sd), trunc = as.numeric(unlist(x$trunc)),
       height = list(mean = as.numeric(pa$height$mean), sd = as.numeric(pa$height$sd), trunc = as.numeric(unlist(pa$height$trunc)))) }
log_or_null <- function(name, seed, extra) if (production) start_run_log(name, master_seed = seed, run_mode = "final", extra = extra) else NULL
say <- function(lf, msg) { cat(msg, "\n"); if (!is.null(lf)) append_run_log(lf, msg) }

if (mode == "individual") {
  N <- as.integer(Sys.getenv("DUPI_ATOPIC_TEST_N", as.character(pr$individual$n_subjects))); stopifnot(as.numeric(pr$individual$dose_mg) == design$dose_mg)   # 환경변수는 시험용
  lf <- log_or_null("atopic_individual", MASTER_SEED, list(variants = VARIANTS, dists = DISTS, n = N))
  CS <- list(); P1 <- list(); WS <- list()
  for (d in DISTS) for (v in VARIANTS) {
    t0 <- Sys.time()
    rv <- resolve_variant(v, design); p <- rv$p
    pop <- run_individual_population(N, p, design, "B0", MASTER_SEED, spec_of(d), jitter = TRUE, sex_ratio_male = as.numeric(pa$sex_ratio_male), model_id = p$model_id, tag = paste0("atopic_", d))
    x <- pop$nca[schedule == "B0"]
    stopifnot(nrow(x) == N, all(c("WT", "BMI", "Rsq_adjusted", "Span_ratio", "AUC_%Extrap_obs", "AUCinf_true", "coverage_true") %in% names(x)))
    for (s_ in names(CRIT_SETS)) x[, (paste0("ok_", s_)) := crit_ok(x, s_)]
    x[, band := weight_band(WT)]
    saveRDS(x[, .(id, WT, HT, BMI, band, AUClast, Cmax, AUCinf, AUCinf_true, coverage_true, pct_extrap_true, pct_extrap, tlast, lambda_ok, Rsq_adjusted, Span_ratio,
                  `AUC_%Extrap_obs`, ok_i, ok_ii, ok_iii, ok_iv)], file.path(out_dir, sprintf("atopic_nca_%s_%s.rds", v, d)))
    tag <- list(variant = v, distribution = d)
    CS[[paste(v, d)]] <- rbind(crit_summary(x)[, band := "all"], crit_summary(x, by = "band"), fill = TRUE)[, names(tag) := tag]
    P1[[paste(v, d)]] <- rbind(pillar1_summary(x)[, band := "all"], pillar1_summary(x, by = "band"), fill = TRUE)[, names(tag) := tag]
    if (v == VARIANTS[1]) WS[[d]] <- data.table(distribution = d, n = N, below_60 = 100 * mean(x$WT < 60), in_60_90 = 100 * mean(x$WT >= 60 & x$WT <= 90), above_90 = 100 * mean(x$WT > 90),
                                                above_100 = 100 * mean(x$WT > 100), wt_mean = mean(x$WT), wt_sd = sd(x$WT), wt_median = median(x$WT), bmi_mean = mean(x$BMI), bmi_sd = sd(x$BMI))
    say(lf, sprintf("%s %s done in %s", v, d, format(Sys.time() - t0)))
  }
  fwrite(rbindlist(CS, use.names = TRUE), file.path(out_dir, "atopic_individual_criteria.csv"))
  fwrite(rbindlist(P1, use.names = TRUE), file.path(out_dir, "atopic_individual_pillar1.csv"))
  fwrite(rbindlist(WS), file.path(out_dir, "atopic_weight_simulated.csv"))
}

if (mode == "truth") {
  N <- 20000L; SCN <- setdiff(unlist(pr$trial$scenarios), "S00"); sc <- load_scenarios()$scenarios
  lf <- log_or_null("atopic_truth", 20260924L, list(variants = VARIANTS, scenarios = SCN, n = N))
  TR <- rbindlist(lapply(VARIANTS, function(v) {
    rv <- resolve_variant(v, design); p <- rv$p
    subj <- with_seed(derive_seed(as.integer(oc$truth_seed), "atopic_truth", v), make_subjects(N, p, spec_of("primary"), as.numeric(pa$sex_ratio_male), 0))
    ref <- truth_metrics(individual_params(p, subj), design$dose_mg, p$model_id, cmax = TRUE)
    r <- rbindlist(lapply(SCN, function(s) { tr <- truth_ratio(ref, truth_metrics(ip_for_mult(p, subj, sc[[s]]$T_multipliers), design$dose_mg, p$model_id, cmax = TRUE))
      data.table(variant = v, scenario = s, auc_ratio = tr$auc_ratio, auc_se_log = tr$auc_se_log, cmax_ratio = tr$cmax_ratio, cmax_se_log = tr$cmax_se_log, n_subjects = tr$n) }))
    say(lf, sprintf("%s truth done", v)); rbind(data.table(variant = v, scenario = "S00", auc_ratio = 1, auc_se_log = 0, cmax_ratio = 1, cmax_se_log = 0, n_subjects = N), r)
  }))
  fwrite(TR, file.path(out_dir, "atopic_truth.csv")); print(TR)
}

if (mode == "trial") {
  v <- args[2]; stopifnot(v %in% VARIANTS)
  cores <- if (length(args) >= 3) as.integer(args[3]) else 2L
  NT <- as.integer(pr$trial$trials); trial_from <- if (length(args) >= 4) as.integer(args[4]) else 1L; trial_to <- if (length(args) >= 5) as.integer(args[5]) else NT
  rv <- resolve_variant(v, design); p <- rv$p; spec <- spec_of(pr$trial$distribution)
  SCN <- unlist(pr$trial$scenarios); scen <- load_scenarios()$scenarios[SCN]; stopifnot(length(scen) == length(SCN))
  AMS <- unlist(pr$trial$analysis_models); EPS <- c(OC_ENDPOINTS_EXT, OC_ENDPOINTS_CRIT)
  n_rows_trial <- length(SCN) * length(EPS) * length(AMS); batch <- 250L
  be_f <- file.path(out_dir, sprintf("atopic_trials_be_%s.csv.gz", v)); dr_f <- file.path(out_dir, sprintf("atopic_trials_drop_%s.csv.gz", v))
  COLS <- c("trial", "scenario", "endpoint", "model", "est", "se", "pass", "n_R", "n_T")
  done <- integer(0)
  if (file.exists(be_f)) { ex <- fread(be_f); cn <- ex[, .N, by = trial]; done <- cn[N == n_rows_trial, trial]
    if (length(setdiff(unique(ex$trial), done))) { ex <- ex[trial %in% done]; tmp <- paste0(be_f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, be_f) }
    if (file.exists(dr_f)) { dr <- fread(dr_f)[trial %in% done]; tmp <- paste0(dr_f, ".rewrite.csv.gz"); fwrite(dr, tmp); file.rename(tmp, dr_f) } }
  append_gz <- function(dt, f) { tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f)); fwrite(dt, tmp, col.names = !file.exists(f))
    if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f) }
  lf <- log_or_null(paste0("atopic_trials_", v), MASTER_SEED, list(variant = v, scenarios = paste(SCN, collapse = ","), trials = sprintf("%d-%d", trial_from, trial_to), cores = cores,
                                                                   distribution = pr$trial$distribution))
  invisible(get_model(p$model_id))
  one <- function(j) run_trial_oc_models(j, p, design, scen, MASTER_SEED, spec, model_id = p$model_id, models = AMS, endpoints = EPS)
  for (b0 in seq(trial_from, trial_to, by = batch)) {
    ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done); if (!length(ids)) next
    t0 <- Sys.time()
    res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
    bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
    if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
    be <- rbindlist(lapply(res, `[[`, "be")); dr <- rbindlist(lapply(res, `[[`, "drop")); st <- rbindlist(lapply(res, `[[`, "strata"))
    stopifnot(nrow(be) == length(ids) * n_rows_trial, identical(sort(unique(st$stratum)), sort(c("40-75", ">75-180"))))
    be[, `:=`(est = signif(est, 8), se = signif(se, 8))]
    append_gz(dr[, !c("lloq", "resid")], dr_f); append_gz(be[, ..COLS], be_f)
    say(lf, sprintf("%s trials %d-%d done in %s", v, min(ids), max(ids), format(Sys.time() - t0)))
  }
  say(lf, sprintf("done: %s", be_f))
}
