#!/usr/bin/env Rscript
# 46_lloq_sensitivity.R — 연구 LLOQ 범위 민감도 (추가 지시 2026-09-26 §2, config/prereg_20260926.yaml section2, config/assay.yaml, D-054).
# 같은 대상자·채혈 시각·잔차(y_raw)를 LLOQ 격자(0.02–0.5 mg/L)마다 다시 검열한다(쌍대). 잔차 변형 fixed(1차)와 scaled(2차: 가산 SD × LLOQ/0.078).
# 모드
#   individual <k2016|k2020>  scripts/10 모집단(20,000명, 60–90 kg, 같은 시드·일정 합집합) → B0 NCA를 LLOQ × 잔차 변형마다. 검사: fixed·0.078이 저장본
#                             (individual_<v>.csv B0 행, nca_<v>_20000.rds B0 대상자 수준)과 같아야 쓴다.
#   cliff                     scripts/34 절벽 모집단(기본 체중, 두 모델)에서 LLOQ별 절벽 끝과 현행 일정의 절벽 채혈점 수. 검사: 0.078 행이 cliff_subjects.rds·cliff_points.csv와 같아야 쓴다.
#   trial [cores] [from] [to] 주 모델(k2016) S00 + 경계(Vmax 양방향, F 하향), 시험 1–5,000, M0·M1. 검사: M0·fixed·0.078이 oc_rejudge_be_k2016과 같아야 쓴다(묶음마다).
# 산출: results/lloq/ (lloq_individual_<model>.csv, lloq_individual_nca_<model>.rds, lloq_individual_check_<model>.csv, lloq_cliff_*.csv, lloq_trials_be_k2016.csv.gz, lloq_trials_drop_k2016.csv.gz)
# 사용법: nice -n 15 Rscript scripts/46_lloq_sensitivity.R <individual k2016|individual k2020|cliff|trial [cores] [from] [to]>
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("individual", "cliff", "trial")) stop("사용법: Rscript scripts/46_lloq_sensitivity.R <individual <model>|cliff|trial [cores] [from] [to]>")
mode <- args[1]
pr <- read_cfg("prereg_20260926.yaml")$section2
design <- read_cfg("trial_design.yaml"); oc <- read_cfg("oc_design.yaml")
LL <- study_lloq_grid(); L0 <- as.numeric(pr$reference_lloq_mg_L)
if (!isTRUE(all.equal(LL, sort(as.numeric(unlist(pr$lloq_grid_mg_L)))))) stop("config/assay.yaml 격자와 사전 등록 격자가 다릅니다")
if (!isTRUE(all.equal(study_lloq(), L0))) stop("연구 LLOQ가 사전 등록 기준 LLOQ와 다릅니다(민감도는 기준 LLOQ에서 저장본을 재현해야 함)")
RES <- unlist(pr$residual_variants); stopifnot(identical(RES, c("fixed", "scaled")))
out_dir <- proj_path(pr$out_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
MASTER_SEED <- 20260923L
rel_eq <- function(a, b, tol) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & (a == b | abs(a / b - 1) <= tol))

if (mode == "individual") {
  model <- args[2]; stopifnot(model %in% unlist(pr$individual$pk_models))
  v <- if (model == "k2016") "base" else "struct2020"
  rv <- resolve_variant(v, design); p <- rv$p
  tag <- paste0("ind_", v); scheds <- design$schedule_analysis; n <- as.integer(design$mc$n_individual)
  logfile <- start_run_log(paste0("lloq_individual_", model), master_seed = MASTER_SEED, run_mode = "final", extra = list(variant = v, n = n, lloqs = LL, resid = RES, tag = tag))
  t0 <- Sys.time()
  pop <- run_individual_population(n, p, design, scheds, MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = tag)
  # 검사 1: 저장된 대상자 수준 B0 NCA
  st <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", v)))[schedule == "B0"]
  mine <- pop$nca[schedule == "B0"]
  setkey(st, id); setkey(mine, id)
  cols <- intersect(c("AUClast", "Cmax", "tlast", "pct_extrap", "reliable", "lambda_ok", "flag_rsq", "flag_extrap", "flag_span", "AUCinf", "AUCinf_true", "coverage_true", "pct_extrap_true"), names(st))
  ok_nca <- nrow(st) == nrow(mine) && all(vapply(cols, function(cn) isTRUE(all.equal(st[[cn]], mine[[cn]], tolerance = 0)), logical(1)))
  r <- lloq_individual(pop, design, p, LL, RES, MASTER_SEED, tag)
  # 검사 2: 요약(fixed, 0.078) = individual_<v>.csv B0 행
  ss <- fread(proj_path("results", "individual", sprintf("individual_%s.csv", v)))[schedule == "B0"]
  me <- r$summary[resid == "fixed" & abs(lloq - L0) < 1e-12]
  num <- setdiff(intersect(names(ss), names(me))[vapply(intersect(names(ss), names(me)), function(cn) is.numeric(ss[[cn]]), logical(1))], c("n_points"))
  ok_sum <- all(vapply(num, function(cn) rel_eq(me[[cn]], ss[[cn]], 1e-9), logical(1)))
  # 검사 3: scaled·0.078 = fixed·0.078(배율 1)
  a <- r$nca[resid == "scaled" & abs(lloq - L0) < 1e-12]; b <- r$nca[resid == "fixed" & abs(lloq - L0) < 1e-12]
  ok_sc <- isTRUE(all.equal(a[, !c("resid")], b[, !c("resid")], tolerance = 0))
  chk <- data.table(model = model, check = c("subject-level B0 NCA equals stored nca rds", "B0 summary equals stored individual csv (relative 1e-9)", "scaled residual at 0.078 equals fixed"),
                    n_items = c(length(cols) * nrow(st), length(num), nrow(a)), pass = c(ok_nca, ok_sum, ok_sc))
  fwrite(chk, file.path(out_dir, sprintf("lloq_individual_check_%s.csv", model)))
  print(chk)
  if (!all(chk$pass)) stop("저장본 재현 실패: 결과를 쓰지 않습니다")
  r$summary[, `:=`(model = model, variant = v, lambda_fail_pct = 100 - lambda_ok_pct)]
  setcolorder(r$summary, c("model", "variant", "resid", "lloq"))
  fwrite(r$summary, file.path(out_dir, sprintf("lloq_individual_%s.csv", model)))
  saveRDS(r$nca[, model := model][], file.path(out_dir, sprintf("lloq_individual_nca_%s.rds", model)))
  append_run_log(logfile, sprintf("done in %s; checks passed (%s)", format(Sys.time() - t0), paste(chk$check, collapse = "; ")))
}

if (mode == "cliff") {
  cf <- oc$cliff; SEED <- as.integer(cf$seed); STEP <- as.numeric(cf$grid_step_day); N <- as.integer(cf$n_subjects); defs <- as.numeric(unlist(cf$definition_days))
  b0 <- get_schedule(design, "B0")
  SCHED <- list(current = b0); SCHED_LABEL <- c(current = "현행")                      # scripts/34와 같은 이름·표기(저장본 대조용)
  min_interval <- vapply(SCHED, function(s) { x <- s[s >= 21]; min(diff(x)) }, numeric(1))
  logfile <- start_run_log("lloq_cliff", master_seed = SEED, run_mode = "final", extra = list(n = N, lloqs = LL))
  st_cs <- readRDS(proj_path("results", "cliff", "cliff_subjects.rds")); st_pts <- fread(proj_path("results", "cliff", "cliff_points.csv"))
  CS <- list(); PT <- list(); CK <- list()
  for (model in c("k2016", "k2020")) {
    t0 <- Sys.time()
    cs <- cliff_subjects(model, "base", N, SEED, STEP, LL, defs, cf$weights$base, design)
    m_ <- model
    a <- cs[abs(lloq - L0) < 1e-12][, lloq := NULL]; b <- st_cs[model == m_ & weight == "base"]
    setkey(a, id); setkey(b, id)
    ok1 <- identical(names(a), names(b)) && isTRUE(all.equal(as.data.frame(a), as.data.frame(b), tolerance = 0))
    pts <- rbindlist(lapply(LL, function(L) rbindlist(lapply(unlist(cf$timing), function(tm)
      cliff_count_points(cs[abs(lloq - L) < 1e-12], tm, model, "base", SCHED, SCHED_LABEL, min_interval, defs, SEED, design)))[, lloq := L]))
    p0 <- pts[abs(lloq - L0) < 1e-12][, lloq := NULL]
    s0 <- st_pts[model == m_ & weight == "base" & schedule == "current"]
    mm <- merge(p0, s0, by = c("model", "weight", "timing", "definition_day", "schedule"), suffixes = c("", ".s"))
    numc <- c("n_subjects", "pct_ge1", "pct_ge2", "pct_ge3", "ge1_lo", "ge1_hi", "n_cliff_in_day29_57")
    ok2 <- nrow(mm) == nrow(s0) && nrow(mm) == 4L && all(vapply(numc, function(cn) all(rel_eq(mm[[cn]], mm[[paste0(cn, ".s")]], 1e-12)), logical(1)))
    CK[[model]] <- data.table(model = model, check = c("cliff subjects at 0.078 equal stored cliff_subjects.rds (base)", "cliff points at 0.078 equal stored cliff_points.csv (current, base)"), pass = c(ok1, ok2))
    CS[[model]] <- cs[, .(n = .N, lloq_never_pct = 100 * mean(is.na(t_lloq)), lloq_studyday_median = median(t_lloq + 1, na.rm = TRUE),
                          lloq_studyday_p05 = quantile(t_lloq + 1, 0.05, na.rm = TRUE, names = FALSE), lloq_studyday_p95 = quantile(t_lloq + 1, 0.95, na.rm = TRUE, names = FALSE),
                          lloq_after_day58_pct = 100 * mean(t_lloq + 1 > 58, na.rm = TRUE), c_start1_below_lloq_pct = 100 * mean(c_start1 < lloq, na.rm = TRUE),
                          len1_median = median(len1, na.rm = TRUE), len1_p95 = quantile(len1, 0.95, na.rm = TRUE, names = FALSE), len2_median = median(len2, na.rm = TRUE)), by = .(model, lloq)]
    PT[[model]] <- pts
    append_run_log(logfile, sprintf("%s done in %s; checks %s", model, format(Sys.time() - t0), paste(CK[[model]]$pass, collapse = ",")))
  }
  ck <- rbindlist(CK); fwrite(ck, file.path(out_dir, "lloq_cliff_check.csv")); print(ck)
  if (!all(ck$pass)) stop("절벽 저장본 재현 실패: 결과를 쓰지 않습니다")
  fwrite(rbindlist(CS), file.path(out_dir, "lloq_cliff_summary.csv")); fwrite(rbindlist(PT), file.path(out_dir, "lloq_cliff_points.csv"))
}

if (mode == "trial") {
  tp <- pr$trial; model <- tp$pk_model; stopifnot(model == "k2016")
  cores <- if (length(args) >= 2) as.integer(args[2]) else 2L
  tr <- as.integer(unlist(tp$trials))
  trial_from <- if (length(args) >= 3) as.integer(args[3]) else tr[1]; trial_to <- if (length(args) >= 4) as.integer(args[4]) else tr[2]
  stopifnot(trial_from >= tr[1], trial_to <= tr[2], trial_to >= trial_from)
  batch <- as.integer(oc$trials$batch); OC_SEED <- as.integer(oc$trials$master_seed)
  rv <- resolve_variant("base", design); p <- rv$p
  sc <- fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", model)))
  SCN <- unlist(tp$scenarios); stopifnot(SCN[1] == "S00")
  scen <- c(list(S00 = list(code = "S00", T_multipliers = list())),
            setNames(lapply(SCN[-1], function(cd) { r <- sc[code == cd]; stopifnot(nrow(r) == 1); list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }), SCN[-1]))
  EPS <- unlist(tp$endpoints); AMS <- unlist(tp$analysis_models)
  n_rows_trial <- length(LL) * length(RES) * length(SCN) * length(EPS) * length(AMS)
  stored <- fread(proj_path("results", "oc", sprintf("oc_rejudge_be_%s.csv.gz", model)))[trial >= trial_from & trial <= trial_to & scenario %in% SCN & endpoint %in% EPS]
  if (nrow(stored) != (trial_to - trial_from + 1L) * length(SCN) * length(EPS)) stop("oc_rejudge_be 저장본이 요청 범위를 모두 갖고 있지 않습니다")
  be_f <- file.path(out_dir, sprintf("lloq_trials_be_%s.csv.gz", model)); dr_f <- file.path(out_dir, sprintf("lloq_trials_drop_%s.csv.gz", model))
  COLS <- c("trial", "lloq", "resid", "scenario", "endpoint", "model", "est", "se", "pass", "n_R", "n_T")
  done <- integer(0)
  if (file.exists(be_f)) {
    ex <- fread(be_f); if (!identical(names(ex), COLS)) stop("기존 파일 열 불일치: ", be_f)
    cn <- ex[, .N, by = trial]; bad <- cn[N != n_rows_trial, trial]
    if (length(bad) || anyDuplicated(ex, by = c("trial", "lloq", "resid", "scenario", "endpoint", "model"))) {
      ex <- unique(ex[!trial %in% bad], by = c("trial", "lloq", "resid", "scenario", "endpoint", "model")); tmp <- paste0(be_f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, be_f) }
    done <- unique(ex$trial); rm(ex)
    if (file.exists(dr_f)) { dr <- fread(dr_f); dr <- dr[trial %in% done]; tmp <- paste0(dr_f, ".rewrite.csv.gz"); fwrite(dr, tmp); file.rename(tmp, dr_f) }
  }
  append_gz <- function(dt, f) { tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f)); fwrite(dt, tmp, col.names = !file.exists(f))
    if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f) }
  logfile <- start_run_log(paste0("lloq_trials_", model), master_seed = OC_SEED, run_mode = "final",
                           extra = list(model = model, scenarios = paste(SCN, collapse = ","), lloqs = LL, resid = RES, trial_from = trial_from, trial_to = trial_to, cores = cores))
  say <- function(msg) { cat(msg, "\n"); append_run_log(logfile, msg) }
  invisible(get_model(p$model_id))
  one <- function(j) run_trial_oc_models(j, p, design, scen, OC_SEED, rv$wt_spec, model_id = p$model_id, lloqs = LL, resid = RES, models = AMS, endpoints = EPS)
  n_chk <- 0L
  for (b0 in seq(trial_from, trial_to, by = batch)) {
    ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done)
    if (!length(ids)) next
    t0 <- Sys.time()
    res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
    bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
    if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
    be <- rbindlist(lapply(res, `[[`, "be")); dr <- rbindlist(lapply(res, `[[`, "drop"))
    stopifnot(nrow(be) == length(ids) * n_rows_trial)
    m0 <- be[model == "M0" & resid == "fixed" & abs(lloq - L0) < 1e-12]
    q <- qt(1 - (1 - design$be$ci_level) / 2, m0$df)
    m0 <- m0[, .(trial, scenario, endpoint, GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se), pass, n_R, n_T)]
    mm <- merge(m0, stored[trial %in% ids], by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
    ok <- nrow(mm) == nrow(m0) && nrow(mm) == nrow(stored[trial %in% ids]) &&
      all((rel_eq(mm$GMR, mm$GMR.s, 1e-6) & rel_eq(mm$CI_lower, mm$CI_lower.s, 1e-6) & rel_eq(mm$CI_upper, mm$CI_upper.s, 1e-6) &
             ((is.na(mm$pass) & is.na(mm$pass.s)) | (mm$pass %in% TRUE & mm$pass.s %in% TRUE) | (mm$pass %in% FALSE & mm$pass.s %in% FALSE)) & mm$n_R == mm$n_R.s & mm$n_T == mm$n_T.s) %in% TRUE)
    if (!ok) { print(head(mm, 5)); stop(sprintf("M0·fixed·%g가 oc_rejudge_be_%s와 다릅니다(시험 %d-%d). 결과를 쓰지 않습니다.", L0, model, min(ids), max(ids))) }
    n_chk <- n_chk + nrow(mm)
    be[, `:=`(est = signif(est, 8), se = signif(se, 8))]
    append_gz(dr, dr_f); append_gz(be[, ..COLS], be_f)
    say(sprintf("trials %d-%d (%d LLOQs x %d residual variants x %d scenarios) done in %s; M0 at %g matched oc_rejudge_be", min(ids), max(ids), length(LL), length(RES), length(SCN), format(Sys.time() - t0), L0))
  }
  say(sprintf("done: %d stored rows matched; output %s", n_chk, be_f))
}
