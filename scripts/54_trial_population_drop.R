#!/usr/bin/env Rscript
# 54_trial_population_drop.R — 시험 모집단(건강인 60–90 kg, 체중 층화) 안의 AUC0-inf 탈락: 층별·처리군별 (정정 지시 2026-09-26 §2-2, §2-3; prereg section6).
# 모드
#   truth  <model>                         F097·KE110·KE120·VM125의 참 AUC0-inf·Cmax 비(200,000명 공통 난수, section1의 F097과 같은 방법)
#   trials <model> [cores] [from] [to]     시험 1–n_trials_regen(2,000) 재생성. 시나리오 S00 + 경계 8개 + F097 + KE110 + KE120 + VM125.
#                                          시험 × 군(REF = 대조군, 시험군은 시나리오) × 배정 층마다 인원: 무작위배정, AUClast 산출, λz 산출, 기준 세트 (i)–(iv) 충족.
# 시드: section1(scripts/44)과 같은 규칙(oc_design.yaml trials$master_seed, simulate_trial_arms). 시험군 대상자는 시나리오 간 같다(공통 난수).
# 검사(묶음마다, 쓰기 전): S00·F097·경계 칸의 군별 인원(층 합계)이 section1 산출(results/oc_models/oc_models_be_<model>.csv.gz, oc_models_crit_be_<model>.csv.gz)의
#   M0 행 n_R·n_T와 같아야 한다: AUCinf_Ai = 세트 (i), AUCinf_A = (ii), AUCinf_Aiii = (iii), AUCinf_Aiv = (iv), AUCinf_B = λz 산출, AUClast = AUClast 산출. 다르면 중단.
# 산출(results/trialpop/): trialpop_counts_<model>.csv.gz, trialpop_truth_<model>.csv
# 사용법: nice -n 19 Rscript scripts/54_trial_population_drop.R <truth|trials> <k2016|k2020> [cores] [from] [to]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || !args[1] %in% c("truth", "trials") || !args[2] %in% c("k2016", "k2020")) stop("사용법: Rscript scripts/54_trial_population_drop.R <truth|trials> <k2016|k2020> [cores] [from] [to]")
mode <- args[1]; model <- args[2]
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml"); pr6 <- read_cfg("prereg_20260926.yaml")$section6; pr1 <- read_cfg("prereg_20260926.yaml")$section1
out_dir <- Sys.getenv("DUPI_TRIALPOP_OUT", proj_path("results", "trialpop")); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)   # 환경변수는 시험용
production <- normalizePath(out_dir, mustWork = FALSE) == normalizePath(proj_path("results", "trialpop"), mustWork = FALSE)
oc_dir <- proj_path("results", "oc"); om_dir <- proj_path(pr1$out_dir)
MASTER_SEED <- as.integer(oc$trials$master_seed); stopifnot(MASTER_SEED == as.integer(pr1$master_seed))
variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p
EXTRA <- c("F097", "KE110", "KE120", "VM125")
sc_yaml <- load_scenarios()$scenarios; stopifnot(all(EXTRA %in% names(sc_yaml)))
log_or_null <- function(name, seed, extra) if (production) start_run_log(name, master_seed = seed, run_mode = "final", extra = extra) else NULL
say <- function(lf, msg) { cat(msg, "\n"); if (!is.null(lf)) append_run_log(lf, msg) }

if (mode == "truth") {
  N <- as.integer(oc$estimand$population$n_subjects); stopifnot(N == 200000L)
  lf <- log_or_null(paste0("trialpop_truth_", model), as.integer(oc$truth_seed), list(model = model, scenarios = EXTRA, n = N))
  subj <- truth_subjects(N, p, design, derive_seed(as.integer(oc$truth_seed), "truth200k", model))
  ref <- readRDS(file.path(oc_dir, sprintf("truth_ref_%s.rds", model)))
  TR <- rbindlist(lapply(EXTRA, function(s) {
    te <- truth_metrics(ip_for_mult(p, subj, sc_yaml[[s]]$T_multipliers), as.numeric(oc$estimand$dose_mg), p$model_id, cmax = TRUE)
    r <- truth_ratio(ref, te); say(lf, sprintf("%s %s truth done", model, s))
    data.table(pk_model = model, scenario = s, auc_ratio = r$auc_ratio, auc_se_log = r$auc_se_log, cmax_ratio = r$cmax_ratio, cmax_se_log = r$cmax_se_log, n_subjects = r$n)
  }))
  f97 <- file.path(om_dir, "truth_F097.csv")                         # section1 요약(scripts/45)이 먼저 만들었으면 같은 값이어야 한다
  if (file.exists(f97)) { a <- fread(f97)[pk_model == model]; if (nrow(a) == 1 && abs(a$auc_ratio / TR[scenario == "F097", auc_ratio] - 1) > 1e-10) stop("F097 참값이 scripts/45와 다릅니다") }
  fwrite(TR, file.path(out_dir, sprintf("trialpop_truth_%s.csv", model))); print(TR)
  say(lf, "done")
}

if (mode == "trials") {
  cores <- if (length(args) >= 3) as.integer(args[3]) else 2L
  NT <- as.integer(pr6$n_trials_regen); trial_from <- if (length(args) >= 4) as.integer(args[4]) else 1L; trial_to <- if (length(args) >= 5) as.integer(args[5]) else NT
  stopifnot(trial_from >= 1, trial_to <= NT, trial_to >= trial_from)
  # 경계 8칸: scripts/44와 같은 규칙(역산 도달 행, 코드 = 기전_방향_목표×100)
  inv <- rbindlist(lapply(list.files(oc_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", model), full.names = TRUE), fread))[reachable == TRUE]
  inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
  bnd <- as.numeric(unlist(oc$boundary_targets)); bcodes <- inv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9, code]; stopifnot(length(bcodes) == 8L)
  bscen <- setNames(lapply(bcodes, function(cd) { r <- inv[code == cd]; list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }), bcodes)
  scen <- c(list(S00 = list(T_multipliers = list())), bscen, lapply(setNames(nm = EXTRA), function(s) list(T_multipliers = sc_yaml[[s]]$T_multipliers)))
  SEC1 <- c("S00", bcodes, "F097")
  # section1 대조 자료: M0 행의 n_R·n_T
  EPMAP <- c(AUClast = "n_auclast", AUCinf_B = "n_lambda", AUCinf_Ai = "n_i", AUCinf_A = "n_ii", AUCinf_Aiii = "n_iii", AUCinf_Aiv = "n_iv")
  s1 <- rbind(fread(file.path(om_dir, sprintf("oc_models_be_%s.csv.gz", model))), fread(file.path(om_dir, sprintf("oc_models_crit_be_%s.csv.gz", model))))[model == "M0" & endpoint %in% names(EPMAP) & trial >= trial_from & trial <= trial_to]
  if (!all(trial_from:trial_to %in% s1$trial)) stop("section1 산출에 대조할 시험이 아직 없습니다: ", paste(head(setdiff(trial_from:trial_to, s1$trial), 5), collapse = ","))
  s1 <- rbind(s1[, .(trial, scenario = "REF", arm = "R", sc_T = scenario, col = unname(EPMAP[endpoint]), n = n_R)], s1[, .(trial, scenario, arm = "T", sc_T = scenario, col = unname(EPMAP[endpoint]), n = n_T)])
  s1_ref <- unique(s1[scenario == "REF", .(trial, col, n)])            # 대조군은 시나리오와 무관하게 같아야 한다
  if (anyDuplicated(s1_ref[, .(trial, col)])) stop("section1: 대조군 인원이 시나리오마다 다릅니다")
  s1 <- rbind(s1_ref[, .(trial, scenario = "REF", col, n)], s1[arm == "T", .(trial, scenario, col, n)])
  counts_f <- file.path(out_dir, sprintf("trialpop_counts_%s.csv.gz", model))
  done <- integer(0)
  if (file.exists(counts_f)) { ex <- fread(counts_f); nrow_trial <- (length(scen) + 1L) * 2L; cn <- ex[, .N, by = trial]; done <- cn[N == nrow_trial, trial]
    if (length(setdiff(unique(ex$trial), done))) { tmp <- paste0(counts_f, ".rewrite.csv.gz"); fwrite(ex[trial %in% done], tmp); file.rename(tmp, counts_f) } }
  append_gz <- function(dt, f) { tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f)); fwrite(dt, tmp, col.names = !file.exists(f))
    if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f) }
  lf <- log_or_null(paste0("trialpop_trials_", model), MASTER_SEED, list(model = model, scenarios = paste(names(scen), collapse = ","), trials = sprintf("%d-%d", trial_from, trial_to), cores = cores))
  invisible(get_model(p$model_id)); sd_days <- get_schedule(design, "B0")
  one <- function(j) {
    sa <- simulate_trial_arms(j, p, design, scen, "B0", MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id)
    nmax <- max(sa$subj$id); arms_l <- c(list(REF = sa$arms$R$S00), setNames(lapply(names(scen), function(sc) sa$arms$T[[sc]]), names(scen)))
    obs <- rbindlist(lapply(seq_along(arms_l), function(k) subset_schedule(arms_l[[k]]$obs, sd_days)[, .(id = id + (k - 1L) * nmax, time, conc)]))
    mp <- rbindlist(lapply(seq_along(arms_l), function(k) data.table(id = arms_l[[k]]$subj$id + (k - 1L) * nmax, sid = arms_l[[k]]$subj$id, scenario = names(arms_l)[k])))
    nca <- mp[run_nca(obs), on = "id"]
    nca <- sa$subj[, .(sid = id, stratum = as.integer(stratum))][nca, on = "sid"]
    for (s_ in names(CRIT_SETS)) nca[, (paste0("ok_", s_)) := crit_ok(.SD, s_)]
    stopifnot(identical(nca$ok_ii, nca$reliable %in% TRUE))
    # 무작위배정 인원(층별)은 NCA에 행이 없어도 세도록 대상자 표에서 만든다
    rnd <- rbindlist(lapply(names(arms_l), function(sc) arms_l[[sc]]$subj[, .(n_rand = .N), by = .(stratum = as.integer(stratum))][, scenario := sc]))
    cnt <- nca[, .(n_auclast = sum(is.finite(AUClast) & AUClast > 0), n_lambda = sum(lambda_ok %in% TRUE), n_i = sum(ok_i), n_ii = sum(ok_ii), n_iii = sum(ok_iii), n_iv = sum(ok_iv)), by = .(scenario, stratum)]
    cnt <- cnt[rnd, on = c("scenario", "stratum")]
    for (cc in c("n_auclast", "n_lambda", "n_i", "n_ii", "n_iii", "n_iv")) set(cnt, which(is.na(cnt[[cc]])), cc, 0L)
    cnt[, `:=`(trial = j, arm = fifelse(scenario == "REF", "R", "T"))]
    setcolorder(cnt, c("trial", "scenario", "arm", "stratum", "n_rand"))[order(scenario, stratum)]
  }
  for (b0 in seq(trial_from, trial_to, by = 250L)) {
    ids <- setdiff(b0:min(b0 + 249L, trial_to), done); if (!length(ids)) next
    t0 <- Sys.time()
    res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
    bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.table(x), logical(1))
    if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
    cnt <- rbindlist(res)
    if (!all(cnt[, .N, by = trial]$N == (length(scen) + 1L) * 2L)) stop("시험마다 (시나리오 + 대조군) × 층 2개 행이어야 합니다")
    # section1 대조
    tot <- melt(cnt[scenario %in% c("REF", SEC1), lapply(.SD, sum), by = .(trial, scenario), .SDcols = unname(EPMAP)], id.vars = c("trial", "scenario"), variable.name = "col", value.name = "n_new")
    tot[, col := as.character(col)]
    mm <- merge(tot, s1[trial %in% ids], by = c("trial", "scenario", "col"), all = TRUE)
    if (anyNA(mm$n_new) || anyNA(mm$n) || any(mm$n_new != mm$n)) { print(head(mm[is.na(n_new) | is.na(n) | n_new != n], 10)); stop("section1 n_R·n_T와 다릅니다. 결과를 쓰지 않습니다.") }
    append_gz(cnt, counts_f)
    say(lf, sprintf("%s trials %d-%d done in %s; per-arm counts equal section1 for %d arm x endpoint rows", model, min(ids), max(ids), format(Sys.time() - t0), nrow(mm)))
  }
  say(lf, sprintf("done: %s", counts_f))
}
