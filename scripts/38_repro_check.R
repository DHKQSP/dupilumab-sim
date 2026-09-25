#!/usr/bin/env Rscript
# 재현성 점검 (config/repro_check.yaml, .github/workflows/repro.yml, README "CI"). 허용 오차는 config에 사전 고정.
# (a) 운용특성 시험 정확 재생성: k2016 시험 1-100, S00·Vmax_down_125를 scripts/31_oc_trials.R와 같은 방식으로 다시 만들어 커밋된 행과 비교
# (b) 참 AUC0-inf 비: Vmax 역산 배율(약 x0.6602)에서 공통 난수 대상자 앞 20,000명 대 커밋된 200,000명 값(MC 표준오차 기반 허용)
# (c) 개인 수준 참 외삽 비율 중앙값(B0): 새 시드 2,000명 대 커밋된 20,000명 값(부트스트랩 표준오차 기반 허용)
# 출력: results/repro/repro_check.csv(item, committed, reproduced, tolerance, pass, ...), results/repro/repro_check_meta.csv.
# 하나라도 불통과면 종료 코드 1. GitHub Actions에서는 작업 요약(GITHUB_STEP_SUMMARY)도 쓴다.
# 사용법: Rscript scripts/38_repro_check.R [코어]
#         Rscript scripts/38_repro_check.R --extract-reference   — 전체 OC 시험 파일에서 (a) 비교 행 발췌본 작성(전체 파일이 커밋 제외일 때 대비)
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
rc <- read_cfg("repro_check.yaml"); oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", "repro"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
A <- rc$oc_trials; B <- rc$truth_ratio; C <- rc$individual_extrap
trial_ids <- seq(as.integer(A$trial_range[[1]]), as.integer(A$trial_range[[2]]))
scen_codes <- as.character(unlist(A$scenarios))
ref_cols <- c("trial", "scenario", "endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T")
full_f <- proj_path(A$reference_file); extract_f <- proj_path(A$reference_extract)

# 커밋된 (a) 비교 행: 전체 파일에서 시험·시나리오 부분집합(중복 행이 있으면 중단)
read_ref_rows <- function(path_) {
  x <- fread(path_, select = ref_cols)
  x <- x[trial %in% trial_ids & scenario %in% scen_codes]
  if (anyDuplicated(x, by = c("trial", "scenario", "endpoint"))) stop("참조 파일에 중복 행: ", path_)
  setorder(x, trial, scenario, endpoint)[]
}

if (length(args) && args[1] == "--extract-reference") {
  if (!file.exists(full_f)) stop("전체 OC 시험 파일이 없습니다: ", full_f)
  x <- read_ref_rows(full_f)
  dir.create(dirname(extract_f), showWarnings = FALSE, recursive = TRUE)
  fwrite(x, extract_f)
  cat(sprintf("발췌본 %d행 → %s\n", nrow(x), extract_f)); quit(save = "no")
}
cores <- if (length(args)) as.integer(args[1]) else as.integer(A$cores)
logfile <- start_run_log("repro_check", master_seed = as.integer(oc$trials$master_seed), run_mode = "final",
                         extra = list(trials = range(trial_ids), scenarios = scen_codes, cores = cores))
t_start <- Sys.time(); items <- list(); meta_x <- list()
add_item <- function(item, committed, reproduced, tolerance, tolerance_type, pass_, rule, detail = "") {
  items[[length(items) + 1]] <<- data.table(item = item, committed = committed, reproduced = reproduced, tolerance = tolerance, pass = isTRUE(pass_),
                                            tolerance_type = tolerance_type, rule = rule, detail = detail)
  cat(sprintf("[%s] %-42s committed %-14s reproduced %-14s tol %-9s %s\n", format(Sys.time(), "%H:%M:%S"), item, format(committed, digits = 8),
              format(reproduced, digits = 8), format(tolerance, digits = 3), if (isTRUE(pass_)) "PASS" else "FAIL"))
  append_run_log(logfile, sprintf("%s: committed %s, reproduced %s, tolerance %s (%s), pass %s", item, format(committed, digits = 15),
                                   format(reproduced, digits = 15), format(tolerance, digits = 6), tolerance_type, isTRUE(pass_)))
}

# ---- (a) 운용특성 시험 정확 재생성 --------------------------------------------------------------
t0 <- Sys.time()
if (file.exists(full_f)) {
  ref_a <- read_ref_rows(full_f); ref_src <- A$reference_file
  if (file.exists(extract_f)) {                                   # 발췌본이 전체 파일과 같은지 확인
    ex <- setorder(fread(extract_f), trial, scenario, endpoint)
    if (!isTRUE(all.equal(as.data.frame(ex), as.data.frame(ref_a), tolerance = 0, check.attributes = FALSE))) stop("발췌본이 전체 파일과 다릅니다: ", extract_f)
  }
} else if (file.exists(extract_f)) {
  ref_a <- setorder(fread(extract_f), trial, scenario, endpoint); ref_src <- A$reference_extract
} else stop("(a) 참조 파일이 없습니다: ", A$reference_file, ", ", A$reference_extract)
expected_rows <- length(trial_ids) * length(scen_codes) * length(OC_ENDPOINTS)
if (nrow(ref_a) != expected_rows) stop(sprintf("(a) 참조 행 수 %d, 기대 %d", nrow(ref_a), expected_rows))
# scripts/31_oc_trials.R와 같은 구성: 변형, 역산 파일의 도달 시나리오 목록(코드 규칙), 시드
model_a <- A$model
rv <- resolve_variant(A$variant, design); p_a <- rv$p
inv <- rbindlist(lapply(list.files(proj_path("results", "oc"), pattern = sprintf("^inversion_%s_.*\\.csv$", model_a), full.names = TRUE), fread))
stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
inv <- inv[reachable == TRUE]
inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
scen_all <- c(list(S00 = list(code = "S00", T_multipliers = list())),
              setNames(lapply(seq_len(nrow(inv)), function(i) list(code = inv$code[i], T_multipliers = setNames(list(inv$multiplier[i]), inv$mechanism[i]))), inv$code))
if (!all(scen_codes %in% names(scen_all))) stop("(a) 시나리오가 역산 목록에 없습니다: ", paste(setdiff(scen_codes, names(scen_all)), collapse = ","))
# 시험 j의 난수(대상자·배정·채혈시각·잔차)는 시나리오와 무관하고 시험군 시나리오는 개체별로 독립 풀이되므로(R/trial.R), 필요한 두 시나리오만 푼다
rr <- run_trials_oc(trial_ids, p_a, design, scen_all[scen_codes], as.integer(oc$trials$master_seed), rv$wt_spec, cores = cores, model_id = p_a$model_id, progress_every = 0)
rep_a <- setorder(rr$be[, ..ref_cols], trial, scenario, endpoint)
m <- merge(ref_a, rep_a, by = c("trial", "scenario", "endpoint"), suffixes = c(".c", ".r"), all = TRUE)
rel <- function(a, b) fifelse(is.na(a) & is.na(b), 0, fifelse(is.na(a) | is.na(b), Inf, abs(b / a - 1)))
cmp_cols <- as.character(unlist(A$compare_cols))
max_rel <- max(vapply(cmp_cols, function(cc) max(rel(m[[paste0(cc, ".c")]], m[[paste0(cc, ".r")]])), numeric(1)))
flag_mis <- sum(xor(is.na(m$pass.c), is.na(m$pass.r)) | ((m$pass.c %in% TRUE) != (m$pass.r %in% TRUE)) |
                  (m$n_R.c != m$n_R.r) %in% TRUE | (m$n_T.c != m$n_T.r) %in% TRUE)
rule_a <- sprintf("max |reproduced/committed - 1| over %d rows x (%s); trials %d-%d, %s, run_trial_oc as scripts/31_oc_trials.R",
                  nrow(m), paste(cmp_cols, collapse = ", "), min(trial_ids), max(trial_ids), paste(scen_codes, collapse = " + "))
add_item("a1_oc_trials_GMR_CI_max_rel_diff", 0, max_rel, as.numeric(A$tolerance_rel), "relative", max_rel <= as.numeric(A$tolerance_rel) && nrow(m) == expected_rows,
         rule_a, sprintf("committed column = ideal difference 0; reference %s", ref_src))
add_item("a2_oc_trials_pass_flag_mismatches", 0, flag_mis, as.numeric(A$pass_flag_mismatch_max), "count", flag_mis <= as.numeric(A$pass_flag_mismatch_max),
         "rows whose pass flag, NA status, n_R or n_T differ", sprintf("%d rows compared", nrow(m)))
cfg <- A$configuration
p2_c <- config_pass(ref_a, oc); p2_r <- config_pass(rep_a, oc)
cfg_col <- paste0("cfg_", cfg)
for (s_ in scen_codes) {
  rc_ <- 100 * mean(p2_c[scenario == s_][[cfg_col]]); rr_ <- 100 * mean(p2_r[scenario == s_][[cfg_col]])
  xc <- sum(p2_c[scenario == s_][[cfg_col]]); xr <- sum(p2_r[scenario == s_][[cfg_col]]); nn <- nrow(p2_r[scenario == s_])
  w <- wilson_ci(xr, nn)
  lab <- if (s_ == "S00") "identical-product pass rate" else "boundary type I error (true AUC0-inf ratio 1.25, subset)"
  add_item(sprintf("a3_%s_pass_rate_%s_pct", cfg, s_), rc_, rr_, 0, "absolute (pct points)", abs(rr_ - rc_) <= 0,
           sprintf("%s: %s pass rate over trials %d-%d", lab, cfg, min(trial_ids), max(trial_ids)),
           sprintf("committed %d/%d, reproduced %d/%d, Wilson 95%% CI %.1f-%.1f%%", xc, nn, xr, nn, w$lo, w$hi))
}
meta_x$runtime_a_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1); meta_x$reference_source_a <- ref_src

# ---- (b) 참 AUC0-inf 비, 앞 20,000명 ----------------------------------------------------------
t0 <- Sys.time()
p_b <- load_params(B$model); dose <- as.numeric(oc$estimand$dose_mg)
inv_b <- fread(proj_path("results", "oc", sprintf("inversion_%s_%s.csv", B$model, B$mechanism)))
row_b <- inv_b[direction == B$direction & abs(target - as.numeric(B$target)) < 1e-9 & reachable == TRUE]
stopifnot(nrow(row_b) == 1)
N_all <- as.integer(oc$estimand$population$n_subjects); N_sub <- as.integer(B$n_subjects)
seed_b <- derive_seed(as.integer(oc$truth_seed), "truth200k", B$model)          # scripts/30_oc_inversion.R와 같은 시드 규칙
subj_b <- truth_subjects(N_all, p_b, design, seed_b)[id <= N_sub]
ref_b <- truth_metrics(individual_params(p_b, subj_b), dose, p_b$model_id, cmax = TRUE)
te_b <- truth_metrics(ip_for_mult(p_b, subj_b, mech_multiplier(B$mechanism, row_b$multiplier)), dose, p_b$model_id, cmax = FALSE)
tr_b <- truth_ratio(ref_b[, .(id, AUCinf)], te_b)
tol_b <- as.numeric(B$z) * row_b$auc_se_log * sqrt(N_all / N_sub - 1)
mult_lab <- sprintf("%s x%.4f", B$mechanism, row_b$multiplier)
add_item(sprintf("b1_true_AUCinf_ratio_%s_%s_first%d", B$mechanism, sprintf("%.4f", row_b$multiplier), N_sub), row_b$auc_ratio, tr_b$auc_ratio, tol_b, "absolute log ratio",
         abs(log(tr_b$auc_ratio / row_b$auc_ratio)) <= tol_b,
         sprintf("%s (multiplier %.15g): first %d of %d CRN subjects vs committed %d-subject value; tolerance = z %g x se200k %.4g x sqrt(%d/%d - 1)",
                 mult_lab, row_b$multiplier, N_sub, N_all, row_b$n_subjects, as.numeric(B$z), row_b$auc_se_log, N_all, N_sub),
         sprintf("reproduced MC SE (log) %.3g; |log difference| %.3g", tr_b$auc_se_log, abs(log(tr_b$auc_ratio / row_b$auc_ratio))))
if (isTRUE(B$screening_check) && row_b$iter_full %in% 0L) {
  tol_s <- as.numeric(oc$inversion$screening_tolerance_rel)
  add_item(sprintf("b2_screening_consistency_first%d", N_sub), as.numeric(B$target), tr_b$auc_ratio, tol_s, "relative",
           abs(tr_b$auc_ratio / as.numeric(B$target) - 1) <= tol_s,
           sprintf("inversion screening bisection used the same first %d subjects and stopped within +-%g of the target (iter_full = 0)", N_sub, tol_s),
           sprintf("relative difference %.3g", abs(tr_b$auc_ratio / as.numeric(B$target) - 1)))
}
meta_x$runtime_b_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

# ---- (c) 개인 수준 참 외삽 비율 중앙값, B0, 2,000명 ------------------------------------------------
t0 <- Sys.time()
rv_c <- resolve_variant(C$variant, design, load_scenarios()); p_c <- rv_c$p
pop <- run_individual_population(as.integer(C$n_subjects), p_c, design, design$schedule_analysis, as.integer(C$master_seed), rv_c$wt_spec,
                                 jitter = TRUE, model_id = p_c$model_id, tag = C$tag)                   # scripts/10_individual_schedules.R와 같은 호출(태그만 다름)
sched_c <- C$schedule
sm <- summarize_individual(pop$nca[schedule == sched_c], last_planned = max(get_schedule(design, sched_c)))
committed_c <- fread(proj_path(C$committed_file))[schedule == sched_c][[C$column]]
stopifnot(length(committed_c) == 1)
add_item(sprintf("c1_%s_%s_%dsubjects", C$column, sched_c, as.integer(C$n_subjects)), committed_c, sm[[C$column]], as.numeric(C$tolerance_abs), "absolute (pct)",
         abs(sm[[C$column]] - committed_c) <= as.numeric(C$tolerance_abs),
         sprintf("median true extrapolated share at %s, %d new subjects (seed %s/%s) vs committed 20,000 (%s); tolerance = 4 x bootstrap SE",
                 sched_c, as.integer(C$n_subjects), C$master_seed, C$tag, C$committed_file),
         sprintf("|difference| %.3g", abs(sm[[C$column]] - committed_c)))
meta_x$runtime_c_s <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

# ---- 출력 --------------------------------------------------------------------------------------
res <- rbindlist(items)
fwrite(res, file.path(out_dir, "repro_check.csv"))
git1 <- function(...) tryCatch(suppressWarnings(system2("git", c("-C", PROJ_ROOT, ...), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
sha <- Sys.getenv("GITHUB_SHA"); sha_src <- "GITHUB_SHA"
if (!nzchar(sha)) { sha <- paste(git1("rev-parse", "HEAD"), collapse = ""); sha_src <- "git rev-parse HEAD" }
dirty <- git1("status", "--porcelain", "--untracked-files=no")
pv <- function(pkg) tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
meta <- c(list(commit = if (nzchar(sha)) sha else NA_character_, commit_source = sha_src,
               tracked_changes_uncommitted = length(dirty), r_version = R.version$version.string,
               rxode2_version = pv("rxode2"), data_table_version = pv("data.table"), testthat_version = pv("testthat"),
               platform = R.version$platform, os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
               cores_detected = parallel::detectCores(), cores_used = cores,
               github_run = if (nzchar(Sys.getenv("GITHUB_RUN_ID"))) sprintf("%s/%s/actions/runs/%s", Sys.getenv("GITHUB_SERVER_URL"), Sys.getenv("GITHUB_REPOSITORY"), Sys.getenv("GITHUB_RUN_ID")) else NA_character_,
               time_start_utc = format(t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), time_end_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
               runtime_total_s = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
               repro_config_sha256 = digest::digest(file = proj_path("config", "repro_check.yaml"), algo = "sha256"),
               n_items = nrow(res), n_pass = sum(res$pass), all_pass = all(res$pass)), meta_x)
fwrite(data.table(field = names(meta), value = unname(vapply(meta, function(v) paste(format(v), collapse = " "), ""))), file.path(out_dir, "repro_check_meta.csv"))
print(res[, .(item, committed, reproduced, tolerance, pass)], digits = 8)
gs <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(gs)) {
  f8 <- function(v) formatC(v, digits = 8, format = "g")
  cat(c(sprintf("## Reproducibility check: %s (%d/%d items pass)", if (all(res$pass)) "PASS" else "FAIL", sum(res$pass), nrow(res)), "",
        sprintf("Commit %s, %s, rxode2 %s, runtime %s s.", substr(meta$commit, 1, 12), R.version$version.string, meta$rxode2_version, meta$runtime_total_s), "",
        "| item | committed | reproduced | tolerance | type | pass |", "|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s |", res$item, f8(res$committed), f8(res$reproduced), f8(res$tolerance), res$tolerance_type, ifelse(res$pass, "yes", "**NO**")), ""),
      file = gs, sep = "\n", append = TRUE)
}
append_run_log(logfile, sprintf("done: %d/%d pass", sum(res$pass), nrow(res)))
if (!all(res$pass)) quit(status = 1)
