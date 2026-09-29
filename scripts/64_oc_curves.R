#!/usr/bin/env Rscript
# 64_oc_curves.R — S9 재설계 운용특성 곡선용 시험 재생성 (지시 2026-09-29 "S9 재설계: 1종·2종 오류", config/prereg_20260929_oc.yaml section8, D-064).
# 경계(0.80, 1.25)와 동일 제품(1.00)이 아닌 칸(기전 F, ke, Vmax, V2, ka × 목표 0.85, 0.90, 0.95, 1.05, 1.11, 1.18 중 도달 가능한 칸)을
# scripts/44와 같은 시드 규칙·역산 배율·함수(R/oc_models.R run_trial_oc_models)로 대상자 수준에서 다시 만들고 M0·M1·M2로 판정한다.
#  시험 1–2,000: 모든 칸. 시험 2,001–5,000: 0.95·1.05 칸만(시나리오 목록이 시험 번호에 따라 줄어든다).
# 검사(쓰기 전, 묶음마다; 다르면 중단):
#  - 시험 1–2,000의 M0를 results/oc/oc_trials_be_<model>.csv.gz(평가변수 6개)와 대조(GMR·CI 상대 차이 ≤ 1e-6, pass·n_R·n_T 동일)
#  - 시작 전 한 번: 시험 2,000을 줄인 목록(0.95·1.05 칸)으로 다시 만들어 전체 목록 결과와 모든 값이 같은지(목록 구성과 무관함 확인)
# 산출(section8 output_dir): oc_curves_be_<model>.csv.gz(OC_ENDPOINTS_EXT × M0·M1·M2), oc_curves_crit_be_<model>.csv.gz(OC_ENDPOINTS_CRIT × M0·M1)
#   열은 scripts/44와 같다: trial, scenario, endpoint, model, est(log GMR, signif 8), se(signif 8), pass, n_R, n_T. 500회 묶음 체크포인트, 재시작 시 완료 시험 건너뜀.
# 사용법: nice -n 10 Rscript scripts/64_oc_curves.R <k2016|k2020> [cores=2] [trial_from=1] [trial_to=5000] [out_dir]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("k2016", "k2020")) stop("사용법: Rscript scripts/64_oc_curves.R <k2016|k2020> [cores] [trial_from] [trial_to] [out_dir]")
model <- args[1]
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml"); s8 <- read_cfg("prereg_20260929_oc.yaml")$section8
cores <- if (length(args) >= 2) as.integer(args[2]) else 2L
n_near <- as.integer(s8$reps$near_one$trials); n_other <- as.integer(s8$reps$other$trials)
trial_from <- if (length(args) >= 3) as.integer(args[3]) else 1L
trial_to <- if (length(args) >= 4) as.integer(args[4]) else n_near
stopifnot(!is.na(cores), cores >= 1, model %in% unlist(s8$pk_models), n_near >= n_other, trial_from >= 1, trial_to <= n_near, trial_to >= trial_from)
default_dir <- proj_path(s8$output_dir)
out_dir <- if (length(args) >= 5) args[5] else default_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
production <- normalizePath(out_dir) == normalizePath(default_dir, mustWork = FALSE)
oc_dir <- proj_path("results", "oc")
batch <- as.integer(oc$trials$batch); MASTER_SEED <- as.integer(oc$trials$master_seed)
stopifnot(MASTER_SEED == 20260923L)
MODELS <- c("M0", "M1", "M2"); MODELS_CRIT <- c("M0", "M1")
EPS <- OC_ENDPOINTS_EXT; EPS_CRIT <- OC_ENDPOINTS_CRIT
COLS <- c("trial", "scenario", "endpoint", "model", "est", "se", "pass", "n_R", "n_T")
variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p

# 칸: 역산 결과의 도달 행(scripts/44와 같은 코드 규칙), 기전·목표는 section8
inv <- rbindlist(lapply(list.files(oc_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", model), full.names = TRUE), fread))
stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
inv <- inv[reachable == TRUE]
inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
s0 <- fread(file.path(oc_dir, sprintf("oc_scenarios_%s.csv", model)))
mm <- merge(inv[, .(code, multiplier)], s0[, .(code, m0 = multiplier)], by = "code", all = TRUE)
if (nrow(mm) != nrow(inv) || anyNA(mm) || any(abs(mm$multiplier / mm$m0 - 1) > 1e-12)) stop("oc_scenarios와 역산 결과의 배율이 다릅니다")
near_t <- as.numeric(unlist(s8$reps$near_one$targets)); other_t <- as.numeric(unlist(s8$reps$other$targets))
on_t <- function(x, tt) vapply(x, function(v) any(abs(v - tt) < 1e-9), logical(1))
cells <- inv[mechanism %in% unlist(s8$cells$mechanisms) & on_t(target, c(near_t, other_t))]
if (anyDuplicated(cells[, .(mechanism, target)])) stop("한 기전·목표에 도달 방향이 둘입니다")
cells[, reps := fifelse(on_t(target, near_t), n_near, n_other)]
setorder(cells, target, mechanism)
mk <- function(cd) { r <- cells[code == cd]; list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }
scen_all <- setNames(lapply(cells$code, mk), cells$code)
scen_for <- function(j) scen_all[cells[reps >= j, code]]
fwrite(cells[, .(pk_model = model, code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio, reps)], file.path(out_dir, sprintf("oc_curves_cells_%s.csv", model)))

# 저장본(M0, 시험 1–2,000)
stored <- fread(file.path(oc_dir, sprintf("oc_trials_be_%s.csv.gz", model)))[scenario %in% cells$code & trial <= n_other]
if (!setequal(unique(stored$scenario), cells$code) || nrow(stored) != n_other * nrow(cells) * length(OC_ENDPOINTS))
  stop("oc_trials_be 저장본이 모든 칸의 시험 1–2,000 × 6개 평가변수를 갖고 있지 않습니다")

rel_ok <- function(a, b, tol) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a / b - 1) <= tol)
same_pass <- function(a, b) (is.na(a) & is.na(b)) | (a %in% TRUE & b %in% TRUE) | (a %in% FALSE & b %in% FALSE)
m0_ci <- function(x) {
  q <- qt(1 - (1 - design$be$ci_level) / 2, x$n_R + x$n_T - 2)
  x[, .(trial, scenario, endpoint, GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se), pass, n_R, n_T)]
}
verify <- function(new, ref, tol, what) {
  m <- merge(new, ref, by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
  if (nrow(m) != nrow(new) || nrow(m) != nrow(ref)) stop(sprintf("%s: 행 수 불일치(재생성 %d, 저장본 %d, 병합 %d)", what, nrow(new), nrow(ref), nrow(m)))
  ok <- (rel_ok(m$GMR, m$GMR.s, tol) & rel_ok(m$CI_lower, m$CI_lower.s, tol) & rel_ok(m$CI_upper, m$CI_upper.s, tol) & same_pass(m$pass, m$pass.s) &
           m$n_R == m$n_R.s & m$n_T == m$n_T.s) %in% TRUE
  if (!all(ok)) { print(head(m[!ok], 10)); stop(sprintf("%s: %d행이 저장본과 다릅니다(시험 %s). 결과를 쓰지 않습니다.", what, sum(!ok), paste(head(unique(m$trial[!ok]), 5), collapse = ","))) }
  nrow(m)
}
append_gz <- function(dt, f) {
  tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f))
  fwrite(dt, tmp, col.names = !file.exists(f))
  if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f)
}
rows_per_trial <- function(j, n_ep, n_mod) length(scen_for(j)) * n_ep * n_mod
done_trials <- function(f, n_ep, n_mod) {                            # 완전한 시험만 남기고 불완전 행은 지운다
  if (!file.exists(f)) return(integer(0))
  ex <- fread(f)
  if (!identical(names(ex), COLS)) stop("기존 파일의 열이 다릅니다: ", f)
  cn <- ex[, .N, by = trial]; cn[, need := vapply(trial, rows_per_trial, 1, n_ep = n_ep, n_mod = n_mod)]
  bad <- cn[N != need]
  if (nrow(bad) || anyDuplicated(ex, by = c("trial", "scenario", "endpoint", "model"))) {
    message(sprintf("기존 파일의 불완전 시험 %d개를 지우고 다시 씁니다: %s", nrow(bad), f))
    ex <- unique(ex[!trial %in% bad$trial], by = c("trial", "scenario", "endpoint", "model"))
    tmp <- paste0(f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, f)
  }
  setdiff(unique(ex$trial), bad$trial)
}
finish_rows <- function(be) {
  stopifnot(identical(unique(be$lloq), study_lloq()), "resid" %in% names(be) && identical(unique(be$resid), "fixed"))
  be[, c("lloq", "resid") := NULL]
  be[, `:=`(est = signif(est, 8), se = signif(se, 8))]
  be[, ..COLS]
}
split_rows <- function(be) list(main = finish_rows(be[endpoint %in% EPS]), crit = finish_rows(be[endpoint %in% EPS_CRIT & model %in% MODELS_CRIT]))
one <- function(j) run_trial_oc_models(j, p, design, scen_for(j), MASTER_SEED, rv$wt_spec, model_id = p$model_id, models = MODELS, endpoints = c(EPS, EPS_CRIT))$be
invisible(get_model(p$model_id))                                     # fork 전에 모델 컴파일

out_f <- file.path(out_dir, sprintf("oc_curves_be_%s.csv.gz", model)); crit_f <- file.path(out_dir, sprintf("oc_curves_crit_be_%s.csv.gz", model))
done <- intersect(done_trials(out_f, length(EPS), length(MODELS)), done_trials(crit_f, length(EPS_CRIT), length(MODELS_CRIT)))
for (f in c(out_f, crit_f)) if (file.exists(f)) { ex <- fread(f); if (any(!ex$trial %in% done)) { tmp <- paste0(f, ".rewrite.csv.gz"); fwrite(ex[trial %in% done], tmp); file.rename(tmp, f) } }
logfile <- if (production) start_run_log(paste0("oc_curves_", model), master_seed = MASTER_SEED, run_mode = "final",
                                          extra = list(model = model, cells = nrow(cells), cells_near_one = sum(cells$reps == n_near), trial_from = trial_from, trial_to = trial_to, cores = cores)) else NULL
say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }
say(sprintf("%s: %d cells (%d with %d trials, %d with %d trials), trials %d-%d, cores %d, output %s", model, nrow(cells), sum(cells$reps == n_near), n_near,
            sum(cells$reps == n_other), n_other, trial_from, trial_to, cores, out_f))

# 목록 구성과 무관함: 시험 n_other를 전체 목록과 줄인 목록으로 만들어 줄인 목록의 칸이 모두 같은지
if (trial_to > n_other) {
  full <- run_trial_oc_models(n_other, p, design, scen_all, MASTER_SEED, rv$wt_spec, model_id = p$model_id, models = MODELS, endpoints = c(EPS, EPS_CRIT))$be
  red <- run_trial_oc_models(n_other, p, design, scen_for(n_other + 1L), MASTER_SEED, rv$wt_spec, model_id = p$model_id, models = MODELS, endpoints = c(EPS, EPS_CRIT))$be
  a <- full[scenario %in% names(scen_for(n_other + 1L))]; setkey(a, scenario, endpoint, model); setkey(red, scenario, endpoint, model)
  cmp <- c("est", "se", "pass", "n_R", "n_T")
  if (nrow(a) != nrow(red) || !isTRUE(all.equal(a[, ..cmp], red[, ..cmp], tolerance = 0, check.attributes = FALSE)))
    stop("시나리오 목록을 줄이면 결과가 달라집니다(시험 ", n_other, ")")
  say(sprintf("check: trial %d gives identical results for the 0.95/1.05 cells with the full and the reduced scenario list (%d rows)", n_other, nrow(red)))
}

n_chk <- 0L
for (b0 in seq(trial_from, trial_to, by = batch)) {
  ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done)
  if (!length(ids)) next
  t0 <- Sys.time()
  res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
  bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.table(x), logical(1))
  if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
  be <- rbindlist(res)
  need <- sum(vapply(ids, function(j) length(scen_for(j)), 1L)) * length(c(EPS, EPS_CRIT)) * length(MODELS)
  stopifnot(nrow(be) == need)
  st_ids <- ids[ids <= n_other]
  if (length(st_ids)) n_chk <- n_chk + verify(m0_ci(be[model == "M0" & endpoint %in% OC_ENDPOINTS & trial %in% st_ids]), stored[trial %in% st_ids], 1e-6, sprintf("%s oc_trials_be 대조", model))
  sp <- split_rows(be); append_gz(sp$crit, crit_f); append_gz(sp$main, out_f)
  say(sprintf("trials %d-%d done in %s (%d scenario-trials)%s", min(ids), max(ids), format(Sys.time() - t0), sum(vapply(ids, function(j) length(scen_for(j)), 1L)),
              if (length(st_ids)) sprintf("; M0 matched the stored oc_trials_be rows for trials %d-%d", min(st_ids), max(st_ids)) else ""))
}
say(sprintf("done: M0 matched %d stored oc_trials_be rows; output %s", n_chk, out_f))
