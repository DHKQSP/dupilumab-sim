#!/usr/bin/env Rscript
# 40_oc_rejudge.R — 경계 시나리오 재생성: AUCinf 규칙 × 플래그 세트 재판정용 (검토 의견 W2 §2, D-040).
# scripts/31_oc_trials.R와 같은 시드·같은 역산 배율로 경계(0.80, 1.25) 시나리오 + S00만 다시 모의하고 run_trial_oc_ext로
# 평가변수 8개(원래 6개 + AUCinf_Ai, AUCinf_Ci = 플래그 세트 (i))를 판정한다.
#  - 31은 시험 1–reps_other를 전체 시나리오 목록으로, 나머지를 경계 목록으로 돌렸다. 여기서는 모든 시험을 경계 목록으로 돌린다
#    (시나리오 목록이 한 시나리오의 결과를 바꾸지 않는다는 가정: 시험군 T는 시나리오별 id 오프셋으로 한 번에 풀고 대상자·채혈·잔차 시드는
#    시나리오와 무관. 아래 저장본 대조가 시험마다 이 가정을 검사한다).
#  - 검사: 재생성한 모든 시험에서 원래 6개 평가변수가 저장본 results/oc/oc_trials_be_<model>.csv.gz와 같아야 한다
#    (GMR·CI 상대 차이 ≤ 1e-6: 저장본은 signif 7; pass·n_R·n_T 동일). 하나라도 다르면 쓰기 전에 중단.
#  - 산출: <out_dir>/oc_rejudge_be_<model>.csv.gz (trial, scenario, endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T; 8개 평가변수, signif 7).
#    500회 묶음 체크포인트: 묶음마다 임시 gz를 만든 뒤 file.append로 붙인다(재시작 시 완료 시험 건너뜀). 같은 출력 파일에 두 실행을 동시에 돌리지 않는다.
# 사용법: nice -n 15 Rscript scripts/40_oc_rejudge.R <k2016|k2020> [cores=3] [trial_from=1] [trial_to=reps_boundary] [out_dir=results/oc]
#   out_dir는 환경변수 DUPI_OC_REJUDGE_DIR로도 지정(인수가 우선). 시험 실행은 임시 디렉터리로(results/oc/에 부분 재판정 파일을 남기지 않는다).
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("k2016", "k2020")) stop("사용법: Rscript scripts/40_oc_rejudge.R <k2016|k2020> [cores] [trial_from] [trial_to] [out_dir]")
model <- args[1]
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
n_bnd <- as.integer(oc$trials$reps_boundary); batch <- as.integer(oc$trials$batch)
trial_from <- if (length(args) >= 3) as.integer(args[3]) else 1L
trial_to <- if (length(args) >= 4) as.integer(args[4]) else n_bnd
stopifnot(!is.na(cores), cores >= 1, !is.na(trial_from), !is.na(trial_to), trial_from >= 1, trial_to >= trial_from, trial_to <= n_bnd)
default_dir <- proj_path("results", "oc")
out_dir <- if (length(args) >= 5) args[5] else Sys.getenv("DUPI_OC_REJUDGE_DIR", default_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
production <- normalizePath(out_dir) == normalizePath(default_dir)

variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p
MASTER_SEED <- as.integer(oc$trials$master_seed)

# 시나리오: 31과 같은 규칙(역산 결과의 도달 행, 코드 = 기전_방향_목표×100)
inv <- rbindlist(lapply(list.files(default_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", model), full.names = TRUE), fread))
stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
inv <- inv[reachable == TRUE]
inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
sc_f <- file.path(default_dir, sprintf("oc_scenarios_%s.csv", model))
if (file.exists(sc_f)) {                                              # 31이 끝나 시나리오 표가 있으면 같은 배율인지 확인
  s0 <- fread(sc_f); m <- merge(inv[, .(code, multiplier)], s0[, .(code, m0 = multiplier)], by = "code", all = TRUE)
  if (nrow(m) != nrow(inv) || any(is.na(m$multiplier)) || any(is.na(m$m0)) || any(abs(m$multiplier / m$m0 - 1) > 1e-12))
    stop(sprintf("%s와 역산 결과의 시나리오·배율이 다릅니다", sc_f))
}
bnd <- as.numeric(unlist(oc$boundary_targets))
bcodes <- inv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9, code]
scen_bnd <- c(list(S00 = list(code = "S00", T_multipliers = list())),
              setNames(lapply(bcodes, function(cd) { r <- inv[code == cd]; list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }), bcodes))
SCN <- names(scen_bnd); n_rows_trial <- length(SCN) * length(OC_ENDPOINTS_EXT)
cat(sprintf("%s: 경계 시나리오 %d개 + S00, 시험 %d–%d, cores %d, 출력 %s\n", model, length(bcodes), trial_from, trial_to, cores, out_dir))

# 저장본(원래 6개 평가변수): 요청 범위의 모든 시험 × 시나리오 × 평가변수가 있어야 한다(31이 아직 돌고 있으면 그 범위는 나중에)
raw_f <- file.path(default_dir, sprintf("oc_trials_be_%s.csv.gz", model))
if (!file.exists(raw_f)) stop("저장본 없음: ", raw_f)
stored <- fread(raw_f)[trial >= trial_from & trial <= trial_to & scenario %in% SCN]
cnt <- stored[, .N, by = trial]
need <- length(SCN) * length(OC_ENDPOINTS)
miss_trials <- setdiff(trial_from:trial_to, cnt[N == need, trial])
if (length(miss_trials))
  stop(sprintf("저장본 %s에 시험 %d개(예: %s)의 경계·S00 행이 완전하지 않습니다(시험당 %d행 필요). 31이 이 범위를 마친 뒤 다시 실행하세요.",
               basename(raw_f), length(miss_trials), paste(head(miss_trials, 5), collapse = ","), need))
if (anyDuplicated(stored, by = c("trial", "scenario", "endpoint"))) stop("저장본에 (trial, scenario, endpoint) 중복")

# 재생성 결과의 원래 6개 평가변수를 저장본과 대조. 불일치 1건이라도 있으면 중단
rel_ok <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a / b - 1) <= 1e-6)
verify_batch <- function(new, ids) {
  n6 <- new[endpoint %in% OC_ENDPOINTS]; ref <- stored[trial %in% ids]
  m <- merge(n6, ref, by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
  if (nrow(m) != nrow(n6) || nrow(m) != nrow(ref)) stop(sprintf("저장본 대조: 행 수 불일치(재생성 %d, 저장본 %d, 병합 %d)", nrow(n6), nrow(ref), nrow(m)))
  ok <- rel_ok(m$GMR, m$GMR.s) & rel_ok(m$CI_lower, m$CI_lower.s) & rel_ok(m$CI_upper, m$CI_upper.s) &
    ((is.na(m$pass) & is.na(m$pass.s)) | (m$pass %in% TRUE & m$pass.s %in% TRUE) | (m$pass %in% FALSE & m$pass.s %in% FALSE)) &
    m$n_R == m$n_R.s & m$n_T == m$n_T.s
  ok <- ok %in% TRUE
  if (!all(ok)) {
    bad <- m[!ok]
    print(head(bad, 10))
    stop(sprintf("저장본 대조 실패: %d행(시험 %s; 시나리오 %s)에서 원래 평가변수가 oc_trials_be_%s.csv.gz와 다릅니다. 재판정 결과를 쓰지 않습니다.",
                 nrow(bad), paste(head(unique(bad$trial), 10), collapse = ","), paste(head(unique(bad$scenario), 5), collapse = ","), model))
  }
  invisible(nrow(m))
}

out_f <- file.path(out_dir, sprintf("oc_rejudge_be_%s.csv.gz", model))
done <- integer(0)
if (file.exists(out_f)) {
  ex <- fread(out_f)
  if (!identical(names(ex), c("trial", "scenario", "endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T"))) stop("기존 재판정 파일의 열이 다릅니다: ", out_f)
  cn <- ex[, .N, by = trial]
  if (any(cn$N != n_rows_trial) || anyDuplicated(ex, by = c("trial", "scenario", "endpoint"))) {
    bad <- cn[N != n_rows_trial, trial]
    message(sprintf("기존 재판정 파일의 불완전 시험 %d개(예: %s)를 지우고 다시 씁니다", length(bad), paste(head(bad, 5), collapse = ",")))
    ex <- unique(ex[!trial %in% bad], by = c("trial", "scenario", "endpoint"))
    tmp <- paste0(out_f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, out_f)
  }
  done <- unique(ex$trial); rm(ex)
}
logfile <- if (production) start_run_log(paste0("oc_rejudge_", model), master_seed = MASTER_SEED, run_mode = "final",
                                          extra = list(model = model, n_scen = length(SCN), trial_from = trial_from, trial_to = trial_to, cores = cores)) else NULL
say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }

invisible(get_model(p$model_id))    # fork 전에 모델 컴파일
one <- function(j) run_trial_oc_ext(j, p, design, scen_bnd, MASTER_SEED, rv$wt_spec, model_id = p$model_id)
n_checked <- 0L
for (b0 in seq(trial_from, trial_to, by = batch)) {
  ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done)
  if (!length(ids)) next
  t0 <- Sys.time()
  res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
  if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
  be <- rbindlist(lapply(res, `[[`, "be"))
  stopifnot(nrow(be) == length(ids) * n_rows_trial)
  n_checked <- n_checked + verify_batch(be, ids)
  be[, `:=`(GMR = signif(GMR, 7), CI_lower = signif(CI_lower, 7), CI_upper = signif(CI_upper, 7))]
  tmp <- tempfile(fileext = ".csv.gz", tmpdir = out_dir)
  fwrite(be, tmp, col.names = !file.exists(out_f))
  if (file.exists(out_f)) { if (!file.append(out_f, tmp)) stop("file.append 실패: ", out_f); unlink(tmp) } else file.rename(tmp, out_f)
  say(sprintf("trials %d-%d (%d scenarios, %d endpoints) done in %s; stored-row check passed", min(ids), max(ids), length(SCN), length(OC_ENDPOINTS_EXT), format(Sys.time() - t0)))
}
say(sprintf("done: %d stored rows matched (6 original endpoints), output %s", n_checked, out_f))
