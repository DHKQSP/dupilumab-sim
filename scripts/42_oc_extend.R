#!/usr/bin/env Rscript
# 42_oc_extend.R — 경계 1종 오류의 적응적 연장 (지시 2026-09-25 §3 마지막 항목, D-040).
# 규칙(모든 모델에 일반 적용; scripts/20_products_5000.R의 5,000회 제품 시나리오와 같은 적응적 상향, R/oc.R oc_extension_select):
#   results/oc/boundary_type1.csv(scripts/33)의 P2 경계 1종 오류 Wilson 95% 구간이 사전 고정 reps_boundary(10,000)회에서 5%를 포함하면
#   (lo ≤ 5 ≤ hi) 그 (모델, 경계 시나리오)를 to(기본 20,000)회로 늘린다. 선택은 표에서 정한다(시나리오를 코드에 적지 않는다).
#   boundary_type1.csv가 이미 연장분을 반영했으면(n_trials > reps_boundary) 그 행의 10,000회 값을 저장본(시험 1–reps_boundary)에서 다시 계산해 판정한다.
# 모의: 선택된 (모델, 시나리오)마다 시험 reps_boundary+1 .. to를 run_trial_oc_ext(평가변수 8개 = 원래 6개 + 세트 (i) AUCinf_Ai·AUCinf_Ci)로,
#   시나리오 목록 = 그 시나리오 하나. 대조군은 시험마다 모의되고 한 시나리오의 결과는 시나리오 목록과 무관하다(scripts/40 머리말).
#   변형·체중 분포·마스터 시드는 scripts/31과 같다(k2016 = base, k2020 = struct2020; oc_design.yaml trials$master_seed).
#   전제 검사(쓰기 전): 시험 reps_boundary(31이 경계 목록으로 돌린 마지막 시험)를 같은 방식(시나리오 하나)으로 다시 모의해 원래 6개 평가변수가
#   저장본 oc_trials_be_<model>.csv.gz와 같은지(GMR·CI 상대 차이 ≤ 1e-6, pass·n 동일) 확인하고, 다르면 중단한다.
# 산출(out_dir, 기본 results/oc; 사전 고정 10,000회 파일 oc_trials_be_<model>.csv.gz는 바꾸지 않는다):
#   extension_decision.csv             모델 × 경계 시나리오: n_before, P2 pass_pct·lo·hi(n_before회), threshold, rule, selected, n_after(목표 시험 수)
#   oc_trials_ext_be_<model>.csv.gz    trial, scenario, endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T (oc_trials_be와 같은 열, 평가변수 8개, signif 7)
#   oc_trials_ext_drop_<model>.csv.gz  run_trial_oc_ext의 탈락 수(REF·T 행, n_reliable_i 포함) + ext_scenario(그 행을 만든 연장 시나리오)
#   500회 묶음 체크포인트: 묶음마다 임시 gz를 만든 뒤 file.append로 붙인다(scripts/40과 같음; drop 먼저, be 나중). 재시작 시 be에 8개 평가변수가
#   모두 있는 (시험, 시나리오)는 건너뛰고 불완전 행·고아 drop 행은 지운다. 같은 출력 파일에 두 실행을 동시에 돌리지 않는다.
# 사용법: nice -n 15 Rscript scripts/42_oc_extend.R [cores=3] [to=20000] [out_dir=results/oc]
#   시험용: out_dir를 임시 디렉터리로, to를 작게(예: 10002). 입력(boundary_type1.csv, 저장본, 시나리오 표)은 results/oc에서 읽는다
#   (경계 표만 환경변수 DUPI_OC_BOUNDARY_TABLE로 바꿀 수 있다: 연장 반영 표의 재계산 경로 시험용).
#   실행 기록(start_run_log: 규칙, 선택된 시나리오)은 운영 실행(out_dir = results/oc)에서만 logs/에 남긴다(scripts/40과 같음).
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
n_bnd <- as.integer(oc$trials$reps_boundary); batch <- as.integer(oc$trials$batch); MASTER_SEED <- as.integer(oc$trials$master_seed)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
n_to <- if (length(args) >= 2) as.integer(args[2]) else 20000L
default_dir <- proj_path("results", "oc")
out_dir <- if (length(args) >= 3) args[3] else default_dir
if (is.na(cores) || cores < 1 || is.na(n_to) || n_to <= n_bnd)
  stop(sprintf("사용법: Rscript scripts/42_oc_extend.R [cores] [to > %d] [out_dir]", n_bnd))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
production <- normalizePath(out_dir) == normalizePath(default_dir)
RULE_CFG <- "P2"; THRESH <- 5
COLS <- c("trial", "scenario", "endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T")
n_rows_trial <- length(OC_ENDPOINTS_EXT)

# ----- 1) 판정: boundary_type1.csv(P2)에서 규칙 적용 ------------------------------------------------------------
bt_f <- Sys.getenv("DUPI_OC_BOUNDARY_TABLE", file.path(default_dir, "boundary_type1.csv"))   # 환경변수는 시험용(기본 results/oc/boundary_type1.csv)
if (!file.exists(bt_f)) stop("경계 표 없음: ", bt_f, " (scripts/33_oc_summary.R를 먼저 실행)")
bt <- fread(bt_f)[config == RULE_CFG]
if (!nrow(bt)) stop("boundary_type1.csv에 P2 행이 없습니다")
stored_cache <- list()
stored_be <- function(mdl) {                                         # 저장본(사전 고정 파일)은 모델당 한 번 읽는다
  if (is.null(stored_cache[[mdl]])) {
    f <- file.path(default_dir, sprintf("oc_trials_be_%s.csv.gz", mdl))
    if (!file.exists(f)) stop("저장본 없음: ", f)
    stored_cache[[mdl]] <<- fread(f)
  }
  stored_cache[[mdl]]
}
if (any(bt$n_trials < n_bnd)) stop(sprintf("경계 시나리오 %d개가 사전 고정 %d회 미만입니다(31 미완료): %s", sum(bt$n_trials < n_bnd), n_bnd,
                                           paste(bt[n_trials < n_bnd, paste(model, scenario)], collapse = ", ")))
re <- which(bt$n_trials > n_bnd)                                      # 이미 연장분이 반영된 표: 10,000회 값을 저장본에서 다시 계산
for (i in re) {
  s <- stored_be(bt$model[i]); sc_ <- bt$scenario[i]
  w <- config_pass(s[scenario == sc_ & trial <= n_bnd & endpoint %in% OC_ENDPOINTS], oc)
  if (nrow(w) != n_bnd) stop(sprintf("저장본 %s %s의 시험 수 %d ≠ %d", bt$model[i], sc_, nrow(w), n_bnd))
  wc <- wilson_ci(sum(w$cfg_P2), nrow(w)); n0 <- bt$n_trials[i]
  set(bt, i, c("n_trials", "pass_pct", "lo", "hi"), list(n_bnd, wc$est, wc$lo, wc$hi))
  message(sprintf("%s %s: boundary_type1.csv는 %d회(연장 반영) → 판정은 저장본 시험 1–%d의 P2 %.2f%% [%.2f, %.2f]", bt$model[i], sc_, n0, n_bnd, wc$est, wc$lo, wc$hi))
}
dec <- oc_extension_select(bt, n_from = n_bnd, n_to = n_to, cfg = RULE_CFG, threshold = THRESH)
fwrite(dec, file.path(out_dir, "extension_decision.csv"))
sel <- dec[selected == TRUE]
cat(sprintf("규칙: %s\n선택: %s\n", dec$rule[1], if (nrow(sel)) paste(sprintf("%s %s (P2 %.2f%% [%.2f, %.2f])", sel$model, sel$scenario, sel$pass_pct, sel$lo, sel$hi), collapse = "; ") else "없음"))
logfile <- if (production) start_run_log("oc_extend", master_seed = MASTER_SEED, run_mode = "final",
                                          extra = list(rule = dec$rule[1], selected = if (nrow(sel)) paste(sel$model, sel$scenario) else "none",
                                                       selected_p2 = if (nrow(sel)) sprintf("%s %s: %.2f%% [%.2f, %.2f] at %d", sel$model, sel$scenario, sel$pass_pct, sel$lo, sel$hi, sel$n_before) else "none",
                                                       trial_from = n_bnd + 1L, trial_to = n_to, cores = cores, out_dir = out_dir)) else NULL
say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }
for (mdl in setdiff(unique(dec$model), unique(sel$model))) {         # 선택 없는 모델에 연장 파일이 남아 있으면 규칙과 어긋난다
  f <- file.path(out_dir, sprintf("oc_trials_ext_be_%s.csv.gz", mdl))
  if (file.exists(f)) stop(sprintf("%s: 규칙이 고른 시나리오가 없는데 연장 파일이 있습니다: %s", mdl, f))
}
if (!nrow(sel)) { say("no boundary scenario selected; nothing to run"); quit(save = "no") }

# ----- 2) 선택된 (모델, 시나리오)별 연장 모의 ---------------------------------------------------------------------
rel_ok <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a / b - 1) <= 1e-6)
same_rows <- function(new, ref) {                                    # 원래 6개 평가변수: 새 모의(signif 전) 대 저장본(signif 7)
  m <- merge(new, ref, by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
  nrow(m) == nrow(new) && nrow(m) == nrow(ref) &&
    all((rel_ok(m$GMR, m$GMR.s) & rel_ok(m$CI_lower, m$CI_lower.s) & rel_ok(m$CI_upper, m$CI_upper.s) &
         ((is.na(m$pass) & is.na(m$pass.s)) | (m$pass %in% TRUE & m$pass.s %in% TRUE) | (m$pass %in% FALSE & m$pass.s %in% FALSE)) &
         m$n_R == m$n_R.s & m$n_T == m$n_T.s) %in% TRUE)
}
append_gz <- function(dt, f) {                                       # 임시 gz를 만든 뒤 붙인다(헤더는 새 파일에만)
  tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f))
  fwrite(dt, tmp, col.names = !file.exists(f))
  if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f)
}
for (mdl in unique(sel$model)) {
  variant <- if (mdl == "k2016") "base" else "struct2020"
  rv <- resolve_variant(variant, design); p <- rv$p
  invisible(get_model(p$model_id))                                  # fork 전에 모델 컴파일
  sc_tab <- fread(file.path(default_dir, sprintf("oc_scenarios_%s.csv", mdl)))
  be_f <- file.path(out_dir, sprintf("oc_trials_ext_be_%s.csv.gz", mdl)); drop_f <- file.path(out_dir, sprintf("oc_trials_ext_drop_%s.csv.gz", mdl))
  scn_m <- sel[model == mdl, scenario]
  # 기존 연장 파일: 열·시나리오 확인, 불완전 (시험, 시나리오)와 고아 drop 행 정리
  done <- data.table(trial = integer(0), scenario = character(0))
  if (file.exists(be_f)) {
    ex <- fread(be_f)
    if (!identical(names(ex), COLS)) stop("기존 연장 파일의 열이 다릅니다: ", be_f)
    if (!all(ex$scenario %in% scn_m)) stop(sprintf("기존 연장 파일 %s에 규칙이 고르지 않은 시나리오가 있습니다: %s", be_f, paste(setdiff(unique(ex$scenario), scn_m), collapse = ",")))
    if (any(ex$trial <= n_bnd)) stop("기존 연장 파일에 사전 고정 범위(≤ reps_boundary)의 시험이 있습니다: ", be_f)
    cn <- ex[, .(N = .N, n_ep = uniqueN(endpoint)), by = .(trial, scenario)]
    bad <- cn[N != n_rows_trial | n_ep != n_rows_trial]
    if (nrow(bad)) {
      message(sprintf("기존 연장 파일의 불완전 (시험, 시나리오) %d개(예: %s)를 지우고 다시 씁니다", nrow(bad), paste(head(paste(bad$trial, bad$scenario), 3), collapse = ", ")))
      ex <- ex[!bad, on = c("trial", "scenario")]
      tmp <- paste0(be_f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, be_f)
    }
    done <- unique(ex[, .(trial, scenario)]); rm(ex)
  }
  if (file.exists(drop_f)) {
    dr <- fread(drop_f)
    keep <- dr[done, on = c(trial = "trial", ext_scenario = "scenario"), nomatch = NULL]
    miss <- done[!dr, on = c(trial = "trial", scenario = "ext_scenario")]
    if (nrow(miss)) stop(sprintf("연장 be 파일의 완료 시험 %d개에 drop 행이 없습니다(예: %s)", nrow(miss), paste(head(paste(miss$trial, miss$scenario), 3), collapse = ", ")))
    if (nrow(keep) != nrow(dr)) {
      message(sprintf("연장 drop 파일의 고아 행 %d개를 지웁니다", nrow(dr) - nrow(keep)))
      tmp <- paste0(drop_f, ".rewrite.csv.gz"); fwrite(keep, tmp); file.rename(tmp, drop_f)
    }
    rm(dr, keep)
  } else if (nrow(done)) stop("연장 be 파일은 있는데 drop 파일이 없습니다: ", drop_f)

  for (sc_ in scn_m) {
    r <- sc_tab[code == sc_]
    if (nrow(r) != 1) stop(sprintf("oc_scenarios_%s.csv에 %s가 없습니다", mdl, sc_))
    d_ <- sel[model == mdl & scenario == sc_]
    if (abs(r$multiplier / d_$multiplier - 1) > 1e-9) stop(sprintf("%s %s: 시나리오 표와 경계 표의 배율이 다릅니다", mdl, sc_))
    scen <- setNames(list(list(code = sc_, T_multipliers = setNames(list(r$multiplier), r$mechanism))), sc_)
    one <- function(j) run_trial_oc_ext(j, p, design, scen, MASTER_SEED, rv$wt_spec, model_id = p$model_id)
    # 사전 고정 파일: 이 시나리오는 시험 1–n_bnd가 모두 있어야 하고(연장은 그 뒤에 붙는다), 시험 n_bnd의 재현으로 시드·변형·배율을 확인
    st <- stored_be(mdl)[scenario == sc_ & endpoint %in% OC_ENDPOINTS]
    tr_ <- sort(unique(st$trial))
    if (!identical(tr_, seq_len(n_bnd)) || nrow(st) != n_bnd * length(OC_ENDPOINTS))
      stop(sprintf("저장본 %s %s: 시험 1–%d의 원래 평가변수가 완전하지 않습니다(시험 %d개, 행 %d개)", mdl, sc_, n_bnd, length(tr_), nrow(st)))
    chk <- one(n_bnd)
    if (!same_rows(chk$be[endpoint %in% OC_ENDPOINTS], st[trial == n_bnd]))
      stop(sprintf("%s %s: 시험 %d 재현이 저장본과 다릅니다(시드·변형·배율 확인). 연장을 쓰지 않습니다.", mdl, sc_, n_bnd))
    say(sprintf("%s %s: trial %d reproduced from the stored file (6 original endpoints); extending trials %d-%d (P2 at %d trials %.2f%% [%.2f, %.2f])",
                mdl, sc_, n_bnd, n_bnd + 1L, n_to, n_bnd, d_$pass_pct, d_$lo, d_$hi))
    for (b0 in seq(n_bnd + 1L, n_to, by = batch)) {
      ids <- setdiff(b0:min(b0 + batch - 1L, n_to), done[scenario == sc_, trial])
      if (!length(ids)) next
      t0 <- Sys.time()
      res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
      bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$be), logical(1))
      if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
      be <- rbindlist(lapply(res, `[[`, "be")); dr <- rbindlist(lapply(res, `[[`, "drop"))
      stopifnot(identical(names(be), COLS), nrow(be) == length(ids) * n_rows_trial, all(be$scenario == sc_), setequal(unique(be$trial), ids),
                all(dr$scenario %in% c("REF", sc_)), nrow(dr) == 2L * length(ids))
      be[, `:=`(GMR = signif(GMR, 7), CI_lower = signif(CI_lower, 7), CI_upper = signif(CI_upper, 7))]
      dr[, ext_scenario := sc_]
      append_gz(dr, drop_f); append_gz(be, be_f)                    # drop 먼저: 중단되면 재시작 시 고아 drop 행을 지운다
      done <- rbind(done, data.table(trial = ids, scenario = sc_))
      say(sprintf("%s %s: trials %d-%d (%d endpoints) done in %s", mdl, sc_, min(ids), max(ids), n_rows_trial, format(Sys.time() - t0)))
    }
  }
}
for (i in seq_len(nrow(sel))) {                                      # 완료 요약: 연장 파일의 시험 수 대 목표
  f <- file.path(out_dir, sprintf("oc_trials_ext_be_%s.csv.gz", sel$model[i]))
  n_ext <- if (file.exists(f)) uniqueN(fread(f)[scenario == sel$scenario[i], trial]) else 0L
  say(sprintf("done: %s %s, %d + %d = %d trials (target %d), output %s", sel$model[i], sel$scenario[i], n_bnd, n_ext, n_bnd + n_ext, sel$n_after[i], f))
}
