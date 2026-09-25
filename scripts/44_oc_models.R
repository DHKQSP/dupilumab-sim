#!/usr/bin/env Rscript
# 44_oc_models.R — 분석 모형 재판정용 시험 재생성 (추가 지시 2026-09-26 §1, config/prereg_20260926.yaml section1, D-053).
# scripts/31·40·42와 같은 시드 규칙(oc_design.yaml trials$master_seed, derive_seed(master, j, "subj"/"alloc"/arm "jitter"/arm "eps"))과
# 같은 역산 배율로 시험을 다시 만들고(R/oc_models.R run_trial_oc_models: 배정 층·체중을 가진 대상자 수준 자료), 시나리오 × 평가변수
# (OC_ENDPOINTS_EXT 8개)마다 M0(pooled t), M1(처리 + 체중 층), M2(처리 + log 체중)로 판정한다.
#  mode main: 시나리오 S00 + 경계 8개 + 제품 F097(시험군 F × 0.97), 시험 1–reps_boundary(10,000).
#  mode ext : 사전 등록 연장 규칙(prereg section1$extension: 사전 고정 10,000회의 P2 Wilson 95% 구간이 5%를 포함하면 20,000회, 분석 모형 M0 또는 M1)으로
#             고른 경계 칸만 시험 10,001–20,000(시나리오 목록 = 그 칸 하나, scripts/42와 같음).
# 검사(쓰기 전, 묶음마다; 하나라도 다르면 중단): 재생성한 M0를 저장본과 대조
#  - S00·경계 8개 평가변수: results/oc/oc_rejudge_be_<model>.csv.gz (signif 7 → 상대 차이 ≤ 1e-6, pass·n_R·n_T 동일; scripts/40과 같음)
#  - 연장: results/oc/oc_trials_ext_be_<model>.csv.gz가 그 칸을 가지면 같은 방식으로 대조. 없으면 시험 10,000을 한 칸 목록으로 다시 만들어 저장본과 대조(scripts/42와 같음)
#  - F097 원래 6개 평가변수: k2016은 results/trials5000/products5000_be_raw_base.csv.gz(시험 1–5,000), k2020은 results/trials/products_be_raw_struct2020.csv.gz
#    (시험 1–500), method pooled_t. 이름 대응 AUCinf_all = AUCinf_B, AUCinf_reliable = AUCinf_A, AUCinf_subC = AUCinf_C. 약 15자리 저장 → 상대 차이 ≤ 1e-12.
# 산출(prereg section1$out_dir, 기본 results/oc_models): oc_models_be_<model>.csv.gz, oc_models_ext_be_<model>.csv.gz
#   열: trial, scenario, endpoint, model, est(log GMR, signif 8), se(signif 8), pass, n_R, n_T. df = n_R + n_T − 2(M0), − 3(M1·M2).
#   GMR = exp(est), 90% CI = exp(est ± qt(0.95, df)·se). pass는 반올림 전 값으로 판정.
#   ext mode는 extension_decision_models.csv(모델 × 경계 칸 × 분석 모형 M0·M1의 10,000회 P2와 선택)를 함께 쓴다.
#   500회 묶음 체크포인트(임시 gz → file.append, 재시작 시 완료 시험 건너뜀).
# 사용법: nice -n 15 Rscript scripts/44_oc_models.R <k2016|k2020> [cores=3] [mode=main|ext] [trial_from] [trial_to] [out_dir]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("k2016", "k2020")) stop("사용법: Rscript scripts/44_oc_models.R <k2016|k2020> [cores] [main|ext] [trial_from] [trial_to] [out_dir]")
model <- args[1]
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml"); pr <- read_cfg("prereg_20260926.yaml")$section1
cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
mode <- if (length(args) >= 3) args[3] else "main"
stopifnot(!is.na(cores), cores >= 1, mode %in% c("main", "ext"), model %in% unlist(pr$pk_models))
n_bnd <- as.integer(oc$trials$reps_boundary); batch <- as.integer(oc$trials$batch); MASTER_SEED <- as.integer(oc$trials$master_seed)
stopifnot(identical(as.integer(unlist(pr$trials)), c(1L, n_bnd)), MASTER_SEED == as.integer(pr$master_seed))
n_to <- as.integer(pr$extension$to)
trial_from <- if (length(args) >= 4) as.integer(args[4]) else if (mode == "main") 1L else n_bnd + 1L
trial_to <- if (length(args) >= 5) as.integer(args[5]) else if (mode == "main") n_bnd else n_to
if (mode == "main") stopifnot(trial_from >= 1, trial_to <= n_bnd, trial_to >= trial_from) else stopifnot(trial_from > n_bnd, trial_to <= n_to, trial_to >= trial_from)
default_dir <- proj_path(pr$out_dir)
out_dir <- if (length(args) >= 6) args[6] else default_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
production <- normalizePath(out_dir) == normalizePath(default_dir, mustWork = FALSE)
oc_dir <- proj_path("results", "oc")
MODELS <- unlist(pr$analysis_models_run); stopifnot(identical(MODELS, c("M0", "M1", "M2")))
EPS <- OC_ENDPOINTS_EXT
# 기준 세트 (iii)·(iv) 평가변수(prereg section4, D-057): 같은 재생성 시험에서 함께 계산해 별도 파일(oc_models_crit_be_<model>)에 M0·M1만 쓴다.
# section1 파일(평가변수 8개 × M0·M1·M2)의 내용·형식은 바뀌지 않는다.
EPS_CRIT <- OC_ENDPOINTS_CRIT; MODELS_CRIT <- unlist(read_cfg("prereg_20260926.yaml")$section4$trial$analysis_models)
stopifnot(identical(EPS_CRIT, unlist(read_cfg("prereg_20260926.yaml")$section4$trial$endpoints_added)), all(MODELS_CRIT %in% MODELS))
COLS <- c("trial", "scenario", "endpoint", "model", "est", "se", "pass", "n_R", "n_T")

variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p

# 시나리오: scripts/40과 같은 규칙(역산 결과의 도달 행, 코드 = 기전_방향_목표×100), 시나리오 표와 배율 대조
inv <- rbindlist(lapply(list.files(oc_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", model), full.names = TRUE), fread))
stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
inv <- inv[reachable == TRUE]
inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
s0 <- fread(file.path(oc_dir, sprintf("oc_scenarios_%s.csv", model)))
mm <- merge(inv[, .(code, multiplier)], s0[, .(code, m0 = multiplier)], by = "code", all = TRUE)
if (nrow(mm) != nrow(inv) || anyNA(mm) || any(abs(mm$multiplier / mm$m0 - 1) > 1e-12)) stop("oc_scenarios와 역산 결과의 배율이 다릅니다")
bnd <- as.numeric(unlist(oc$boundary_targets))
bcodes <- inv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9, code]
stopifnot(length(bcodes) == 8L)
bscen <- setNames(lapply(bcodes, function(cd) { r <- inv[code == cd]; list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }), bcodes)
prod_scen <- lapply(pr$scenarios$products, function(m) list(T_multipliers = m))
stopifnot(identical(names(prod_scen), "F097"), identical(prod_scen$F097$T_multipliers, list(F = 0.97)))
sc_yaml <- load_scenarios()$scenarios$F097$T_multipliers
if (!identical(as.numeric(sc_yaml$F), 0.97)) stop("config/scenarios.yaml F097이 F × 0.97이 아닙니다")

rel_ok <- function(a, b, tol) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a / b - 1) <= tol)
same_pass <- function(a, b) (is.na(a) & is.na(b)) | (a %in% TRUE & b %in% TRUE) | (a %in% FALSE & b %in% FALSE)
m0_ci <- function(x) {                                                # M0 행의 GMR·CI(반올림 전 est·se에서)
  q <- qt(1 - (1 - design$be$ci_level) / 2, x$n_R + x$n_T - 2)
  x[, .(trial, scenario, endpoint, GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se), pass, n_R, n_T)]
}
verify <- function(new, ref, tol, what) {                            # new, ref: trial, scenario, endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T
  m <- merge(new, ref, by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
  if (nrow(m) != nrow(new) || nrow(m) != nrow(ref)) stop(sprintf("%s: 행 수 불일치(재생성 %d, 저장본 %d, 병합 %d)", what, nrow(new), nrow(ref), nrow(m)))
  ok <- (rel_ok(m$GMR, m$GMR.s, tol) & rel_ok(m$CI_lower, m$CI_lower.s, tol) & rel_ok(m$CI_upper, m$CI_upper.s, tol) & same_pass(m$pass, m$pass.s) &
           m$n_R == m$n_R.s & m$n_T == m$n_T.s) %in% TRUE
  if (!all(ok)) { print(head(m[!ok], 10)); stop(sprintf("%s: %d행이 저장본과 다릅니다(시험 %s). 결과를 쓰지 않습니다.", what, sum(!ok), paste(head(unique(m$trial[!ok]), 5), collapse = ","))) }
  nrow(m)
}
PROD_MAP <- c(Cmax = "Cmax", AUClast = "AUClast", AUCinf_all = "AUCinf_B", AUCinf_reliable = "AUCinf_A", AUCinf_true = "AUCinf_true", AUCinf_subC = "AUCinf_C")
read_products_f097 <- function() {
  f <- if (model == "k2016") proj_path("results", "trials5000", "products5000_be_raw_base.csv.gz") else proj_path("results", "trials", "products_be_raw_struct2020.csv.gz")
  x <- fread(f)[scenario == "F097" & method == "pooled_t" & schedule == "B0"]
  if (anyDuplicated(x, by = c("trial", "endpoint"))) stop("제품 저장본 F097 중복: ", f)
  x[, endpoint := unname(PROD_MAP[endpoint])]
  list(file = f, rows = x[, .(trial, scenario, endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T)])
}
append_gz <- function(dt, f) {
  tmp <- tempfile(fileext = ".csv.gz", tmpdir = dirname(f))
  fwrite(dt, tmp, col.names = !file.exists(f))
  if (file.exists(f)) { if (!file.append(f, tmp)) stop("file.append 실패: ", f); unlink(tmp) } else if (!file.rename(tmp, f)) stop("file.rename 실패: ", f)
}
keep_only <- function(f, keep, by_scen = FALSE) {                     # 두 파일에 모두 완전한 것만 남긴다(한쪽만 붙은 묶음 제거)
  if (!file.exists(f)) return(invisible(NULL)); ex <- fread(f); keys <- if (by_scen) c("trial", "scenario") else "trial"
  n0 <- nrow(ex); ex <- if (by_scen) ex[keep, on = keys, nomatch = NULL] else ex[trial %in% keep]
  if (nrow(ex) != n0) { message(sprintf("%s: 다른 파일과 맞지 않는 행 %d개를 지웁니다", basename(f), n0 - nrow(ex))); tmp <- paste0(f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, f) }
  invisible(NULL)
}
done_trials <- function(f, n_rows_trial, by_scen = FALSE) {         # 완전한 시험(또는 (시험, 시나리오))만 남기고 불완전 행은 지운다
  if (!file.exists(f)) return(if (by_scen) data.table(trial = integer(0), scenario = character(0)) else integer(0))
  ex <- fread(f)
  if (!identical(names(ex), COLS)) stop("기존 파일의 열이 다릅니다: ", f)
  keys <- if (by_scen) c("trial", "scenario") else "trial"
  cn <- ex[, .N, by = keys]; bad <- cn[N != n_rows_trial]
  if (nrow(bad) || anyDuplicated(ex, by = c("trial", "scenario", "endpoint", "model"))) {
    message(sprintf("기존 파일의 불완전 묶음 %d개를 지우고 다시 씁니다: %s", nrow(bad), f))
    ex <- unique(ex[!bad, on = keys], by = c("trial", "scenario", "endpoint", "model"))
    tmp <- paste0(f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, f)
  }
  if (by_scen) unique(ex[, .(trial, scenario)]) else unique(ex$trial)
}
finish_rows <- function(be) {
  stopifnot(identical(unique(be$lloq), study_lloq()), "resid" %in% names(be) && identical(unique(be$resid), "fixed"))
  be[, c("lloq", "resid") := NULL]
  be[, `:=`(est = signif(est, 8), se = signif(se, 8))]
  be[, ..COLS]
}
split_rows <- function(be) list(main = finish_rows(be[endpoint %in% EPS]), crit = finish_rows(be[endpoint %in% EPS_CRIT & model %in% MODELS_CRIT]))
invisible(get_model(p$model_id))                                     # fork 전에 모델 컴파일

if (mode == "main") {
  scen <- c(list(S00 = list(code = "S00", T_multipliers = list())), bscen, prod_scen)
  SCN <- names(scen); n_rows_trial <- length(SCN) * length(EPS) * length(MODELS); n_rows_crit <- length(SCN) * length(EPS_CRIT) * length(MODELS_CRIT)
  stored <- fread(file.path(oc_dir, sprintf("oc_rejudge_be_%s.csv.gz", model)))[trial >= trial_from & trial <= trial_to]
  if (!setequal(unique(stored$scenario), setdiff(SCN, "F097")) || nrow(stored) != (trial_to - trial_from + 1L) * (length(SCN) - 1L) * length(EPS))
    stop("oc_rejudge_be 저장본이 요청 범위의 S00·경계 × 8개 평가변수를 모두 갖고 있지 않습니다")
  pf <- read_products_f097(); prod_rows <- pf$rows[trial >= trial_from & trial <= trial_to]
  out_f <- file.path(out_dir, sprintf("oc_models_be_%s.csv.gz", model)); crit_f <- file.path(out_dir, sprintf("oc_models_crit_be_%s.csv.gz", model))
  done <- intersect(done_trials(out_f, n_rows_trial), done_trials(crit_f, n_rows_crit))
  keep_only(out_f, done); keep_only(crit_f, done)
  logfile <- if (production) start_run_log(paste0("oc_models_", model), master_seed = MASTER_SEED, run_mode = "final",
                                            extra = list(model = model, mode = mode, scenarios = paste(SCN, collapse = ","), trial_from = trial_from, trial_to = trial_to,
                                                         cores = cores, verify_f097 = basename(pf$file), verify_f097_trials = if (nrow(prod_rows)) paste(range(prod_rows$trial), collapse = "-") else "none")) else NULL
  say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }
  say(sprintf("%s main: %d scenarios (%s), trials %d-%d, cores %d, output %s", model, length(SCN), paste(SCN, collapse = ","), trial_from, trial_to, cores, out_f))
  one <- function(j) run_trial_oc_models(j, p, design, scen, MASTER_SEED, rv$wt_spec, model_id = p$model_id, models = MODELS, endpoints = c(EPS, EPS_CRIT))$be
  n_chk <- c(rejudge = 0L, products = 0L)
  for (b0 in seq(trial_from, trial_to, by = batch)) {
    ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done)
    if (!length(ids)) next
    t0 <- Sys.time()
    res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
    bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.table(x), logical(1))
    if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
    be <- rbindlist(res)
    stopifnot(nrow(be) == length(ids) * length(SCN) * length(c(EPS, EPS_CRIT)) * length(MODELS))
    m0 <- m0_ci(be[model == "M0" & endpoint %in% EPS])
    n_chk["rejudge"] <- n_chk["rejudge"] + verify(m0[scenario != "F097"], stored[trial %in% ids], 1e-6, sprintf("%s oc_rejudge_be 대조", model))
    pr_ids <- intersect(ids, prod_rows$trial)
    if (length(pr_ids)) n_chk["products"] <- n_chk["products"] + verify(m0[scenario == "F097" & trial %in% pr_ids & endpoint %in% PROD_MAP], prod_rows[trial %in% pr_ids], 1e-12, sprintf("%s F097 제품 저장본 대조", model))
    sp <- split_rows(be); stopifnot(nrow(sp$main) == length(ids) * n_rows_trial, nrow(sp$crit) == length(ids) * n_rows_crit)
    append_gz(sp$crit, crit_f); append_gz(sp$main, out_f)
    say(sprintf("trials %d-%d (%d scenarios x %d endpoints x %d models; criteria sets iii-iv in the crit file) done in %s; M0 matched the stored rows (oc_rejudge_be%s)", min(ids), max(ids), length(SCN), length(EPS), length(MODELS),
                format(Sys.time() - t0), if (length(pr_ids)) sprintf(", F097 products trials %d-%d", min(pr_ids), max(pr_ids)) else ""))
  }
  say(sprintf("done: M0 matched %d stored oc_rejudge_be rows and %d stored F097 product rows; output %s", n_chk["rejudge"], n_chk["products"], out_f))
} else {
  # 연장 판정: 사전 고정 10,000회(main 산출)의 P2, 분석 모형 prereg extension$analysis_models 각각. 하나라도 규칙에 걸리면 그 칸을 연장
  main_f <- file.path(default_dir, sprintf("oc_models_be_%s.csv.gz", model))
  if (!file.exists(main_f)) stop("main 산출 없음: ", main_f)
  cfg_eps <- OC_REJUDGE_CONFIGS[[pr$extension$config]]
  mb <- fread(main_f)[scenario %in% bcodes & model %in% unlist(pr$extension$analysis_models) & endpoint %in% cfg_eps]
  setnames(mb, "model", "analysis_model")
  wide <- dcast(mb[, .(trial, scenario, analysis_model, endpoint, ok = pass %in% TRUE)], trial + scenario + analysis_model ~ endpoint, value.var = "ok")
  wide[, cfg := Reduce(`&`, .SD), .SDcols = cfg_eps]
  dec <- wide[, .(n_before = .N, n_pass = sum(cfg)), by = .(scenario, analysis_model)]
  if (any(dec$n_before != n_bnd) || nrow(dec) != length(bcodes) * length(unlist(pr$extension$analysis_models))) stop("main 산출이 경계 칸 × 분석 모형마다 10,000회가 아닙니다")
  dec[, c("pass_pct", "lo", "hi") := wilson_ci(n_pass, n_before)]
  thr <- as.numeric(pr$extension$threshold_pct)
  dec[, triggers := lo <= thr & hi >= thr]
  dec[, pk_model := model]
  selc <- dec[, .(selected = any(triggers)), by = scenario]
  dec <- merge(dec, selc, by = "scenario"); dec[, n_after := fifelse(selected, n_to, n_bnd)]
  dec[, rule := pr$extension$rule]
  setcolorder(dec, c("pk_model", "scenario", "analysis_model", "n_before", "n_pass", "pass_pct", "lo", "hi", "triggers", "selected", "n_after", "rule"))
  dec_f <- file.path(out_dir, "extension_decision_models.csv")
  old <- if (file.exists(dec_f)) fread(dec_f)[pk_model != model] else NULL
  fwrite(rbind(old, dec, fill = TRUE), dec_f)
  sel <- selc[selected == TRUE, scenario]
  logfile <- if (production) start_run_log(paste0("oc_models_ext_", model), master_seed = MASTER_SEED, run_mode = "final",
                                            extra = list(model = model, mode = mode, rule = pr$extension$rule, selected = if (length(sel)) paste(sel, collapse = ",") else "none",
                                                         trial_from = trial_from, trial_to = trial_to, cores = cores)) else NULL
  say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }
  for (i in seq_len(nrow(dec))) say(sprintf("%s %s %s: P2 %.2f%% [%.2f, %.2f] at %d trials -> %s", model, dec$scenario[i], dec$analysis_model[i], dec$pass_pct[i], dec$lo[i], dec$hi[i], dec$n_before[i],
                                            if (dec$triggers[i]) "rule met" else "rule not met"))
  if (!length(sel)) { say("no boundary cell selected; nothing to run"); quit(save = "no") }
  out_f <- file.path(out_dir, sprintf("oc_models_ext_be_%s.csv.gz", model)); crit_f <- file.path(out_dir, sprintf("oc_models_ext_crit_be_%s.csv.gz", model))
  n_rows_trial <- length(EPS) * length(MODELS); n_rows_crit <- length(EPS_CRIT) * length(MODELS_CRIT)
  done <- fintersect(done_trials(out_f, n_rows_trial, by_scen = TRUE), done_trials(crit_f, n_rows_crit, by_scen = TRUE))
  keep_only(out_f, done, TRUE); keep_only(crit_f, done, TRUE)
  ext_f <- file.path(oc_dir, sprintf("oc_trials_ext_be_%s.csv.gz", model))
  ext_st <- if (file.exists(ext_f)) fread(ext_f) else NULL
  rej <- fread(file.path(oc_dir, sprintf("oc_rejudge_be_%s.csv.gz", model)))[trial == n_bnd]
  for (sc_ in sel) {
    scen <- bscen[sc_]
    one <- function(j) run_trial_oc_models(j, p, design, scen, MASTER_SEED, rv$wt_spec, model_id = p$model_id, models = MODELS, endpoints = c(EPS, EPS_CRIT))$be
    chk <- one(n_bnd)                                                  # 한 칸 목록으로 시험 10,000을 다시 만들어 저장본과 대조(scripts/42와 같은 전제 검사)
    verify(m0_ci(chk[model == "M0" & endpoint %in% EPS]), rej[scenario == sc_], 1e-6, sprintf("%s %s 시험 %d 재현", model, sc_, n_bnd))
    st_sc <- if (!is.null(ext_st)) ext_st[scenario == sc_] else NULL
    say(sprintf("%s %s: trial %d reproduced; extending trials %d-%d%s", model, sc_, n_bnd, trial_from, trial_to,
                if (!is.null(st_sc) && nrow(st_sc)) sprintf(" (M0 checked against %s, stored trials %d-%d)", basename(ext_f), min(st_sc$trial), max(st_sc$trial)) else " (no stored extension rows for M0)"))
    for (b0 in seq(trial_from, trial_to, by = batch)) {
      ids <- setdiff(b0:min(b0 + batch - 1L, trial_to), done[scenario == sc_, trial])
      if (!length(ids)) next
      t0 <- Sys.time()
      res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
      bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.table(x), logical(1))
      if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","), " — ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = "; "))
      be <- rbindlist(res)
      stopifnot(nrow(be) == length(ids) * length(c(EPS, EPS_CRIT)) * length(MODELS), all(be$scenario == sc_))
      st_ids <- if (!is.null(st_sc)) intersect(ids, st_sc$trial) else integer(0)
      if (length(st_ids)) verify(m0_ci(be[model == "M0" & endpoint %in% EPS & trial %in% st_ids]), st_sc[trial %in% st_ids], 1e-6, sprintf("%s %s 연장 저장본 대조", model, sc_))
      sp <- split_rows(be); append_gz(sp$crit, crit_f); append_gz(sp$main, out_f)
      done <- rbind(done, data.table(trial = ids, scenario = sc_))
      say(sprintf("%s %s: trials %d-%d done in %s%s", model, sc_, min(ids), max(ids), format(Sys.time() - t0), if (length(st_ids)) "; M0 matched the stored extension rows" else ""))
    }
  }
  say(sprintf("done: extension output %s", out_f))
}
