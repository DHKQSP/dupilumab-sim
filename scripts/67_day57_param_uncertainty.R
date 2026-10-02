#!/usr/bin/env Rscript
# 67_day57_param_uncertainty.R — 파라미터 불확실성(2절), 시험 수준 중첩(3절 상한), ①② 모델 불확실성 범위(B)
#   (사전 등록 config/prereg_20261002_day57.yaml section2·section3 nested·addendum_B, D-066). 입력: config/param_uncertainty.yaml(논문 추정값과 RSE).
# 입력이 완전하지 않으면(status pending, 표본을 뽑을 항목에 rse도 not_reported도 없음) 멈춘다. 시험용 입력은 환경변수 DUPI_PU_CONFIG로만 바꿀 수 있고,
# 그때 산출은 DUPI_PU_OUT(필수)에 쓴다(정식 결과 폴더에는 쓰지 않음).
# 모드
#   sets <k2016|k2020> [from] [to]        세트 from..to(기본 1..1000; 예비 실행 1..200). 세트 0 = 논문 추정값(기준). 25세트 묶음 파일(있으면 건너뜀)
#   nested <k2016|k2020> [from] [to]      세트 from..to(기본 1..200) x 시험 200회(3절 설계), 세트별 파일(있으면 건너뜀)
#   summary [max_set]                     세트 간 분포(중앙값, 2.5·97.5 백분위)와 중첩 예측 분포
# 사용법: nice -n 19 Rscript scripts/67_day57_param_uncertainty.R <mode> ...
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("sets", "nested", "summary")) stop("사용법: Rscript scripts/67_day57_param_uncertainty.R <sets|nested|summary> ...")
mode <- args[1]
PR <- read_cfg("prereg_20261002_day57.yaml"); design <- read_cfg("trial_design.yaml")
test_cfg <- Sys.getenv("DUPI_PU_CONFIG", ""); test_out <- Sys.getenv("DUPI_PU_OUT", "")
if (nzchar(test_cfg) && !nzchar(test_out)) stop("시험용 입력(DUPI_PU_CONFIG)에는 DUPI_PU_OUT이 필요합니다")
PU <- if (nzchar(test_cfg)) yaml::read_yaml(test_cfg) else read_cfg("param_uncertainty.yaml")
out_dir <- if (nzchar(test_out)) test_out else proj_path(PR$output_dir)
for (d in c("psets", "nested")) dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
MS <- 20261002L; LLOQ <- study_lloq(); stopifnot(isTRUE(all.equal(LLOQ, 0.078)))
DAYS <- as.numeric(unlist(PR$timepoints_study_day)); EL <- DAYS - 1
VAR <- c(k2016 = "base", k2020 = "struct2020")
NSUB <- as.integer(PR$section2_parameter_uncertainty$per_set$n_subjects); stopifnot(NSUB == 2000L)

# ---- 입력 검사: 표본을 뽑을 항목 목록(순서 고정 = 설정 파일 순서) ----
pu_entries <- function(model) {
  m <- PU$models[[model]]; if (is.null(m)) stop("입력 없음: ", model)
  grp <- intersect(c("theta", "omega2", "omega_sd", "sigma"), names(m))
  rbindlist(lapply(grp, function(g) rbindlist(lapply(names(m[[g]]), function(nm) { e <- m[[g]][[nm]]
    data.table(group = g, name = nm, estimate = as.numeric(e$estimate), fixed = isTRUE(e$fixed), not_reported = isTRUE(e$not_reported),
               rse = if (is.null(e$rse)) NA_real_ else as.numeric(e$rse)) }))))
}
check_inputs <- function(model, p) {
  if (identical(PU$status, "pending")) stop("config/param_uncertainty.yaml status가 pending입니다: 논문의 RSE를 넣고 커밋한 뒤 실행하십시오(D-066)")
  E <- pu_entries(model)
  miss <- E[!fixed & !not_reported & !is.finite(rse)]
  if (nrow(miss)) stop("RSE가 없는 항목(rse 또는 not_reported 필요): ", paste(model, miss$group, miss$name, collapse = "; "))
  if (any(E$rse[is.finite(E$rse)] <= 0 | E$rse[is.finite(E$rse)] > 2)) stop("RSE 범위 밖(0 < rse <= 2, 비율로 입력)")
  # 추정값이 모델 config와 같은지
  cur <- function(g, nm) switch(g,
    theta = if (nm == "theta_WT") p$cov$theta_WT else if (nm == "Mpc") p$theta[["k12"]] / p$theta[["k21"]] else p$theta[[nm]],
    omega2 = p$omega[[nm]]^2, omega_sd = p$omega[[nm]], sigma = p$sigma[[nm]])
  for (i in seq_len(nrow(E))) { v <- cur(E$group[i], E$name[i]); tol <- if (E$name[i] == "Mpc") 2e-3 else 1e-3
    if (abs(v / E$estimate[i] - 1) > tol) stop(sprintf("추정값 불일치 %s %s %s: config %g, 입력 %g", model, E$group[i], E$name[i], v, E$estimate[i])) }
  E
}
# 세트 s의 파라미터: P_s = P_hat * exp(z * sdlog), sdlog = sqrt(log(1 + RSE^2)); 독립; 고정·미보고 항목은 그대로. 세트 0 = 논문 추정값
param_set <- function(p, E, model, s) {
  ps <- p; samp <- E[!fixed & !not_reported]
  z <- if (s == 0) rep(0, nrow(samp)) else with_seed(derive_seed(MS, "pset", model, s), rnorm(nrow(samp)))
  val <- samp$estimate * exp(z * sqrt(log(1 + samp$rse^2)))
  for (i in seq_len(nrow(samp))) { g <- samp$group[i]; nm <- samp$name[i]; v <- val[i]
    if (g == "theta") { if (nm == "theta_WT") ps$cov$theta_WT <- v else if (nm == "Mpc") NULL else ps$theta[[nm]] <- v }
    if (g == "omega2") ps$omega[[nm]] <- sqrt(v)
    if (g == "omega_sd") ps$omega[[nm]] <- v
    if (g == "sigma") ps$sigma[[nm]] <- v }
  if (model == "k2020" && "Mpc" %in% E$name) {                          # k21 = k12 / Mpc (논문의 모수화)
    mpc <- if ("Mpc" %in% samp$name) val[samp$name == "Mpc"] else E[name == "Mpc", estimate]
    ps$theta[["k21"]] <- ps$theta[["k12"]] / mpc }
  if (!(ps$theta[["F"]] > 0 && ps$theta[["F"]] < 1)) stop(sprintf("세트 %d: F = %g가 (0, 1) 밖입니다(사전 등록에 처리 규칙이 없어 멈춤)", s, ps$theta[["F"]]))
  list(p = ps, values = as.list(setNames(val, paste(samp$group, samp$name, sep = "."))))
}
params_of <- function(model) { rv <- resolve_variant(VAR[[model]], design); stopifnot(rv$p$ada$fraction == 0); rv }

if (mode == "sets") {
  model <- args[2]; stopifnot(model %in% names(VAR))
  NS <- as.integer(PR$section2_parameter_uncertainty$sampling$n_sets); stopifnot(NS == 1000L)
  from <- if (length(args) >= 3) as.integer(args[3]) else 1L; to <- if (length(args) >= 4) as.integer(args[4]) else NS
  stopifnot(from >= 0, to <= NS, to >= from)
  rv <- params_of(model); p <- rv$p; E <- check_inputs(model, p)
  lf <- if (!nzchar(test_out)) start_run_log(sprintf("day57_psets_%s_%d_%d", model, from, to), master_seed = MS, run_mode = "final", extra = list(model = model, from = from, to = to)) else NULL
  B <- 25L                                                              # 묶음: 세트 0 단독, 그 뒤 1-25, 26-50, ...(경계 고정 = 이어서 실행해도 같은 파일)
  blocks <- c(if (from == 0) list(0L), lapply(seq((max(from, 1L) - 1L) %/% B * B + 1L, to, by = B), function(b0) max(b0, from, 1L):min(b0 + B - 1L, to)))
  for (ss in blocks) {
    if (!length(ss) || max(ss) < min(ss)) next
    f <- file.path(out_dir, "psets", sprintf("day57_psets_%s_%04d_%04d.csv", model, min(ss), max(ss)))
    if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
    t0 <- Sys.time()
    R <- rbindlist(lapply(ss, function(s) {
      st <- param_set(p, E, model, s); ps <- st$p
      pop <- run_individual_population(NSUB, ps, design, "B0", MS, rv$wt_spec, jitter = TRUE, model_id = ps$model_id, tag = paste("pu", model, s))
      x <- pop$nca; stopifnot(nrow(x) == NSUB)
      ip <- individual_params(ps, pop$subj)
      sol <- solve_model(ip, CJ(id = ip$id, time = EL), design$dose_mg, model_id = ps$model_id)
      sh <- sol[, .(pct = 100 * mean(C >= LLOQ)), by = time][order(time)]
      row <- data.table(pk_model = model, set = s, n = NSUB, fail_iii_pct = 100 * mean(!crit_ok(x, "iii")), fail_i_pct = 100 * mean(!crit_ok(x, "i")),
                        cov_true_p05 = quantile(x$coverage_true, 0.05, names = FALSE), cov_true_median = median(x$coverage_true))
      for (k in seq_len(nrow(sh))) set(row, j = sprintf("above_day%d_pct", sh$time[k] + 1), value = sh$pct[k])
      cbind(row, as.data.table(st$values))
    }), fill = TRUE)
    tmp <- file.path(dirname(f), paste0("tmp_", basename(f))); fwrite(R, tmp); file.rename(tmp, f)
    msg <- sprintf("%s sets %d-%d in %s", model, min(ss), max(ss), format(Sys.time() - t0)); cat(msg, "\n"); if (!is.null(lf)) append_run_log(lf, msg)
  }
}

if (mode == "nested") {
  model <- args[2]; stopifnot(model %in% names(VAR))
  from <- if (length(args) >= 3) as.integer(args[3]) else 1L; to <- if (length(args) >= 4) as.integer(args[4]) else 200L
  stopifnot(from >= 1, to <= 200, to >= from)
  rv <- params_of(model); p <- rv$p; E <- check_inputs(model, p); wt <- rv$wt_spec; n_arm <- as.integer(design$n_per_arm); NTR <- 200L
  st <- strata_from_weight_spec(wt, as.numeric(design$stratification$split_kg$value))
  EL3 <- sort(unique(c(57, as.numeric(unlist(PR$additional_sampling_study_day)) - 1)))
  lf <- if (!nzchar(test_out)) start_run_log(sprintf("day57_nested_%s_%d_%d", model, from, to), master_seed = MS, run_mode = "final", extra = list(model = model, from = from, to = to)) else NULL
  for (s in from:to) {
    f <- file.path(out_dir, "nested", sprintf("day57_nested_%s_set%04d.csv.gz", model, s))
    if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
    t0 <- Sys.time(); ps <- param_set(p, E, model, s)$p
    lst <- lapply(seq_len(NTR), function(j) {
      sj <- with_seed(derive_seed(MS, "trial3n", model, s, j, "subj"), make_subjects(2 * n_arm, ps, wt, design$weight$sex_ratio_male$value, 0))
      sj <- with_seed(derive_seed(MS, "trial3n", model, s, j, "alloc"), assign_arms_stratified(sj, st$breaks, st$labels))
      e <- rbindlist(lapply(c("R", "T"), function(a) with_seed(derive_seed(MS, "trial3n", model, s, j, a, "eps"), draw_eps(sj[arm == a]$id, EL3, ps$sigma))))
      uid0 <- (j - 1L) * (2L * n_arm)
      list(s = sj[, `:=`(trial = j, uid = uid0 + id)], e = e[, `:=`(trial = j, uid = uid0 + id)])
    })
    subj <- rbindlist(lapply(lst, `[[`, "s")); eps <- rbindlist(lapply(lst, `[[`, "e"))
    ip <- individual_params(ps, copy(subj)[, id := uid])
    sol <- solve_model(ip, CJ(id = ip$id, time = EL3), design$dose_mg, model_id = ps$model_id)
    d <- merge(sol[, .(uid = id, planned = time, C)], subj[, .(uid, trial, arm)], by = "uid")
    d <- merge(d, eps[, .(uid, planned, eps_p, eps_a)], by = c("uid", "planned"))
    d[, y := C * (1 + eps_p) + eps_a]
    res <- d[, .(n = .N, n_true = sum(C >= LLOQ), n_obs = sum(y >= LLOQ)), by = .(trial, arm, study_day = planned + 1)][, `:=`(pk_model = model, set = s)]
    stopifnot(all(res$n == n_arm), nrow(res) == NTR * 2L * length(EL3))
    tmp <- file.path(dirname(f), paste0("tmp_", basename(f))); fwrite(res, tmp); file.rename(tmp, f)
    msg <- sprintf("%s nested set %d in %s", model, s, format(Sys.time() - t0)); cat(msg, "\n"); if (!is.null(lf)) append_run_log(lf, msg)
  }
}

if (mode == "summary") {
  max_set <- if (length(args) >= 2) as.integer(args[2]) else 1000L
  P <- rbindlist(lapply(list.files(file.path(out_dir, "psets"), pattern = "^day57_psets_.*\\.csv$", full.names = TRUE), fread), fill = TRUE)[set <= max_set]
  stopifnot(nrow(P) == uniqueN(P[, .(pk_model, set)]))
  cols <- c(grep("^above_day", names(P), value = TRUE), "fail_iii_pct", "cov_true_p05")
  q <- function(v, pr) quantile(v, pr, names = FALSE, type = 7)
  S <- rbindlist(lapply(unique(P$pk_model), function(m) { x <- P[pk_model == m & set >= 1]; ref <- P[pk_model == m & set == 0]
    rbindlist(lapply(cols, function(cl) { v <- x[[cl]]
      mcvar <- if (grepl("pct$", cl)) mean(v / 100 * (1 - v / 100) / NSUB) * 1e4 else NA_real_   # 세트 하나의 이항 몬테카를로 분산(%^2)
      data.table(pk_model = m, metric = cl, n_sets = length(v), reference_set0 = if (nrow(ref)) ref[[cl]] else NA_real_, median = median(v), lo = q(v, 0.025), hi = q(v, 0.975),
                 sd_between = sd(v), sd_mc_binomial = sqrt(mcvar)) })) }))
  fwrite(S, file.path(out_dir, sprintf("day57_param_uncertainty_summary_%d.csv", max_set)))
  NF <- list.files(file.path(out_dir, "nested"), pattern = "^day57_nested_.*\\.csv\\.gz$", full.names = TRUE)
  if (length(NF)) {
    N <- rbindlist(lapply(NF, fread))
    NS <- N[, .(n_sets = uniqueN(set), n_arm_trials = .N, median = median(n_true), p95 = q(n_true, 0.95), mean = mean(n_true), share_ge1_pct = 100 * mean(n_true >= 1),
                median_obs = median(n_obs), p95_obs = q(n_obs, 0.95)), by = .(pk_model, study_day)]
    fwrite(NS, file.path(out_dir, "day57_nested_summary.csv")); print(NS)
  }
  print(S)
}
