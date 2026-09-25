#!/usr/bin/env Rscript
# 47_sample_size.R — 표본 수 표 (추가 지시 2026-09-26 §3, config/prereg_20260926.yaml section3, D-055).
# 목표: P2(AUClast와 Cmax 모두 90% CI 80.00–125.00%) 동시 검정력. 분석 모형 M0(pooled t)와 M1(처리 + 체중 층).
# 자료 생성 모형: log Y = tau·T + beta·(log WT − 평균) + e (AUClast, Cmax), WT ~ N(75, 9) 60–90 절단, 층 = WT > 75, 층 안 1:1(홀수 잔여 교대).
#   입력(모의에서 추정, 저장 20,000명 B0): beta(log 체중 기울기), Cmax 총 log SD, AUClast–Cmax 총 상관, Cmax/AUClast SD 비. 2016 모델(1차), Model 1(민감도).
#   격자 CV = AUClast 총 개체 간 CV(sigma^2 = log(1 + CV^2)). e의 분산 = 총분산 − beta^2·Var(log WT), 공분산은 총 상관을 맞춘다.
# 모드
#   stat        입력 추정(ss_inputs.csv) → 분석식 검정력(격자, 필요 n) → 통계 모의 5,000회 × 72칸 × 입력 2벌(M0·M1 같은 시험)
#   pk <model>  PK 모델 시험 n = 117: F_down_090, F_down_095(참 AUC0-inf 비 0.90, 0.95) 5,000회, M0·M1·M2(운용특성 시드). M0는 oc_trials_be 시험 1–2,000과 대조
# 산출: results/sample_size/ (ss_inputs.csv, ss_weight_moments.csv, ss_power_analytic.csv, ss_n_needed.csv, ss_power_mc.csv, ss_pk_trials_<model>.csv.gz)
# 사용법: Rscript scripts/47_sample_size.R stat [cores]   |   nice -n 15 Rscript scripts/47_sample_size.R pk <k2016|k2020> [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("stat", "pk")) stop("사용법: Rscript scripts/47_sample_size.R <stat [cores]|pk <model> [cores]>")
mode <- args[1]
pr <- read_cfg("prereg_20260926.yaml")$section3
design <- read_cfg("trial_design.yaml"); oc <- read_cfg("oc_design.yaml")
out_dir <- Sys.getenv("DUPI_SS_OUT", proj_path(pr$out_dir)); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)   # 환경변수는 시험용
CL <- design$be$ci_level; LIMS <- as.numeric(unlist(design$be$limits))
G <- pr$grid; CVS <- as.numeric(unlist(G$auclast_cv_pct)); GMRS <- as.numeric(unlist(G$true_gmr)); NS <- as.integer(unlist(G$n_evaluable_per_arm))
SS_SEED <- 20260926L; B_MC <- as.integer(Sys.getenv("DUPI_SS_B", "5000"))                     # 시험용 축소(운영 실행은 설정하지 않음)
INPUT_VARIANTS <- c(k2016 = "base", k2020 = "struct2020")
EVAL_FRAC <- as.numeric(design$n_per_arm) / as.numeric(design$n_randomized_per_arm)   # 117/130

# ---- 체중 분포 모멘트(log WT): 전체·층 간·층 안 분산, 수치 적분 -------------------------------------------------------------
weight_moments <- function(wt, split) {
  mu <- as.numeric(wt$mean); s <- as.numeric(wt$sd); lo <- as.numeric(wt$trunc[1]); hi <- as.numeric(wt$trunc[2])
  f <- function(w) dnorm(w, mu, s)
  I <- function(a, b, g) integrate(function(w) g(w) * f(w), a, b, rel.tol = 1e-12)$value
  Z <- I(lo, hi, function(w) 1); p1 <- I(lo, split, function(w) 1) / Z; p2 <- 1 - p1
  m1 <- I(lo, split, log) / (p1 * Z); m2 <- I(split, hi, log) / (p2 * Z)
  v1 <- I(lo, split, function(w) log(w)^2) / (p1 * Z) - m1^2; v2 <- I(split, hi, function(w) log(w)^2) / (p2 * Z) - m2^2
  ma <- p1 * m1 + p2 * m2
  data.table(p_stratum2 = p2, mean_lwt = ma, var_between = p1 * (m1 - ma)^2 + p2 * (m2 - ma)^2, var_within = p1 * v1 + p2 * v2)[, var_total := var_between + var_within][]
}
WT <- weight_spec_from_design(design, "base"); SPLIT <- as.numeric(design$stratification$split_kg$value)
WM <- weight_moments(WT, SPLIT)

# ---- 입력 추정(저장 20,000명 B0) ----------------------------------------------------------------------------------------------
est_inputs <- function(m) {
  v <- INPUT_VARIANTS[[m]]
  x <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", v)))[schedule == "B0" & is.finite(AUClast) & AUClast > 0 & is.finite(Cmax) & Cmax > 0]
  la <- log(x$AUClast); lc <- log(x$Cmax); lw <- log(x$WT)
  ba <- unname(coef(lm(la ~ lw))[2]); bc <- unname(coef(lm(lc ~ lw))[2])
  data.table(input_model = m, source = sprintf("results/individual/nca_%s_20000.rds, B0", v), n_subjects = nrow(x), sd_auc = sd(la), sd_cmax = sd(lc),
             cv_auc_pct = 100 * sqrt(exp(sd(la)^2) - 1), cv_cmax_pct = 100 * sqrt(exp(sd(lc)^2) - 1), rho_total = cor(la, lc), beta_auc = ba, beta_cmax = bc,
             ratio_sd_cmax_auc = sd(lc) / sd(la), var_lwt_sample = var(lw), r2_stratum_auc = summary(lm(la ~ I(x$WT > SPLIT)))$r.squared)
}

# 칸의 분산 구조: AUClast CV와 Cmax 규칙(fixed = 모의 Cmax SD, ratio = SD 비 × AUClast SD)
cell_cov <- function(inp, cv, cmax_rule) {
  sA <- sqrt(log(1 + (cv / 100)^2)); sC <- if (cmax_rule == "fixed") inp$sd_cmax else inp$ratio_sd_cmax_auc * sA
  cov_t <- inp$rho_total * sA * sC
  vb <- WM$var_between; vt <- WM$var_total
  list(sA = sA, sC = sC, cov_t = cov_t,
       vwA = sA^2 - inp$beta_auc^2 * vb, vwC = sC^2 - inp$beta_cmax^2 * vb, cw = cov_t - inp$beta_auc * inp$beta_cmax * vb,                 # 층 안
       veA = sA^2 - inp$beta_auc^2 * vt, veC = sC^2 - inp$beta_cmax^2 * vt, ce = cov_t - inp$beta_auc * inp$beta_cmax * vt)                # 체중 보정 잔차(e)
}

# 분석식: 두 처리 추정량의 이변량 정규 직사각형(SE는 기댓값에 고정). M1: 추정량·SE 모두 층 안 분산(df 2n − 3). M0: 추정량 층 안, SE 총분산(df 2n − 2)
p2_power_analytic <- function(n, cv, gmr_a, gmr_c, inp, am, cmax_rule = "fixed") {
  cc <- cell_cov(inp, cv, cmax_rule); k <- 2 / n
  sdA <- sqrt(cc$vwA * k); sdC <- sqrt(cc$vwC * k); rho <- cc$cw / sqrt(cc$vwA * cc$vwC)
  if (am == "M1") { seA <- sdA; seC <- sdC; df <- 2 * n - 3 } else { seA <- sqrt(cc$sA^2 * k); seC <- sqrt(cc$sC^2 * k); df <- 2 * n - 2 }
  q <- qt(1 - (1 - CL) / 2, df)
  lA <- (log(LIMS[1]) + q * seA - log(gmr_a)) / sdA; uA <- (log(LIMS[2]) - q * seA - log(gmr_a)) / sdA
  lC <- (log(LIMS[1]) + q * seC - log(gmr_c)) / sdC; uC <- (log(LIMS[2]) - q * seC - log(gmr_c)) / sdC
  if (lA >= uA || lC >= uC) return(0)
  r1 <- sqrt(1 - rho^2)
  integrate(function(z) dnorm(z) * (pnorm((uC - rho * z) / r1) - pnorm((lC - rho * z) / r1)), lA, uA, rel.tol = 1e-10)$value
}
n_needed <- function(target, cv, gmr_a, gmr_c, inp, am, cmax_rule) {
  for (n in 20:600) if (p2_power_analytic(n, cv, gmr_a, gmr_c, inp, am, cmax_rule) >= target) return(n)
  NA_integer_
}
cmax_gmr <- function(g, rule, k_mech) if (rule == "equal") g else exp(k_mech * log(g))

if (mode == "stat") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
  logfile <- start_run_log("sample_size_stat", master_seed = SS_SEED, run_mode = "final", extra = list(cvs = CVS, gmrs = GMRS, ns = NS, B = B_MC))
  INP <- rbindlist(lapply(names(INPUT_VARIANTS), est_inputs))
  fwrite(INP, file.path(out_dir, "ss_inputs.csv")); fwrite(WM, file.path(out_dir, "ss_weight_moments.csv"))
  # 기전 Cmax 비(F 하향 역산): log Cmax 비 / log AUC 비
  kk <- rbindlist(lapply(names(INPUT_VARIANTS), function(m) fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", m)))[code %in% c("F_down_090", "F_down_095"),
                  .(input_model = m, code, auc_ratio, cmax_ratio, k = log(cmax_ratio) / log(auc_ratio))]))
  K_MECH <- kk[, .(k_mech = mean(k)), by = input_model]
  fwrite(kk, file.path(out_dir, "ss_cmax_gmr_mechanistic.csv"))
  # 분석식 격자
  AN <- CJ(input_model = names(INPUT_VARIANTS), cmax_cv_rule = c("fixed", "ratio"), cmax_gmr_rule = c("equal", "mechanistic"), cv = CVS, gmr = GMRS, n = NS, analysis_model = c("M0", "M1"))
  AN[, power_pct := 100 * mapply(function(m, cr, gr, cv, g, n, am) {
    inp <- INP[input_model == m]; k <- K_MECH[input_model == m, k_mech]
    p2_power_analytic(n, cv, g, cmax_gmr(g, gr, k), inp, am, cr) }, input_model, cmax_cv_rule, cmax_gmr_rule, cv, gmr, n, analysis_model)]
  fwrite(AN, file.path(out_dir, "ss_power_analytic.csv"))
  NN <- CJ(input_model = names(INPUT_VARIANTS), cmax_cv_rule = c("fixed", "ratio"), cmax_gmr_rule = c("equal", "mechanistic"), cv = CVS, gmr = GMRS, analysis_model = c("M0", "M1"), target_pct = c(90, 85))
  NN[, n_evaluable_per_arm := mapply(function(m, cr, gr, cv, g, am, tg) {
    inp <- INP[input_model == m]; k <- K_MECH[input_model == m, k_mech]
    n_needed(tg / 100, cv, g, cmax_gmr(g, gr, k), inp, am, cr) }, input_model, cmax_cv_rule, cmax_gmr_rule, cv, gmr, analysis_model, target_pct)]
  NN[, n_randomized_per_arm := as.integer(ceiling(n_evaluable_per_arm / EVAL_FRAC - 1e-9))]
  fwrite(NN, file.path(out_dir, "ss_n_needed.csv"))
  append_run_log(logfile, "analytic grid done")
  # 통계 모의: 입력 2벌 × 72칸(1차 규칙: Cmax CV 고정, Cmax GMR = AUClast GMR), 5,000회, M0·M1 같은 시험
  sim_cell <- function(m, cv, g, n) {
    inp <- INP[input_model == m]; cc <- cell_cov(inp, cv, "fixed")
    seed <- derive_seed(SS_SEED, "samplesize", m, cv, g, n)
    with_seed(seed, {
      N <- 2L * n; B <- B_MC
      wt <- rtrunc_norm(B * N, as.numeric(WT$mean), as.numeric(WT$sd), as.numeric(WT$trunc[1]), as.numeric(WT$trunc[2]))
      d <- data.table(tid = rep(seq_len(B), each = N), stratum = as.integer(wt > SPLIT), lwt = log(wt), key_ = runif(B * N))
      d[, r_ := frank(key_, ties.method = "first"), by = .(tid, stratum)]
      d[, ns_ := .N, by = .(tid, stratum)]
      d[, arm := fifelse(r_ <= ns_ %/% 2L, "R", fifelse(r_ <= 2L * (ns_ %/% 2L), "T", NA_character_))]
      d[is.na(arm), arm := { o <- order(key_); a <- character(.N); a[o] <- rep(c("R", "T"), length.out = .N); a }, by = tid]   # 홀수 층 잔여(시험당 0 또는 2명) 교대
      Se <- matrix(c(cc$veA, cc$ce, cc$ce, cc$veC), 2)
      z <- matrix(rnorm(2 * B * N), ncol = 2) %*% chol(Se)
      tt <- as.numeric(d$arm == "T"); lwc <- d$lwt - WM$mean_lwt
      d[, `:=`(yA = exp(log(g) * tt + inp$beta_auc * lwc + z[, 1]), yC = exp(log(g) * tt + inp$beta_cmax * lwc + z[, 2]))]
      stopifnot(all(d[, .N, by = .(tid, arm)]$N == n))
      res <- list()
      for (ep in c("A", "C")) {
        dd <- d[, .(tid, arm, y = if (ep == "A") yA else yC, stratum, lwt)]
        r0 <- be_m0_fast(copy(dd), "tid", CL, LIMS); r1 <- be_models_fast(copy(dd), "tid", CL, LIMS)[model == "M1"]
        res[[ep]] <- rbind(r0[, .(tid, model, pass, est, se)], r1[, .(tid, model, pass, est, se)])[, endpoint := ep]
      }
      w <- dcast(rbindlist(res), tid + model ~ endpoint, value.var = "pass")
      m_ <- m; cv_ <- cv; g_ <- g; n_ <- n                               # 인자를 지역 변수로(data.table 스코프 lint, D-022)
      w[, .(input_model = m_, cv = cv_, gmr = g_, n = n_, n_trials = .N, n_pass = sum(A %in% TRUE & C %in% TRUE), n_pass_auc = sum(A %in% TRUE), n_pass_cmax = sum(C %in% TRUE)), by = .(analysis_model = model)]
    })
  }
  cells <- CJ(m = names(INPUT_VARIANTS), cv = CVS, g = GMRS, n = NS)
  t0 <- Sys.time()
  MC <- rbindlist(parallel::mclapply(seq_len(nrow(cells)), function(i) sim_cell(cells$m[i], cells$cv[i], cells$g[i], cells$n[i]), mc.cores = cores, mc.preschedule = FALSE))
  MC[, c("power_pct", "lo", "hi") := wilson_ci(n_pass, n_trials)]
  MC <- merge(MC, AN[cmax_cv_rule == "fixed" & cmax_gmr_rule == "equal", .(input_model, cv, gmr, n, analysis_model, analytic_pct = power_pct)], by = c("input_model", "cv", "gmr", "n", "analysis_model"))
  MC[, `:=`(diff_pp = power_pct - analytic_pct, analytic_in_ci = analytic_pct >= lo & analytic_pct <= hi)]
  fwrite(MC, file.path(out_dir, "ss_power_mc.csv"))
  append_run_log(logfile, sprintf("statistical MC done in %s: %d cells, analytic within the Wilson 95%% interval in %d", format(Sys.time() - t0), nrow(MC), sum(MC$analytic_in_ci)))
  print(MC[input_model == "k2016" & gmr == 0.95, .(cv, n, analysis_model, power_pct, analytic_pct, diff_pp)])
}

if (mode == "pk") {
  model <- args[2]; stopifnot(model %in% names(INPUT_VARIANTS))
  cores <- if (length(args) >= 3) as.integer(args[3]) else 2L
  rv <- resolve_variant(INPUT_VARIANTS[[model]], design); p <- rv$p
  OC_SEED <- as.integer(oc$trials$master_seed); batch <- as.integer(oc$trials$batch)
  sc <- fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", model)))
  SCN <- c("F_down_090", "F_down_095")
  scen <- setNames(lapply(SCN, function(cd) { r <- sc[code == cd]; stopifnot(nrow(r) == 1); list(code = cd, T_multipliers = setNames(list(r$multiplier), r$mechanism)) }), SCN)
  EPS <- c("Cmax", "AUClast"); AMS <- c("M0", "M1", "M2"); NTR <- 5000L
  n_rows_trial <- length(SCN) * length(EPS) * length(AMS)
  stored <- fread(proj_path("results", "oc", sprintf("oc_trials_be_%s.csv.gz", model)))[scenario %in% SCN & endpoint %in% EPS]
  st_tr <- sort(unique(stored$trial))
  be_f <- file.path(out_dir, sprintf("ss_pk_trials_%s.csv.gz", model)); COLS <- c("trial", "scenario", "endpoint", "model", "est", "se", "pass", "n_R", "n_T")
  done <- if (file.exists(be_f)) { ex <- fread(be_f); cn <- ex[, .N, by = trial]; unique(cn[N == n_rows_trial, trial]) } else integer(0)
  if (length(done) && file.exists(be_f)) { ex <- ex[trial %in% done]; tmp <- paste0(be_f, ".rewrite.csv.gz"); fwrite(ex, tmp); file.rename(tmp, be_f); rm(ex) }
  logfile <- start_run_log(paste0("sample_size_pk_", model), master_seed = OC_SEED, run_mode = "final", extra = list(model = model, scenarios = paste(SCN, collapse = ","), trials = NTR, cores = cores,
                                                                                                         verify = sprintf("oc_trials_be_%s trials %d-%d", model, min(st_tr), max(st_tr))))
  say <- function(msg) { cat(msg, "\n"); append_run_log(logfile, msg) }
  invisible(get_model(p$model_id))
  one <- function(j) run_trial_oc_models(j, p, design, scen, OC_SEED, rv$wt_spec, model_id = p$model_id, models = AMS, endpoints = EPS)$be
  rel_eq <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a / b - 1) <= 1e-6)
  for (b0 in seq(1L, NTR, by = batch)) {
    ids <- setdiff(b0:min(b0 + batch - 1L, NTR), done)
    if (!length(ids)) next
    t0 <- Sys.time()
    res <- if (cores > 1) parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE) else lapply(ids, one)
    bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.table(x), logical(1))
    if (any(bad)) stop("실패한 시험 반복: ", paste(ids[bad], collapse = ","))
    be <- rbindlist(res); stopifnot(nrow(be) == length(ids) * n_rows_trial)
    vid <- intersect(ids, st_tr); msg <- ""
    if (length(vid)) {
      m0 <- be[model == "M0" & trial %in% vid]; q <- qt(1 - (1 - CL) / 2, m0$df)
      m0 <- m0[, .(trial, scenario, endpoint, GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se), pass, n_R, n_T)]
      mm <- merge(m0, stored[trial %in% vid], by = c("trial", "scenario", "endpoint"), all = TRUE, suffixes = c("", ".s"))
      ok <- nrow(mm) == nrow(m0) && nrow(mm) == nrow(stored[trial %in% vid]) &&
        all((rel_eq(mm$GMR, mm$GMR.s) & rel_eq(mm$CI_lower, mm$CI_lower.s) & rel_eq(mm$CI_upper, mm$CI_upper.s) & mm$pass == mm$pass.s & mm$n_R == mm$n_R.s & mm$n_T == mm$n_T.s) %in% TRUE)
      if (!ok) stop(sprintf("M0가 oc_trials_be_%s와 다릅니다(시험 %d-%d)", model, min(vid), max(vid)))
      msg <- sprintf("; M0 matched oc_trials_be for trials %d-%d", min(vid), max(vid))
    }
    be[, `:=`(est = signif(est, 8), se = signif(se, 8))]
    tmp <- tempfile(fileext = ".csv.gz", tmpdir = out_dir); fwrite(be[, ..COLS], tmp, col.names = !file.exists(be_f))
    if (file.exists(be_f)) { if (!file.append(be_f, tmp)) stop("file.append 실패"); unlink(tmp) } else file.rename(tmp, be_f)
    say(sprintf("%s trials %d-%d done in %s%s", model, min(ids), max(ids), format(Sys.time() - t0), msg))
  }
  say(sprintf("done: %s", be_f))
}
