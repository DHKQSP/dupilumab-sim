#!/usr/bin/env Rscript
# 75_m15_summary.R — [9] 그림, [10] 요약 표·핵심 결론(5줄 이내, 전제 점검)·가정 목록, 그리고 파생 표(OC 통과 확률, 짝 차이 Tango 구간,
#   기전 OC와 50% 교차(부트스트랩), AUROC(시험 단위 군집 부트스트랩)). 계획서 plans/M15_analysis_plan_20261020.md 2.1·8·10절, D-068.
# 모드
#   tables [cores]   results/m15/summary/ — oc_pass.csv, oc_consumer_risk.csv, oc_paired_diff.csv, mech_pass.csv, mech_crossing.csv, mech_crossing_diff.csv,
#                    auroc.csv, auroc_diff.csv, pauc_zero_share.csv
#   figures          results/m15/figures/ — fig1~fig7 PNG(영문) + 같은 값의 표(CSV)
#   summary          results/m15/summary/ — m15_summary_table.csv/.md(영문), m15_key_conclusions_en.md, m15_assumptions_register.csv, m15_summary_ko.md
#   pilot            results/m15/trials/pilot/pilot_report.csv — 11절 예비 실행 보고 항목(실행 시간, 층 항 제외 시험 수, 평가변수별 결측 시험 수)
# 사용법: nice -n 5 Rscript scripts/75_m15_summary.R <tables|figures|summary> [cores]
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("tables", "figures", "summary", "pilot")) stop("사용법: Rscript scripts/75_m15_summary.R <tables|figures|summary|pilot> [cores]")
mode <- args[1]; cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
RD <- proj_path(CFG$output_dir); SD <- file.path(RD, "summary"); FD <- file.path(RD, "figures")
dir.create(SD, showWarnings = FALSE, recursive = TRUE); dir.create(FD, showWarnings = FALSE, recursive = TRUE)
NBASE <- as.integer(CFG$mechanisms$base_n_randomized_per_arm); B_BOOT <- 2000L
PRIMARY <- CFG$models$primary
COND <- m15_conditional_label(RD, PRIMARY)                          # 조건부 표시(계획서 6.2절): 결과 표마다 열로, 그림마다 설명으로
combo_lab <- function(dd) { p <- strsplit(basename(dd), "_")[[1]]; list(pk_model = p[1], population = paste(p[2:(length(p) - 1)], collapse = "_"), dose_mg = as.numeric(p[length(p)])) }
read_dir <- function(dd, pat) { f <- list.files(dd, pattern = pat, full.names = TRUE); if (!length(f)) return(NULL); rbindlist(lapply(f, fread), fill = TRUE) }
# 완전성 점검: 키마다 한 행, 등록된 시험 번호가 빠짐없이(중복·누락이면 중단; 묶음 파일을 지우고 그 범위만 다시 실행)
check_trials <- function(x, keys, expect, what) {
  if (anyDuplicated(x, by = c("trial", keys))) stop(what, ": 중복 행(같은 시험이 두 번)")
  got <- sort(unique(x$trial)); if (!identical(as.integer(got), as.integer(expect))) stop(what, sprintf(": 시험 %d개(기대 %d개), 누락 예: %s", length(got), length(expect), paste(head(setdiff(expect, got), 5), collapse = ",")))
  cn <- x[, .N, by = keys]; if (any(cn$N != length(expect))) stop(what, sprintf(": 칸 %d개에서 시험 수가 %d개가 아님", sum(cn$N != length(expect)), length(expect)))   # 칸마다 모든 시험
  invisible(TRUE)
}
NT_MAIN <- as.integer(CFG$oc$n_trials); NT_EXT <- as.integer(CFG$oc$extension$n_trials)
bit <- function(mask, v) bitwAnd(as.integer(mask), as.integer(2^(v - 1))) > 0
wil <- function(k, n) { w <- wilson_ci(k, n); list(pct = w$est, lo = w$lo, hi = w$hi) }

# ================================================================ 예비 실행 보고(11절) ================================================================
if (mode == "pilot") {
  PR <- list(c("k2020_S1_300", 1L, 500L), c("k2020_S2_200", 1L, 50L))
  out <- rbindlist(lapply(PR, function(z) {
    oc <- read_dir(file.path(RD, "trials", "oc", z[1]), "^oc_\\d+_\\d+\\.csv\\.gz$")[trial >= as.integer(z[2]) & trial <= as.integer(z[3])]
    me <- read_dir(file.path(RD, "trials", "mech", z[1]), "^mech_\\d+_\\d+\\.csv\\.gz$")[trial >= as.integer(z[2]) & trial <= as.integer(z[3])]
    a <- oc[model == "M1", .(item = "studies with the stratum term dropped (M1, AUC0-last)", value = sum(m1_dropped > 0), n_studies = uniqueN(trial)), by = .(boundary, schedule, n_rand)]
    ec <- grep("^est_", names(me), value = TRUE)
    b <- rbindlist(lapply(ec, function(cc) me[, .(item = sprintf("studies with %s missing at 130 per arm (fewer than 2 analysable subjects in an arm)", sub("^est_", "", cc)),
                                                value = uniqueN(trial[!is.finite(get(cc))]), n_studies = uniqueN(trial)), by = .(schedule)]))   # 판 하나라도 결측인 시험 수[, `:=`(boundary = "fixed75", n_rand = NBASE)]
    rbind(a, b, fill = TRUE)[, combo := z[1]]
  }), fill = TRUE)
  pl <- readLines(file.path(RD, "run_progress.log")); rs <- grep("pilot restarted", pl); if (length(rs)) pl <- pl[(max(rs) + 1):length(pl)]   # 시드 수정 뒤 예비 실행만(D-068 이탈 1)
  dn <- grep("done scripts/73_m15_trials.R (oc|mech) k2020 (S1 300 1 500|S2 200 1 50)", pl, value = TRUE)
  tm <- data.table(combo = NA_character_, item = sub(".*done scripts/73_m15_trials.R ", "run time: ", dn), value = as.numeric(sub(".* in ([0-9]+)s$", "\\1", dn)), n_studies = NA_integer_)
  out <- rbind(out, tm, fill = TRUE); setcolorder(out, c("combo", "item")); fwrite(out, file.path(RD, "trials", "pilot", "pilot_report.csv")); print(out[value > 0 | grepl("run time", item)])
}

# ================================================================ 파생 표 ================================================================
if (mode == "tables") {
  lf <- start_run_log("m15_summary_tables", master_seed = MS, run_mode = "final", extra = list(cores = cores, B = B_BOOT))
  # ---------------------------------------------------------------- [3] OC 통과 확률 ----------------------------------------------------------------
  oc_dirs <- list.dirs(file.path(RD, "trials", "oc"), recursive = FALSE)
  OCX <- list()                                                    # 시험별 통과(짝 차이용): 키 = 묶음
  OCP <- rbindlist(lapply(oc_dirs, function(dd) {
    x <- read_dir(dd, "^oc_\\d+_\\d+\\.csv\\.gz$"); if (is.null(x)) return(NULL); v <- fread(file.path(dd, "versions.csv")); lb <- combo_lab(dd)
    check_trials(x, c("boundary", "schedule", "n_rand", "model"), seq_len(NT_MAIN), basename(dd))
    long <- rbindlist(lapply(seq_len(nrow(v)), function(i) x[, .(trial, boundary, schedule, n_rand, model, g = v$g[i], pass = bit(mask, v$vi[i]), pass_auc = bit(mask_auc, v$vi[i]),
                                                                pass_cmax = bit(mask_cmax, v$vi[i]), m1_dropped, n_stratum2)]))
    OCX[[basename(dd)]] <<- long[n_rand == NBASE][, `:=`(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg)]
    long[, .(n_trials = .N, n_pass = sum(pass), n_pass_auc = sum(pass_auc), n_pass_cmax = sum(pass_cmax), m1_dropped_studies = sum(m1_dropped > 0),
             n_stratum2_median = as.numeric(median(n_stratum2))), by = .(boundary, schedule, n_rand, model, g)][, `:=`(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg)]
  }))
  OCP[model == "M0", m1_dropped_studies := NA_integer_]             # M0에는 층 항이 없다(73은 M1 값을 두 모형 행에 붙인다)
  OCP[, c("pass_pct", "pass_lo", "pass_hi") := wil(n_pass, n_trials)]
  OCP[, c("auc_only_pct", "auc_only_lo", "auc_only_hi") := wil(n_pass_auc, n_trials)]
  OCP[, c("cmax_only_pct", "cmax_only_lo", "cmax_only_hi") := wil(n_pass_cmax, n_trials)]
  setcolorder(OCP, c("pk_model", "population", "dose_mg", "boundary", "schedule", "model", "n_rand", "g")); m15_fwrite_cond(OCP, file.path(SD, "oc_pass.csv"), COND)
  # 소비자 위험(10,000회 = 본 실행 1~1,000 + 확장 1,001~10,000; k2020)
  CR <- rbindlist(lapply(list.dirs(file.path(RD, "trials", "ocext"), recursive = FALSE), function(dd) {
    x <- read_dir(dd, "^ocext_\\d+_\\d+\\.csv\\.gz$"); if (is.null(x)) return(NULL); v <- fread(file.path(dd, "versions.csv")); lb <- combo_lab(dd)
    check_trials(x, c("boundary", "schedule", "n_rand", "model"), (NT_MAIN + 1L):NT_EXT, paste("ocext", basename(dd)))
    rbindlist(lapply(seq_len(nrow(v)), function(i) x[, .(n_trials = .N, n_pass = sum(bit(mask, v$vi[i])), n_pass_auc = sum(bit(mask_auc, v$vi[i])), n_pass_cmax = sum(bit(mask_cmax, v$vi[i]))),
                                                     by = .(boundary, schedule, n_rand, model)][, g := v$g[i]]))[, `:=`(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg)]
  }))
  if (nrow(CR)) {
    m <- merge(CR, OCP[, .(pk_model, population, dose_mg, boundary, schedule, n_rand, model, g, m_n = n_trials, m_p = n_pass, m_a = n_pass_auc, m_c = n_pass_cmax)],
               by = c("pk_model", "population", "dose_mg", "boundary", "schedule", "n_rand", "model", "g"))
    m[, `:=`(n_trials = n_trials + m_n, n_pass = n_pass + m_p, n_pass_auc = n_pass_auc + m_a, n_pass_cmax = n_pass_cmax + m_c)][, c("m_n", "m_p", "m_a", "m_c") := NULL]
    stopifnot(all(m$n_trials == NT_EXT))
    m[, c("pass_pct", "pass_lo", "pass_hi") := wil(n_pass, n_trials)]; m[, c("auc_only_pct", "auc_only_lo", "auc_only_hi") := wil(n_pass_auc, n_trials)]
    m[, c("cmax_only_pct", "cmax_only_lo", "cmax_only_hi") := wil(n_pass_cmax, n_trials)]
    m15_fwrite_cond(m, file.path(SD, "oc_consumer_risk.csv"), COND)
  }
  # 균등 모양 민감도
  SHP <- rbindlist(lapply(list.dirs(file.path(RD, "trials", "ocshape"), recursive = FALSE), function(dd) {
    x <- read_dir(dd, "^ocshape_\\d+_\\d+\\.csv\\.gz$"); if (is.null(x)) return(NULL); v <- fread(file.path(dd, "versions.csv")); lb <- combo_lab(dd)
    check_trials(x, c("boundary", "schedule", "n_rand", "model"), seq_len(as.integer(CFG$oc$weight_shape_sensitivity$n_trials)), paste("ocshape", basename(dd)))
    rbindlist(lapply(seq_len(nrow(v)), function(i) x[, .(n_trials = .N, n_pass = sum(bit(mask, v$vi[i]))), by = .(boundary, schedule, n_rand, model)][, g := v$g[i]]))[,
      `:=`(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg)] }))
  if (nrow(SHP)) { SHP[, c("pass_pct", "pass_lo", "pass_hi") := wil(n_pass, n_trials)]; m15_fwrite_cond(SHP, file.path(SD, "oc_weight_shape.csv"), COND) }
  # 짝 차이(Tango): 130명, g = 0.80, 0.95, 1.00, 1.25
  X <- rbindlist(OCX)[abs(g - 0.80) < 1e-9 | abs(g - 0.95) < 1e-9 | abs(g - 1) < 1e-9 | abs(g - 1.25) < 1e-9]
  pdiff <- function(a, b, by_cols, ct_lab, lab_a, lab_b) {
    z <- merge(a, b, by = c(by_cols, "trial"), suffixes = c("_a", "_b")); if (!nrow(z)) return(NULL)
    z[, { bb <- sum(pass_a & !pass_b); cc <- sum(!pass_a & pass_b); r <- m15_tango_ci(bb, cc, .N)
      .(contrast = ct_lab, design_a = lab_a, design_b = lab_b, n_trials = .N, pass_a_pct = 100 * mean(pass_a), pass_b_pct = 100 * mean(pass_b),
        diff_pp = 100 * r[["diff"]], lo_pp = 100 * r[["lo"]], hi_pp = 100 * r[["hi"]], discordant_a_only = bb, discordant_b_only = cc) }, by = by_cols]
  }
  K <- X[, .(trial, pk_model, population, dose_mg, boundary, schedule, model, g, pass)]
  PD <- rbind(
    pdiff(K[model == "M1"][, model := NULL], K[model == "M0"][, model := NULL], c("pk_model", "population", "dose_mg", "boundary", "schedule", "g"), "stratum term", "M1", "M0"),
    pdiff(K[boundary == "median"][, boundary := NULL], K[boundary == "fixed75"][, boundary := NULL], c("pk_model", "population", "dose_mg", "schedule", "model", "g"), "boundary", "median", "fixed75"),
    pdiff(K[schedule == "RM36_50"][, schedule := NULL], K[schedule == "B0"][, schedule := NULL], c("pk_model", "population", "dose_mg", "boundary", "model", "g"), "schedule", "RM36_50", "B0"),
    pdiff(K[schedule == "ADD64_71"][, schedule := NULL], K[schedule == "B0"][, schedule := NULL], c("pk_model", "population", "dose_mg", "boundary", "model", "g"), "schedule", "ADD64_71", "B0"),
    pdiff(K[population == "S2"][, population := NULL], K[population == "S1"][, population := NULL], c("pk_model", "dose_mg", "boundary", "schedule", "model", "g"), "population", "S2", "S1"),
    pdiff(K[dose_mg == 200][, dose_mg := NULL], K[dose_mg == 300][, dose_mg := NULL], c("pk_model", "population", "boundary", "schedule", "model", "g"), "dose", "200", "300"), fill = TRUE)
  rd <- as.numeric(sub(".*at least ([0-9.]+) percentage points.*", "\\1", CFG$interpretation$design_difference)); stopifnot(is.finite(rd))
  PD[, changes_oc := (lo_pp > 0 | hi_pp < 0) & abs(diff_pp) >= rd]
  PD[, n_rand := NBASE]                                              # 짝 차이는 군당 130명, g 0.80·0.95·1.00·1.25에서만(위 OCX, X)
  m15_fwrite_cond(PD, file.path(SD, "oc_paired_diff.csv"), COND)
  append_run_log(lf, sprintf("OC tables: %d cells, %d paired differences", nrow(OCP), nrow(PD)))

  # ---------------------------------------------------------------- [4] 기전 OC와 50% 교차 ----------------------------------------------------------------
  TG <- fread(file.path(RD, "truth", "truth_gmr.csv"))
  mech_dirs <- list.dirs(file.path(RD, "trials", "mech"), recursive = FALSE)
  rd_mech <- function(kind, dd) { d2 <- file.path(RD, "trials", kind, basename(dd)); if (!dir.exists(d2)) return(NULL)
    r <- list(cnt = read_dir(d2, sprintf("^cnt_%s_\\d+_\\d+\\.csv\\.gz$", kind)), base = read_dir(d2, sprintf("^%s_\\d+_\\d+\\.csv\\.gz$", kind)), v = fread(file.path(d2, "versions.csv")))
    nt_ <- if (kind == "mechsupp") as.integer(CFG$mechanisms$supplementary$n_trials) else as.integer(CFG$mechanisms$n_trials)
    if (!is.null(r$base)) check_trials(r$base, c("schedule", "vi"), seq_len(nt_), paste(kind, basename(dd)))
    if (!is.null(r$cnt)) stopifnot(all(r$cnt[, .(n = sum(n_trials)), by = .(schedule, vi, n_rand, model)]$n == nt_))
    r }
  MP <- rbindlist(lapply(c("mech", "mechsupp"), function(kind) rbindlist(lapply(list.dirs(file.path(RD, "trials", kind), recursive = FALSE), function(dd) {
    r <- rd_mech(kind, dd); if (is.null(r$cnt)) return(NULL); lb <- combo_lab(dd)
    merge(r$cnt[, .(n_trials = sum(n_trials), n_pass = sum(n_pass), n_pass_auc = sum(n_pass_auc), n_pass_cmax = sum(n_pass_cmax)), by = .(schedule, vi, n_rand, model)], r$v, by = "vi")[,
      `:=`(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg, set = kind)] }))))
  MP[, c("pass_pct", "pass_lo", "pass_hi") := wil(n_pass, n_trials)]
  setcolorder(MP, c("pk_model", "population", "dose_mg", "set", "schedule", "mechanism", "multiplier", "n_rand", "model")); m15_fwrite_cond(MP, file.path(SD, "mech_pass.csv"), COND)
  # 시험별 통과 행렬(130명, M1): 기전 본 실행 + 보충(같은 시험 j, 항등은 같다)
  BIDX <- m15_boot_index(as.integer(CFG$mechanisms$n_trials), B_BOOT, derive_seed(MS, "crossing_boot"))
  cross_one <- function(mult, p) {                                  # 1.00에서 바깥쪽 첫 50% 미만 교차(log 배율 선형 보간)
    o <- order(abs(log(mult))); mult <- mult[o]; p <- p[o]
    if (p[1] < 0.5) return(1)
    for (k in 2:length(mult)) if (p[k] < 0.5) return(exp(log(mult[k - 1]) + (0.5 - p[k - 1]) / (p[k] - p[k - 1]) * (log(mult[k]) - log(mult[k - 1]))))
    NA_real_
  }
  CRS <- list(); CRB <- list()
  for (dd in mech_dirs) {
    lb <- combo_lab(dd); r <- rd_mech("mech", dd); if (is.null(r$base)) next
    rs <- rd_mech("mechsupp", dd)
    for (sn in unique(r$base$schedule)) {
      b <- r$base[schedule == sn, .(trial, vi, pass = pass_M1)][r$v, on = "vi"]
      if (!is.null(rs) && !is.null(rs$base) && sn %in% rs$base$schedule) {
        bs <- rs$base[schedule == sn & vi > 1, .(trial, vi, pass = pass_M1)][rs$v, on = "vi", nomatch = NULL]
        b <- rbind(b, bs[, names(b), with = FALSE])
      }
      for (mech in setdiff(unique(b$mechanism), "identity")) for (side in c("down", "up")) {
        bb <- b[mechanism %in% c(mech, "identity") & (if (side == "down") multiplier <= 1 else multiplier >= 1)]
        W <- dcast(bb, trial ~ multiplier, value.var = "pass", fun.aggregate = function(z) as.integer(z[1]))
        mult <- as.numeric(names(W)[-1]); P <- as.matrix(W[, -1]); nt <- nrow(P)
        est <- cross_one(mult, colMeans(P))
        bt <- vapply(BIDX, function(i) { r_ <- match(i, W$trial); r_ <- r_[!is.na(r_)]; cross_one(mult, colMeans(P[r_, , drop = FALSE])) }, 1)
        key <- paste(lb$pk_model, lb$population, lb$dose_mg, sn, mech, side, sep = "|")
        CRB[[key]] <- bt
        # 미도달 재표집은 격자 끝 너머로 검열(버리지 않음): 아래쪽은 0, 위쪽은 Inf로 두고 백분위수를 구한다
        btc <- ifelse(is.finite(bt), bt, if (side == "down") 0 else Inf); qb <- quantile(btc, c(0.025, 0.975), names = FALSE, type = 7)
        nr_ <- sum(!is.finite(bt)); beyond_outer <- nr_ >= 0.025 * length(bt)       # 바깥쪽 한계: 미도달 재표집이 2.5% 이상이면 격자 끝 너머
        beyond <- function(q) beyond_outer && ((side == "down" && identical(q, qb[1])) || (side == "up" && identical(q, qb[2])))
        # 유도 참 GMR(교차점에서 log 배율 선형 보간)
        tg <- TG[pk_model == lb$pk_model & population == lb$population & dose_mg == lb$dose_mg & (mechanism == mech | mechanism == "identity") &
                 metric %in% c("AUCinf", "Cmax", paste0("AUClast_dense_", sn))][(if (side == "down") multiplier <= 1 else multiplier >= 1)]
        ig <- function(mt) { z <- tg[metric == mt][order(multiplier)]; if (!is.finite(est) || nrow(z) < 2 || est < min(z$multiplier) || est > max(z$multiplier)) return(NA_real_)
          exp(approx(log(z$multiplier), log(z$gmr), xout = log(est))$y) }
        CRS[[key]] <- data.table(pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg, schedule = sn, mechanism = mech, side = side, model = "M1", n_rand = NBASE,
                                 grid_end = if (side == "down") min(mult) else max(mult), p_at_grid_end = colMeans(P)[which.max(abs(log(mult)))], reached = is.finite(est),
                                 crossing_multiplier = est, crossing_pct = 100 * (est - 1),
                                 lo_pct = if (!is.finite(est) || beyond(qb[1])) NA_real_ else 100 * (qb[1] - 1), hi_pct = if (!is.finite(est) || beyond(qb[2])) NA_real_ else 100 * (qb[2] - 1),
                                 outer_bound_beyond_grid = beyond(if (side == "down") qb[1] else qb[2]), boot_not_reached_pct = 100 * mean(!is.finite(bt)),
                                 induced_gmr_AUClast = ig(paste0("AUClast_dense_", sn)), induced_gmr_AUCinf = ig("AUCinf"), induced_gmr_Cmax = ig("Cmax"), n_trials = nt)
      }
    }
  }
  CRT <- rbindlist(CRS); m15_fwrite_cond(CRT, file.path(SD, "mech_crossing.csv"), COND)
  # 교차 차이(같은 재표집 색인): S2 - S1, 200 - 300, 일정 - B0 (%p 단위 차이; 어느 쪽이든 미도달인 재표집은 제외, 비율 보고)
  CRT[, ckey := paste(pk_model, population, dose_mg, schedule, mechanism, side, sep = "|")]
  cdiff <- function(key_a, key_b, ct_lab) { a <- CRB[[key_a]]; b <- CRB[[key_b]]; if (is.null(a) || is.null(b)) return(NULL)
    sd_ <- if (grepl("\\|down$", key_a)) -Inf else Inf                                      # 미도달 = 격자 끝 너머(아래쪽 -Inf, 위쪽 Inf)로 검열: 그 쪽 한계는 정할 수 없음(NA)
    ac <- ifelse(is.finite(a), a, sd_); bc_ <- ifelse(is.finite(b), b, sd_)
    d <- 100 * (ac - bc_); d[!is.finite(a) & !is.finite(b)] <- NA_real_                       # 둘 다 미도달이면 차이를 정할 수 없다
    ea <- CRT$crossing_multiplier[match(key_a, CRT$ckey)]; eb <- CRT$crossing_multiplier[match(key_b, CRT$ckey)]
    na_ <- 100 * mean(is.na(d)); q_ <- if (na_ <= 2.5) quantile(d, c(0.025, 0.975), names = FALSE, na.rm = TRUE) else c(NA_real_, NA_real_)
    q_[!is.finite(q_)] <- NA_real_                                                            # 무한대에 놓인 한계 = 정할 수 없음
    data.table(contrast = ct_lab, a = key_a, b = key_b, diff_pp = 100 * (ea - eb), lo_pp = if (is.finite(ea) && is.finite(eb)) q_[1] else NA_real_,
               hi_pp = if (is.finite(ea) && is.finite(eb)) q_[2] else NA_real_, boot_not_reached_pct = 100 * mean(!is.finite(a) | !is.finite(b)), boot_both_not_reached_pct = na_) }
  CD <- rbindlist(lapply(seq_len(nrow(CRT)), function(i) { r <- CRT[i]; k <- function(pop = r$population, dose = r$dose_mg, sn = r$schedule) paste(r$pk_model, pop, dose, sn, r$mechanism, r$side, sep = "|")
    rbind(if (r$population == "S2") cdiff(k(), k(pop = "S1"), "population S2 - S1"), if (r$dose_mg == 200) cdiff(k(), k(dose = 300), "dose 200 - 300"),
          if (r$schedule != "B0") cdiff(k(), k(sn = "B0"), sprintf("schedule %s - B0", r$schedule))) }))
  m15_fwrite_cond(CD, file.path(SD, "mech_crossing_diff.csv"), COND)
  append_run_log(lf, sprintf("mechanism tables: %d crossings", nrow(CRT)))

  # ---------------------------------------------------------------- [5] AUROC ----------------------------------------------------------------
  AIDX <- m15_boot_index(as.integer(CFG$mechanisms$n_trials), B_BOOT, derive_seed(MS, "auroc_boot"))
  PA <- c("pAUC_0_14", "pAUC_14_28", "pAUC_28_56", "pAUC_28_56_ext")
  stat_cols <- function(b) {                                         # 통계량 열: 주(AUCinf_A M1, 부분 AUC 산술 평균비) + 민감도
    out <- b[, .(trial, schedule, vi, AUCinf_A = abs(est_AUCinf_A), AUCinf_B = abs(est_AUCinf_B))]
    for (e in PA) {
      mr <- b[[paste0("mean_R_", e)]]; mt <- b[[paste0("mean_T_", e)]]; ok <- b[[paste0("n_R_", e)]] >= 2 & b[[paste0("n_T_", e)]] >= 2 & mr > 0 & mt > 0
      set(out, j = paste0(e, "_meanratio"), value = fifelse(ok, abs(log(mt / mr)), NA_real_))
      set(out, j = paste0(e, "_logc"), value = abs(b[[paste0("est_", e)]]))
    }
    out
  }
  SIZES <- list(`5` = c(0.95, 1.05), `10` = c(0.90, 1.10), `15` = c(0.85, 1.15))
  au_jobs <- rbindlist(lapply(mech_dirs, function(dd) { lb <- combo_lab(dd); data.table(dir = dd, pk_model = lb$pk_model, population = lb$population, dose_mg = lb$dose_mg) }))
  AU <- parallel::mclapply(seq_len(nrow(au_jobs)), function(jj) {
    dd <- au_jobs$dir[jj]; r <- rd_mech("mech", dd); if (is.null(r$base)) return(NULL)
    s <- stat_cols(r$base); v <- r$v; ST <- setdiff(names(s), c("trial", "schedule", "vi"))
    zs <- r$base[vi == 1, c("schedule", grep("^(zero_R|n_R|zero_T|n_T)_", names(r$base), value = TRUE)), with = FALSE]
    ZS <- rbindlist(lapply(PA, function(e) zs[, .(endpoint = e, zero_share_R = sum(get(paste0("zero_R_", e)), na.rm = TRUE) / sum(get(paste0("n_R_", e)), na.rm = TRUE),
                                                  zero_share_T_identity = sum(get(paste0("zero_T_", e)), na.rm = TRUE) / sum(get(paste0("n_T_", e)), na.rm = TRUE)), by = schedule]))
    # 양성(차이를 넣은 판)의 시험군 0 비율: 판마다
    ZV <- rbindlist(lapply(PA, function(e) r$base[, .(endpoint = e, zT = sum(get(paste0("zero_T_", e)), na.rm = TRUE), nT = sum(get(paste0("n_T_", e)), na.rm = TRUE)), by = .(schedule, vi)]))
    rows <- list(); draws <- list()
    for (sn in unique(s$schedule)) {
      w <- s[schedule == sn]; trials <- sort(unique(w$trial))
      idx_all <- lapply(AIDX, function(i) { pi_ <- match(i, trials); pi_[!is.na(pi_)] })   # 재표집 색인을 시험 위치로(일정마다 한 번)
      col_at <- function(x, st_) { z <- w[vi == x]; z[[st_]][match(trials, z$trial)] }
      for (st in ST) {
        neg <- col_at(1L, st)
        for (mech in unique(v$mechanism[v$mechanism != "identity"])) for (sz in names(SIZES)) for (dir_ in c("down", "up", "pooled")) {
          ms <- SIZES[[sz]]; ms <- if (dir_ == "down") ms[1] else if (dir_ == "up") ms[2] else ms
          vis <- v[mechanism == mech & abs(multiplier - ms[1]) < 1e-9 | mechanism == mech & abs(multiplier - ms[length(ms)]) < 1e-9, vi]
          pos <- matrix(sapply(vis, function(x) col_at(x, st)), ncol = length(vis))
          keep <- is.finite(neg) & apply(is.finite(pos), 1, all)
          pk <- pos[keep, , drop = FALSE]; nk <- neg[keep]
          ep_ <- sub("_(meanratio|logc)$", "", st); zT_ <- if (ep_ %in% PA) { zz <- ZV[schedule == sn & endpoint == ep_ & vi %in% vis]; sum(zz$zT) / sum(zz$nT) } else NA_real_
          if (sum(keep) < 10) { rows[[length(rows) + 1]] <- data.table(schedule = sn, statistic = st, mechanism = mech, size_pct = as.numeric(sz), direction = dir_, n_studies = sum(keep), auroc = NA_real_,
                                                                        zero_share_T_pos = zT_); next }
          a <- .m15_mw(as.vector(pk), nk); dl <- m15_auroc(as.vector(pk), nk)
          # 같은 색인(전체 시험 기준)에서 제외된 시험을 뺀다
          idx <- if (all(keep)) idx_all else { map_ <- cumsum(keep); lapply(idx_all, function(pi_) map_[pi_[keep[pi_]]]) }
          bt <- m15_auroc_draws_w(pk, nk, idx)                                               # 가중 Mann-Whitney(재표집 자료 직접 계산과 같은 값, 단위 시험)
          key <- paste(sn, st, mech, sz, dir_, sep = "|"); draws[[key]] <- bt
          rows[[length(rows) + 1]] <- data.table(schedule = sn, statistic = st, mechanism = mech, size_pct = as.numeric(sz), direction = dir_, n_studies = sum(keep), auroc = a,
                                                 lo = quantile(bt, 0.025, names = FALSE, na.rm = TRUE), hi = quantile(bt, 0.975, names = FALSE, na.rm = TRUE),
                                                 delong_lo = dl$lo, delong_hi = dl$hi, zero_share_T_pos = zT_)
        }
      }
    }
    list(rows = rbindlist(rows, fill = TRUE)[, `:=`(pk_model = au_jobs$pk_model[jj], population = au_jobs$population[jj], dose_mg = au_jobs$dose_mg[jj])],
         zs = ZS[, `:=`(pk_model = au_jobs$pk_model[jj], population = au_jobs$population[jj], dose_mg = au_jobs$dose_mg[jj])], draws = draws)
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(AU, function(x) is.null(x) || inherits(x, "try-error"), TRUE); if (length(AU) != nrow(au_jobs) || any(bad)) stop("AUROC 실패 또는 작업자 결과 없음: ", which(bad))
  AUR <- rbindlist(lapply(AU, `[[`, "rows"), fill = TRUE); ZS <- rbindlist(lapply(AU, `[[`, "zs"))
  ZS[, not_interpretable := zero_share_R > 0.5]; m15_fwrite_cond(ZS, file.path(SD, "pauc_zero_share.csv"), COND)
  AUR[, endpoint := sub("_(meanratio|logc)$", "", statistic)]
  AUR <- merge(AUR, ZS[, .(pk_model, population, dose_mg, schedule, endpoint, zero_share_R, not_interpretable)], by = c("pk_model", "population", "dose_mg", "schedule", "endpoint"), all.x = TRUE)
  AUR[, role := fifelse(statistic %in% c("AUCinf_A", paste0(c("pAUC_0_14", "pAUC_14_28", "pAUC_28_56"), "_meanratio")), "primary", "sensitivity")]
  m15_fwrite_cond(AUR, file.path(SD, "auroc.csv"), COND)
  # 차이(같은 부트스트랩 색인): 일정 - B0, 부분 AUC - AUCinf_A, S2 - S1, 200 - 300 (양방향 합침)
  DR <- list(); names(AU) <- paste(au_jobs$pk_model, au_jobs$population, au_jobs$dose_mg, sep = "|")
  est_of <- function(cb, key) { p <- strsplit(cb, "\\|")[[1]]; k <- strsplit(key, "\\|")[[1]]
    AUR[pk_model == p[1] & population == p[2] & dose_mg == as.numeric(p[3]) & schedule == k[1] & statistic == k[2] & mechanism == k[3] & size_pct == as.numeric(k[4]) & direction == k[5], auroc] }
  ni_of <- function(cb, key) { p <- strsplit(cb, "\\|")[[1]]; k <- strsplit(key, "\\|")[[1]]
    z <- AUR[pk_model == p[1] & population == p[2] & dose_mg == as.numeric(p[3]) & schedule == k[1] & statistic == k[2] & mechanism == k[3] & size_pct == as.numeric(k[4]) & direction == k[5], not_interpretable]
    isTRUE(z[1]) }
  adiff <- function(ca, ka, cbb, kb, ct_lab) { a <- AU[[ca]]$draws[[ka]]; b <- AU[[cbb]]$draws[[kb]]; if (is.null(a) || is.null(b)) return(NULL); d <- a - b
    data.table(contrast = ct_lab, a = paste(ca, ka, sep = "|"), b = paste(cbb, kb, sep = "|"), diff = est_of(ca, ka) - est_of(cbb, kb),
               lo = quantile(d, 0.025, names = FALSE, na.rm = TRUE), hi = quantile(d, 0.975, names = FALSE, na.rm = TRUE),
               not_interpretable_a = ni_of(ca, ka), not_interpretable_b = ni_of(cbb, kb))[, not_interpretable := not_interpretable_a | not_interpretable_b] }
  for (cb in names(AU)) for (key in names(AU[[cb]]$draws)) {
    k <- strsplit(key, "\\|")[[1]]; if (k[5] != "pooled") next; p <- strsplit(cb, "\\|")[[1]]
    if (k[1] != "B0") DR[[length(DR) + 1]] <- adiff(cb, key, cb, paste(c("B0", k[-1]), collapse = "|"), sprintf("schedule %s - B0", k[1]))
    if (k[2] != "AUCinf_A") DR[[length(DR) + 1]] <- adiff(cb, key, cb, paste(c(k[1], "AUCinf_A", k[3:5]), collapse = "|"), sprintf("endpoint %s - AUCinf_A", k[2]))
    if (p[2] == "S2") DR[[length(DR) + 1]] <- adiff(cb, key, paste(p[1], "S1", p[3], sep = "|"), key, "population S2 - S1")
    if (p[3] == "200") DR[[length(DR) + 1]] <- adiff(cb, key, paste(p[1], p[2], "300", sep = "|"), key, "dose 200 - 300")
  }
  m15_fwrite_cond(rbindlist(DR), file.path(SD, "auroc_diff.csv"), COND)
  append_run_log(lf, sprintf("AUROC tables: %d cells, %d differences", nrow(AUR), length(DR)))
  append_run_log(lf, "done")
}

# ================================================================ 그림 ================================================================
fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
POPL <- c(S1 = "Global healthy-volunteer population", S2 = "Korean healthy-volunteer population")
SCHL <- c(B0 = "Confirmed schedule", RM36_50 = "Reduced (Days 36 to 50 removed)", ADD64_71 = "Extended (Days 64 and 71 added)")
MODL <- c(k2020 = "Kovalenko 2020 Model 1 (primary)", k2016 = "Kovalenko 2016 (sensitivity)")
MODS <- c(k2020 = "Kovalenko 2020 Model 1\n(primary)", k2016 = "Kovalenko 2016\n(sensitivity)")   # 그림 띠 표기(두 줄, 잘림 방지)
SCHS <- c(B0 = "Confirmed", RM36_50 = "Reduced", ADD64_71 = "Extended")
COLS3 <- c(VIZ$s1, VIZ$s2, VIZ$s3)
save_fig <- function(g, name, w, h, tab, lim200 = FALSE) {
  cap <- c(if (nzchar(COND)) COND, if (lim200) paste("200 mg panels:", M15_LIM200)); if (length(cap)) g <- g + labs(caption = paste(cap, collapse = "\n"))
  ggsave(file.path(FD, paste0(name, ".png")), g, width = w, height = h, dpi = 150, bg = "white"); m15_fwrite_cond(tab, file.path(FD, paste0(name, ".csv")), COND) }
if (mode == "figures") {
  lf <- start_run_log("m15_summary_figures", master_seed = MS, run_mode = "final")
  # fig1: [1] 예측 대 관측(판정 단위, 고정 체중 가정)
  U <- fread(file.path(RD, "validation", "validation_units.csv"))[wmode %in% c("primary", "wunc") & metric %in% c("AUClast_mean", "AUClast_cv_within", "Cmax_mean") & is.finite(observed)]
  U[, `:=`(model_lab = factor(MODS[pk_model], MODS), weight_lab = fifelse(wmode == "primary", "Fixed weight assumption", "Uncertain-mean weight assumption"),
           metric_lab = c(AUClast_mean = "AUC0-last mean (mg x day/L)", AUClast_cv_within = "AUC0-last within-arm CV (%)", Cmax_mean = "Cmax mean (mg/L)")[metric])]
  ULAB <- c(tdu12265_300 = "TDU12265 300 mg", tdu12265_600 = "TDU12265 600 mg", PKM12350 = "PKM12350", PKM14161 = "PKM14161", `R668-HV-1108` = "R668-HV-1108"); stopifnot(all(U$unit %in% names(ULAB)))
  U[, unit_lab := paste0(ULAB[unit], fifelse(dev_status == "external", " (out-of-sample)", ""))]   # 개발 자료 여부는 모델마다 달라 행마다 y 눈금을 따로 둔다(개정 1 표기, 계획서 16.5절)
  U[, unit_lab := factor(unit_lab, unique(unit_lab[order(match(unit, rev(names(ULAB))))]))]
  g1 <- ggplot(U, aes(y = unit_lab)) + geom_errorbarh(aes(xmin = pi_lo, xmax = pi_hi, colour = weight_lab), height = 0.25, position = position_dodge(width = 0.5)) +
    geom_point(aes(x = pred_median, colour = weight_lab), position = position_dodge(width = 0.5), size = 1.8) + geom_point(aes(x = observed), shape = 4, size = 2.6, stroke = 1, colour = VIZ$ink) +
    facet_grid(model_lab ~ metric_lab, scales = "free") + scale_colour_manual(values = COLS3[1:2], name = NULL) + labs(x = "Predicted median and 95% prediction interval; cross = observed", y = NULL,
    title = "Posterior predictive check (Li 2020 healthy-volunteer studies)",
    subtitle = "Five rule units; parallel-group studies pool both arms; rule R4 judges out-of-sample units against the 99% interval") + theme_dupi()
  save_fig(g1, "fig1_validation", 12, 6.5, U)
  # fig2: [3] OC 곡선(130명, M1, 75 kg 경계)
  OCP <- fread(file.path(SD, "oc_pass.csv"))[n_rand == NBASE & model == "M1" & boundary == "fixed75"]
  OCP[, `:=`(sched_lab = factor(SCHL[schedule], SCHL), panel = sprintf("%s, %d mg", c(S1 = "Global", S2 = "Korean")[population], dose_mg), model_lab = factor(MODS[pk_model], MODS))]
  g2 <- ggplot(OCP, aes(g, pass_pct / 100, colour = sched_lab, linetype = sched_lab)) + geom_ribbon(aes(ymin = pass_lo / 100, ymax = pass_hi / 100, fill = sched_lab, group = sched_lab), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) + scale_linetype_manual(values = c("solid", "solid", "22"), name = NULL) +   # 확장 일정(점선)이 확정 일정 위에 그려져 겹쳐도 둘 다 보인다
    geom_vline(xintercept = c(0.8, 1.25), linetype = 2, colour = VIZ$muted) + facet_grid(model_lab ~ panel) + scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25)) +
    scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = COLS3, name = NULL) + scale_fill_manual(values = COLS3, guide = "none") +
    labs(x = "True GMR (log scale)", y = "Probability of concluding similarity", title = "Operating characteristic curves", subtitle = "130 randomized per arm, stratified analysis, 75 kg boundary; band = Wilson 95% CI") + theme_dupi() +
    theme(panel.spacing.x = unit(1.4, "lines"))
  save_fig(g2, "fig2_oc_curves", 12, 6.5, OCP, lim200 = TRUE)
  # fig3: [4] 기전 OC(130명, M1, 확정 일정)
  MP <- fread(file.path(SD, "mech_pass.csv"))[n_rand == NBASE & model == "M1" & schedule == "B0"]
  idn <- MP[mechanism == "identity" & set == "mech"]; MPm <- rbind(MP[mechanism != "identity"], rbindlist(lapply(c("F", "ka", "ke", "Vmax"), function(m) copy(idn)[, mechanism := m])))
  MPm[, `:=`(diff_pct = 100 * (multiplier - 1), dose_lab = paste(dose_mg, "mg"), panel = sprintf("%s\n%s", c(k2020 = "Model 1", k2016 = "Kovalenko 2016")[pk_model], c(S1 = "Global", S2 = "Korean")[population]))]
  MPm[, panel := factor(panel, c("Model 1\nGlobal", "Model 1\nKorean", "Kovalenko 2016\nGlobal", "Kovalenko 2016\nKorean"))]; stopifnot(!anyNA(MPm$panel))
  g3 <- ggplot(MPm[order(multiplier)], aes(diff_pct, pass_pct / 100, colour = dose_lab)) + geom_line(aes(group = interaction(dose_lab, panel, mechanism)), linewidth = 0.6) + geom_point(aes(shape = set), size = 1.4) +
    geom_hline(yintercept = 0.5, linetype = 3, colour = VIZ$muted) + facet_grid(panel ~ mechanism) + scale_y_continuous(labels = scales::percent) +
    scale_colour_manual(values = COLS3[1:2], name = NULL) + scale_shape_manual(values = c(mech = 16, mechsupp = 1), labels = c(mech = "Main range", mechsupp = "Supplementary"), name = NULL) +
    labs(x = "Test-arm parameter difference (%)", y = "Probability of concluding similarity", title = "Mechanism-specific operating characteristics",
         subtitle = "Confirmed schedule, 130 randomized per arm, stratified analysis; dotted line = 50%\nModel 1 = Kovalenko 2020 Model 1 (primary model); Kovalenko 2016 = sensitivity model") + theme_dupi()
  save_fig(g3, "fig3_mechanism_oc", 12, 8, MPm, lim200 = TRUE)
  # fig4: [5] AUROC 대 차이 크기(양방향 합침, 주 통계량)
  AUR <- fread(file.path(SD, "auroc.csv"))[direction == "pooled" & role == "primary" & pk_model == PRIMARY & dose_mg == as.numeric(CFG$doses_mg$primary)]
  AUR[, `:=`(sched_lab = factor(SCHL[schedule], SCHL), ep_lab = c(AUCinf = "AUC0-inf", pAUC_0_14 = "pAUC 0 to 14 d", pAUC_14_28 = "pAUC 14 to 28 d", pAUC_28_56 = "pAUC 28 to 56 d")[endpoint],
             panel = paste(c(S1 = "Global", S2 = "Korean")[population], mechanism))]
  AUR[endpoint == "AUCinf_A", ep_lab := "AUC0-inf"]
  g4 <- ggplot(AUR, aes(size_pct, auroc, colour = sched_lab, group = sched_lab)) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.6, position = position_dodge(width = 1.6)) +
    geom_point(aes(shape = not_interpretable %in% TRUE), position = position_dodge(width = 1.6), size = 1.6) + geom_line(position = position_dodge(width = 1.6)) +
    facet_grid(ep_lab ~ panel) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "Interpretable", `TRUE` = "Zero share above 50%"), name = NULL) +
    scale_colour_manual(values = COLS3, name = NULL) + scale_x_continuous(breaks = c(5, 10, 15)) + coord_cartesian(ylim = c(0.4, 1)) +
    labs(x = "Injected difference (%)", y = "AUROC (95% cluster bootstrap CI)", title = "Secondary discrimination of endpoints without a margin",
         subtitle = sprintf("%s, %d mg, 130 randomized per arm; both directions pooled", MODL[[PRIMARY]], as.numeric(CFG$doses_mg$primary))) + theme_dupi()
  save_fig(g4, "fig4_auroc", 14, 9, AUR)
  # fig5: [6] 편향·정밀도 상자 그림(주 모델, 두 용량, 두 집단)
  IS <- readRDS(file.path(RD, "individual", "individual_subjects.rds"))[pk_model == PRIMARY]
  E <- rbind(IS[is.finite(AUClast) & AUClast > 0, .(pk_model, population, dose_mg, schedule, endpoint = "AUC0-last (dense truth)", err = 100 * (AUClast / AUClast_true - 1))],
             IS[lambda_ok & reliable & is.finite(AUCINF_obs), .(pk_model, population, dose_mg, schedule, endpoint = "AUC0-inf (reliable set)", err = 100 * (AUCINF_obs / AUCinf_true - 1))],
             IS[is.finite(Cmax), .(pk_model, population, dose_mg, schedule, endpoint = "Cmax", err = 100 * (Cmax / Cmax_true - 1))])
  E[, `:=`(sched_lab = factor(SCHL[schedule], SCHL), panel = sprintf("%s, %d mg", c(S1 = "Global", S2 = "Korean")[population], dose_mg))]
  Es <- E[, .(n = .N, q05 = quantile(err, 0.05), q25 = quantile(err, 0.25), q50 = median(err), q75 = quantile(err, 0.75), q95 = quantile(err, 0.95), rmse = sqrt(mean(err^2))), by = .(pk_model, population, dose_mg, schedule, endpoint)]
  Es[, `:=`(sched_lab = factor(SCHL[schedule], SCHL), panel = sprintf("%s, %d mg", c(S1 = "Global", S2 = "Korean")[population], dose_mg))]
  g5 <- ggplot(Es, aes(x = sched_lab, fill = sched_lab)) + geom_boxplot(aes(ymin = q05, lower = q25, middle = q50, upper = q75, ymax = q95), stat = "identity", width = 0.6, colour = VIZ$ink2) +
    geom_hline(yintercept = 0, colour = VIZ$muted) + facet_grid(endpoint ~ panel, scales = "free_y") + scale_fill_manual(values = COLS3, name = NULL) + scale_x_discrete(labels = function(x) SCHS[match(x, SCHL)]) +
    labs(x = NULL, y = "Percent error against the model truth", title = "Bias and precision of NCA parameters by sampling schedule",
         subtitle = { na <- range(Es[endpoint != "AUC0-inf (reliable set)", n]); nr <- range(Es[endpoint == "AUC0-inf (reliable set)", n])
           sprintf("%s, %s simulated subjects per cell (AUC0-last, Cmax); box = quartiles, whiskers = 5th and 95th percentiles\nAUC0-inf uses the reliable-terminal-phase subset (%s to %s subjects per box)",
                   MODL[[PRIMARY]], if (na[1] == na[2]) format(na[1], big.mark = ",") else paste(format(na, big.mark = ","), collapse = " to "), format(nr[1], big.mark = ","), format(nr[2], big.mark = ",")) }) + theme_dupi()
  save_fig(g5, "fig5_bias_precision", 13, 8, Es, lim200 = TRUE)
  # fig6: [7] 검정력 대 n
  PSIM <- fread(file.path(RD, "samplesize", "power_sim.csv"))[schedule == "B0" & boundary == "fixed75" & model == "M1"]
  PAN <- fread(file.path(RD, "samplesize", "power_analytic.csv"))[schedule == "B0" & boundary == "fixed75" & model == "M1"]
  PSIM[, panel := sprintf("%s, %d mg", c(S1 = "Global", S2 = "Korean")[population], dose_mg)]
  PAN[, panel := sprintf("%s, %d mg", c(S1 = "Global", S2 = "Korean")[population], dose_mg)]
  PL <- rbind(PSIM[, .(pk_model, panel, n_rand, true_gmr, what = "Simulated (Wilson 95% CI)", value = power_pct, lo, hi)],
              PAN[, .(pk_model, panel, n_rand, true_gmr, what = "Analytic, model CV", value = power_model_cv_pct, lo = NA_real_, hi = NA_real_)],
              PAN[, .(pk_model, panel, n_rand, true_gmr, what = "Assurance, Li 2020 300 mg CV (registered)", value = assurance_cv_pct, lo = NA_real_, hi = NA_real_)],
              PAN[, .(pk_model, panel, n_rand, true_gmr, what = "Assurance, CV shifted to the cell (post hoc)", value = assurance_cv_cellshift_pct, lo = NA_real_, hi = NA_real_)])
  WH <- c("Analytic, model CV", "Assurance, Li 2020 300 mg CV (registered)", "Assurance, CV shifted to the cell (post hoc)", "Simulated (Wilson 95% CI)")
  stopifnot(all(PL$what %in% WH)); PL[, what := factor(what, WH)]
  RL <- CJ(pk_model = c("k2020", "k2016"), true_gmr = sort(unique(PL$true_gmr)), sorted = FALSE)[, gl := sprintf("%s\nTrue GMR %.2f", c(k2020 = "Model 1", k2016 = "Kovalenko 2016")[pk_model], true_gmr)]
  PL <- merge(PL, RL, by = c("pk_model", "true_gmr")); PL[, gl := factor(gl, RL$gl)]
  g6 <- ggplot(PL, aes(n_rand, value / 100, colour = what, linetype = what)) + geom_errorbar(aes(ymin = lo / 100, ymax = hi / 100), width = 2, na.rm = TRUE) + geom_line() + geom_point(size = 1.2) +
    facet_grid(gl ~ panel) + scale_y_continuous(labels = scales::percent) + scale_colour_manual(values = c(VIZ$s1, VIZ$s2, VIZ$s2, VIZ$s3), name = NULL) +
    scale_linetype_manual(values = c("solid", "solid", "22", "solid"), name = NULL) + guides(colour = guide_legend(nrow = 2)) + geom_vline(xintercept = NBASE, linetype = 3, colour = VIZ$muted) +
    labs(x = "Randomized subjects per arm (10% withdrawal)", y = "Probability of concluding similarity", title = "Power and assurance against sample size",
         subtitle = "Confirmed schedule, 75 kg boundary, stratified analysis; dotted line = 130 per arm\nThe registered assurance draws its CV from the Li 2020 300 mg arms at both doses; the post hoc assurance shifts it to the cell by the model (plan section 13.4)\nModel 1 = Kovalenko 2020 Model 1 (primary model); Kovalenko 2016 = sensitivity model") + theme_dupi()
  save_fig(g6, "fig6_power_n", 13, 11, PL, lim200 = TRUE)
  # fig7: [8] forest plot(300 mg, 두 모델, S1·S2)
  WR <- fread(file.path(RD, "weight", "weight_ratio.csv"))[population %in% c("S1", "S2") & dose_mg == as.numeric(CFG$doses_mg$primary)]
  WR[, `:=`(lab = sprintf("%s P%d (%.1f kg)", c(S1 = "Global", S2 = "Korean")[population], round(100 * prob), WT), met_lab = c(AUCinf = "AUC0-inf", Cmax = "Cmax", AUClast_dense_B0 = "AUC0-last (confirmed schedule)")[metric],
            model_lab = factor(MODL[pk_model], MODL))]
  WR[, lab := factor(lab, rev(unique(lab[order(population, prob)])))]
  g7 <- ggplot(WR, aes(y = lab)) + geom_vline(xintercept = 1, colour = VIZ$muted) + geom_errorbarh(aes(xmin = pi_q05, xmax = pi_q95), height = 0, colour = VIZ$grid, linewidth = 2.2) +
    geom_errorbarh(aes(xmin = paired_q05, xmax = paired_q95), height = 0.3, colour = VIZ$s1) + geom_point(aes(x = typical_ratio), colour = VIZ$s1, size = 1.8) +
    facet_grid(model_lab ~ met_lab) + scale_x_log10() +
    labs(x = "Ratio to 75 kg (log scale): point = typical subject; bar = 90% interval of the paired ratio\ngrey band = 90% prediction interval with between-subject variability", y = NULL,
         title = "Body-weight effect on exposure at enrolled-weight percentiles", subtitle = sprintf("%d mg, simulated through the weight covariate only", as.numeric(CFG$doses_mg$primary))) + theme_dupi()
  save_fig(g7, "fig7_weight_forest", 13, 7, WR)
  append_run_log(lf, "done")
}

# ================================================================ 요약 표·결론·가정 ================================================================
if (mode == "summary") {
  lf <- start_run_log("m15_summary", master_seed = MS, run_mode = "final")
  RULE <- fread(file.path(RD, "validation", "validation_rule.csv"))
  # 개정 1(계획서 16.5절): '실패 규칙' -> 'misfit rule' 표기. 등록 산출물(validation/)의 결과 범주 문자열은 그대로 두고 여기서 바꿔 쓴다
  oc_word <- function(x) sub("^passes under both weight assumptions$", "does not meet the misfit rule under either weight assumption",
                             sub("^fails under ", "meets the misfit rule under ", x))
  OCP <- fread(file.path(SD, "oc_pass.csv")); CR <- fread(file.path(SD, "oc_consumer_risk.csv")); PD <- fread(file.path(SD, "oc_paired_diff.csv"), na.strings = c("", "NA"))
  CRT <- fread(file.path(SD, "mech_crossing.csv")); AUR <- fread(file.path(SD, "auroc.csv"))
  SY <- fread(file.path(RD, "samplesize", "syneos_comparison.csv")); PAN <- fread(file.path(RD, "samplesize", "power_analytic.csv"))
  WR <- fread(file.path(RD, "weight", "weight_ratio.csv"))
  base <- list(population = "S1", dose_mg = as.numeric(CFG$doses_mg$primary), schedule = "B0", boundary = "fixed75", n_rand = NBASE, model = "M1")
  pick <- function(d, pm, extra = list()) { s <- c(base, extra); z <- d[d$pk_model == pm]; for (nm in intersect(names(s), names(z))) z <- z[z[[nm]] == s[[nm]]]; z }
  gsel <- function(z, gg) z[abs(z$g - gg) < 1e-9]
  rows <- list(); addr <- function(item, k2020, k2016, src) rows[[length(rows) + 1]] <<- data.table(item = item, model_1 = k2020, kovalenko_2016 = k2016, source = src)
  vo <- function(pm) { r <- RULE[pk_model == pm & wmode == "primary" & rule == "any"]; sprintf("%s (fixed weight); %s", if (r$failure) "misfit rule met" else "misfit rule not met", oc_word(RULE[pk_model == pm, outcome][1])) }
  addr("Posterior predictive check outcome (Li 2020 healthy-volunteer studies; section 6.2)", vo("k2020"), vo("k2016"), "results/m15/validation/validation_rule.csv")
  cell <- function(pm, gg, col = "pass", d = OCP) { z <- gsel(pick(d, pm), gg); if (!nrow(z)) return("-"); sprintf("%s%% (%s to %s)", fmt(z[[paste0(col, "_pct")]]), fmt(z[[paste0(col, "_lo")]]), fmt(z[[paste0(col, "_hi")]])) }
  for (gg in c(0.95, 1.00)) addr(sprintf("Power at true GMR %.2f, primary scenario", gg), cell("k2020", gg), cell("k2016", gg), "results/m15/summary/oc_pass.csv")
  NCR <- unique(CR[pk_model == PRIMARY, n_trials]); stopifnot(length(NCR) == 1, NCR == NT_EXT)       # 소비자 위험 연구 수는 파일에서
  for (gg in c(0.80, 1.25)) for (col in c("pass", "auc_only", "cmax_only")) {
    lab <- c(pass = "joint, both endpoints at the margin", auc_only = "AUC0-last-only pass", cmax_only = "Cmax-only pass")[[col]]
    addr(sprintf("Consumer risk at true GMR %.2f (%s)", gg, lab), cell("k2020", gg, col, if (nrow(CR)) CR else OCP), cell("k2016", gg, col),
         sprintf("results/m15/summary/oc_consumer_risk.csv (Model 1, %s studies); oc_pass.csv (Kovalenko 2016, %s)", fmt(NCR, 0), fmt(NT_MAIN, 0)))
  }
  pdc <- function(pm, ct_, a_, gg = 0.95) { z <- PD[PD$pk_model == pm & PD$contrast == ct_ & PD$design_a == a_ & abs(PD$g - gg) < 1e-9]
    for (nm in c("population", "dose_mg", "schedule", "boundary", "model")) if (nm %in% names(z) && !is.na(match(nm, names(base)))) { v <- base[[nm]]; if (!all(is.na(z[[nm]]))) z <- z[is.na(z[[nm]]) | z[[nm]] == v] }
    if (nrow(z) != 1) { if (pm == PRIMARY) stop(sprintf("요약 표 %s/%s: 기본 시나리오 칸 %d개", ct_, a_, nrow(z))); return("-") }
    sprintf("%s pp (%s to %s)", fmt(z$diff_pp), fmt(z$lo_pp), fmt(z$hi_pp)) }
  for (ct in list(c("stratum term", "M1", "Effect of the stratum term (stratified minus unstratified)"), c("boundary", "median", "Effect of the median boundary (median minus 75 kg)"),
                  c("population", "S2", "Effect of the population (Korean minus global)"), c("dose", "200", "Effect of the dose (200 minus 300 mg; 200 mg results limited by the known AUC0-last under-prediction for the 200 mg presentation)"),
                  c("schedule", "RM36_50", "Effect of the reduced schedule (minus confirmed)"), c("schedule", "ADD64_71", "Effect of the extended schedule (minus confirmed)")))
    addr(sprintf("%s on power at true GMR 0.95", ct[3]), pdc("k2020", ct[1], ct[2]), pdc("k2016", ct[1], ct[2]), "results/m15/summary/oc_paired_diff.csv")
  for (mech in c("F", "ka", "ke", "Vmax")) {
    cc <- function(pm) { z <- CRT[pk_model == pm & population == "S1" & dose_mg == base$dose_mg & schedule == "B0" & mechanism == mech]
      paste(vapply(c("down", "up"), function(sd_) { r <- z[side == sd_]; if (!nrow(r)) return("-"); if (!r$reached) return(sprintf("not reached (%s%% at %s%%)", fmt(100 * r$p_at_grid_end), fmt(100 * (r$grid_end - 1), 0)))
        bl <- function(x) if (is.na(x)) "beyond the grid end" else fmt(x)
        sprintf("%s%% (%s to %s; AUC0-last GMR %s, Cmax GMR %s)", fmt(r$crossing_pct), bl(r$lo_pct), bl(r$hi_pct), fmt(r$induced_gmr_AUClast, 3), fmt(r$induced_gmr_Cmax, 3)) }, ""), collapse = "; ") }
    addr(sprintf("50%% crossing, %s (lower; upper)", mech), cc("k2020"), cc("k2016"), "results/m15/summary/mech_crossing.csv")
  }
  for (ep in c("AUCinf_A", "pAUC_0_14_meanratio", "pAUC_14_28_meanratio", "pAUC_28_56_meanratio")) {
    ac <- function(pm) { z <- AUR[pk_model == pm & population == "S1" & dose_mg == base$dose_mg & schedule == "B0" & statistic == ep & size_pct == 10 & direction == "pooled"]
      if (!nrow(z)) return("-"); paste(sprintf("%s %s (%s to %s)%s", z$mechanism, fmt(z$auroc, 3), fmt(z$lo, 3), fmt(z$hi, 3), fifelse(z$not_interpretable %in% TRUE, " [zero share above 50%]", "")), collapse = "; ") }
    addr(sprintf("AUROC at a 10%% difference (both directions pooled), %s", c(AUCinf_A = "AUC0-inf (absolute log GMR, stratified analysis)", pAUC_0_14_meanratio = "partial AUC 0 to 14 days (absolute log ratio of arithmetic means)",
                                                    pAUC_14_28_meanratio = "partial AUC 14 to 28 days (absolute log ratio of arithmetic means)",
                                                    pAUC_28_56_meanratio = "partial AUC 28 to 56 days (absolute log ratio of arithmetic means)")[[ep]]), ac("k2020"), ac("k2016"), "results/m15/summary/auroc.csv")
  }
  ss <- function(pm, col) { z <- SY[pk_model == pm & cv_source == "model, enrolled global population" & true_gmr == 0.95 & model == "M1"]; if (!nrow(z)) return("-")
    zr <- z[cmax_rule == "ratio"][[col]]; zf <- z[cmax_rule == "fixed"][[col]]; stopifnot(length(zr) == 1, length(zf) == 1)
    if (zr == zf) as.character(zr) else sprintf("%s (fixed Cmax rule %s)", zr, zf) }
  addr("Randomized per arm for 85% power at true GMR 0.95 (model CV; Cmax SD by the ratio rule)", ss("k2020", "n_rand_85"), ss("k2016", "n_rand_85"), "results/m15/samplesize/syneos_comparison.csv")
  addr("Randomized per arm for 90% power at true GMR 0.95 (model CV; Cmax SD by the ratio rule)", ss("k2020", "n_rand_90"), ss("k2016", "n_rand_90"), "results/m15/samplesize/syneos_comparison.csv")
  asr <- function(pm) { z <- pick(PAN, pm, list())[true_gmr == 0.95]; if (!nrow(z)) "-" else sprintf("%s%%", fmt(z$assurance_cv_pct)) }
  addr("Assurance over CV uncertainty at 130 per arm, true GMR 0.95", asr("k2020"), asr("k2016"), "results/m15/samplesize/power_analytic.csv")
  syz <- function(pm) { z <- SY[pk_model == pm & true_gmr == 0.95 & model == "M1"]; zr <- z[cmax_rule == "ratio"]; zf <- z[cmax_rule == "fixed"][match(zr$cv_source, cv_source)]
    stopifnot(nrow(zr) == nrow(zf), !anyNA(zf$cv_source)); list(r = zr, f = zf) }
  sy <- function(pm) { z <- syz(pm); paste(sprintf("%s: CV %s%%, power %s%% (ratio rule) or %s%% (fixed rule)", z$r$cv_source, fmt(z$r$cv_pct), fmt(z$r$power_n117_eval_pct), fmt(z$f$power_n117_eval_pct)), collapse = "; ") }
  sn_ <- function(pm) { z <- syz(pm); paste(sprintf("%s: %s and %s (ratio rule); %s and %s (fixed rule)", z$r$cv_source, z$r$n_rand_85, z$r$n_rand_90, z$f$n_rand_85, z$f$n_rand_90), collapse = "; ") }
  addr("Protocol comparison at 117 evaluable per arm, true GMR 0.95, by CV source (Cmax SD by the ratio rule or the fixed rule of the stored protocol figures)", sy("k2020"), sy("k2016"), "results/m15/samplesize/syneos_comparison.csv")
  addr("Randomized per arm for 85% and 90% power at true GMR 0.95, by CV source", sn_("k2020"), sn_("k2016"), "results/m15/samplesize/syneos_comparison.csv")
  wr <- function(pm, pr) { z <- WR[pk_model == pm & population == "S1" & dose_mg == base$dose_mg & abs(prob - pr) < 1e-9]
    ml <- c(AUCinf = "AUC0-inf", Cmax = "Cmax", AUClast_dense_B0 = "AUC0-last (dense truth, confirmed schedule)"); stopifnot(all(z$metric %in% names(ml)))
    paste(sprintf("%s %s (paired 90%%: %s to %s)", ml[z$metric], fmt(z$typical_ratio, 3), fmt(z$paired_q05, 3), fmt(z$paired_q95, 3)), collapse = "; ") }
  for (pr in c(0.05, 0.95)) addr(sprintf("Weight effect at the %dth percentile (global population), ratio to 75 kg", round(100 * pr)), wr("k2020", pr), wr("k2016", pr), "results/m15/weight/weight_ratio.csv")
  # 적용 가능성(계획서 6.3절): 개발 자료 전체의 체중 범위는 저장소에 없다(논문 표 접근 불가). 알려진 내부 연구의 체중 범위(S1 범위 = R668-AS-0907, 두 모델 내부;
  #   S2 범위 = TDU12265, Model 1 내부; 계획서 6.2절 내부·외부 구분) 밖의 등록 체중 비율을 외삽으로 보고한다(등록 집단 200,000명, 74와 같은 시드)
  PP <- CFG$populations; KNOWN <- list(k2020 = list(`R668-AS-0907` = PP$S1$weight, TDU12265 = PP$S2$weight), k2016 = list(`R668-AS-0907` = PP$S1$weight))
  APL <- rbindlist(lapply(c("S1", "S2"), function(pv) {
    wt <- m15_subjects(200000L, PP[[pv]], CFG$inclusion, as.numeric(PP$weight_height_rho), as.numeric(PP$sex_ratio_male), derive_seed(MS, "popq"))$WT
    rbindlist(lapply(names(KNOWN), function(pm) { kr <- KNOWN[[pm]]; lo <- min(vapply(kr, function(z) as.numeric(z$lo), 1)); hi <- max(vapply(kr, function(z) as.numeric(z$hi), 1))
      stopifnot(length(kr) == 1 || (as.numeric(kr[[2]]$lo) <= as.numeric(kr[[1]]$hi)))                  # 두 범위가 겹쳐 하나의 구간
      data.table(pk_model = pm, population = pv, enrolled_min_kg = min(wt), enrolled_max_kg = max(wt), internal_studies_with_known_range = paste(names(kr), collapse = " and "),
                 known_range_lo_kg = lo, known_range_hi_kg = hi, outside_known_pct = 100 * mean(wt < lo | wt > hi),
                 development_data_full_range = "not recorded in the repository (publication tables not accessible)") })) }))
  m15_fwrite_cond(APL, file.path(SD, "m15_applicability.csv"), COND)
  ap <- function(pm) { z <- APL[pk_model == pm]; paste(sprintf("%s population %s to %s kg: %s%% outside %s (%s to %s kg)", c(S1 = "global", S2 = "Korean")[z$population], fmt(z$enrolled_min_kg), fmt(z$enrolled_max_kg),
                                                       fmt(z$outside_known_pct), z$internal_studies_with_known_range, fmt(z$known_range_lo_kg), fmt(z$known_range_hi_kg)), collapse = "; ") }
  addr("Applicability: enrolled weights outside the documented weight ranges of development studies (reported as outside the documented range; the full development-data range is not recorded; amendment section 16.9)", ap("k2020"), ap("k2016"),
       "results/m15/summary/m15_applicability.csv")
  if (RULE[pk_model == PRIMARY & wmode == "primary", failure][1]) addr("Conditionality", "All results conditional on a primary model that met the pre-specified misfit rule of the posterior predictive check (plan section 6.2; cause under diagnosis, amendment section 16.7)", "",
                                                                         "results/m15/validation/validation_rule.csv")
  TAB <- rbindlist(rows); fwrite(TAB, file.path(SD, "m15_summary_table.csv"))
  # 약어: 문서에 나오는 것만, 처음 나오는 순서로(본문 앞에 한 줄)
  ABBR <- c(`AUC0-inf` = "area under the concentration-time curve from time zero to infinity", `AUC0-last` = "area under the concentration-time curve from time zero to the last quantifiable concentration",
            AUC = "area under the concentration-time curve", AUROC = "area under the receiver operating characteristic curve", CI = "confidence interval", Cmax = "maximum concentration",
            CV = "coefficient of variation", F = "bioavailability", GMR = "geometric mean ratio", ka = "absorption rate constant", ke = "linear elimination rate constant",
            Vmax = "maximum Michaelis-Menten elimination rate", pp = "percentage points", RSE = "relative standard error", SD = "standard deviation", PK = "pharmacokinetic", NCA = "non-compartmental analysis",
            SE = "standard error", MC = "Monte Carlo", OC = "operating characteristic", UTC = "Coordinated Universal Time", lambda_z = "terminal elimination rate constant", pAUC = "partial AUC",
            QC = "quality control", CSV = "comma-separated values", ASCII = "American Standard Code for Information Interchange", AI = "artificial intelligence", ADA = "anti-drug antibodies", BLA = "Biologics License Application")
  abbr_line <- function(txt) { txt <- paste(txt, collapse = "\n")
    pos <- vapply(names(ABBR), function(a) as.integer(regexpr(paste0("(?<![A-Za-z0-9_])", a, "(?![A-Za-z0-9_])"), txt, perl = TRUE)), 1L); pos <- sort(pos[pos > 0])
    if (!length(pos)) return(character(0)); c(paste0("Abbreviations: ", paste(sprintf("%s, %s", names(pos), ABBR[names(pos)]), collapse = "; "), "."), "") }
  body_ <- c("Primary scenario: global healthy-volunteer population, 300 mg, confirmed schedule, 75 kg boundary, 130 randomized per arm, stratified analysis. Every value is read from the file in the last column.",
             "Intervals: Wilson 95% confidence interval (CI) for probabilities; paired Tango 95% CI for differences in power; study-level cluster bootstrap 95% CI for 50% crossings and AUROC; 5th to 95th percentile of the paired individual ratio for the weight effect.", "",
             "| Item | Kovalenko 2020 Model 1 (primary) | Kovalenko 2016 (sensitivity) | Source |", "|---|---|---|---|",
             sprintf("| %s | %s | %s | %s |", TAB$item, TAB$model_1, TAB$kovalenko_2016, TAB$source))
  md <- c("# Results summary table (model analysis plan 2026-10-20)", "", abbr_line(body_), body_)
  if (RULE[pk_model == PRIMARY & wmode == "primary", failure][1]) md <- c(md, "", "All results are conditional on a primary model that met the pre-specified misfit rule of the posterior predictive check (plan section 6.2; cause under diagnosis, amendment section 16.7).")
  md <- c(md, "", "Parameter uncertainty is not propagated in this table; every value is conditional on the published parameter point estimates. Parameter uncertainty is evaluated separately (plan section 16.8, results/m15/amend1/param_assurance).")
  # 이탈·사후 분석 표(계획서 13.1절): 계획서 13.3·13.4절의 번호 항목 첫 문장을 자동으로 옮긴다
  sec_items <- function(h0, h1, type) { i0 <- grep(h0, pl_); i1 <- grep(h1, pl_); stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0)
    x <- pl_[(i0 + 1):(i1 - 1)]; it <- grep("^[0-9]+\\. ", x, value = TRUE)
    first2 <- vapply(strsplit(sub("^[0-9]+\\. ", "", it), "(?<=\\.)\\s+(?=[A-Z])", perl = TRUE), function(z) paste(head(z, 2), collapse = " "), "")   # 첫 두 문장
    data.table(number = as.integer(sub("^([0-9]+)\\..*", "\\1", it)), type = type, summary = first2, dlog = vapply(regmatches(it, regexec("\\((?:amendment|addendum) 1, (D-0[0-9]{2})\\)", it)), function(m) if (length(m) == 2) m[2] else "D-068", "")) }
  pl_ <- readLines(proj_path("plans", "M15_analysis_plan_20261020.md"), warn = FALSE)
  DEV <- rbind(sec_items("^### 13\\.3 ", "^### 13\\.4 ", "deviation after registration"), sec_items("^### 13\\.4 ", "^## 14\\. ", "post hoc analysis"))
  DEV[, `:=`(plan_section = fifelse(type == "post hoc analysis", "13.4", "13.3"), decision_log = dlog)][, dlog := NULL]   # 개정 1 항목은 본문의 결정 기록 번호(D-069, 정정은 D-070), 추록 1은 D-071
  m15_fwrite_cond(DEV, file.path(SD, "m15_deviations.csv"), COND)
  md <- c(md, "", "## Deviations and post hoc analyses (plan sections 13.3 and 13.4; full text in the plan)", "", "| Section | Number | Type | Summary |", "|---|---|---|---|",
          sprintf("| %s | %d | %s | %s |", DEV$plan_section, DEV$number, DEV$type, DEV$summary))
  ia <- grep("^Abbreviations: ", md); if (length(ia) == 1) md[ia] <- abbr_line(md[-ia])[1]          # 약어 줄을 표 전체 기준으로 다시
  writeLines(md, file.path(SD, "m15_summary_table.md"))
  # 핵심 결론(최대 5줄, 2.1절 해석 기준, 전제 점검)
  pz <- function(gg) gsel(pick(OCP, PRIMARY), gg); crz <- function(gg) gsel(pick(if (nrow(CR)) CR else OCP, PRIMARY), gg)
  p95 <- pz(0.95); c80 <- crz(0.80); c125 <- crz(1.25); stopifnot(nrow(p95) == 1, nrow(c80) == 1, nrow(c125) == 1)
  adequate <- p95$pass_pct >= 90 & max(c80$pass_hi, c125$pass_hi) <= 5
  rv <- RULE[pk_model == PRIMARY & wmode == "primary" & rule == "any"]
  # 결론 1(개정 1, 계획서 16.5절): 2.1절 3번(검출 확률에 따른 'reproduced'/'not contradicted')은 16.7절 완료까지 중지. 중간 문구는 결과 파일에서 생성
  MLAB <- c(k2020 = "Model 1", k2016 = "Kovalenko 2016"); RMET <- c(AUClast_mean = "AUC0-last mean", AUClast_cv_within = "within-arm CV of AUC0-last", Cmax_mean = "Cmax mean")
  ULB <- c(tdu12265_300 = "TDU12265 300 mg", tdu12265_600 = "TDU12265 600 mg")
  met_txt <- function(pm, md) { r <- RULE[pk_model == pm & wmode == md & rule != "any" & met == TRUE]
    if (!nrow(r)) return("")
    paste(sprintf("rule %s: the %s %s %s outside %s %d%% prediction interval", r$rule, vapply(strsplit(r$units_outside, ";"), function(u) paste(fifelse(u %in% names(ULB), ULB[u], u), collapse = " and "), ""),
                  RMET[r$metric], fifelse(r$n_outside > 1, "are", "is"), fifelse(r$n_outside > 1, "their", "its"), round(100 * r$pi_level)), collapse = "; ") }
  cat_txt <- function(pm, ko = FALSE) { f1 <- RULE[pk_model == pm & wmode == "primary", failure][1]; f2 <- RULE[pk_model == pm & wmode == "wunc", failure][1]
    if (ko) return(if (f1 && f2) "두 체중 가정 모두에서 충족" else if (f1) "고정 체중 가정에서만 충족" else if (f2) "평균 불확실 체중 가정에서만 충족" else "두 가정 모두에서 미충족")
    if (f1 && f2) sprintf("met the misfit rule under both weight assumptions (%s)", met_txt(pm, "primary")) else
      if (f1) sprintf("met the misfit rule under the fixed weight assumption only (%s) and did not meet it under the uncertain-mean weight assumption", met_txt(pm, "primary")) else
      if (f2) sprintf("met the misfit rule under the uncertain-mean weight assumption only (%s)", met_txt(pm, "wunc")) else "did not meet it under either weight assumption" }
  OTHER <- setdiff(unique(RULE$pk_model), PRIMARY); stopifnot(length(OTHER) == 1)
  # 결론 3: 설계 요인별로 기준을 충족한 칸 수(주 모델, 참 GMR 0.95, 130명)와 기본 시나리오 칸의 충족 여부
  pd95 <- PD[pk_model == PRIMARY & abs(g - 0.95) < 1e-9]
  is_base <- function(z) (is.na(z$population) | z$population == base$population) & (is.na(z$dose_mg) | z$dose_mg == base$dose_mg) & (is.na(z$schedule) | z$schedule == base$schedule) &
    (is.na(z$boundary) | z$boundary == base$boundary) & (is.na(z$model) | z$model == base$model)
  pd95[, base_cell := is_base(pd95)]
  stopifnot(!anyNA(pd95$changes_oc))
  sg <- pd95[, .(k = sum(changes_oc), n = .N, in_base = any(changes_oc & base_cell), d_lo = if (any(changes_oc)) min(diff_pp[changes_oc]) else NA_real_,
                 d_hi = if (any(changes_oc)) max(diff_pp[changes_oc]) else NA_real_), by = .(contrast, design_a, design_b)][k > 0]
  DLAB <- c(RM36_50 = "reduced schedule", ADD64_71 = "extended schedule", B0 = "confirmed schedule", M1 = "stratified analysis", M0 = "unstratified analysis",
            median = "median boundary", fixed75 = "75 kg boundary", S2 = "Korean population", S1 = "global population", `200` = "200 mg", `300` = "300 mg")
  DLAB_KO <- c(RM36_50 = "축소 일정", ADD64_71 = "확장 일정", B0 = "확정 일정", M1 = "층 포함 분석", M0 = "층 없는 분석", median = "중앙값 경계", fixed75 = "75 kg 경계",
               S2 = "한국인 집단", S1 = "전체 집단", `200` = "200 mg", `300` = "300 mg")
  stopifnot(all(c(sg$design_a, sg$design_b) %in% names(DLAB)))
  sg[, rng := fifelse(d_lo == d_hi, sprintf("%s pp", fmt(d_lo)), sprintf("%s to %s pp", fmt(d_lo), fmt(d_hi)))]
  sg[, rng_ko := fifelse(d_lo == d_hi, sprintf("%s%%p", fmt(d_lo)), sprintf("%s~%s%%p", fmt(d_lo), fmt(d_hi)))]
  CLAB_KO <- c(`stratum term` = "층 항", boundary = "경계", schedule = "일정", population = "집단", dose = "용량")
  sig <- sg
  sig_txt <- if (nrow(sg)) paste(sprintf("%s (%s minus %s): %d of %d cells, %s%s%s", sg$contrast, DLAB[sg$design_a], DLAB[sg$design_b], sg$k, sg$n, sg$rng, fifelse(sg$in_base, ", including the primary scenario", ""), fifelse(sg$contrast == "dose", " (200 mg results limited by the known AUC0-last under-prediction for the 200 mg presentation)", "")), collapse = "; ") else "none"
  sig_txt_ko <- if (nrow(sg)) paste(sprintf("%s(%s 빼기 %s): %d/%d칸, %s%s%s", CLAB_KO[sg$contrast], DLAB_KO[sg$design_a], DLAB_KO[sg$design_b], sg$k, sg$n, sg$rng_ko, fifelse(sg$in_base, ", 기본 시나리오 포함", ""), fifelse(sg$contrast == "dose", "(200 mg 결과는 200 mg 제형의 알려진 AUC0-last 과소 예측으로 제한)", "")), collapse = "; ") else "없음"
  side_txt <- function(z, ko = FALSE) paste(sprintf("%s %s %s", z$mechanism, if (ko) fifelse(z$side == "down", "아래쪽", "위쪽") else fifelse(z$side == "down", "lower", "upper"),
    fifelse(z$reached, paste0(fmt(z$crossing_pct), "%"), sprintf(if (ko) "격자 안에서 도달 안 함(%2$s%% 차이에서 %1$s%%)" else "not reached within the grid (%s%% at %s%%)", fmt(100 * z$p_at_grid_end), fmt(100 * (z$grid_end - 1), 0)))),
    collapse = ", ")
  bc <- CRT[pk_model == PRIMARY & population == "S1" & dose_mg == base$dose_mg & schedule == "B0"]
  # 결론 5: 축소 일정의 검정력 차이와 측정 품질(주 모델, 기본 시나리오의 집단·용량)
  BP <- fread(file.path(RD, "individual", "bias_precision.csv")); SSI5 <- fread(file.path(RD, "individual", "ss_inputs_m15.csv"))
  bsel <- function(d, sc) { z <- d[d$pk_model == PRIMARY & d$population == base$population & d$dose_mg == base$dose_mg & d$schedule == sc]; z }
  q5 <- lapply(c(B0 = "B0", RM36_50 = "RM36_50"), function(sc) {
    b <- bsel(BP, sc); a <- bsel(AUR, sc)[statistic == "AUCinf_A" & mechanism == "F" & size_pct == 10 & direction == "pooled"]; v <- bsel(SSI5, sc)
    stopifnot(nrow(b[endpoint == "AUClast"]) == 1, nrow(b[endpoint == "AUCinf_ruleA"]) == 1, nrow(a) == 1, nrow(v) == 1)
    list(bias = b[endpoint == "AUClast", bias_median_pct], rel = 100 * b[endpoint == "AUCinf_ruleA", n] / b[endpoint == "AUClast", n], auroc = a$auroc, cv = v$cv_auc_pct) })
  rs <- pd95[contrast == "schedule" & design_a == "RM36_50" & base_cell == TRUE]; stopifnot(nrow(rs) == 1)
  shorter <- isTRUE(rs$changes_oc) && rs$diff_pp > 0 && q5$RM36_50$bias < q5$B0$bias && q5$RM36_50$auroc < q5$B0$auroc && q5$RM36_50$cv < q5$B0$cv
  c5 <- sprintf("5. The reduced schedule changes power in the primary scenario by %s pp (paired 95%% CI %s to %s) against the confirmed schedule, with model AUC0-last CV %s%% vs %s%%, median AUC0-last error against the dense truth %s%% vs %s%%, a reliable AUC0-inf in %s%% vs %s%% of subjects, and an AUC0-inf AUROC for a 10%% bioavailability difference (both directions pooled) of %s vs %s%s.",
                fmt(rs$diff_pp), fmt(rs$lo_pp), fmt(rs$hi_pp), fmt(q5$RM36_50$cv), fmt(q5$B0$cv), fmt(q5$RM36_50$bias), fmt(q5$B0$bias), fmt(q5$RM36_50$rel, 0), fmt(q5$B0$rel, 0),
                fmt(q5$RM36_50$auroc, 3), fmt(q5$B0$auroc, 3), if (shorter) "; the power gain reflects a shorter, less variable AUC0-last, not better measurement" else "")
  c5_ko <- sprintf("5. 축소 일정은 기본 시나리오 검정력을 확정 일정 대비 %s%%p(짝 95%% 구간 %s~%s) 바꾼다. 모델 AUC0-last CV %s%% 대 %s%%, 조밀 참값 대비 AUC0-last 중앙 오차 %s%% 대 %s%%, 신뢰할 AUC0-inf 비율 %s%% 대 %s%%, 생체이용률 10%% 차이(양방향 합침)의 AUC0-inf AUROC %s 대 %s%s.",
                   fmt(rs$diff_pp), fmt(rs$lo_pp), fmt(rs$hi_pp), fmt(q5$RM36_50$cv), fmt(q5$B0$cv), fmt(q5$RM36_50$bias), fmt(q5$B0$bias), fmt(q5$RM36_50$rel, 0), fmt(q5$B0$rel, 0),
                   fmt(q5$RM36_50$auroc, 3), fmt(q5$B0$auroc, 3), if (shorter) ". 검정력 증가는 더 짧고 덜 변하는 AUC0-last 때문이지 더 나은 측정 때문이 아니다" else "")
  HDR <- if (rv$failure) "All results are conditional on a primary model that met the pre-specified misfit rule of the posterior predictive check (plan section 6.2; cause under diagnosis, amendment section 16.7) and on the published parameter point estimates (parameter uncertainty is evaluated separately in amendment section 16.8)." else
    "All results are conditional on the published parameter point estimates (parameter uncertainty is evaluated separately in amendment section 16.8)."
  HDR_KO <- if (rv$failure) "주의: 모든 결과는 사후 예측 점검(Li 2020 건강인 시험)의 사전 등록 misfit 규칙을 충족한 주 모델(Kovalenko 2020 Model 1, 고정 체중 가정; 계획서 6.2절, 원인 진단 중 16.7절)과 공개 파라미터 점추정값에 조건부다(파라미터 불확실성은 개정 16.8절에서 따로 평가)." else
    "주의: 모든 결과는 공개 파라미터 점추정값에 조건부다(파라미터 불확실성은 개정 16.8절에서 따로 평가)."
  C <- c(sprintf("1. Posterior predictive check (Li 2020 healthy-volunteer studies): %s %s; %s %s. The cause is under diagnosis (amendment section 16.7).",
                 MLAB[[PRIMARY]], cat_txt(PRIMARY), MLAB[[OTHER]], cat_txt(OTHER)),
         sprintf("2. At 130 randomized per arm (primary scenario) the simulated power at a true GMR of 0.95 is %s%% and the joint consumer risk is %s%% at 0.80 and %s%% at 1.25 (single-endpoint maxima %s%% and %s%%); the design %s the pre-specified adequacy reading.",
                 fmt(p95$pass_pct), fmt(c80$pass_pct), fmt(c125$pass_pct), fmt(max(c80$auc_only_pct, c80$cmax_only_pct)), fmt(max(c125$auc_only_pct, c125$cmax_only_pct)), if (adequate) "meets" else "does not meet"),
         sprintf("3. Design factors that change power at a true GMR of 0.95 by at least 2 percentage points (pp) with a paired 95%% CI excluding 0 (primary model, 130 per arm; cells counted over populations, doses, schedules, boundaries and analysis models): %s.", sig_txt),
         sprintf("4. The conclusion probability falls to 50%% at test-arm differences of %s (primary scenario).", side_txt(bc)),
         c5)
  stopifnot(length(C) <= 5)
  writeLines(c("# Key conclusions (generated from the result files)", "", HDR, "", abbr_line(c(HDR, C)), C),
             file.path(SD, "m15_key_conclusions_en.md"))
  # 한국어 요약(내부용): 같은 수치
  K <- c("# M15 분석 결과 요약 (내부용, 수치는 results/m15의 결과 파일에서 자동 생성)", "",
         HDR_KO, "",
         sprintf("1. 사후 예측 점검(Li 2020 건강인 시험): 주 모델은 misfit 규칙을 %s, 민감도 모델은 %s. 원인은 진단 중(개정 16.7절).", cat_txt(PRIMARY, TRUE), cat_txt(OTHER, TRUE)),
         sprintf("2. 기본 시나리오(S1, 300 mg, 확정 일정, 75 kg 경계, 군당 130명, 층 포함 분석): 참 GMR 0.95 검정력 %s%%, 결합 소비자 위험 0.80에서 %s%%, 1.25에서 %s%%(단일 평가변수 최대 %s%%, %s%%) → 사전 해석 기준 %s.",
                 fmt(p95$pass_pct), fmt(c80$pass_pct), fmt(c125$pass_pct), fmt(max(c80$auc_only_pct, c80$cmax_only_pct)), fmt(max(c125$auc_only_pct, c125$cmax_only_pct)), if (adequate) "충족" else "미충족"),
         sprintf("3. 참 GMR 0.95 검정력을 2%%p 이상 바꾸고 짝 95%% 구간이 0을 포함하지 않는 설계 요인(주 모델, 군당 130명; 집단·용량·일정·경계·분석 모형 칸 수): %s.", sig_txt_ko),
         sprintf("4. 결론 확률 50%% 교차(기본 시나리오): %s.", side_txt(bc, TRUE)),
         c5_ko, "",
         "요약 표: results/m15/summary/m15_summary_table.md(영문, 이탈·사후 분석 표 포함; results/m15/summary/m15_deviations.csv). 적용 가능성: results/m15/summary/m15_applicability.csv. 그림: results/m15/figures/.")
  writeLines(K, file.path(SD, "m15_summary_ko.md"))
  # 가정 목록: 계획서 12절 표를 CSV로
  pl <- readLines(proj_path("plans", "M15_analysis_plan_20261020.md"), warn = FALSE)
  i0 <- grep("^## 12\\. Assumptions", pl); i1 <- grep("^## 13\\.", pl); tl <- pl[i0:i1]; tl <- tl[grepl("^\\| B[0-9]+ ", tl)]
  AS <- rbindlist(lapply(strsplit(sub("^\\| ", "", sub(" \\|$", "", tl)), " \\| "), function(x) as.data.table(as.list(x))))
  setnames(AS, c("ID", "Assumption", "Value used", "Basis", "Impact evidence", "Action before submission")); fwrite(AS, file.path(SD, "m15_assumptions_register.csv"))
  if (nzchar(COND)) writeLines(c("# Notice on the M15 results", "", COND, "",
                                 "This applies to every result under results/m15 except the validation folder: the trial, truth, individual, population, weight, sample-size and summary outputs and the figures.",
                                 "Result tables written by the summary, sample-size and weight steps carry it in a column named conditionality."), file.path(RD, "NOTICE_conditional.md"))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE)
  en_files <- list.files(RD, "\\.(csv|md)$", recursive = TRUE, full.names = TRUE); en_files <- en_files[!grepl("_ko\\.md$", en_files)]   # results/m15의 모든 영문 csv·md
  for (f in en_files) check_english(f)
  append_run_log(lf, sprintf("English text check: %d files", length(en_files)))
  print(TAB); cat(C, sep = "\n")
  append_run_log(lf, "done")
}
