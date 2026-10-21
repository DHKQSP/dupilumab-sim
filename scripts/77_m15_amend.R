#!/usr/bin/env Rscript
# 77_m15_amend.R — 개정 1(계획서 16절, D-069): 일치 점검, 추출, 진단, 교차, 부록, robustness. 등록 산출물은 읽기만 하고 results/m15/amend1/에 쓴다.
# 모드
#   mechcheck [model] [pop]   mechsuppx 일치 점검(16.3절 6): (i) Model 1 B0, 판 1~16 = 저장 mechsupp(모든 열), 묶음별 개수 행, versions.csv 1~16행;
#                             (ii) 두 모델, 세 일정, 항등 판(vi 1) = 저장 mech의 vi 1 행(모든 열). 있는 시험만 비교. 불일치면 종료 코드 3
#   mechx                     기전 OC 확장 표(16.3절 6): mech + mechsuppx 판, 참값 + 참값 확장 → amend1/mech_pass_ext.csv, mech_crossing_ext.csv, mech_crossing_diff_ext.csv
#                             (등록 부트스트랩 색인; 저장 교차와 일치 점검)
#   weightcheck               체중 효과 점검(16.4절 1): weight_ratio.csv 모든 행 재계산(1e-6) → amend1/weight_check.csv; 진단 모델 체중 효과 → amend1/weight_ratio_diagnostic.csv
#   extract [cores]           P1 추출(16.3절 1~4, 16.4절 2·4): amend1/p1/ — sampling_schedule.csv(검정력, lambda_z 성공률, 신뢰 말단상 몫, 오차, AUROC, 짝 대상자 부트스트랩 차이),
#                             crossings.csv, coverage.csv, ada.csv, cell_detail.csv, population_contrast.csv
#   ppcstored                 사후 예측 점검 표(16.7절 1·2, 저장 반복만): amend1/ppc/arm_table.csv, seven_arm_bias.csv
# 사용법: nice -n 5 Rscript scripts/77_m15_amend.R <mode> ...
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("mechcheck", "mechx", "weightcheck", "extract", "ppcstored")
if (!length(args) || !args[1] %in% MODES) stop("사용법: Rscript scripts/77_m15_amend.R <", paste(MODES, collapse = "|"), "> ...")
mode <- args[1]
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
AM <- CFG$amendment_20261021; stopifnot(!is.null(AM))
RD <- proj_path(CFG$output_dir); AD <- file.path(RD, "amend1"); dir.create(AD, showWarnings = FALSE, recursive = TRUE)
TD <- file.path(RD, "trials")

# 두 표가 열·행 모두 같은지(키로 정렬; 숫자는 상대 1e-12, 나머지는 동일). 반환: 최대 상대 차이(행·열 불일치면 Inf)
same_tab <- function(a, b, keys) {
  if (nrow(a) != nrow(b) || !setequal(names(a), names(b))) return(Inf)
  setcolorder(b, names(a)); setkeyv(a, keys); setkeyv(b, keys)
  if (!isTRUE(all.equal(a[, ..keys], b[, ..keys], check.attributes = FALSE))) return(Inf)
  max(vapply(setdiff(names(a), keys), function(cc) { x <- a[[cc]]; y <- b[[cc]]
    if (!identical(is.na(x), is.na(y))) return(Inf)
    z <- !is.na(x); if (!any(z)) return(0)
    if (is.numeric(x) && is.numeric(y)) max(abs(x[z] - y[z]) / pmax(abs(y[z]), 1e-300)) else if (identical(as.character(x[z]), as.character(y[z]))) 0 else Inf }, 1))
}
rd_files <- function(d, pat) { f <- sort(list.files(d, pattern = pat, full.names = TRUE)); if (!length(f)) return(NULL); setNames(lapply(f, fread), basename(f)) }

if (mode == "mechcheck") {
  X <- AM$mechsuppx; ms <- if (length(args) >= 2) args[2] else unlist(X$models); ps <- if (length(args) >= 3) args[3] else unlist(X$populations)
  lf <- start_run_log("m15_amend_mechcheck", master_seed = MS, run_mode = "final", extra = list(models = ms, populations = ps))
  R <- list()
  for (m in ms) for (pv in ps) {
    cell <- sprintf("%s_%s_%d", m, pv, as.integer(X$dose_mg)); dx <- file.path(TD, "mechsuppx", cell); if (!dir.exists(dx)) next
    bx <- rd_files(dx, "^mechsuppx_\\d+_\\d+\\.csv\\.gz$"); cx <- rd_files(dx, "^cnt_mechsuppx_\\d+_\\d+\\.csv\\.gz$"); if (is.null(bx)) next
    BX <- rbindlist(bx); trials <- sort(unique(BX$trial))
    # (ii) 항등 판: 세 일정, 저장 mech의 vi 1 행
    dm <- file.path(TD, "mech", cell); BM <- rbindlist(rd_files(dm, "^mech_\\d+_\\d+\\.csv\\.gz$"))
    a <- BX[vi == 1]; b <- BM[vi == 1 & trial %in% trials]
    R[[length(R) + 1]] <- data.table(cell = cell, check = "identity rows (vi 1) vs mech, all schedules", n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), copy(b), c("schedule", "trial")))
    # (i) Model 1 B0: 판 1~16 = 저장 mechsupp
    if (m == CFG$models$primary) {
      ds <- file.path(TD, "mechsupp", cell); BS <- rbindlist(rd_files(ds, "^mechsupp_\\d+_\\d+\\.csv\\.gz$"))
      a <- BX[schedule == "B0" & vi <= 16]; b <- BS[trial %in% trials]
      R[[length(R) + 1]] <- data.table(cell = cell, check = "B0 versions 1 to 16 vs mechsupp", n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), copy(b), c("schedule", "trial", "vi")))
      cs <- rd_files(ds, "^cnt_mechsupp_\\d+_\\d+\\.csv\\.gz$")
      for (nm in names(cx)) { s_ <- sub("^cnt_mechsuppx_", "cnt_mechsupp_", nm); if (is.null(cs[[s_]])) { R[[length(R) + 1]] <- data.table(cell = cell, check = paste("counts", nm), n_rows = NA_integer_, n_stored = 0L, max_rel_diff = Inf); next }
        a <- cx[[nm]][schedule == "B0" & vi <= 16]; b <- copy(cs[[s_]])
        R[[length(R) + 1]] <- data.table(cell = cell, check = paste("counts", nm), n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), b, c("schedule", "vi", "n_rand", "model"))) }
      vx <- fread(file.path(dx, "versions.csv")); vs <- fread(file.path(ds, "versions.csv"))
      R[[length(R) + 1]] <- data.table(cell = cell, check = "versions.csv rows 1 to 16", n_rows = nrow(vx[vi <= 16]), n_stored = nrow(vs), max_rel_diff = same_tab(copy(vx[vi <= 16]), copy(vs), "vi"))
    }
  }
  CK <- rbindlist(R); stopifnot(nrow(CK) > 0); CK[, pass := is.finite(max_rel_diff) & max_rel_diff <= as.numeric(AM$identity_check$csv_relative) & n_rows == n_stored]
  fwrite(CK, file.path(AD, "mechsuppx_identity_check.csv")); print(CK)
  append_run_log(lf, sprintf("mechcheck: %d checks, all pass %s", nrow(CK), all(CK$pass)))
  if (!all(CK$pass)) quit(status = 3)
}

NBASE <- as.integer(CFG$mechanisms$base_n_randomized_per_arm); B_BOOT <- 2000L; PRIMARY <- CFG$models$primary
COND <- m15_conditional_label(RD, PRIMARY)

# ---------------------------------------------------------------- mechx: 기전 OC 확장(16.3절 6) ----------------------------------------------------------------
if (mode == "mechx") {
  lf <- start_run_log("m15_amend_mechx", master_seed = MS, run_mode = "final", extra = list(B = B_BOOT))
  X <- AM$mechsuppx; DOSE <- as.numeric(X$dose_mg)
  ck <- fread(file.path(AD, "mechsuppx_identity_check.csv")); stopifnot(nrow(ck) > 0, all(ck$pass))                    # 일치 점검 통과가 전제
  TG <- rbind(fread(file.path(RD, "truth", "truth_gmr.csv")), fread(file.path(RD, "truth", "truth_gmr_ext.csv")), fill = TRUE)
  stopifnot(!anyDuplicated(TG[, .(pk_model, population, dose_mg, mechanism, multiplier, metric)]))
  rdx <- function(kind, cell) { d <- file.path(TD, kind, cell); if (!dir.exists(d)) return(NULL)
    b <- rbindlist(rd_files(d, sprintf("^%s_\\d+_\\d+\\.csv\\.gz$", kind))); cn <- rbindlist(rd_files(d, sprintf("^cnt_%s_\\d+_\\d+\\.csv\\.gz$", kind)))
    nt <- if (kind == "mech") as.integer(CFG$mechanisms$n_trials) else as.integer(X$n_trials)
    stopifnot(all(b[, .N, by = .(schedule, vi)]$N == nt), all(b[, uniqueN(trial), by = .(schedule, vi)]$V1 == nt))     # 모든 시험, 중복 없음
    list(base = b, cnt = cn, v = fread(file.path(d, "versions.csv"))) }
  wil <- function(k, n) { w <- wilson_ci(k, n); list(pct = w$est, lo = w$lo, hi = w$hi) }
  BIDX <- m15_boot_index(as.integer(CFG$mechanisms$n_trials), B_BOOT, derive_seed(MS, "crossing_boot"))
  cross_one <- function(mult, p) {                                  # scripts/75와 같은 규칙: 1.00에서 바깥쪽 첫 50% 미만 교차(log 배율 선형 보간)
    o <- order(abs(log(mult))); mult <- mult[o]; p <- p[o]
    if (p[1] < 0.5) return(1)
    for (k in 2:length(mult)) if (p[k] < 0.5) return(exp(log(mult[k - 1]) + (0.5 - p[k - 1]) / (p[k] - p[k - 1]) * (log(mult[k]) - log(mult[k - 1]))))
    NA_real_
  }
  MPX <- list(); CRS <- list(); CRB <- list()
  for (m in unlist(X$models)) for (pv in unlist(X$populations)) {
    cell <- sprintf("%s_%s_%d", m, pv, as.integer(DOSE)); rm_ <- rdx("mech", cell); rx <- rdx("mechsuppx", cell); stopifnot(!is.null(rm_), !is.null(rx))
    MPX[[cell]] <- merge(rx$cnt[, .(n_trials = sum(n_trials), n_pass = sum(n_pass), n_pass_auc = sum(n_pass_auc), n_pass_cmax = sum(n_pass_cmax)), by = .(schedule, vi, n_rand, model)], rx$v, by = "vi")[,
      `:=`(pk_model = m, population = pv, dose_mg = DOSE, set = "mechsuppx")]
    for (sn in unique(rm_$base$schedule)) {
      b <- rm_$base[schedule == sn, .(trial, vi, pass = pass_M1)][rm_$v, on = "vi"]
      bs <- rx$base[schedule == sn & vi > 1, .(trial, vi, pass = pass_M1)][rx$v, on = "vi", nomatch = NULL]
      b <- rbind(b, bs[, names(b), with = FALSE])
      for (mech in setdiff(unique(b$mechanism), "identity")) for (side in c("down", "up")) {
        bb <- b[mechanism %in% c(mech, "identity") & (if (side == "down") multiplier <= 1 else multiplier >= 1)]
        W <- dcast(bb, trial ~ multiplier, value.var = "pass", fun.aggregate = function(z) as.integer(z[1]))
        mult <- as.numeric(names(W)[-1]); P <- as.matrix(W[, -1]); nt <- nrow(P)
        est <- cross_one(mult, colMeans(P))
        bt <- vapply(BIDX, function(i) { r_ <- match(i, W$trial); r_ <- r_[!is.na(r_)]; cross_one(mult, colMeans(P[r_, , drop = FALSE])) }, 1)
        key <- paste(m, pv, DOSE, sn, mech, side, sep = "|"); CRB[[key]] <- bt
        btc <- ifelse(is.finite(bt), bt, if (side == "down") 0 else Inf); qb <- quantile(btc, c(0.025, 0.975), names = FALSE, type = 7)
        nr_ <- sum(!is.finite(bt)); beyond_outer <- nr_ >= 0.025 * length(bt)
        beyond <- function(q) beyond_outer && ((side == "down" && identical(q, qb[1])) || (side == "up" && identical(q, qb[2])))
        tg <- TG[pk_model == m & population == pv & dose_mg == DOSE & (mechanism == mech | mechanism == "identity") &
                 metric %in% c("AUCinf", "Cmax", paste0("AUClast_dense_", sn))][(if (side == "down") multiplier <= 1 else multiplier >= 1)]
        ig <- function(mt) { z <- tg[metric == mt][order(multiplier)]; if (!is.finite(est) || nrow(z) < 2 || est < min(z$multiplier) || est > max(z$multiplier)) return(NA_real_)
          exp(approx(log(z$multiplier), log(z$gmr), xout = log(est))$y) }
        CRS[[key]] <- data.table(pk_model = m, population = pv, dose_mg = DOSE, schedule = sn, mechanism = mech, side = side, model = "M1", n_rand = NBASE,
                                 grid_end = if (side == "down") min(mult) else max(mult), p_at_grid_end = colMeans(P)[which.max(abs(log(mult)))], reached = is.finite(est),
                                 crossing_multiplier = est, crossing_pct = 100 * (est - 1),
                                 lo_pct = if (!is.finite(est) || beyond(qb[1])) NA_real_ else 100 * (qb[1] - 1), hi_pct = if (!is.finite(est) || beyond(qb[2])) NA_real_ else 100 * (qb[2] - 1),
                                 outer_bound_beyond_grid = beyond(if (side == "down") qb[1] else qb[2]), boot_not_reached_pct = 100 * mean(!is.finite(bt)),
                                 induced_gmr_AUClast = ig(paste0("AUClast_dense_", sn)), induced_gmr_AUCinf = ig("AUCinf"), induced_gmr_Cmax = ig("Cmax"), n_trials = nt)
      }
    }
  }
  MP <- rbindlist(MPX); MP[, c("pass_pct", "pass_lo", "pass_hi") := wil(n_pass, n_trials)]
  setcolorder(MP, c("pk_model", "population", "dose_mg", "set", "schedule", "mechanism", "multiplier", "n_rand", "model")); m15_fwrite_cond(MP, file.path(AD, "mech_pass_ext.csv"), COND)
  CRT <- rbindlist(CRS); CRT[, ckey := paste(pk_model, population, dose_mg, schedule, mechanism, side, sep = "|")]
  # 저장 교차와 일치: 저장 격자 안에서 도달하고 미도달 재표집이 없던 칸은 같아야 한다(그 칸은 새 판이 결과를 바꾸지 않는다)
  ST <- fread(file.path(RD, "summary", "mech_crossing.csv"))[dose_mg == DOSE]; ST[, ckey := paste(pk_model, population, dose_mg, schedule, mechanism, side, sep = "|")]
  cmp <- merge(ST[reached == TRUE & boot_not_reached_pct == 0, .(ckey, s_est = crossing_multiplier, s_lo = lo_pct, s_hi = hi_pct)], CRT[, .(ckey, est = crossing_multiplier, lo = lo_pct, hi = hi_pct)], by = "ckey")
  cmp[, max_abs_diff := pmax(abs(est - s_est), abs(lo - s_lo), abs(hi - s_hi), na.rm = TRUE)]
  fwrite(cmp, file.path(AD, "mech_crossing_ext_check.csv")); append_run_log(lf, sprintf("stored crossings compared: %d, max abs diff %g", nrow(cmp), max(cmp$max_abs_diff)))
  if (!nrow(cmp) || max(cmp$max_abs_diff) > 1e-9) stop("저장 교차와 불일치")
  CRT[, source := fifelse(ckey %in% cmp$ckey, "unchanged from stored (reached within the stored grid)", "extended grid (mechsuppx)")]
  m15_fwrite_cond(CRT[, !"ckey"], file.path(AD, "mech_crossing_ext.csv"), COND)
  cdiff <- function(key_a, key_b, ct_lab) { a <- CRB[[key_a]]; b <- CRB[[key_b]]; if (is.null(a) || is.null(b)) return(NULL)      # scripts/75와 같은 규칙
    sd_ <- if (grepl("\\|down$", key_a)) -Inf else Inf
    ac <- ifelse(is.finite(a), a, sd_); bc_ <- ifelse(is.finite(b), b, sd_)
    d <- 100 * (ac - bc_); d[!is.finite(a) & !is.finite(b)] <- NA_real_
    ea <- CRT$crossing_multiplier[match(key_a, CRT$ckey)]; eb <- CRT$crossing_multiplier[match(key_b, CRT$ckey)]
    na_ <- 100 * mean(is.na(d)); q_ <- if (na_ <= 2.5) quantile(d, c(0.025, 0.975), names = FALSE, na.rm = TRUE) else c(NA_real_, NA_real_)
    q_[!is.finite(q_)] <- NA_real_
    data.table(contrast = ct_lab, a = key_a, b = key_b, diff_pp = 100 * (ea - eb), lo_pp = if (is.finite(ea) && is.finite(eb)) q_[1] else NA_real_,
               hi_pp = if (is.finite(ea) && is.finite(eb)) q_[2] else NA_real_, boot_not_reached_pct = 100 * mean(!is.finite(a) | !is.finite(b)), boot_both_not_reached_pct = na_) }
  CD <- rbindlist(lapply(seq_len(nrow(CRT)), function(i) { r <- CRT[i]; k <- function(pop = r$population, sn = r$schedule) paste(r$pk_model, pop, r$dose_mg, sn, r$mechanism, r$side, sep = "|")
    rbind(if (r$population == "S2") cdiff(k(), k(pop = "S1"), "population S2 - S1"), if (r$schedule != "B0") cdiff(k(), k(sn = "B0"), sprintf("schedule %s - B0", r$schedule))) }))
  m15_fwrite_cond(CD, file.path(AD, "mech_crossing_diff_ext.csv"), COND)
  print(CRT[, .(pk_model, population, schedule, mechanism, side, reached, crossing_pct, lo_pct, hi_pct, p_at_grid_end, grid_end)])
  append_run_log(lf, sprintf("done: %d crossings, %d differences", nrow(CRT), nrow(CD)))
}

# ---------------------------------------------------------------- weightcheck: 체중 효과 점검(16.4절 1) ----------------------------------------------------------------
if (mode == "weightcheck") {
  lf <- start_run_log("m15_amend_weightcheck", master_seed = MS, run_mode = "final")
  WR <- fread(file.path(RD, "weight", "weight_ratio.csv")); TOL <- as.numeric(AM$weight_check$tolerance); REF <- as.numeric(CFG$weight_effect$reference_kg)
  # 독립 계산: 전형 대상자(η = 0) AUC0-inf = 4,000일까지 ODE 누적 AUC, Cmax = 0.005일 격자 최댓값의 세 점 포물선 보정(등록 함수는 0.01일 격자 + 같은 보정)
  indep <- function(p, wts, dose) {
    ip <- individual_params(p, data.table(id = seq_along(wts), WT = wts))
    s <- solve_model(ip, CJ(id = seq_along(wts), time = c(seq(0.005, 60, by = 0.005), 400, 4000)), dose, model_id = p$model_id)
    a <- s[time == 4000, .(id, AUCinf = auc)]
    cm <- s[time <= 60][order(id, time)][, { i <- which.max(C); v <- C[i]
      if (i > 1 && i < .N) { tt <- time[(i - 1):(i + 1)]; cc <- C[(i - 1):(i + 1)]; A <- cbind(tt^2, tt, 1); co <- solve(A, cc); if (co[1] < 0) v <- co[3] - co[2]^2 / (4 * co[1]) }
      .(Cmax = v) }, by = id]
    merge(a, cm, by = "id")[order(id)] }
  CK <- rbindlist(lapply(split(WR, by = c("pk_model", "dose_mg", "population")), function(z) {
    p <- load_params(z$pk_model[1]); d <- z$dose_mg[1]; w <- sort(unique(z$WT)); I <- indep(p, c(w, REF), d); r0 <- I[.N]
    T_ <- m15_typical_metrics(p, c(w, REF), d, last_day = max(as.numeric(CFG$schedules$B0))); t0 <- T_[.N]
    rbindlist(lapply(seq_len(nrow(z)), function(i) { j <- match(z$WT[i], w); mt <- z$metric[i]
      rc <- switch(mt, AUCinf = I$AUCinf[j] / r0$AUCinf, Cmax = I$Cmax[j] / r0$Cmax, AUClast_dense_B0 = T_$AUClast_dense[j] / t0$AUClast_dense)
      data.table(pk_model = z$pk_model[i], dose_mg = d, population = z$population[i], prob = z$prob[i], WT = z$WT[i], metric = mt, theta_WT = p$cov$theta_WT,
                 stored_typical_ratio = z$typical_ratio[i], recomputed = rc, method = if (mt == "AUClast_dense_B0") "registered function" else "independent solver call",
                 rel_diff = abs(rc / z$typical_ratio[i] - 1)) })) }))
  CK[, pass := rel_diff <= TOL]
  m15_fwrite_cond(CK, file.path(AD, "weight_check.csv"), COND)
  # 진단 모델(16.4절 1): 공변량 모델(BMI 인자 1; 평균 키 BMI; 조건부 중앙 키 BMI), 기본 모델. 300 mg, 인종 인자 1. BMI_ref는 고정 키에서 비에 영향 없음
  WB <- AM$weight_check$covariate_model_bmi; DOSE <- 300; PCT <- as.numeric(unlist(CFG$weight_effect$percentiles)) / 100
  WQ <- fread(file.path(RD, "populations", "population_quantiles.csv"))
  POPS <- CFG$populations
  tm_bmi <- function(p, wts, bmi, dose) {                            # m15_typical_metrics와 같은 풀이(0.01일 격자, 400/4,000일), eta 표에 BMI 열
    pp <- p; ip <- individual_params(pp, data.table(id = seq_along(wts), WT = wts, BMI = bmi, ada = 0))
    s <- solve_model(ip, CJ(id = seq_along(wts), time = c(seq(0.01, 60, by = 0.01), 400, 4000)), dose, model_id = pp$model_id)
    a <- s[time == 4000, .(id, AUCinf = auc)]
    cm <- s[time <= 60][order(id, time)][, { i <- which.max(C); v <- C[i]
      if (i > 1 && i < .N) { tt <- time[(i - 1):(i + 1)]; cc <- C[(i - 1):(i + 1)]; A <- cbind(tt^2, tt, 1); co <- solve(A, cc); if (co[1] < 0) v <- co[3] - co[2]^2 / (4 * co[1]) }
      .(Cmax = v) }, by = id]
    merge(a, cm, by = "id")[order(id)] }
  DG <- rbindlist(lapply(c("S1", "S2"), function(pv) {
    w <- WQ[population == pv][match(PCT, round(prob, 6)), WT]; stopifnot(length(w) == length(PCT), !anyNA(w)); wts <- c(w, REF)
    sub_ <- m15_subjects(200000L, POPS[[pv]], CFG$inclusion, as.numeric(POPS$weight_height_rho), as.numeric(POPS$sex_ratio_male), derive_seed(MS, "popq"))
    hcond <- vapply(wts, function(x) median(sub_$HT[abs(sub_$WT - x) <= 1]), 1)                  # 조건부 중앙 키(코풀라, 200,000명, ±1 kg)
    hmean <- as.numeric(WB$height_cm[[pv]])
    pc <- load_params("k2020cov"); pc$cov$ke_bmi_exp <- as.numeric(WB$exponent); pc$cov$BMI_ref <- as.numeric(WB$bmi_ref_kg_m2)
    v <- list(list(lab = "covariate model, BMI factor 1", p = load_params("k2020cov"), bmi = rep(NA_real_, length(wts))),
              list(lab = sprintf("covariate model, ke x (BMI/%g)^%g, mean height %g cm", WB$bmi_ref_kg_m2, WB$exponent, hmean), p = pc, bmi = wts / (hmean / 100)^2),
              list(lab = sprintf("covariate model, ke x (BMI/%g)^%g, conditional median height", WB$bmi_ref_kg_m2, WB$exponent), p = pc, bmi = wts / (hcond / 100)^2),
              list(lab = "base model", p = load_params("k2020base"), bmi = rep(NA_real_, length(wts))),
              list(lab = "Model 1 (registered)", p = load_params("k2020"), bmi = rep(NA_real_, length(wts))),
              list(lab = "Kovalenko 2016 (registered)", p = load_params("k2016"), bmi = rep(NA_real_, length(wts))))
    rbindlist(lapply(v, function(x) { b_ <- if (all(is.na(x$bmi))) rep(25, length(wts)) else x$bmi; r <- tm_bmi(x$p, wts, b_, DOSE); r0 <- r[.N]
      data.table(population = pv, variant = x$lab, prob = PCT, WT = w, height_cm = if (grepl("conditional", x$lab)) hcond[seq_along(w)] else if (grepl("mean height", x$lab)) hmean else NA_real_,
                 ratio_AUCinf = r$AUCinf[seq_along(w)] / r0$AUCinf, ratio_Cmax = r$Cmax[seq_along(w)] / r0$Cmax,
                 role = if (grepl("registered", x$lab)) "registered model (same solver, for comparison)" else "diagnostic model, not registered") })) }))
  DG[, dose_mg := DOSE]
  m15_fwrite_cond(DG, file.path(AD, "weight_ratio_diagnostic.csv"), COND)
  print(CK[metric == "AUCinf" & dose_mg == 300 & prob %in% c(0.05, 0.95), .(pk_model, population, prob, WT, stored_typical_ratio, recomputed, rel_diff)]); print(DG[prob %in% c(0.05, 0.95)])
  append_run_log(lf, sprintf("weight check: %d rows, all pass %s, max rel diff %g", nrow(CK), all(CK$pass), max(CK$rel_diff)))
  if (!all(CK$pass)) quit(status = 3)
}

# ---------------------------------------------------------------- extract: P1 표(16.3절 1~4, 16.4절 2·4) ----------------------------------------------------------------
if (mode == "extract") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
  lf <- start_run_log("m15_amend_extract", master_seed = MS, run_mode = "final", extra = list(cores = cores))
  PD1 <- file.path(AD, "p1"); dir.create(PD1, showWarnings = FALSE)
  SD <- file.path(RD, "summary"); DOSE <- 300; G <- 0.95; MODS <- c(CFG$models$primary, CFG$models$sensitivity); SCH <- names(CFG$schedules)
  OCP <- fread(file.path(SD, "oc_pass.csv")); PDF <- fread(file.path(SD, "oc_paired_diff.csv"), na.strings = c("", "NA"))
  AUR <- fread(file.path(SD, "auroc.csv")); AUD <- fread(file.path(SD, "auroc_diff.csv"))
  EXT <- fread(file.path(RD, "individual", "extrapolation_lambda_z.csv")); BP <- fread(file.path(RD, "individual", "bias_precision.csv"))
  SSI <- fread(file.path(RD, "individual", "ss_inputs_m15.csv")); WM <- fread(file.path(RD, "samplesize", "weight_moments.csv"))
  CRX <- fread(file.path(AD, "mech_crossing_ext.csv")); CDX <- fread(file.path(AD, "mech_crossing_diff_ext.csv"))
  gsel <- function(z) z[abs(z$g - G) < 1e-9]
  # 1) 짝 대상자 부트스트랩(16.3절 1): 2,000회, 20,000명, 모든 일정·두 모델에 같은 색인
  RDS <- readRDS(file.path(RD, "individual", "individual_subjects.rds")); RDS <- RDS[RDS$dose_mg == DOSE]
  NB <- as.integer(AM$sampling$paired_bootstrap$n_resamples); NS <- as.integer(CFG$individual$n_subj)
  BI <- with_seed(derive_seed(MS, "amend_boot"), lapply(seq_len(NB), function(b) sample.int(NS, NS, replace = TRUE)))
  stat_sub <- function(x) {                                         # 등록 정의(scripts/72 summ_schedule)와 같은 집합
    ok_last <- is.finite(x$AUClast) & x$AUClast > 0; e_last <- 100 * (x$AUClast / x$AUClast_true - 1); e_last[!ok_last] <- NA
    okA <- x$lambda_ok & is.finite(x$AUCINF_obs) & x$reliable; eA <- 100 * (x$AUCINF_obs / x$AUCinf_true - 1); eA[!okA] <- NA
    list(success = x$lambda_ok & x$pct_extrap < 20, reliable = x$reliable, e_last = e_last, eA = eA) }
  boot_cell <- function(m, pv) {
    S <- lapply(SCH, function(sn) { x <- RDS[RDS$pk_model == m & RDS$population == pv & RDS$schedule == sn][order(id)]; stopifnot(nrow(x) == NS); stat_sub(x) }); names(S) <- SCH
    f_ <- function(st, i) { el <- st$e_last[i]; ea <- st$eA[i]
      c(success_pct = 100 * mean(st$success[i], na.rm = TRUE), reliable_pct = 100 * mean(st$reliable[i]), bias_auclast_pct = median(el, na.rm = TRUE),
        rmse_log_auclast = sqrt(mean(log(1 + el / 100)^2, na.rm = TRUE)), bias_aucinfA_pct = median(ea, na.rm = TRUE), rmse_log_aucinfA = sqrt(mean(log(1 + ea / 100)^2, na.rm = TRUE))) }
    est <- lapply(S, function(st) f_(st, seq_len(NS)))
    bt <- parallel::mclapply(BI, function(i) lapply(S, function(st) f_(st, i)), mc.cores = cores)
    rbindlist(lapply(setdiff(SCH, "B0"), function(sn) { D <- t(vapply(bt, function(b) b[[sn]] - b$B0, est$B0))
      data.table(pk_model = m, population = pv, schedule = sn, metric = names(est$B0), value = est[[sn]], value_B0 = est$B0, diff = est[[sn]] - est$B0,
                 lo = apply(D, 2, quantile, 0.025, names = FALSE), hi = apply(D, 2, quantile, 0.975, names = FALSE)) })) }
  BOOT <- rbindlist(lapply(MODS, function(m) rbindlist(lapply(c("S1", "S2"), function(pv) boot_cell(m, pv)))))
  # 등록 저장 값과 일치(점추정): extrapolation_lambda_z.csv, bias_precision.csv
  chk <- merge(BOOT[metric == "success_pct", .(pk_model, population, schedule, v = value)], EXT[dose_mg == DOSE, .(pk_model, population, schedule, s = success_pct)], by = c("pk_model", "population", "schedule"))
  stopifnot(nrow(chk) > 0, max(abs(chk$v - chk$s)) < 1e-9)
  chk <- merge(BOOT[metric == "bias_auclast_pct", .(pk_model, population, schedule, v = value)], BP[dose_mg == DOSE & endpoint == "AUClast", .(pk_model, population, schedule, s = bias_median_pct)], by = c("pk_model", "population", "schedule"))
  stopifnot(nrow(chk) > 0, max(abs(chk$v - chk$s)) < 1e-9)
  m15_fwrite_cond(BOOT[, dose_mg := DOSE], file.path(PD1, "schedule_paired_bootstrap.csv"), COND)
  # 2) 일정 표(16.3절 1·2)
  SS <- rbindlist(lapply(MODS, function(m) rbindlist(lapply(c("S1", "S2"), function(pv) rbindlist(lapply(SCH, function(sn) {
    o <- gsel(OCP[pk_model == m & population == pv & dose_mg == DOSE & schedule == sn & boundary == "fixed75" & model == "M1" & n_rand == NBASE]); stopifnot(nrow(o) == 1)
    d <- if (sn == "B0") NULL else gsel(PDF[pk_model == m & population == pv & dose_mg == DOSE & contrast == "schedule" & design_a == sn & design_b == "B0" & boundary == "fixed75" & model == "M1"])
    if (!is.null(d)) stopifnot(nrow(d) == 1)
    e <- EXT[pk_model == m & population == pv & dose_mg == DOSE & schedule == sn]; stopifnot(nrow(e) == 1)
    bp <- function(ep, col) { z <- BP[pk_model == m & population == pv & dose_mg == DOSE & schedule == sn & endpoint == ep]; stopifnot(nrow(z) == 1); z[[col]] }
    a <- AUR[pk_model == m & population == pv & dose_mg == DOSE & schedule == sn & statistic == "AUCinf_A" & mechanism == "F" & size_pct == 10 & direction == "pooled"]; stopifnot(nrow(a) == 1)
    ad <- if (sn == "B0") NULL else AUD[contrast == sprintf("schedule %s - B0", sn) & a == paste(m, pv, DOSE, sn, "AUCinf_A", "F", 10, "pooled", sep = "|")]
    if (!is.null(ad)) stopifnot(nrow(ad) == 1)
    data.table(pk_model = m, population = pv, dose_mg = DOSE, schedule = sn, power_pct = o$pass_pct, power_lo = o$pass_lo, power_hi = o$pass_hi,
               power_diff_vs_B0_pp = if (is.null(d)) NA_real_ else d$diff_pp, power_diff_lo = if (is.null(d)) NA_real_ else d$lo_pp, power_diff_hi = if (is.null(d)) NA_real_ else d$hi_pp,
               lambda_z_success_pct = e$success_pct, lambda_z_success_lo = e$success_lo, lambda_z_success_hi = e$success_hi, reliable_terminal_phase_pct = e$reliable_pct,
               auclast_median_error_pct = bp("AUClast", "bias_median_pct"), auclast_rmse_log = bp("AUClast", "rmse_log"),
               aucinf_A_median_error_pct = bp("AUCinf_ruleA", "bias_median_pct"), aucinf_A_rmse_log = bp("AUCinf_ruleA", "rmse_log"), aucinf_A_n = bp("AUCinf_ruleA", "n"),
               aucinf_computable_median_error_pct = bp("AUCinf_computable", "bias_median_pct"),
               auroc_F10_pooled = a$auroc, auroc_lo = a$lo, auroc_hi = a$hi, auroc_diff_vs_B0 = if (is.null(ad)) NA_real_ else ad$diff,
               auroc_diff_lo = if (is.null(ad)) NA_real_ else ad$lo, auroc_diff_hi = if (is.null(ad)) NA_real_ else ad$hi) }))))))
  m15_fwrite_cond(SS, file.path(PD1, "sampling_schedule.csv"), COND)
  # 3) 교차(16.3절 1, 16.4절 4): F(저장 격자), ke·Vmax·ka(확장 격자)
  m15_fwrite_cond(CRX[, .(pk_model, population, dose_mg, schedule, mechanism, side, reached, crossing_pct, lo_pct, hi_pct, outer_bound_beyond_grid, boot_not_reached_pct,
                          grid_end, p_at_grid_end, induced_gmr_AUClast, induced_gmr_AUCinf, induced_gmr_Cmax, source)], file.path(PD1, "crossings.csv"), COND)
  # 4) 채혈 커버리지(16.3절 3): Day 57·64·71 세 정의 + lambda_z 성공률, 신뢰 말단상 몫
  CV <- fread(file.path(AD, "coverage_quantifiable.csv"))[, conditionality := NULL]
  CV <- merge(CV, EXT[dose_mg == DOSE, .(pk_model, population, schedule, lambda_z_success_pct = success_pct, reliable_terminal_phase_pct = reliable_pct)], by = c("pk_model", "population", "schedule"))
  m15_fwrite_cond(CV, file.path(PD1, "coverage.csv"), COND)
  AD_ <- fread(file.path(AD, "ada_concentrations.csv"))[, conditionality := NULL]; setcolorder(AD_, c("pk_model", "population", "dose_mg", "study_day", "planned"))
  m15_fwrite_cond(AD_, file.path(PD1, "ada.csv"), COND)
  # 5) 칸 상세(16.4절 2): 검정력(Wilson; 모의 오차만), 총 CV, 분석 모형 잔차 CV, 위 층 몫·평가 가능 117명 중 기대 인원
  NE <- as.integer(AM$cell_detail$n_evaluable_per_arm)
  CD <- rbindlist(lapply(MODS, function(m) rbindlist(lapply(c("S1", "S2"), function(pv) rbindlist(lapply(c("fixed75", "median"), function(bd) rbindlist(lapply(c("M1", "M0"), function(am) {
    o <- gsel(OCP[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == bd & model == am & n_rand == NBASE]); stopifnot(nrow(o) == 1)
    si <- SSI[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0"]; w <- WM[population == pv & boundary == bd]; stopifnot(nrow(si) == 1, nrow(w) == 1)
    s2a <- si$sd_auc^2 - if (am == "M1") si$beta_auc^2 * w$var_between else 0; s2c <- si$sd_cmax^2 - if (am == "M1") si$beta_cmax^2 * w$var_between else 0
    data.table(pk_model = m, population = pv, dose_mg = DOSE, schedule = "B0", boundary = bd, boundary_kg = w$boundary_kg, analysis_model = am, n_rand = NBASE, true_gmr = G,
               power_pct = o$pass_pct, power_mc_lo = o$pass_lo, power_mc_hi = o$pass_hi, interval_note = "Wilson 95% interval: simulation error of 1,000 studies only",
               cv_auclast_total_pct = si$cv_auc_pct, cv_cmax_total_pct = si$cv_cmax_pct, cv_auclast_residual_pct = 100 * sqrt(exp(s2a) - 1), cv_cmax_residual_pct = 100 * sqrt(exp(s2c) - 1),
               share_upper_stratum_pct = 100 * w$p_stratum2, expected_upper_of_117 = NE * w$p_stratum2, n_studies = o$n_trials) }))))))))
  m15_fwrite_cond(CD, file.path(PD1, "cell_detail.csv"), COND)
  # 6) 집단 대비(16.4절 4): S2 - S1, 300 mg, B0
  PC <- rbindlist(lapply(MODS, function(m) {
    pw <- gsel(PDF[pk_model == m & dose_mg == DOSE & contrast == "population" & design_a == "S2" & design_b == "S1" & schedule == "B0" & boundary == "fixed75" & model == "M1"]); stopifnot(nrow(pw) == 1)
    o <- function(pv) gsel(OCP[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE])
    a <- function(pv) AUR[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & statistic == "AUCinf_A" & mechanism == "F" & size_pct == 10 & direction == "pooled"]
    ad <- AUD[contrast == "population S2 - S1" & a == paste(m, "S2", DOSE, "B0", "AUCinf_A", "F", 10, "pooled", sep = "|")]; stopifnot(nrow(ad) == 1)
    r1 <- data.table(pk_model = m, quantity = c("power at true GMR 0.95 (%)", "AUROC, AUC0-inf, F 10% (directions pooled)"),
                     S1 = c(o("S1")$pass_pct, a("S1")$auroc), S1_lo = c(o("S1")$pass_lo, a("S1")$lo), S1_hi = c(o("S1")$pass_hi, a("S1")$hi),
                     S2 = c(o("S2")$pass_pct, a("S2")$auroc), S2_lo = c(o("S2")$pass_lo, a("S2")$lo), S2_hi = c(o("S2")$pass_hi, a("S2")$hi),
                     diff_S2_minus_S1 = c(pw$diff_pp, ad$diff), diff_lo = c(pw$lo_pp, ad$lo), diff_hi = c(pw$hi_pp, ad$hi), interval = c("Wilson; paired Tango", "cluster bootstrap (shared indices)"))
    r2 <- rbindlist(lapply(c("F", "ke", "Vmax"), function(mech) rbindlist(lapply(c("down", "up"), function(sd_) {
      g1 <- CRX[pk_model == m & population == "S1" & schedule == "B0" & mechanism == mech & side == sd_]; g2 <- CRX[pk_model == m & population == "S2" & schedule == "B0" & mechanism == mech & side == sd_]
      dd <- CDX[contrast == "population S2 - S1" & a == paste(m, "S2", DOSE, "B0", mech, sd_, sep = "|")]; stopifnot(nrow(g1) == 1, nrow(g2) == 1, nrow(dd) == 1)
      data.table(pk_model = m, quantity = sprintf("50%% crossing, %s %s (percent difference)", mech, if (sd_ == "down") "lower" else "upper"),
                 S1 = g1$crossing_pct, S1_lo = g1$lo_pct, S1_hi = g1$hi_pct, S2 = g2$crossing_pct, S2_lo = g2$lo_pct, S2_hi = g2$hi_pct,
                 diff_S2_minus_S1 = dd$diff_pp, diff_lo = dd$lo_pp, diff_hi = dd$hi_pp, interval = "study bootstrap (shared indices)") })))) 
    rbind(r1, r2) }))
  m15_fwrite_cond(PC[, `:=`(dose_mg = DOSE, schedule = "B0")], file.path(PD1, "population_contrast.csv"), COND)
  print(SS[, .(pk_model, population, schedule, power_pct, power_diff_vs_B0_pp, lambda_z_success_pct, reliable_terminal_phase_pct, auclast_median_error_pct, auroc_F10_pooled)])
  print(BOOT); print(CD); print(PC)
  append_run_log(lf, "done")
}

# ---------------------------------------------------------------- ppcstored: 사후 예측 점검 표(16.7절 1·2, 저장 반복) ----------------------------------------------------------------
if (mode == "ppcstored") {
  lf <- start_run_log("m15_amend_ppcstored", master_seed = MS, run_mode = "final")
  PD2 <- file.path(AD, "ppc"); dir.create(PD2, showWarnings = FALSE)
  V <- CFG$validation; REP <- fread(file.path(RD, "validation", "validation_reps.csv.gz")); L <- AM$labels$in_sample
  ARMS <- fread(file.path(RD, "validation", "validation_arms.csv"))
  q_ <- function(x, p) quantile(x, p, names = FALSE, type = 7, na.rm = TRUE)
  ob <- function(aid, e) { o <- V$arms[[aid]]$observed; if (is.null(o) || is.null(o[[paste0(e, "_mean")]])) NA_real_ else as.numeric(o[[paste0(e, "_mean")]]) }
  AT <- rbindlist(lapply(split(REP[wmode %in% c("primary", "wunc")], by = c("pk_model", "arm", "wmode")), function(z) {
    a <- V$arms[[z$arm[1]]]; m <- z$pk_model[1]
    rbindlist(lapply(c("Cmax", "AUClast"), function(e) { x <- z[[paste0(e, "_mean")]]; o <- ob(z$arm[1], e)
      w <- if (is.finite(o)) wilson_ci(sum(x <= o), length(x)) else list(est = NA_real_)
      data.table(pk_model = m, arm = z$arm[1], wmode = z$wmode[1], study = a$study, dose_mg = a$dose_mg, n_subj = a$n_subj, endpoint = e,
                 in_sample_status = L[[m]][[a$study]], observed_mean = o, pred_median = median(x), pi95_lo = q_(x, 0.025), pi95_hi = q_(x, 0.975), pi99_lo = q_(x, 0.005), pi99_hi = q_(x, 0.995),
                 ratio_pred_obs = median(x) / o, pct_rank = w$est, outside_pi95 = if (is.finite(o)) o < q_(x, 0.025) | o > q_(x, 0.975) else NA,
                 outside_pi99 = if (is.finite(o)) o < q_(x, 0.005) | o > q_(x, 0.995) else NA, n_rep = length(x)) })) }))
  AT[, verdict := fifelse(is.na(observed_mean), "no observed value in the repository", fifelse(outside_pi99, "outside the 99% interval", fifelse(outside_pi95, "outside the 95% interval only", "inside the 95% interval")))]
  # 저장 표(95% 예측구간)와 일치
  ck <- merge(AT[, .(pk_model, arm, wmode, metric = paste0(endpoint, "_mean"), pred_median, pi95_lo, pi95_hi)], ARMS[, .(pk_model, arm, wmode, metric, s_med = pred_median, s_lo = pi_lo, s_hi = pi_hi)], by = c("pk_model", "arm", "wmode", "metric"))
  stopifnot(nrow(ck) == nrow(AT), max(abs(c(ck$pred_median - ck$s_med, ck$pi95_lo - ck$s_lo, ck$pi95_hi - ck$s_hi))) < 1e-9)
  setorder(AT, pk_model, wmode, study, arm, endpoint); m15_fwrite_cond(AT, file.path(PD2, "arm_table.csv"), "")
  # 7개 군(요청의 약 80 kg 300 mg): 관측 6개 군 n 가중 평균 AUC0-last 대 반복별 n 가중 예측 평균
  SA <- AM$ppc$seven_arms; arms <- unlist(SA$arms); oa <- arms[vapply(arms, function(aid) is.finite(ob(aid, "AUClast")), TRUE)]
  nn <- vapply(oa, function(aid) as.numeric(V$arms[[aid]]$n_subj), 1); stopifnot(length(oa) == as.integer(SA$observed_arms), sum(nn) == as.numeric(SA$observed_n))
  obs7 <- sum(nn * vapply(oa, ob, 1, e = "AUClast")) / sum(nn)
  B7 <- rbindlist(lapply(c("k2020", "k2016"), function(m) rbindlist(lapply(c("primary", "wunc", "w70", "w86"), function(md) rbindlist(lapply(list(obs6 = oa, all7 = arms), function(set_) {
    z <- REP[pk_model == m & wmode == md & arm %in% set_]; if (!nrow(z)) return(NULL); stopifnot(uniqueN(z$arm) == length(set_))
    pr <- z[, .(v = sum(n * AUClast_mean) / sum(n), wt = sum(n * WT_mean) / sum(n)), by = rep]$v
    wt <- z[, .(wt = sum(n * WT_mean) / sum(n)), by = rep]$wt
    data.table(pk_model = m, wmode = md, arm_set = if (length(set_) == length(oa)) "six arms with observed AUC0-last (n 102)" else "all seven arms (prediction only)", observed_mean_6arms = obs7,
               pred_median = median(pr), ratio_median = median(pr) / obs7, ratio_pi95_lo = q_(pr, 0.025) / obs7, ratio_pi95_hi = q_(pr, 0.975) / obs7,
               ratio_pi99_lo = q_(pr, 0.005) / obs7, ratio_pi99_hi = q_(pr, 0.995) / obs7, pct_rank_observed = 100 * mean(pr <= obs7), realized_mean_weight_kg = mean(wt), n_rep = length(pr)) }))))))
  B7[, interval_note := "prediction intervals of the ratio (not confidence intervals)"]
  m15_fwrite_cond(B7, file.path(PD2, "seven_arm_bias.csv"), "")
  print(AT[wmode == "primary" & is.finite(observed_mean)]); print(B7)
  append_run_log(lf, "done")
}
