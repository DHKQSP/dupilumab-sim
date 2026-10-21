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
#   p1report                  P1 보고(16.12절): robustness 분류(16.10절, Model 1 대 Kovalenko 2016), P1 결과 표 1개, 핵심 결론 5줄 → amend1/p1/
#   appendix                  브리핑 부록 초안(16.6절, 영문, ICH M15 구성, 표 위주) → amend1/appendix/M15_MS_appendix_draft.md, .docx(pandoc markdown-smart); 영문·ASCII 점검
#   postcheck <before_dir>    표기 반영 사후 단계 재실행 점검(16.1절): summary, samplesize, weight, figures의 CSV를 재실행 전 사본과 열마다 비교(숫자 열 동일, 바뀐 표기 열만 허용)
#   ppc                       P2 사후 예측 점검 진단(16.7절 4~7): 대안 매핑 군 표·7개 군 편향, 판정 규칙 오경보율(두 NCA 규칙, 연구 간 인자 CV 0·10%), 검출 확률,
#                             연구별 규칙, 넓힌 한계로 본 관측 판정, 결론 1 고정 서식 → amend1/ppc/
# 사용법: nice -n 5 Rscript scripts/77_m15_amend.R <mode> ...
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("mechcheck", "mechx", "weightcheck", "extract", "ppcstored", "p1report", "appendix", "ppc", "postcheck")
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
  # 독립 계산(solve_model 따로 호출): 전형 대상자(η = 0) AUC0-inf = 등록 정의대로 400일까지 ODE 누적 AUC(남은 양 < 1e-7 x 용량 확인;
  #   transit 구조 모델은 약 2,700일 이후 수치 언더플로로 NaN이 되어 4,000일 값을 쓰지 않는다),
  # Cmax = 등록 정의(0.01일 격자 최댓값의 세 점 포물선 보정, 8.6절). 참고 열: 0.001일 격자 Cmax 비(수치 정밀도)
  TT <- c(seq(0.001, 60, by = 0.001), seq(61, 400, by = 1))
  pmax3 <- function(time, C) { i <- which.max(C); v <- C[i]
    if (i > 1 && i < length(C)) { tt <- time[(i - 1):(i + 1)]; cc <- C[(i - 1):(i + 1)]; A <- cbind(tt^2, tt, 1); co <- solve(A, cc); if (co[1] < 0) v <- co[3] - co[2]^2 / (4 * co[1]) }
    v }
  indep <- function(p, wts, dose) {
    ip <- individual_params(p, data.table(id = seq_along(wts), WT = wts))
    s <- solve_model(ip, CJ(id = seq_along(wts), time = TT), dose, model_id = p$model_id); stopifnot(!anyNA(s$auc), !anyNA(s$C))
    stc <- setdiff(names(s), c("id", "time", "C", "auc")); a <- s[abs(time - 400) < 1e-9][, resid := rowSums(.SD), .SDcols = stc]
    stopifnot(all(a$resid / dose < 1e-7)); a <- a[, .(id, AUCinf = auc)]
    g01 <- s[time <= 60 + 1e-9 & abs(time * 100 - round(time * 100)) < 1e-6][order(id, time)]
    cm <- g01[, .(Cmax = pmax3(time, C)), by = id]; cf <- s[time <= 60 + 1e-9][order(id, time)][, .(Cmax_fine = pmax3(time, C)), by = id]
    merge(merge(a, cm, by = "id"), cf, by = "id")[order(id)] }
  CK <- rbindlist(lapply(split(WR, by = c("pk_model", "dose_mg", "population")), function(z) {
    p <- load_params(z$pk_model[1]); d <- z$dose_mg[1]; w <- sort(unique(z$WT)); I <- indep(p, c(w, REF), d); r0 <- I[.N]
    T_ <- m15_typical_metrics(p, c(w, REF), d, last_day = max(as.numeric(CFG$schedules$B0))); t0 <- T_[.N]
    rbindlist(lapply(seq_len(nrow(z)), function(i) { j <- match(z$WT[i], w); mt <- z$metric[i]
      rc <- switch(mt, AUCinf = I$AUCinf[j] / r0$AUCinf, Cmax = I$Cmax[j] / r0$Cmax, AUClast_dense_B0 = T_$AUClast_dense[j] / t0$AUClast_dense)
      data.table(pk_model = z$pk_model[i], dose_mg = d, population = z$population[i], prob = z$prob[i], WT = z$WT[i], metric = mt, theta_WT = p$cov$theta_WT,
                 stored_typical_ratio = z$typical_ratio[i], recomputed = rc, method = if (mt == "AUClast_dense_B0") "registered function" else "independent solver call",
                 rel_diff = abs(rc / z$typical_ratio[i] - 1), cmax_ratio_fine_grid = if (mt == "Cmax") I$Cmax_fine[j] / r0$Cmax_fine else NA_real_) })) }))
  CK[, pass := rel_diff <= TOL]
  CK[, within_registered_solver_tolerance := rel_diff <= 5e-6]       # 등록 참값 풀이 허용오차(atol 1e-8, rtol 1e-6; D-042: 개인별 상대 차이 최대 5e-6)
  CK[, note := fifelse(pass, "", "registered tolerance 1e-6 not met; the stored value was solved with the registered truth tolerances (D-042), the independent value is unchanged at tighter tolerances")]
  m15_fwrite_cond(CK, file.path(AD, "weight_check.csv"), COND)
  # 진단 모델(16.4절 1): 공변량 모델(BMI 인자 1; 평균 키 BMI; 조건부 중앙 키 BMI), 기본 모델. 300 mg, 인종 인자 1. BMI_ref는 고정 키에서 비에 영향 없음
  WB <- AM$weight_check$covariate_model_bmi; DOSE <- 300; PCT <- as.numeric(unlist(CFG$weight_effect$percentiles)) / 100
  WQ <- fread(file.path(RD, "populations", "population_quantiles.csv"))
  POPS <- CFG$populations
  tm_bmi <- function(p, wts, bmi, dose) {                            # m15_typical_metrics와 같은 풀이(0.01일 격자, 400/4,000일), eta 표에 BMI 열
    pp <- p; ip <- individual_params(pp, data.table(id = seq_along(wts), WT = wts, BMI = bmi, ada = 0))
    s <- solve_model(ip, CJ(id = seq_along(wts), time = c(seq(0.01, 60, by = 0.01), seq(61, 400, by = 1))), dose, model_id = pp$model_id); stopifnot(!anyNA(s$auc))
    stc <- setdiff(names(s), c("id", "time", "C", "auc")); a <- s[abs(time - 400) < 1e-9][, resid := rowSums(.SD), .SDcols = stc]
    stopifnot(all(a$resid / dose < 1e-7)); a <- a[, .(id, AUCinf = auc)]
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
  append_run_log(lf, sprintf("weight check: %d rows, %d not within 1e-6 (max rel diff %g), all within the registered solver tolerance %s", nrow(CK), sum(!CK$pass), max(CK$rel_diff), all(CK$within_registered_solver_tolerance)))
  if (!all(CK$within_registered_solver_tolerance)) quit(status = 3)          # 풀이 허용오차를 넘으면 구현 문제로 보고 멈춘다(1e-6 미충족 행은 보고하고 이탈로 기록)
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
    nr_ <- REP[pk_model == m & wmode == md, max(rep)]; if (!is.finite(nr_)) return(NULL)
    z <- rbind(REP[pk_model == m & wmode == md & arm %in% set_], REP[pk_model == m & wmode == "primary" & arm %in% setdiff(set_, REP[pk_model == m & wmode == md, unique(arm)]) & rep <= nr_])   # 범위만 알려진 군(R668-AS-0907)은 모드와 무관(등록 규칙과 같음)
    stopifnot(uniqueN(z$arm) == length(set_), all(z[, .N, by = arm]$N == nr_))
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

# ---------------------------------------------------------------- 공통 서식 ----------------------------------------------------------------
fmt <- function(x, d = 1) ifelse(is.na(x), "-", formatC(x, format = "f", digits = d, big.mark = ","))
ci_ <- function(x, lo, hi, d = 1, unit = "") ifelse(is.na(x), "-", sprintf("%s%s (%s to %s)", fmt(x, d), unit, fmt(lo, d), fmt(hi, d)))
MLAB <- c(k2020 = "Model 1", k2016 = "Kovalenko 2016"); PLAB <- c(S1 = "global", S2 = "Korean")
SLAB <- c(B0 = "confirmed schedule", RM36_50 = "reduced schedule (Days 36 to 50 removed)", ADD64_71 = "extended schedule (Days 64 and 71 added)")

# ---------------------------------------------------------------- p1report: robustness 분류(16.10절), P1 표, 결론 5줄 ----------------------------------------------------------------
if (mode == "p1report") {
  lf <- start_run_log("m15_amend_p1report", master_seed = MS, run_mode = "final")
  PD1 <- file.path(AD, "p1"); SD <- file.path(RD, "summary"); DOSE <- 300; G <- 0.95
  RB <- AM$robustness; TPP <- as.numeric(RB$thresholds$pp); TAU <- as.numeric(RB$thresholds$auroc); TRM <- as.numeric(RB$thresholds$rmse_log)
  PDF <- fread(file.path(SD, "oc_paired_diff.csv"), na.strings = c("", "NA"))[abs(g - G) < 1e-9 & n_rand == NBASE]
  AUD <- fread(file.path(SD, "auroc_diff.csv")); CDX <- fread(file.path(AD, "mech_crossing_diff_ext.csv")); BOOT <- fread(file.path(PD1, "schedule_paired_bootstrap.csv"))
  OCP <- fread(file.path(SD, "oc_pass.csv")); PAN <- fread(file.path(RD, "samplesize", "power_analytic.csv"))
  sgn <- function(d, lo, hi, thr) ifelse(is.na(d) | is.na(lo) | is.na(hi), NA_integer_, ifelse((lo > 0 | hi < 0) & abs(d) >= thr, as.integer(sign(d)), 0L))
  cls <- function(s1, s2, d1, d2) { if (is.na(s1) || is.na(s2)) return("no interval, descriptive")
    if (s1 != 0 && s1 == s2) return("1 same direction"); if (s1 == 0 && s2 == 0) return("2 no material effect in either model")
    if (s1 != 0 && s2 != 0) return("4 opposite directions")
    other <- if (s1 != 0) d2 else d1; nz <- if (s1 != 0) s1 else s2
    sprintf("3 material in %s only (other model's point estimate %s sign)", if (s1 != 0) "Model 1" else "Kovalenko 2016", if (sign(other) == nz) "same" else "opposite") }
  row_ <- function(id, contrast, metric, cell, x1, x2, thr, unit) { s1 <- sgn(x1[1], x1[2], x1[3], thr); s2 <- sgn(x2[1], x2[2], x2[3], thr)
    data.table(row_id = id, contrast = contrast, metric = metric, cell = cell, unit = unit, threshold = thr,
               m1_diff = x1[1], m1_lo = x1[2], m1_hi = x1[3], m1_sign = s1, k16_diff = x2[1], k16_lo = x2[2], k16_hi = x2[3], k16_sign = s2, class = cls(s1, s2, x1[1], x2[1])) }
  pw_ <- function(m, ct, a, b, pv = "S1", bd = "fixed75", sn = "B0", am = "M1", dose_ = DOSE) {
    z <- PDF[pk_model == m & contrast == ct & design_a == a & design_b == b & dose_mg == (if (ct == "dose") 200 else dose_)]
    if (ct != "population") z <- z[population == pv]; if (ct != "boundary") z <- z[boundary == bd]; if (ct != "schedule") z <- z[schedule == sn]; if (ct != "stratum term") z <- z[model == am]
    if (ct == "dose") z <- PDF[pk_model == m & contrast == "dose" & design_a == "200" & design_b == "300" & population == pv & boundary == bd & schedule == sn & model == am]
    stopifnot(nrow(z) == 1); c(z$diff_pp, z$lo_pp, z$hi_pp) }
  au_ <- function(m, ct, pv, sn) { z <- AUD[contrast == ct & a == paste(m, pv, DOSE, sn, "AUCinf_A", "F", 10, "pooled", sep = "|")]; stopifnot(nrow(z) == 1); c(z$diff, z$lo, z$hi) }
  cr_ <- function(m, ct, pv, sn, mech, sd_) { z <- CDX[contrast == ct & a == paste(m, pv, DOSE, sn, mech, sd_, sep = "|")]; stopifnot(nrow(z) == 1); c(z$diff_pp, z$lo_pp, z$hi_pp) }
  bs_ <- function(m, pv, sn, mt) { z <- BOOT[pk_model == m & population == pv & schedule == sn & metric == mt]; stopifnot(nrow(z) == 1); c(z$diff, z$lo, z$hi) }
  R <- list()
  for (pv in c("S1", "S2")) {
    cell <- sprintf("%s population, 300 mg, 75 kg boundary, stratified analysis, 130 randomized per arm", PLAB[[pv]])
    for (sn in c("ADD64_71", "RM36_50")) { ct <- sprintf("schedule: %s minus confirmed", sub(" \\(.*", "", SLAB[[sn]]))
      R[[length(R) + 1]] <- row_(paste0("sched_power_", sn, "_", pv), ct, "power at true GMR 0.95 (Tango 95% CI)", cell, pw_("k2020", "schedule", sn, "B0", pv), pw_("k2016", "schedule", sn, "B0", pv), TPP, "pp")
      R[[length(R) + 1]] <- row_(paste0("sched_auroc_", sn, "_", pv), ct, "AUROC of AUC0-inf, F 10%, directions pooled (bootstrap 95% CI)", cell, au_("k2020", sprintf("schedule %s - B0", sn), pv, sn), au_("k2016", sprintf("schedule %s - B0", sn), pv, sn), TAU, "AUROC")
      for (sd_ in c("down", "up")) R[[length(R) + 1]] <- row_(paste0("sched_Fcross_", sd_, "_", sn, "_", pv), ct, sprintf("50%% crossing, F %s (bootstrap 95%% CI)", if (sd_ == "down") "lower" else "upper"), cell,
                                                              cr_("k2020", sprintf("schedule %s - B0", sn), pv, sn, "F", sd_), cr_("k2016", sprintf("schedule %s - B0", sn), pv, sn, "F", sd_), TPP, "pp")
      for (mt in c("success_pct", "reliable_pct", "bias_auclast_pct", "rmse_log_auclast")) {
        lab <- c(success_pct = "lambda_z success (AUC0-inf estimable with extrapolation below 20%)", reliable_pct = "reliable-terminal-phase share", bias_auclast_pct = "median percent error of AUC0-last against the dense truth",
                 rmse_log_auclast = "RMSE of the log ratio of AUC0-last against the dense truth")[[mt]]
        R[[length(R) + 1]] <- row_(paste0("sched_", mt, "_", sn, "_", pv), ct, paste(lab, "(paired subject bootstrap 95% CI)"), cell, bs_("k2020", pv, sn, mt), bs_("k2016", pv, sn, mt),
                                   if (mt == "rmse_log_auclast") TRM else TPP, if (mt == "rmse_log_auclast") "log" else "pp") } }
    R[[length(R) + 1]] <- row_(paste0("strat_term_", pv), "stratification: stratified minus unstratified analysis", "power at true GMR 0.95 (Tango 95% CI)", sub(", stratified analysis", "", cell), pw_("k2020", "stratum term", "M1", "M0", pv), pw_("k2016", "stratum term", "M1", "M0", pv), TPP, "pp")
    R[[length(R) + 1]] <- row_(paste0("strat_boundary_", pv), "stratification: median minus 75 kg boundary", "power at true GMR 0.95 (Tango 95% CI)", sub("75 kg boundary, ", "", cell), pw_("k2020", "boundary", "median", "fixed75", pv), pw_("k2016", "boundary", "median", "fixed75", pv), TPP, "pp")
    R[[length(R) + 1]] <- row_(paste0("dose_", pv), "dose: 200 minus 300 mg (limited by the 200 mg presentation; not used in the appendix)", "power at true GMR 0.95 (Tango 95% CI)", sub("300 mg, ", "", cell), pw_("k2020", "dose", "200", "300", pv), pw_("k2016", "dose", "200", "300", pv), TPP, "pp")
  }
  cellP <- "Korean minus global population, 300 mg, confirmed schedule, 75 kg boundary, stratified analysis, 130 randomized per arm"
  R[[length(R) + 1]] <- row_("pop_power", "population: Korean minus global", "power at true GMR 0.95 (Tango 95% CI)", cellP, pw_("k2020", "population", "S2", "S1"), pw_("k2016", "population", "S2", "S1"), TPP, "pp")
  R[[length(R) + 1]] <- row_("pop_auroc", "population: Korean minus global", "AUROC of AUC0-inf, F 10%, directions pooled (bootstrap 95% CI)", cellP, au_("k2020", "population S2 - S1", "S2", "B0"), au_("k2016", "population S2 - S1", "S2", "B0"), TAU, "AUROC")
  for (mech in c("F", "ke", "Vmax")) for (sd_ in c("down", "up")) R[[length(R) + 1]] <- row_(paste0("pop_cross_", mech, "_", sd_), "population: Korean minus global", sprintf("50%% crossing, %s %s (bootstrap 95%% CI)", mech, if (sd_ == "down") "lower" else "upper"), cellP,
                                                                                              cr_("k2020", "population S2 - S1", "S2", "B0", mech, sd_), cr_("k2016", "population S2 - S1", "S2", "B0", mech, sd_), TPP, "pp")
  ROB <- rbindlist(R)
  # 다른 칸의 같은 범주 수(검정력 대비): 같은 대비, 300 mg, 다른 집단·경계·일정·분석 모형
  other_cells <- function(ct, a, b) { z <- PDF[contrast == ct & design_a == a & design_b == b & (dose_mg == DOSE | ct == "dose")]
    k <- c("population", "boundary", "schedule", "model"); k <- intersect(k, names(z)); w <- dcast(z, as.formula(paste(paste(k, collapse = " + "), "~ pk_model")), value.var = c("diff_pp", "lo_pp", "hi_pp"))
    w[, cl := mapply(function(d1, l1, h1, d2, l2, h2) cls(sgn(d1, l1, h1, TPP), sgn(d2, l2, h2, TPP), d1, d2), diff_pp_k2020, lo_pp_k2020, hi_pp_k2020, diff_pp_k2016, lo_pp_k2016, hi_pp_k2016)]; w }
  ROB[, other_cells_same_class := NA_character_]
  for (i in which(grepl("power", ROB$metric))) { r <- ROB[i]
    ct <- if (grepl("^schedule", r$contrast)) "schedule" else if (grepl("stratified minus", r$contrast)) "stratum term" else if (grepl("boundary", r$contrast)) "boundary" else if (grepl("^dose", r$contrast)) "dose" else "population"
    ab <- switch(ct, schedule = c(if (grepl("extended", r$contrast)) "ADD64_71" else "RM36_50", "B0"), `stratum term` = c("M1", "M0"), boundary = c("median", "fixed75"), dose = c("200", "300"), population = c("S2", "S1"))
    w <- other_cells(ct, ab[1], ab[2]); set(ROB, i, "other_cells_same_class", sprintf("%d of %d cells", sum(substr(w$cl, 1, 1) == substr(r$class, 1, 1)), nrow(w))) }
  # 표본 수 행: 130명 모의 검정력의 Wilson 하한 >= 90%(2.1절 1번의 검정력 부분), 분석식 90% 필요 n
  ss <- rbindlist(lapply(c("S1", "S2"), function(pv) rbindlist(lapply(c("k2020", "k2016"), function(m) {
    o <- OCP[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE & abs(g - G) < 1e-9]
    pa <- PAN[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE & abs(true_gmr - G) < 1e-9]; stopifnot(nrow(o) == 1, nrow(pa) == 1)
    data.table(population = pv, pk_model = m, power = o$pass_pct, lo = o$pass_lo, hi = o$pass_hi, meets = o$pass_lo >= 90, n_rand_90 = pa$n_rand_90, n_rand_85 = pa$n_rand_85) }))))
  SSR <- ss[, .(row_id = paste0("samplesize_", population), contrast = "sample size", metric = "simulated power at 130 randomized per arm, lower Wilson 95% bound at least 90%",
                cell = sprintf("%s population, 300 mg, confirmed schedule, 75 kg boundary, stratified analysis", PLAB[population]),
                class = { mm <- meets[pk_model == "k2020"]; kk <- meets[pk_model == "k2016"]; if (mm && kk) "both meet" else if (!mm && !kk) "neither meets" else "one meets" },
                m1_text = sprintf("%s%% (%s to %s); n for 90%%: %d", fmt(power[pk_model == "k2020"]), fmt(lo[pk_model == "k2020"]), fmt(hi[pk_model == "k2020"]), n_rand_90[pk_model == "k2020"]),
                k16_text = sprintf("%s%% (%s to %s); n for 90%%: %d", fmt(power[pk_model == "k2016"]), fmt(lo[pk_model == "k2016"]), fmt(hi[pk_model == "k2016"]), n_rand_90[pk_model == "k2016"])), by = population][, population := NULL]
  ROB[, `:=`(m1_text = ci_(m1_diff, m1_lo, m1_hi, if (unit == "pp") 1 else 3), k16_text = ci_(k16_diff, k16_lo, k16_hi, if (unit == "pp") 1 else 3))]
  ROB <- rbind(ROB, SSR, fill = TRUE)
  ROB[, corrected_mapping := as.character(RB$corrected_mapping_column)]
  m15_fwrite_cond(ROB, file.path(PD1, "robustness_p1.csv"), COND)
  print(ROB[, .(row_id, class, m1_text, k16_text, other_cells_same_class)])
  append_run_log(lf, sprintf("robustness rows: %d; classes %s", nrow(ROB), paste(names(table(substr(ROB$class, 1, 1))), table(substr(ROB$class, 1, 1)), sep = ":", collapse = ", ")))

  # ---------------- P1 결과 표 1개(두 모델 나란히; 값은 amend1/p1/*.csv에서) ----------------
  SS <- fread(file.path(PD1, "sampling_schedule.csv")); CR <- fread(file.path(PD1, "crossings.csv")); CV <- fread(file.path(PD1, "coverage.csv")); AD_ <- fread(file.path(PD1, "ada.csv"))
  CD <- fread(file.path(PD1, "cell_detail.csv")); PC <- fread(file.path(PD1, "population_contrast.csv")); VD <- fread(file.path(AD, "variance_decomposition.csv"))
  WC <- fread(file.path(AD, "weight_check.csv")); WD <- fread(file.path(AD, "weight_ratio_diagnostic.csv"))
  T <- list(); addt <- function(sec, item, v1, v2, src) T[[length(T) + 1]] <<- data.table(section = sec, item = item, model_1 = v1, kovalenko_2016 = v2, source = src)
  bypop <- function(f) paste(sprintf("%s: %s", PLAB, vapply(c("S1", "S2"), f, "")), collapse = "; ")
  for (sn in c("B0", "ADD64_71", "RM36_50")) {
    sv <- function(m, pv) SS[pk_model == m & population == pv & schedule == sn]
    lab <- SLAB[[sn]]
    addt("Sampling", sprintf("Power at true GMR 0.95, %s (%%; Wilson 95%% CI%s)", lab, if (sn == "B0") "" else "; difference from the confirmed schedule in pp, Tango 95% CI"),
         bypop(function(pv) { z <- sv("k2020", pv); paste0(ci_(z$power_pct, z$power_lo, z$power_hi), if (sn != "B0") paste0("; ", ci_(z$power_diff_vs_B0_pp, z$power_diff_lo, z$power_diff_hi)) else "") }),
         bypop(function(pv) { z <- sv("k2016", pv); paste0(ci_(z$power_pct, z$power_lo, z$power_hi), if (sn != "B0") paste0("; ", ci_(z$power_diff_vs_B0_pp, z$power_diff_lo, z$power_diff_hi)) else "") }), "amend1/p1/sampling_schedule.csv")
    addt("Sampling", sprintf("AUC0-inf estimable with extrapolation below 20%% (lambda_z success), %s (%%)", lab),
         bypop(function(pv) { z <- sv("k2020", pv); ci_(z$lambda_z_success_pct, z$lambda_z_success_lo, z$lambda_z_success_hi) }), bypop(function(pv) { z <- sv("k2016", pv); ci_(z$lambda_z_success_pct, z$lambda_z_success_lo, z$lambda_z_success_hi) }), "amend1/p1/sampling_schedule.csv")
    addt("Sampling", sprintf("Reliable-terminal-phase share, %s (%%)", lab), bypop(function(pv) fmt(sv("k2020", pv)$reliable_terminal_phase_pct)), bypop(function(pv) fmt(sv("k2016", pv)$reliable_terminal_phase_pct)), "amend1/p1/sampling_schedule.csv")
    addt("Sampling", sprintf("AUC0-last error against the dense truth, %s (median %%; RMSE of the log ratio)", lab),
         bypop(function(pv) { z <- sv("k2020", pv); sprintf("%s; %s", fmt(z$auclast_median_error_pct, 2), fmt(z$auclast_rmse_log, 3)) }), bypop(function(pv) { z <- sv("k2016", pv); sprintf("%s; %s", fmt(z$auclast_median_error_pct, 2), fmt(z$auclast_rmse_log, 3)) }), "amend1/p1/sampling_schedule.csv")
    addt("Sampling", sprintf("AUC0-inf error, reliable-terminal-phase set, %s (median %%; RMSE of the log ratio)", lab),
         bypop(function(pv) { z <- sv("k2020", pv); sprintf("%s; %s", fmt(z$aucinf_A_median_error_pct, 2), fmt(z$aucinf_A_rmse_log, 3)) }), bypop(function(pv) { z <- sv("k2016", pv); sprintf("%s; %s", fmt(z$aucinf_A_median_error_pct, 2), fmt(z$aucinf_A_rmse_log, 3)) }), "amend1/p1/sampling_schedule.csv")
    addt("Sampling", sprintf("AUROC of AUC0-inf for a 10%% bioavailability difference, directions pooled, %s%s", lab, if (sn == "B0") "" else " (difference from the confirmed schedule)"),
         bypop(function(pv) { z <- sv("k2020", pv); paste0(ci_(z$auroc_F10_pooled, z$auroc_lo, z$auroc_hi, 3), if (sn != "B0") paste0("; ", ci_(z$auroc_diff_vs_B0, z$auroc_diff_lo, z$auroc_diff_hi, 3)) else "") }),
         bypop(function(pv) { z <- sv("k2016", pv); paste0(ci_(z$auroc_F10_pooled, z$auroc_lo, z$auroc_hi, 3), if (sn != "B0") paste0("; ", ci_(z$auroc_diff_vs_B0, z$auroc_diff_lo, z$auroc_diff_hi, 3)) else "") }), "amend1/p1/sampling_schedule.csv")
    cx <- function(m, pv) paste(vapply(c("F", "ke", "Vmax", "ka"), function(mech) paste(vapply(c("down", "up"), function(sd_) { z <- CR[pk_model == m & population == pv & schedule == sn & mechanism == mech & side == sd_]; stopifnot(nrow(z) == 1)
      if (isTRUE(z$reached)) sprintf("%s%%", fmt(z$crossing_pct)) else sprintf("not reached (%s%% at %s%%)", fmt(100 * z$p_at_grid_end), fmt(100 * (z$grid_end - 1), 0)) }, ""), collapse = " / ") |> (\(x) sprintf("%s %s", mech, x))(), ""), collapse = "; ")
    addt("Sampling", sprintf("50%% crossing of the conclusion probability, %s (lower / upper, percent difference)", lab), bypop(function(pv) cx("k2020", pv)), bypop(function(pv) cx("k2016", pv)), "amend1/p1/crossings.csv")
  }
  cvv <- function(m, dd, sn = "B0") bypop(function(pv) { z <- CV[pk_model == m & population == pv & schedule == sn & planned_elapsed_day == dd]; if (!nrow(z)) "-" else sprintf("raw %s, true %s, retained %s", fmt(z$raw_pct), fmt(z$true_pct), fmt(z$used_pct)) })
  addt("Coverage", "Day 57 sample quantifiable, confirmed schedule (% of subjects)", cvv("k2020", 56), cvv("k2016", 56), "amend1/p1/coverage.csv")
  for (dd in c(63, 70)) addt("Coverage", sprintf("Day %d sample quantifiable, extended schedule (%% of subjects)", dd + 1), cvv("k2020", dd, "ADD64_71"), cvv("k2016", dd, "ADD64_71"), "amend1/p1/coverage.csv")
  for (sd_ in unique(AD_$study_day)) addt("Anti-drug antibody sampling", sprintf("Day %d true serum concentration (mg/L): median (5th to 95th percentile); below LLOQ %%", sd_),
    bypop(function(pv) { z <- AD_[pk_model == "k2020" & population == pv & study_day == sd_]; sprintf("%s (%s to %s); %s", fmt(z$true_median, 2), fmt(z$true_p05, 2), fmt(z$true_p95, 2), fmt(z$true_below_lloq_pct)) }),
    bypop(function(pv) { z <- AD_[pk_model == "k2016" & population == pv & study_day == sd_]; sprintf("%s (%s to %s); %s", fmt(z$true_median, 2), fmt(z$true_p05, 2), fmt(z$true_p95, 2), fmt(z$true_below_lloq_pct)) }), "amend1/p1/ada.csv")
  wr <- function(m, pv) { z <- WC[pk_model == m & dose_mg == DOSE & population == pv & metric == "AUCinf" & prob %in% c(0.05, 0.95)][order(prob)]; sprintf("%s / %s (exponent %s)", fmt(z$recomputed[1], 6), fmt(z$recomputed[2], 6), fmt(z$theta_WT[1], 3)) }
  addt("Body weight", "Typical AUC0-inf ratio to 75 kg at the 5th / 95th percentile weight (independent recomputation of the stored values)", bypop(function(pv) wr("k2020", pv)), bypop(function(pv) wr("k2016", pv)), "amend1/weight_check.csv")
  addt("Body weight", "Stored weight ratios reproduced (all doses, populations, metrics, percentiles)", sprintf("%d of %d rows within 1e-6; all within the registered solver tolerance", WC[pk_model == "k2020", sum(pass)], WC[pk_model == "k2020", .N]),
       sprintf("%d of %d rows within 1e-6; all within the registered solver tolerance", WC[pk_model == "k2016", sum(pass)], WC[pk_model == "k2016", .N]), "amend1/weight_check.csv")
  dg <- function(pv) paste(WD[population == pv & role != "registered model (same solver, for comparison)" & prob == 0.05, sprintf("%s %s", variant, fmt(ratio_AUCinf, 3))], collapse = "; ")
  addt("Body weight", "Diagnostic models (not registered): typical AUC0-inf ratio at the 5th percentile weight", bypop(dg), "", "amend1/weight_ratio_diagnostic.csv")
  cdt <- function(m, pv) paste(vapply(c("fixed75", "median"), function(bd) { z <- CD[pk_model == m & population == pv & boundary == bd]
    sprintf("%s kg: stratified %s, unstratified %s; residual CV %s%% vs %s%%; upper stratum %s%%", fmt(z$boundary_kg[1], 1), ci_(z[analysis_model == "M1", power_pct], z[analysis_model == "M1", power_mc_lo], z[analysis_model == "M1", power_mc_hi]),
            ci_(z[analysis_model == "M0", power_pct], z[analysis_model == "M0", power_mc_lo], z[analysis_model == "M0", power_mc_hi]), fmt(z[analysis_model == "M1", cv_auclast_residual_pct]), fmt(z[analysis_model == "M0", cv_auclast_residual_pct]),
            fmt(z$share_upper_stratum_pct[1])) }, ""), collapse = " | ")
  addt("Body weight", "Cell detail at 130 randomized per arm, true GMR 0.95: power (%, Monte Carlo interval), residual CV of AUC0-last, upper-stratum share", bypop(function(pv) cdt("k2020", pv)), bypop(function(pv) cdt("k2016", pv)), "amend1/p1/cell_detail.csv")
  vdt <- function(m, pv) { z <- VD[pk_model == m & population == pv]; sprintf("log weight R2 %s; stratum R2 %s (75 kg), %s (median); absorbed %s%% and %s%%", fmt(z$r2_log_weight[1], 3), fmt(z[boundary == "fixed75", r2_stratum], 3),
                                                                             fmt(z[boundary == "median", r2_stratum], 3), fmt(100 * z[boundary == "fixed75", stratum_share_of_weight_r2], 0), fmt(100 * z[boundary == "median", stratum_share_of_weight_r2], 0)) }
  addt("Body weight", "Variance of log AUC0-last explained by log weight and by the stratum indicator (confirmed schedule)", bypop(function(pv) vdt("k2020", pv)), bypop(function(pv) vdt("k2016", pv)), "amend1/variance_decomposition.csv")
  for (q in unique(PC$quantity)) { z <- function(m) PC[pk_model == m & quantity == q]
    d_ <- if (grepl("AUROC", q)) 3 else 1
    addt("Korean against global population", sprintf("%s: global; Korean; Korean minus global", q), sprintf("%s; %s; %s", ci_(z("k2020")$S1, z("k2020")$S1_lo, z("k2020")$S1_hi, d_), ci_(z("k2020")$S2, z("k2020")$S2_lo, z("k2020")$S2_hi, d_), ci_(z("k2020")$diff_S2_minus_S1, z("k2020")$diff_lo, z("k2020")$diff_hi, d_)),
         sprintf("%s; %s; %s", ci_(z("k2016")$S1, z("k2016")$S1_lo, z("k2016")$S1_hi, d_), ci_(z("k2016")$S2, z("k2016")$S2_lo, z("k2016")$S2_hi, d_), ci_(z("k2016")$diff_S2_minus_S1, z("k2016")$diff_lo, z("k2016")$diff_hi, d_)), "amend1/p1/population_contrast.csv") }
  TAB <- rbindlist(T); m15_fwrite_cond(TAB, file.path(PD1, "p1_results_table.csv"), COND)

  # ---------------- 핵심 결론 5줄(16.3절 4 보수적 쪽 규칙, 16.6절 범주) ----------------
  rb_ <- function(rid) { z <- ROB[row_id == rid]; stopifnot(nrow(z) == 1); z }
  cls_txt <- function(ids) { z <- ROB[row_id %in% ids]; z[, k := substr(class, 1, 1)]
    one <- function(r) { lab <- sub(" \\(.*", "", r$metric)
      if (r$k == "1") sprintf("%s changes in the same direction in both models (smaller difference %s)", lab, if (abs(r$m1_diff) <= abs(r$k16_diff)) r$m1_text else r$k16_text)
      else if (r$k == "2") sprintf("%s shows no material change in either model (larger difference %s)", lab, if (abs(r$m1_diff) >= abs(r$k16_diff)) r$m1_text else r$k16_text)
      else if (r$k == "3") sprintf("%s changes materially in %s only (Model 1 %s; Kovalenko 2016 %s)", lab, if (r$m1_sign != 0) "Model 1" else "Kovalenko 2016", r$m1_text, r$k16_text)
      else sprintf("%s changes in opposite directions (Model 1 %s; Kovalenko 2016 %s)", lab, r$m1_text, r$k16_text) }
    paste(vapply(seq_len(nrow(z)), function(i) one(z[i]), ""), collapse = "; ") }
  ids_s <- function(sn, pv = "S1") c(paste0("sched_power_", sn, "_", pv), paste0("sched_success_pct_", sn, "_", pv), paste0("sched_reliable_pct_", sn, "_", pv), paste0("sched_bias_auclast_pct_", sn, "_", pv),
                                     paste0("sched_rmse_log_auclast_", sn, "_", pv), paste0("sched_auroc_", sn, "_", pv), paste0("sched_Fcross_down_", sn, "_", pv), paste0("sched_Fcross_up_", sn, "_", pv))
  worst <- function(col, m_better = "higher", sn = "B0", pv = "S1") { v <- SS[population == pv & schedule == sn][[col]]; names(v) <- SS[population == pv & schedule == sn, pk_model]
    if (m_better == "higher") v[which.min(v)] else v[which.max(abs(v))] }
  d57 <- CV[schedule == "B0" & planned_elapsed_day == 56]
  C <- c(
    sprintf("1. Extended schedule (Days 64 and 71 added) against the confirmed schedule, global population, 300 mg: %s.", cls_txt(ids_s("ADD64_71"))),
    sprintf("2. Reduced schedule (Days 36 to 50 removed) against the confirmed schedule, global population, 300 mg: %s.", cls_txt(ids_s("RM36_50"))),
    sprintf("3. Coverage on the confirmed schedule (less favourable model quoted): the Day 57 sample is retained as quantifiable in the NCA for %s%% of global and %s%% of Korean subjects; AUC0-inf is estimable with extrapolation below 20%% for %s%% and %s%%, and %s%% and %s%% of subjects have a reliable terminal phase. For anti-drug antibody samples, the median true serum concentration across models and populations is %s mg/L at Day 29, %s mg/L at Day 57 and %s mg/L at Day 85, with up to %s%% of Day 85 values below the LLOQ.",
            fmt(min(d57[population == "S1", used_pct])), fmt(min(d57[population == "S2", used_pct])), fmt(min(SS[population == "S1" & schedule == "B0", lambda_z_success_pct])), fmt(min(SS[population == "S2" & schedule == "B0", lambda_z_success_pct])),
            fmt(min(SS[population == "S1" & schedule == "B0", reliable_terminal_phase_pct])), fmt(min(SS[population == "S2" & schedule == "B0", reliable_terminal_phase_pct])),
            paste(fmt(range(AD_[study_day == 29, true_median]), 1), collapse = " to "), paste(fmt(range(AD_[study_day == 57, true_median]), 2), collapse = " to "), paste(fmt(range(AD_[study_day == 85, true_median]), 2), collapse = " to "),
            fmt(max(AD_[study_day == 85, true_below_lloq_pct]))),
    sprintf("4. Body weight: each model runs with its own weight exponent (typical AUC0-inf ratio at the global 5th percentile weight %s for Model 1 and %s for Kovalenko 2016; near-equal by coincidence of exponent and sensitivity, amendment section 16.4 item 1). Log weight explains %s to %s of the variance of log AUC0-last; the 75 kg stratum term absorbs %s to %s%% of that share in the global and %s to %s%% in the Korean population. Stratified against unstratified analysis: %s; median against 75 kg boundary: %s.",
            fmt(WC[pk_model == "k2020" & dose_mg == DOSE & population == "S1" & metric == "AUCinf" & prob == 0.05, recomputed], 6), fmt(WC[pk_model == "k2016" & dose_mg == DOSE & population == "S1" & metric == "AUCinf" & prob == 0.05, recomputed], 6),
            fmt(min(VD$r2_log_weight), 3), fmt(max(VD$r2_log_weight), 3), fmt(100 * min(VD[population == "S1" & boundary == "fixed75", stratum_share_of_weight_r2]), 0), fmt(100 * max(VD[population == "S1" & boundary == "fixed75", stratum_share_of_weight_r2]), 0),
            fmt(100 * min(VD[population == "S2" & boundary == "fixed75", stratum_share_of_weight_r2]), 0), fmt(100 * max(VD[population == "S2" & boundary == "fixed75", stratum_share_of_weight_r2]), 0),
            cls_txt("strat_term_S1"), cls_txt("strat_boundary_S1")),
    sprintf("5. Korean against global population (300 mg, confirmed schedule): %s.", cls_txt(c("pop_power", "pop_auroc", paste0("pop_cross_", rep(c("F", "ke", "Vmax"), each = 2), "_", c("down", "up"))))))
  stopifnot(length(C) <= 5)
  HDRQ <- "Question to the requester (amendment section 16.2): the request describes Model 1 as the Kovalenko 2020 covariate model (Vc proportional to weight^0.817, ke proportional to BMI^0.368). The registered Model 1 is Kovalenko 2020 Table 1 Model 1 (weight on Vc only, exponent 0.711). Please confirm which model is Model 1 before the appendix is finalized; until then, Model 1 means the registered model, and the covariate and base models are diagnostic only."
  md <- c("# P1 report (amendment 1, plan section 16)", "", HDRQ, "", if (nzchar(COND)) c(COND, "") else NULL,
          "## Key conclusions (P1; conservative-side rule of section 16.3 item 4; classes of section 16.10)", "", C, "",
          "## P1 results table (both models side by side; every value from the file in the last column)", "",
          "| Section | Item | Model 1 (primary) | Kovalenko 2016 (sensitivity) | Source |", "|---|---|---|---|---|", sprintf("| %s | %s | %s | %s | results/m15/%s |", TAB$section, TAB$item, TAB$model_1, TAB$kovalenko_2016, TAB$source), "",
          "## Robustness classification (section 16.10; primary scenario cell and the Korean cell)", "",
          "| Contrast | Metric | Cell | Model 1 | Kovalenko 2016 | Class | Other cells in the same class |", "|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | %s |", ROB$contrast, ROB$metric, ROB$cell, ROB$m1_text, ROB$k16_text, ROB$class, fifelse(is.na(ROB$other_cells_same_class), "-", ROB$other_cells_same_class)))
  pl <- readLines(proj_path("plans", "M15_analysis_plan_20261020.md"), warn = FALSE); i0 <- grep("^### 16\\.13 ", pl); tl <- pl[(i0 + 1):length(pl)]; tl <- tl[seq_len(which(!grepl("^\\||^$", tl))[1] - 1)]
  md <- c(md, "", "## Departures from the request (plan section 16.13)", tl[nzchar(tl)])
  writeLines(md, file.path(PD1, "p1_report.md")); writeLines(c("# P1 key conclusions", "", HDRQ, "", C), file.path(PD1, "p1_key_conclusions.md"))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE); for (ff in c(file.path(PD1, "p1_report.md"), file.path(PD1, "p1_key_conclusions.md"), file.path(PD1, "p1_results_table.csv"), file.path(PD1, "robustness_p1.csv"))) check_english(ff)
  cat(C, sep = "\n"); append_run_log(lf, "P1 report written")
}

# ---------------------------------------------------------------- appendix: 브리핑 부록 초안(16.6절) ----------------------------------------------------------------
if (mode == "appendix") {
  lf <- start_run_log("m15_amend_appendix", master_seed = MS, run_mode = "final")
  APD <- file.path(AD, "appendix"); dir.create(APD, showWarnings = FALSE); PD1 <- file.path(AD, "p1"); SD <- file.path(RD, "summary"); DOSE <- 300; G <- 0.95
  SS <- fread(file.path(PD1, "sampling_schedule.csv")); CR <- fread(file.path(PD1, "crossings.csv")); CV <- fread(file.path(PD1, "coverage.csv")); ADA <- fread(file.path(PD1, "ada.csv"))
  CD <- fread(file.path(PD1, "cell_detail.csv")); ROB <- fread(file.path(PD1, "robustness_p1.csv")); VD <- fread(file.path(AD, "variance_decomposition.csv")); WC <- fread(file.path(AD, "weight_check.csv"))
  AT <- fread(file.path(AD, "ppc", "arm_table.csv")); B7 <- fread(file.path(AD, "ppc", "seven_arm_bias.csv")); RULE <- fread(file.path(RD, "validation", "validation_rule.csv"))
  OCP <- fread(file.path(SD, "oc_pass.csv")); OCR <- fread(file.path(SD, "oc_consumer_risk.csv")); PAN <- fread(file.path(RD, "samplesize", "power_analytic.csv"))
  KC <- readLines(file.path(SD, "m15_key_conclusions_en.md"), warn = FALSE); c1 <- grep("^1\\. Posterior predictive check", KC, value = TRUE); stopifnot(length(c1) == 1)
  gh <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
  tbl <- function(h, rows) c(paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", paste(rep("---", length(h)), collapse = "|"), "|"), vapply(rows, function(r) paste0("| ", paste(r, collapse = " | "), " |"), ""))
  pick <- function(d, m, pv, extra = TRUE) d[pk_model == m & population == pv & dose_mg == DOSE & extra]
  L <- c("# Modelling and simulation appendix: design of the phase 1 pharmacokinetic similarity study", "",
         sprintf("Draft for the European Medicines Agency (EMA) scientific advice briefing document. Structured as in the International Council for Harmonisation (ICH) M15 guideline. Generated on %s from the result files at commit %s (plan `plans/M15_analysis_plan_20261020.md`, amendment 1, section 16). Every number is read from a result file.", format(Sys.Date()), gh), "",
         "Open question: the requester describes Model 1 as the Kovalenko 2020 covariate model. In this appendix Model 1 is the registered Kovalenko 2020 Table 1 Model 1 (body weight on the central volume only). This will be confirmed before the appendix is finalized.", "",
         if (nzchar(COND)) c(paste("Conditionality:", COND), "") else NULL)
  # 1-3 질문, 사용 맥락, 위험
  L <- c(L, "## 1. Question of interest, context of use and model risk", "",
         tbl(c("Element", "Statement"), list(
           c("Question of interest", "Does the planned phase 1 design (single 300 mg subcutaneous dose in healthy adults, parallel groups, weight-stratified randomization, sampling schedule, sample size) detect pharmacokinetic (PK) differences between the proposed biosimilar and the reference product while concluding similarity for similar products?"),
           c("Context of use", "Design support only. The similarity conclusion of the study rests on the non-compartmental analysis (NCA) of the observed data: the 90% confidence interval (CI) of the geometric mean ratio (GMR) within 0.80 to 1.25 for AUC0-last and Cmax. Model results are not used in that conclusion."),
           c("Model influence", "Medium: one of several inputs to the design choices, with statistical sample-size methods, precedent and feasibility."),
           c("Consequence of a wrong decision", "Medium: an underpowered study (a similar product fails) or a design less able to detect a true PK difference; the design can be revised before dosing, and the observed-data NCA remains the basis of the conclusion."),
           c("Model risk", "Medium (medium influence, medium consequence). Credibility activities: code verification, posterior predictive check against the original healthy-volunteer studies, applicability assessment, a structural sensitivity model for every main result, and the uncertainty that is not propagated stated explicitly."))), "")
  # 4 모델
  L <- c(L, "## 2. Models and sources", "",
         tbl(c("Role", "Model", "Structure", "Covariates", "Recorded development data"), list(
           c("Primary", "Kovalenko 2020, Table 1 Model 1", "Two compartments, parallel linear and Michaelis-Menten elimination, three transit compartments and first-order absorption", "Body weight on the central volume (exponent 0.711, reference 75 kg); no other covariate", "R668-AS-0907, TDU12265, PKM14161, R668-AD-1117 (repository record)"),
           c("Sensitivity", "Kovalenko 2016, Table 2 (data below the LLOQ included)", "Two compartments, parallel linear and Michaelis-Menten elimination, first-order absorption", "Body weight on the central volume (exponent 0.705, reference 75 kg)", "R668-AS-0907, R668-HV-1108 and patient studies (Kovalenko 2016 Table 1)"),
           c("Diagnostic only", "Kovalenko 2020 phase 3 covariate and base models (values supplied by the requester)", "As Model 1 with fixed absorption and elimination values; between-subject variability borrowed from Model 1", "Covariate model: weight 0.817 and albumin on the central volume, BMI, EASI and race on linear elimination", "Phase 3 data (repository record); not used for the design grid"))), "",
         "Estimates and structure of the registered models are used as published; nothing is re-estimated. The publications could not be opened from the analysis environment, so the repository records are the source.", "")
  # 5 사후 예측 점검
  pr <- AT[wmode == "primary" & is.finite(observed_mean) & dose_mg >= 300]
  rows5 <- lapply(split(pr, by = c("study", "arm", "endpoint")), function(z) { a <- z[pk_model == "k2020"]; b <- z[pk_model == "k2016"]
    c(sprintf("%s (%s, %d mg, n %d)", a$study, a$arm, a$dose_mg, a$n_subj), sub("_mean$", "", a$endpoint), fmt(a$observed_mean, 1),
      sprintf("%s (%s to %s); %s; %s", fmt(a$pred_median, 1), fmt(a$pi99_lo, 1), fmt(a$pi99_hi, 1), fmt(a$ratio_pred_obs, 2), a$in_sample_status),
      sprintf("%s (%s to %s); %s; %s", fmt(b$pred_median, 1), fmt(b$pi99_lo, 1), fmt(b$pi99_hi, 1), fmt(b$ratio_pred_obs, 2), b$in_sample_status)) })
  b7 <- B7[wmode == "primary" & grepl("^six", arm_set)]
  L <- c(L, "## 3. Posterior predictive check (Li 2020 healthy-volunteer studies)", "",
         "Arm means simulated 5,000 times with each study's design, dose, schedule and LLOQ (0.078 mg/L) under the registered fixed weight assumption (normal 78 kg, SD 10 kg, where weights are not reported). Each cell: predicted median (99% prediction interval); predicted/observed ratio; in-sample status.", "",
         tbl(c("Study arm", "Endpoint", "Observed mean", "Model 1", "Kovalenko 2016"), unname(rows5)), "",
         sprintf("Pooled 300 mg arms with an observed AUC0-last (n 102, observed n-weighted mean %s mg x day/L): predicted/observed ratio %s (95%% prediction interval %s to %s) for Model 1 and %s (%s to %s) for Kovalenko 2016; the simulated mean weight was %s kg. Model 1 over-predicts AUC0-last in these healthy volunteers; Kovalenko 2016 does not.",
                 fmt(b7$observed_mean_6arms[1], 0), fmt(b7[pk_model == "k2020", ratio_median], 3), fmt(b7[pk_model == "k2020", ratio_pi95_lo], 3), fmt(b7[pk_model == "k2020", ratio_pi95_hi], 3),
                 fmt(b7[pk_model == "k2016", ratio_median], 3), fmt(b7[pk_model == "k2016", ratio_pi95_lo], 3), fmt(b7[pk_model == "k2016", ratio_pi95_hi], 3), fmt(b7[pk_model == "k2020", realized_mean_weight_kg], 1)), "",
         sub("^1\\. ", "", c1), "", "The diagnosis of this misfit (alternative healthy-volunteer mappings, false-positive rate of the rule with and without between-study variation) is in progress and is not reported here.", "")
  # 6 설계
  sch <- vapply(names(CFG$schedules), function(sn) sprintf("%s: elapsed days %s", sub(" \\(.*", "", SLAB[[sn]]), paste(CFG$schedules[[sn]], collapse = ", ")), "")
  L <- c(L, "## 4. Simulation design", "",
         tbl(c("Factor", "Setting"), list(
           c("Populations", sprintf("Global: normal %s kg (SD %s), truncated %s to %s kg; Korean: normal %s kg (SD %s), truncated %s to %s kg; inclusion men 60 to 90 kg, women 50 to 90 kg", fmt(as.numeric(CFG$populations$S1$weight$lo) / 2 + as.numeric(CFG$populations$S1$weight$hi) / 2, 2),
                                    fmt((as.numeric(CFG$populations$S1$weight$hi) - as.numeric(CFG$populations$S1$weight$lo)) / 4, 3), CFG$populations$S1$weight$lo, CFG$populations$S1$weight$hi,
                                    fmt(as.numeric(CFG$populations$S2$weight$lo) / 2 + as.numeric(CFG$populations$S2$weight$hi) / 2, 2), fmt((as.numeric(CFG$populations$S2$weight$hi) - as.numeric(CFG$populations$S2$weight$lo)) / 4, 3),
                                    CFG$populations$S2$weight$lo, CFG$populations$S2$weight$hi)),
           c("Dose", "300 mg single subcutaneous dose (200 mg results are not used in this appendix)"),
           c("Sampling schedules", paste(sch, collapse = "; ")),
           c("Sample size", "80 to 160 randomized per arm (base 130); 10% withdrawal completely at random; permuted blocks of two within weight strata"),
           c("Strata and analysis", "Two weight strata at 75 kg or at the enrolled median; stratified (treatment + stratum) and unstratified analysis of log AUC0-last and log Cmax"),
           c("Product differences", "True GMR grid 0.80 to 1.25 (21 points) applied to test-arm concentrations; mechanism versions with test-arm F, ka, ke or Vmax multiplied by 0.60 to 2.00"),
           c("Replicates and seeds", "1,000 simulated studies per cell (10,000 at true GMRs 0.80 and 1.25 for Model 1); common random numbers across designs; 20,000 subjects for the sampling evaluation"),
           c("Software", "R 4.3.3, rxode2 5.1.7, project NCA engine (linear-up log-down; Phoenix best-fit lambda_z rules)"))), "")
  # 7 결과
  cellr <- function(m, pv, gg) { z <- OCP[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE & abs(g - gg) < 1e-9]; stopifnot(nrow(z) == 1); z }
  crr <- function(m, pv, gg) { if (m == PRIMARY) { z <- OCR[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE & abs(g - gg) < 1e-9]; stopifnot(nrow(z) == 1); z } else cellr(m, pv, gg) }
  oc_rows <- lapply(c("S1", "S2"), function(pv) lapply(c("k2020", "k2016"), function(m) { p95 <- cellr(m, pv, 0.95); p1 <- cellr(m, pv, 1.00); c8 <- crr(m, pv, 0.80); c12 <- crr(m, pv, 1.25)
    c(PLAB[[pv]], MLAB[[m]], ci_(p95$pass_pct, p95$pass_lo, p95$pass_hi), ci_(p1$pass_pct, p1$pass_lo, p1$pass_hi), sprintf("%s; %s", ci_(c8$pass_pct, c8$pass_lo, c8$pass_hi), fmt(max(c8$auc_only_pct, c8$cmax_only_pct))),
      sprintf("%s; %s", ci_(c12$pass_pct, c12$pass_lo, c12$pass_hi), fmt(max(c12$auc_only_pct, c12$cmax_only_pct))), format(c8$n_trials, big.mark = ",")) }))
  L <- c(L, "## 5. Results", "", "### 5.1 Operating characteristics (300 mg, confirmed schedule, 75 kg boundary, stratified analysis, 130 randomized per arm)", "",
         tbl(c("Population", "Model", "Power at GMR 0.95 (%)", "Power at GMR 1.00 (%)", "Consumer risk at 0.80 (%): joint; largest single endpoint", "Consumer risk at 1.25 (%): joint; largest single endpoint", "Studies at the margins"), do.call(c, oc_rows)), "",
         "Intervals are Wilson 95% intervals of the simulation (Monte Carlo) error only; they do not include parameter or model uncertainty.", "")
  sr <- lapply(c("S1", "S2"), function(pv) lapply(c("B0", "ADD64_71", "RM36_50"), function(sn) { a <- SS[pk_model == "k2020" & population == pv & schedule == sn]; b <- SS[pk_model == "k2016" & population == pv & schedule == sn]
    c(PLAB[[pv]], sub(" \\(.*", "", SLAB[[sn]]), sprintf("%s / %s", fmt(a$power_pct), fmt(b$power_pct)), sprintf("%s / %s", fmt(a$lambda_z_success_pct), fmt(b$lambda_z_success_pct)),
      sprintf("%s / %s", fmt(a$reliable_terminal_phase_pct), fmt(b$reliable_terminal_phase_pct)), sprintf("%s / %s", fmt(a$auclast_median_error_pct, 1), fmt(b$auclast_median_error_pct, 1)),
      sprintf("%s / %s", fmt(a$auroc_F10_pooled, 3), fmt(b$auroc_F10_pooled, 3))) }))
  d57 <- CV[schedule == "B0" & planned_elapsed_day == 56]
  ada_rows <- lapply(sort(unique(ADA$study_day)), function(sd_) { z <- ADA[study_day == sd_]
    c(sprintf("Day %d", sd_), sprintf("%s to %s", fmt(min(z$true_median), 2), fmt(max(z$true_median), 2)), sprintf("%s to %s", fmt(min(z$true_p05), 2), fmt(max(z$true_p95), 2)), sprintf("%s to %s", fmt(min(z$true_below_lloq_pct)), fmt(max(z$true_below_lloq_pct)))) })
  L <- c(L, "### 5.2 Sampling evaluation (300 mg; Model 1 / Kovalenko 2016)", "",
         tbl(c("Population", "Schedule", "Power at GMR 0.95 (%)", "AUC0-inf estimable, extrapolation below 20% (%)", "Reliable terminal phase (%)", "AUC0-last median error vs dense truth (%)", "AUROC, 10% F difference"), do.call(c, sr)), "",
         sprintf("Day 57 sample on the confirmed schedule: quantifiable and retained in the NCA for %s to %s%% of subjects (both models, both populations). Rows with a confidence interval are compared between models in section 6.", fmt(min(d57$used_pct)), fmt(max(d57$used_pct))), "",
         "Predicted true serum concentrations at candidate anti-drug antibody sampling days (single 300 mg dose; range over both models and both populations; 20,000 subjects each):", "",
         tbl(c("Study day", "Median (mg/L)", "5th to 95th percentile (mg/L)", "Below the LLOQ (%)"), ada_rows), "")
  ssr <- lapply(c("S1", "S2"), function(pv) lapply(c("k2020", "k2016"), function(m) { pa <- PAN[pk_model == m & population == pv & dose_mg == DOSE & schedule == "B0" & boundary == "fixed75" & model == "M1" & n_rand == NBASE & abs(true_gmr - G) < 1e-9]; p95 <- cellr(m, pv, 0.95)
    c(PLAB[[pv]], MLAB[[m]], fmt(pa$cv_model_pct), as.character(pa$n_rand_85), as.character(pa$n_rand_90), ci_(p95$pass_pct, p95$pass_lo, p95$pass_hi), fmt(pa$assurance_cv_pct)) }))
  L <- c(L, "### 5.3 Sample size (true GMR 0.95, stratified analysis, 10% withdrawal)", "",
         tbl(c("Population", "Model", "Model CV of AUC0-last (%)", "Randomized per arm for 85% power", "Randomized per arm for 90% power", "Simulated power at 130 per arm (%)", "Assurance over CV uncertainty at 130 (%)"), do.call(c, ssr)), "")
  wrr <- function(m, pv, pp) WC[pk_model == m & dose_mg == DOSE & population == pv & metric == "AUCinf" & abs(prob - pp) < 1e-9, recomputed]
  L <- c(L, "### 5.4 Body weight and population", "",
         tbl(c("Quantity", "Model 1", "Kovalenko 2016"), list(
           c("Typical AUC0-inf ratio to 75 kg, global 5th / 95th percentile weight", sprintf("%s / %s", fmt(wrr("k2020", "S1", 0.05), 3), fmt(wrr("k2020", "S1", 0.95), 3)), sprintf("%s / %s", fmt(wrr("k2016", "S1", 0.05), 3), fmt(wrr("k2016", "S1", 0.95), 3))),
           c("Typical AUC0-inf ratio to 75 kg, Korean 5th / 95th percentile weight", sprintf("%s / %s", fmt(wrr("k2020", "S2", 0.05), 3), fmt(wrr("k2020", "S2", 0.95), 3)), sprintf("%s / %s", fmt(wrr("k2016", "S2", 0.05), 3), fmt(wrr("k2016", "S2", 0.95), 3))),
           c("Share of log AUC0-last variance explained by log weight (global; Korean)", sprintf("%s; %s", fmt(VD[pk_model == "k2020" & population == "S1", r2_log_weight][1], 3), fmt(VD[pk_model == "k2020" & population == "S2", r2_log_weight][1], 3)), sprintf("%s; %s", fmt(VD[pk_model == "k2016" & population == "S1", r2_log_weight][1], 3), fmt(VD[pk_model == "k2016" & population == "S2", r2_log_weight][1], 3))),
           c("Part of it absorbed by a 75 kg stratum term (global; Korean)", sprintf("%s%%; %s%%", fmt(100 * VD[pk_model == "k2020" & population == "S1" & boundary == "fixed75", stratum_share_of_weight_r2], 0), fmt(100 * VD[pk_model == "k2020" & population == "S2" & boundary == "fixed75", stratum_share_of_weight_r2], 0)),
             sprintf("%s%%; %s%%", fmt(100 * VD[pk_model == "k2016" & population == "S1" & boundary == "fixed75", stratum_share_of_weight_r2], 0), fmt(100 * VD[pk_model == "k2016" & population == "S2" & boundary == "fixed75", stratum_share_of_weight_r2], 0))),
           c("Share of enrolled subjects above 75 kg (global; Korean)", sprintf("%s%%; %s%%", fmt(CD[pk_model == "k2020" & population == "S1" & boundary == "fixed75", share_upper_stratum_pct][1]), fmt(CD[pk_model == "k2020" & population == "S2" & boundary == "fixed75", share_upper_stratum_pct][1])), "same populations"))), "",
         "The two models' weight ratios are near-equal by coincidence: weight acts only on the central volume in both, and the weight exponent times the sensitivity of exposure to bioavailability is 1.110 for Model 1 and 1.107 for Kovalenko 2016. Each model was confirmed to run with its own exponent.", "")
  # 8 비교(robustness) + 본문 결론
  rb <- ROB[!grepl("^dose", contrast) & !grepl("^samplesize", row_id)]
  L <- c(L, "## 6. Comparison between the two models", "",
         "Each design contrast is classified by its sign in each model: a sign counts when the 95% interval excludes 0 and the difference is at least 2 percentage points (0.02 for AUROC and for the root mean square error (RMSE) of the log ratio). Class 1: same direction in both models; class 2: no material effect in either model; class 3: material in one model only; class 4: opposite directions.", "",
         tbl(c("Contrast", "Metric", "Model 1", "Kovalenko 2016", "Class"), lapply(seq_len(nrow(rb)), function(i) c(sub(" \\(.*", "", rb$contrast[i]), sub(" \\(.*", "", rb$metric[i]), rb$m1_text[i], rb$k16_text[i], substr(rb$class[i], 1, 1)))), "")
  body <- rb[substr(class, 1, 1) %in% c("1", "2")]
  bl <- vapply(split(body, by = "contrast"), function(z) sprintf("- %s: %s.", z$contrast[1], paste(vapply(seq_len(nrow(z)), function(i) sprintf("%s %s", sub(" \\(.*", "", z$metric[i]), if (substr(z$class[i], 1, 1) == "1") sprintf("changes in the same direction in both models (Model 1 %s; Kovalenko 2016 %s)", z$m1_text[i], z$k16_text[i]) else "shows no material change in either model"), ""), collapse = "; ")), "")
  ssc <- ROB[grepl("^samplesize", row_id)]
  L <- c(L, "Conclusions supported by both models (classes 1 and 2 only):", "", bl,
         sprintf("- Sample size: simulated power at 130 randomized per arm has a lower 95%% bound of at least 90%%: %s.", paste(sprintf("%s population %s", c(samplesize_S1 = "global", samplesize_S2 = "Korean")[ssc$row_id], ssc$class), collapse = "; ")), "")
  # 9 가정·한계
  L <- c(L, "## 7. Assumptions and limitations", "",
         tbl(c("Item", "Statement"), list(
           c("Model identity", "Model 1 is the registered Kovalenko 2020 Table 1 Model 1; the requester's description matches the phase 3 covariate model (to be confirmed)."),
           c("Posterior predictive check", "Model 1 met the pre-specified misfit rule under the fixed weight assumption (one unit outside its 99% interval) and over-predicts AUC0-last in the healthy volunteers; the cause is under diagnosis. Kovalenko 2016 did not meet the rule. Several studies are in-sample for each model."),
           c("Parameter uncertainty", "Not propagated in these results; published point estimates are used. Propagation from the supplied standard errors is in progress."),
           c("Populations", "The Korean population differs from the global one by body weight only (Japanese TDU12265 weight range in place of national survey data); no race effect is in the registered models."),
           c("Weights of literature studies", "Studies without reported weights assume a normal 78 kg (SD 10 kg) mean weight; results depend on it."),
           c("Assay and windows", "LLOQ 0.078 mg/L and visit windows are placeholders until the assay and the protocol are final."),
           c("Immunogenicity", "Anti-drug antibodies are not modelled; the concentrations in section 5.2 are inputs for setting the assay drug tolerance."),
           c("200 mg presentation", "Excluded from this appendix (known AUC0-last under-prediction for that presentation)."))), "")
  abl <- "Abbreviations: AUC0-inf, area under the concentration-time curve from time zero to infinity; AUC0-last, area to the last quantifiable concentration; AUROC, area under the receiver operating characteristic curve; BMI, body mass index; CI, confidence interval; Cmax, maximum concentration; CV, coefficient of variation; EASI, Eczema Area and Severity Index; F, bioavailability; GMR, geometric mean ratio; ka, absorption rate constant; ke, linear elimination rate constant; lambda_z, terminal elimination rate constant; LLOQ, lower limit of quantification; NCA, non-compartmental analysis; PK, pharmacokinetic; RMSE, root mean square error; SD, standard deviation; Vmax, maximum Michaelis-Menten elimination rate."
  L <- c(L, abl)
  fm <- file.path(APD, "M15_MS_appendix_draft.md"); writeLines(L, fm)
  fd <- file.path(APD, "M15_MS_appendix_draft.docx")
  rc <- system2("pandoc", c("-f", "markdown-smart", "-t", "docx", "-o", shQuote(fd), shQuote(fm)), stdout = TRUE, stderr = TRUE); stopifnot(file.exists(fd))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE); check_english(fm)
  bx <- readBin(fm, "raw", file.info(fm)$size); stopifnot(!any(as.integer(bx) > 127L))
  dt <- officer::docx_summary(officer::read_docx(fd)); tx <- paste(dt$text, collapse = "\n"); ft <- tempfile(fileext = ".txt"); writeLines(tx, ft); check_english(ft)
  stopifnot(!any(utf8ToInt(tx) > 127L))
  nw <- length(strsplit(paste(L, collapse = " "), "\\s+")[[1]])
  append_run_log(lf, sprintf("appendix written: %d words, English and ASCII checks passed (markdown and Word)", nw)); cat("words", nw, "\n")
}

# ---------------------------------------------------------------- ppc: P2 사후 예측 점검 진단(16.7절 4~7) ----------------------------------------------------------------
if (mode == "ppc") {
  lf <- start_run_log("m15_amend_ppc", master_seed = MS, run_mode = "final")
  PD2 <- file.path(AD, "ppc"); dir.create(PD2, showWarnings = FALSE); VA <- file.path(AD, "validation_alt")
  V <- CFG$validation; PPC <- AM$ppc; MAPS <- setNames(PPC$mappings, vapply(PPC$mappings, `[[`, "", "label")); dz <- read_cfg("design_clot2021.yaml")
  REPS <- as.integer(V$replicates); REPS_S <- as.integer(V$replicates_sensitivity); UNITS <- V$rule_units; fr <- V$failure_rule; DET <- as.numeric(unlist(fr$detection$factors)); PI <- as.numeric(V$pi_level)
  q_ <- function(x, p) quantile(x, p, names = FALSE, type = 7)
  RS <- list(registered = rbind(fread(file.path(VA, "validation_reps.csv.gz")), fread(file.path(RD, "validation", "validation_reps.csv.gz"))[wmode %in% c("w70", "w86")][, pk_model := c(k2020 = "a_model1", k2016 = "a_k2016")[pk_model]]),
             request = fread(file.path(VA, "validation_reps_blq0.csv.gz")))
  ST <- fread(file.path(RD, "validation", "validation_rule.csv"))
  stopifnot(all(names(MAPS) %in% unique(RS$registered$pk_model)), all(names(MAPS) %in% unique(RS$request$pk_model)))
  # 단위 통계(합침): 평균 = sum n mean / N, 군 안 CV = sqrt(sum (n - 1) SD^2 / (N - k)) / 평균
  unit_rep <- function(d) d[, { sdA <- AUClast_cv / 100 * AUClast_mean; N <- sum(n); k <- .N; mA <- sum(n * AUClast_mean) / N
    .(AUClast_mean = mA, Cmax_mean = sum(n * Cmax_mean) / N, AUClast_cv_within = 100 * sqrt(sum((n - 1) * sdA^2) / (N - k)) / mA) }, by = rep]
  UN <- c(UNITS, list(PKM14271 = list(arms = c("pkm14271_T", "pkm14271_R"))))                 # 연구별 규칙 두 번째 판(200 mg)
  ustudy <- vapply(names(UN), function(u) V$arms[[unlist(UN[[u]]$arms)[1]]]$study, "")
  ob_arm <- function(aid) { o <- V$arms[[aid]]$observed; if (is.null(o)) return(NULL); g <- function(nm) if (is.null(o[[nm]])) NA_real_ else as.numeric(o[[nm]])
    list(n = as.numeric(V$arms[[aid]]$n_subj), AUClast_mean = g("AUClast_mean"), AUClast_sd = g("AUClast_sd"), Cmax_mean = g("Cmax_mean"), Cmax_sd = g("Cmax_sd")) }
  ob_unit <- function(u) { lst <- lapply(unlist(UN[[u]]$arms), ob_arm); n <- vapply(lst, `[[`, 1, "n"); N <- sum(n)
    mA <- sum(n * vapply(lst, `[[`, 1, "AUClast_mean")) / N; sA <- vapply(lst, `[[`, 1, "AUClast_sd"); mC <- sum(n * vapply(lst, `[[`, 1, "Cmax_mean")) / N
    c(AUClast_mean = mA, Cmax_mean = mC, AUClast_cv_within = 100 * sqrt(sum((n - 1) * sA^2) / (N - length(n))) / mA) }
  OBS <- rbindlist(lapply(names(UN), function(u) { o <- ob_unit(u); data.table(unit = u, metric = names(o), obs = unname(o)) }))
  sig <- sqrt(log(1 + (as.numeric(PPC$between_study_factor$cv_pct) / 100)^2))
  BSF <- rbindlist(lapply(unique(ustudy), function(st) data.table(study = st, rep = seq_len(REPS), g = with_seed(derive_seed(MS, "val", st, "bsf"), exp(sig * rnorm(REPS))))))
  GU <- merge(data.table(unit = names(ustudy), study = ustudy), BSF, by = "study", allow.cartesian = TRUE)[, .(unit, rep, g)]
  ext_units <- function(rm) names(UNITS)[!vapply(names(UNITS), function(u) ustudy[[u]] %in% unlist(dz$model_development_data[[rm]]), TRUE)]
  UR_of <- function(REP, lab, md) { nmax <- REP[pk_model == lab & wmode == md, max(rep)]
    rbindlist(lapply(names(UN), function(u) { arms <- unlist(UN[[u]]$arms); md_ <- if (nrow(REP[pk_model == lab & wmode == md & arm %in% arms])) md else "primary"
      z <- REP[pk_model == lab & wmode == md_ & arm %in% arms & rep <= nmax]; stopifnot(uniqueN(z$arm) == length(arms)); unit_rep(z)[, unit := u] }))[, nmax := nmax][] }
  # 규칙 평가: 단위 집합·수준·임계, 관측 판정, 반복마다 '관측'으로 놓은 오경보율(연구 간 인자 선택), 검출 확률(배율 x 인자)
  eval_rules <- function(UR, rm, use_g, defs) { nmax <- UR$nmax[1]; ur <- merge(UR, GU, by = c("unit", "rep")); farep <- list(); detrep <- lapply(DET, function(f) list()); rows <- list()
    for (rn in names(defs)) { d <- defs[[rn]]; us <- if (identical(d$units, "all")) names(UNITS) else if (identical(d$units, "external")) ext_units(rm) else unlist(d$units)
      thr <- d$outside_pi_min; thr <- as.integer(if (is.list(thr)) thr[[rm]] else thr); lev <- if (is.null(d$pi_level)) PI else as.numeric(d$pi_level); ql <- (1 - lev) / 2
      z <- ur[unit %in% us, c("unit", "rep", "g", d$metric), with = FALSE]; setnames(z, d$metric, "v"); if (!use_g) z[, g := 1]
      z[, `:=`(lo = q_(v, ql), hi = q_(v, 1 - ql)), by = unit]
      im <- grepl("_mean$", d$metric); o <- merge(unique(z[, .(unit, lo, hi)]), OBS[metric == d$metric], by = "unit")[, out := obs < lo | obs > hi]
      zw <- copy(z)[, vg := if (im) g * v else v][, `:=`(lw = q_(vg, ql), hw = q_(vg, 1 - ql)), by = unit]                            # 인자를 곱한 반복의 한계(넓힌 예측구간)
      ow <- merge(unique(zw[, .(unit, lw, hw)]), OBS[metric == d$metric], by = "unit")[, outw := obs < lw | obs > hw]
      kf <- function(f) z[, .(k = sum((if (im) f * g * v else v) < lo | (if (im) f * g * v else v) > hi)), by = rep]
      k1 <- kf(1); farep[[rn]] <- k1[k >= thr, rep]
      for (j in seq_along(DET)) { kk <- kf(DET[j]); detrep[[j]][[rn]] <- kk[k >= thr, rep] }
      rows[[rn]] <- data.table(rule = rn, metric = d$metric, pi_level = lev, units = paste(us, collapse = ";"), threshold = thr, n_outside = sum(o$out), units_outside = paste(o[out == TRUE, unit], collapse = ";"),
                               met = sum(o$out) >= thr, met_widened_limits = sum(ow$outw) >= thr, false_positive_pct = 100 * mean(k1$k >= thr),
                               t(setNames(vapply(seq_along(DET), function(j) 100 * length(detrep[[j]][[rn]]) / nmax, 1), sprintf("detect_x%.2f_pct", DET)))) }
    R <- rbindlist(rows); list(R = R, farep = farep, detrep = detrep, nmax = nmax) }
  union_row <- function(ev, rules, lab) { nmax <- ev$nmax
    data.table(rule = lab, metric = "", pi_level = NA_real_, units = paste(rules, collapse = "+"), threshold = NA_integer_, n_outside = NA_integer_, units_outside = "",
               met = any(ev$R[rule %in% rules, met]), met_widened_limits = any(ev$R[rule %in% rules, met_widened_limits]), false_positive_pct = 100 * length(unique(unlist(ev$farep[rules]))) / nmax,
               t(setNames(vapply(seq_along(DET), function(j) 100 * length(unique(unlist(ev$detrep[[j]][rules]))) / nmax, 1), sprintf("detect_x%.2f_pct", DET)))) }
  DEFS <- list(R1 = fr$R1, R2 = fr$R2, R4 = fr$R4, R5 = fr$R5)
  PS <- PPC$per_study_rule; mk <- function(u, metric) list(metric = metric, units = list(u), outside_pi_min = 1, pi_level = as.numeric(PS$level))
  ps_defs <- function(with14271) { a <- unlist(PS$main$auc); cc <- unlist(PS$main$cmax); if (with14271) { a <- c(a, unlist(PS$with_pkm14271$adds$auc)); cc <- c(cc, unlist(PS$with_pkm14271$adds$cmax)) }
    c(setNames(lapply(a, mk, metric = "AUClast_mean"), paste0("S_auc_", a)), setNames(lapply(cc, mk, metric = "Cmax_mean"), paste0("S_cmax_", cc))) }
  RES <- list()
  for (nr in names(RS)) for (lab in names(MAPS)) for (md in unique(RS[[nr]][pk_model == lab, wmode])) for (ug in c(FALSE, TRUE)) {
    rm <- MAPS[[lab]]$rule_model; UR <- UR_of(RS[[nr]], lab, md)
    ev <- eval_rules(UR, rm, ug, DEFS); dx <- !(MAPS[[lab]]$model %in% c(CFG$models$primary, CFG$models$sensitivity))
    R <- rbind(ev$R, union_row(ev, c("R1", "R2", "R4", "R5"), "any (registered rule)"), union_row(ev, c("R1", "R2", "R5"), "any of R1, R2, R5 (diagnostic-mapping rule)"))
    for (w14 in c(FALSE, TRUE)) { e2 <- eval_rules(UR, rm, ug, ps_defs(w14)); R <- rbind(R, union_row(e2, names(ps_defs(w14)), if (w14) "per-study rule with PKM14271 (9 comparisons)" else "per-study rule (7 comparisons)")) }
    R[, decides := fifelse(dx, rule == "any of R1, R2, R5 (diagnostic-mapping rule)", rule == "any (registered rule)")]
    R[rule == "R4" & dx, rule := "R4 (Model 1 units, descriptive)"]
    RES[[length(RES) + 1]] <- R[, `:=`(mapping = lab, model = MAPS[[lab]]$model, role = MAPS[[lab]]$role, wmode = md, nca_rule = nr, between_study_cv_pct = if (ug) as.numeric(PPC$between_study_factor$cv_pct) else 0, n_rep = ev$nmax)] }
  RT <- rbindlist(RES, fill = TRUE); setcolorder(RT, c("mapping", "model", "role", "wmode", "nca_rule", "between_study_cv_pct", "rule"))
  # 등록 모델, 등록 NCA 규칙, 인자 없음: 저장된 판정·오경보율과 같아야 한다
  ck <- merge(RT[mapping %in% c("a_model1", "a_k2016") & nca_rule == "registered" & between_study_cv_pct == 0 & rule %in% c("R1", "R2", "R4", "R5", "any (registered rule)") & wmode %in% c("primary", "wunc")][,
               `:=`(pk_model = c(a_model1 = "k2020", a_k2016 = "k2016")[mapping], rule2 = sub("any \\(registered rule\\)", "any", rule))],
              ST[, .(pk_model, wmode, rule2 = rule, met_s = met, fa_s = false_alarm_sim_pct)], by = c("pk_model", "wmode", "rule2"))
  stopifnot(nrow(ck) == 20, all(ck$met == ck$met_s), max(abs(ck$false_positive_pct - ck$fa_s)) < 1e-9)
  m15_fwrite_cond(RT, file.path(PD2, "rule_false_positive.csv"), "")
  # 대안 매핑 군 표, 7개 군 편향(같은 정의)
  L_ <- AM$labels$in_sample
  AT <- rbindlist(lapply(split(RS$registered[wmode %in% c("primary", "wunc")], by = c("pk_model", "arm", "wmode")), function(z) { a <- V$arms[[z$arm[1]]]; lab <- z$pk_model[1]; m <- MAPS[[lab]]$model
    rbindlist(lapply(c("Cmax", "AUClast"), function(e) { x <- z[[paste0(e, "_mean")]]; o <- ob_arm(z$arm[1]); o <- if (is.null(o)) NA_real_ else o[[paste0(e, "_mean")]]
      data.table(mapping = lab, model = m, arm = z$arm[1], wmode = z$wmode[1], study = a$study, dose_mg = a$dose_mg, n_subj = a$n_subj, endpoint = e,
                 in_sample_status = if (m %in% names(L_)) L_[[m]][[a$study]] else as.character(L_$diagnostic), observed_mean = o, pred_median = median(x), pi99_lo = q_(x, 0.005), pi99_hi = q_(x, 0.995),
                 ratio_pred_obs = median(x) / o, outside_pi99 = if (is.finite(o)) o < q_(x, 0.005) | o > q_(x, 0.995) else NA) })) }))
  m15_fwrite_cond(AT, file.path(PD2, "arm_table_mappings.csv"), "")
  SA <- PPC$seven_arms; arms <- unlist(SA$arms); oa <- arms[vapply(arms, function(aid) { o <- ob_arm(aid); !is.null(o) && is.finite(o$AUClast_mean) }, TRUE)]
  nn <- vapply(oa, function(aid) as.numeric(V$arms[[aid]]$n_subj), 1); obs7 <- sum(nn * vapply(oa, function(aid) ob_arm(aid)$AUClast_mean, 1)) / sum(nn)
  B7 <- rbindlist(lapply(names(RS), function(nr) rbindlist(lapply(names(MAPS), function(lab) rbindlist(lapply(unique(RS[[nr]][pk_model == lab, wmode]), function(md) {
    z <- RS[[nr]][pk_model == lab & wmode == md & arm %in% oa]; stopifnot(uniqueN(z$arm) == length(oa))
    pr <- z[, .(v = sum(n * AUClast_mean) / sum(n)), by = rep]$v; wt <- z[, .(w = sum(n * WT_mean) / sum(n)), by = rep]$w
    data.table(mapping = lab, model = MAPS[[lab]]$model, wt_ref_kg = if (is.null(MAPS[[lab]]$wt_ref)) 75 else as.numeric(MAPS[[lab]]$wt_ref), wmode = md, nca_rule = nr, observed_mean_6arms = obs7,
               ratio_median = median(pr) / obs7, ratio_pi95_lo = q_(pr, 0.025) / obs7, ratio_pi95_hi = q_(pr, 0.975) / obs7, ratio_pi99_lo = q_(pr, 0.005) / obs7, ratio_pi99_hi = q_(pr, 0.995) / obs7,
               pct_rank_observed = 100 * mean(pr <= obs7), realized_mean_weight_kg = mean(wt), n_rep = length(pr)) }))))))
  m15_fwrite_cond(B7, file.path(PD2, "seven_arm_bias_mappings.csv"), "")
  if (file.exists(file.path(VA, "nca_rule_difference.csv"))) file.copy(file.path(VA, "nca_rule_difference.csv"), file.path(PD2, "nca_rule_difference.csv"), overwrite = TRUE)
  # 결론 1 고정 서식(16.7절 7): 고정 체중 가정, 등록 NCA 규칙
  TH <- as.numeric(PPC$conclusion_template$threshold_fp10_pct)
  pick_ <- function(lab, ug, rl) RT[mapping == lab & wmode == "primary" & nca_rule == "registered" & between_study_cv_pct == ug & rule == rl]
  dev_k <- function(m) sum(vapply(names(UNITS), function(u) ustudy[[u]] %in% unlist(dz$model_development_data[[m]]), TRUE))
  OC <- c(k2020 = "a_model1", k2016 = "a_k2016"); outc <- ST[wmode == "primary" & rule == "any", .(pk_model, outcome)]
  dxmaps <- names(MAPS)[vapply(MAPS, function(x) !(x$model %in% c("k2020", "k2016")) && is.null(x$wt_ref), TRUE)]
  dx_txt <- paste(vapply(dxmaps, function(lab) { r <- pick_(lab, 0, "any of R1, R2, R5 (diagnostic-mapping rule)"); b <- B7[mapping == lab & wmode == "primary" & nca_rule == "registered"]
    sprintf("%s %s the rule (pooled 300 mg AUC0-last ratio %s)", MAPS[[lab]]$role, if (r$met) "met" else "did not meet", fmt(b$ratio_median, 3)) }, ""), collapse = "; ")
  C1 <- vapply(c("k2020", "k2016"), function(m) { lab <- OC[[m]]; r0 <- pick_(lab, 0, "any (registered rule)"); r10 <- pick_(lab, as.numeric(PPC$between_study_factor$cv_pct), "any (registered rule)")
    rr <- RT[mapping == lab & wmode == "primary" & nca_rule == "registered" & between_study_cv_pct == 0 & rule %in% c("R1", "R2", "R4", "R5") & met == TRUE]
    thr_txt <- if (nrow(rr)) paste(sprintf("%s (%s)", rr$rule, gsub(";", " and ", rr$units_outside)), collapse = ", ") else "no rule"
    cat_ <- if (r0$met && r10$false_positive_pct <= TH) "The misfit is unlikely to be explained by chance or by 10% between-study variation." else if (r0$met) "The misfit is compatible with 10% between-study variation and is not taken as evidence of model error." else "The check does not contradict the model."
    sprintf("Posterior predictive check (Li 2020 healthy-volunteer studies; %d of 5 rule units in the recorded development data of %s): %s %s the registered misfit rule (%s) through %s. When %s is correct, the rule is met by chance with probability %s%% without and %s%% with 10%% between-study exposure variation; the probability of detecting a 25%% bias is %s%% and %s%% (the smaller of factors 0.75 and 1.25). %s%s",
            dev_k(m), MLAB[[m]], MLAB[[m]], if (r0$met) "met" else "did not meet", sub("^fails under", "met under", sub("^passes under", "not met under", outc[pk_model == m, outcome])), thr_txt, MLAB[[m]], fmt(r0$false_positive_pct), fmt(r10$false_positive_pct),
            fmt(min(r0$detect_x0.75_pct, r0$detect_x1.25_pct), 0), fmt(min(r10$detect_x0.75_pct, r10$detect_x1.25_pct), 0), cat_,
            if (m == PRIMARY) sprintf(" Mappings (b)+(c), (c) and (d), which are other published models with borrowed variability and not corrections of Model 1: %s.", dx_txt) else "") }, "")
  writeLines(c("# Conclusion 1 (fixed template of plan section 16.7 item 7)", "", paste0("- ", C1)), file.path(PD2, "conclusion1.md"))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE); for (ff in list.files(PD2, "\\.(csv|md)$", full.names = TRUE)) check_english(ff)
  print(RT[decides == TRUE & between_study_cv_pct %in% c(0, 10) & wmode == "primary", .(mapping, nca_rule, between_study_cv_pct, met, false_positive_pct, detect_x0.75_pct, detect_x1.25_pct)]); cat(C1, sep = "\n")
  append_run_log(lf, "done")
}

# ---------------------------------------------------------------- postcheck: 사후 단계 재실행 점검(16.1절) ----------------------------------------------------------------
if (mode == "postcheck") {
  bd <- args[2]; stopifnot(dir.exists(bd)); lf <- start_run_log("m15_amend_postcheck", master_seed = MS, run_mode = "final", extra = list(before = bd))
  ROWS_OK <- c("summary/m15_assumptions_register.csv", "summary/m15_deviations.csv", "summary/m15_summary_table.csv")       # 계획서 표에서 만든 표: 행이 바뀔 수 있다(가정·이탈 추가)
  fl <- list.files(bd, "\\.csv$", recursive = TRUE)
  R <- rbindlist(lapply(fl, function(rel) { a <- fread(file.path(bd, rel), na.strings = c("", "NA")); fb <- file.path(RD, rel)
    if (!file.exists(fb)) return(data.table(file = rel, status = "missing after re-run", numeric_max_rel_diff = NA_real_, text_columns_changed = ""))
    b <- fread(fb, na.strings = c("", "NA"))
    if (!identical(names(a), names(b))) return(data.table(file = rel, status = if (rel %in% ROWS_OK) "columns differ (plan-derived table)" else "columns differ", numeric_max_rel_diff = NA_real_, text_columns_changed = paste(setdiff(union(names(a), names(b)), intersect(names(a), names(b))), collapse = ";")))
    if (nrow(a) != nrow(b)) return(data.table(file = rel, status = if (rel %in% ROWS_OK) "rows differ (plan-derived table)" else "rows differ", numeric_max_rel_diff = NA_real_, text_columns_changed = ""))
    num <- names(a)[vapply(a, is.numeric, TRUE) & vapply(b, is.numeric, TRUE)]
    dn <- if (length(num)) max(vapply(num, function(cc) { x <- a[[cc]]; y <- b[[cc]]; if (!identical(is.na(x), is.na(y))) return(Inf); z <- !is.na(x); if (!any(z)) 0 else max(abs(x[z] - y[z]) / pmax(abs(x[z]), 1e-300)) }, 1)) else 0
    txt <- setdiff(names(a), num); ch <- txt[vapply(txt, function(cc) !identical(as.character(a[[cc]]), as.character(b[[cc]])), TRUE)]
    data.table(file = rel, status = if (dn <= 1e-12 && !length(ch)) "identical" else if (dn <= 1e-12) "text differs only" else "numeric difference", numeric_max_rel_diff = dn, text_columns_changed = paste(ch, collapse = ";")) }))
  ALLOWED_TXT <- c("conditionality", "item", "model_1", "kovalenko_2016", "unit_lab", "summary", "decision_log", "Assumption", "Value used", "Basis", "Impact evidence", "Action before submission", "ID")
  R[, ok := status %in% c("identical", "columns differ (plan-derived table)", "rows differ (plan-derived table)") | (status == "text differs only" & vapply(strsplit(text_columns_changed, ";"), function(x) all(x %in% ALLOWED_TXT), TRUE))]
  fwrite(R, file.path(AD, "post_rerun_check.csv")); print(R[status != "identical"]); append_run_log(lf, sprintf("post re-run check: %d files, %d identical, all ok %s", nrow(R), sum(R$status == "identical"), all(R$ok)))
  if (!all(R$ok)) quit(status = 3)
}
