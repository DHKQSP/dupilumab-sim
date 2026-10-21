#!/usr/bin/env Rscript
# 78_m15_tm_fraction.R — 추록 1(계획서 17절, config addendum_tm_20261021, D-071): 중심구획 target-mediated 소실 비중 f_TM(t).
# 모델 구조·추정값은 바꾸지 않는다(Model 1 주, Kovalenko 2016 민감도). 대상자 1,000명(자기 IIV, 잔차 없음, ADA 없음)은 모든 시나리오·체중·모델·파라미터 변화가 공유한다.
# 모드
#   run [cores]   두 모델 x 시나리오 A(300 mg 단회)·B(600 mg 뒤 300 mg 2주 간격) x 체중 60·75·90 kg x (기준 + Vm·ke·F x 0.8·1.2)
#                 → results/m15/tm/tm_subjects.csv.gz(대상자별 지표), tm_typical.csv(전형 대상자), tm_curves.csv(시각별 f_TM 분위수), tm_checks.csv(17.7절 점검)
#   report        표 1(17.5절 2·3·5), 표 2(17.5절 4), 그림 1(f_TM 곡선과 확정 채혈일), 그림 2(탄력도 막대), 요약 → results/m15/tm/
# 사용법: nice -n 5 Rscript scripts/78_m15_tm_fraction.R run 4; Rscript scripts/78_m15_tm_fraction.R report
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("run", "report")
if (!length(args) || !args[1] %in% MODES) stop("사용법: Rscript scripts/78_m15_tm_fraction.R <", paste(MODES, collapse = "|"), "> ...")
mode <- args[1]
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
TM <- CFG$addendum_tm_20261021; stopifnot(!is.null(TM), identical(TM$decision, "D-071"))
OD <- proj_path(TM$output_dir); dir.create(OD, showWarnings = FALSE, recursive = TRUE)
N <- as.integer(TM$n_subjects); XS <- as.numeric(unlist(TM$thresholds)); WTS <- as.numeric(unlist(TM$weights_kg)); WT_HEAD <- as.numeric(TM$headline_weight_kg)
MODELS <- c(TM$models$primary, TM$models$sensitivity); stopifnot(identical(MODELS, c("k2020", "k2016")))
LLOQ <- study_lloq(); stopifnot(abs(LLOQ - as.numeric(TM$lloq_mg_L)) < 1e-12)
SA <- TM$scenarios$A; SB <- TM$scenarios$B
seq3 <- function(v) { v <- as.numeric(unlist(v)); seq(v[1], v[2], by = v[3]) }
GRID_A <- sort(unique(round(c(seq3(SA$grid$fine), seq3(SA$grid$coarse)), 8)))
T_END <- as.numeric(SA$auc_last_end); T_MAX_A <- max(GRID_A); TINF_A <- as.numeric(SA$aucinf_time); TFB_A <- as.numeric(SA$aucinf_fallback)
SS <- as.numeric(unlist(SB$ss_interval)); PREV <- as.numeric(unlist(SB$previous_interval)); TINF_B <- as.numeric(SB$aucinf_time); TFB_B <- as.numeric(SB$aucinf_fallback)
GRID_B <- sort(unique(round(c(PREV[1], seq(SS[1], SS[2], by = as.numeric(SB$grid_step))), 8)))
DOSE <- list(A = list(t = as.numeric(unlist(SA$dose_times)), a = as.numeric(unlist(SA$dose_mg))), B = list(t = as.numeric(unlist(SB$dose_times)), a = as.numeric(unlist(SB$dose_mg))))
stopifnot(length(DOSE$A$t) == 1L, DOSE$A$a == 300, DOSE$B$a[1] == 600, all(DOSE$B$a[-1] == 300), all(diff(DOSE$B$t) == 14), max(DOSE$B$t) == SS[1], diff(SS) == 14, PREV[2] == SS[1])
PLOT_A <- round(c(0, seq(0.01, 2, by = 0.01), seq(2.05, T_END, by = 0.05)), 8); PLOT_B <- round(seq(SS[1], SS[2], by = 0.05), 8)
stopifnot(all(PLOT_A[-1] %in% GRID_A), all(PLOT_B %in% GRID_B))
# 확정 채혈일(17.5절 1): config trial_design B0 + 투여 전 Day 0, 연구일 d = 경과일 d - 1(Day 1은 0.25일 자리표시자, 가정 B02)
SD_ <- TM$sampling_days; SDAY <- as.integer(unlist(SD_$study_day)); SEL <- as.numeric(unlist(SD_$elapsed))
stopifnot(identical(SEL, c(0, as.numeric(unlist(CFG$schedules$B0)))), length(SDAY) == length(SEL), all(SEL[-(1:2)] == SDAY[-(1:2)] - 1), SDAY[1] == 0L, SDAY[2] == 1L)
PCT <- as.numeric(unlist(TM$percentiles)); QT <- as.integer(TM$quantile_type)
MODL <- c(k2020 = "Model 1 (primary)", k2016 = "Kovalenko 2016 (sensitivity)")
qtl <- function(x, pr) as.numeric(quantile(x, pr, type = QT, names = FALSE, na.rm = FALSE))

# ---------------------------------------------------------------- run ----------------------------------------------------------------
if (mode == "run") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
  SEED <- derive_seed(MS, "tm")
  lf <- start_run_log("m15_tm_run", master_seed = MS, run_mode = "final", extra = list(seed_tm = SEED, n = N, cores = cores))
  Z <- with_seed(SEED, matrix(rnorm(N * length(M15_ETA_NAMES)), N)); colnames(Z) <- paste0("z_", M15_ETA_NAMES)
  EL <- TM$elasticity; MULT <- as.numeric(unlist(EL$multipliers)); PARS <- as.character(unlist(EL$parameters)); stopifnot(identical(PARS, c("Vmax", "ke", "F")), identical(MULT, c(0.8, 1.2)))
  JOBS <- rbind(CJ(model = MODELS, scenario = c("A", "B"), WT = WTS, param = "none", mult = 1, sorted = FALSE),
                CJ(model = MODELS, scenario = c("A", "B"), WT = WTS, param = PARS, mult = MULT, sorted = FALSE))
  append_run_log(lf, sprintf("jobs %d; seed %d; grid A %d points (to %g d), grid B %d points", nrow(JOBS), SEED, length(GRID_A), T_MAX_A, length(GRID_B)))
  ip_of <- function(p, w) {
    subj <- cbind(data.table(id = seq_len(N), WT = w), as.data.table(Z))
    rbind(individual_params(p, data.table(id = 0L, WT = w, ada = 0)), individual_params(p, m15_eta_table(subj, p)))   # id 0 = 전형 대상자(η = 0)
  }
  # AUC0-inf: t_inf에서 몸 안 총량 / 마지막 투여량 > 1e-7이면 fallback 시각까지 다시 푼다(transit 모델은 약 2,700일 뒤 언더플로라 4,000일을 쓰지 않는다)
  auc_inf <- function(s, ipc, p, sc, t_inf, t_fb) {
    d_last <- DOSE[[sc]]$a[length(DOSE[[sc]]$a)]; ti_ <- t_inf
    a <- s[abs(time - ti_) < 1e-6, .(id, auc_inf = auc, resid)]
    long <- a[resid / d_last > 1e-7, id]
    if (length(long)) { s2 <- m15_solve_regimen(ipc[id %in% long], p$model_id, DOSE[[sc]]$t, DOSE[[sc]]$a, t_fb); a[s2[abs(time - t_fb) < 1e-6, .(id, a2 = auc, r2 = resid)], on = "id", `:=`(auc_inf = i.a2, resid = i.r2)] }
    a[, fallback := id %in% long][]
  }
  one <- function(k) {
    k_ <- k; jb <- JOBS[k_]; p <- load_params(jb$model); if (jb$param != "none") p <- apply_multipliers(p, setNames(list(jb$mult), jb$param))
    ip <- ip_of(p, jb$WT); sc <- jb$scenario
    chunks <- split(ip, ceiling(seq_len(nrow(ip)) / 251))
    met <- list(); cur <- list()
    for (ipc in chunks) {
      pr <- ipc[, .(id, Vm = Vmax, ke, Km)]
      if (sc == "A") {
        s <- m15_solve_regimen(ipc, p$model_id, DOSE$A$t, DOSE$A$a, c(GRID_A, TINF_A))
        m <- m15_tm_single(s[time <= T_MAX_A + 1e-9], pr, XS, LLOQ, T_END, T_MAX_A)
        a <- auc_inf(s, ipc, p, "A", TINF_A, TFB_A); m <- merge(m, a[, .(id, AUCinf = auc_inf, resid_inf = resid, fallback)], by = "id")
        if (jb$param == "none") { g <- merge(rbind(data.table(id = ipc$id, time = 0, C = 0), s[round(time, 8) %in% PLOT_A, .(id, time, C)]), pr, by = "id"); cur[[length(cur) + 1L]] <- g[, .(id, time = round(time, 8), C, f = m15_ftm(C, Vm, ke, Km))] }
      } else {
        s <- m15_solve_regimen(ipc, p$model_id, DOSE$B$t, DOSE$B$a, c(GRID_B, TINF_B))
        m <- m15_tm_ss(s, pr, SS, PREV)
        a <- auc_inf(s, ipc, p, "B", TINF_B, TFB_B); m <- merge(m, a[, .(id, AUC_w16_inf_total = auc_inf, resid_inf = resid, fallback)], by = "id")
        m[, `:=`(AUC_w16_inf = AUC_w16_inf_total - auc_at_start)][, c("AUC_w16_inf_total", "auc_at_start") := NULL]
        if (jb$param == "none") { g <- merge(s[round(time, 8) %in% PLOT_B, .(id, time, C)], pr, by = "id"); cur[[length(cur) + 1L]] <- g[, .(id, time = round(time, 8), C, f = m15_ftm(C, Vm, ke, Km))] }
      }
      met[[length(met) + 1L]] <- merge(m, ipc[, .(id, Vm = Vmax, ke_i = ke, Km, Vc_i, ka, F)], by = "id")
    }
    M <- rbindlist(met, fill = TRUE)[, `:=`(model = jb$model, scenario = sc, WT = jb$WT, param = jb$param, mult = jb$mult)]
    CU <- NULL
    if (length(cur)) {
      cc <- rbindlist(cur); ll_ <- LLOQ
      # pq = 참 농도가 LLOQ 이상인 대상자 몫(17.5절 1; 그림 1의 정량 한계 표시와 채혈일별 정량 몫)
      CU <- merge(cc[id > 0, .(q05 = qtl(f, PCT[1]), q50 = qtl(f, PCT[2]), q95 = qtl(f, PCT[3]), pq = mean(C >= ll_), n = .N), by = time], cc[id == 0, .(time, typical = f, typical_C = C)], by = "time")
      CU[, `:=`(model = jb$model, scenario = sc, WT = jb$WT)]
    }
    list(M = M, CU = CU)
  }
  for (m_ in MODELS) invisible(get_model(load_params(m_)$model_id))   # fork 전에 컴파일(자식 프로세스 간 경합 방지, R/mc.R와 같은 규칙)
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(nrow(JOBS)), function(k) tryCatch(one(k), error = function(e) structure(conditionMessage(e), class = "job_error")), mc.cores = cores, mc.preschedule = FALSE)
  bad <- which(vapply(res, function(r) inherits(r, "job_error") || is.null(r), TRUE))
  if (length(bad)) {
    msg <- vapply(res[bad], function(r) if (is.null(r)) "no result (worker died)" else paste(as.character(r), collapse = " "), "")
    append_run_log(lf, paste("job errors:", paste(bad, msg, collapse = " | "))); stop("job errors: ", paste(bad, msg, collapse = " | "))
  }
  append_run_log(lf, sprintf("jobs done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  M <- rbindlist(lapply(res, `[[`, "M"), fill = TRUE); CU <- rbindlist(lapply(res, `[[`, "CU"))
  setcolorder(M, c("model", "scenario", "WT", "param", "mult", "id"))
  TY <- M[id == 0]; SU <- M[id > 0]
  # ---------------- 점검(17.7절) ----------------
  CK <- list()
  # (a) 단회 투여 일정 풀이 = 등록 단회 풀이(.solve_param_table), 같은 격자·허용오차
  p1 <- load_params("k2020"); ip5 <- ip_of(p1, WT_HEAD)[id %in% 1:5]
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", "ktr")
  prm5 <- as.data.frame(ip5[, c(pcols, "F"), with = FALSE]); names(prm5)[names(prm5) == "F"] <- "Fbio"
  r1 <- m15_solve_regimen(ip5, p1$model_id, 0, 300, c(GRID_A, TINF_A)); r2 <- .solve_param_table(get_model(p1$model_id), prm5, c(GRID_A, TINF_A), 300)
  r12 <- merge(r1[, .(id, time, C, auc)], r2[, .(id, time, C2 = C, auc2 = auc)], by = c("id", "time"))
  d_a <- max(abs(r12$C - r12$C2) / pmax(abs(r12$C2), 1e-300), abs(r12$auc - r12$auc2) / pmax(abs(r12$auc2), 1e-300))
  CK[[1]] <- data.table(check = "a_regimen_equals_registered_single_dose_solver", value = d_a, limit = 1e-12, pass = d_a <= 1e-12 && nrow(r12) == 5 * length(c(GRID_A, TINF_A)))
  # (b) 전형 대상자 = m15_typical_metrics(8.6절), 75 kg, 두 모델
  for (m_ in MODELS) {
    ref <- m15_typical_metrics(load_params(m_), WT_HEAD, 300, last_day = T_END); ty <- TY[model == m_ & scenario == "A" & WT == WT_HEAD & param == "none"]
    d_b <- max(abs(ty$AUCinf / ref$AUCinf - 1), abs(ty$Cmax / ref$Cmax - 1), abs(ty$AUClast / ref$AUClast_dense - 1))
    CK[[length(CK) + 1L]] <- data.table(check = paste0("b_typical_equals_m15_typical_metrics_", m_), value = d_b, limit = 1e-5, pass = is.finite(d_b) && d_b <= 1e-5)
  }
  # (c) f_TM(C_x) = x, 모든 대상자
  A0 <- M[scenario == "A"]
  d_c <- max(vapply(XS, function(x) { tag <- sprintf("%02d", round(100 * x)); cx <- A0[[paste0("C_", tag)]]; max(abs(m15_ftm(cx, A0$Vm, A0$ke_i, A0$Km) - x)) }, 1), abs(m15_ftm(A0$C_star, A0$Vm, A0$ke_i, A0$Km) - 0.5))
  CK[[length(CK) + 1L]] <- data.table(check = "c_ftm_at_threshold_concentration_equals_threshold", value = d_c, limit = 1e-12, pass = d_c <= 1e-12)
  # (d) 단봉: 문턱마다 상승·하강 교차 각각 1번 이하(시나리오 A 전체)
  nmx <- max(unlist(lapply(XS, function(x) { tag <- sprintf("%02d", round(100 * x)); c(A0[[paste0("nup_", tag)]], A0[[paste0("ndn_", tag)]]) })))
  CK[[length(CK) + 1L]] <- data.table(check = "d_at_most_one_rising_and_one_falling_crossing", value = nmx, limit = 1, pass = nmx <= 1)
  # (e) 교차 시각에서 다시 푼 농도 = C_x(Model 1 75 kg, 문턱마다 교차한 대상자 앞 20명; registered solve_model)
  for (x in XS) {
    tag <- sprintf("%02d", round(100 * x)); z <- SU[model == "k2020" & scenario == "A" & WT == WT_HEAD & param == "none" & get(paste0("status_", tag)) == "crossed"][order(id)][seq_len(min(20, .N))]
    ipz <- ip_of(p1, WT_HEAD)[match(z$id, id)]
    cz <- vapply(seq_len(nrow(z)), function(i) { ii <- ipz[i]; tx_ <- z[[paste0("t_", tag)]][i]; sol <- solve_model(ii, data.table(id = ii$id, time = tx_), 300, model_id = p1$model_id, atol = TRUTH_ATOL, rtol = TRUTH_RTOL); sol$C[1] }, 1)
    d_e <- if (nrow(z)) max(abs(cz / z[[paste0("C_", tag)]] - 1)) else NA_real_
    CK[[length(CK) + 1L]] <- data.table(check = paste0("e_resolved_concentration_at_crossing_", tag), value = d_e, limit = 1e-4, pass = is.finite(d_e) && d_e <= 1e-4)
  }
  # (f) 유한값, 잔여량
  FIN <- list(A = c("Cmax", "AUClast", "AUCinf", "ftm_cmax", "C_star", paste0("C_", sprintf("%02d", round(100 * XS))), paste0("share_", sprintf("%02d", round(100 * XS))), "amt_tm_share"),
              B = c("AUCtau_ss", "AUCtau_prev", "Cmax_ss", "Ctrough", "Cmin_ss", "Cavg_ss", "AUC_w16_inf", "ftm_trough", "ftm_cavg", "ftm_cmax", "ftm_cmin", "ftm_tmean", "amt_tm_share", "ss_ratio_auc", "ss_ratio_trough"))
  FIN$A <- c(FIN$A, paste0("share_post_", sprintf("%02d", round(100 * XS))), "t_true_last")
  stopifnot(all(unlist(FIN) %in% names(M)))                                    # 이름이 빠지면 점검이 조용히 통과하지 않도록
  nf <- sum(vapply(c("A", "B"), function(sc_) sum(vapply(FIN[[sc_]], function(cc) sum(!is.finite(M[scenario == sc_][[cc]])), 1)), 1))
  CK[[length(CK) + 1L]] <- data.table(check = "f_all_reported_values_finite", value = nf, limit = 0, pass = nf == 0)
  rmax <- max(M$resid_inf / 300); nfb <- sum(M$fallback)                       # 마지막 투여량은 두 시나리오 모두 300 mg
  CK[[length(CK) + 1L]] <- data.table(check = "f_amount_remaining_at_auc_inf_time_over_last_dose", value = rmax, limit = 1e-7, pass = rmax <= 1e-7)
  CK[[length(CK) + 1L]] <- data.table(check = "f_subjects_using_fallback_time (count, descriptive)", value = nfb, limit = NA_real_, pass = TRUE)
  # (g) 다회 투여(시나리오 B): 일정 풀이 = 독립적으로 만든 data.frame 사건표의 rxSolve(Model 1 75 kg 앞 5명; 98, 112, 119, 126일 C와 auc, 상대 1e-10),
  #     그리고 두 번째 이후 투여마다 직전·직후(± 1e-6일) 몸 안 총량 증가 = 투여량 x F(상대 1e-4)
  tB <- c(98, 112, 119, 126); nd <- length(DOSE$B$t)
  sB <- m15_solve_regimen(ip5, p1$model_id, DOSE$B$t, DOSE$B$a, tB)
  evd <- data.frame(time = c(DOSE$B$t, tB), amt = c(DOSE$B$a, rep(NA_real_, length(tB))), evid = c(rep(1L, nd), rep(0L, length(tB))), cmt = c(rep("depot", nd), rep(NA_character_, length(tB))))
  evd <- evd[order(evd$time, -evd$evid), ]
  sI <- rxode2::rxSolve(get_model(p1$model_id), params = prm5, events = evd, atol = TRUTH_ATOL, rtol = TRUTH_RTOL, maxsteps = 2000000L, cores = 1L, returnType = "data.table", addDosing = FALSE)
  if ("sim.id" %in% names(sI)) sI[, id := prm5$id[sim.id]] else sI[, id := prm5$id[as.integer(id)]]
  gI <- merge(sB[, .(id, time, C, auc)], sI[, .(id, time, C2 = C, auc2 = auc)], by = c("id", "time"))
  d_g1 <- if (nrow(gI) == 5 * length(tB)) max(abs(gI$C / gI$C2 - 1), abs(gI$auc / gI$auc2 - 1)) else Inf
  CK[[length(CK) + 1L]] <- data.table(check = "g_regimen_equals_independent_event_table", value = d_g1, limit = 1e-10, pass = d_g1 <= 1e-10)
  dj <- DOSE$B$t[-1]; sJ <- m15_solve_regimen(ip5, p1$model_id, DOSE$B$t, DOSE$B$a, sort(c(dj - 1e-6, dj + 1e-6)))
  sJ <- merge(sJ[, .(id, time, resid)], ip5[, .(id, F)], by = "id"); setorder(sJ, id, time)
  jmp <- sJ[, .(jump = resid[seq(2, .N, 2)] - resid[seq(1, .N, 2)], amt = DOSE$B$a[-1] * F[1]), by = id]
  d_g2 <- max(abs(jmp$jump / jmp$amt - 1))
  CK[[length(CK) + 1L]] <- data.table(check = "g_amount_jump_at_each_repeat_dose_equals_dose_times_F", value = d_g2, limit = 1e-4, pass = nrow(jmp) == 5 * (nd - 1) && d_g2 <= 1e-4)
  CK <- rbindlist(CK)
  fwrite(SU, file.path(OD, "tm_subjects.csv.gz")); fwrite(TY, file.path(OD, "tm_typical.csv")); fwrite(CU, file.path(OD, "tm_curves.csv")); fwrite(CK, file.path(OD, "tm_checks.csv"))
  print(CK)
  append_run_log(lf, sprintf("checks: %d of %d pass; subjects rows %d; fallback %d", sum(CK$pass), nrow(CK), nrow(SU), nfb))
  if (!all(CK$pass)) quit(status = 3)
}

# ---------------------------------------------------------------- report ----------------------------------------------------------------
if (mode == "report") {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  lf <- start_run_log("m15_tm_report", master_seed = MS, run_mode = "final")
  CK <- fread(file.path(OD, "tm_checks.csv")); stopifnot(all(CK$pass))
  SU <- fread(file.path(OD, "tm_subjects.csv.gz")); TY <- fread(file.path(OD, "tm_typical.csv")); CU <- fread(file.path(OD, "tm_curves.csv"))
  DEC <- as.integer(TM$display_decimals); stopifnot(DEC == 1L)
  f1 <- function(x) { s <- formatC(round(x, DEC), format = "f", digits = DEC); s[s == "-0.0"] <- "0.0"; s[is.infinite(x) & x > 0] <- "not reached"; s[is.na(x)] <- "NA"; s }
  day <- function(t) t + 1                                     # 연구일 = 투여 후 경과일 + 1(config day_convention)
  # Monte Carlo 오차(17.6절): 대상자 부트스트랩, 모든 행이 같은 재표집 색인(derive_seed(20261020, "tm_boot"))
  NB <- as.integer(TM$mc_error$n_boot); BIDX <- with_seed(derive_seed(MS, "tm_boot"), matrix(sample.int(N, N * NB, replace = TRUE), N))
  B0 <- SU[param == "none"]; T0 <- TY[param == "none"]
  rows <- list()
  # x: 대상자 순서(id 1..N) 벡터. cond = TRUE이면 NA는 그 행의 대상에서 빠진 대상자(조건부 요약), 아니면 NA가 없어야 한다
  add <- function(section, quantity, unit, model, w, typ, x, note = "", cond = FALSE) {
    stopifnot(length(x) == N); if (!cond) stopifnot(!anyNA(x))
    v <- x[!is.na(x)]; q <- if (length(v)) qtl(v, PCT) else rep(NA_real_, 3)
    se <- if (length(v) >= 2) { X <- matrix(x[BIDX], N); bq <- apply(X, 2, function(z) { z <- z[!is.na(z)]; if (length(z)) qtl(z, PCT) else rep(NA_real_, 3) }); apply(bq, 1, sd) } else rep(NA_real_, 3)
    rows[[length(rows) + 1L]] <<- data.table(section = section, quantity = quantity, unit = unit, model = model, WT = w, typical = typ, median = q[2], p05 = q[1], p95 = q[3],
      n_subjects = length(v), se_median = se[2], se_p05 = se[1], se_p95 = se[3], share_pct = NA_real_, note = note)
  }
  addshare <- function(section, quantity, model, w, flag, note = "") rows[[length(rows) + 1L]] <<- data.table(section = section, quantity = quantity, unit = "% of subjects", model = model, WT = w,
    typical = NA_real_, median = NA_real_, p05 = NA_real_, p95 = NA_real_, n_subjects = N, se_median = NA_real_, se_p05 = NA_real_, se_p95 = NA_real_, share_pct = 100 * mean(flag),
    note = sprintf("%s%sMonte Carlo SE %s percentage points", note, if (nzchar(note)) "; " else "", formatC(100 * sqrt(mean(flag) * (1 - mean(flag)) / N), format = "f", digits = 2)))
  SX <- "Threshold concentrations (do not depend on weight)"; SA_ <- sprintf("A: %g mg single SC dose", DOSE$A$a); SB_ <- sprintf("B: 600 mg, then 300 mg every 2 weeks; steady-state interval (days %g to %g)", SS[1], SS[2])
  for (m_ in MODELS) for (w in WTS) {
    a <- B0[model == m_ & scenario == "A" & WT == w][order(id)]; ta <- T0[model == m_ & scenario == "A" & WT == w]
    b <- B0[model == m_ & scenario == "B" & WT == w][order(id)]; tb <- T0[model == m_ & scenario == "B" & WT == w]
    stopifnot(nrow(a) == N, nrow(b) == N, nrow(ta) == 1, nrow(tb) == 1, identical(a$id, seq_len(N)), identical(b$id, seq_len(N)))
    add(SX, "C* = Vm/ke - Km (f_TM = 50%)", "mg/L", m_, w, ta$C_star, a$C_star)
    for (x in XS[-1]) { tag <- sprintf("%02d", round(100 * x)); add(SX, sprintf("Concentration at f_TM = %d%%", round(100 * x)), "mg/L", m_, w, ta[[paste0("C_", tag)]], a[[paste0("C_", tag)]]) }
    add(SA_, "Cmax", "mg/L", m_, w, ta$Cmax, a$Cmax)
    add(SA_, "Minimum f_TM (at Cmax)", "%", m_, w, 100 * ta$ftm_cmax, 100 * a$ftm_cmax)
    add(SA_, "True Tlast at the LLOQ", "study Day", m_, w, day(ta$t_true_last), day(a$t_true_last), note = sprintf("concentrations after it are below the LLOQ (%s mg/L)", format(LLOQ)))
    for (x in XS) {
      tag <- sprintf("%02d", round(100 * x)); st <- a[[paste0("status_", tag)]]; tx <- a[[paste0("t_", tag)]]; tx[st %in% c("censored", "never")] <- Inf   # 도달하지 않음(분위수에서 +Inf로)
      sty <- ta[[paste0("status_", tag)]]; txy <- ta[[paste0("t_", tag)]]; if (sty %in% c("censored", "never")) txy <- Inf
      add(SA_, sprintf("Study Day from which f_TM stays at or above %d%%, all subjects", round(100 * x)), "study Day", m_, w, day(txy), day(tx),
          note = "study Day = days after the dose + 1; Day 1.0 = the dose time (f_TM at or above the threshold from the dose onward), not the Day 1 sample")
      xc <- ifelse(st == "throughout", NA_real_, tx)
      add(SA_, sprintf("Study Day from which f_TM stays at or above %d%%, subjects whose f_TM falls below it", round(100 * x)), "study Day", m_, w, if (sty == "throughout") NA_real_ else day(txy), day(xc), cond = TRUE)
      cc_ <- ifelse(st == "crossed", a[[paste0("C_", tag)]], NA_real_)
      add(SA_, sprintf("Concentration at that time, subjects whose f_TM falls below %d%% (= concentration at f_TM = %d%%)", round(100 * x), round(100 * x)), "mg/L", m_, w,
          if (sty == "crossed") ta[[paste0("C_", tag)]] else NA_real_, cc_, note = "equals the threshold concentration by construction", cond = TRUE)
      addshare(SA_, sprintf("Subjects with f_TM at or above %d%% from the dose onward (Cmax not above the threshold)", round(100 * x)), m_, w, st == "throughout")
      addshare(SA_, sprintf("Subjects not reaching f_TM %d%% by day %g", round(100 * x), T_MAX_A), m_, w, st %in% c("censored", "never"))
    }
    add(SA_, "Share of AUC0-last accumulated while f_TM is at or above 50%", "%", m_, w, 100 * ta$share_50, 100 * a$share_50)
    add(SA_, "Of which after Cmax (falling phase; supplementary)", "%", m_, w, 100 * ta$share_post_50, 100 * a$share_post_50)
    for (x in XS[-1]) { tag <- sprintf("%02d", round(100 * x))
      add(SA_, sprintf("Share of AUC0-last accumulated while f_TM is at or above %d%% (supplementary)", round(100 * x)), "%", m_, w, 100 * ta[[paste0("share_", tag)]], 100 * a[[paste0("share_", tag)]]) }
    add(SB_, "f_TM at Ctrough,ss", "%", m_, w, 100 * tb$ftm_trough, 100 * b$ftm_trough)
    add(SB_, "f_TM at Cavg,ss", "%", m_, w, 100 * tb$ftm_cavg, 100 * b$ftm_cavg)
    add(SB_, "f_TM at Cmax,ss (lowest f_TM in the interval)", "%", m_, w, 100 * tb$ftm_cmax, 100 * b$ftm_cmax)
    add(SB_, "f_TM at Cmin,ss (highest f_TM in the interval)", "%", m_, w, 100 * tb$ftm_cmin, 100 * b$ftm_cmin)
    add(SB_, "Time-averaged f_TM over the interval (supplementary)", "%", m_, w, 100 * tb$ftm_tmean, 100 * b$ftm_tmean)
    add(SB_, "Ctrough,ss", "mg/L", m_, w, tb$Ctrough, b$Ctrough)
    add(SB_, "Cavg,ss", "mg/L", m_, w, tb$Cavg_ss, b$Cavg_ss)
    add(SB_, "Cmax,ss", "mg/L", m_, w, tb$Cmax_ss, b$Cmax_ss)
    add(SB_, "Cmin,ss", "mg/L", m_, w, tb$Cmin_ss, b$Cmin_ss)
    add(SB_, sprintf("Steady-state check: AUCtau (days %g to %g) / AUCtau (days %g to %g)", SS[1], SS[2], PREV[1], PREV[2]), "ratio x 100", m_, w, 100 * tb$ss_ratio_auc, 100 * b$ss_ratio_auc)
    add(SB_, sprintf("Steady-state check: C(%g) / C(%g)", SS[2], SS[1]), "ratio x 100", m_, w, 100 * tb$ss_ratio_trough, 100 * b$ss_ratio_trough)
    add(SB_, sprintf("AUCtau,ss / single-dose AUC0-inf of scenario A, same subject (supplementary; 100 under linear PK)"), "ratio x 100", m_, w, 100 * tb$AUCtau_ss / ta$AUCinf, 100 * b$AUCtau_ss / a$AUCinf)
  }
  T1 <- rbindlist(rows)[, model_label := MODL[model]]
  fwrite(T1, file.path(OD, "tm_table1.csv"))
  T1[, txt := ""][is.na(share_pct), txt := sprintf("%s; %s (%s to %s)%s", f1(typical), f1(median), f1(p05), f1(p95), fifelse(n_subjects < N, sprintf(" [n = %d]", n_subjects), ""))][!is.na(share_pct), txt := f1(share_pct)]
  wd <- dcast(T1, section + quantity + unit + model ~ WT, value.var = "txt", fun.aggregate = function(z) z[1], fill = "")
  ord <- unique(T1[, .(section, quantity, model)])[, o := .I]; wd <- merge(wd, ord, by = c("section", "quantity", "model"))[order(match(model, MODELS), o)]
  wc <- as.character(WTS)
  md1 <- c(sprintf("# Table 1. Target-mediated share of central elimination after a single %g mg SC dose (A) and at steady state (B)", DOSE$A$a), "",
           sprintf("Plan section 17 (addendum 1, D-071). f_TM = [Vm/(Km + C)] / [Vm/(Km + C) + ke], C = central concentration. %s subjects per scenario with between-subject variability (the same subjects in every scenario, weight and model), no residual error, no anti-drug antibodies. Each cell: typical subject; median (5th to 95th percentile) over subjects. Rows in %% of subjects give shares of the %s subjects. Model 1 = Kovalenko 2020 Table 1 Model 1 (primary model); Kovalenko 2016 = sensitivity model. Study Day = days after the (first) dose + 1.",
                   format(N, big.mark = ","), format(N, big.mark = ",")), "")
  for (m_ in MODELS) {
    md1 <- c(md1, sprintf("## %s", MODL[[m_]]), "", sprintf("| Section | Quantity | Unit | %s |", paste(paste(wc, "kg"), collapse = " | ")), paste0("|", paste(rep("---", 3 + length(wc)), collapse = "|"), "|"))
    z <- wd[model == m_]
    for (i in seq_len(nrow(z))) md1 <- c(md1, sprintf("| %s | %s | %s | %s |", z$section[i], z$quantity[i], z$unit[i], paste(vapply(wc, function(cc) z[[cc]][i], ""), collapse = " | ")))
    md1 <- c(md1, "")
  }
  md1 <- c(md1, "Notes.",
           "- Study Day 1.0 in a crossing row is the dose time: f_TM is at or above the threshold from the dose onward because Cmax does not exceed the threshold concentration. The rows for subjects whose f_TM falls below the threshold leave those subjects out. 'not reached' = f_TM still below the threshold at day 140, or the percentile includes such subjects.",
           "- In the rows for subjects whose f_TM falls below the threshold, [n = ...] is the number of such subjects, and NA in the typical-subject position means that the typical subject does not fall below the threshold.",
           "- The concentration at the crossing equals the threshold concentration by construction; it depends only on Vm, ke and Km, not on weight. Below the LLOQ, f_TM approaches Vm / (Vm + ke Km), which depends on the fixed Km.",
           sprintf("- AUC0-last = model AUC from the dose to the earlier of the true Tlast at the LLOQ (%s mg/L) and day %g after the dose (study Day %g), the dense-sampling truth of plan section 8.2(a). The share counts every time with f_TM at or above the threshold: the absorption phase before C first exceeds the threshold concentration and the phase after the crossing ('of which after Cmax').", format(LLOQ), T_END, T_END + 1),
           sprintf("- Numbers are shown to one decimal as requested; Monte Carlo standard errors of the percentiles (subject bootstrap, %s resamples) are in tm_table1.csv.", format(NB, big.mark = ",")))
  writeLines(md1, file.path(OD, "tm_table1.md"))
  # 채혈일별 참 농도 >= LLOQ 몫(보충, 17.5절 1): 그림 1의 정량 한계 표시와 같은 자료
  QD <- CU[scenario == "A" & round(time, 8) %in% round(SEL, 8)][, `:=`(elapsed = time, study_day = SDAY[match(round(time, 8), round(SEL, 8))], quantifiable_pct = 100 * pq)]
  QD <- QD[, .(model, model_label = MODL[model], WT, study_day, elapsed, quantifiable_pct, n)][order(match(model, MODELS), WT, elapsed)]
  stopifnot(nrow(QD) == length(MODELS) * length(WTS) * length(SEL))
  fwrite(QD, file.path(OD, "tm_quantifiable_by_day.csv"))
  # ---------------- 표 2(17.5절 4): 같은 대상자의 기하평균 % 변화 ----------------
  MA <- c("AUClast", "Cmax", "AUCinf"); MB <- c("AUCtau_ss", "Cmax_ss", "Ctrough", "AUC_w16_inf")
  MLAB <- c(AUClast = sprintf("A: AUC0-last (sampling to Day %d)", T_END + 1), Cmax = "A: Cmax", AUCinf = "A: AUC0-inf (reference)", AUCtau_ss = "B: AUCtau,ss", Cmax_ss = "B: Cmax,ss",
            Ctrough = TM$elasticity$labels$Ctrough, AUC_w16_inf = TM$elasticity$labels$AUC_w16_inf)
  base <- SU[param == "none"]; pert <- SU[param != "none"]
  el <- rbindlist(lapply(c(MA, MB), function(mt) {
    sc_ <- if (mt %in% MA) "A" else "B"
    j <- merge(pert[scenario == sc_, c("model", "scenario", "WT", "param", "mult", "id", mt), with = FALSE], base[scenario == sc_, c("model", "WT", "id", mt), with = FALSE], by = c("model", "WT", "id"), suffixes = c("", "_base"))
    lr <- log(j[[mt]] / j[[paste0(mt, "_base")]]); j[, lr := lr]
    tt <- merge(TY[param != "none" & scenario == sc_, c("model", "WT", "param", "mult", mt), with = FALSE], TY[param == "none" & scenario == sc_, c("model", "WT", mt), with = FALSE], by = c("model", "WT"), suffixes = c("", "_base"))
    tt[, typ := 100 * (get(mt) / get(paste0(mt, "_base")) - 1)]
    o <- j[, .(gm_pct = 100 * (exp(mean(lr)) - 1), se_gm_pct = 100 * exp(mean(lr)) * sd(lr) / sqrt(.N), median_pct = 100 * (exp(median(lr)) - 1), p05_pct = 100 * (exp(qtl(lr, PCT[1])) - 1), p95_pct = 100 * (exp(qtl(lr, PCT[3])) - 1), n = .N),
           by = .(model, scenario, WT, param, mult)]
    o <- merge(o, tt[, .(model, WT, param, mult, typical_pct = typ)], by = c("model", "WT", "param", "mult"))
    o[, `:=`(metric = mt, metric_label = MLAB[[mt]])]
  }))
  stopifnot(all(el$n == N))
  el[, change := sprintf("%s %s%d%%", fifelse(param == "Vmax", "Vm", param), fifelse(mult > 1, "+", "-"), round(abs(mult - 1) * 100))]
  fwrite(el[, .(model, model_label = MODL[model], scenario, WT, param, mult, change, metric, metric_label, gm_pct, se_gm_pct, median_pct, p05_pct, p95_pct, typical_pct, n)], file.path(OD, "tm_elasticity_full.csv"))
  CH <- c("Vm -20%", "Vm +20%", "ke -20%", "ke +20%", "F -20%", "F +20%")
  T2 <- dcast(el[WT == WT_HEAD], model + change ~ metric, value.var = "gm_pct")[order(match(model, MODELS), match(change, CH))]
  setcolorder(T2, c("model", "change", MA, MB)); T2[, model_label := MODL[model]]
  fwrite(T2, file.path(OD, "tm_table2.csv"))
  md2 <- c(sprintf("# Table 2. Percentage change in exposure when Vm, ke or F is changed by 20%% (%g kg)", WT_HEAD), "",
           sprintf("Plan section 17.5 item 4 (addendum 1, D-071). Change in the geometric mean over the same %s subjects (between-subject variability included, no residual error) when the typical value of one parameter is multiplied by 0.8 or 1.2 in both scenarios. A = %g mg single SC dose; B = 600 mg, then 300 mg every 2 weeks, steady-state interval days %g to %g after the first dose (study Days %g to %g).",
                   format(N, big.mark = ","), DOSE$A$a, SS[1], SS[2], SS[1] + 1, SS[2] + 1), "",
           sprintf("| Model | Change | %s |", paste(MLAB[c(MA, MB)], collapse = " | ")), paste0("|", paste(rep("---", 2 + length(c(MA, MB))), collapse = "|"), "|"))
  for (i in seq_len(nrow(T2))) md2 <- c(md2, sprintf("| %s | %s | %s |", T2$model_label[i], T2$change[i], paste(vapply(c(MA, MB), function(cc) f1(T2[[cc]][i]), ""), collapse = " | ")))
  md2 <- c(md2, "", "Notes.",
           sprintf("- AUC0-last = model AUC from the dose to the earlier of the true Tlast at the LLOQ and study Day %g. In scenario B, AUC0-last and Cmax of the request are read as AUCtau,ss and Cmax,ss; Ctrough,ss is added.", T_END + 1),
           "- The reference columns are not study endpoints. The B reference is the AUC from the Week 16 dose to infinity with dosing stopped; it includes the washout of drug from all earlier doses, so it is not the steady-state counterpart of the single-dose AUC0-inf (that counterpart is AUCtau,ss).",
           "- Other weights, the Monte Carlo standard error, the median individual change with its 5th and 95th percentiles, and the typical subject: tm_elasticity_full.csv.")
  writeLines(md2, file.path(OD, "tm_table2.md"))
  # ---------------- 그림 1(17.5절 1) ----------------
  PAN <- c(sprintf("A: %g kg", WTS), sprintf("B: %g kg, steady state", WT_HEAD))
  cu <- CU[(scenario == "A") | (scenario == "B" & WT == WT_HEAD)][, panel := factor(fifelse(scenario == "A", sprintf("A: %g kg", WT), sprintf("B: %g kg, steady state", WT)), levels = PAN)]
  cu[, model_lab := factor(MODL[model], levels = MODL[MODELS])]
  sl <- CJ(model_lab = factor(MODL[MODELS], levels = MODL[MODELS]), panel = factor(PAN[seq_along(WTS)], levels = PAN), k = seq_along(SEL))
  sl[, `:=`(time = SEL[k], lab = as.character(SDAY[k]), y = c(104, 110, 116)[(k - 1) %% 3 + 1])]
  tq <- B0[scenario == "A", .(lo = qtl(t_true_last, PCT[1]), md = qtl(t_true_last, PCT[2]), hi = qtl(t_true_last, PCT[3])), by = .(model, WT)]
  tq[, `:=`(model_lab = factor(MODL[model], levels = MODL[MODELS]), panel = factor(sprintf("A: %g kg", WT), levels = PAN), hi_draw = pmin(hi, T_END), y = 6)]
  rl <- data.table(y = 100 * XS)
  g1 <- ggplot(cu, aes(time)) +
    geom_hline(data = rl, aes(yintercept = y), colour = VIZ$muted, linewidth = 0.3, linetype = "dashed") +
    geom_vline(data = sl, aes(xintercept = time), colour = VIZ$muted, linewidth = 0.25, alpha = 0.7) +
    geom_text(data = sl, aes(x = time, y = y, label = lab), size = 2.2, colour = VIZ$ink2) +
    geom_ribbon(aes(ymin = 100 * q05, ymax = 100 * q95), fill = VIZ$s1, alpha = 0.22) +
    geom_line(aes(y = 100 * q50), colour = VIZ$s1, linewidth = 0.7) +
    geom_segment(data = tq, aes(x = lo, xend = hi_draw, y = y, yend = y), colour = VIZ$ink2, linewidth = 1.1, inherit.aes = FALSE) +
    geom_point(data = tq, aes(x = md, y = y), colour = VIZ$ink, size = 1.8, inherit.aes = FALSE) +
    facet_grid(model_lab ~ panel, scales = "free_x") +
    scale_y_continuous(breaks = c(0, 25, 50, 80, 95), limits = c(0, 119), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(breaks = function(l) if (l[2] - l[1] > 30) seq(0, 56, by = 7) else seq(SS[1], SS[2], by = 7), expand = expansion(mult = c(0.01, 0.02)),
                       sec.axis = sec_axis(~ . + 1, name = "Study Day", breaks = function(l) if (l[2] - l[1] > 30) seq(1, T_END + 1, by = 7) else seq(SS[1] + 1, SS[2] + 1, by = 7))) +
    labs(x = "Time after the first dose (days)", y = "f_TM (%)",
         title = "Target-mediated share of central elimination, f_TM",
         subtitle = sprintf("Line = median, band = 5th to 95th percentile of %s subjects; dashed lines = 50%%, 80%%, 95%%. Grey lines and numbers = confirmed sampling days (study Day; Day 0 = pre-dose, Day 1 = %g hours\nafter the dose (placeholder, assumption B02), Day d = d - 1 days after the dose). Black bar near the bottom = 5th to 95th percentile of the true Tlast at the LLOQ (%s mg/L; point = median):\nlater f_TM values are at concentrations below the LLOQ. A = %g mg single SC dose; B = 600 mg, then 300 mg every 2 weeks, days %g to %g.",
                            format(N, big.mark = ","), 24 * SEL[2], format(LLOQ), DOSE$A$a, SS[1], SS[2])) +
    theme_dupi() + theme(panel.spacing.x = unit(1.4, "lines"), plot.subtitle = element_text(size = 8))
  ggsave(file.path(OD, "tm_fig1_ftm_curves.png"), g1, width = 12, height = 6.8, dpi = 150, bg = "white")
  # ---------------- 그림 2(17.5절 4) ----------------
  PAIR <- c(AUClast = "AUC", AUCtau_ss = "AUC", Cmax = "Cmax", Cmax_ss = "Cmax", AUCinf = "Reference", AUC_w16_inf = "Reference")
  PL <- c(AUC = sprintf("AUC (A: AUC0-last, sampling to Day %d; B: AUCtau,ss)", T_END + 1), Cmax = "Cmax (A: Cmax; B: Cmax,ss)",
          Reference = "Reference (A: AUC0-inf; B: AUC from the Week 16 dose to infinity, dosing stopped)")
  e2 <- el[WT == WT_HEAD & metric %in% names(PAIR)][, `:=`(pair = factor(PL[PAIR[metric]], levels = PL), scen = factor(fifelse(scenario == "A", "A: 300 mg single dose", "B: steady state"), levels = c("A: 300 mg single dose", "B: steady state")),
                                                   change = factor(change, levels = CH), model_lab = factor(MODL[model], levels = MODL[MODELS]))]
  e2[, `:=`(lab = f1(gm_pct), vj = fifelse(gm_pct >= 0, -0.35, 1.35))]
  g2 <- ggplot(e2, aes(change, gm_pct, fill = scen)) + geom_hline(yintercept = 0, colour = VIZ$ink2, linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, colour = VIZ$surface, linewidth = 0.5) +
    geom_text(aes(label = lab, vjust = vj), position = position_dodge(width = 0.78), size = 2.3, colour = VIZ$ink2) +
    facet_grid(pair ~ model_lab, labeller = labeller(pair = label_wrap_gen(26))) + scale_fill_manual(values = c(VIZ$s1, VIZ$s2), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = NULL, y = "Change in the geometric mean (%)", title = sprintf("Exposure change when Vm, ke or F is changed by 20%% (%g kg)", WT_HEAD),
         subtitle = sprintf("Same %s subjects, between-subject variability included; one parameter changed at a time in both scenarios", format(N, big.mark = ",")),
         caption = sprintf("AUC0-last = AUC from the dose to the earlier of the true Tlast at the LLOQ and study Day %g. The B reference includes the washout of drug from all earlier doses and is not an NCA endpoint;\nthe steady-state counterpart of the single-dose AUC0-inf is AUCtau,ss.", T_END + 1)) +
    theme_dupi() + theme(strip.text.y = element_text(size = 8), plot.subtitle = element_text(size = 8.5), plot.caption = element_text(size = 7.5, colour = VIZ$ink2, hjust = 0),
                         plot.margin = margin(6, 14, 6, 6))
  ggsave(file.path(OD, "tm_fig2_elasticity.png"), g2, width = 10.5, height = 8, dpi = 150, bg = "white")
  fwrite(e2[, .(model, scenario, WT, change, metric, pair, gm_pct)], file.path(OD, "tm_fig2_data.csv"))
  # ---------------- 요약(영문, 결과 파일에서 생성) ----------------
  g <- function(m_, q, w = WT_HEAD) { r <- T1[model == m_ & quantity == q & WT == w]; stopifnot(nrow(r) == 1); r }
  ftxt <- function(r, u = "") sprintf("%s%s (typical %s; 5th to 95th percentile %s to %s)", f1(r$median), u, f1(r$typical), f1(r$p05), f1(r$p95))
  qd <- function(x) sprintf("Study Day from which f_TM stays at or above %d%%, all subjects", x)
  qs <- "Share of AUC0-last accumulated while f_TM is at or above 50%"
  rep_ <- c("# Target-mediated elimination share (addendum 1, plan section 17, D-071)", "",
            sprintf("Descriptive model-based quantities; they change no registered result. Model 1 (Kovalenko 2020 Table 1 Model 1) is the primary model and Kovalenko 2016 the sensitivity model; structure and estimates are as published. Medians over the simulated subjects at %g kg unless stated; study Day = days after the dose + 1.", WT_HEAD), "")
  for (m_ in MODELS) rep_ <- c(rep_, sprintf("- %s: C* %s. After %g mg the lowest f_TM (at Cmax) is %s. f_TM stays at or above 50%% from study Day %s, 80%% from Day %s and 95%% from Day %s (Day 1.0 = from the dose). The share of AUC0-last accrued while f_TM is at or above 50%% is %s, of which %s after Cmax. At steady state of regimen B, f_TM is %s at Ctrough,ss, %s at Cavg,ss and %s at Cmax,ss.",
      MODL[[m_]], ftxt(g(m_, "C* = Vm/ke - Km (f_TM = 50%)"), " mg/L"), DOSE$A$a, ftxt(g(m_, "Minimum f_TM (at Cmax)"), "%"), ftxt(g(m_, qd(50))), ftxt(g(m_, qd(80))), ftxt(g(m_, qd(95))),
      ftxt(g(m_, qs), "%"), ftxt(g(m_, "Of which after Cmax (falling phase; supplementary)"), "%"), ftxt(g(m_, "f_TM at Ctrough,ss"), "%"), ftxt(g(m_, "f_TM at Cavg,ss"), "%"), ftxt(g(m_, "f_TM at Cmax,ss (lowest f_TM in the interval)"), "%")))
  tl_ <- function(m_) g(m_, "True Tlast at the LLOQ")
  rep_ <- c(rep_, "", sprintf("True Tlast at the LLOQ after %g mg (%g kg; coverage statement quoted Kovalenko 2016 first, plan section 16.3 item 4): Kovalenko 2016 median study Day %s (Model 1 %s); later f_TM values in Figure 1 are at concentrations below the LLOQ, where f_TM approaches Vm / (Vm + ke Km) and depends on the fixed Km. Share of subjects with a quantifiable true concentration at each sampling day: tm_quantifiable_by_day.csv.",
                              DOSE$A$a, WT_HEAD, f1(tl_("k2016")$median), f1(tl_("k2020")$median)))
  sup <- B0[WT == WT_HEAD, .(amt = 100 * median(amt_tm_share), tmean = 100 * median(ftm_tmean)), by = .(model, scenario)]
  nl <- vapply(MODELS, function(m_) f1(g(m_, "AUCtau,ss / single-dose AUC0-inf of scenario A, same subject (supplementary; 100 under linear PK)")$median), "")
  rep_ <- c(rep_, "", sprintf("Supplementary (not requested), %g kg medians: share of the eliminated amount through the target-mediated pathway %s; time-averaged f_TM over the B interval %s; AUCtau,ss as a percentage of the same subject's single-dose AUC0-inf (100 under linear PK) %s.", WT_HEAD,
                              paste(sprintf("%s %s %s%%", MODL[sup$model], fifelse(sup$scenario == "A", sprintf("A days 0 to %g", T_END), sprintf("B days %g to %g", SS[1], SS[2])), f1(sup$amt)), collapse = "; "),
                              paste(sprintf("%s %s%%", MODL[sup[scenario == "B"]$model], f1(sup[scenario == "B"]$tmean)), collapse = "; "), paste(sprintf("%s %s%%", MODL[MODELS], nl), collapse = "; ")),
            sprintf("Numbers are shown to one decimal as requested; Monte Carlo standard errors of the percentiles are in tm_table1.csv (largest for a median: %s study Days among the crossing-day rows, %s percentage points among the f_TM and share rows).",
                    f1(max(T1[unit == "study Day" & is.finite(se_median)]$se_median)), f1(max(T1[unit == "%" & is.finite(se_median)]$se_median))),
            sprintf("Checks (plan section 17.7): %d of %d passed (tm_checks.csv).", sum(CK$pass), nrow(CK)), "",
            "Outputs: tm_table1.md/.csv (items 2, 3 and 5), tm_table2.md/.csv (item 4), tm_fig1_ftm_curves.png (item 1), tm_fig2_elasticity.png (item 4), tm_elasticity_full.csv, tm_quantifiable_by_day.csv, tm_subjects.csv.gz, tm_typical.csv, tm_curves.csv.")
  writeLines(rep_, file.path(OD, "tm_report.md"))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE)
  for (ff in file.path(OD, c("tm_table1.md", "tm_table2.md", "tm_report.md", "tm_table1.csv", "tm_table2.csv", "tm_elasticity_full.csv", "tm_quantifiable_by_day.csv"))) check_english(ff)
  append_run_log(lf, "report done")
}
