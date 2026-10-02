#!/usr/bin/env Rscript
# 66_day57_mc.R — Day 57 이후 정량 가능 비율: 개인차 몬테카를로(1절), 시험 수준(3절), 라벨 비교(4절), 기존 결과 재현(A)와 요약
#   (사전 등록 config/prereg_20261002_day57.yaml, D-066). 모델 구조와 추정값은 그대로 쓴다(재조정 없음).
# 모드
#   individual <k2016|k2020>              1절 + A: 절벽 분석(scripts/34)과 같은 호출·시드로 20,000명을 다시 만들고 저장 결과와 대조, 시점별 정량한계 위 비율
#   trials <k2016|k2020> [from] [to]      3절: 시험 from..to(기본 1..10,000), 250회 묶음 파일(이미 있는 묶음은 건너뜀 = 이어서 실행)
#   label <k2016|k2020> <Q2W|QW>          4절: 아토피 가정 체중, 600 mg 부하 뒤 300 mg Q2W 또는 QW를 day 364까지, 마지막 투여 뒤 정량한계 미만까지 시간
#   summary                               모든 요약 표와 한국어 요약(results/day57/)
# 사용법: nice -n 19 Rscript scripts/66_day57_mc.R <mode> ...
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("individual", "trials", "label", "summary")) stop("사용법: Rscript scripts/66_day57_mc.R <individual|trials|label|summary> ...")
mode <- args[1]
PR <- read_cfg("prereg_20261002_day57.yaml"); design <- read_cfg("trial_design.yaml"); oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff
out_dir <- proj_path(PR$output_dir); dir.create(file.path(out_dir, "trials"), showWarnings = FALSE, recursive = TRUE)
MS <- 20261002L                                                        # 사전 등록 시드 근간(3·4절)
LLOQ <- study_lloq(); stopifnot(isTRUE(all.equal(LLOQ, 0.078)))
DAYS <- as.numeric(unlist(PR$timepoints_study_day)); EL <- DAYS - 1     # 연구일 d = 경과일 d - 1
VAR <- c(k2016 = "base", k2020 = "struct2020")
params_of <- function(model) {                                         # 절벽 분석(load_params)과 시험 엔진(resolve_variant)의 파라미터가 같아야 한다
  rv <- resolve_variant(VAR[[model]], design); p0 <- load_params(model)
  stopifnot(isTRUE(all.equal(rv$p$theta, p0$theta)), isTRUE(all.equal(rv$p$omega, p0$omega)), isTRUE(all.equal(rv$p$sigma, p0$sigma)), rv$p$ada$fraction == 0)
  rv
}

# ---------------------------------------------------------------- 1절 + A ----------------------------------------------------------------
if (mode == "individual") {
  model <- args[2]; stopifnot(model %in% names(VAR))
  N <- as.integer(PR$section1_monte_carlo$n_subjects_per_model); SEED <- as.integer(cf$seed); STEP <- as.numeric(cf$grid_step_day); defs <- as.numeric(unlist(cf$definition_days))
  stopifnot(N == as.integer(cf$n_subjects), SEED == 20260928L)
  lf <- start_run_log(paste0("day57_individual_", model), master_seed = SEED, run_mode = "final", extra = list(n = N, step = STEP))
  rv <- params_of(model); p <- rv$p
  t0 <- Sys.time()
  cs <- cliff_subjects(model, "base", N, SEED, STEP, LLOQ, defs, cf$weights[["base"]], design)   # scripts/34와 같은 호출
  cs[, lloq := NULL]
  append_run_log(lf, sprintf("cliff_subjects done in %s", format(Sys.time() - t0)))
  # A: 저장 결과와 대조(대상자별 t_lloq, 요약)
  mdl <- model; st <- readRDS(proj_path("results", "cliff", "cliff_subjects.rds"))[model == mdl & weight == "base"]
  m <- merge(cs[, .(id, t_new = t_lloq, tmax_new = tmax, len1_new = len1)], st[, .(id, t_old = t_lloq, tmax_old = tmax, len1_old = len1)], by = "id")
  stopifnot(nrow(m) == N)
  dif <- function(a, b) max(abs(a - b), na.rm = TRUE)
  ss <- fread(proj_path("results", "cliff", "cliff_summary.csv"))[model == mdl & weight == "base"]; stopifnot(nrow(ss) == 1)
  q <- function(x, pr) quantile(x + 1, pr, na.rm = TRUE, names = FALSE)
  rep_rows <- rbind(
    data.table(item = "per-subject t_lloq, max abs difference (days)", stored = NA_real_, new = dif(m$t_new, m$t_old), identical = identical(m$t_new, m$t_old)),
    data.table(item = "per-subject tmax, max abs difference (days)", stored = NA_real_, new = dif(m$tmax_new, m$tmax_old), identical = identical(m$tmax_new, m$tmax_old)),
    data.table(item = "LLOQ reach study day, median", stored = ss$lloq_studyday_median, new = q(cs$t_lloq, 0.5), identical = identical(ss$lloq_studyday_median, q(cs$t_lloq, 0.5))),
    data.table(item = "LLOQ reach study day, 5th percentile", stored = ss$lloq_studyday_p05, new = q(cs$t_lloq, 0.05), identical = identical(ss$lloq_studyday_p05, q(cs$t_lloq, 0.05))),
    data.table(item = "LLOQ reach study day, 95th percentile", stored = ss$lloq_studyday_p95, new = q(cs$t_lloq, 0.95), identical = identical(ss$lloq_studyday_p95, q(cs$t_lloq, 0.95))),
    data.table(item = "share reaching the LLOQ after study Day 58 (%)", stored = ss$lloq_after_day58_pct, new = 100 * mean(cs$t_lloq + 1 > 58, na.rm = TRUE),
               identical = identical(ss$lloq_after_day58_pct, 100 * mean(cs$t_lloq + 1 > 58, na.rm = TRUE))))
  rep_rows[, `:=`(pk_model = model, diff = new - stored)]
  rep_rows[, agree_csv_precision := fifelse(is.na(stored), identical, abs(diff) <= 1e-12 * pmax(1, abs(stored)))]   # 저장 요약은 CSV(유효숫자 15자리) → 1e-12 상대 허용
  # 1절: 시점별 정량한계 위 비율(도달 시각 지표) + 직접 풀이 지표 대조
  subj <- with_seed(derive_seed(SEED, model, "base", "subj"), make_subjects(N, p, cliff_weight_spec(design, "base"), 0.5, 0))
  ip <- individual_params(p, subj)
  sol <- solve_model(ip, CJ(id = ip$id, time = EL), design$dose_mg, model_id = p$model_id)
  dd <- merge(sol[, .(id, el = time, C)], cs[, .(id, t_lloq)], by = "id")
  dd[, `:=`(above_direct = C >= LLOQ, above_reach = is.na(t_lloq) | t_lloq > el)]
  SH <- dd[, { w <- wilson_ci(sum(above_reach), .N); .(n = .N, n_above = sum(above_reach), pct = w$est, lo = w$lo, hi = w$hi,   # wilson_ci는 %로 반환
                                                      n_disagree_direct = sum(above_direct != above_reach)) }, by = el][, `:=`(pk_model = model, study_day = el + 1)][]
  RP <- data.table(pk_model = model, n = N, n_never = sum(is.na(cs$t_lloq)), p50 = q(cs$t_lloq, 0.5), p90 = q(cs$t_lloq, 0.9), p95 = q(cs$t_lloq, 0.95), p99 = q(cs$t_lloq, 0.99),
                   max = max(cs$t_lloq + 1, na.rm = TRUE))
  fwrite(SH, file.path(out_dir, sprintf("day57_share_by_day_%s.csv", model))); fwrite(RP, file.path(out_dir, sprintf("day57_reach_percentiles_%s.csv", model)))
  fwrite(rep_rows, file.path(out_dir, sprintf("day57_repro_check_%s.csv", model)))
  fwrite(cs[, .(id, WT, tmax, t_lloq)], file.path(out_dir, sprintf("day57_reach_subjects_%s.csv.gz", model)))
  print(SH); print(RP); print(rep_rows)
  append_run_log(lf, sprintf("done: disagreements direct vs reach %d; t_lloq identical %s", sum(SH$n_disagree_direct), rep_rows$identical[1]))
}

# ---------------------------------------------------------------- 3절 ----------------------------------------------------------------
if (mode == "trials") {
  model <- args[2]; stopifnot(model %in% names(VAR))
  S3 <- PR$section3_trial_level; NT <- as.integer(S3$n_trials); stopifnot(NT == 10000L)
  from <- if (length(args) >= 3) as.integer(args[3]) else 1L; to <- if (length(args) >= 4) as.integer(args[4]) else NT
  stopifnot(from >= 1, to <= NT, to >= from)
  rv <- params_of(model); p <- rv$p; wt <- rv$wt_spec; n_arm <- as.integer(design$n_per_arm); stopifnot(n_arm == 117L)
  st <- strata_from_weight_spec(wt, as.numeric(design$stratification$split_kg$value))
  EL3 <- sort(unique(c(57, as.numeric(unlist(PR$additional_sampling_study_day)) - 1)))    # Day 58(참고), 64, 71, 85
  lf <- start_run_log(sprintf("day57_trials_%s_%d_%d", model, from, to), master_seed = MS, run_mode = "final", extra = list(model = model, from = from, to = to))
  B <- 250L
  for (b0 in seq((from - 1) %/% B * B + 1, to, by = B)) {
    js <- max(b0, from):min(b0 + B - 1, to)
    f <- file.path(out_dir, "trials", sprintf("day57_trials_%s_%05d_%05d.csv.gz", model, min(js), max(js)))
    if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
    t0 <- Sys.time()
    lst <- lapply(js, function(j) {
      s <- with_seed(derive_seed(MS, "trial3", model, j, "subj"), make_subjects(2 * n_arm, p, wt, design$weight$sex_ratio_male$value, p$ada$fraction))
      s <- with_seed(derive_seed(MS, "trial3", model, j, "alloc"), assign_arms_stratified(s, st$breaks, st$labels))
      e <- rbindlist(lapply(c("R", "T"), function(a) with_seed(derive_seed(MS, "trial3", model, j, a, "eps"), draw_eps(s[arm == a]$id, EL3, p$sigma))))
      uid0 <- (j - 1L) * (2L * n_arm)
      list(s = s[, `:=`(trial = j, uid = uid0 + id)], e = e[, `:=`(trial = j, uid = uid0 + id)])
    })
    subj <- rbindlist(lapply(lst, `[[`, "s")); eps <- rbindlist(lapply(lst, `[[`, "e"))
    ip <- individual_params(p, copy(subj)[, id := uid])
    sol <- solve_model(ip, CJ(id = ip$id, time = EL3), design$dose_mg, model_id = p$model_id)
    d <- merge(sol[, .(uid = id, planned = time, C)], subj[, .(uid, trial, arm, stratum)], by = "uid")
    d <- merge(d, eps[, .(uid, planned, eps_p, eps_a)], by = c("uid", "planned"))
    d[, y := C * (1 + eps_p) + eps_a]
    res <- d[, .(n = .N, n_true = sum(C >= LLOQ), n_obs = sum(y >= LLOQ)), by = .(trial, arm, study_day = planned + 1)]
    stopifnot(all(res$n == n_arm), nrow(res) == length(js) * 2L * length(EL3))
    tmp <- file.path(dirname(f), paste0("tmp_", basename(f))); fwrite(res[, pk_model := model], tmp); file.rename(tmp, f)   # 완성된 묶음만 정식 이름
    msg <- sprintf("%s trials %d-%d written in %s", model, min(js), max(js), format(Sys.time() - t0)); cat(msg, "\n"); append_run_log(lf, msg)
    rm(lst, subj, eps, ip, sol, d); gc(FALSE)
  }
  append_run_log(lf, "done")
}

# ---------------------------------------------------------------- 4절 ----------------------------------------------------------------
if (mode == "label") {
  model <- args[2]; reg <- args[3]; stopifnot(model %in% names(VAR), reg %in% c("Q2W", "QW"))
  S4 <- PR$section4_label; N <- as.integer(S4$simulation$n_subjects); stopifnot(N == 10000L)
  pa <- read_cfg("population_atopic.yaml")$distributions$primary
  wt <- list(dist = "lognormal", mean = as.numeric(pa$mean), sd = as.numeric(pa$sd), trunc = as.numeric(unlist(pa$trunc)))
  stopifnot(wt$mean == 78, wt$sd == 19, identical(wt$trunc, c(40, 180)))
  rv <- params_of(model); p <- rv$p; mod <- get_model(p$model_id)
  tau <- c(Q2W = 14, QW = 7)[[reg]]; t_last <- 364
  dts <- c(0, seq(tau, t_last, by = tau)); amts <- c(600, rep(300, length(dts) - 1)); stopifnot(max(dts) == t_last)
  grid <- t_last + seq(0.1, 300, by = 0.1); trough_t <- c(t_last - tau - 1e-3, t_last - 1e-3)
  lf <- start_run_log(sprintf("day57_label_%s_%s", model, reg), master_seed = MS, run_mode = "final", extra = list(model = model, regimen = reg, n = N, doses = length(dts)))
  subj <- with_seed(derive_seed(MS, "label", model, reg), make_subjects(N, p, wt, 0.5, 0))
  ip <- individual_params(p, subj)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  out <- list(); chunk <- 1000L; t0 <- Sys.time()
  for (b in seq_len(ceiling(N / chunk))) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, N)
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    ev <- rxode2::add.sampling(rxode2::et(amt = amts, time = dts, cmt = "depot"), c(trough_t, grid))
    s <- rxode2::rxSolve(mod, params = prm, events = ev, atol = 1e-10, rtol = 1e-8, maxsteps = 2000000L, cores = 1L, returnType = "data.table", addDosing = FALSE)
    s[, id := prm$id[sim.id]]
    out[[b]] <- s[, {
      tr <- vapply(trough_t, function(x) C[which.min(abs(time - x))], numeric(1)); g <- time > t_last; tt <- time[g] - t_last; cc <- C[g]
      k <- which(cc < LLOQ)[1]
      tl <- if (is.na(k)) NA_real_ else if (k == 1) tt[1] else tt[k - 1] + (log(LLOQ) - log(cc[k - 1])) / (log(cc[k]) - log(cc[k - 1])) * (tt[k] - tt[k - 1])
      .(trough_prev = tr[1], trough_last = tr[2], c_last_plus = cc[1], t_below_days = tl)
    }, by = id]
  }
  R <- rbindlist(out)[, `:=`(pk_model = model, regimen = reg, WT = subj$WT[match(id, subj$id)])]
  fwrite(R, file.path(out_dir, sprintf("day57_label_subjects_%s_%s.csv.gz", model, reg)))
  msg <- sprintf("%s %s: n %d, median %.2f weeks, NA %d, trough ratio median %.4f, in %s", model, reg, nrow(R), median(R$t_below_days, na.rm = TRUE) / 7, sum(is.na(R$t_below_days)),
                 median(R$trough_last / R$trough_prev), format(Sys.time() - t0))
  cat(msg, "\n"); append_run_log(lf, msg)
}

# ---------------------------------------------------------------- 요약 ----------------------------------------------------------------
if (mode == "summary") {
  PK <- unlist(PR$pk_models)
  rd <- function(pat) rbindlist(lapply(PK, function(m) { f <- file.path(out_dir, sprintf(pat, m)); if (file.exists(f)) fread(f) else stop("없음: ", f) }))
  SH <- rd("day57_share_by_day_%s.csv"); RP <- rd("day57_reach_percentiles_%s.csv"); RC <- rd("day57_repro_check_%s.csv")
  stopifnot(all(SH$n_disagree_direct == 0))
  fwrite(SH[, .(pk_model, study_day, n, n_above, pct, lo, hi, n_disagree_direct)], file.path(out_dir, "day57_share_by_day.csv"))
  fwrite(RP, file.path(out_dir, "day57_reach_percentiles.csv")); fwrite(RC[, .(pk_model, item, stored, new, diff, identical, agree_csv_precision)], file.path(out_dir, "day57_repro_check.csv"))
  stopifnot(all(RC$agree_csv_precision))
  # 3절
  TF <- list.files(file.path(out_dir, "trials"), pattern = "^day57_trials_.*\\.csv\\.gz$", full.names = TRUE)
  TR <- rbindlist(lapply(TF, fread))
  chk <- TR[, .(n_trials = uniqueN(trial)), by = pk_model]; print(chk)
  stopifnot(all(chk$n_trials == as.integer(PR$section3_trial_level$n_trials)), nrow(TR) == uniqueN(TR[, .(pk_model, trial, arm, study_day)]))
  smry <- function(x, unit, meas) data.table(unit = unit, measure = meas, n = length(x), median = median(x), p95 = quantile(x, 0.95, names = FALSE, type = 7),
                                             mean = mean(x), share_ge1_pct = 100 * mean(x >= 1), max = max(x))
  TS <- rbindlist(lapply(PK, function(m) rbindlist(lapply(sort(unique(TR$study_day)), function(dy) {
    a <- TR[pk_model == m & study_day == dy]; tr <- a[, .(t_true = sum(n_true), t_obs = sum(n_obs)), by = trial]
    rbind(smry(a$n_true, "per arm", "true"), smry(a$n_obs, "per arm", "observed"), smry(tr$t_true, "per trial (234)", "true"), smry(tr$t_obs, "per trial (234)", "observed"))[, `:=`(pk_model = m, study_day = dy)] }))))
  setcolorder(TS, c("pk_model", "study_day", "unit", "measure"))
  fwrite(TS, file.path(out_dir, "day57_trials_summary.csv"))
  # 4절
  LB <- rbindlist(lapply(PK, function(m) rbindlist(lapply(c("Q2W", "QW"), function(r) fread(file.path(out_dir, sprintf("day57_label_subjects_%s_%s.csv.gz", m, r)))))))
  lab <- PR$section4_label$label; LW <- c(Q2W = as.numeric(lab$q2w_300mg_weeks), QW = as.numeric(lab$qw_300mg_weeks))
  LS <- LB[, .(n = .N, n_never = sum(is.na(t_below_days)), median_weeks = median(t_below_days, na.rm = TRUE) / 7, p25_weeks = quantile(t_below_days, 0.25, na.rm = TRUE, names = FALSE) / 7,
               p75_weeks = quantile(t_below_days, 0.75, na.rm = TRUE, names = FALSE) / 7, trough_ratio_median = median(trough_last / trough_prev)), by = .(pk_model, regimen)]
  LS[, `:=`(label_weeks = LW[regimen], diff_pct = 100 * (median_weeks - LW[regimen]) / LW[regimen])]
  LS[, interpretation := fifelse(abs(diff_pct) <= 20, "within_20pct", fifelse(diff_pct < 0, "model_shorter", "model_longer"))]
  stopifnot(all(abs(LS$trough_ratio_median - 1) <= 0.01))           # 정상상태 확인(사전 등록)
  fwrite(LS, file.path(out_dir, "day57_label_comparison.csv"))
  print(SH[study_day %in% c(58, 64, 71, 85)]); print(RP); print(RC); print(TS[measure == "true"]); print(LS)
  # 한국어 요약(수치는 위 표에서 자동 생성). 2절(모델 불확실성)은 입력 대기라 자리만 둔다.
  f1 <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
  ML <- c(k2020 = "2020 모델(주 모델)", k2016 = "2016 모델(민감도 모델)")
  L <- c("# Day 57 이후 정량 가능 비율 (D-066, config/prereg_20261002_day57.yaml)", "",
         "모델 구조와 추정값은 그대로 썼다(재조정 없음). 수치는 results/day57/*.csv에서 자동 생성했다.", "",
         "## 1절 몬테카를로(모델당 20,000명, 시험 모집단 60~90 kg, 300 mg 단회 피하)", "")
  for (m in PK) { x <- SH[pk_model == m][order(study_day)]; r <- RP[pk_model == m]
    L <- c(L, sprintf("- %s: 참 농도가 정량한계 위인 비율 %s", ML[[m]], paste(sprintf("Day %d %s%%(%s~%s)", x$study_day, f1(x$pct, 2), f1(x$lo, 2), f1(x$hi, 2)), collapse = ", ")),
           sprintf("  - 정량한계 도달 연구일: 중앙값 Day %s, 90번째 백분위수 %s, 95번째 %s, 99번째 %s", f1(r$p50), f1(r$p90), f1(r$p95), f1(r$p99))) }
  L <- c(L, "", "## A. 기존 결과 재현", "")
  for (m in PK) { r <- RC[pk_model == m]
    L <- c(L, sprintf("- %s: 대상자별 정량한계 도달 시각이 저장 결과(results/cliff/cliff_subjects.rds)와 완전히 같다(최대 차이 %s일). Day 58 이후 비율 %s%%(저장값 %s%%). 요약 분위수의 차이(최대 %.1e일)는 저장 요약이 CSV(유효숫자 15자리)라서 생긴 표기 차이다.",
                      ML[[m]], format(r[grepl("t_lloq", item)]$new), format(r[grepl("Day 58", item)]$new), format(r[grepl("Day 58", item)]$stored), max(abs(r$diff), na.rm = TRUE))) }
  L <- c(L, "", "## 3절 시험 수준(arm당 평가 가능 117명, 체중 층 배정, 시험 10,000회; 참 농도 기준)", "")
  for (m in PK) for (dy in sort(unique(TS$study_day))) { a <- TS[pk_model == m & study_day == dy & measure == "true"]
    L <- c(L, sprintf("- %s Day %d: arm당 중앙값 %s명(95번째 백분위수 %s명, 최대 %s명), 시험 전체(234명) 중앙값 %s명(95번째 %s명); 1명 이상인 arm %s%%",
                      ML[[m]], dy, f1(a[unit == "per arm"]$median, 0), f1(a[unit == "per arm"]$p95, 0), f1(a[unit == "per arm"]$max, 0),
                      f1(a[unit != "per arm"]$median, 0), f1(a[unit != "per arm"]$p95, 0), f1(a[unit == "per arm"]$share_ge1_pct, 1))) }
  L <- c(L, "", "## 4절 라벨 비교(아토피 가정 체중, 600 mg 부하 뒤 day 364까지 반복 투여, 질환 공변량 없음)", "")
  for (i in seq_len(nrow(LS))) { r <- LS[i]
    L <- c(L, sprintf("- %s %s: 마지막 투여 뒤 정량한계(78 ng/mL) 미만까지 중앙값 %s주(사분위 %s~%s주) 대 라벨 %s주, 차이 %s%% → %s",
                      ML[[r$pk_model]], r$regimen, f1(r$median_weeks), f1(r$p25_weeks), f1(r$p75_weeks), f1(r$label_weeks, 0), f1(r$diff_pct),
                      c(within_20pct = "±20% 이내: 정량한계 근처 꼬리 예측이 외부 자료와 맞는다", model_shorter = "모델이 짧다: Day 58 이후 비율이 과소추정일 수 있다(한계로만 적음)",
                        model_longer = "모델이 길다: 과대추정 방향(보수적)")[[r$interpretation]])) }
  L <- c(L, "", "## 2절·3절 중첩·B (입력 대기)", "", "- 논문의 상대 표준오차(config/param_uncertainty.yaml)가 없어 실행하지 않았다(D-066).")
  writeLines(L, file.path(out_dir, "day57_summary_ko.md"))
}
