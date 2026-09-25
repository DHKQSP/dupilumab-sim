#!/usr/bin/env Rscript
# 55_trial_population_summary.R — 시험 모집단(건강인 60–90 kg, 체중 층화, B0, 300 mg) 안에서 AUC0-inf 제외 논리의 근거 (정정 지시 2026-09-26 §2, §3; prereg section6, D-059).
# 입력: 저장된 20,000명 NCA(results/individual/nca_base_20000.rds = 2016 모델, nca_struct2020_20000.rds = Model 1, nca_resid12_20000.rds = 2016 모델 비례 잔차 12%),
#       section1 재생성 시험(results/oc_models/oc_models_be_<model>.csv.gz, oc_models_crit_be_<model>.csv.gz; M0 행의 n_R·n_T, AUClast GMR),
#       scripts/54 산출(results/trialpop/trialpop_counts_<model>.csv.gz, trialpop_truth_<model>.csv), 역산 참값(results/oc/inversion_all.csv).
# 산출(results/trialpop/):
#   tp_failure_by_set.csv          §2-1 모델 × 세트: 미달 %, λz 산출 불가 %, 산출되나 미달 %(Wilson 95%)
#   tp_failure_reasons.csv         §2-1 모델 × 세트 × 사유: 위계 분해(산출 불가 → adjusted R² → 외삽 → span)와 비배타 플래그 비율
#   tp_retained_per_arm.csv        §2-1 모델 × 세트: section1 S00 시험 1–10,000의 arm당 잔류 인원(중앙값, 5·95백분위, 평균)
#   tp_residual_sensitivity.csv    §2-1 잔차 민감성: 2016 모델(비례 24.2%), 2016 모델 비례 12%, Model 1(15.0%)
#   tp_strata_individual.csv       §2-2 모델 × 세트 × 층: 미달 %(Wilson), 무거운 층 − 가벼운 층 차이(Newcombe 95%)
#   tp_strata_composition.csv      §2-2 모델 × 시나리오 × 분석군(AUClast, 규칙 A 세트 i–iv): 무거운 층 비율의 무작위배정 대비 편차와 arm 간 차이(시험 간 분포)
#   tp_arm_difference.csv          §2-3 모델 × 시나리오 × 세트: 시험군 − 대조군 미달 %p(평균, 95% 구간, 5·95백분위), 참 AUC0-inf 비
#   tp_arm_difference_fit.csv      §2-3 모델 × 세트: 평균 차이 ~ log(참 비) 선형 적합(기술 통계)
#   tp_characteristics.csv         §2-4 모델 × 세트: 미달 − 잔류 체중 차이, 참 AUC0-inf·AUClast 기하평균비(Welch 95%)
#   tp_coverage_individual.csv     §3 모델 × 지표: 채혈 구간 커버리지, 관측 대 참 비(AUClast, AUCinf 규칙 B, 규칙 A 세트 i), 외삽 항
#   tp_gmr_agreement.csv           §3 section1 M0: AUClast GMR 기하평균 대 참 AUC0-inf 비(경계 16칸, S00, F097)
#   tp_identity_check.csv          재생성 시험의 군별 인원이 section1 n_R·n_T와 같은지(모델별 비교 행 수)
#   fig_tp_failure_by_set.png, fig_tp_arm_difference.png, fig_tp_strata.png, tp_conclusion_en.md, tp_conclusion_ko.md
# 사용법: Rscript scripts/55_trial_population_summary.R   (환경변수 DUPI_TRIALPOP_SUMMARY_OUT, DUPI_TRIALPOP_ALLOW_PARTIAL=1은 시험용)
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
pr6 <- read_cfg("prereg_20260926.yaml")$section6; pr1 <- read_cfg("prereg_20260926.yaml")$section1; oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
tp_dir <- proj_path("results", "trialpop"); out_dir <- Sys.getenv("DUPI_TRIALPOP_SUMMARY_OUT", tp_dir); tp_in <- Sys.getenv("DUPI_TRIALPOP_IN", tp_dir)   # 환경변수는 시험용; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
allow_partial <- nzchar(Sys.getenv("DUPI_TRIALPOP_ALLOW_PARTIAL", ""))
om_dir <- proj_path(pr1$out_dir); PK <- c("k2016", "k2020"); SETS <- names(CRIT_SETS)
MODEL_EN <- c(k2016 = "2016 model", k2020 = "Model 1"); SET_EN <- c(i = "(i)", ii = "(ii)", iii = "(iii)", iv = "(iv)")
NMAP <- c(i = "n_i", ii = "n_ii", iii = "n_iii", iv = "n_iv"); EPA <- c(i = "AUCinf_Ai", ii = "AUCinf_A", iii = "AUCinf_Aiii", iv = "AUCinf_Aiv")
split_kg <- as.numeric(design$stratification$split_kg$value); NT1 <- as.integer(oc$trials$reps_boundary); NTR <- as.integer(pr6$n_trials_regen)
f1 <- function(x) formatC(round(x + sign(x) * 1e-9, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x + sign(x) * 1e-9, 2) + 0, format = "f", digits = 2)
f3 <- function(x) formatC(round(x + sign(x) * 1e-9, 3) + 0, format = "f", digits = 3)
rg <- function(x, f = f1, u = "%", sep = " to ") if (identical(f(min(x)), f(max(x)))) sprintf("%s%s", f(min(x)), u) else sprintf("%s%s%s%s%s", f(min(x)), u, sep, f(max(x)), u)
wil <- function(k, n) { w <- wilson_ci(k, n); c(est = w$est, lo = w$lo, hi = w$hi) }
newcombe <- function(k1, n1, k2, n2) { a <- wil(k1, n1); b <- wil(k2, n2); d <- b[["est"]] - a[["est"]]      # 차이 = 2 − 1 (%p)
  c(diff = d, lo = d - sqrt((b[["est"]] - b[["lo"]])^2 + (a[["hi"]] - a[["est"]])^2), hi = d + sqrt((b[["hi"]] - b[["est"]])^2 + (a[["est"]] - a[["lo"]])^2)) }
welch <- function(x, y) { t <- t.test(x, y); c(diff = unname(t$estimate[1] - t$estimate[2]), lo = t$conf.int[1], hi = t$conf.int[2]) }
IND <- list(k2016 = readRDS(proj_path("results", "individual", "nca_base_20000.rds"))[schedule == "B0"], k2020 = readRDS(proj_path("results", "individual", "nca_struct2020_20000.rds"))[schedule == "B0"])
R12 <- readRDS(proj_path("results", "individual", "nca_resid12_20000.rds"))[schedule == "B0"]
for (m in PK) { x <- IND[[m]]; stopifnot(nrow(x) == 20000L, min(x$WT) >= 60, max(x$WT) <= 90); x[, stratum := fifelse(WT <= split_kg, 1L, 2L)] }
sig <- c(k2016 = resolve_variant("base", design)$p$sigma[["prop"]], k2020 = resolve_variant("struct2020", design)$p$sigma[["prop"]], resid12 = resolve_variant("resid12", design)$p$sigma[["prop"]])

# ---- §2-1 기준 세트별 미달 --------------------------------------------------------------------------------------------------
thr <- function(s_) CRIT_SETS[[s_]]
FB <- rbindlist(lapply(PK, function(m) { x <- IND[[m]]; rbindlist(lapply(SETS, function(s_) { ok <- crit_ok(x, s_); lz <- !(x$lambda_ok %in% TRUE); N <- nrow(x)
  a <- wil(sum(!ok), N); b <- wil(sum(lz), N); c_ <- wil(sum(!ok & !lz), N)
  data.table(pk_model = m, set = s_, n = N, fail_pct = a[["est"]], fail_lo = a[["lo"]], fail_hi = a[["hi"]], lz_pct = b[["est"]], lz_lo = b[["lo"]], lz_hi = b[["hi"]],
             est_fail_pct = c_[["est"]], est_fail_lo = c_[["lo"]], est_fail_hi = c_[["hi"]], fail_per_arm = 117 * mean(!ok)) })) }))
fwrite(FB, file.path(out_dir, "tp_failure_by_set.csv"))
RS <- rbindlist(lapply(PK, function(m) { x <- IND[[m]]; lz <- x$lambda_ok %in% TRUE
  rbindlist(lapply(SETS, function(s_) { t_ <- thr(s_); r2 <- lz & !(is.finite(x$Rsq_adjusted) & x$Rsq_adjusted >= t_$adj_r2_min)
    ex <- lz & !(is.finite(x[["AUC_%Extrap_obs"]]) & x[["AUC_%Extrap_obs"]] <= t_$extrap_max_pct)
    sp <- if (is.na(t_$span_ratio_min)) rep(FALSE, nrow(x)) else lz & !(is.finite(x$Span_ratio) & x$Span_ratio >= t_$span_ratio_min)
    stopifnot(identical(!crit_ok(x, s_), !lz | r2 | ex | sp))
    data.table(pk_model = m, set = s_, reason = c("lambda-z not estimable", "adjusted R-squared below threshold", "extrapolation above 20%", "span ratio below threshold"),
               hierarchical_pct = 100 * c(mean(!lz), mean(r2), mean(ex & !r2), mean(sp & !r2 & !ex)), any_pct = 100 * c(mean(!lz), mean(r2), mean(ex), mean(sp))) })) }))
fwrite(RS, file.path(out_dir, "tp_failure_reasons.csv"))
S1 <- rbindlist(lapply(PK, function(m) rbind(fread(file.path(om_dir, sprintf("oc_models_be_%s.csv.gz", m))), fread(file.path(om_dir, sprintf("oc_models_crit_be_%s.csv.gz", m))))[model == "M0"][, pk_model := m]))
nt1 <- S1[scenario == "S00" & endpoint == "AUClast", uniqueN(trial), by = pk_model]
if (any(nt1$V1 < NT1) && !allow_partial) stop("section1 S00 시험 수 부족: ", paste(nt1[V1 < NT1, sprintf("%s %d", pk_model, V1)], collapse = ", "))
RP <- rbindlist(lapply(PK, function(m) rbindlist(lapply(c(SETS, "lambda"), function(s_) { ep <- if (s_ == "lambda") "AUCinf_B" else EPA[[s_]]
  x <- S1[pk_model == m & scenario == "S00" & endpoint == ep & trial <= NT1]; v <- c(x$n_R, x$n_T)
  data.table(pk_model = m, set = s_, n_trials = nrow(x), n_arms = length(v), retained_median = median(v), retained_p05 = quantile(v, 0.05, names = FALSE), retained_p95 = quantile(v, 0.95, names = FALSE),
             retained_mean = mean(v), retained_min = min(v), retained_max = max(v)) }))))
fwrite(RP, file.path(out_dir, "tp_retained_per_arm.csv"))
RSD <- rbindlist(lapply(list(k2016 = IND$k2016, resid12 = R12, k2020 = IND$k2020), function(x) rbindlist(lapply(SETS, function(s_) { a <- wil(sum(!crit_ok(x, s_)), nrow(x))
  data.table(set = s_, fail_pct = a[["est"]], fail_lo = a[["lo"]], fail_hi = a[["hi"]], adj_r2_median = median(x$Rsq_adjusted, na.rm = TRUE)) }))), idcol = "variant")
RSD[, `:=`(sigma_prop_pct = 100 * sig[variant], label = c(k2016 = "2016 model", resid12 = "2016 model, proportional residual 12%", k2020 = "Model 1")[variant])]
fwrite(RSD, file.path(out_dir, "tp_residual_sensitivity.csv"))

# ---- §2-2 층별 ------------------------------------------------------------------------------------------------------------------
SI <- rbindlist(lapply(PK, function(m) { x <- IND[[m]]; rbindlist(lapply(SETS, function(s_) { f <- !crit_ok(x, s_); k1 <- sum(f[x$stratum == 1]); n1 <- sum(x$stratum == 1); k2 <- sum(f[x$stratum == 2]); n2 <- sum(x$stratum == 2)
  a <- wil(k1, n1); b <- wil(k2, n2); d <- newcombe(k1, n1, k2, n2)
  data.table(pk_model = m, set = s_, n_light = n1, fail_light_pct = a[["est"]], fail_light_lo = a[["lo"]], fail_light_hi = a[["hi"]], n_heavy = n2, fail_heavy_pct = b[["est"]], fail_heavy_lo = b[["lo"]], fail_heavy_hi = b[["hi"]],
             diff_pp = d[["diff"]], diff_lo = d[["lo"]], diff_hi = d[["hi"]]) })) }))
fwrite(SI, file.path(out_dir, "tp_strata_individual.csv"))
CNT <- rbindlist(lapply(PK, function(m) { f <- file.path(tp_in, sprintf("trialpop_counts_%s.csv.gz", m)); if (!file.exists(f)) { if (allow_partial) return(NULL) else stop("없음: ", f) }; fread(f)[, pk_model := m] }))
ntr <- CNT[, uniqueN(trial), by = pk_model]
if ((nrow(ntr) < 2 || any(ntr$V1 < NTR)) && !allow_partial) stop("scripts/54 시험 수 부족")
# 재검사: 재생성 인원(층 합계)이 section1 M0 행의 n_R·n_T와 같은가(scripts/54가 묶음마다 검사한 것을 산출 파일 전체에서 다시 확인)
EPM <- c(AUClast = "n_auclast", AUCinf_B = "n_lambda", AUCinf_Ai = "n_i", AUCinf_A = "n_ii", AUCinf_Aiii = "n_iii", AUCinf_Aiv = "n_iv")
if (nrow(CNT)) {
  tt <- melt(CNT[, lapply(.SD, sum), by = .(pk_model, trial, scenario), .SDcols = unname(EPM)], id.vars = c("pk_model", "trial", "scenario"), variable.name = "col", value.name = "n_new")[, col := as.character(col)]
  s1c <- S1[endpoint %in% names(EPM) & trial %in% unique(CNT$trial)]
  s1c <- rbind(unique(s1c[, .(pk_model, trial, scenario = "REF", col = unname(EPM[endpoint]), n = n_R)]), s1c[, .(pk_model, trial, scenario, col = unname(EPM[endpoint]), n = n_T)])
  if (anyDuplicated(s1c[, .(pk_model, trial, scenario, col)])) stop("section1 대조군 인원이 시나리오마다 다릅니다")
  mm <- merge(tt, s1c, by = c("pk_model", "trial", "scenario", "col"))
  IDC <- mm[, .(rows_compared = .N, identical_rows = sum(n_new == n), trials = uniqueN(trial), scenarios = uniqueN(scenario)), by = pk_model][, ok := rows_compared == identical_rows]
  fwrite(IDC, file.path(out_dir, "tp_identity_check.csv")); if (!all(IDC$ok)) stop("재생성 인원이 section1과 다릅니다")
}
SETSX <- c(auclast = "n_auclast", setNames(NMAP, SETS))
# 시험 × 시나리오 × 분석군: arm별 무거운 층 비율
comp <- function(d, col) { w <- dcast(d, pk_model + trial + scenario ~ stratum, value.var = col); setnames(w, c("1", "2"), c("a1", "a2")); w[, h := 100 * a2 / (a1 + a2)][, .(pk_model, trial, scenario, h)] }
REFC <- CNT[scenario == "REF"]; TC <- CNT[scenario != "REF"]
CP <- rbindlist(lapply(names(SETSX), function(k) {
  hr <- comp(REFC, SETSX[[k]]); hrr <- comp(REFC, "n_rand"); ht <- comp(TC, SETSX[[k]]); htr <- comp(TC, "n_rand")
  z <- ht[hr[, .(pk_model, trial, h_R = h)], on = c("pk_model", "trial")][htr[, .(pk_model, trial, scenario, h_T_rand = h)], on = c("pk_model", "trial", "scenario")][hrr[, .(pk_model, trial, h_R_rand = h)], on = c("pk_model", "trial")]
  z[, .(analysis_set = k, n_trials = .N, dev_R_median = median(h_R - h_R_rand), dev_R_p05 = quantile(h_R - h_R_rand, 0.05, names = FALSE), dev_R_p95 = quantile(h_R - h_R_rand, 0.95, names = FALSE),
        dev_T_median = median(h - h_T_rand), dev_T_p05 = quantile(h - h_T_rand, 0.05, names = FALSE), dev_T_p95 = quantile(h - h_T_rand, 0.95, names = FALSE),
        armdiff_median = median(h - h_R), armdiff_p05 = quantile(h - h_R, 0.05, names = FALSE), armdiff_p95 = quantile(h - h_R, 0.95, names = FALSE),
        armdiff_mean = mean(h - h_R), armdiff_abs_gt5_pct = 100 * mean(abs(h - h_R) > 5), rand_armdiff_abs_max = max(abs(h_T_rand - h_R_rand))), by = .(pk_model, scenario)]
}))
fwrite(CP, file.path(out_dir, "tp_strata_composition.csv"))

# ---- §2-3 처리군 의존 탈락 ----------------------------------------------------------------------------------------------------
inv <- fread(proj_path("results", "oc", "inversion_all.csv"))[reachable %in% TRUE]; inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
TRU <- rbind(inv[, .(pk_model = model, scenario = code, mechanism, direction, auc_ratio, cmax_ratio)],
             rbindlist(lapply(PK, function(m) { f <- file.path(tp_in, sprintf("trialpop_truth_%s.csv", m)); if (!file.exists(f)) { if (allow_partial) return(NULL) else stop("없음: ", f) }
               fread(f)[, .(pk_model, scenario, mechanism = c(F097 = "F", KE110 = "ke", KE120 = "ke", VM125 = "Vmax")[scenario], direction = c(F097 = "down", KE110 = "up", KE120 = "up", VM125 = "up")[scenario], auc_ratio, cmax_ratio)] })),
             data.table(pk_model = PK, scenario = "S00", mechanism = "-", direction = "-", auc_ratio = 1, cmax_ratio = 1))
tot <- CNT[, lapply(.SD, sum), by = .(pk_model, trial, scenario), .SDcols = c("n_rand", unname(NMAP))]
AD <- rbindlist(lapply(SETS, function(s_) { col <- NMAP[[s_]]
  r <- tot[scenario == "REF", .(pk_model, trial, fR = 100 * (1 - get(col) / n_rand))]; t_ <- tot[scenario != "REF", .(pk_model, trial, scenario, fT = 100 * (1 - get(col) / n_rand))]
  z <- r[t_, on = c("pk_model", "trial")]
  z[, .(set = s_, n_trials = .N, fail_R_mean = mean(fR), fail_T_mean = mean(fT), diff_mean = mean(fT - fR), diff_lo = mean(fT - fR) - 1.96 * sd(fT - fR) / sqrt(.N), diff_hi = mean(fT - fR) + 1.96 * sd(fT - fR) / sqrt(.N),
        diff_p05 = quantile(fT - fR, 0.05, names = FALSE), diff_p95 = quantile(fT - fR, 0.95, names = FALSE)), by = .(pk_model, scenario)] }))
AD <- TRU[AD, on = c("pk_model", "scenario")]
if (anyNA(AD$auc_ratio) && !allow_partial) stop("참값이 없는 시나리오: ", paste(unique(AD[is.na(auc_ratio), scenario]), collapse = ", "))
fwrite(AD, file.path(out_dir, "tp_arm_difference.csv"))
AF <- AD[is.finite(auc_ratio), { fit <- lm(diff_mean ~ log(auc_ratio)); ci <- confint(fit)[2, ]
  .(n_scenarios = .N, slope_pp_per_log = coef(fit)[[2]], slope_lo = ci[[1]], slope_hi = ci[[2]], pp_per_10pct_lower = coef(fit)[[2]] * log(0.9), r_squared = summary(fit)$r.squared) }, by = .(pk_model, set)]
fwrite(AF, file.path(out_dir, "tp_arm_difference_fit.csv"))

# ---- §2-4 미달자 특성 --------------------------------------------------------------------------------------------------------
CH <- rbindlist(lapply(PK, function(m) { x <- IND[[m]]; rbindlist(lapply(SETS, function(s_) { f <- !crit_ok(x, s_)
  w <- welch(x$WT[f], x$WT[!f]); a <- welch(log(x$AUCinf_true[f]), log(x$AUCinf_true[!f])); l <- welch(log(x$AUClast[f & x$AUClast > 0]), log(x$AUClast[!f & x$AUClast > 0]))
  data.table(pk_model = m, set = s_, n_fail = sum(f), n_retained = sum(!f), wt_diff_kg = w[["diff"]], wt_diff_lo = w[["lo"]], wt_diff_hi = w[["hi"]],
             true_aucinf_gmr = exp(a[["diff"]]), true_aucinf_gmr_lo = exp(a[["lo"]]), true_aucinf_gmr_hi = exp(a[["hi"]]), auclast_gmr = exp(l[["diff"]]), auclast_gmr_lo = exp(l[["lo"]]), auclast_gmr_hi = exp(l[["hi"]])) })) }))
fwrite(CH, file.path(out_dir, "tp_characteristics.csv"))

# ---- §3 커버리지 정의 분리 ------------------------------------------------------------------------------------------------------
qs <- function(v, what) { v <- v[is.finite(v)]; data.table(metric = what, n = length(v), median = median(v), p05 = quantile(v, 0.05, names = FALSE), p95 = quantile(v, 0.95, names = FALSE), min = min(v), max = max(v), sd_log = sd(log(v))) }
CV <- rbindlist(lapply(PK, function(m) { x <- IND[[m]]; lz <- x$lambda_ok %in% TRUE; ai <- crit_ok(x, "i")
  rbind(qs(x$coverage_true, "window coverage (true AUC0-tlast / true AUC0-inf)"), qs(x$AUClast / x$AUCinf_true, "observed-to-true, AUClast (all subjects)"),
        qs(x$AUClast[lz] / x$AUCinf_true[lz], "observed-to-true, AUClast (lambda-z estimable)"), qs(x$AUCinf[lz] / x$AUCinf_true[lz], "observed-to-true, AUCinf rule B (lambda-z estimable)"),
        qs(x$AUClast[ai] / x$AUCinf_true[ai], "observed-to-true, AUClast (meeting set (i))"), qs(x$AUCinf[ai] / x$AUCinf_true[ai], "observed-to-true, AUCinf rule A (meeting set (i))"),
        qs(x$AUCinf[lz] / x$AUClast[lz], "extrapolation factor AUCinf / AUClast (lambda-z estimable)"))[, pk_model := m] }))
fwrite(CV, file.path(out_dir, "tp_coverage_individual.csv"))
sc1 <- rbindlist(lapply(PK, function(m) fread(proj_path("results", "oc", sprintf("oc_scenarios_%s.csv", m)))[, .(pk_model = model, scenario = code, target, auc_ratio)]))
sc1 <- sc1[abs(target - 0.80) < 1e-9 | abs(target - 1.25) < 1e-9]
tf97 <- rbindlist(lapply(PK, function(m) { f <- file.path(tp_in, sprintf("trialpop_truth_%s.csv", m)); if (file.exists(f)) fread(f)[scenario == "F097", .(pk_model, scenario, auc_ratio)] else NULL }))
GA <- S1[endpoint == "AUClast" & trial <= NT1 & is.finite(est), .(n_trials = .N, gm_gmr = exp(mean(est))), by = .(pk_model, scenario)]
GA <- rbind(sc1[, .(pk_model, scenario, auc_ratio)], tf97, data.table(pk_model = PK, scenario = "S00", auc_ratio = 1))[GA, on = c("pk_model", "scenario"), nomatch = NULL]
GA[, `:=`(abs_diff = abs(gm_gmr - auc_ratio), rel_diff_pct = 100 * (gm_gmr / auc_ratio - 1))]
fwrite(GA, file.path(out_dir, "tp_gmr_agreement.csv"))

# ---- 그림 --------------------------------------------------------------------------------------------------------------------
fb <- melt(FB[, .(pk_model, set, lz = lz_pct, est = est_fail_pct)], id.vars = c("pk_model", "set"), variable.name = "reason", value.name = "pct")
fb[, `:=`(reason = factor(fifelse(reason == "lz", "lambda-z not estimable", "estimable, criteria not met"), c("estimable, criteria not met", "lambda-z not estimable")), ml = factor(MODEL_EN[pk_model], MODEL_EN))]
g1 <- ggplot(fb, aes(set, pct, fill = reason)) + geom_col(width = 0.65) + geom_errorbar(data = FB[, .(set, ml = factor(MODEL_EN[pk_model], MODEL_EN), lo = fail_lo, hi = fail_hi)], aes(set, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.2, colour = "#52514e") +
  facet_wrap(~ml) + scale_fill_manual(values = c("estimable, criteria not met" = "#2a78d6", "lambda-z not estimable" = "#eb6834"), name = NULL) +
  labs(x = "Criteria set: (i) adj. R-squared 0.80; (ii) (i) + span 2; (iii) adj. R-squared 0.90; (iv) (iii) + span 3 (all with extrapolation at most 20%)", y = "Subjects without a reliable AUC0-inf (%)",
       title = "Study population (healthy, 60 to 90 kg): subjects without a reliable AUC0-inf by criteria set", subtitle = "Planned schedule, 300 mg, 20,000 subjects per model; bars: Wilson 95% interval of the total") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_tp_failure_by_set.png"), g1, width = 9, height = 4.8, dpi = 150, bg = "white")
ad <- AD[set %in% c("i", "iii") & is.finite(auc_ratio)][, `:=`(ml = factor(MODEL_EN[pk_model], MODEL_EN), sl = factor(paste("criteria set", SET_EN[set]), paste("criteria set", SET_EN[c("i", "iii")])),
                                                              grp = fifelse(scenario %in% c("S00", "F097", "KE110", "KE120", "VM125"), "product scenarios", "boundary scenarios"))]
g2 <- ggplot(ad, aes(auc_ratio, diff_mean, colour = grp)) + geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.3) + geom_vline(xintercept = 1, colour = "#b5b3ad", linewidth = 0.3, linetype = "dashed") +
  geom_pointrange(aes(ymin = diff_lo, ymax = diff_hi), size = 0.25) + geom_text(aes(label = scenario), size = 2.2, vjust = -0.9, show.legend = FALSE) +
  facet_grid(sl ~ ml) + scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25)) + scale_colour_manual(values = c("product scenarios" = "#2a78d6", "boundary scenarios" = "#c4502a"), name = NULL) +
  labs(x = "True AUC0-inf ratio, test to reference (log scale)", y = "Test minus reference: subjects failing (percentage points)",
       title = "Study population: the share failing the criteria differs between arms when the products differ", subtitle = sprintf("117 subjects per arm, %s trials per scenario; bars: 95%% interval of the mean", format(NTR, big.mark = ","))) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_tp_arm_difference.png"), g2, width = 10, height = 6.5, dpi = 150, bg = "white")
si <- melt(SI[, .(pk_model, set, light = fail_light_pct, heavy = fail_heavy_pct, llo = fail_light_lo, lhi = fail_light_hi, hlo = fail_heavy_lo, hhi = fail_heavy_hi)], id.vars = c("pk_model", "set"), measure.vars = list(c("light", "heavy"), c("llo", "hlo"), c("lhi", "hhi")),
           variable.name = "stratum", value.name = c("pct", "lo", "hi"))
si[, `:=`(stratum = factor(c("60 to 75 kg", "above 75 to 90 kg")[as.integer(stratum)], c("60 to 75 kg", "above 75 to 90 kg")), ml = factor(MODEL_EN[pk_model], MODEL_EN))]
g3 <- ggplot(si, aes(set, pct, fill = stratum)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) + geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.75), width = 0.2, colour = "#52514e") +
  facet_wrap(~ml) + scale_fill_manual(values = c("60 to 75 kg" = "#8fb8ea", "above 75 to 90 kg" = "#1c5fb0"), name = "Randomization stratum") +
  labs(x = "Criteria set", y = "Subjects without a reliable AUC0-inf (%)", title = "Study population: subjects without a reliable AUC0-inf by randomization stratum",
       subtitle = "Planned schedule, 300 mg, 20,000 subjects per model; bars: Wilson 95% intervals") + theme_minimal(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_tp_strata.png"), g3, width = 9, height = 4.5, dpi = 150, bg = "white")

# ---- 결론 문안 -----------------------------------------------------------------------------------------------------------------
fbr <- function(s_, col = "fail_pct", f = f1, u = "%") rg(FB[set == s_][[col]], f, u)
rpr <- function(s_) { x <- RP[set == s_]; sprintf("median %s (5th to 95th percentile %s to %s)", rg(x$retained_median, function(v) formatC(v, format = "f", digits = 0), ""), min(x$retained_p05), max(x$retained_p95)) }
sir <- function(s_, col, f = f1) rg(SI[set == s_][[col]], f, "")
adr <- function(sc, s_, col = "diff_mean") rg(AD[scenario == sc & set == s_][[col]], f2, "")
en <- c("# Trial population (healthy adults, 60 to 90 kg, weight-stratified randomization): why AUC0-inf is not a reliable primary endpoint", "",
  "## Failure by criteria set (planned schedule, 300 mg, 20,000 subjects per model; two models)", "",
  sprintf("- Set (i) adjusted R-squared at least 0.80: %s fail (lambda-z not estimable %s; estimable but failing %s). Set (ii) with span ratio at least 2: %s. Set (iii) adjusted R-squared at least 0.90: %s. Set (iv) with span at least 3: %s.",
          fbr("i"), fbr("i", "lz_pct", f2), fbr("i", "est_fail_pct"), fbr("ii"), fbr("iii"), fbr("iv")),
  sprintf("- Retained subjects per arm of 117 under rule A in the trials (identical products, %s trials per model): set (i) %s; set (iii) %s; set (iv) %s.", format(NT1, big.mark = ","), rpr("i"), rpr("iii"), rpr("iv")),
  sprintf("- Residual sensitivity: with the 2016 model's proportional residual of %s%% the set (iii) failure share is %s; with the same model and a residual of %s%% it is %s; Model 1 (%s%%) gives %s. A larger proportional residual scatters the terminal log concentrations around the regression line and lowers the adjusted R-squared, so a stricter threshold removes more subjects.",
          f1(100 * sig[["k2016"]]), sprintf("%s%%", f1(RSD[variant == "k2016" & set == "iii", fail_pct])), f1(100 * sig[["resid12"]]), sprintf("%s%%", f1(RSD[variant == "resid12" & set == "iii", fail_pct])), f1(100 * sig[["k2020"]]), sprintf("%s%%", f1(RSD[variant == "k2020" & set == "iii", fail_pct]))), "",
  "## Failure by randomization stratum", "",
  sprintf("- Heavier stratum (above 75 to 90 kg) minus lighter stratum (60 to 75 kg), percentage points: set (i) %s; set (iii) %s (Newcombe 95%% intervals in the table).", sir("i", "diff_pp", f2), sir("iii", "diff_pp", f2)),
  { x <- CP[scenario == "S00" & analysis_set %in% c("auclast", "i", "iii")]; if (nrow(x)) sprintf("- Between-arm difference in the share of the heavier stratum, identical products (median; 5th to 95th percentile over trials, percentage points): AUClast analysis set %s; rule A set (i) %s; rule A set (iii) %s.",
      paste(sprintf("%s (%s to %s)", f1(x[analysis_set == "auclast", armdiff_median]), f1(x[analysis_set == "auclast", armdiff_p05]), f1(x[analysis_set == "auclast", armdiff_p95])), collapse = "; "),
      paste(sprintf("%s (%s to %s)", f1(x[analysis_set == "i", armdiff_median]), f1(x[analysis_set == "i", armdiff_p05]), f1(x[analysis_set == "i", armdiff_p95])), collapse = "; "),
      paste(sprintf("%s (%s to %s)", f1(x[analysis_set == "iii", armdiff_median]), f1(x[analysis_set == "iii", armdiff_p05]), f1(x[analysis_set == "iii", armdiff_p95])), collapse = "; ")) else "" }, "",
  "## Treatment-dependent failure (test minus reference, percentage points, mean over trials)", "",
  if (nrow(AD)) sprintf("- Set (i): identical products %s; F x0.97 %s; ke x1.10 %s; ke x1.20 %s; Vmax x1.25 %s. Set (iii): %s; %s; %s; %s; %s.", adr("S00", "i"), adr("F097", "i"), adr("KE110", "i"), adr("KE120", "i"), adr("VM125", "i"),
                        adr("S00", "iii"), adr("F097", "iii"), adr("KE110", "iii"), adr("KE120", "iii"), adr("VM125", "iii")) else "",
  if (nrow(AF)) sprintf("- Across the scenarios, each 10%% lower true AUC0-inf ratio goes with %s percentage points more test-arm failures under set (i) and %s under set (iii) (descriptive linear fit; R-squared %s and %s).",
                        rg(AF[set == "i", pp_per_10pct_lower], f2, ""), rg(AF[set == "iii", pp_per_10pct_lower], f2, ""), rg(AF[set == "i", r_squared], f2, ""), rg(AF[set == "iii", r_squared], f2, "")) else "", "",
  "## Failing versus retained subjects", "",
  sprintf("- Body weight difference, failing minus retained: set (i) %s kg, set (iv) %s kg; true AUC0-inf geometric mean ratio, failing to retained: set (i) %s, set (iii) %s (ranges over the two models).",
          rg(CH[set == "i", wt_diff_kg], f2, ""), rg(CH[set == "iv", wt_diff_kg], f2, ""), rg(CH[set == "i", true_aucinf_gmr], f3, ""), rg(CH[set == "iii", true_aucinf_gmr], f3, "")), "",
  "## Coverage and the observed-to-true ratio", "",
  sprintf("- Window coverage (true AUC0-tlast / true AUC0-inf): median %s, minimum %s.", rg(100 * CV[metric == "window coverage (true AUC0-tlast / true AUC0-inf)", median], f1), rg(100 * CV[metric == "window coverage (true AUC0-tlast / true AUC0-inf)", min], f1)),
  sprintf("- Observed-to-true ratio, subjects with an estimable lambda-z: AUClast median %s, 5th percentile %s, 95th percentile %s, SD of log %s; AUCinf rule B median %s, 5th percentile %s, 95th percentile %s, SD of log %s (ranges over the two models).",
          rg(100 * CV[metric == "observed-to-true, AUClast (lambda-z estimable)", median], f1), rg(100 * CV[metric == "observed-to-true, AUClast (lambda-z estimable)", p05], f1), rg(100 * CV[metric == "observed-to-true, AUClast (lambda-z estimable)", p95], f1),
          rg(CV[metric == "observed-to-true, AUClast (lambda-z estimable)", sd_log], f3, ""), rg(100 * CV[metric == "observed-to-true, AUCinf rule B (lambda-z estimable)", median], f1),
          rg(100 * CV[metric == "observed-to-true, AUCinf rule B (lambda-z estimable)", p05], f1), rg(100 * CV[metric == "observed-to-true, AUCinf rule B (lambda-z estimable)", p95], f1), rg(CV[metric == "observed-to-true, AUCinf rule B (lambda-z estimable)", sd_log], f3, "")),
  if (nrow(GA)) sprintf("- Geometric mean of the trial AUClast GMRs against the true AUC0-inf ratio (%d scenarios, pooled t-test): largest absolute difference %s.", nrow(GA), f3(max(GA$abs_diff))) else "")
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", paste(en, collapse = "\n"))) stop("영문 결론 규칙 위반")
writeLines(en, file.path(out_dir, "tp_conclusion_en.md"))
ko <- c("# 시험 모집단(건강인 60–90 kg, 체중 층화) 안의 AUC0-inf 탈락", "",
  sprintf("- 미달(두 모델): (i) %s, (ii) %s, (iii) %s, (iv) %s. arm당 잔류(규칙 A, 동일 제품): (i) %s, (iii) %s.", fbr("i"), fbr("ii"), fbr("iii"), fbr("iv"), rpr("i"), rpr("iii")),
  sprintf("- 층 차이(무거운 − 가벼운, %%p): (i) %s, (iii) %s.", sir("i", "diff_pp", f2), sir("iii", "diff_pp", f2)),
  if (nrow(AD)) sprintf("- 처리군 의존(시험군 − 대조군, %%p, 세트 (i)): S00 %s, F097 %s, KE110 %s, KE120 %s, VM125 %s.", adr("S00", "i"), adr("F097", "i"), adr("KE110", "i"), adr("KE120", "i"), adr("VM125", "i")) else "",
  sprintf("- 미달자: 체중 차이 (i) %s kg, 참 AUC0-inf 기하평균비 (i) %s, (iii) %s.", rg(CH[set == "i", wt_diff_kg], f2, ""), rg(CH[set == "i", true_aucinf_gmr], f3, ""), rg(CH[set == "iii", true_aucinf_gmr], f3, "")))
writeLines(ko, file.path(out_dir, "tp_conclusion_ko.md"))
cat(en, sep = "\n")
