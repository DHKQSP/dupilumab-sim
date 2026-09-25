#!/usr/bin/env Rscript
# 45_oc_models_summary.R — 분석 모형 재판정 요약 (추가 지시 2026-09-26 §1, config/prereg_20260926.yaml section1, D-053).
# 입력: results/oc_models/oc_models_be_<model>.csv.gz(시험 1–10,000), oc_models_ext_be_<model>.csv.gz(연장 칸, 있으면),
#       results/oc/oc_scenarios_<model>.csv(배율·참값 비), results/oc/oc_rejudge_be_<model>.csv.gz·oc_trials_ext_be_<model>.csv.gz(M0 재확인),
#       results/oc/truth_ref_<model>.rds(F097 참값 계산용 200,000명 대조 참값).
# 산출(results/oc_models/):
#   m0_reverification.csv        재생성 M0 대 저장본(쓰기 때 검사와 별도로 전체 파일을 다시 대조)
#   truth_F097.csv               F097 참 AUC0-inf·Cmax 비(200,000명 공통 난수, 역산과 같은 방법)
#   type1_models.csv             경계 16칸 × 분석 모형 × 구성: 통과율, Wilson 95% 구간, 분류(최종 시험 수)와 사전 고정 10,000회 값
#   type1_paired_change.csv      같은 시험에서 M1 − M0, M2 − M0 (%p, 쌍대 95% 구간, 불일치 시험 수)
#   classification_change.csv    칸 × 구성별 M0·M1·M2 분류
#   sd_se_models.csv             칸(S00·F097 포함) × 평가변수 × 분석 모형: 편향, 시험 간 SD, 시험 안 SE 중앙값, SD/SE, 이론값·정규 근사 예측
#   power_models.csv             S00·F097 × 분석 모형 × 구성: 검정력(Wilson), M0 대비 쌍대 차이
#   p2_decomposition_models.csv  P2 − 5 = (AUCinf_true 단독 − 5) + (AUClast 단독 − AUCinf_true 단독) + (P2 − AUClast 단독)
#   not_estimable.csv            분석 모형·평가변수별 판정 불가(pass NA) 시험 수
#   expectations_check.csv       사전 등록한 사용자 예상의 대조
#   fig_type1_models.png         P2와 불편 참조(AUCinf_true 단독)의 경계 1종 오류, M0·M1·M2
#   oc_models_conclusion_en.md, oc_models_conclusion_ko.md
source("R/00_setup.R"); source_project()
pr <- read_cfg("prereg_20260926.yaml")$section1
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
out_dir <- Sys.getenv("DUPI_OC_MODELS_DIR", proj_path(pr$out_dir)); oc_dir <- proj_path("results", "oc")   # 환경변수는 시험용(작은 산출 디렉터리)
CFG <- lapply(pr$configuration_endpoints, unlist)
stopifnot(identical(names(CFG), unlist(pr$configurations)), all(unlist(CFG) %in% OC_ENDPOINTS_EXT))
for (cf in intersect(names(CFG), names(OC_REJUDGE_CONFIGS))) stopifnot(setequal(CFG[[cf]], OC_REJUDGE_CONFIGS[[cf]]))
for (cf in c("F3A", "F3B", "F3C", "AUClast_only")) stopifnot(setequal(CFG[[cf]], unlist(oc$configurations[[cf]]$endpoints)))
ALPHA <- 5; CL <- design$be$ci_level; LIMS <- as.numeric(unlist(design$be$limits)); Z <- qnorm(0.975)
n_bnd <- as.integer(Sys.getenv("DUPI_OC_MODELS_TEST_N", oc$trials$reps_boundary)); AM <- unlist(pr$analysis_models_run)   # 시험용 축소 시험 수(운영 실행은 설정하지 않음)
PK <- unlist(pr$pk_models)
MODEL_EN <- c(k2016 = "Kovalenko 2016 (primary)", k2020 = "Kovalenko 2020 Model 1"); MODEL_KO <- c(k2016 = "Kovalenko 2016 (주)", k2020 = "Kovalenko 2020 Model 1")
MODEL_EN_S <- c(k2016 = "2016 model", k2020 = "Model 1")
MECH_EN <- c(F = "bioavailability (F)", ka = "absorption rate (ka)", ke = "linear elimination (ke)", Vmax = "target-mediated elimination capacity (Vmax)", Km = "binding constant (Km)", V2 = "peripheral volume (V2)")
MECH_KO <- c(F = "F(흡수량)", ka = "ka(흡수 속도)", ke = "ke(선형 소실)", Vmax = "Vmax(표적 매개 소실)", Km = "Km(결합)", V2 = "V2(말초 분포용적)")
AM_EN <- c(M0 = "M0 (pooled t)", M1 = "M1 (ANOVA, weight stratum)", M2 = "M2 (ANCOVA, log weight)")
CFG_EN <- c(P2 = "P2", G2_Aii = "G2-A(ii)", G2_Ai = "G2-A(i)", G2_B = "G2-B", G2_Cii = "G2-C(ii)", G2_Ci = "G2-C(i)", F3A = "F3-A", F3B = "F3-B", F3C = "F3-C",
            AUClast_only = "AUClast only", AUCinf_true_only = "AUCinf_true only")
CLASS_KO <- c(conservative = "보수적", nominal = "명목", exceeding = "초과")
f2 <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2); s2 <- function(x) sprintf("%+.2f", round(x, 2) + 0)
f3 <- function(x) formatC(round(x, 3) + 0, format = "f", digits = 3); f4 <- function(x) formatC(round(x, 4) + 0, format = "f", digits = 4)
fint <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
rng <- function(x, f = f2, sep = " to ") if (length(x)) paste0(f(min(x)), sep, f(max(x))) else "none"

# ---- 1) 자료 읽기와 완전성 --------------------------------------------------------------------------------------------
sc_meta <- rbindlist(lapply(PK, function(m) fread(file.path(oc_dir, sprintf("oc_scenarios_%s.csv", m)))[, .(pk_model = model, scenario = code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio)]))
bnd <- as.numeric(unlist(oc$boundary_targets))
sc_meta <- sc_meta[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9]
stopifnot(sc_meta[, .N, by = pk_model]$N == 8L)
read_be <- function(m) {
  x <- fread(file.path(out_dir, sprintf("oc_models_be_%s.csv.gz", m)))
  scn <- c("S00", sc_meta[pk_model == m, scenario], names(pr$scenarios$products))
  if (!setequal(unique(x$scenario), scn)) stop(m, ": 시나리오 집합이 사전 등록과 다릅니다")
  cnt <- x[, .N, by = trial]
  if (!identical(sort(cnt$trial), seq_len(n_bnd)) || any(cnt$N != length(scn) * length(OC_ENDPOINTS_EXT) * length(AM))) stop(m, ": main 산출이 시험 1-", n_bnd, " 전체가 아닙니다")
  fe <- file.path(out_dir, sprintf("oc_models_ext_be_%s.csv.gz", m))
  if (file.exists(fe)) {
    e <- fread(fe)
    if (any(e$trial <= n_bnd) || !all(e$scenario %in% sc_meta[pk_model == m, scenario])) stop(m, ": 연장 파일 범위 오류")
    ce <- e[, .(n = uniqueN(trial), mx = max(trial), N = .N), by = scenario]
    if (any(ce$N != ce$n * length(OC_ENDPOINTS_EXT) * length(AM)) || any(ce$mx != n_bnd + ce$n)) stop(m, ": 연장 파일이 연속된 완전한 시험이 아닙니다")
    x <- rbind(x, e)
  }
  if (anyDuplicated(x, by = c("trial", "scenario", "endpoint", "model"))) stop(m, ": 중복 행")
  x[, pk_model := m]
  x[, df := n_R + n_T - fifelse(model == "M0", 2L, 3L)]
  x[]
}
BE <- rbindlist(lapply(PK, read_be))
dec_f <- file.path(out_dir, "extension_decision_models.csv")
DEC <- if (file.exists(dec_f)) fread(dec_f) else NULL
for (m in PK) {                                                        # 연장 판정 파일과 연장 산출이 맞는지
  ext_sc <- unique(BE[pk_model == m & trial > n_bnd, scenario])
  sel <- if (!is.null(DEC) && m %in% DEC$pk_model) unique(DEC[pk_model == m & selected == TRUE, scenario]) else character(0)
  if (!setequal(ext_sc, sel)) stop(sprintf("%s: 연장 산출(%s)과 연장 판정(%s)이 다릅니다", m, paste(ext_sc, collapse = ","), paste(sel, collapse = ",")))
  if (length(sel) && any(BE[pk_model == m & scenario %in% sel, max(trial), by = scenario]$V1 != as.integer(pr$extension$to))) stop(m, ": 연장 칸이 목표 시험 수에 도달하지 않았습니다")
}

# ---- 2) M0 재확인(전체 파일) --------------------------------------------------------------------------------------------
rel <- function(a, b) ifelse(is.na(a) & is.na(b), 0, abs(a / b - 1))
reverify <- function(m) {
  x <- BE[pk_model == m & model == "M0"]
  q <- qt(1 - (1 - CL) / 2, x$df)
  x[, `:=`(GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se))]
  out <- list()
  st <- fread(file.path(oc_dir, sprintf("oc_rejudge_be_%s.csv.gz", m)))[trial %in% unique(x$trial)]
  mm <- merge(x[trial <= n_bnd & scenario != "F097"], st, by = c("trial", "scenario", "endpoint"), suffixes = c("", ".s"))
  out[[1]] <- data.table(pk_model = m, stored_file = sprintf("results/oc/oc_rejudge_be_%s.csv.gz", m), n_rows = nrow(mm), n_expected = nrow(st),
                         max_rel_diff = max(rel(mm$GMR, mm$GMR.s), rel(mm$CI_lower, mm$CI_lower.s), rel(mm$CI_upper, mm$CI_upper.s), na.rm = TRUE),
                         n_pass_differ = sum(!((is.na(mm$pass) & is.na(mm$pass.s)) | (mm$pass %in% TRUE & mm$pass.s %in% TRUE) | (mm$pass %in% FALSE & mm$pass.s %in% FALSE))),
                         n_n_differ = sum(mm$n_R != mm$n_R.s | mm$n_T != mm$n_T.s), tolerance = 1e-6, note = "est and se stored at 8 significant digits, stored file at 7")
  fe <- file.path(oc_dir, sprintf("oc_trials_ext_be_%s.csv.gz", m))
  if (file.exists(fe) && any(x$trial > n_bnd)) {
    se_ <- fread(fe)[trial %in% unique(x$trial)]; mm <- merge(x[trial > n_bnd], se_, by = c("trial", "scenario", "endpoint"), suffixes = c("", ".s"))
    out[[2]] <- data.table(pk_model = m, stored_file = sprintf("results/oc/oc_trials_ext_be_%s.csv.gz", m), n_rows = nrow(mm), n_expected = nrow(se_[scenario %in% unique(x[trial > n_bnd, scenario])]),
                           max_rel_diff = max(rel(mm$GMR, mm$GMR.s), rel(mm$CI_lower, mm$CI_lower.s), rel(mm$CI_upper, mm$CI_upper.s), na.rm = TRUE),
                           n_pass_differ = sum(!((is.na(mm$pass) & is.na(mm$pass.s)) | (mm$pass %in% TRUE & mm$pass.s %in% TRUE) | (mm$pass %in% FALSE & mm$pass.s %in% FALSE))),
                           n_n_differ = sum(mm$n_R != mm$n_R.s | mm$n_T != mm$n_T.s), tolerance = 1e-6, note = "extension trials")
  }
  PM <- c(Cmax = "Cmax", AUClast = "AUClast", AUCinf_all = "AUCinf_B", AUCinf_reliable = "AUCinf_A", AUCinf_true = "AUCinf_true", AUCinf_subC = "AUCinf_C")
  pf <- if (m == "k2016") "results/trials5000/products5000_be_raw_base.csv.gz" else "results/trials/products_be_raw_struct2020.csv.gz"
  pp <- fread(proj_path(pf))[scenario == "F097" & method == "pooled_t" & schedule == "B0" & trial %in% unique(x$trial)]; pp[, endpoint := unname(PM[endpoint])]
  mm <- merge(x[scenario == "F097"], pp[, .(trial, scenario, endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T)], by = c("trial", "scenario", "endpoint"), suffixes = c("", ".s"))
  out[[3]] <- data.table(pk_model = m, stored_file = pf, n_rows = nrow(mm), n_expected = nrow(pp),
                         max_rel_diff = max(rel(mm$GMR, mm$GMR.s), rel(mm$CI_lower, mm$CI_lower.s), rel(mm$CI_upper, mm$CI_upper.s), na.rm = TRUE),
                         n_pass_differ = sum(!((is.na(mm$pass) & is.na(mm$pass.s)) | (mm$pass %in% TRUE & mm$pass.s %in% TRUE) | (mm$pass %in% FALSE & mm$pass.s %in% FALSE))),
                         n_n_differ = sum(mm$n_R != mm$n_R.s | mm$n_T != mm$n_T.s), tolerance = 1e-6, note = sprintf("F097, trials %d-%d, 6 original endpoints", min(pp$trial), max(pp$trial)))
  rbindlist(out)
}
RV <- rbindlist(lapply(PK, reverify))
RV[, ok := n_rows == n_expected & max_rel_diff <= tolerance & n_pass_differ == 0 & n_n_differ == 0]
fwrite(RV, file.path(out_dir, "m0_reverification.csv"))
if (!all(RV$ok)) { print(RV); stop("M0 재확인 실패") }

# ---- 3) F097 참값(200,000명, 역산과 같은 방법) ---------------------------------------------------------------------------
tf <- file.path(out_dir, "truth_F097.csv")
if (!file.exists(tf)) {
  TR <- rbindlist(lapply(PK, function(m) {
    p <- load_params(m); N <- as.integer(oc$estimand$population$n_subjects)
    subj <- truth_subjects(N, p, design, derive_seed(as.integer(oc$truth_seed), "truth200k", m))
    ref <- readRDS(file.path(oc_dir, sprintf("truth_ref_%s.rds", m)))
    te <- truth_metrics(ip_for_mult(p, subj, list(F = 0.97)), as.numeric(oc$estimand$dose_mg), p$model_id, cmax = TRUE)
    r <- truth_ratio(ref, te)
    data.table(pk_model = m, scenario = "F097", multiplier_F = 0.97, auc_ratio = r$auc_ratio, auc_se_log = r$auc_se_log, cmax_ratio = r$cmax_ratio, cmax_se_log = r$cmax_se_log, n_subjects = r$n)
  }))
  fwrite(TR, tf)
}
TRF <- fread(tf)
truth_tab <- rbind(sc_meta[, .(pk_model, scenario, auc_ratio, cmax_ratio)], TRF[, .(pk_model, scenario, auc_ratio, cmax_ratio)],
                   data.table(pk_model = PK, scenario = "S00", auc_ratio = 1, cmax_ratio = 1))

# ---- 4) 구성 통과 ------------------------------------------------------------------------------------------------------
W <- dcast(BE[, .(pk_model, trial, scenario, model, endpoint, ok = pass %in% TRUE)], pk_model + trial + scenario + model ~ endpoint, value.var = "ok")
for (cf in names(CFG)) W[, (cf) := Reduce(`&`, .SD), .SDcols = CFG[[cf]]]
NE <- BE[, .(n_trials = uniqueN(trial), n_not_estimable = sum(is.na(pass))), by = .(pk_model, model, endpoint)]
fwrite(NE, file.path(out_dir, "not_estimable.csv"))

rate <- function(x) { w <- wilson_ci(sum(x), length(x)); list(n_trials = length(x), n_pass = sum(x), pass_pct = w$est, lo = w$lo, hi = w$hi) }
long <- melt(W, id.vars = c("pk_model", "trial", "scenario", "model"), measure.vars = names(CFG), variable.name = "config", value.name = "cfg_pass")
long[, config := as.character(config)]
T1 <- long[scenario %in% sc_meta$scenario, rate(cfg_pass), by = .(pk_model, scenario, model, config)]
T1_10k <- long[scenario %in% sc_meta$scenario & trial <= n_bnd, rate(cfg_pass), by = .(pk_model, scenario, model, config)]
setnames(T1_10k, c("n_trials", "n_pass", "pass_pct", "lo", "hi"), paste0(c("n_trials", "n_pass", "pass_pct", "lo", "hi"), "_10k"))
T1 <- merge(T1, T1_10k, by = c("pk_model", "scenario", "model", "config"))
T1 <- sc_meta[T1, on = c("pk_model", "scenario"), nomatch = NULL]
T1[, class := classify_type1(lo, hi, ALPHA)][, class_10k := classify_type1(lo_10k, hi_10k, ALPHA)]
setnames(T1, "model", "analysis_model")
setcolorder(T1, c("pk_model", "scenario", "mechanism", "direction", "target", "multiplier", "auc_ratio", "cmax_ratio", "analysis_model", "config"))
setorder(T1, pk_model, mechanism, direction, analysis_model, config)
fwrite(T1, file.path(out_dir, "type1_models.csv"))

paired <- function(sc_set) {
  a <- dcast(long[scenario %in% sc_set], pk_model + scenario + config + trial ~ model, value.var = "cfg_pass")
  rbindlist(lapply(setdiff(AM, "M0"), function(mm) a[, { d <- paired_prop_diff_ci(get(mm), M0)
    .(comparison = paste(mm, "minus M0"), diff_pp = d$est, lo = d$lo, hi = d$hi, n_trials = d$n, n_gain = sum(get(mm) & !M0), n_loss = sum(!get(mm) & M0)) }, by = .(pk_model, scenario, config)]))
}
PC <- paired(sc_meta$scenario)
fwrite(PC, file.path(out_dir, "type1_paired_change.csv"))
CC <- dcast(T1[, .(pk_model, scenario, config, analysis_model, class)], pk_model + scenario + config ~ analysis_model, value.var = "class")
CC[, changed_M1 := M1 != M0][, changed_M2 := M2 != M0]
fwrite(CC, file.path(out_dir, "classification_change.csv"))

# ---- 5) 편향, 시험 간 SD / 시험 안 SE -----------------------------------------------------------------------------------
SS <- BE[is.finite(est) & is.finite(se), .(n_trials = .N, mean_log_gmr = mean(est), sd_log_gmr = sd(est), se_median = median(se), df_median = median(df)),
         by = .(pk_model, scenario, endpoint, analysis_model = model)]
SS <- truth_tab[SS, on = c("pk_model", "scenario")]
SS[, true_ratio := fifelse(endpoint == "Cmax", cmax_ratio, auc_ratio)][, c("auc_ratio", "cmax_ratio") := NULL]
SS[, `:=`(bias_log = mean_log_gmr - log(true_ratio), bias_mc_se = sd_log_gmr / sqrt(n_trials), sd_se_ratio = sd_log_gmr / se_median)]
SS[, `:=`(bias_lo = bias_log - Z * bias_mc_se, bias_hi = bias_log + Z * bias_mc_se, bias_pct = 100 * (exp(bias_log) - 1))]
SS[, sd_se_lo := sd_se_ratio * sqrt(qchisq(0.025, n_trials - 1) / (n_trials - 1))][, sd_se_hi := sd_se_ratio * sqrt(qchisq(0.975, n_trials - 1) / (n_trials - 1))]   # SD의 카이제곱 구간(se 중앙값 고정)
SS[, pred_pct := 100 * be_pass_prob_normal(mean_log_gmr, sd_log_gmr, se_median, df_median, CL, LIMS)]
SS[scenario %in% sc_meta$scenario, theory_boundary_pct := 100 * be_boundary_type1_theory(se_median, df_median, true_ratio, CL, LIMS)$p_pass]
SS[scenario %in% sc_meta$scenario, bias_dir := bias_direction(bias_log, true_ratio, bias_lo, bias_hi)]
setorder(SS, pk_model, scenario, endpoint, analysis_model)
fwrite(SS, file.path(out_dir, "sd_se_models.csv"))

# ---- 6) 검정력 ---------------------------------------------------------------------------------------------------------
PWS <- c("S00", names(pr$scenarios$products))
PW <- long[scenario %in% PWS, rate(cfg_pass), by = .(pk_model, scenario, analysis_model = model, config)]
PWd <- paired(PWS)[, .(pk_model, scenario, config, analysis_model = sub(" minus M0", "", comparison), diff_vs_M0_pp = diff_pp, diff_lo = lo, diff_hi = hi)]
PW <- merge(PW, PWd, by = c("pk_model", "scenario", "config", "analysis_model"), all.x = TRUE)
PW <- truth_tab[PW, on = c("pk_model", "scenario")]
setorder(PW, pk_model, scenario, config, analysis_model)
fwrite(PW, file.path(out_dir, "power_models.csv"))

# ---- 7) P2 분해 --------------------------------------------------------------------------------------------------------
g <- function(cf) T1[config == cf, .(pk_model, scenario, analysis_model, v = pass_pct)]
DE <- Reduce(function(a, b) merge(a, b, by = c("pk_model", "scenario", "analysis_model")),
             list(setnames(g("P2"), "v", "p2_pct"), setnames(g("AUClast_only"), "v", "auclast_pct"), setnames(g("AUCinf_true_only"), "v", "ref_pct")))
DE[, `:=`(ref_minus_5_pp = ref_pct - ALPHA, auclast_minus_ref_pp = auclast_pct - ref_pct, p2_minus_auclast_pp = p2_pct - auclast_pct, p2_minus_5_pp = p2_pct - ALPHA)]
stopifnot(all(abs(DE$ref_minus_5_pp + DE$auclast_minus_ref_pp + DE$p2_minus_auclast_pp - DE$p2_minus_5_pp) < 1e-9))
DE <- merge(DE, SS[endpoint == "AUClast", .(pk_model, scenario, analysis_model, auclast_bias_pct = bias_pct, auclast_bias_lo_pct = 100 * (exp(bias_lo) - 1),
                                              auclast_bias_hi_pct = 100 * (exp(bias_hi) - 1), auclast_bias_dir = bias_dir)], by = c("pk_model", "scenario", "analysis_model"))
DE <- merge(DE, SS[endpoint == "AUCinf_true", .(pk_model, scenario, analysis_model, ref_sd_se = sd_se_ratio)], by = c("pk_model", "scenario", "analysis_model"))
DE <- merge(DE, T1[config == "P2", .(pk_model, scenario, analysis_model, n_trials, lo, hi, class)], by = c("pk_model", "scenario", "analysis_model"))
DE <- sc_meta[DE, on = c("pk_model", "scenario")]
setorder(DE, pk_model, mechanism, direction, analysis_model)
fwrite(DE, file.path(out_dir, "p2_decomposition_models.csv"))

# ---- 8) 사용자 예상 대조(사전 등록 section1$expectation_checks의 기준) ----------------------------------------------------
EC <- pr$expectation_checks
ref1 <- T1[config == "AUCinf_true_only" & analysis_model == "M1"]; ss1 <- SS[endpoint == "AUCinf_true" & analysis_model == "M1" & scenario %in% sc_meta$scenario]
ss0 <- SS[endpoint == "AUCinf_true" & analysis_model == "M0" & scenario %in% sc_meta$scenario]
pool1 <- wilson_ci(sum(ref1$n_pass), sum(ref1$n_trials))
e1a_ok <- pool1$lo <= ALPHA && pool1$hi >= ALPHA && all(ref1$pass_pct >= EC$e1_cell_range_pct[[1]] & ref1$pass_pct <= EC$e1_cell_range_pct[[2]])
e1b_ok <- all(ss1$sd_se_ratio >= EC$e1_sd_se_range[[1]] & ss1$sd_se_ratio <= EC$e1_sd_se_range[[2]])
pc2 <- PC[config == "P2" & comparison == "M1 minus M0"]
e2_ok <- sum(pc2$diff_pp > 0) >= EC$e2_min_cells_up && mean(pc2$diff_pp) > 0 && mean(pc2$diff_pp) <= EC$e2_max_mean_rise_pp
v2 <- T1[pk_model == "k2020" & mechanism == "V2" & config == "P2" & analysis_model == "M1"]
e3_ok <- abs(v2$pass_pct - EC$e3_expected_pct) <= EC$e3_tolerance_pp
EX <- data.table(
  expectation = c("M1: unbiased reference (AUCinf_true only) close to 5%", "M1: between-trial SD / within-trial SE close to 1", "P2 rises slightly under M1",
                  "Model 1 V2 x2.28 cell about 5.6% under M1 (exceeding)"),
  criterion = c(sprintf("pooled 16-cell Wilson 95%% interval includes 5%% and every cell within %s-%s%%", EC$e1_cell_range_pct[[1]], EC$e1_cell_range_pct[[2]]),
                sprintf("every cell within %s-%s (AUCinf_true)", EC$e1_sd_se_range[[1]], EC$e1_sd_se_range[[2]]),
                sprintf("M1 minus M0 > 0 in at least %d of 16 cells and mean rise above 0 and at most %s points", EC$e2_min_cells_up, EC$e2_max_mean_rise_pp),
                sprintf("point estimate within %s +/- %s%%", EC$e3_expected_pct, EC$e3_tolerance_pp)),
  observed = c(sprintf("pooled %s%% [%s, %s]; cells %s%% (M0 %s%%)", f2(pool1$est), f2(pool1$lo), f2(pool1$hi), rng(ref1$pass_pct), rng(T1[config == "AUCinf_true_only" & analysis_model == "M0", pass_pct])),
               sprintf("M1 %s (M0 %s)", rng(ss1$sd_se_ratio, f3), rng(ss0$sd_se_ratio, f3)),
               sprintf("M1 minus M0 > 0 in %d of %d cells; mean %s points (range %s)", sum(pc2$diff_pp > 0), nrow(pc2), s2(mean(pc2$diff_pp)), rng(pc2$diff_pp, s2)),
               sprintf("%s%% [%s, %s], %s, %s trials (M0 %s%% [%s, %s], %s)", f2(v2$pass_pct), f2(v2$lo), f2(v2$hi), v2$class, fint(v2$n_trials),
                       f2(T1[pk_model == "k2020" & mechanism == "V2" & config == "P2" & analysis_model == "M0", pass_pct]),
                       f2(T1[pk_model == "k2020" & mechanism == "V2" & config == "P2" & analysis_model == "M0", lo]),
                       f2(T1[pk_model == "k2020" & mechanism == "V2" & config == "P2" & analysis_model == "M0", hi]),
                       T1[pk_model == "k2020" & mechanism == "V2" & config == "P2" & analysis_model == "M0", class])),
  consistent = c(e1a_ok, e1b_ok, e2_ok, e3_ok))
fwrite(EX, file.path(out_dir, "expectations_check.csv"))

# ---- 9) 그림: P2와 불편 참조, 분석 모형별 ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
COLS <- c(M0 = "#2a78d6", M1 = "#eb6834", M2 = "#1baf7a")           # dataviz 기준 팔레트 1-3번(3개는 전쌍 CVD 검증 통과)
fd <- T1[config %in% c("P2", "AUCinf_true_only")]
fd[, cell := sprintf("%s %s %s x%.3g", c(k2016 = "2016", k2020 = "2020")[pk_model], mechanism, direction, multiplier)]   # "M1"은 분석 모형과 헷갈리므로 연도로
fd[, panel := factor(fifelse(config == "P2", "P2 (AUClast + Cmax)", "AUCinf_true only (unbiased reference)"), c("P2 (AUClast + Cmax)", "AUCinf_true only (unbiased reference)"))]
fd[, cell := factor(cell, rev(unique(fd[order(pk_model, mechanism, direction), cell])))]
pl <- ggplot(fd, aes(x = pass_pct, y = cell, colour = analysis_model)) +
  geom_vline(xintercept = ALPHA, colour = "#52514e", linewidth = 0.4, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  facet_wrap(~panel, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = COLS, labels = AM_EN, name = "Analysis model") +
  labs(x = "Boundary type I error, % (Wilson 95% interval)", y = NULL,
       title = "Boundary type I error by analysis model", subtitle = "16 boundary cells (true AUC0-inf ratio 0.80 or 1.25; 2016 = primary PK model, 2020 = Kovalenko 2020 Model 1); dashed line 5%") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "fig_type1_models.png"), pl, width = 9, height = 6.5, dpi = 150, bg = "white")

# ---- 10) 결론 문안 -----------------------------------------------------------------------------------------------------
cell_en <- function(r) sprintf("%s, %s %s (multiplier x%s, true AUC0-inf ratio %s)", MODEL_EN[[r$pk_model]], MECH_EN[[r$mechanism]], r$direction, formatC(r$multiplier, format = "fg", digits = 3), f4(r$auc_ratio))
cell_ko <- function(r) sprintf("%s · %s %s(배율 ×%s, 참 AUC0-inf 비 %s)", MODEL_KO[[r$pk_model]], MECH_KO[[r$mechanism]], c(up = "↑", down = "↓")[[r$direction]], formatC(r$multiplier, format = "fg", digits = 3), f4(r$auc_ratio))
over <- DE[p2_pct > ALPHA]
p2c <- function(m, a, cf = "P2") T1[pk_model == m & analysis_model == a & config == cf]
first_en <- if (nrow(over)) {
  cells <- unique(over[, .(pk_model, scenario)])
  paste(vapply(seq_len(nrow(cells)), function(i) {
    rr <- DE[pk_model == cells$pk_model[i] & scenario == cells$scenario[i]]; r1 <- rr[1]
    by_m <- paste(vapply(AM, function(a) { x <- rr[analysis_model == a]; sprintf("%s %s%% [%s, %s], %s", a, f2(x$p2_pct), f2(x$lo), f2(x$hi), x$class) }, ""), collapse = "; ")
    ov <- rr[p2_pct > ALPHA]
    cause <- paste(vapply(seq_len(nrow(ov)), function(k) { x <- ov[k]
      sprintf("under %s the AUClast GMR is biased %s by %s%% (95%% interval %s to %s) relative to the true AUC0-inf ratio, which moves the pass rate from the unbiased reference %s%% (%s points versus 5%%) to %s%% for AUClast alone (%s points); Cmax failures change it by %s points to the P2 value",
              x$analysis_model, if (identical(x$auclast_bias_dir, "toward_1")) "toward 1" else if (identical(x$auclast_bias_dir, "away_from_1")) "away from 1" else "with an interval including 0",
              s2(x$auclast_bias_pct), s2(x$auclast_bias_lo_pct), s2(x$auclast_bias_hi_pct), f2(x$ref_pct), s2(x$ref_minus_5_pp), f2(x$auclast_pct), s2(x$auclast_minus_ref_pp), s2(x$p2_minus_auclast_pp)) }, ""), collapse = "; ")
    plaus <- if (r1$mechanism == "V2") sprintf("The cell corresponds to a proposed biosimilar with the same amino acid sequence as the reference IgG4 antibody whose peripheral volume of distribution is %s times that of the reference, all other parameters being equal (in the 2016 model the same true ratio requires %s times).",
                                               formatC(r1$multiplier, format = "fg", digits = 3), formatC(sc_meta[pk_model == "k2016" & mechanism == "V2" & direction == r1$direction & abs(target - r1$target) < 1e-9, multiplier], format = "fg", digits = 3))
             else sprintf("The cell corresponds to a proposed biosimilar with the same amino acid sequence as the reference IgG4 antibody whose %s is %s times that of the reference, all other parameters being equal.", MECH_EN[[r1$mechanism]], formatC(r1$multiplier, format = "fg", digits = 3))
    sprintf("%s: P2 boundary type I error %s. Quantitative cause: %s. %s", cell_en(r1), by_m, cause, plaus) }, ""), collapse = " ")
} else "No P2 boundary cell has a point estimate above 5% under any analysis model."
first_ko <- if (nrow(over)) {
  cells <- unique(over[, .(pk_model, scenario)])
  paste(vapply(seq_len(nrow(cells)), function(i) {
    rr <- DE[pk_model == cells$pk_model[i] & scenario == cells$scenario[i]]; r1 <- rr[1]
    by_m <- paste(vapply(AM, function(a) { x <- rr[analysis_model == a]; sprintf("%s %s%% [%s, %s] %s", a, f2(x$p2_pct), f2(x$lo), f2(x$hi), CLASS_KO[[x$class]]) }, ""), collapse = "; ")
    ov <- rr[p2_pct > ALPHA]
    cause <- paste(vapply(seq_len(nrow(ov)), function(k) { x <- ov[k]
      sprintf("%s에서 AUClast GMR이 참 AUC0-inf 비 대비 %s%%(95%% 구간 %s ~ %s) %s 편향되어, 불편 참조 %s%%(5%% 대비 %s%%p)에서 AUClast 단독 %s%%(%s%%p)로 올라가고, Cmax 불통과가 %s%%p를 더해 P2가 된다",
              x$analysis_model, s2(x$auclast_bias_pct), s2(x$auclast_bias_lo_pct), s2(x$auclast_bias_hi_pct),
              if (identical(x$auclast_bias_dir, "toward_1")) "1 쪽으로" else if (identical(x$auclast_bias_dir, "away_from_1")) "1에서 멀어지는 쪽으로" else "(0 포함)",
              f2(x$ref_pct), s2(x$ref_minus_5_pp), f2(x$auclast_pct), s2(x$auclast_minus_ref_pp), s2(x$p2_minus_auclast_pp)) }, ""), collapse = "; ")
    plaus <- if (r1$mechanism == "V2") sprintf("이 칸은 대조약(IgG4)과 아미노산 서열이 같은 시험약의 말초 분포용적이 대조약의 %s배이고 나머지 파라미터가 같은 경우다(2016 모델에서는 같은 참값 비에 %s배가 필요).",
                                               formatC(r1$multiplier, format = "fg", digits = 3), formatC(sc_meta[pk_model == "k2016" & mechanism == "V2" & direction == r1$direction & abs(target - r1$target) < 1e-9, multiplier], format = "fg", digits = 3))
             else sprintf("이 칸은 대조약(IgG4)과 아미노산 서열이 같은 시험약의 %s가 대조약의 %s배이고 나머지 파라미터가 같은 경우다.", MECH_KO[[r1$mechanism]], formatC(r1$multiplier, format = "fg", digits = 3))
    sprintf("%s: P2 경계 1종 오류 %s. 정량 원인: %s. %s", cell_ko(r1), by_m, cause, plaus) }, ""), collapse = " ")
} else "어느 분석 모형에서도 P2 점추정이 5%를 넘는 경계 칸은 없다."

cls_count <- function(a, cf) { x <- T1[analysis_model == a & config == cf]; k <- table(factor(x$class, OC_TYPE1_CLASSES)); k }
cls_txt_en <- function(a, cf) { k <- cls_count(a, cf); sprintf("conservative %d, nominal %d, exceeding %d", k[["conservative"]], k[["nominal"]], k[["exceeding"]]) }
cls_txt_ko <- function(a, cf) { k <- cls_count(a, cf); sprintf("보수적 %d, 명목 %d, 초과 %d", k[["conservative"]], k[["nominal"]], k[["exceeding"]]) }
fam <- function(a, cf) { x <- T1[analysis_model == a & config == cf]; list(n_over = sum(x$pass_pct > ALPHA), n_lo = sum(x$lo > ALPHA), mx = max(x$pass_pct), cell = x[which.max(pass_pct)]) }
chg <- CC[config == "P2" & changed_M1 == TRUE]
pwr <- function(m, sc, a, cf = "P2") PW[pk_model == m & scenario == sc & analysis_model == a & config == cf]
rs <- function(a) SS[endpoint == "AUCinf_true" & analysis_model == a & scenario %in% sc_meta$scenario, sd_se_ratio]
ra <- function(a) SS[endpoint == "AUClast" & analysis_model == a & scenario %in% sc_meta$scenario, sd_se_ratio]
ref_rng <- function(a) T1[config == "AUCinf_true_only" & analysis_model == a, pass_pct]
p2_rng <- function(a) T1[config == "P2" & analysis_model == a, pass_pct]
d12 <- PC[config == "P2" & comparison == "M1 minus M0"]; d22 <- PC[config == "P2" & comparison == "M2 minus M0"]

n_zero <- T1[config == "P2" & analysis_model == "M1" & n_pass == 0, .N]
en <- c("# Analysis-model re-judgement of the boundary type I error (M0 pooled t, M1 weight-stratum ANOVA, M2 log-weight ANCOVA)", "",
  sprintf("Generated by `scripts/45_oc_models_summary.R` from `results/oc_models/` (pre-registered in `config/prereg_20260926.yaml`, section1). Trials: %s per cell and model (the extended cells at %s).",
          fint(n_bnd), fint(as.integer(pr$extension$to))), "",
  "## Conclusions", "",
  paste("1.", first_en), "",
  sprintf("2. Classification of P2 over the 16 boundary cells: M0 %s; M1 %s; M2 %s. P2 ranges %s%% (M0), %s%% (M1) and %s%% (M2). On the same trials M1 minus M0 is %s points (mean %s; positive in %d of 16 cells) and M2 minus M0 is %s points (mean %s). Cells whose P2 class changes from M0 to M1: %s.",
          cls_txt_en("M0", "P2"), cls_txt_en("M1", "P2"), cls_txt_en("M2", "P2"), rng(p2_rng("M0")), rng(p2_rng("M1")), rng(p2_rng("M2")),
          rng(d12$diff_pp, s2), s2(mean(d12$diff_pp)), sum(d12$diff_pp > 0), rng(d22$diff_pp, s2), s2(mean(d22$diff_pp)),
          if (nrow(chg)) paste(sprintf("%s %s (%s to %s)", MODEL_EN_S[chg$pk_model], chg$scenario, chg$M0, chg$M1), collapse = "; ") else "none"), "",
  sprintf("3. Unbiased reference (AUCinf_true only, the true AUC0-inf of each subject analysed like an endpoint): %s%% under M0, %s%% under M1 and %s%% under M2. The between-trial SD of the log GMR divided by the median within-trial SE is %s (M0), %s (M1) and %s (M2) for AUCinf_true, and %s (M0), %s (M1) and %s (M2) for AUClast. The SE of M0 ignores the weight stratum used in the randomization, so it overstates the sampling variability of the treatment difference; M1 includes the stratum.",
          rng(ref_rng("M0")), rng(ref_rng("M1")), rng(ref_rng("M2")), rng(rs("M0"), f3), rng(rs("M1"), f3), rng(rs("M2"), f3), rng(ra("M0"), f3), rng(ra("M1"), f3), rng(ra("M2"), f3)), "",
  sprintf("4. Configurations with AUC0-inf under M1 (cells with point estimate above 5%% / with Wilson lower bound above 5%%, maximum): %s.",
          paste(vapply(setdiff(names(CFG), c("P2", "AUClast_only", "AUCinf_true_only")), function(cf) { f <- fam("M1", cf)
            sprintf("%s %d/16 (%d), %s%%", CFG_EN[[cf]], f$n_over, f$n_lo, f2(f$mx)) }, ""), collapse = "; ")),
  sprintf("   Under M0 the same counts are: %s.", paste(vapply(setdiff(names(CFG), c("P2", "AUClast_only", "AUCinf_true_only")), function(cf) { f <- fam("M0", cf)
            sprintf("%s %d/16 (%d), %s%%", CFG_EN[[cf]], f$n_over, f$n_lo, f2(f$mx)) }, ""), collapse = "; ")), "",
  sprintf("5. Power of P2: identical products (S00) %s%% (M0), %s%% (M1), %s%% (M2) for the 2016 model and %s%%, %s%%, %s%% for Model 1; F x0.97 (true AUC0-inf ratio %s and %s, Cmax ratio %s and %s) %s%%, %s%%, %s%% (2016 model) and %s%%, %s%%, %s%% (Model 1).",
          f2(pwr("k2016", "S00", "M0")$pass_pct), f2(pwr("k2016", "S00", "M1")$pass_pct), f2(pwr("k2016", "S00", "M2")$pass_pct),
          f2(pwr("k2020", "S00", "M0")$pass_pct), f2(pwr("k2020", "S00", "M1")$pass_pct), f2(pwr("k2020", "S00", "M2")$pass_pct),
          f4(TRF[pk_model == "k2016", auc_ratio]), f4(TRF[pk_model == "k2020", auc_ratio]), f4(TRF[pk_model == "k2016", cmax_ratio]), f4(TRF[pk_model == "k2020", cmax_ratio]),
          f2(pwr("k2016", "F097", "M0")$pass_pct), f2(pwr("k2016", "F097", "M1")$pass_pct), f2(pwr("k2016", "F097", "M2")$pass_pct),
          f2(pwr("k2020", "F097", "M0")$pass_pct), f2(pwr("k2020", "F097", "M1")$pass_pct), f2(pwr("k2020", "F097", "M2")$pass_pct)), "",
  "6. Expectations recorded before the results:", "",
  paste0("   - ", EX$expectation, ": ", EX$observed, " (criterion: ", EX$criterion, "): ", ifelse(EX$consistent, "consistent", "not consistent")),
  sprintf("   - Note: in %d cells (absorption rate down) Cmax fails in every trial, so P2 passes in no trial under any model and cannot rise; the registered criterion counts them as not rising.", n_zero), "",
  sprintf("7. Verification: the regenerated M0 equals the stored results in %s rows (%s); maximum relative difference %s.",
          fint(sum(RV$n_rows)), paste(sprintf("%s %s rows", basename(RV$stored_file), fint(RV$n_rows)), collapse = ", "), formatC(max(RV$max_rel_diff), format = "e", digits = 1)), "",
  sprintf("8. Not estimable (pass counted as not passing): M1 %s trial-endpoint results, M2 %s, M0 %s, out of %s each.",
          fint(sum(NE[model == "M1", n_not_estimable])), fint(sum(NE[model == "M2", n_not_estimable])), fint(sum(NE[model == "M0", n_not_estimable])), fint(sum(NE[model == "M0", n_trials]))),
  "", "The choice of the primary analysis model is the sponsor's; M0 and M1 are reported side by side.")
ko <- c("# 분석 모형별 경계 1종 오류 재판정 (M0 pooled t, M1 체중 층 ANOVA, M2 log 체중 ANCOVA)", "",
  sprintf("`scripts/45_oc_models_summary.R`가 `results/oc_models/`에서 생성(사전 등록 `config/prereg_20260926.yaml` section1). 칸·모형당 시험 %s회(연장 칸 %s회).", fint(n_bnd), fint(as.integer(pr$extension$to))), "",
  "## 결론", "",
  paste("1.", first_ko), "",
  sprintf("2. P2 분류(경계 16칸): M0 %s; M1 %s; M2 %s. P2 범위 %s%%(M0), %s%%(M1), %s%%(M2). 같은 시험의 M1 − M0 %s%%p(평균 %s, 16칸 중 %d칸 양수), M2 − M0 %s%%p(평균 %s). M0 → M1 분류가 바뀐 칸: %s.",
          cls_txt_ko("M0", "P2"), cls_txt_ko("M1", "P2"), cls_txt_ko("M2", "P2"), rng(p2_rng("M0"), f2, " ~ "), rng(p2_rng("M1"), f2, " ~ "), rng(p2_rng("M2"), f2, " ~ "),
          rng(d12$diff_pp, s2, " ~ "), s2(mean(d12$diff_pp)), sum(d12$diff_pp > 0), rng(d22$diff_pp, s2, " ~ "), s2(mean(d22$diff_pp)),
          if (nrow(chg)) paste(sprintf("%s %s(%s → %s)", MODEL_KO[chg$pk_model], chg$scenario, CLASS_KO[chg$M0], CLASS_KO[chg$M1]), collapse = "; ") else "없음"), "",
  sprintf("3. 불편 참조(AUCinf_true 단독): M0 %s%%, M1 %s%%, M2 %s%%. 시험 간 log GMR SD / 시험 안 SE 중앙값: AUCinf_true M0 %s, M1 %s, M2 %s; AUClast M0 %s, M1 %s, M2 %s. M0의 SE는 배정에 쓴 체중 층을 반영하지 않아 처리 차이의 표본 변동을 과대 추정하고, M1은 층을 넣는다.",
          rng(ref_rng("M0"), f2, " ~ "), rng(ref_rng("M1"), f2, " ~ "), rng(ref_rng("M2"), f2, " ~ "), rng(rs("M0"), f3, " ~ "), rng(rs("M1"), f3, " ~ "), rng(rs("M2"), f3, " ~ "),
          rng(ra("M0"), f3, " ~ "), rng(ra("M1"), f3, " ~ "), rng(ra("M2"), f3, " ~ ")), "",
  sprintf("4. M1에서 AUC0-inf를 넣은 구성(점추정 > 5%% 칸 수/16 (Wilson 하한 > 5%% 칸 수), 최대): %s.",
          paste(vapply(setdiff(names(CFG), c("P2", "AUClast_only", "AUCinf_true_only")), function(cf) { f <- fam("M1", cf); sprintf("%s %d/16 (%d), %s%%", CFG_EN[[cf]], f$n_over, f$n_lo, f2(f$mx)) }, ""), collapse = "; ")),
  sprintf("   M0: %s.", paste(vapply(setdiff(names(CFG), c("P2", "AUClast_only", "AUCinf_true_only")), function(cf) { f <- fam("M0", cf); sprintf("%s %d/16 (%d), %s%%", CFG_EN[[cf]], f$n_over, f$n_lo, f2(f$mx)) }, ""), collapse = "; ")), "",
  sprintf("5. P2 검정력: 동일 제품(S00) 2016 %s%%(M0), %s%%(M1), %s%%(M2), Model 1 %s%%, %s%%, %s%%; F ×0.97(참 AUC0-inf 비 %s·%s, Cmax 비 %s·%s) 2016 %s%%, %s%%, %s%%, Model 1 %s%%, %s%%, %s%%.",
          f2(pwr("k2016", "S00", "M0")$pass_pct), f2(pwr("k2016", "S00", "M1")$pass_pct), f2(pwr("k2016", "S00", "M2")$pass_pct),
          f2(pwr("k2020", "S00", "M0")$pass_pct), f2(pwr("k2020", "S00", "M1")$pass_pct), f2(pwr("k2020", "S00", "M2")$pass_pct),
          f4(TRF[pk_model == "k2016", auc_ratio]), f4(TRF[pk_model == "k2020", auc_ratio]), f4(TRF[pk_model == "k2016", cmax_ratio]), f4(TRF[pk_model == "k2020", cmax_ratio]),
          f2(pwr("k2016", "F097", "M0")$pass_pct), f2(pwr("k2016", "F097", "M1")$pass_pct), f2(pwr("k2016", "F097", "M2")$pass_pct),
          f2(pwr("k2020", "F097", "M0")$pass_pct), f2(pwr("k2020", "F097", "M1")$pass_pct), f2(pwr("k2020", "F097", "M2")$pass_pct)), "",
  "6. 결과 전 기록한 예상과의 대조:", "",
  paste0("   - ", EX$expectation, ": ", EX$observed, " (기준: ", EX$criterion, "): ", ifelse(EX$consistent, "일치", "불일치")),
  sprintf("   - 주: %d칸(ka 하향)은 모든 시험에서 Cmax가 불통과라 어느 모형에서도 P2 통과가 0이어서 오를 수 없다. 등록 기준은 이 칸을 '오르지 않음'으로 센다.", n_zero), "",
  sprintf("7. 검증: 재생성 M0가 저장본 %s행과 일치(최대 상대 차이 %s).", fint(sum(RV$n_rows)), formatC(max(RV$max_rel_diff), format = "e", digits = 1)), "",
  sprintf("8. 판정 불가(불통과로 셈): M1 %s, M2 %s, M0 %s (각 %s개 시험 × 평가변수 중).", fint(sum(NE[model == "M1", n_not_estimable])), fint(sum(NE[model == "M2", n_not_estimable])),
          fint(sum(NE[model == "M0", n_not_estimable])), fint(sum(NE[model == "M0", n_trials]))), "",
  "주분석 모형 선택은 스폰서 결정이며, 보고서는 M0·M1을 나란히 제시한다.")
writeLines(en, file.path(out_dir, "oc_models_conclusion_en.md")); writeLines(ko, file.path(out_dir, "oc_models_conclusion_ko.md"))
txt <- paste(en, collapse = "\n")
if (grepl("[ㄱ-ㆎ가-힣]|–|—|−", txt)) stop("영문 결론에 한글·en/em dash·U+2212가 있습니다")
cat(en, sep = "\n")
