#!/usr/bin/env Rscript
# 43_p2_interpretation.R — P2 경계 1종 오류의 해석 (지시 2026-09-25 §3). 새 모의 없음.
# (1) 기준: 추정량이 참값 비에 대해 불편이면 두 단측 검정(90% CI가 80.00–125.00% 안)의 경계 1종 오류는 명목 5%에 가깝다(t 기반, df = 2 × 117 − 2 = 232).
#     (a) 이론(R/oc_interpret.R be_boundary_type1_theory): 참값이 정확히 경계면 P(통과) = 5% − P(반대쪽 한계 불통과).
#         se = 저장 CI에서 역산한 (log 상한 − log 하한) / (2·qt(0.95, df_시험))의 시험 중앙값(시나리오·평가변수별), df = 232
#     (b) 경험: AUCinf_true(시험 대상자 개인 모델 적분 AUC0-inf, 참값에 불편) 단독의 경계 통과율(Wilson 95% 구간)
# (2) 분류(모델 × 경계 시나리오 × 구성 P2·AUClast 단독·AUCinf_true 단독; P2·AUClast 단독은 R/oc.R config_pass):
#     Wilson 상한 < 5% 보수적, 하한 > 5% 초과, 그 밖 명목(R/oc_interpret.R classify_type1)
# (3) 원인: 같은 시험 항등식 P2 − 5 = (AUCinf_true − 5) + (AUClast 단독 − AUCinf_true) + (P2 − AUClast 단독)(뒤 두 항은 쌍대 차이, 95% 구간).
#     편향 = mean(log GMR) − log(참값 비)(AUC 계열은 참 AUC0-inf 비, Cmax는 참 Cmax 비), MC 표준오차 = SD / √n.
#     방향: 1에서 멀어지면(0.80에서 음, 1.25에서 양) 1종 오류 감소, 1 쪽이면 증가. scripts/41 p2_bias_boundary.csv와 같은지 확인.
# (4) 문구(ko/en): 수치는 모두 위 표에서, 고정 문구의 전제는 stopifnot(어긋나면 생성 실패). P2 최종 점추정이 5%를 넘는 곳은 이름·크기·구간·시험 수로
#     적고 scripts/33 결론 첫 문단(oc_conclusion_ko.md, oc_conclusion_en.md)·boundary_type1.csv와 대조한다.
# 입력: results/oc/oc_trials_be_<model>.csv.gz(사전 고정 reps_boundary회), 있으면 oc_trials_ext_be_<model>.csv.gz(scripts/42 적응적 연장, 시험 > reps_boundary),
#       oc_scenarios_<model>.csv. 연장 파일이 있으면 합친 값이 최종값이고 사전 고정 reps_boundary회 값도 함께 적는다(열 *_prereg).
# 산출(out_dir, 기본 results/oc): p2_interpretation.csv, p2_interpretation_ko.md, p2_interpretation_en.md
# 사용법: nice -n 15 Rscript scripts/43_p2_interpretation.R [out_dir=results/oc] [in_dir=results/oc]
#   시험용: in_dir를 임시 디렉터리(입력 파일의 링크 + 시험용 연장 파일)로 주면 연장 경로를 운영 파일을 건드리지 않고 확인할 수 있다.
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
in_dir <- if (length(args) >= 2) args[2] else proj_path("results", "oc")
out_dir <- if (length(args) >= 1) args[1] else proj_path("results", "oc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
bnd <- as.numeric(unlist(oc$boundary_targets)); N_BND <- as.integer(oc$trials$reps_boundary); N_ARM <- as.integer(oc$trials$n_per_arm)
CL <- as.numeric(design$be$ci_level); LIMS <- as.numeric(unlist(design$be$limits))
ALPHA <- 100 * (1 - CL) / 2; DF_NOM <- 2L * N_ARM - 2L; Z <- qnorm(0.975)
# 고정 문구("90% CI가 80.00–125.00% 안", "명목 5%", "df = 2 × 117 − 2 = 232")의 전제
stopifnot(abs(CL - 0.90) < 1e-12, isTRUE(all.equal(LIMS, c(0.80, 1.25))), abs(ALPHA - 5) < 1e-9, N_ARM == 117L, DF_NOM == 232L,
          isTRUE(all.equal(sort(bnd), LIMS)), identical(oc$trials$method, "pooled_t"))
MODELS <- intersect(c("k2016", "k2020"), sub("^oc_trials_be_(.*)\\.csv\\.gz$", "\\1", list.files(in_dir, pattern = "^oc_trials_be_.*\\.csv\\.gz$")))
stopifnot(length(MODELS) > 0)
MODEL_KO <- c(k2016 = "Kovalenko 2016 (주)", k2020 = "Kovalenko 2020 Model 1"); MODEL_EN <- c(k2016 = "Kovalenko 2016 (primary)", k2020 = "Kovalenko 2020 Model 1")
MECH_KO <- c(F = "F(흡수량)", ka = "ka(흡수 속도)", ke = "ke(선형 소실)", Vmax = "Vmax(표적 매개 소실)", Km = "Km(결합)", V2 = "V2(분포)")      # scripts/33과 같은 표기
MECH_EN <- c(F = "bioavailability (F)", ka = "absorption rate (ka)", ke = "linear elimination (ke)", Vmax = "target-mediated elimination capacity (Vmax)",
             Km = "binding constant (Km)", V2 = "peripheral volume (V2)")
# 해석 대상 구성: P2·AUClast 단독은 사전 고정 정의(config_pass), AUCinf_true 단독은 불편 기준(참고)
CFG3 <- data.table(config = c("P2", "AUClast_only", "AUCinf_true_only"), col = c("cfg_P2", "cfg_AUClast_only", "ref"), auc_endpoint = c("AUClast", "AUClast", "AUCinf_true"),
                   label_ko = c("P2", "AUClast 단독", "AUCinf_true 단독"), label_en = c("P2", "AUClast only", "AUCinf_true only"))
stopifnot(identical(unlist(oc$configurations$P2$endpoints), c("AUClast", "Cmax")), identical(unlist(oc$configurations$AUClast_only$endpoints), "AUClast"))
CLASS_KO <- c(conservative = "보수적", nominal = "명목", exceeding = "초과"); CLASS_EN <- c(conservative = "conservative", nominal = "nominal", exceeding = "exceeding")
f2 <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2); s2 <- function(x) sprintf("%+.2f", round(x, 2) + 0)
f3 <- function(x) formatC(round(x, 3) + 0, format = "f", digits = 3); f4 <- function(x) formatC(round(x, 4) + 0, format = "f", digits = 4)
fint <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
fe <- function(x) formatC(x, format = "e", digits = 1)
rng <- function(x, fmt = f2, sep = "–") if (isTRUE(all.equal(min(x), max(x)))) fmt(min(x)) else paste0(fmt(min(x)), sep, fmt(max(x)))
fb <- function(x) ifelse(abs(round(x, 2) - 5) < 1e-9 & x != 5, f3(x), f2(x))   # Wilson 끝점: 반올림하면 5.00이 되는 값은 셋째 자리까지(분류가 보이게)
ci_k <- function(p, l, h) sprintf("%s%% [%s, %s]", f2(p), fb(l), fb(h)); ci_e <- function(p, l, h) sprintf("%s%% (95%% CI %s to %s)", f2(p), fb(l), fb(h))

# ---- 입력 ----------------------------------------------------------------------------------------------------------
COLS <- c("trial", "scenario", "endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T")
scen_info <- function(mdl) {
  s <- fread(file.path(in_dir, sprintf("oc_scenarios_%s.csv", mdl)))
  s <- s[s$reachable %in% TRUE & (abs(s$target - bnd[1]) < 1e-9 | abs(s$target - bnd[2]) < 1e-9)]
  stopifnot(nrow(s) > 0, !anyDuplicated(s$code), all(s$within_tol %in% TRUE), all(abs(s$auc_ratio / s$target - 1) <= as.numeric(oc$inversion$tolerance_rel) + 1e-12))
  s[, .(scenario = code, mechanism, direction, target, multiplier, auc_ratio, auc_se_log, cmax_ratio, cmax_se_log)]
}
# 경계 시나리오의 원래 6개 평가변수. 6개가 모두 있는 (시험, 시나리오)만 쓴다(추가 중인 파일의 불완전 묶음 제외, 수를 기록)
read_be <- function(f, scn) {
  x <- fread(f)
  if (!all(COLS %in% names(x))) stop("열이 없습니다: ", f, " (", paste(setdiff(COLS, names(x)), collapse = ", "), ")")
  x <- x[x$scenario %in% scn & x$endpoint %in% OC_ENDPOINTS, COLS, with = FALSE]
  if (anyDuplicated(x, by = c("trial", "scenario", "endpoint"))) stop("(trial, scenario, endpoint) 중복: ", f)
  k <- x[, .N, by = .(trial, scenario)]
  list(be = x[k[k$N == length(OC_ENDPOINTS), .(trial, scenario)], on = c("trial", "scenario")], n_incomplete = sum(k$N != length(OC_ENDPOINTS)))
}

# ---- 분석(한 자료 집합) ------------------------------------------------------------------------------------------------
analyze <- function(be, info) {
  be <- copy(be)
  be[, df := n_R + n_T - 2]
  be[, se := (log(CI_upper) - log(CI_lower)) / (2 * qt(1 - (1 - CL) / 2, df))]      # 저장 CI(pooled t)에서 역산한 log 척도 표준오차
  w <- config_pass(be, oc)
  stopifnot(identical(w$cfg_P2, (w$AUClast %in% TRUE) & (w$Cmax %in% TRUE)), identical(w$cfg_AUClast_only, w$AUClast %in% TRUE))
  w[, `:=`(ref = AUCinf_true %in% TRUE, cmx = Cmax %in% TRUE)]
  stopifnot(!any(w$cfg_P2 & !w$cfg_AUClast_only), !any(w$cfg_P2 & !w$cmx))              # 구성 정의상 P2 통과 ⇒ AUClast·Cmax 통과
  rt <- rbindlist(lapply(seq_len(nrow(CFG3)), function(i) {
    cc <- CFG3$col[i]
    w[, { x <- get(cc); wc <- wilson_ci(sum(x), .N); .(n_trials = .N, n_pass = sum(x), pass_pct = wc$est, lo = wc$lo, hi = wc$hi) }, by = "scenario"][, config := CFG3$config[i]][]
  }))
  rt[, class := classify_type1(lo, hi, ALPHA)]
  pd <- w[, {
    a <- paired_prop_diff_ci(cfg_AUClast_only, ref); b <- paired_prop_diff_ci(cfg_P2, cfg_AUClast_only)
    .(n_w = .N, p2_pct = 100 * mean(cfg_P2), al_pct = 100 * mean(cfg_AUClast_only), ref_pct = 100 * mean(ref), cmax_pass_pct = 100 * mean(cmx),
      auclast_minus_ref_pp = a$est, auclast_minus_ref_lo = a$lo, auclast_minus_ref_hi = a$hi,
      p2_minus_auclast_pp = b$est, p2_minus_auclast_lo = b$lo, p2_minus_auclast_hi = b$hi, n_auclast_pass_cmax_fail = sum(cfg_AUClast_only & !cmx)) }, by = "scenario"]
  pd[, ref_minus_5_pp := ref_pct - ALPHA]
  # 항등식(같은 시험 집합에서 정확히): P2 − 5 = (AUCinf_true − 5) + (AUClast 단독 − AUCinf_true) + (P2 − AUClast 단독)
  stopifnot(all(abs((pd$p2_pct - ALPHA) - (pd$ref_minus_5_pp + pd$auclast_minus_ref_pp + pd$p2_minus_auclast_pp)) < 1e-9),
            all(pd$p2_minus_auclast_pp <= 1e-12), all(abs(pd$p2_minus_auclast_pp + 100 * pd$n_auclast_pass_cmax_fail / pd$n_w) < 1e-9))
  bb <- be[is.finite(GMR) & GMR > 0, .(n_bias = .N, mean_log_gmr = mean(log(GMR)), sd_log_gmr = sd(log(GMR)), se_med = median(se), df_med = median(df),
                                       df_nom_share = mean(df == DF_NOM), cor_se_lgmr = cor(se, log(GMR))), by = .(scenario, endpoint)]
  bb <- merge(info, bb, by = "scenario")
  bb[, true_ratio := fifelse(endpoint == "Cmax", cmax_ratio, auc_ratio)]
  bb[, `:=`(bias_log = mean_log_gmr - log(true_ratio), bias_mc_se = sd_log_gmr / sqrt(n_bias))]
  bb[, `:=`(bias_lo = bias_log - Z * bias_mc_se, bias_hi = bias_log + Z * bias_mc_se, bias_pct = 100 * (exp(bias_log) - 1), sd_se_ratio = sd_log_gmr / se_med)]
  bb[, bias_dir := bias_direction(bias_log, true_ratio, bias_lo, bias_hi)]
  # 이론·정규 근사(AUC 계열: 참 AUC0-inf 비가 경계. Cmax는 참값이 경계가 아니므로 제외). 이론은 참값이 정확히 경계(목표)일 때와 시나리오의 실제 참값(역산 허용 오차 ±0.1%)일 때
  ia <- which(bb$endpoint != "Cmax")
  th_b <- be_boundary_type1_theory(bb$se_med[ia], bb$df_med[ia], bb$target[ia], CL, LIMS)
  th_t <- be_boundary_type1_theory(bb$se_med[ia], bb$df_med[ia], bb$auc_ratio[ia], CL, LIMS)
  stopifnot(all(abs(th_b$p_near - (1 - CL) / 2) < 1e-12), all(abs(th_b$p_pass - (th_b$p_near - th_b$p_far)) < 1e-15))   # 경계: 가까운 한계 5%, 통과 = 5% − 먼 한계
  bb[ia, `:=`(theory_boundary_pct = 100 * th_b$p_pass, theory_far_prob = th_b$p_far, theory_at_truth_pct = 100 * th_t$p_pass,
              pred_unbiased_pct = 100 * be_pass_prob_normal(log(true_ratio), sd_log_gmr, se_med, df_med, CL, LIMS),
              pred_pct = 100 * be_pass_prob_normal(mean_log_gmr, sd_log_gmr, se_med, df_med, CL, LIMS))]
  list(w = w, rt = rt, pd = pd, bb = bb)
}

# ---- 모델별 실행 ------------------------------------------------------------------------------------------------------
A_FIN <- list(); A_PRE <- list(); ST <- list(); INFO <- list()
for (mdl in MODELS) {
  info <- scen_info(mdl); INFO[[mdl]] <- info
  mn <- read_be(file.path(in_dir, sprintf("oc_trials_be_%s.csv.gz", mdl)), info$scenario)
  stopifnot(setequal(unique(mn$be$scenario), info$scenario))
  pre <- mn$be[trial <= N_BND]                                                          # 사전 고정 시험 1–reps_boundary
  ext_f <- file.path(in_dir, sprintf("oc_trials_ext_be_%s.csv.gz", mdl)); has_ext <- file.exists(ext_f); n_ext_inc <- 0L
  fin <- pre
  if (has_ext) {
    ex <- read_be(ext_f, info$scenario); n_ext_inc <- ex$n_incomplete
    stopifnot(all(ex$be$trial > N_BND), !nrow(ex$be[pre, on = c("trial", "scenario", "endpoint"), nomatch = NULL]))   # 연장분은 사전 고정 범위 밖(scripts/42)
    for (sc_ in unique(ex$be$scenario)) {                               # scripts/33·41과 같은 조건: 사전 고정 시험 1..reps_boundary 뒤에 빈틈 없이 이어진다
      tp <- sort(unique(pre[scenario == sc_, trial])); te <- sort(unique(ex$be[scenario == sc_, trial]))
      if (!identical(tp, seq_len(N_BND)) || !identical(te, N_BND + seq_along(te))) stop(sprintf("%s %s: 연장 시험이 사전 고정 시험 1-%d 뒤에 빈틈 없이 이어지지 않습니다", mdl, sc_, N_BND))
    }
    fin <- rbind(pre, ex$be)
  }
  stopifnot(nrow(mn$be[trial > N_BND]) == 0)                                           # 사전 고정 파일에 reps_boundary 초과 시험 없음
  A_PRE[[mdl]] <- analyze(pre, info)
  A_FIN[[mdl]] <- if (has_ext) analyze(fin, info) else A_PRE[[mdl]]
  n_pre <- A_PRE[[mdl]]$pd[, .(scenario, n_pre = n_w)]; n_fin <- A_FIN[[mdl]]$pd[, .(scenario, n_fin = n_w)]
  nn <- merge(n_pre, n_fin, by = "scenario")
  ST[[mdl]] <- data.table(model = mdl, n_scen = nrow(info), n_pre_min = min(nn$n_pre), n_pre_max = max(nn$n_pre), complete = min(nn$n_pre) >= N_BND,
                          has_ext = has_ext, ext_scen = paste(nn[n_fin > n_pre, sprintf("%s:%d", scenario, n_fin - n_pre)], collapse = ";"),
                          n_incomplete_main = mn$n_incomplete, n_incomplete_ext = n_ext_inc)
}
ST <- rbindlist(ST)

# ---- scripts/41 p2_bias_boundary.csv 교차 확인: 같은 시험 수면 값이 같아야 하고(1e-9), 시험 수가 다르면(41이 부분 파일로 실행) 앞 n회로 다시 계산해 비교 --------
b41_f <- file.path(in_dir, "p2_bias_boundary.csv")
B41 <- if (file.exists(b41_f)) fread(b41_f)[endpoint %in% OC_ENDPOINTS] else NULL
bias_on <- function(be, info) {
  x <- be[is.finite(GMR) & GMR > 0, .(n = .N, m = mean(log(GMR)), s = sd(log(GMR))), by = .(scenario, endpoint)]
  x <- merge(info[, .(scenario, auc_ratio, cmax_ratio)], x, by = "scenario")
  x[, .(scenario, endpoint, n, mean_log_gmr = m, bias_log = m - log(fifelse(endpoint == "Cmax", cmax_ratio, auc_ratio)), bias_mc_se = s / sqrt(n))]
}
CHK41 <- rbindlist(lapply(MODELS, function(mdl) {
  sc <- INFO[[mdl]]$scenario
  if (is.null(B41) || !nrow(B41[model == mdl])) return(data.table(model = mdl, scenario = sc, bias41_status = "absent", bias41_n_trials = NA_integer_))
  mine <- A_FIN[[mdl]]$bb[, .(scenario, endpoint, n = n_bias, mean_log_gmr, bias_log, bias_mc_se)]
  b <- merge(B41[model == mdl, .(scenario, endpoint, n41 = n_trials, m41 = mean_log_gmr, b41 = bias_log, s41 = bias_mc_se)], mine, by = c("scenario", "endpoint"), all = TRUE)
  if (anyNA(b$n41) || anyNA(b$n)) stop("p2_bias_boundary.csv와 경계 시나리오·평가변수 집합이 다릅니다(", mdl, ")")
  same <- b[n41 == n]
  stopifnot(all(abs(same$m41 - same$mean_log_gmr) < 1e-9), all(abs(same$b41 - same$bias_log) < 1e-9), all(abs(same$s41 - same$bias_mc_se) < 1e-9))
  pre_be <- NULL
  out <- b[, {
    if (all(n41 == n)) .(bias41_status = "equal", bias41_n_trials = as.integer(n[1]))
    else {
      if (is.null(pre_be)) { f <- file.path(in_dir, sprintf("oc_trials_be_%s.csv.gz", mdl)); pre_be <<- read_be(f, sc)$be }
      k <- unique(n41); stopifnot(length(k) == 1)
      sub_ <- bias_on(pre_be[pre_be$scenario == .BY$scenario & pre_be$trial <= k], INFO[[mdl]])
      e <- merge(.SD[, .(endpoint, m41, b41, s41)], sub_, by = "endpoint")
      ok <- nrow(e) == .N && all(e$n == k) && all(abs(e$m41 - e$mean_log_gmr) < 1e-9) && all(abs(e$s41 - e$bias_mc_se) < 1e-9)
      .(bias41_status = if (ok) "stale_equal_on_subset" else "stale_differs", bias41_n_trials = as.integer(k))
    } }, by = "scenario"]
  out[, model := mdl][]
}))
b41_txt <- function(st, n, lang) {
  if (lang == "ko") switch(st, equal = sprintf("같음(시험 %s회)", fint(n)), stale_equal_on_subset = sprintf("41은 시험 %s회로 실행(시험 1-%s에서 같음)", fint(n), fint(n)),
                           stale_differs = sprintf("41은 시험 %s회로 실행(시험 1-%s에서도 다름)", fint(n), fint(n)), absent = "41 파일 없음")
  else switch(st, equal = sprintf("equal (%s trials)", fint(n)), stale_equal_on_subset = sprintf("stale: scripts/41 used %s trials (equal on trials 1-%s); rerun scripts/41", fint(n), fint(n)),
              stale_differs = sprintf("stale: scripts/41 used %s trials (NOT equal on trials 1-%s); rerun scripts/41", fint(n), fint(n)), absent = "p2_bias_boundary.csv absent")
}
CHK41[, bias41_check := mapply(b41_txt, bias41_status, bias41_n_trials, MoreArgs = list(lang = "en"))]

# ---- 결과 표(모델 × 경계 시나리오 × 구성) --------------------------------------------------------------------------------
mk_table <- function(a, a_pre, info, mdl) {
  rt <- merge(info, a$rt, by = "scenario")
  rt[, auc_endpoint := CFG3$auc_endpoint[match(config, CFG3$config)]]
  b <- a$bb[endpoint %in% c("AUClast", "AUCinf_true"), .(scenario, auc_endpoint = endpoint, se_median = se_med, df_median = df_med, df232_share = df_nom_share, sd_log_gmr,
            sd_se_ratio, cor_se_log_gmr = cor_se_lgmr, theory_boundary_pct, theory_far_prob, theory_at_truth_pct, pred_unbiased_pct, pred_pct,
            bias_log, bias_mc_se, bias_lo, bias_hi, bias_pct, bias_dir)]
  rt <- merge(rt, b, by = c("scenario", "auc_endpoint"))
  rf <- a$rt[config == "AUCinf_true_only", .(scenario, ref_lo = lo, ref_hi = hi, ref_class = class)]
  cm <- a$bb[endpoint == "Cmax", .(scenario, cmax_bias_log = bias_log, cmax_bias_mc_se = bias_mc_se, cmax_bias_pct = bias_pct, cmax_bias_dir = bias_dir)]
  rt <- merge(merge(merge(rt, a$pd, by = "scenario"), rf, by = "scenario"), cm, by = "scenario")
  wb <- dcast(a$bb, scenario ~ endpoint, value.var = c("bias_pct", "bias_mc_se"))
  setnames(wb, sub("^bias_mc_se_(.*)$", "bias_mc_se_pct_\\1", names(wb)))
  for (cc in grep("^bias_mc_se_pct_", names(wb), value = TRUE)) set(wb, j = cc, value = 100 * wb[[cc]])   # MC 표준오차(100 × log)
  rt <- merge(rt, wb, by = "scenario")
  pr <- a_pre$rt[, .(scenario, config, n_trials_prereg = n_trials, pass_pct_prereg = pass_pct, lo_prereg = lo, hi_prereg = hi, class_prereg = class)]
  rt <- merge(rt, pr, by = c("scenario", "config"))
  rt[, `:=`(model = mdl, gap_pp = pass_pct - ALPHA)]
  rt[]
}
TAB <- rbindlist(lapply(MODELS, function(mdl) mk_table(A_FIN[[mdl]], A_PRE[[mdl]], INFO[[mdl]], mdl)), use.names = TRUE)
TAB <- merge(TAB, CHK41, by = c("model", "scenario"))
# 표 값 = 자료: 시험 수는 시나리오별 행 수, 분류는 Wilson 구간에서 다시 계산
for (mdl in MODELS) {
  nw <- A_FIN[[mdl]]$w[, .N, by = "scenario"]
  x <- TAB[model == mdl][nw, on = "scenario"]
  stopifnot(all(x$n_trials == x$N), all(x$n_trials == x$n_w))
}
stopifnot(identical(TAB$class, classify_type1(TAB$lo, TAB$hi, ALPHA)), all(TAB$n_trials >= TAB$n_trials_prereg),
          all(TAB[class == "conservative", hi < ALPHA]), all(TAB[class == "exceeding", lo > ALPHA]), all(TAB[class == "nominal", lo <= ALPHA & hi >= ALPHA]),
          all(TAB[auc_endpoint %in% c("AUClast", "AUCinf_true"), df_median == DF_NOM]))
# 편향 방향 = 부호(참값 < 1이면 음이 1에서 멀어짐)
stopifnot(all(TAB[bias_dir == "away_from_1", sign(bias_log) == sign(log(target))]), all(TAB[bias_dir == "toward_1", sign(bias_log) == -sign(log(target))]),
          all(TAB[bias_dir == "none", bias_lo <= 0 & bias_hi >= 0]))

# ---- 한 줄 원인 ------------------------------------------------------------------------------------------------------
# 요소: 기준(AUCinf_true − 5, Wilson 구간이 5% 제외면 유의), AUClast(AUClast 단독 − AUCinf_true, 쌍대 95% 구간이 0 제외면 유의), Cmax(P2 − AUClast 단독, 상한 < 0이면 유의).
# 유의한 요소를 크기 순으로 적고, AUClast 요소는 편향 방향과 부호가 맞을 때만 "편향 때문" 문구를 쓴다. Cmax가 차이(P2 − 5)의 절반 이상이면 "대부분".
cause_one <- function(r, lang) {
  ko <- lang == "ko"; comps <- list()
  gap <- r$pass_pct - ALPHA
  if (r$ref_class != "nominal") {
    stopifnot((r$ref_class == "conservative") == (r$ref_minus_5_pp < 0))
    comps$ref <- list(v = r$ref_minus_5_pp, txt = if (ko) sprintf("편향 없는 기준도 5%% %s(AUCinf_true %s%%, %s%%p)", if (r$ref_minus_5_pp < 0) "아래" else "위", f2(r$ref_pct), s2(r$ref_minus_5_pp))
                                              else sprintf("the unbiased reference is itself %s 5%% (AUCinf_true %s%%, %s points)", if (r$ref_minus_5_pp < 0) "below" else "above", f2(r$ref_pct), s2(r$ref_minus_5_pp)))
  }
  le <- r$auclast_minus_ref_pp
  if (r$auclast_minus_ref_lo > 0 || r$auclast_minus_ref_hi < 0) {
    d <- r$bias_dir; cons <- (d == "away_from_1" && le < 0) || (d == "toward_1" && le > 0)
    txt <- if (d == "none") {
      if (ko) sprintf("AUClast 편향은 0과 구별되지 않고(편향 %s%%) AUClast 단독 − AUCinf_true %s%%p", s2(r$bias_pct), s2(le))
      else sprintf("AUClast bias not distinguishable from 0 (bias %s%%), AUClast only minus AUCinf_true %s points", s2(r$bias_pct), s2(le))
    } else if (cons) {
      if (ko) sprintf("AUClast 편향이 %s(편향 %s%%; AUClast 단독 − AUCinf_true %s%%p)", if (d == "away_from_1") "1에서 멀어지는 방향" else "1 쪽", s2(r$bias_pct), s2(le))
      else sprintf("AUClast bias %s (bias %s%%; AUClast only minus AUCinf_true %s points)", if (d == "away_from_1") "away from 1" else "toward 1", s2(r$bias_pct), s2(le))
    } else {
      if (ko) sprintf("AUClast 단독 − AUCinf_true %s%%p(AUClast 편향 %s%%는 %s: 방향 불일치, 편향 외 요인)", s2(le), s2(r$bias_pct), if (d == "away_from_1") "1에서 멀어지는 방향" else "1 쪽")
      else sprintf("AUClast only minus AUCinf_true %s points (opposite to the AUClast bias of %s%%, which is %s: a factor other than bias)", s2(le), s2(r$bias_pct), if (d == "away_from_1") "away from 1" else "toward 1")
    }
    comps$last <- list(v = le, txt = txt, cons = cons, dir = d)
  }
  if (r$config == "P2" && r$p2_minus_auclast_hi < 0) {
    dom <- gap < 0 && abs(r$p2_minus_auclast_pp) >= 0.5 * abs(gap)
    comps$cmax <- list(v = r$p2_minus_auclast_pp, dom = dom,
                       txt = if (ko) sprintf("Cmax 불통과%s(Cmax 통과 %s%%; P2 − AUClast 단독 %s%%p)", if (dom) "가 대부분" else "", f2(r$cmax_pass_pct), s2(r$p2_minus_auclast_pp))
                             else sprintf("Cmax failures%s (Cmax pass %s%%; P2 minus AUClast only %s points)", if (dom) " dominate" else "", f2(r$cmax_pass_pct), s2(r$p2_minus_auclast_pp)))
  }
  out <- switch(r$class,
    conservative = if (ko) "보수적" else "conservative",
    nominal = if (r$pass_pct > ALPHA) { if (ko) "점추정 5% 초과, Wilson 구간 5% 포함(명목)" else "point estimate above 5% with the Wilson interval including 5% (nominal)" } else { if (ko) "명목" else "nominal" },
    exceeding = if (ko) "5% 초과(Wilson 하한 > 5%)" else "above 5% (Wilson lower bound above 5%)")
  lab <- if (ko) CFG3[config == r$config, label_ko] else CFG3[config == r$config, label_en]
  o <- if (length(comps)) order(-abs(vapply(comps, `[[`, 0, "v"))) else integer(0)
  head_ <- if (length(comps)) paste(vapply(comps[o], `[[`, "", "txt"), collapse = " + ") else if (ko) "유의한 요인 없음" else "no component distinguishable from 0"
  line <- if (ko) sprintf("%s → %s %s%%: %s", head_, lab, f2(r$pass_pct), out) else sprintf("%s; hence %s %s%%: %s", head_, lab, f2(r$pass_pct), out)
  list(line = line, primary = if (length(comps)) names(comps)[o[1]] else "none", last_cons = if (!is.null(comps$last)) comps$last$cons else NA,
       cmax_dom = if (!is.null(comps$cmax)) comps$cmax$dom else FALSE)
}
ref_note <- function(r, lang) {
  inc0 <- r$bias_lo <= 0 && r$bias_hi >= 0
  if (lang == "ko") sprintf("편향 없는 기준: 편향 %s%% [%s, %s]%s, 시험 간 SD / 시험 안 SE %s, 정규 근사 %s%%(관측 %s%%)", s2(r$bias_pct), s2(100 * r$bias_lo), s2(100 * r$bias_hi),
                            if (inc0) "(0 포함)" else "(0 제외)", f3(r$sd_se_ratio), f2(r$pred_pct), f2(r$pass_pct))
  else sprintf("unbiased reference: bias %s%% (95%% CI %s to %s, %s), between-trial SD / within-trial SE %s, normal approximation %s%% (observed %s%%)", s2(r$bias_pct), s2(100 * r$bias_lo),
               s2(100 * r$bias_hi), if (inc0) "includes 0" else "excludes 0", f3(r$sd_se_ratio), f2(r$pred_pct), f2(r$pass_pct))
}
TAB[, `:=`(cause_ko = NA_character_, cause_en = NA_character_, cause_primary = NA_character_, auclast_bias_consistent = NA, cmax_dominant = NA)]
for (i in seq_len(nrow(TAB))) {
  r <- TAB[i]
  if (r$config == "AUCinf_true_only") { set(TAB, i, c("cause_ko", "cause_en"), list(ref_note(r, "ko"), ref_note(r, "en"))); next }
  ck <- cause_one(r, "ko"); ce <- cause_one(r, "en")
  stopifnot(identical(ck$primary, ce$primary), identical(ck$last_cons, ce$last_cons))
  set(TAB, i, c("cause_ko", "cause_en", "cause_primary", "auclast_bias_consistent", "cmax_dominant"), list(ck$line, ce$line, ck$primary, ck$last_cons, ck$cmax_dom))
}
# 원인 문구의 전제: "1에서 멀어지는 방향"/"1 쪽" 문구는 부호 일치 행에서만, "대부분"은 Cmax가 차이의 절반 이상일 때만
stopifnot(all(TAB[grepl("AUClast 편향이 1에서 멀어지는 방향", cause_ko), bias_dir == "away_from_1" & auclast_minus_ref_pp < 0]),
          all(TAB[grepl("AUClast 편향이 1 쪽", cause_ko), bias_dir == "toward_1" & auclast_minus_ref_pp > 0]),
          all(TAB[grepl("Cmax 불통과가 대부분", cause_ko), config == "P2" & gap_pp < 0 & abs(p2_minus_auclast_pp) >= 0.5 * abs(gap_pp) & p2_minus_auclast_hi < 0]),
          all(TAB[grepl("%: 보수적$", cause_ko), class == "conservative"]), all(TAB[config != "AUCinf_true_only" & class == "conservative", grepl("%: 보수적$", cause_ko)]), all(TAB[grepl("점추정 5% 초과", cause_ko), class == "nominal" & pass_pct > ALPHA]))

ord_cols <- c("model", "scenario", "mechanism", "direction", "target", "multiplier", "auc_ratio", "cmax_ratio", "config", "auc_endpoint",
              "n_trials", "n_pass", "pass_pct", "lo", "hi", "class", "gap_pp", "n_trials_prereg", "pass_pct_prereg", "lo_prereg", "hi_prereg", "class_prereg",
              "se_median", "df_median", "df232_share", "theory_boundary_pct", "theory_far_prob", "theory_at_truth_pct", "sd_log_gmr", "sd_se_ratio", "cor_se_log_gmr",
              "pred_unbiased_pct", "pred_pct", "bias_log", "bias_mc_se", "bias_lo", "bias_hi", "bias_pct", "bias_dir",
              "ref_pct", "ref_lo", "ref_hi", "ref_class", "ref_minus_5_pp", "al_pct", "auclast_minus_ref_pp", "auclast_minus_ref_lo", "auclast_minus_ref_hi",
              "p2_pct", "p2_minus_auclast_pp", "p2_minus_auclast_lo", "p2_minus_auclast_hi", "n_auclast_pass_cmax_fail", "cmax_pass_pct",
              "cmax_bias_log", "cmax_bias_mc_se", "cmax_bias_pct", "cmax_bias_dir",
              paste0("bias_pct_", OC_ENDPOINTS), paste0("bias_mc_se_pct_", OC_ENDPOINTS),
              "cause_primary", "auclast_bias_consistent", "cmax_dominant", "cause_ko", "cause_en", "bias41_status", "bias41_n_trials", "bias41_check")
stopifnot(all(ord_cols %in% names(TAB)))
TAB <- TAB[, ord_cols, with = FALSE][order(match(model, MODELS), target, mechanism, direction, match(config, CFG3$config))]

# ---- scripts/33 대조: boundary_type1.csv(P2·AUClast 단독)와 결론 첫 문단 ------------------------------------------------------
# 33의 시험 수가 최종(연장 포함) 또는 사전 고정 값과 같은 기준을 찾아 그 기준으로 정확히 대조한다. 어느 쪽과도 다르면(33이 오래됨) 문구에 적고 대조는 생략한다.
bt_f <- file.path(in_dir, "boundary_type1.csv"); cc_ko_f <- file.path(in_dir, "oc_conclusion_ko.md"); cc_en_f <- file.path(in_dir, "oc_conclusion_en.md")
C33 <- list(status = "absent", basis = NA_character_, note = "")
if (file.exists(bt_f)) {
  bt <- fread(bt_f)[config %in% c("P2", "AUClast_only") & model %in% MODELS]
  pick <- function(pfx) {
    x <- if (pfx == "final") TAB[config %in% c("P2", "AUClast_only"), .(model, scenario, config, n = n_trials, p = pass_pct, l = lo, h = hi)]
         else TAB[config %in% c("P2", "AUClast_only"), .(model, scenario, config, n = n_trials_prereg, p = pass_pct_prereg, l = lo_prereg, h = hi_prereg)]
    m <- merge(x, bt[, .(model, scenario, config, n33 = n_trials, p33 = pass_pct, l33 = lo, h33 = hi)], by = c("model", "scenario", "config"), all.x = TRUE)
    if (anyNA(m$n33) || any(m$n33 != m$n)) return(NULL)
    stopifnot(all(abs(m$p33 - m$p) < 1e-9), all(abs(m$l33 - m$l) < 1e-9), all(abs(m$h33 - m$h) < 1e-9))   # 같은 시험이면 같은 값
    x[config == "P2"]
  }
  P33 <- pick("final"); basis <- "final"
  if (is.null(P33)) { P33 <- pick("prereg"); basis <- "prereg" }
  if (is.null(P33)) {
    C33 <- list(status = "stale", basis = NA_character_, note = sprintf("boundary_type1.csv trial counts (%s) differ from these data", paste(unique(bt[config == "P2", n_trials]), collapse = "/")))
    warning("scripts/33 결과(boundary_type1.csv)의 시험 수가 이 자료와 다릅니다. scripts/33을 다시 실행하세요.")
  } else {
    C33 <- list(status = "checked", basis = basis, note = "")
    ex33 <- P33[p > ALPHA][order(-p)]
    info_all <- rbindlist(lapply(MODELS, function(mdl) INFO[[mdl]][, model := mdl]))
    ex33 <- merge(ex33, info_all, by = c("model", "scenario"))[order(-p)]
    # 33 문단의 수치 표기 자릿수(소수 1자리 또는 2자리)를 읽어 같은 자릿수로 비교한다. 시험 수가 적혀 있으면("시험 N회", "N trials") 그것도 비교
    same_num <- function(txt, x) { d <- nchar(sub("^[^.]*\\.?", "", txt)); identical(txt, formatC(x, format = "f", digits = d)) }
    chk_para <- function(f, lang) {
      if (!file.exists(f)) return(FALSE)
      para <- readLines(f, encoding = "UTF-8", warn = FALSE)[1]
      if (lang == "ko") {
        rx <- "참값 ([0-9.]+)\\(배율 ([^)]+)\\): ([0-9.]+)% \\[([0-9.]+), ([0-9.]+)\\](?:, 시험 ([0-9,]+)회)?"
        rx_max <- "최대 ([0-9.]+)% \\[([0-9.]+), ([0-9.]+)\\]"
        yes <- grepl("5%를 넘는 경우가 있다", para, fixed = TRUE); no <- grepl("5% 이하다", para, fixed = TRUE)
      } else {
        rx <- "true ratio ([0-9.]+) \\(multiplier ([^)]+)\\): ([0-9.]+)% \\(95% CI ([0-9.]+) to ([0-9.]+)(?:, ([0-9,]+) trials)?\\)"
        rx_max <- "the largest is ([0-9.]+)% \\(95% CI ([0-9.]+) to ([0-9.]+)\\)"
        yes <- grepl("exceeds 5% in some cases", para, fixed = TRUE); no <- grepl("is at most 5%", para, fixed = TRUE)
      }
      stopifnot(xor(yes, no), yes == (nrow(ex33) > 0))                                     # 첫 문단의 5% 초과 여부 = 이 자료
      if (nrow(ex33)) {
        segs <- strsplit(para, "; ", fixed = TRUE)[[1]]
        mm <- regmatches(para, gregexpr(rx, para, perl = TRUE))[[1]]
        stopifnot(length(mm) == nrow(ex33))
        v <- regmatches(mm, regexec(rx, mm, perl = TRUE))
        for (i in seq_len(nrow(ex33))) {
          r <- ex33[i]; g <- v[[i]]
          seg <- segs[grepl(mm[i], segs, fixed = TRUE)]
          lab_m <- if (lang == "ko") MODEL_KO[[r$model]] else c(k2016 = "2016 model", k2020 = "Model 1")[[r$model]]
          lab_c <- if (lang == "ko") MECH_KO[[r$mechanism]] else MECH_EN[[r$mechanism]]
          stopifnot(length(seg) == 1, grepl(lab_m, seg, fixed = TRUE), grepl(lab_c, seg, fixed = TRUE),
                    g[2] == sprintf("%.2f", r$target), g[3] == sprintf("%.3g", r$multiplier), same_num(g[4], r$p), same_num(g[5], r$l), same_num(g[6], r$h),
                    length(g) < 7 || !nzchar(g[7]) || as.integer(gsub(",", "", g[7])) == r$n)
        }
      } else {
        mx <- P33[which.max(p)]; g <- regmatches(para, regexec(rx_max, para))[[1]]
        stopifnot(length(g) == 4, same_num(g[2], mx$p), same_num(g[3], mx$l), same_num(g[4], mx$h))
      }
      TRUE
    }
    C33$ko <- chk_para(cc_ko_f, "ko"); C33$en <- chk_para(cc_en_f, "en")
  }
}

# ---- 문구 --------------------------------------------------------------------------------------------------------------
B <- TAB; P2 <- B[config == "P2"]; AL <- B[config == "AUClast_only"]; RF <- B[config == "AUCinf_true_only"]
scn_lab <- function(r, lang) if (lang == "ko") sprintf("%s(참값 %.2f, ×%.3g)", r$scenario, r$target, r$multiplier) else sprintf("%s (true %.2f, x%.3g)", r$scenario, r$target, r$multiplier)
full_lab <- function(r, lang) {
  if (lang == "ko") sprintf("%s · %s %s(배율 ×%.3g) · 참값 %.2f", MODEL_KO[[r$model]], MECH_KO[[r$mechanism]], r$direction, r$multiplier, r$target)
  else sprintf("%s, %s, %s (multiplier x%.3g), true ratio %.2f", MODEL_EN[[r$model]], MECH_EN[[r$mechanism]], r$direction, r$multiplier, r$target)
}
trials_txt <- function(mdl, lang) {
  s <- ST[model == mdl]; x <- B[model == mdl & config == "P2"]
  n <- if (min(x$n_trials) == max(x$n_trials)) fint(min(x$n_trials)) else sprintf("%s-%s", fint(min(x$n_trials)), fint(max(x$n_trials)))
  ex <- x[n_trials > n_trials_prereg]
  if (lang == "ko") sprintf("%s: 경계 시나리오 %d개, 시나리오당 %s회%s%s%s", MODEL_KO[[mdl]], s$n_scen, n, if (s$complete) "(사전 고정 reps_boundary 완료)" else "(사전 고정 파일이 아직 부분: 완료 후 같은 코드로 다시 생성)",
                            if (nrow(ex)) sprintf("; 적응적 연장(scripts/42): %s", paste(sprintf("%s %s회(사전 고정 %s + 연장 %s)", ex$scenario, fint(ex$n_trials), fint(ex$n_trials_prereg), fint(ex$n_trials - ex$n_trials_prereg)), collapse = ", ")) else "",
                            if (s$n_incomplete_main + s$n_incomplete_ext > 0) sprintf("; 평가변수 6개가 다 있지 않은 (시험, 시나리오) %d개 제외(추가 중인 묶음)", s$n_incomplete_main + s$n_incomplete_ext) else "")
  else sprintf("%s: %d boundary scenarios, %s trials per scenario%s%s%s", MODEL_EN[[mdl]], s$n_scen, n, if (s$complete) " (prespecified reps_boundary complete)" else " (prespecified file still partial; regenerate with the same code when complete)",
               if (nrow(ex)) sprintf("; adaptive extension (scripts/42): %s", paste(sprintf("%s %s trials (prespecified %s plus extension %s)", ex$scenario, fint(ex$n_trials), fint(ex$n_trials_prereg), fint(ex$n_trials - ex$n_trials_prereg)), collapse = ", ")) else "",
               if (s$n_incomplete_main + s$n_incomplete_ext > 0) sprintf("; %d (trial, scenario) pairs without all 6 endpoints excluded (batch being appended)", s$n_incomplete_main + s$n_incomplete_ext) else "")
}
# (1) 기준 문장
thB <- B[config %in% c("P2", "AUCinf_true_only")]                     # AUClast(P2 행), AUCinf_true(기준 행): 평가변수별 한 번
th_min <- min(thB$theory_boundary_pct); far_max <- max(thB$theory_far_prob)
# "5%에서 무시할 만큼 작은 확률을 뺀 값, 5% 미만": 반대쪽 한계 확률 > 0(수학적으로 5% 미만; 표시값은 부동소수 오차 1e-12 이내), 최대 1e-3 미만
stopifnot(all(thB$theory_far_prob > 0), all(thB$theory_boundary_pct <= ALPHA + 1e-12), all(abs(thB$theory_boundary_pct - (ALPHA - 100 * thB$theory_far_prob)) < 1e-10), far_max < 1e-3)
th_dig <- max(4, min(9, ceiling(-log10(ALPHA - th_min)) + 1))
fth <- function(x) formatC(x, format = "f", digits = th_dig)
se_al <- P2$se_median; se_tr <- RF$se_median
base_ko <- sprintf("추정량이 참값 비에 대해 불편하면 두 단측 검정(90%% CI가 80.00–125.00%% 안)의 경계 1종 오류는 명목 5%%에 가깝다(t 기반, df = 2 × %d − 2 = %d). 참값이 정확히 0.80(또는 1.25)이면 가까운 한계 쪽 확률이 정확히 5%%이고, 통과 확률은 5%%에서 반대쪽 한계에서 불통과할 확률(0.80이면 상한 > 125%%, 1.25면 하한 < 80%%)을 뺀 값이다. 저장 CI에서 역산한 시나리오별 표준오차(log 척도, 시험 중앙값: AUClast %s, AUCinf_true %s)에서 그 확률은 최대 %s(= %s%%p)이므로 이론값은 %s%% 이상, 5%% 미만이다.",
                   N_ARM, DF_NOM, rng(se_al, f4), rng(se_tr, f4), fe(far_max), fe(100 * far_max), fth(th_min))
base_en <- sprintf("If the estimator is unbiased for the true ratio, the boundary type I error of the two one-sided tests (90%% CI within 80.00 to 125.00%%) is close to the nominal 5%% (t-based, df = 2 x %d - 2 = %d). When the true ratio is exactly 0.80 (or 1.25), the probability on the near limit is exactly 5%%, and the pass probability is 5%% minus the probability of failing the far limit (upper bound above 125%% at 0.80, lower bound below 80%% at 1.25). With the per-scenario standard errors implied by the stored CIs (log scale, median over trials: AUClast %s, AUCinf_true %s), that probability is at most %s (%s percentage points), so the theoretical value is at least %s%% and below 5%%.",
                   N_ARM, DF_NOM, rng(se_al, f4, " to "), rng(se_tr, f4, " to "), fe(far_max), fe(100 * far_max), fth(th_min))
tol_ko <- sprintf("역산 허용 오차(±%s%%) 때문에 시나리오의 실제 참 AUC0-inf 비는 %s이고, 그 값에서의 이론값은 AUClast %s%%, AUCinf_true %s%%다.", format(100 * as.numeric(oc$inversion$tolerance_rel)),
                  rng(RF$auc_ratio, f4), rng(P2$theory_at_truth_pct), rng(RF$theory_at_truth_pct))
tol_en <- sprintf("Because of the inversion tolerance (+/-%s%%), the actual true AUC0-inf ratios of the scenarios are %s, where the theoretical values are %s%% for AUClast and %s%% for AUCinf_true.", format(100 * as.numeric(oc$inversion$tolerance_rel)),
                  rng(RF$auc_ratio, f4, " to "), rng(P2$theory_at_truth_pct, f2, " to "), rng(RF$theory_at_truth_pct, f2, " to "))
# (b) AUCinf_true 경험: 편향 0 포함 여부, 분류, SD/SE 비, 정규 근사
rf_inc0 <- RF[bias_lo <= 0 & bias_hi >= 0]; rf_exc0 <- RF[!(bias_lo <= 0 & bias_hi >= 0)]
n_lt1 <- sum(RF$sd_se_ratio < 1)
cnt_txt <- function(x, lang) { k <- table(factor(x$class, levels = OC_TYPE1_CLASSES)); stopifnot(sum(k) == nrow(x))
  if (lang == "ko") sprintf("보수적 %d, 명목 %d, 초과 %d", k[["conservative"]], k[["nominal"]], k[["exceeding"]]) else sprintf("conservative %d, nominal %d, exceeding %d", k[["conservative"]], k[["nominal"]], k[["exceeding"]]) }
# 관측이 정규 근사 예측보다 이항 95% 폭 밖으로 큰 시나리오: se와 log GMR의 상관이 통과를 늘리는 방향(0.80에서 음, 1.25에서 양)인지 함께 적는다
RF[, `:=`(pred_z = (pass_pct - pred_pct) / (100 * sqrt(pred_pct / 100 * (1 - pred_pct / 100) / n_trials)), cor_up = cor_se_log_gmr * sign(log(target)) > 0)]
rf_hi <- RF[pred_z > Z]; rf_lo <- RF[pred_z < -Z]
n_below <- sum(RF$pass_pct < ALPHA)
if (n_below <= nrow(RF) / 2) stop("전제 불일치: \"AUCinf_true 통과율이 5%보다 약간 낮은 쪽으로 기운다\"는 문장의 전제(과반이 5% 미만)가 결과와 다르다")
emp_ko <- sprintf("경험 기준(AUCinf_true 단독, 시험 대상자의 개인 모델 적분 AUC0-inf): 경계 통과율 %s%%(%d개 모델 × 경계 시나리오 중 %d개가 5%% 미만). 편향은 %s%%이고 %s. 분류는 %s. 5%%보다 약간 낮은 쪽으로 기우는 이유: 시험 간 log GMR 표준편차가 시험 안 표준오차(중앙값)보다 작다(비 %s, %d개 중 %d개 < 1). 즉 90%% CI가 추정량의 실제 변동보다 넓어 경계 통과가 줄고, 정규 근사(편향 0, 시험 간 SD, se 중앙값)의 예측은 %s%%다. 이 비가 1보다 작은 원인은 이 분석에서 검정하지 않았다(후보: 체중 층 안 1:1 배정; pooled t는 층을 모형에 넣지 않는다).%s",
                  rng(RF$pass_pct), nrow(RF), n_below, rng(RF$bias_pct, s2, " ~ "), if (!nrow(rf_exc0)) "모두 95% 구간이 0을 포함한다(불편과 부합)" else sprintf("%d개는 95%% 구간이 0을 제외한다(%s)", nrow(rf_exc0), paste(sprintf("%s %s %s%%", MODEL_KO[rf_exc0$model], rf_exc0$scenario, s2(rf_exc0$bias_pct)), collapse = ", ")),
                  paste(vapply(MODELS, function(mdl) sprintf("%s %s", MODEL_KO[[mdl]], cnt_txt(RF[model == mdl], "ko")), ""), collapse = "; "),
                  rng(RF$sd_se_ratio, f3), nrow(RF), n_lt1, rng(RF$pred_unbiased_pct),
                  if (nrow(rf_hi) || nrow(rf_lo)) sprintf(" 관측이 정규 근사 예측(편향 포함)과 이항 95%% 폭 밖으로 다른 곳: %s.", paste(c(
                    if (nrow(rf_hi)) sprintf("%s %s 관측 %s%% > 예측 %s%%(se와 log GMR의 상관 %s%s)", MODEL_KO[rf_hi$model], rf_hi$scenario, f2(rf_hi$pass_pct), f2(rf_hi$pred_pct), f2(rf_hi$cor_se_log_gmr),
                                             ifelse(rf_hi$cor_up, ": 통과 쪽 시험의 se가 작아 통과를 늘리는 방향", ": 통과를 늘리는 방향이 아님, 원인 미확인")),
                    if (nrow(rf_lo)) sprintf("%s %s 관측 %s%% < 예측 %s%%", MODEL_KO[rf_lo$model], rf_lo$scenario, f2(rf_lo$pass_pct), f2(rf_lo$pred_pct))), collapse = "; ")) else "")
emp_en <- sprintf("Empirical reference (AUCinf_true alone, the individual model-integrated AUC0-inf of the trial subjects): boundary pass rate %s%% (%d model by boundary scenario cells, %d of them below 5%%). The bias is %s%%, and %s. Classification: %s. Why it leans slightly below 5%%: the between-trial SD of log GMR is smaller than the within-trial standard error (median) (ratio %s; below 1 in %d of %d). The 90%% CI is therefore wider than the actual variation of the estimator, which lowers the boundary pass rate; the normal approximation (bias 0, between-trial SD, median se) predicts %s%%. The reason for a ratio below 1 was not tested in this analysis (candidate: 1:1 allocation within weight strata, which the pooled t does not model).%s",
                  rng(RF$pass_pct, f2, " to "), nrow(RF), n_below, rng(RF$bias_pct, s2, " to "), if (!nrow(rf_exc0)) "every 95% CI includes 0 (consistent with unbiasedness)" else sprintf("the 95%% CI excludes 0 in %d (%s)", nrow(rf_exc0), paste(sprintf("%s %s %s%%", MODEL_EN[rf_exc0$model], rf_exc0$scenario, s2(rf_exc0$bias_pct)), collapse = ", ")),
                  paste(vapply(MODELS, function(mdl) sprintf("%s %s", MODEL_EN[[mdl]], cnt_txt(RF[model == mdl], "en")), ""), collapse = "; "),
                  rng(RF$sd_se_ratio, f3, " to "), n_lt1, nrow(RF), rng(RF$pred_unbiased_pct, f2, " to "),
                  if (nrow(rf_hi) || nrow(rf_lo)) sprintf(" Observed and normal-approximation (with bias) values differ by more than the binomial 95%% margin in: %s.", paste(c(
                    if (nrow(rf_hi)) sprintf("%s %s observed %s%% above predicted %s%% (correlation of se with log GMR %s%s)", MODEL_EN[rf_hi$model], rf_hi$scenario, f2(rf_hi$pass_pct), f2(rf_hi$pred_pct), f2(rf_hi$cor_se_log_gmr),
                                             ifelse(rf_hi$cor_up, ", in the direction that gives passing trials a smaller se", ", not in the direction that raises passing; cause not identified")),
                    if (nrow(rf_lo)) sprintf("%s %s observed %s%% below predicted %s%%", MODEL_EN[rf_lo$model], rf_lo$scenario, f2(rf_lo$pass_pct), f2(rf_lo$pred_pct))), collapse = "; ")) else "")
stopifnot(n_lt1 == sum(RF$sd_se_ratio < 1), all(RF[pred_z > Z, pass_pct > pred_pct]))
if (n_lt1 < nrow(RF) / 2) stop("전제 불일치: AUCinf_true의 시험 간 SD / 시험 안 SE가 대부분 1 미만이라는 설명 문장의 전제가 결과와 다르다")

# (2) 분류 개수·명목/초과 사례
class_lines <- function(lang) unlist(lapply(MODELS, function(mdl) {
  x <- B[model == mdl]
  sprintf(if (lang == "ko") "  - %s, 경계 %d개: %s" else "  - %s, %d boundary scenarios: %s", if (lang == "ko") MODEL_KO[[mdl]] else MODEL_EN[[mdl]], uniqueN(x$scenario),
          paste(vapply(CFG3$config, function(cf) sprintf("%s %s", CFG3[config == cf, if (lang == "ko") label_ko else label_en], cnt_txt(x[config == cf], lang)), ""), collapse = "; "))
}))
named_cases <- function(lang) {
  x <- B[class != "conservative"]
  if (!nrow(x)) return(if (lang == "ko") "  - 없음(모든 칸이 Wilson 상한 < 5%)." else "  - None (every cell has a Wilson upper bound below 5%).")
  unlist(lapply(CFG3$config, function(cf) {
    y <- x[config == cf]; if (!nrow(y)) return(NULL)
    lab <- CFG3[config == cf, if (lang == "ko") label_ko else label_en]
    sprintf("  - %s: %s.", lab, paste(vapply(seq_len(nrow(y)), function(i) { r <- y[i]
      if (lang == "ko") sprintf("%s %s %s%% [%s, %s](%s, 시험 %s회)", MODEL_KO[[r$model]], r$scenario, f2(r$pass_pct), fb(r$lo), fb(r$hi), CLASS_KO[[r$class]], fint(r$n_trials))
      else sprintf("%s %s %s%% (95%% CI %s to %s; %s, %s trials)", MODEL_EN[[r$model]], r$scenario, f2(r$pass_pct), fb(r$lo), fb(r$hi), CLASS_EN[[r$class]], fint(r$n_trials)) }, ""), collapse = "; "))
  }))
}
# (3) P2 최종 점추정 5% 초과: 이름·방향·크기·구간·시험 수. 시험 수는 자료의 시험 수(위 검사)
EX <- P2[pass_pct > ALPHA][order(-pass_pct)]
for (i in seq_len(nrow(EX))) stopifnot(EX$n_trials[i] == nrow(A_FIN[[EX$model[i]]]$w[scenario == EX$scenario[i]]))
ex_ko <- if (nrow(EX)) sprintf("P2의 최종 점추정이 5%%를 넘는 곳: %s.", paste(vapply(seq_len(nrow(EX)), function(i) { r <- EX[i]
    sprintf("%s: %s, 시험 %s회%s%s", full_lab(r, "ko"), ci_k(r$pass_pct, r$lo, r$hi), fint(r$n_trials),
            if (r$class == "nominal") "(Wilson 구간이 5%를 포함해 명목으로 분류)" else "(Wilson 하한 > 5%: 초과)",
            if (r$n_trials > r$n_trials_prereg) sprintf("; 사전 고정 %s회 값 %s", fint(r$n_trials_prereg), ci_k(r$pass_pct_prereg, r$lo_prereg, r$hi_prereg)) else "") }, ""), collapse = "; ")) else
  sprintf("P2의 최종 점추정은 모든 경계 시나리오에서 5%% 이하다(최대 %s, %s, 시험 %s회).", ci_k(max(P2$pass_pct), P2[which.max(pass_pct), lo], P2[which.max(pass_pct), hi]), full_lab(P2[which.max(pass_pct)], "ko"), fint(P2[which.max(pass_pct), n_trials]))
ex_en <- if (nrow(EX)) sprintf("Where the final point estimate of P2 exceeds 5%%: %s.", paste(vapply(seq_len(nrow(EX)), function(i) { r <- EX[i]
    sprintf("%s: %s, %s trials%s%s", full_lab(r, "en"), ci_e(r$pass_pct, r$lo, r$hi), fint(r$n_trials),
            if (r$class == "nominal") " (the Wilson interval includes 5%, so classified nominal)" else " (Wilson lower bound above 5%: exceeding)",
            if (r$n_trials > r$n_trials_prereg) sprintf("; prespecified %s-trial value %s", fint(r$n_trials_prereg), ci_e(r$pass_pct_prereg, r$lo_prereg, r$hi_prereg)) else "") }, ""), collapse = "; ")) else
  sprintf("The final point estimate of P2 is at most 5%% in every boundary scenario (highest %s, %s, %s trials).", ci_e(max(P2$pass_pct), P2[which.max(pass_pct), lo], P2[which.max(pass_pct), hi]), full_lab(P2[which.max(pass_pct)], "en"), fint(P2[which.max(pass_pct), n_trials]))
c33_ko <- switch(C33$status,
  checked = sprintf("scripts/33 결론 첫 문단(`oc_conclusion_ko.md`%s)과 `boundary_type1.csv`의 P2·AUClast 단독 값이 이 표의 %s 값과 같다(사례, 배율, 통과율·구간을 33 문단의 표기 자릿수로, 시험 수가 적혀 있으면 시험 수도 대조; stopifnot).",
                    if (isTRUE(C33$en)) ", `oc_conclusion_en.md`" else "", if (C33$basis == "final") "최종" else sprintf("사전 고정 %s회", fint(N_BND))),
  stale = sprintf("주의: scripts/33 결과가 이 자료와 시험 수가 달라 대조하지 못했다(%s). scripts/33을 다시 실행한 뒤 이 스크립트를 다시 실행한다.", C33$note),
  absent = "scripts/33 결과(`boundary_type1.csv`)가 없어 대조하지 않았다.")
c33_en <- switch(C33$status,
  checked = sprintf("The first paragraph of the scripts/33 conclusion (`oc_conclusion_ko.md`%s) and the P2 and AUClast-only values in `boundary_type1.csv` agree with the %s values here (cases, multipliers, pass rates and intervals compared at the precision printed by scripts/33, and trial counts when printed; stopifnot).",
                    if (isTRUE(C33$en)) " and `oc_conclusion_en.md`" else "", if (C33$basis == "final") "final" else sprintf("prespecified %s-trial", fint(N_BND))),
  stale = sprintf("Caution: the scripts/33 results have different trial counts from these data, so they were not cross-checked (%s). Rerun scripts/33 and then this script.", C33$note),
  absent = "No scripts/33 results (`boundary_type1.csv`), so no cross-check.")
if (C33$status == "checked") stopifnot(isTRUE(C33$ko))
if (C33$status == "checked" && C33$basis == "prereg" && any(TAB$n_trials > TAB$n_trials_prereg)) {       # 33이 연장분을 아직 반영하지 않음
  c33_ko <- paste(c33_ko, "scripts/33은 적응적 연장분을 아직 반영하지 않았다(연장 파일 생성 뒤 scripts/33을 다시 실행하면 최종 값으로 대조된다).")
  c33_en <- paste(c33_en, "scripts/33 does not yet include the adaptive extension rows (rerun scripts/33 after the extension file is written to cross-check the final values).")
}
# (4) 원인 요약: 편향 방향과 AUClast 효과 부호의 일치
dir_counts <- function(x, lang) {
  sig <- x[auclast_minus_ref_lo > 0 | auclast_minus_ref_hi < 0]
  k_aw <- sig[bias_dir == "away_from_1"]; k_to <- sig[bias_dir == "toward_1"]
  stopifnot(all(sig$auclast_bias_consistent == ((sig$bias_dir == "away_from_1" & sig$auclast_minus_ref_pp < 0) | (sig$bias_dir == "toward_1" & sig$auclast_minus_ref_pp > 0))))
  n_cons <- sum(sig$auclast_bias_consistent, na.rm = TRUE)
  aw <- x[bias_dir == "away_from_1"]; to <- x[bias_dir == "toward_1"]; no <- x[bias_dir == "none"]
  awr <- function(sep) paste(vapply(sort(unique(aw$target)), function(tg) { y <- aw[abs(target - tg) < 1e-9]
    sprintf(if (sep == "ko") "참값 %.2f에서 %d개 %s%%" else "%d at true %.2f, %s%%", if (sep == "ko") tg else nrow(y), if (sep == "ko") nrow(y) else tg, rng(y$bias_pct, s2, if (sep == "ko") " ~ " else " to ")) }, ""), collapse = ", ")
  if (lang == "ko") sprintf("AUClast 편향(참 AUC0-inf 비 대비): 1에서 멀어지는 방향 %d개(%s), 1 쪽 %d개%s, 0 포함 %d개. AUClast 단독 − AUCinf_true가 0과 구별되는 %d개 중 %d개에서 부호가 편향 방향과 맞는다(멀어지면 음, 1 쪽이면 양).",
                            nrow(aw), if (nrow(aw)) awr("ko") else "-", nrow(to), if (nrow(to)) sprintf("(%s)", paste(sprintf("%s %s %s%%", MODEL_KO[to$model], to$scenario, s2(to$bias_pct)), collapse = ", ")) else "", nrow(no), nrow(sig), n_cons)
  else sprintf("AUClast bias (versus the true AUC0-inf ratio): away from 1 in %d (%s), toward 1 in %d%s, CI including 0 in %d. Of the %d cells where AUClast only minus AUCinf_true is distinguishable from 0, the sign agrees with the bias direction in %d (negative when away from 1, positive when toward 1).",
               nrow(aw), if (nrow(aw)) awr("en") else "-", nrow(to), if (nrow(to)) sprintf(" (%s)", paste(sprintf("%s %s %s%%", MODEL_EN[to$model], to$scenario, s2(to$bias_pct)), collapse = ", ")) else "", nrow(no), nrow(sig), n_cons)
}

# 원인 요약(P2): 각 칸의 가장 큰 요소별 칸 수, Cmax가 대부분인 칸, 보수적이 아닌 칸의 한 줄 원인
cause_summary <- function(lang) {
  k <- P2[, .N, by = "cause_primary"][order(-N)]
  lab <- if (lang == "ko") c(last = "AUClast 효과(AUClast 단독 − AUCinf_true)", cmax = "Cmax 불통과", ref = "편향 없는 기준과 5%의 차이", none = "유의한 요소 없음")
         else c(last = "the AUClast effect (AUClast only minus AUCinf_true)", cmax = "Cmax failures", ref = "the gap between the unbiased reference and 5%", none = "no component distinguishable from 0")
  cd <- P2[cmax_dominant %in% TRUE]; nc <- P2[class != "conservative"]
  stopifnot(sum(k$N) == nrow(P2), all(cd$cause_primary == "cmax"))
  if (lang == "ko") c(sprintf("- 원인(P2, 칸마다 가장 큰 요소): %s.%s", paste(sprintf("%s %d칸", lab[k$cause_primary], k$N), collapse = ", "),
                              if (nrow(cd)) sprintf(" Cmax 불통과가 대부분인 칸: %s.", paste(sprintf("%s %s(Cmax 통과 %s%%, P2 %s%%)", MODEL_KO[cd$model], cd$scenario, f2(cd$cmax_pass_pct), f2(cd$pass_pct)), collapse = ", ")) else ""),
                      if (nrow(nc)) sprintf("  - 보수적이 아닌 P2 칸: %s", sprintf("%s %s: %s", MODEL_KO[nc$model], nc$scenario, nc$cause_ko)))
  else c(sprintf("- Causes (P2, largest component per cell): %s.%s", paste(sprintf("%s in %d", lab[k$cause_primary], k$N), collapse = ", "),
                 if (nrow(cd)) sprintf(" Cells where Cmax failures dominate: %s.", paste(sprintf("%s %s (Cmax pass %s%%, P2 %s%%)", MODEL_EN[cd$model], cd$scenario, f2(cd$cmax_pass_pct), f2(cd$pass_pct)), collapse = ", ")) else ""),
         if (nrow(nc)) sprintf("  - P2 cells that are not conservative: %s", sprintf("%s %s: %s", MODEL_EN[nc$model], nc$scenario, nc$cause_en)))
}

# 표
cell <- function(p, l, h, cls, lang) sprintf(if (lang == "ko") "%s [%s, %s] %s" else "%s [%s to %s] %s", f2(p), fb(l), fb(h), if (lang == "ko") CLASS_KO[cls] else CLASS_EN[cls])
class_table <- function(mdl, lang) {
  x <- B[model == mdl]; sc <- unique(x$scenario)
  c(sprintf("| %s | %s |", if (lang == "ko") "경계 시나리오" else "Boundary scenario", paste(if (lang == "ko") CFG3$label_ko else CFG3$label_en, collapse = " | ")), paste0("|", strrep("---|", nrow(CFG3) + 1)),
    vapply(sc, function(s) { y <- x[scenario == s][match(CFG3$config, config)]
      paste0("| ", scn_lab(y[1], lang), " | ", paste(cell(y$pass_pct, y$lo, y$hi, y$class, lang), collapse = " | "), " |") }, ""))
}
base_table <- function(mdl, lang) {
  x <- B[model == mdl]; sc <- unique(x$scenario)
  h <- if (lang == "ko") c("경계 시나리오", "참 AUC0-inf 비", "se AUClast / AUCinf_true", "반대쪽 한계 확률 AUClast / AUCinf_true", "이론(실제 참값) AUClast / AUCinf_true %", "AUCinf_true 통과율 % [Wilson]", "AUCinf_true 편향 % (MC SE)", "SD / SE", "정규 근사 %")
       else c("Boundary scenario", "True AUC0-inf ratio", "se AUClast / AUCinf_true", "Far-limit probability AUClast / AUCinf_true", "Theory (actual truth) AUClast / AUCinf_true %", "AUCinf_true pass % [Wilson]", "AUCinf_true bias % (MC SE)", "SD / SE", "Normal approx. %")
  c(paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", strrep("---|", length(h))),
    vapply(sc, function(s) { a <- x[scenario == s & config == "P2"]; r <- x[scenario == s & config == "AUCinf_true_only"]
      paste0("| ", paste(c(scn_lab(r, lang), f4(r$auc_ratio), sprintf("%s / %s", f4(a$se_median), f4(r$se_median)), sprintf("%s / %s", fe(a$theory_far_prob), fe(r$theory_far_prob)),
                           sprintf("%s / %s", f2(a$theory_at_truth_pct), f2(r$theory_at_truth_pct)), cell(r$pass_pct, r$lo, r$hi, r$class, lang),
                           sprintf("%s (%s)", s2(r$bias_pct), f2(100 * r$bias_mc_se)), f3(r$sd_se_ratio), f2(r$pred_pct)), collapse = " | "), " |") }, ""))
}
cause_table <- function(mdl, lang) {
  x <- B[model == mdl]; sc <- unique(x$scenario)
  h <- if (lang == "ko") c("경계 시나리오", "P2 %", "P2 − 5", "= 기준 − 5", "+ AUClast 단독 − AUCinf_true [95%]", "+ P2 − AUClast 단독 [95%]", "Cmax 통과 %", "AUClast 편향 % (MC SE)", "Cmax 편향 % (MC SE)")
       else c("Boundary scenario", "P2 %", "P2 minus 5", "= reference minus 5", "+ AUClast only minus AUCinf_true [95%]", "+ P2 minus AUClast only [95%]", "Cmax pass %", "AUClast bias % (MC SE)", "Cmax bias % (MC SE)")
  dc <- function(e, l, u) sprintf(if (lang == "ko") "%s [%s, %s]" else "%s [%s to %s]", s2(e), s2(l), s2(u))
  c(paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", strrep("---|", length(h))),
    vapply(sc, function(s) { r <- x[scenario == s & config == "P2"]
      paste0("| ", paste(c(scn_lab(r, lang), f2(r$pass_pct), s2(r$gap_pp), s2(r$ref_minus_5_pp), dc(r$auclast_minus_ref_pp, r$auclast_minus_ref_lo, r$auclast_minus_ref_hi),
                           dc(r$p2_minus_auclast_pp, r$p2_minus_auclast_lo, r$p2_minus_auclast_hi), f2(r$cmax_pass_pct), sprintf("%s (%s)", s2(r$bias_pct), f2(100 * r$bias_mc_se)),
                           sprintf("%s (%s)", s2(r$cmax_bias_pct), f2(100 * r$cmax_bias_mc_se))), collapse = " | "), " |") }, ""))
}
bias_table <- function(mdl, lang) {
  x <- B[model == mdl & config == "P2"]
  h <- c(if (lang == "ko") "경계 시나리오" else "Boundary scenario", OC_ENDPOINTS, if (lang == "ko") "scripts/41 대조" else "scripts/41 check")
  c(paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", strrep("---|", length(h))),
    vapply(seq_len(nrow(x)), function(i) { r <- x[i]
      paste0("| ", paste(c(scn_lab(r, lang), vapply(OC_ENDPOINTS, function(e) sprintf("%s (%s)", s2(r[[paste0("bias_pct_", e)]]), f2(r[[paste0("bias_mc_se_pct_", e)]])), ""),
                           b41_txt(r$bias41_status, r$bias41_n_trials, lang)), collapse = " | "), " |") }, ""))
}
b41_sum <- unique(TAB[, .(model, bias41_status, bias41_n_trials)])
b41_ko <- if (all(TAB$bias41_status == "equal")) "`p2_bias_boundary.csv`(scripts/41)의 원래 6개 평가변수 값과 같다(1e-9, 같은 시험 수)." else
  sprintf("`p2_bias_boundary.csv`(scripts/41) 대조: %s. 시험 수가 같은 칸은 값이 같고(1e-9), 시험 수가 다른 칸은 41이 그보다 앞선(부분) 파일로 실행된 것이다(scripts/41을 다시 실행하면 이 표와 같아진다).",
          paste(sprintf("%s %s", MODEL_KO[b41_sum$model], mapply(b41_txt, b41_sum$bias41_status, b41_sum$bias41_n_trials, MoreArgs = list(lang = "ko"))), collapse = "; "))
b41_en <- if (all(TAB$bias41_status == "equal")) "Equal to the 6 original endpoints in `p2_bias_boundary.csv` (scripts/41) to 1e-9, with the same trial counts." else
  sprintf("Check against `p2_bias_boundary.csv` (scripts/41): %s. Cells with the same trial count are equal (1e-9); a different trial count means scripts/41 ran on an earlier (partial) file, and rerunning scripts/41 makes it equal to this table.",
          paste(sprintf("%s %s", MODEL_EN[b41_sum$model], mapply(b41_txt, b41_sum$bias41_status, b41_sum$bias41_n_trials, MoreArgs = list(lang = "en"))), collapse = "; "))
if (any(TAB$bias41_status == "stale_differs")) warning("p2_bias_boundary.csv(scripts/41)가 같은 앞 시험에서도 이 표와 다릅니다")

ko <- c("# P2 경계 1종 오류의 해석 (지시 2026-09-25 §3)", "",
  "자료: 저장된 시험 수준 결과(`results/oc/oc_trials_be_<모델>.csv.gz`, scripts/31: arm당 117명, B0, pooled t, 90% CI 80.00–125.00%)와 있으면 적응적 연장분(`oc_trials_ext_be_<모델>.csv.gz`, scripts/42). 새 모의 없음. 생성: `scripts/43_p2_interpretation.R`. 수치: `p2_interpretation.csv`(모델 × 경계 시나리오 × 구성).", "",
  paste0("- 시험 수: ", vapply(MODELS, trials_txt, "", lang = "ko"), "."), "",
  "## 요약", "",
  paste0("- 기준(이론): ", base_ko, " ", tol_ko),
  paste0("- ", emp_ko),
  "- 분류(Wilson 95% 구간 대 5%: 상한 < 5% 보수적, 하한 > 5% 초과, 그 밖 명목):", class_lines("ko"),
  "- 명목·초과 사례(크기, Wilson 95% 구간):", named_cases("ko"), cause_summary("ko"),
  paste0("- ", ex_ko, " ", c33_ko), "",
  "## 1. 기준: 불편 추정량의 경계 1종 오류", "",
  "- 이론값은 R/oc_interpret.R `be_boundary_type1_theory`(t 기반, se 고정). 참값이 정확히 0.80 또는 1.25이면 이론값 = 5% − 반대쪽 한계 확률(표의 확률, 0–1 척도), 이론(실제 참값) = 시나리오의 역산 참값에서. 정규 근사 = 시험 간 log GMR의 평균·SD와 se 중앙값으로 계산한 AUCinf_true 통과율(`be_pass_prob_normal`). SD / SE = 시험 간 log GMR 표준편차 / 시험 안 se 중앙값.", "")
for (mdl in MODELS) ko <- c(ko, sprintf("### %s", MODEL_KO[[mdl]]), "", base_table(mdl, "ko"), "")
ko <- c(ko, "## 2. 분류 (통과율 %, Wilson 95% 구간)", "")
for (mdl in MODELS) ko <- c(ko, sprintf("### %s", MODEL_KO[[mdl]]), "", class_table(mdl, "ko"), "")
ko <- c(ko, "## 3. 원인: 편향과 Cmax", "",
  "- 방향: 추정량이 참값 비보다 1에서 멀어지는 쪽(0.80에서 음, 1.25에서 양)으로 치우치면 CI가 범위 밖으로 밀려 경계 1종 오류가 낮아지고, 1 쪽으로 치우치면 높아진다.",
  "- 항등식(같은 시험, 정확히 성립; stopifnot): P2 − 5 = (AUCinf_true − 5) + (AUClast 단독 − AUCinf_true) + (P2 − AUClast 단독). 첫 항은 편향 없는 기준과 5%의 차이, 둘째 항은 AUClast(편향 + 변동)의 효과, 셋째 항은 Cmax 불통과의 효과(항상 ≤ 0).",
  paste0("- ", dir_counts(P2, "ko")), "")
for (mdl in MODELS) {
  ko <- c(ko, sprintf("### %s", MODEL_KO[[mdl]]), "", cause_table(mdl, "ko"), "", "한 줄 원인(P2; 95% 구간이 0 또는 5%를 제외하는 요소만, 크기 순):", "",
          sprintf("- %s: %s", P2[model == mdl, scenario], P2[model == mdl, cause_ko]), "", "AUClast 단독:", "", sprintf("- %s: %s", AL[model == mdl, scenario], AL[model == mdl, cause_ko]), "")
}
ko <- c(ko, "### 평가변수별 편향 (%, MC 표준오차; 참 AUC0-inf 비 대비, Cmax는 참 Cmax 비 대비)", "", paste0("- ", b41_ko), "")
for (mdl in MODELS) ko <- c(ko, sprintf("#### %s", MODEL_KO[[mdl]]), "", bias_table(mdl, "ko"), "")
ko <- c(ko, "## 전제 검사(stopifnot)", "",
  "- 고정 문구의 설계값: 90% CI, 한계 0.80–1.25, arm당 117명(df 232), 명목 5%, pooled t. AUClast·AUCinf_true의 시험별 df 중앙값 = 232.",
  "- P2·AUClast 단독 = config_pass(사전 고정 정의), P2 통과 ⇒ AUClast·Cmax 통과. 분류 = Wilson 구간에서 다시 계산한 classify_type1. 보고한 시험 수 = 자료의 시나리오별 시험 수.",
  "- 이론: 경계에서 가까운 한계 확률 = 5%, 통과 = 5% − 반대쪽 한계 확률. 항등식 P2 − 5 = 세 항의 합. 편향 방향 문구 = 부호(0.80에서 음 = 1에서 멀어짐), \"편향 때문\" 문구는 AUClast 효과의 부호가 편향 방향과 맞을 때만, \"Cmax 불통과가 대부분\"은 Cmax 효과가 차이의 절반 이상일 때만.",
  "- scripts/41 p2_bias_boundary.csv(같은 시험 수면 1e-9), scripts/33 boundary_type1.csv·결론 첫 문단(같은 시험 기준으로 정확히)과 대조. 영문 파일에 한글·em dash·en dash·U+2212 없음.")

en <- c("# Interpretation of the P2 boundary type I error (directive 2026-09-25 section 3)", "",
  "Data: saved trial-level results (`results/oc/oc_trials_be_<model>.csv.gz`, scripts/31: 117 per arm, B0, pooled t, 90% CI within 80.00 to 125.00%) plus the adaptive extension rows when present (`oc_trials_ext_be_<model>.csv.gz`, scripts/42). No new simulation. Generated by `scripts/43_p2_interpretation.R`. Numbers: `p2_interpretation.csv` (model by boundary scenario by configuration).", "",
  paste0("- Trials: ", vapply(MODELS, trials_txt, "", lang = "en"), "."), "",
  "## Summary", "",
  paste0("- Baseline (theory): ", base_en, " ", tol_en),
  paste0("- ", emp_en),
  "- Classification (Wilson 95% CI versus 5%: upper bound below 5% conservative, lower bound above 5% exceeding, otherwise nominal):", class_lines("en"),
  "- Nominal and exceeding cells (size, Wilson 95% CI):", named_cases("en"), cause_summary("en"),
  paste0("- ", ex_en, " ", c33_en), "",
  "## 1. Baseline: boundary type I error of an unbiased estimator", "",
  "- Theory from R/oc_interpret.R `be_boundary_type1_theory` (t-based, se fixed). At a true ratio of exactly 0.80 or 1.25 the theoretical value is 5% minus the far-limit probability (shown in the table on the 0 to 1 scale); theory (actual truth) = at the inverted true ratio of the scenario. Normal approximation = AUCinf_true pass rate from the between-trial mean and SD of log GMR and the median se (`be_pass_prob_normal`). SD / SE = between-trial SD of log GMR divided by the median within-trial se.", "")
for (mdl in MODELS) en <- c(en, sprintf("### %s", MODEL_EN[[mdl]]), "", base_table(mdl, "en"), "")
en <- c(en, "## 2. Classification (pass rate %, Wilson 95% CI)", "")
for (mdl in MODELS) en <- c(en, sprintf("### %s", MODEL_EN[[mdl]]), "", class_table(mdl, "en"), "")
en <- c(en, "## 3. Causes: bias and Cmax", "",
  "- Direction: an estimator biased away from 1 relative to the true ratio (negative at 0.80, positive at 1.25) pushes the CI out of the limits and lowers the boundary type I error; a bias toward 1 raises it.",
  "- Identity (same trials, exact; stopifnot): P2 minus 5 = (AUCinf_true minus 5) + (AUClast only minus AUCinf_true) + (P2 minus AUClast only). The first term is the gap between the unbiased reference and 5%, the second the AUClast effect (bias and variation), the third the effect of Cmax failures (never positive).",
  paste0("- ", dir_counts(P2, "en")), "")
for (mdl in MODELS) {
  en <- c(en, sprintf("### %s", MODEL_EN[[mdl]]), "", cause_table(mdl, "en"), "", "One-line cause (P2; only components whose 95% CI excludes 0 or 5%, largest first):", "",
          sprintf("- %s: %s", P2[model == mdl, scenario], P2[model == mdl, cause_en]), "", "AUClast only:", "", sprintf("- %s: %s", AL[model == mdl, scenario], AL[model == mdl, cause_en]), "")
}
en <- c(en, "### Bias by endpoint (%, MC SE; versus the true AUC0-inf ratio, Cmax versus the true Cmax ratio)", "", paste0("- ", b41_en), "")
for (mdl in MODELS) en <- c(en, sprintf("#### %s", MODEL_EN[[mdl]]), "", bias_table(mdl, "en"), "")
en <- c(en, "## Premise checks (stopifnot)", "",
  "- Design values behind the fixed wording: 90% CI, limits 0.80 to 1.25, 117 per arm (df 232), nominal 5%, pooled t. The median per-trial df of AUClast and AUCinf_true is 232.",
  "- P2 and AUClast only come from config_pass (prespecified definitions), and a P2 pass implies AUClast and Cmax passes. Classes are recomputed from the Wilson intervals (classify_type1). Reported trial counts equal the per-scenario trial counts of the data.",
  "- Theory: at the boundary the near-limit probability is 5% and the pass probability is 5% minus the far-limit probability. The identity P2 minus 5 = sum of the three terms holds. Bias direction wording follows the sign (negative at 0.80 = away from 1); the \"bias\" wording is used only when the sign of the AUClast effect agrees with the bias direction, and \"Cmax failures dominate\" only when the Cmax effect is at least half of the gap.",
  "- Cross-checked against scripts/41 p2_bias_boundary.csv (1e-9 when the trial counts match) and against scripts/33 boundary_type1.csv and the first conclusion paragraph (exactly, on the same trials). No Korean, em dash, en dash or U+2212 in the English file.")
# 영문: 한글·em dash·en dash·U+2212 없음
stopifnot(!any(grepl("[ᄀ-ᇿ㄰-㆏가-힣]", en)), !any(grepl("[—–−]", en)))
fwrite(TAB, file.path(out_dir, "p2_interpretation.csv"))
writeLines(ko, file.path(out_dir, "p2_interpretation_ko.md")); writeLines(en, file.path(out_dir, "p2_interpretation_en.md"))
print(ST)
print(TAB[, .(model, scenario, config, n_trials, pass_pct = round(pass_pct, 2), lo = round(lo, 2), hi = round(hi, 2), class, bias_pct = round(bias_pct, 3), bias_dir)])
cat("theory at boundary:", fth(th_min), "-", fth(max(thB$theory_boundary_pct)), " far max", fe(far_max), "\n")
cat("scripts/33 check:", C33$status, C33$basis, "\n")
cat("written:", file.path(out_dir, c("p2_interpretation.csv", "p2_interpretation_ko.md", "p2_interpretation_en.md")), sep = "\n  ")
