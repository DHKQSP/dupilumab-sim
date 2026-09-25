#!/usr/bin/env Rscript
# 41_oc_rules_flags.R — G2(AUCinf + Cmax)와 AUCinf 단독의 경계 1종 오류: AUCinf 처리 규칙 A·B·C × 플래그 세트 (i)·(ii) (검토 의견 W2 §2, D-040).
# 새 모의 없음. 입력: results/oc/oc_trials_be_<model>.csv.gz(scripts/31, 원래 6개 평가변수), 있으면 <rejudge_dir>/oc_rejudge_be_<model>.csv.gz
# (scripts/40, 같은 시드 재생성: 세트 (i) 평가변수 AUCinf_Ai·AUCinf_Ci). 31이 아직 돌고 있으면 부분 파일의 완전한 (시험, 시나리오)만 쓰고 시험 수를 적는다.
#  규칙 A: 플래그 대상 제외 / B: λz 산출 가능 전원 포함(플래그 무관 → (i) = (ii)) / C: 플래그 대상은 AUClast 대입
#  세트 (i): λz 산출 가능 & Rsq_adjusted ≥ 0.80 & AUC_%Extrap_obs ≤ 20% / 세트 (ii): (i) & Span_ratio ≥ 2 (사전 고정 구성의 AUCinf_A·AUCinf_C)
# 산출(out_dir, 기본 results/oc):
#  g2_rules_flags.csv          모델 × 시나리오(경계 + S00) × 구성: 통과율(Wilson 95% 구간), 5% 초과(점추정, Wilson 하한)
#  g2_decomposition.csv        같은 시험 쌍대 차이(%p, 95% 구간): G2-B − G2-A(세트별, 탈락·선택 효과), G2-B − P2(외삽 효과), G2-C − P2, G2-A − P2(합), 세트 (i) − (ii)
#  p2_bias_boundary.csv        경계 시나리오 × 평가변수: mean(log GMR) − log(참 AUC0-inf 비)(Cmax는 참 Cmax 비), MC 표준오차, 시험 수
#  g2_rules_conclusion_ko.md / g2_rules_conclusion_en.md   수치로 생성, 고정 문구의 전제는 stopifnot으로 검사(어긋나면 생성 실패)
# 적응적 연장(scripts/42): <ext_dir>/oc_trials_ext_be_<model>.csv.gz가 있으면 연장 시험(평가변수 8개: 세트 (i) 포함)을 붙인다(규칙 × 플래그 표, 분해, 치우침 모두).
#  <ext_dir>/extension_decision.csv의 선택 시나리오만, 사전 고정 시험 1..n_before 뒤에 빈틈 없이 이어져야 한다. 없으면 사전 고정 시험만 쓴다.
# 사용법: nice -n 15 Rscript scripts/41_oc_rules_flags.R [out_dir=results/oc] [rejudge_dir=results/oc] [ext_dir=results/oc]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE)
oc <- read_cfg("oc_design.yaml")
in_dir <- proj_path("results", "oc")
out_dir <- if (length(args) >= 1) args[1] else in_dir
rj_dir <- if (length(args) >= 2) args[2] else in_dir
ext_dir <- if (length(args) >= 3) args[3] else in_dir
dec_f <- file.path(ext_dir, "extension_decision.csv"); DEC <- if (file.exists(dec_f)) fread(dec_f) else NULL
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
bnd <- as.numeric(unlist(oc$boundary_targets)); N_BND <- as.integer(oc$trials$reps_boundary)
MODELS <- intersect(c("k2016", "k2020"), sub("^oc_trials_be_(.*)\\.csv\\.gz$", "\\1", list.files(in_dir, pattern = "^oc_trials_be_.*\\.csv\\.gz$")))
stopifnot(length(MODELS) > 0)
MODEL_KO <- c(k2016 = "Kovalenko 2016 (주)", k2020 = "Kovalenko 2020 Model 1"); MODEL_EN <- c(k2016 = "Kovalenko 2016 (primary)", k2020 = "Kovalenko 2020 Model 1")
EP_I <- c("AUCinf_Ai", "AUCinf_Ci")
LIM <- 5                                                   # 명목 유의수준(%)
# 구성 메타데이터(정의는 R/oc.R OC_REJUDGE_CONFIGS)
CFG <- data.table(config = names(OC_REJUDGE_CONFIGS),
                  family = c("P2", rep("G2", 5), rep("AUCinf_only", 5)),
                  rule = c("AUClast", "A", "B", "C", "A", "C", "A", "B", "C", "A", "C"),
                  flag_set = c("-", "(ii)", "(i)=(ii)", "(ii)", "(i)", "(i)", "(ii)", "(i)=(ii)", "(ii)", "(i)", "(i)"),
                  label = c("P2", "G2-A(ii)", "G2-B", "G2-C(ii)", "G2-A(i)", "G2-C(i)", "AUCinf-A(ii)", "AUCinf-B", "AUCinf-C(ii)", "AUCinf-A(i)", "AUCinf-C(i)"))
CFG[, needs_i := vapply(OC_REJUDGE_CONFIGS[config], function(e) any(e %in% EP_I), logical(1))]
stopifnot(identical(CFG$config, names(OC_REJUDGE_CONFIGS)))
f2 <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2); s2 <- function(x) sprintf("%+.2f", round(x, 2) + 0)
fint <- function(x) format(x, big.mark = ",")

# 시나리오 정보: 31과 같은 규칙으로 역산 결과에서 만든다. oc_scenarios_<model>.csv가 있으면 같은지 확인
scen_info <- function(mdl) {
  inv <- rbindlist(lapply(list.files(in_dir, pattern = sprintf("^inversion_%s_.*\\.csv$", mdl), full.names = TRUE), fread))
  stopifnot(nrow(inv) > 0, all(inv[reachable == TRUE, within_tol]))
  inv <- inv[reachable == TRUE]
  inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
  f <- file.path(in_dir, sprintf("oc_scenarios_%s.csv", mdl))
  if (file.exists(f)) { s0 <- fread(f); m <- merge(inv[, .(code, multiplier, auc_ratio)], s0[, .(code, m0 = multiplier, a0 = auc_ratio)], by = "code", all = TRUE)
    stopifnot(nrow(m) == nrow(inv), !anyNA(m), all(abs(m$multiplier / m$m0 - 1) < 1e-12), all(abs(m$auc_ratio / m$a0 - 1) < 1e-12)) }
  b <- inv[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9]
  rbind(b[, .(scenario = code, region = "boundary", mechanism, direction, target, multiplier, auc_ratio, auc_se_log, cmax_ratio, cmax_se_log)],
        data.table(scenario = "S00", region = "identical", mechanism = "-", direction = "-", target = 1, multiplier = 1, auc_ratio = 1, auc_se_log = 0, cmax_ratio = 1, cmax_se_log = 0))
}
# 완전한 (시험, 시나리오)만: 평가변수 집합이 모두 있는 행
complete_rows <- function(be, eps) {
  be <- unique(be, by = c("trial", "scenario", "endpoint"))
  k <- be[endpoint %in% eps, .N, by = .(trial, scenario)][N == length(eps)]
  be[k[, .(trial, scenario)], on = c("trial", "scenario")]
}

rates <- list(); decomp <- list(); bias <- list(); status <- list(); W <- list(); NT <- list(); EXTN <- list()
for (mdl in MODELS) {
  info <- scen_info(mdl); SCN <- info$scenario
  be <- fread(file.path(in_dir, sprintf("oc_trials_be_%s.csv.gz", mdl)))[scenario %in% SCN]
  be <- complete_rows(be[endpoint %in% OC_ENDPOINTS], OC_ENDPOINTS)
  stopifnot(all(info$scenario %in% be$scenario))
  # 재판정 파일(세트 (i)): 원래 6개 평가변수가 저장본과 같은지 확인한 뒤 AUCinf_Ai·AUCinf_Ci만 붙인다
  rj_f <- file.path(rj_dir, sprintf("oc_rejudge_be_%s.csv.gz", mdl)); has_rj <- file.exists(rj_f); n_rj <- 0L
  if (has_rj) {
    rj <- complete_rows(fread(rj_f)[scenario %in% SCN], OC_ENDPOINTS_EXT)
    m <- merge(rj[endpoint %in% OC_ENDPOINTS], be, by = c("trial", "scenario", "endpoint"), suffixes = c("", ".s"))
    stopifnot(nrow(m) == nrow(rj[endpoint %in% OC_ENDPOINTS]),                                   # 재판정 시험은 모두 저장본에 있다
              identical(m$pass, m$pass.s), all(m$n_R == m$n_R.s), all(m$n_T == m$n_T.s),
              all((is.na(m$GMR) & is.na(m$GMR.s)) | abs(m$GMR / m$GMR.s - 1) <= 1e-6))
    # 세트 (i) ⊇ 세트 (ii): 규칙 A 포함 수는 (i)에서 줄지 않는다
    nn <- dcast(rj[endpoint %in% c("AUCinf_A", "AUCinf_Ai")], trial + scenario ~ endpoint, value.var = c("n_R", "n_T"))
    stopifnot(all(nn$n_R_AUCinf_Ai >= nn$n_R_AUCinf_A), all(nn$n_T_AUCinf_Ai >= nn$n_T_AUCinf_A))
    be <- rbind(be, rj[endpoint %in% EP_I]); n_rj <- uniqueN(rj$trial)
  }
  n_pre <- be[endpoint == "Cmax", .(N = uniqueN(trial)), by = scenario]   # 사전 고정 파일의 시나리오별 시험 수(연장 전)
  # 적응적 연장(scripts/42): 평가변수 8개 전부. 선택된 시나리오만, 사전 고정 시험 1..n_before 뒤에 이어지는 완전한 시험만 허용
  ext_f <- file.path(ext_dir, sprintf("oc_trials_ext_be_%s.csv.gz", mdl)); has_ext <- file.exists(ext_f); ext_n <- data.table(scenario = character(0), n_ext = integer(0))
  if (has_ext) {
    if (is.null(DEC)) stop("연장 파일이 있는데 extension_decision.csv가 없습니다: ", ext_f)
    ext <- fread(ext_f)
    if (!identical(names(ext), c("trial", "scenario", "endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T"))) stop("연장 파일의 열이 다릅니다: ", ext_f)
    sel_m <- DEC[model == mdl & selected == TRUE]
    if (!all(unique(ext$scenario) %in% intersect(sel_m$scenario, SCN))) stop(sprintf("%s: 연장 파일에 규칙이 고르지 않은 시나리오가 있습니다: %s", mdl, paste(setdiff(unique(ext$scenario), sel_m$scenario), collapse = ",")))
    for (sc_ in unique(ext$scenario)) {
      nb <- sel_m[scenario == sc_, n_before]
      et <- ext[scenario == sc_, .(n = .N, n_ep = uniqueN(endpoint), all8 = all(OC_ENDPOINTS_EXT %in% endpoint)), by = trial]
      if (n_pre[scenario == sc_, N] != nb || max(be[scenario == sc_, trial]) != nb || !identical(sort(et$trial), nb + seq_len(nrow(et))) ||
          any(et$n != length(OC_ENDPOINTS_EXT)) || !all(et$all8))
        stop(sprintf("%s %s: 연장 시험(%d개)이 사전 고정 시험 1-%d 뒤에 평가변수 %d개로 완전하게 이어지지 않습니다", mdl, sc_, nrow(et), nb, length(OC_ENDPOINTS_EXT)))
      ext_n <- rbind(ext_n, data.table(scenario = sc_, n_ext = nrow(et)))
    }
    nn <- dcast(ext[endpoint %in% c("AUCinf_A", "AUCinf_Ai")], trial + scenario ~ endpoint, value.var = c("n_R", "n_T"))   # 세트 (i) ⊇ 세트 (ii)
    stopifnot(all(nn$n_R_AUCinf_Ai >= nn$n_R_AUCinf_A), all(nn$n_T_AUCinf_Ai >= nn$n_T_AUCinf_A))
    be <- rbind(be, ext); rm(ext)
  }
  w <- config_pass_ext(be)
  # 사전 고정 구성 정의(config_pass)와 같은지: P2, G2 = G2-A(ii), AUCinf_only = AUCinf-A(ii)
  wc <- config_pass(be[endpoint %in% OC_ENDPOINTS], oc)
  mm <- merge(w[, .(trial, scenario, cfg_P2, cfg_G2_Aii, cfg_AUCinf_Aii)], wc[, .(trial, scenario, P2 = cfg_P2, G2 = cfg_G2, AO = cfg_AUCinf_only)], by = c("trial", "scenario"))
  stopifnot(nrow(mm) == nrow(w), identical(mm$cfg_P2, mm$P2), identical(mm$cfg_G2_Aii, mm$G2), identical(mm$cfg_AUCinf_Aii, mm$AO))
  w[, Cmax_pass := Cmax == 1L]
  # 구성 정의상 G2·P2 통과 ⇒ Cmax 통과(시험 단위, 어떤 시험 부분집합에서도 G2·P2 통과율 ≤ Cmax 통과율)
  for (cf in CFG[family %in% c("P2", "G2"), config]) stopifnot(!any(w[[paste0("cfg_", cf)]] %in% TRUE & !w$Cmax_pass))
  W[[mdl]] <- w
  n_tr <- w[, .N, by = scenario]
  status[[mdl]] <- data.table(model = mdl, n_trials_min = min(n_pre$N), n_trials_max = max(n_pre$N), complete = min(n_pre$N) >= N_BND,
                              rejudge = has_rj, n_rejudge_trials = n_rj, rejudge_complete = has_rj && n_rj >= max(n_pre$N),
                              extension = has_ext, ext_scenarios = paste(ext_n$scenario, collapse = ","), n_ext_trials = paste(ext_n$n_ext, collapse = ","))
  EXTN[[mdl]] <- ext_n[, model := mdl]
  NT[[mdl]] <- merge(info[, .(scenario, mechanism, direction, target)], n_tr[, .(scenario, n_trials = N)], by = "scenario")[, model := mdl]

  # (1) 통과율
  rt <- rbindlist(lapply(seq_len(nrow(CFG)), function(i) {
    col <- paste0("cfg_", CFG$config[i])
    w[, { x <- get(col); x <- x[!is.na(x)]; k <- length(x)
          wc_ <- if (k) wilson_ci(sum(x), k) else data.table(est = NA_real_, lo = NA_real_, hi = NA_real_)
          .(n_trials = k, n_pass = sum(x), pass_pct = wc_$est, lo = wc_$lo, hi = wc_$hi) }, by = scenario][, config := CFG$config[i]]
  }))
  rt <- merge(CFG[, .(config, family, rule, flag_set, label, needs_i)], rt, by = "config")
  rt[, source := fifelse(needs_i, "oc_rejudge_be", "oc_trials_be")]
  rt[scenario %in% ext_n$scenario, source := paste0(source, "+oc_trials_ext_be")]   # 적응적 연장 시나리오: 연장 파일의 시험 포함
  rt[, note := fifelse(needs_i & !has_rj & n_trials == 0, "rejudge file absent: run scripts/40_oc_rejudge.R",
                       fifelse(needs_i & n_trials < n_tr[match(rt$scenario, scenario), N], if (has_rj) "rejudge covers fewer trials than stored" else "rejudge file absent: extension trials only", ""))]
  rt[n_trials == 0, `:=`(pass_pct = NA_real_, lo = NA_real_, hi = NA_real_, n_pass = NA_integer_)]
  rt[, `:=`(gt5_point = pass_pct > LIM, gt5_lower = lo > LIM)]
  rt <- merge(info, rt, by = "scenario"); rt[, model := mdl]
  rates[[mdl]] <- rt

  # (2) 분해: 같은 시험 쌍대 차이(a − b). 비교에 필요한 열이 모두 있는 시험만
  CON <- data.table(contrast = c("sel_ii", "sel_i", "extrap", "subst_ii", "subst_i", "total_ii", "total_i", "flag_A", "flag_C"),
                    a = c("G2_B", "G2_B", "G2_B", "G2_Cii", "G2_Ci", "G2_Aii", "G2_Ai", "G2_Ai", "G2_Ci"),
                    b = c("G2_Aii", "G2_Ai", "P2", "P2", "P2", "P2", "P2", "G2_Aii", "G2_Cii"),
                    flag_set = c("(ii)", "(i)", "-", "(ii)", "(i)", "(ii)", "(i)", "(i) vs (ii)", "(i) vs (ii)"),
                    meaning = c("dropout/selection effect (rule B minus rule A)", "dropout/selection effect (rule B minus rule A)", "extrapolation effect (rule B minus AUClast)",
                                "rule C minus P2", "rule C minus P2", "total G2-A minus P2 = extrapolation minus selection", "total G2-A minus P2 = extrapolation minus selection",
                                "flag set (i) minus (ii), rule A", "flag set (i) minus (ii), rule C"))
  dc <- rbindlist(lapply(seq_len(nrow(CON)), function(i) {
    ca <- paste0("cfg_", CON$a[i]); cb <- paste0("cfg_", CON$b[i])
    w[, { ok <- !is.na(get(ca)) & !is.na(get(cb))
          if (any(ok)) { d <- paired_prop_diff_ci(get(ca)[ok], get(cb)[ok]); .(n_trials = sum(ok), est = d$est, lo = d$lo, hi = d$hi, pct_a = 100 * mean(get(ca)[ok]), pct_b = 100 * mean(get(cb)[ok])) }
          else .(n_trials = 0L, est = NA_real_, lo = NA_real_, hi = NA_real_, pct_a = NA_real_, pct_b = NA_real_) }, by = scenario][, contrast := CON$contrast[i]]
  }))
  dc <- merge(CON, dc, by = "contrast")
  # 합의 항등식: (G2-A − P2) = (G2-B − P2) − (G2-B − G2-A), 같은 시험 집합에서 정확히
  for (fs in c("ii", "i")) {
    a_ <- paste0("cfg_G2_A", fs); ok <- !is.na(w[[a_]])
    if (!any(ok)) next
    chk <- w[ok, .(t = 100 * mean(get(a_) - cfg_P2), e = 100 * mean(cfg_G2_B - cfg_P2), s = 100 * mean(cfg_G2_B - get(a_))), by = scenario]
    stopifnot(all(abs(chk$t - (chk$e - chk$s)) < 1e-9))
    if (fs == "ii") stopifnot(all(abs(dc[contrast == "total_ii"][chk, on = "scenario", est] - chk$t) < 1e-9))
  }
  dc[, note := fifelse(n_trials == 0, "rejudge file absent: run scripts/40_oc_rejudge.R", "")]
  dc <- merge(info, dc, by = "scenario"); dc[, model := mdl]
  decomp[[mdl]] <- dc

  # (3) 치우침: mean(log GMR) − log(참값 비). AUC 계열은 참 AUC0-inf 비, Cmax는 참 Cmax 비. 경계 시나리오만
  eps_b <- c("Cmax", "AUClast", "AUCinf_A", "AUCinf_B", "AUCinf_C", "AUCinf_true", EP_I)
  bb <- be[endpoint %in% eps_b & scenario %in% info[region == "boundary", scenario] & is.finite(GMR) & GMR > 0,
           .(n_trials = .N, mean_log_gmr = mean(log(GMR)), sd_log_gmr = sd(log(GMR))), by = .(scenario, endpoint)]
  bb <- CJ(scenario = info[region == "boundary", scenario], endpoint = eps_b)[bb, on = c("scenario", "endpoint"), `:=`(n_trials = i.n_trials, mean_log_gmr = i.mean_log_gmr, sd_log_gmr = i.sd_log_gmr)]
  bb[is.na(n_trials), n_trials := 0L]
  bb <- merge(info, bb, by = "scenario")
  bb[, `:=`(truth_metric = fifelse(endpoint == "Cmax", "Cmax", "AUC0-inf"), true_ratio = fifelse(endpoint == "Cmax", cmax_ratio, auc_ratio),
            truth_se_log = fifelse(endpoint == "Cmax", cmax_se_log, auc_se_log))]
  bb[, `:=`(bias_log = mean_log_gmr - log(true_ratio), bias_mc_se = sd_log_gmr / sqrt(n_trials))]
  bb[, `:=`(bias_lo = bias_log - qnorm(0.975) * bias_mc_se, bias_hi = bias_log + qnorm(0.975) * bias_mc_se, bias_pct = 100 * (exp(bias_log) - 1),
            toward_one = sign(bias_log) == sign(-log(true_ratio)), bias_toward_one_log = bias_log * sign(-log(true_ratio)))]   # 양수 = 그 평가변수의 참값 비에서 1 쪽으로 치우침
  bb[, note := fifelse(n_trials == 0, "rejudge file absent: run scripts/40_oc_rejudge.R", fifelse(endpoint == "AUCinf_true", "reference: individual model AUC0-inf of trial subjects (no NCA)", ""))]
  bb[, model := mdl]
  bias[[mdl]] <- bb
}
RT <- rbindlist(rates); DC <- rbindlist(decomp); BI <- rbindlist(bias); ST <- rbindlist(status)
ord <- function(x) x[order(match(model, MODELS), region != "boundary", target, mechanism, direction)]
RT <- ord(RT); BI <- ord(BI)
DC <- DC[order(match(model, MODELS), region != "boundary", target, mechanism, direction, match(contrast, c("extrap", "sel_ii", "sel_i", "total_ii", "total_i", "subst_ii", "subst_i", "flag_A", "flag_C")))]
setcolorder(RT, c("model", "scenario", "region", "mechanism", "direction", "target", "multiplier", "auc_ratio", "cmax_ratio", "config", "family", "rule", "flag_set", "label",
                  "n_trials", "n_pass", "pass_pct", "lo", "hi", "gt5_point", "gt5_lower", "source", "note"))
RT[, c("auc_se_log", "cmax_se_log", "needs_i") := NULL]
RT <- RT[order(match(model, MODELS), region != "boundary", target, mechanism, direction, match(config, CFG$config))]
setcolorder(DC, c("model", "scenario", "region", "mechanism", "direction", "target", "multiplier", "auc_ratio", "contrast", "a", "b", "flag_set", "meaning", "n_trials", "pct_a", "pct_b", "est", "lo", "hi", "note"))
DC[, c("auc_se_log", "cmax_ratio", "cmax_se_log") := NULL]
BI <- BI[, .(model, scenario, mechanism, direction, target, multiplier, endpoint, truth_metric, true_ratio, truth_se_log, n_trials, mean_log_gmr, bias_log, bias_mc_se, bias_lo, bias_hi, bias_pct, toward_one, bias_toward_one_log, note)]
BI <- BI[order(match(model, MODELS), target, mechanism, direction, match(endpoint, c("Cmax", "AUClast", "AUCinf_A", "AUCinf_B", "AUCinf_C", "AUCinf_Ai", "AUCinf_Ci", "AUCinf_true")))]
fwrite(RT, file.path(out_dir, "g2_rules_flags.csv")); fwrite(DC, file.path(out_dir, "g2_decomposition.csv")); fwrite(BI, file.path(out_dir, "p2_bias_boundary.csv"))

# =============================================================================================
# 자동 문구. 모든 수치는 위 표에서, 고정 문구의 전제는 stopifnot으로 검사
B <- RT[region == "boundary"]
has_i <- function(mdl) ST[model == mdl, rejudge]
# (a) 구성별 5% 초과 요약
EX <- B[!is.na(pass_pct), .(n_scen = .N, n_gt5 = sum(gt5_point), n_lo_gt5 = sum(gt5_lower), max_pct = max(pass_pct), max_scn = scenario[which.max(pass_pct)],
                           max_lo = lo[which.max(pass_pct)], max_hi = hi[which.max(pass_pct)], min_pct = min(pass_pct), n_trials = max(n_trials)), by = .(model, config, label, family)]
EX <- EX[order(match(model, MODELS), match(config, CFG$config))]
stopifnot(all(EX$n_lo_gt5 <= EX$n_gt5))                                         # Wilson 하한 > 5%면 점추정 > 5%
# (b) 규칙 C가 P2에 가까워지는가: |G2-C − P2| < |G2-A − P2| (세트별, 경계 시나리오 전부)
CP <- rbindlist(lapply(c("ii", "i"), function(fs) {
  x <- DC[region == "boundary" & contrast %in% paste0(c("subst_", "total_"), fs) & n_trials > 0]
  if (!nrow(x)) return(NULL)
  y <- dcast(x, model + scenario ~ contrast, value.var = "est")
  setnames(y, paste0(c("subst_", "total_"), fs), c("C_P2", "A_P2"))
  y[, flag := fs][]
}), use.names = TRUE)
CP[, closer := abs(C_P2) < abs(A_P2)]
cmax_pass <- rbindlist(lapply(MODELS, function(mdl) W[[mdl]][scenario %in% B$scenario, .(model = mdl, cmax_pct = 100 * mean(Cmax_pass), n_trials = .N), by = scenario]))
CP <- merge(CP, cmax_pass, by = c("model", "scenario"))
# 동률(두 차이가 모두 0)이면서 Cmax가 한 번도 통과하지 않은 시나리오: G2·P2가 규칙과 무관하게 0%(구성 정의상 G2, P2 ≤ Cmax 통과)
CP[, tie_cmax0 := !closer & abs(C_P2) < 1e-12 & abs(A_P2) < 1e-12 & cmax_pct == 0]
# (c) 치우침 방향: G2-A(ii)가 5%를 넘는 시나리오에서 AUCinf_A GMR이 1 쪽으로 치우치는가, AUClast는?
BX <- merge(BI[endpoint %in% c("AUClast", "AUCinf_A", "AUCinf_B", "AUCinf_C")], B[config == "G2_Aii", .(model, scenario, g2a_gt5 = gt5_point)], by = c("model", "scenario"))

scn_lab <- function(s, lang, mdl) {
  r <- B[model == mdl & scenario == s][1]
  if (lang == "ko") sprintf("%s(참값 %.2f, ×%.3g)", s, r$target, r$multiplier) else sprintf("%s (true %.2f, x%.3g)", s, r$target, r$multiplier)
}
ci_k <- function(p, l, h) sprintf("%s%% [%s, %s]", f2(p), f2(l), f2(h)); ci_e <- function(p, l, h) sprintf("%s%% (95%% CI %s to %s)", f2(p), f2(l), f2(h))
cell_k <- function(p, l, h) ifelse(is.na(p), "NA", sprintf("%s [%s, %s]", f2(p), f2(l), f2(h))); cell_e <- function(p, l, h) ifelse(is.na(p), "NA", sprintf("%s [%s to %s]", f2(p), f2(l), f2(h)))
dcell_k <- function(e, l, h) ifelse(is.na(e), "NA", sprintf("%s [%s, %s]", s2(e), s2(l), s2(h))); dcell_e <- function(e, l, h) ifelse(is.na(e), "NA", sprintf("%s [%s to %s]", s2(e), s2(l), s2(h)))
trials_txt <- function(mdl, lang) {
  s <- ST[model == mdl]
  n <- if (s$n_trials_min == s$n_trials_max) fint(s$n_trials_min) else sprintf("%s-%s", fint(s$n_trials_min), fint(s$n_trials_max))   # 사전 고정 파일의 시험 수
  e <- if (s$extension) oc_n_trials_text(NT[[mdl]], lang, DEC, exceptions_only = TRUE) else ""                                  # 적응적 연장: 실제 수(R/oc.R)
  e <- if (nzchar(e)) paste0("; ", e) else ""
  if (lang == "ko") sprintf("%s: 시나리오당 %s회%s%s", MODEL_KO[[mdl]], n, if (s$complete) "(완료)" else "(scripts/31 진행 중인 부분 파일; 완료 후 같은 코드로 다시 생성)", e)
  else sprintf("%s: %s trials per scenario%s%s", MODEL_EN[[mdl]], n, if (s$complete) " (complete)" else " (partial file while scripts/31 is still running; regenerate with the same code when complete)", e)
}
ext_txt <- function(lang) {                                          # 연장 파일이 있을 때만 머리말에 한 줄
  x <- rbindlist(EXTN)
  if (!nrow(x)) return(character(0))
  if (lang == "ko") sprintf("- 적응적 연장: %s. `oc_trials_ext_be_<모델>.csv.gz`(scripts/42: 사전 고정 %s회에서 P2 경계 1종 오류의 Wilson 95%% 구간이 5%%를 포함한 경계 시나리오를 같은 시드 규칙으로 이어서 모의, 평가변수 8개)의 시험을 모든 표에 더했다.",
                            paste(sprintf("%s %s(연장 시험 %s회)", MODEL_KO[x$model], x$scenario, format(x$n_ext, big.mark = ",", trim = TRUE)), collapse = "; "), fint(N_BND))
  else sprintf("- Adaptive extension: %s. The trials of `oc_trials_ext_be_<model>.csv.gz` (scripts/42: boundary scenarios whose P2 boundary type I error Wilson 95%% CI at the pre-registered %s trials includes 5%%, continued with the same seed rule, all 8 endpoints) are added to every table.",
               paste(sprintf("%s %s (%s extension trials)", MODEL_EN[x$model], x$scenario, format(x$n_ext, big.mark = ",", trim = TRUE)), collapse = "; "), fint(N_BND))
}
rj_txt <- function(mdl, lang) {
  s <- ST[model == mdl]
  if (!s$rejudge) return(if (lang == "ko") sprintf("%s: 재판정 파일 없음 → 세트 (i) 구성(G2-A(i), G2-C(i), AUCinf-A(i), AUCinf-C(i))은 NA%s. `scripts/40_oc_rejudge.R %s` 실행 후 다시 생성.", MODEL_KO[[mdl]],
                                                        if (s$extension) "(적응적 연장 시나리오는 연장 시험만으로 계산)" else "", mdl)
                               else sprintf("%s: no rejudge file, so the flag set (i) configurations (G2-A(i), G2-C(i), AUCinf-A(i), AUCinf-C(i)) are NA%s. Run `scripts/40_oc_rejudge.R %s` and regenerate.", MODEL_EN[[mdl]],
                                            if (s$extension) " (the adaptively extended scenarios use the extension trials only)" else "", mdl))
  if (lang == "ko") sprintf("%s: 재판정 파일의 시험 %s회(원래 6개 평가변수가 저장본과 일치)%s.", MODEL_KO[[mdl]], fint(s$n_rejudge_trials), if (s$rejudge_complete) "" else ", 저장본보다 적어 세트 (i) 구성은 그 시험들로만 계산")
  else sprintf("%s: rejudge file with %s trials (the 6 original endpoints match the stored rows)%s.", MODEL_EN[[mdl]], fint(s$n_rejudge_trials), if (s$rejudge_complete) "" else "; fewer than stored, so the flag set (i) configurations use those trials only")
}
# 5% 초과 문장
exceed_txt <- function(mdl, fam, lang) {
  x <- EX[model == mdl & family %in% fam]
  yes <- x[n_gt5 > 0]; no <- x[n_gt5 == 0]
  na_ <- setdiff(CFG[family %in% fam, config], x$config)
  lab <- function(r) sprintf("%s %d/%d", r$label, r$n_gt5, r$n_scen)
  if (lang == "ko") {
    a <- if (nrow(yes)) paste(vapply(seq_len(nrow(yes)), function(i) { r <- yes[i]
      sprintf("%s개 시나리오(Wilson 하한 > 5%%: %d개; 최대 %s, %s)", lab(r), r$n_lo_gt5, ci_k(r$max_pct, r$max_lo, r$max_hi), r$max_scn) }, ""), collapse = "; ") else "없음"
    b <- if (nrow(no)) paste(vapply(seq_len(nrow(no)), function(i) { r <- no[i]; sprintf("%s(최대 %s, %s)", r$label, ci_k(r$max_pct, r$max_lo, r$max_hi), r$max_scn) }, ""), collapse = "; ") else "없음"
    c(sprintf("  - 점추정 > 5%%: %s.", a), sprintf("  - 모든 경계 시나리오에서 ≤ 5%%: %s.", b),
      if (length(na_)) sprintf("  - 미산출(재판정 파일 없음): %s.", paste(CFG[config %in% na_, label], collapse = ", ")))
  } else {
    a <- if (nrow(yes)) paste(vapply(seq_len(nrow(yes)), function(i) { r <- yes[i]
      sprintf("%s scenarios (Wilson lower bound above 5%%: %d; highest %s, %s)", lab(r), r$n_lo_gt5, ci_e(r$max_pct, r$max_lo, r$max_hi), r$max_scn) }, ""), collapse = "; ") else "none"
    b <- if (nrow(no)) paste(vapply(seq_len(nrow(no)), function(i) { r <- no[i]; sprintf("%s (highest %s, %s)", r$label, ci_e(r$max_pct, r$max_lo, r$max_hi), r$max_scn) }, ""), collapse = "; ") else "none"
    c(sprintf("  - Point estimate above 5%%: %s.", a), sprintf("  - At or below 5%% in every boundary scenario: %s.", b),
      if (length(na_)) sprintf("  - Not computed (no rejudge file): %s.", paste(CFG[config %in% na_, label], collapse = ", ")))
  }
}
# 시나리오별 표
rate_table <- function(mdl, fam, lang) {
  cf <- CFG[family %in% fam, config]; x <- B[model == mdl & config %in% cf]
  hdr <- c(sprintf("| %s | %s |", if (lang == "ko") "경계 시나리오" else "Boundary scenario", paste(CFG[config %in% cf, label], collapse = " | ")), paste0("|", strrep("---|", length(cf) + 1)))
  rows <- vapply(unique(x$scenario), function(s) {
    y <- x[scenario == s][match(cf, config)]
    paste0("| ", scn_lab(s, lang, mdl), " | ", paste(if (lang == "ko") cell_k(y$pass_pct, y$lo, y$hi) else cell_e(y$pass_pct, y$lo, y$hi), collapse = " | "), " |") }, "")
  c(hdr, rows)
}
DCOLS <- c("sel_ii", "sel_i", "extrap", "subst_ii", "subst_i", "total_ii", "total_i")
DLAB <- c(sel_ii = "G2-B − G2-A(ii)", sel_i = "G2-B − G2-A(i)", extrap = "G2-B − P2", subst_ii = "G2-C(ii) − P2", subst_i = "G2-C(i) − P2", total_ii = "G2-A(ii) − P2", total_i = "G2-A(i) − P2")
DLAB_EN <- setNames(gsub(" − ", " minus ", DLAB), names(DLAB))
dec_table <- function(mdl, lang) {
  x <- DC[model == mdl & region == "boundary" & contrast %in% DCOLS]
  lab <- if (lang == "ko") DLAB else DLAB_EN
  hdr <- c(sprintf("| %s | %s |", if (lang == "ko") "경계 시나리오" else "Boundary scenario", paste(lab[DCOLS], collapse = " | ")), paste0("|", strrep("---|", length(DCOLS) + 1)))
  rows <- vapply(unique(x$scenario), function(s) {
    y <- x[scenario == s][match(DCOLS, contrast)]
    paste0("| ", scn_lab(s, lang, mdl), " | ", paste(if (lang == "ko") dcell_k(y$est, y$lo, y$hi) else dcell_e(y$est, y$lo, y$hi), collapse = " | "), " |") }, "")
  c(hdr, rows)
}
# 분해 문장(모델별): 외삽 효과 부호, 세트별 선택 효과 부호, 합이 가장 큰 시나리오의 분해(합 = 외삽 − 선택이 같은 시험 집합에서 성립할 때만)
dec_txt <- function(mdl, lang) {
  x <- DC[model == mdl & region == "boundary" & n_trials > 0]
  e <- x[contrast == "extrap"]; e_pos <- sum(e$lo > 0); e_neg <- sum(e$hi < 0); e0 <- nrow(e) - e_pos - e_neg
  out <- if (lang == "ko") sprintf("- 외삽 효과(G2-B − P2): 경계 %d개 중 %d개에서 양(95%% 구간 하한 > 0), %d개에서 음(상한 < 0), %d개는 0을 포함(범위 %s ~ %s%%p).", nrow(e), e_pos, e_neg, e0, s2(min(e$est)), s2(max(e$est)))
         else sprintf("- Extrapolation effect (G2-B minus P2): positive in %d of %d boundary scenarios (95%% CI lower bound above 0), negative in %d (upper bound below 0), CI including 0 in %d (range %s to %s points).", e_pos, nrow(e), e_neg, e0, s2(min(e$est)), s2(max(e$est)))
  for (fs in c("ii", "i")) {
    s <- x[contrast == paste0("sel_", fs)]; if (!nrow(s)) next
    t_ <- x[contrast == paste0("total_", fs)]
    n_neg <- sum(s$hi < 0); n_pos <- sum(s$lo > 0); n0 <- nrow(s) - n_neg - n_pos
    k <- t_[which.max(abs(est)), scenario]; tk <- t_[scenario == k, est]; ek <- e[scenario == k, est]; sk <- s[scenario == k, est]
    same_set <- t_[scenario == k, n_trials] == e[scenario == k, n_trials]
    if (same_set) stopifnot(abs(tk - (ek - sk)) < 1e-9)                  # 같은 시험 집합이면 항등식이 정확히 성립
    fl <- sprintf("(%s)", fs)
    out <- c(out, if (lang == "ko")
      sprintf("- 세트 %s 선택 효과(G2-B − G2-A%s): 경계 %d개 중 %d개에서 음(95%% 구간 상한 < 0; 플래그 대상 제외가 통과율을 높임), %d개에서 양(하한 > 0), %d개는 0을 포함(범위 %s ~ %s%%p). 합(G2-A%s − P2)이 가장 큰 시나리오는 %s(%s%%p)%s.",
              fl, fl, nrow(s), n_neg, n_pos, n0, s2(min(s$est)), s2(max(s$est)), fl, k, s2(tk), if (same_set) sprintf(" = 외삽 %s%%p − 선택 %s%%p", s2(ek), s2(sk)) else "(세트 (i) 시험 수가 저장본보다 적어 분해는 생략)")
      else sprintf("- Flag set %s selection effect (G2-B minus G2-A%s): negative in %d of %d boundary scenarios (95%% CI upper bound below 0; excluding flagged subjects raises the pass rate), positive in %d (lower bound above 0), CI including 0 in %d (range %s to %s points). The largest total (G2-A%s minus P2) is in %s (%s points)%s.",
                   fl, fl, n_neg, nrow(s), n_pos, n0, s2(min(s$est)), s2(max(s$est)), fl, k, s2(tk), if (same_set) sprintf(" = extrapolation %s points minus selection %s points", s2(ek), s2(sk)) else " (decomposition omitted because flag set (i) covers fewer trials than stored)"))
  }
  out
}
# 규칙 C 해석 문장: 모든 경계 시나리오에서 |G2-C − P2| < |G2-A − P2|일 때만 일반 해석, 아니면 예외 명시
c_txt <- function(mdl, fs, lang) {
  x <- CP[model == mdl & flag == fs]
  if (!nrow(x)) return(if (lang == "ko") sprintf("- 세트 (%s): 재판정 파일 없음(미산출).", fs) else sprintf("- Flag set (%s): no rejudge file (not computed).", fs))
  fl <- sprintf("(%s)", fs); ex <- x[!(closer)]
  rng_c <- sprintf("%s ~ %s", s2(min(x$C_P2)), s2(max(x$C_P2))); rng_a <- sprintf("%s ~ %s", s2(min(x$A_P2)), s2(max(x$A_P2)))
  rng_ce <- sprintf("%s to %s", s2(min(x$C_P2)), s2(max(x$C_P2))); rng_ae <- sprintf("%s to %s", s2(min(x$A_P2)), s2(max(x$A_P2)))
  if (!nrow(ex)) {
    stopifnot(all(x$closer))
    return(if (lang == "ko") sprintf("- 세트 %s: G2-C%s는 플래그 대상자에 AUClast를 대입하므로 P2에 가까워진다. 경계 %d개 모두에서 |G2-C%s − P2| < |G2-A%s − P2| (G2-C − P2 %s%%p, G2-A − P2 %s%%p).", fl, fl, nrow(x), fl, fl, rng_c, rng_a)
           else sprintf("- Flag set %s: G2-C%s substitutes AUClast for flagged subjects and therefore approaches P2. In all %d boundary scenarios |G2-C%s minus P2| < |G2-A%s minus P2| (G2-C minus P2 %s points, G2-A minus P2 %s points).", fl, fl, nrow(x), fl, fl, rng_ce, rng_ae))
  }
  ok <- x[(closer)]
  exk <- vapply(seq_len(nrow(ex)), function(i) { r <- ex[i]
    if (r$tie_cmax0) sprintf("%s(두 차이 모두 0: Cmax가 %s회 중 한 번도 통과하지 않아 G2와 P2가 규칙과 무관하게 0%%)", r$scenario, fint(r$n_trials))
    else sprintf("%s(|G2-C − P2| %s%%p ≥ |G2-A − P2| %s%%p)", r$scenario, f2(abs(r$C_P2)), f2(abs(r$A_P2))) }, "")
  exe <- vapply(seq_len(nrow(ex)), function(i) { r <- ex[i]
    if (r$tie_cmax0) sprintf("%s (both differences are 0: Cmax passed in none of %s trials, so G2 and P2 are 0%% whatever the rule)", r$scenario, fint(r$n_trials))
    else sprintf("%s (|G2-C minus P2| %s points is not below |G2-A minus P2| %s points)", r$scenario, f2(abs(r$C_P2)), f2(abs(r$A_P2))) }, "")
  if (lang == "ko") sprintf("- 세트 %s: \"G2-C는 AUClast 대입으로 P2에 가까워진다\"는 경계 %d개 중 %d개에서 성립한다(|G2-C − P2| < |G2-A − P2|; 그 시나리오들의 G2-C − P2 %s%%p). 예외: %s.",
                            fl, nrow(x), nrow(ok), if (nrow(ok)) sprintf("%s ~ %s", s2(min(ok$C_P2)), s2(max(ok$C_P2))) else "NA", paste(exk, collapse = "; "))
  else sprintf("- Flag set %s: \"G2-C approaches P2 because it substitutes AUClast\" holds in %d of %d boundary scenarios (|G2-C minus P2| < |G2-A minus P2|; G2-C minus P2 %s points in those). Exceptions: %s.",
               fl, nrow(ok), nrow(x), if (nrow(ok)) sprintf("%s to %s", s2(min(ok$C_P2)), s2(max(ok$C_P2))) else "NA", paste(exe, collapse = "; "))
}
# 치우침 문장. 크기는 100 × log 차이, "1 쪽" 부호(bias_toward_one_log > 0 = 범위 안쪽)로 적는다
bias_txt <- function(mdl, lang) {
  x <- BX[model == mdl]; ga <- x[endpoint == "AUCinf_A" & g2a_gt5 == TRUE]; al <- x[endpoint == "AUClast"]
  tr <- BI[model == mdl & endpoint == "AUCinf_true"]
  h <- function(v) f2(100 * v); rg <- function(v, lang) if (lang == "ko") sprintf("%s ~ %s", h(min(v)), h(max(v))) else sprintf("%s to %s", h(min(v)), h(max(v)))
  out <- character(0)
  if (nrow(ga)) {
    gi <- ga[(toward_one)]; go <- ga[!(toward_one)]
    stopifnot(all(gi$bias_toward_one_log > 0), all(go$bias_toward_one_log <= 0))
    out <- c(out, if (lang == "ko")
      sprintf("- G2-A(ii)가 5%%를 넘는 경계 %d개 중 %d개에서 AUCinf_A(규칙 A, 세트 (ii))의 GMR이 참 AUC0-inf 비보다 1 쪽(범위 안쪽)으로 치우친다%s%s.", nrow(ga), nrow(gi),
              if (nrow(gi)) sprintf("(크기 %s, 100 × log 차이)", rg(gi$bias_toward_one_log, "ko")) else "", if (nrow(go)) sprintf("; 예외 %s", paste(go$scenario, collapse = ", ")) else "")
      else sprintf("- In %d of the %d boundary scenarios where G2-A(ii) exceeds 5%%, the AUCinf_A (rule A, flag set (ii)) GMR is biased toward 1, that is into the acceptance range, relative to the true AUC0-inf ratio%s%s.", nrow(gi), nrow(ga),
                   if (nrow(gi)) sprintf(" (by %s, 100 x log difference)", rg(gi$bias_toward_one_log, "en")) else "", if (nrow(go)) sprintf("; exceptions %s", paste(go$scenario, collapse = ", ")) else ""))
  }
  ao <- al[!(toward_one)]; ai <- al[(toward_one)]
  stopifnot(all(ao$bias_toward_one_log <= 0), all(ai$bias_toward_one_log > 0))
  out <- c(out, if (lang == "ko")
    sprintf("- AUClast의 GMR은 참 AUC0-inf 비 대비 경계 %d개 중 %d개에서 1 반대쪽(범위 바깥쪽)으로 치우친다%s%s. 참고: 시험 대상자의 개인 모델 AUC0-inf(AUCinf_true)의 치우침은 %s(MC 표준오차 최대 %s).",
            nrow(al), nrow(ao), if (nrow(ao)) sprintf("(크기 %s)", rg(-ao$bias_toward_one_log, "ko")) else "",
            if (nrow(ai)) sprintf("; 안쪽: %s", paste(sprintf("%s %s", ai$scenario, h(ai$bias_toward_one_log)), collapse = ", ")) else "", rg(tr$bias_log, "ko"), h(max(tr$bias_mc_se)))
    else sprintf("- The AUClast GMR is biased away from 1 (out of the acceptance range) relative to the true AUC0-inf ratio in %d of %d boundary scenarios%s%s. Reference: the bias of the individual model AUC0-inf of the trial subjects (AUCinf_true) is %s (largest MC SE %s).",
                 nrow(ao), nrow(al), if (nrow(ao)) sprintf(" (by %s)", rg(-ao$bias_toward_one_log, "en")) else "",
                 if (nrow(ai)) sprintf("; toward 1 in %s", paste(sprintf("%s (%s)", ai$scenario, h(ai$bias_toward_one_log)), collapse = ", ")) else "", rg(tr$bias_log, "en"), h(max(tr$bias_mc_se))))
  out
}

# Cmax가 한 번도 통과하지 않은 경계 시나리오: G2·P2는 AUC 규칙과 무관하게 0%. 규칙 효과는 AUCinf 단독에서만 보인다
cmax_txt <- function(mdl, lang) {
  x <- cmax_pass[model == mdl & scenario %in% B[region == "boundary", scenario] & cmax_pct == 0]
  if (!nrow(x)) return(character(0))
  vapply(seq_len(nrow(x)), function(i) { r <- x[i]
    g <- B[model == mdl & scenario == r$scenario & family %in% c("P2", "G2") & !is.na(pass_pct)]; stopifnot(all(g$pass_pct == 0))
    ao <- B[model == mdl & scenario == r$scenario & family == "AUCinf_only" & !is.na(pass_pct)]
    if (lang == "ko") sprintf("- %s %s: Cmax가 시험 %s회 중 한 번도 통과하지 않아 G2와 P2는 AUC 규칙과 무관하게 0%%다. 규칙 차이는 AUCinf 단독에서만 보인다(%s).", MODEL_KO[[mdl]], r$scenario, fint(r$n_trials),
                              paste(sprintf("%s %s%%", ao$label, f2(ao$pass_pct)), collapse = ", "))
    else sprintf("- %s, %s: Cmax passed in none of %s trials, so G2 and P2 are 0%% whatever the AUC rule. The rule differences show only for AUCinf alone (%s).", MODEL_EN[[mdl]], r$scenario, fint(r$n_trials),
                 paste(sprintf("%s %s%%", ao$label, f2(ao$pass_pct)), collapse = ", ")) }, "")
}

hdr_ko <- c("# G2 경계 1종 오류: AUCinf 처리 규칙 × 플래그 세트 (검토 의견 W2 §2)", "",
  "자료: 저장된 시험 수준 결과(`results/oc/oc_trials_be_<모델>.csv.gz`, scripts/31: arm당 117명, B0, pooled t, 90% CI 80.00–125.00%), 세트 (i) 평가변수는 `oc_rejudge_be_<모델>.csv.gz`(scripts/40: 같은 시드로 경계 시나리오 + S00 재생성). 새 모의 없음. 생성: `scripts/41_oc_rules_flags.R`. 수치: `g2_rules_flags.csv`, `g2_decomposition.csv`, `p2_bias_boundary.csv`.", "",
  paste0("- 시험 수: ", vapply(MODELS, trials_txt, "", lang = "ko"), "."), paste0("- ", vapply(MODELS, rj_txt, "", lang = "ko")), ext_txt("ko"), "",
  "## 정의", "",
  "- 규칙 A: 플래그 대상(또는 λz 산출 불가) 제외. 규칙 B: λz 산출 가능 전원 포함(플래그와 무관하므로 (i) = (ii)). 규칙 C: 플래그 대상(또는 λz 산출 불가)은 AUCinf 자리에 AUClast 대입.",
  "- 플래그 세트 (i): λz 산출 가능 & adj R² ≥ 0.80 & AUC_%Extrap_obs ≤ 20%. 세트 (ii): (i) & span ratio ≥ 2(사전 고정 구성 G2·F3A·F3C의 규칙 A·C, D-039).",
  "- 구성: P2 = AUClast + Cmax, G2-x = AUCinf(규칙 x) + Cmax, AUCinf-x = AUCinf(규칙 x) 단독. 경계 시나리오 = 참 AUC0-inf 비 0.80 또는 1.25(200,000명 공통 난수 역산, 기전·방향별). 5% 초과 판단은 점추정과 Wilson 95% 하한 두 가지로 적는다.", "")
hdr_en <- c("# G2 boundary type I error: AUCinf handling rule by reliability flag set (review W2 section 2)", "",
  "Data: saved trial-level results (`results/oc/oc_trials_be_<model>.csv.gz`, scripts/31: 117 per arm, B0, pooled t, 90% CI within 80.00 to 125.00%); the flag set (i) endpoints come from `oc_rejudge_be_<model>.csv.gz` (scripts/40: boundary scenarios and S00 regenerated with the same seeds). No new simulation. Generated by `scripts/41_oc_rules_flags.R`. Numbers: `g2_rules_flags.csv`, `g2_decomposition.csv`, `p2_bias_boundary.csv`.", "",
  paste0("- Trials: ", vapply(MODELS, trials_txt, "", lang = "en"), "."), paste0("- ", vapply(MODELS, rj_txt, "", lang = "en")), ext_txt("en"), "",
  "## Definitions", "",
  "- Rule A: exclude flagged subjects (or lambda-z not estimable). Rule B: include every subject with an estimable lambda-z (does not depend on the flags, so (i) = (ii)). Rule C: use AUClast in place of AUCinf for flagged subjects (or lambda-z not estimable).",
  "- Flag set (i): lambda-z estimable, adjusted R2 at least 0.80 and AUC_%Extrap_obs at most 20%. Flag set (ii): (i) plus span ratio at least 2 (rules A and C of the prespecified configurations G2, F3A and F3C, D-039).",
  "- Configurations: P2 = AUClast + Cmax, G2-x = AUCinf (rule x) + Cmax, AUCinf-x = AUCinf (rule x) alone. Boundary scenario = true AUC0-inf ratio 0.80 or 1.25 (inverted with 200,000 CRN subjects, per mechanism and direction). Exceeding 5% is reported both for the point estimate and for the Wilson 95% lower bound.", "")
ko <- hdr_ko; en <- hdr_en
ko <- c(ko, "## 1. 5% 초과 여부 (경계 시나리오)", ""); en <- c(en, "## 1. Which combinations exceed 5% (boundary scenarios)", "")
for (mdl in MODELS) {
  ko <- c(ko, sprintf("- %s", MODEL_KO[[mdl]]), "  - G2 계열:", sub("^  ", "    ", exceed_txt(mdl, c("P2", "G2"), "ko")), "  - AUCinf 단독:", sub("^  ", "    ", exceed_txt(mdl, "AUCinf_only", "ko")))
  en <- c(en, sprintf("- %s", MODEL_EN[[mdl]]), "  - G2 family:", sub("^  ", "    ", exceed_txt(mdl, c("P2", "G2"), "en")), "  - AUCinf alone:", sub("^  ", "    ", exceed_txt(mdl, "AUCinf_only", "en")))
}
ko <- c(ko, "", unlist(lapply(MODELS, cmax_txt, lang = "ko")), ""); en <- c(en, "", unlist(lapply(MODELS, cmax_txt, lang = "en")), "")
ko <- c(ko, "## 2. 시나리오별 통과율 (%, Wilson 95% 구간)", ""); en <- c(en, "## 2. Pass rate by scenario (%, Wilson 95% CI)", "")
for (mdl in MODELS) {
  ko <- c(ko, sprintf("### %s", MODEL_KO[[mdl]]), "", rate_table(mdl, c("P2", "G2"), "ko"), "", rate_table(mdl, "AUCinf_only", "ko"), "")
  en <- c(en, sprintf("### %s", MODEL_EN[[mdl]]), "", rate_table(mdl, c("P2", "G2"), "en"), "", rate_table(mdl, "AUCinf_only", "en"), "")
}
ko <- c(ko, "## 3. 분해 (같은 시험의 쌍대 차이, %p, 95% 구간)", "",
        "- 항등식(같은 시험 집합에서 정확히 성립, stopifnot): G2-A − P2 = (G2-B − P2) − (G2-B − G2-A) = 외삽 효과 − 선택 효과. 선택 효과(G2-B − G2-A)가 음이면 플래그 대상 제외(규칙 A)가 경계 통과율을 높인다는 뜻이다.", "")
en <- c(en, "## 3. Decomposition (paired differences within trials, percentage points, 95% CI)", "",
        "- Identity (exact on the same trials, checked with stopifnot): G2-A minus P2 = (G2-B minus P2) minus (G2-B minus G2-A) = extrapolation effect minus selection effect. A negative selection effect (G2-B minus G2-A) means that excluding flagged subjects (rule A) raises the boundary pass rate.", "")
for (mdl in MODELS) {
  dk <- dec_txt(mdl, "ko"); de <- dec_txt(mdl, "en")
  ko <- c(ko, sprintf("### %s", MODEL_KO[[mdl]]), "", dk, "", dec_table(mdl, "ko"), "")
  en <- c(en, sprintf("### %s", MODEL_EN[[mdl]]), "", de, "", dec_table(mdl, "en"), "")
}
ko <- c(ko, "## 4. 규칙 C와 P2", ""); en <- c(en, "## 4. Rule C and P2", "")
for (mdl in MODELS) {
  ko <- c(ko, sprintf("- %s", MODEL_KO[[mdl]]), sub("^- ", "  - ", c(c_txt(mdl, "ii", "ko"), c_txt(mdl, "i", "ko"))))
  en <- c(en, sprintf("- %s", MODEL_EN[[mdl]]), sub("^- ", "  - ", c(c_txt(mdl, "ii", "en"), c_txt(mdl, "i", "en"))))
}
ko <- c(ko, "", "## 5. 치우침 (시험 GMR 대 참값 비, `p2_bias_boundary.csv`)", ""); en <- c(en, "", "## 5. Bias (trial GMR versus true ratio, `p2_bias_boundary.csv`)", "")
for (mdl in MODELS) {
  ko <- c(ko, sprintf("- %s", MODEL_KO[[mdl]]), sub("^- ", "  - ", bias_txt(mdl, "ko")))
  en <- c(en, sprintf("- %s", MODEL_EN[[mdl]]), sub("^- ", "  - ", bias_txt(mdl, "en")))
}
ko <- c(ko, "", "## 전제 검사(stopifnot)", "",
  "- 재판정 파일의 원래 6개 평가변수 = 저장본(pass·n 동일, GMR 상대 차이 ≤ 1e-6), 세트 (i)의 규칙 A 포함 수 ≥ 세트 (ii).",
  "- P2·G2-A(ii)·AUCinf-A(ii) = 사전 고정 구성 P2·G2·AUCinf_only(config_pass)와 시험마다 같음. G2·P2 통과율 ≤ Cmax 통과율.",
  "- 분해 항등식, Wilson 하한 > 5%인 칸은 점추정도 > 5%. 규칙 C 해석 문장은 |G2-C − P2| < |G2-A − P2|가 모든 경계 시나리오에서 성립할 때만 일반형으로 쓰고, 아니면 예외를 적는다. 동률 예외는 Cmax 통과 0회일 때만 그 이유로 적는다.")
en <- c(en, "", "## Premise checks (stopifnot)", "",
  "- The 6 original endpoints in the rejudge file equal the stored rows (pass and n identical, GMR relative difference at most 1e-6), and rule A under flag set (i) includes at least as many subjects as under (ii).",
  "- P2, G2-A(ii) and AUCinf-A(ii) equal the prespecified configurations P2, G2 and AUCinf_only (config_pass) in every trial. G2 and P2 pass rates are at most the Cmax pass rate.",
  "- The decomposition identity holds, and every cell whose Wilson lower bound exceeds 5% also has a point estimate above 5%. The rule C interpretation is written in its general form only when |G2-C minus P2| < |G2-A minus P2| in every boundary scenario; otherwise the exceptions are listed. A tie is attributed to Cmax only when Cmax passed in no trial.")
# 영문: 한글·em dash·en dash 없음
stopifnot(!any(grepl("[ᄀ-ᇿ㄰-㆏가-힣]", en)), !any(grepl("—|–|−", en)))
writeLines(ko, file.path(out_dir, "g2_rules_conclusion_ko.md")); writeLines(en, file.path(out_dir, "g2_rules_conclusion_en.md"))
print(ST)
print(EX[, .(model, label, n_scen, n_gt5, n_lo_gt5, max_pct = round(max_pct, 2), max_scn, n_trials)])
cat("written:", file.path(out_dir, c("g2_rules_flags.csv", "g2_decomposition.csv", "p2_bias_boundary.csv", "g2_rules_conclusion_ko.md", "g2_rules_conclusion_en.md")), sep = "\n  ")
