#!/usr/bin/env Rscript
# 39_reliability_flags.R — 신뢰 플래그 두 세트와 "후기 채혈이 촘촘해지면 신뢰 충족률이 낮아지는" 기전(검토 의견 W2 §1, D-039, D-044).
# 저장된 개인 수준 NCA(results/individual/nca_<variant>_20000.rds, 20,000명 × 일정 B0, Bminus, D1–D4)만 쓴다. 모의 없음.
#  플래그 세트 (i):  λz 산출 가능 & Rsq_adjusted ≥ 0.80 & AUC_%Extrap_obs ≤ 20%   (D-039 이전 정의)
#  플래그 세트 (ii): (i) & Span_ratio ≥ 2                                          (현행 reliable 열, D-039; 같음을 검사)
# 산출(results/reliability/):
#  1. reliability_by_schedule.csv          변형 × 일정: λz 산출 %, (i)·(ii) 충족률(Wilson 95% 구간), span 단독 손실 = (i) − (ii)
#  2. reliability_paired_vs_B0.csv         같은 대상자 B0 대비 증감(%p, 쌍대 95% CI), 기준 (c)의 세트별 판정, 기준 (d) 값(플래그와 무관).
#                                           (ii) 증감은 results/trials/schedule_decision_<variant>.csv의 c_reliable_gain_pp와 0.01%p 이내여야 한다(아니면 중단).
#  3. dropout_reasons_by_schedule.csv       탈락 사유(중복 포함)와 배타적 조합(열), dropout_reason_combinations_by_schedule.csv(긴 형식, nca_engine B0 표와 같은 꼴)
#  4. reliability_lz_window_by_schedule.csv λz 창(하한·상한·폭), 점 수, 반감기, span ratio, adj R² 분포(λz 산출 가능 대상자)
#     reliability_span_transition.csv       B0 → 추가 일정 대상자 이동: B0 신뢰 충족 → span 플래그(손실), 계속 충족, B0 span 플래그 → 충족(회복)
#     reliability_span_transition_windows.csv 손실 대상자의 명목 창(연구일) 변화 상위 유형
#     reliability_short_window_conversion.csv B0 (ii) 충족자 중 추가 일정에서 최단 명목 3점 창(Day 22 이후)으로 바뀐 대상자와 그중 span 플래그 비율
#     mechanism_summary_ko.md / mechanism_summary_en.md 수치로 생성하고 주장 방향은 stopifnot으로 검사(전제가 어긋나면 생성 실패)
#  5. fig_reliability_mechanism.png(한국어) / fig_reliability_mechanism_en.png(영문)
#  6. reliability_two_flag_sets_summary.csv 변형별 B0 (i)·(ii)와 본문에 쓰는 두 모델 범위(results/rationale/pillar1_two_model_range.csv와 대조)
# 사용법: nice -n 15 Rscript scripts/39_reliability_flags.R
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios(); rel <- get_nca_rules()$standard$reliability; dr <- design$decision_rule
out_dir <- proj_path("results", "reliability"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("reliability_flags", master_seed = NA, run_mode = "final", extra = list(input = "results/individual/nca_<variant>_20000.rds"))

VARIANTS <- c("base", "struct2020", "vmax080_both", "vmax125_both", "iiv150", "resid12", "weight_alt", "ada10")
MECH_VARIANTS <- c("base", "struct2020", "vmax080_both", "vmax125_both")     # 기전 분석 대상(검토 의견 W2 §1-4)
TWO_MODELS <- c("base", "struct2020")                                        # 본문 두 모델 범위
SCHED <- design$schedule_analysis                                            # B0, Bminus, D1, D2, D3, D4
DENSE <- c("D1", "D2", "D3", "D4")                                           # 후기 채혈 추가 일정
EXTRAP <- "AUC_%Extrap_obs"
EPS <- 1e-9
LABEL_EN <- c(base = "Kovalenko 2016 (primary)", struct2020 = "Kovalenko 2020 Model 1", vmax080_both = "Vmax x0.8 (both arms)",
              vmax125_both = "Vmax x1.25 (both arms)", iiv150 = "IIV omega2 x1.5", resid12 = "Proportional residual 12%",
              weight_alt = "Weight mean 72, SD 10, 50 to 90 kg", ada10 = "ADA-like subgroup 10%, ke x2 after Day 14")
SHORT_KO <- c(base = "2016 모델", struct2020 = "Model 1", vmax080_both = "Vmax ×0.8", vmax125_both = "Vmax ×1.25")
SHORT_EN <- c(base = "2016 model", struct2020 = "Model 1", vmax080_both = "Vmax x0.8", vmax125_both = "Vmax x1.25")
label_ko <- function(v) sc$variants[[v]]$label
days_of <- function(sh) get_schedule(design, sh)
B0_DAYS <- days_of("B0")
added_days <- function(sh) setdiff(days_of(sh), B0_DAYS)
# 반올림은 R/summarize.R의 round_half_away(0에서 먼 쪽, 음의 0 제거): 20,000명 비율은 0.005의 배수라 round()로는 같은 값이
# summary_en·key_numbers(fmt_num)와 0.01 다르게 찍힌다(예: -1.935 → -1.93 대 -1.94). 같은 문서에 두 표기가 나오지 않게 맞춘다
f1 <- function(x) formatC(round_half_away(x, 1), format = "f", digits = 1); f2 <- function(x) formatC(round_half_away(x, 2), format = "f", digits = 2)
s2 <- function(x) sprintf("%+.2f", round_half_away(x, 2))
q05 <- function(x) quantile(x, 0.05, na.rm = TRUE, names = FALSE); q95 <- function(x) quantile(x, 0.95, na.rm = TRUE, names = FALSE)

# 명목 채혈일(투여 후 일)로 되돌리기: 허용창 ±1일, Day 22 이후 간격 ≥ 3일이라 가장 가까운 명목일이 유일하다(아래 일치율로 확인)
nearest_nominal <- function(t, days) { days <- sort(days); i <- findInterval(t, (days[-1] + days[-length(days)]) / 2) + 1L; days[i] }
# Day 22(투여 후 21일) 이후 명목일로 만들 수 있는 가장 짧은 3점 창(일): 3점 창이 span ≥ 2를 얻으려면 반감기 ≤ 이 값/2
min3_span <- function(days) { d <- sort(days[days >= 21]); min(d[-(1:2)] - d[seq_len(length(d) - 2)]) }

# ---------------------------------------------------------------------------------------------
# 입력: 변형별 NCA 저장본. 두 플래그 세트를 원 열에서 다시 계산하고 저장된 플래그와 같은지 검사
load_nca <- function(v) {
  f <- proj_path("results", "individual", sprintf("nca_%s_20000.rds", v))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  x <- x[schedule %in% SCHED]
  x[, rel_i := (lambda_ok & Rsq_adjusted >= rel$adj_r2_min & get(EXTRAP) <= rel$extrap_max_pct) %in% TRUE]
  x[, rel_ii := (rel_i & Span_ratio >= rel$span_ratio_min) %in% TRUE]
  x[, window := Lambda_z_upper - Lambda_z_lower]
  stopifnot(all(x$rel_ii == (x$reliable %in% TRUE)),                                             # (ii) = reliable 열
            all(x$rel_i == (x$lambda_ok & !(x$flag_rsq %in% TRUE) & !(x$flag_extrap %in% TRUE))),   # (i) = 플래그 두 개 기준
            all(x[lambda_ok == TRUE, abs(Span_ratio - window / HL_Lambda_z) < 1e-9]),
            x[, uniqueN(id), by = schedule][, all(V1 == nrow(x) / uniqueN(x$schedule))])
  x[, variant := v][]
}
NCA <- list()
for (v in VARIANTS) { x <- load_nca(v); if (!is.null(x)) NCA[[v]] <- x else message("없음: ", v) }
stopifnot(all(TWO_MODELS %in% names(NCA)), all(MECH_VARIANTS %in% names(NCA)))
append_run_log(logfile, sprintf("loaded variants: %s", paste(names(NCA), collapse = ",")))
sched_levels <- function(s) factor(s, levels = SCHED)

# ---------------------------------------------------------------------------------------------
# 1. 일정별 신뢰 충족률(두 세트)
by_sched <- rbindlist(lapply(names(NCA), function(v) {
  NCA[[v]][, {
    wi <- wilson_ci(sum(rel_i), .N); wii <- wilson_ci(sum(rel_ii), .N)
    .(n = .N, n_points = length(days_of(.BY$schedule)), lambda_ok_pct = 100 * mean(lambda_ok),
      reliable_i_pct = wi$est, reliable_i_lo = wi$lo, reliable_i_hi = wi$hi,
      reliable_ii_pct = wii$est, reliable_ii_lo = wii$lo, reliable_ii_hi = wii$hi,
      span_only_loss_pct = wi$est - wii$est, flag_span_pct = 100 * mean(flag_span %in% TRUE))
  }, by = schedule][, `:=`(variant = v, model_label = label_ko(v), model_label_en = LABEL_EN[[v]])]
}))
setcolorder(by_sched, c("variant", "model_label", "model_label_en", "schedule"))
by_sched <- by_sched[order(match(variant, VARIANTS), sched_levels(schedule))]
stopifnot(all(by_sched$span_only_loss_pct >= -1e-9))                    # (ii) ⊂ (i)
# 이미 보고된 개인 수준 표와 대조: (ii) = reliable_pct, (i) = reliable_no_span_pct, λz 산출 %
for (v in names(NCA)) {
  f <- proj_path("results", "individual", sprintf("individual_%s.csv", v)); if (!file.exists(f)) next
  ind <- fread(f); m <- merge(by_sched[variant == v], ind[, .(schedule, reliable_pct, reliable_no_span_pct, lambda_ok_pct_ind = lambda_ok_pct)], by = "schedule")
  stopifnot(nrow(m) == nrow(by_sched[variant == v]), all(abs(m$reliable_ii_pct - m$reliable_pct) < 1e-9),
            all(abs(m$reliable_i_pct - m$reliable_no_span_pct) < 1e-9), all(abs(m$lambda_ok_pct - m$lambda_ok_pct_ind) < 1e-9))
}
fwrite(by_sched, file.path(out_dir, "reliability_by_schedule.csv"))

# ---------------------------------------------------------------------------------------------
# 2. B0 대비 쌍대 증감(같은 대상자). 기준 (c)는 두 세트로, 기준 (d)는 플래그와 무관(비구획 외삽 > 20%)이라 그대로
paired <- rbindlist(lapply(names(NCA), function(v) {
  x <- NCA[[v]]
  b0 <- x[schedule == "B0", .(id, i0 = rel_i, ii0 = rel_ii, fs0 = flag_span %in% TRUE, x0 = !is.na(pct_extrap) & pct_extrap > 20)]
  sd_f <- proj_path("results", "trials", sprintf("schedule_decision_%s.csv", v))
  sdc <- if (file.exists(sd_f)) fread(sd_f) else NULL
  rbindlist(lapply(setdiff(SCHED, "B0"), function(sh) {
    m <- merge(b0, x[schedule == sh, .(id, i1 = rel_i, ii1 = rel_ii, fs1 = flag_span %in% TRUE, x1 = !is.na(pct_extrap) & pct_extrap > 20)], by = "id")
    gi <- paired_prop_diff_ci(m$i1, m$i0); gii <- paired_prop_diff_ci(m$ii1, m$ii0)
    d_b0 <- x[schedule == "B0", 100 * mean(pct_extrap > 20, na.rm = TRUE)]; d_s <- x[schedule == sh, 100 * mean(pct_extrap > 20, na.rm = TRUE)]
    out <- data.table(variant = v, model_label = label_ko(v), model_label_en = LABEL_EN[[v]], schedule = sh, n = nrow(m),
      added_points = length(days_of(sh)) - length(B0_DAYS),
      reliable_i_B0_pct = 100 * mean(m$i0), reliable_i_pct = 100 * mean(m$i1), gain_i_pp = gi$est, gain_i_lo = gi$lo, gain_i_hi = gi$hi,
      n_i_gained = sum(!m$i0 & m$i1), n_i_lost = sum(m$i0 & !m$i1),
      reliable_ii_B0_pct = 100 * mean(m$ii0), reliable_ii_pct = 100 * mean(m$ii1), gain_ii_pp = gii$est, gain_ii_lo = gii$lo, gain_ii_hi = gii$hi,
      n_ii_gained = sum(!m$ii0 & m$ii1), n_ii_lost = sum(m$ii0 & !m$ii1),
      span_only_loss_B0_pct = 100 * (mean(m$i0) - mean(m$ii0)), span_only_loss_pct = 100 * (mean(m$i1) - mean(m$ii1)),
      flag_span_B0_pct = 100 * mean(m$fs0), flag_span_pct = 100 * mean(m$fs1), flag_span_change_pp = 100 * (mean(m$fs1) - mean(m$fs0)),
      n_reliable_to_span = sum(m$ii0 & m$fs1), n_span_to_reliable = sum(m$fs0 & m$ii1),
      crit_c_threshold_pp = dr$reliability_gain_pp_min, crit_c_i = gi$est >= dr$reliability_gain_pp_min, crit_c_ii = gii$est >= dr$reliability_gain_pp_min,
      d_extrap20_B0 = d_b0, d_extrap20 = d_s, d_extrap20_ratio = d_s / d_b0, crit_d = d_s / d_b0 <= dr$extrap_gt20_ratio_max)
    if (!is.null(sdc) && sh %in% sdc$schedule) {
      s <- sdc[schedule == sh]
      # (ii) 증감이 판정표와 0.01%p 이내가 아니면 중단(같은 저장본에서 나온 값이어야 한다)
      if (abs(s$c_reliable_gain_pp - gii$est) > 0.01) stop(sprintf("%s %s: (ii) 증감 %.4f ≠ schedule_decision c_reliable_gain_pp %.4f", v, sh, gii$est, s$c_reliable_gain_pp))
      stopifnot(abs(s$d_extrap20_ratio - out$d_extrap20_ratio) < 1e-9, abs(s$c_reliable_B0 - out$reliable_ii_B0_pct) < 1e-9)
      out[, `:=`(decision_c_reliable_gain_pp = s$c_reliable_gain_pp, decision_match_0.01pp = TRUE,
                 crit_a = s$crit_a, crit_b = s$crit_b, recommend_decision = s$recommend,
                 recommend_ii = (s$crit_a %in% TRUE) | (s$crit_b %in% TRUE) | crit_c_ii | (crit_d %in% TRUE),
                 recommend_i = (s$crit_a %in% TRUE) | (s$crit_b %in% TRUE) | crit_c_i | (crit_d %in% TRUE))]
      stopifnot(out$recommend_ii == out$recommend_decision)
    }
    out
  }), fill = TRUE)                                   # 판정표가 D2만 있는 변형(iiv150, resid12, weight_alt, ada10)은 나머지 일정의 판정 열이 NA
}), fill = TRUE)
paired <- paired[order(match(variant, VARIANTS), sched_levels(schedule))]
stopifnot(all(abs(paired$gain_ii_pp - paired$gain_i_pp + paired$span_only_loss_pct - paired$span_only_loss_B0_pct) < 1e-9))
fwrite(paired, file.path(out_dir, "reliability_paired_vs_B0.csv"))

# ---------------------------------------------------------------------------------------------
# 3. 탈락 사유: 중복 포함(λz 산출 불가, adj R² < 0.80, 외삽 > 20%, span < 2)과 배타적 조합
EXCL <- c("reliable", "lambda_fail", "rsq", "extrap", "span", "rsq+extrap", "rsq+span", "extrap+span", "rsq+extrap+span")
all_nca <- rbindlist(NCA, use.names = TRUE, fill = TRUE)
all_nca[, `:=`(fr = flag_rsq %in% TRUE, fe = flag_extrap %in% TRUE, fsp = flag_span %in% TRUE)]
all_nca[, combo := fifelse(!lambda_ok, "lambda_fail", fifelse(rel_ii, "reliable",
                    gsub("^\\+|\\+$", "", gsub("\\+\\++", "+", paste(fifelse(fr, "rsq", ""), fifelse(fe, "extrap", ""), fifelse(fsp, "span", ""), sep = "+")))))]
stopifnot(all(all_nca$combo %in% EXCL))
drop_w <- all_nca[, {
  r <- .(n = .N, n_quant0 = sum(n_quant == 0), lambda_fail_pct = 100 * mean(!lambda_ok), flag_rsq_pct = 100 * mean(fr), flag_extrap_pct = 100 * mean(fe),
         flag_span_pct = 100 * mean(fsp), any_flag_or_fail_pct = 100 * mean(!rel_ii))
  ex <- setNames(lapply(EXCL, function(k) 100 * mean(combo == k)), paste0("excl_", gsub("\\+", "_", EXCL), "_pct"))
  c(r, ex)
}, by = .(variant, schedule)]
drop_w[, `:=`(model_label = vapply(variant, label_ko, ""), model_label_en = unname(LABEL_EN[variant]))]
setcolorder(drop_w, c("variant", "model_label", "model_label_en", "schedule"))
excl_cols <- grep("^excl_", names(drop_w), value = TRUE)
stopifnot(all(abs(rowSums(drop_w[, ..excl_cols]) - 100) < 1e-9), all(abs(drop_w$excl_reliable_pct + drop_w$any_flag_or_fail_pct - 100) < 1e-9))
drop_w <- drop_w[order(match(variant, VARIANTS), sched_levels(schedule))]
fwrite(drop_w, file.path(out_dir, "dropout_reasons_by_schedule.csv"))
drop_l <- all_nca[, .N, by = .(variant, schedule, lambda_fail = !lambda_ok, rsq = fr, extrap = fe, span = fsp)][, pct := 100 * N / sum(N), by = .(variant, schedule)]
drop_l <- drop_l[order(match(variant, VARIANTS), sched_levels(schedule), -pct)]
setnames(drop_l, "N", "n")
fwrite(drop_l, file.path(out_dir, "dropout_reason_combinations_by_schedule.csv"))

# ---------------------------------------------------------------------------------------------
# 4-1. λz 창 분포(λz 산출 가능 대상자): 하한·상한·폭(상한 − 하한), 점 수(3, 4, 5, 6+), 반감기, span ratio, adj R²
dist4 <- function(x, nm) { x <- as.numeric(x); setNames(list(median(x), q05(x), q95(x), mean(x)), paste0(nm, c("_median", "_p05", "_p95", "_mean"))) }
lz <- rbindlist(lapply(MECH_VARIANTS, function(v) {
  NCA[[v]][lambda_ok == TRUE, {
    c(list(n = uniqueN(NCA[[v]][schedule == .BY$schedule, id]), n_lambda_ok = .N, span_flag_pct_of_lambda_ok = 100 * mean(Span_ratio < rel$span_ratio_min),
           min_3pt_nominal_window_days = min3_span(days_of(.BY$schedule))),
      dist4(Lambda_z_lower, "lower"), dist4(Lambda_z_upper, "upper"), dist4(window, "window"),
      dist4(No_points_lambda_z, "npts"),
      list(npts3_pct = 100 * mean(No_points_lambda_z == 3), npts4_pct = 100 * mean(No_points_lambda_z == 4),
           npts5_pct = 100 * mean(No_points_lambda_z == 5), npts6plus_pct = 100 * mean(No_points_lambda_z >= 6)),
      dist4(HL_Lambda_z, "HL"), dist4(Span_ratio, "span"), dist4(Rsq_adjusted, "rsq"))
  }, by = schedule][, `:=`(variant = v, model_label = label_ko(v), model_label_en = LABEL_EN[[v]])]
}))
setcolorder(lz, c("variant", "model_label", "model_label_en", "schedule"))
lz <- lz[order(match(variant, VARIANTS), sched_levels(schedule))]
stopifnot(all(abs(lz$npts3_pct + lz$npts4_pct + lz$npts5_pct + lz$npts6plus_pct - 100) < 1e-9))
fwrite(lz, file.path(out_dir, "reliability_lz_window_by_schedule.csv"))

# 4-2. 대상자 이동 B0 → 일정 s(같은 대상자, 공통 관측치: 추가 일정은 B0 관측치에 점만 더한다)
#  lost_span:        B0 (ii) 충족 → s에서 span 플래그(다른 플래그 동반 가능)
#  stayed_reliable:  B0·s 모두 (ii) 충족
#  gained_from_span: B0 span 플래그 → s에서 (ii) 충족
grp_stats <- function(d) {
  d[, .(n = .N,
    window_B0_median = median(W0), window_median = median(W1), HL_B0_median = median(HL0), HL_median = median(HL1),
    lower_B0_median = median(lo0), lower_median = median(lo1), upper_B0_median = median(up0), upper_median = median(up1),
    npts_B0_median = as.numeric(median(n0)), npts_median = as.numeric(median(n1)), span_B0_median = median(SR0), span_median = median(SR1),
    rsq_B0_median = median(R0), rsq_median = median(R1),
    d_window_median = median(W1 - W0), d_HL_median = median(HL1 - HL0), d_lower_median = median(lo1 - lo0),
    d_upper_median = median(up1 - up0), d_npts_median = as.numeric(median(n1 - n0)), d_rsq_median = median(R1 - R0),
    window_unchanged_pct = 100 * mean(abs(lo1 - lo0) < EPS & abs(up1 - up0) < EPS & n1 == n0),
    window_shorter_pct = 100 * mean(W1 < W0 - EPS), window_longer_pct = 100 * mean(W1 > W0 + EPS),
    HL_longer_pct = 100 * mean(HL1 > HL0 * (1 + EPS)), HL_shorter_pct = 100 * mean(HL1 < HL0 * (1 - EPS)),
    lower_later_pct = 100 * mean(lo1 > lo0 + EPS), upper_later_pct = 100 * mean(up1 > up0 + EPS),
    npts_fewer_pct = 100 * mean(n1 < n0), npts3_B0_pct = 100 * mean(n0 == 3), npts3_pct = 100 * mean(n1 == 3),
    rsq_higher_pct = 100 * mean(R1 > R0 + EPS),
    mean_dlog_window = mean(log(W1 / W0)), mean_dlog_HL = mean(log(HL1 / HL0)), mean_dlog_span = mean(log(SR1 / SR0)),
    window_with_added_day_pct = 100 * mean(k_add1 >= 1), HL_lt1_pct = 100 * mean(HL1 < 1),
    other_flag_pct = 100 * mean(fr1 | fe1),
    nominal_map_consistent_pct = 100 * mean(k_nom1 == n1 & k_nom0 == n0))]
}
count_in <- function(lo, up, days) if (!length(days)) integer(length(lo)) else as.integer(rowSums(outer(lo, days, "<=") & outer(up, days, ">=")))
trans_list <- list(); win_list <- list(); fig_pairs <- list(); conv_list <- list()
for (v in MECH_VARIANTS) {
  x <- NCA[[v]]
  b0 <- x[schedule == "B0", .(id, ii0 = rel_ii, fs0 = flag_span %in% TRUE, l0 = lambda_ok, W0 = window, HL0 = HL_Lambda_z, lo0 = Lambda_z_lower, up0 = Lambda_z_upper,
                             n0 = No_points_lambda_z, SR0 = Span_ratio, R0 = Rsq_adjusted)]
  b0[, `:=`(nlo0 = nearest_nominal(lo0, B0_DAYS), nup0 = nearest_nominal(up0, B0_DAYS))]
  b0[, k_nom0 := count_in(nlo0, nup0, B0_DAYS)]
  for (sh in setdiff(SCHED, "B0")) {
    dsh <- days_of(sh); add <- added_days(sh)
    m <- merge(b0, x[schedule == sh, .(id, ii1 = rel_ii, fs1 = flag_span %in% TRUE, fr1 = flag_rsq %in% TRUE, fe1 = flag_extrap %in% TRUE, l1 = lambda_ok,
                                      W1 = window, HL1 = HL_Lambda_z, lo1 = Lambda_z_lower, up1 = Lambda_z_upper, n1 = No_points_lambda_z, SR1 = Span_ratio, R1 = Rsq_adjusted)], by = "id")
    m[, `:=`(nlo1 = nearest_nominal(lo1, dsh), nup1 = nearest_nominal(up1, dsh))]
    m[, `:=`(k_nom1 = count_in(nlo1, nup1, dsh), k_add1 = count_in(nlo1, nup1, add))]
    m[, grp := fifelse(ii0 & fs1, "lost_span", fifelse(ii0 & ii1, "stayed_reliable", fifelse(fs0 & ii1, "gained_from_span", NA_character_)))]
    g <- m[!is.na(grp), grp_stats(.SD), by = grp]
    g[, `:=`(variant = v, model_label = label_ko(v), model_label_en = LABEL_EN[[v]], schedule = sh, n_subjects = nrow(m),
             pct_of_subjects = 100 * n / nrow(m), pct_of_B0_reliable = fifelse(grp == "gained_from_span", NA_real_, 100 * n / sum(m$ii0)),
             window_share_of_dlog_span = fifelse(abs(mean_dlog_span) > 0.05, mean_dlog_window / mean_dlog_span, NA_real_))]
    trans_list[[paste(v, sh)]] <- g
    # 반증 검사: 짧은 창이 곧 span 플래그인가? B0 (ii) 충족자 중 추가 일정에서 최단 명목 3점 창(Day 22 이후)으로 바뀐 대상자의 span 플래그 비율과
    # 그 창의 적합 반감기. 창 폭만으로 플래그가 정해진다면 이 비율이 높아야 한다(반감기가 짧을 때 주로 그 창이 선택되면 낮다)
    if (sh %in% DENSE) {
      m3s <- min3_span(dsh)
      sw <- m[ii0 & (l1 %in% TRUE) & nlo1 >= 21 & (nup1 - nlo1) <= m3s]
      conv_list[[paste(v, sh)]] <- data.table(variant = v, model_label = label_ko(v), model_label_en = LABEL_EN[[v]], schedule = sh, min3_window_days = m3s, n_B0_reliable = sum(m$ii0),
        n_short = nrow(sw), short_pct_of_B0_reliable = 100 * nrow(sw) / sum(m$ii0), span_flag_pct = 100 * mean(sw$fs1),
        HL_B0_median = median(sw$HL0), HL_median = median(sw$HL1), npts3_pct = 100 * mean(sw$n1 == 3),
        n_lost = sum(m$ii0 & m$fs1), lost_in_short_pct = 100 * sum(sw$fs1) / sum(m$ii0 & m$fs1))
    }
    # 손실 대상자의 명목 창 변화(연구일 = 투여 후 일 + 1)
    lz_ <- m[grp == "lost_span"]
    if (nrow(lz_)) {
      w <- lz_[, .N, by = .(window_B0 = sprintf("Day %g-%g", nlo0 + 1, nup0 + 1), npts_B0 = n0, window_sched = sprintf("Day %g-%g", nlo1 + 1, nup1 + 1), npts_sched = n1)]
      w[, `:=`(pct_of_lost = 100 * N / sum(N), variant = v, schedule = sh)]
      win_list[[paste(v, sh)]] <- head(w[order(-N)], 5)
    }
    if (v %in% TWO_MODELS && sh %in% DENSE) fig_pairs[[paste(v, sh)]] <- g[grp == "lost_span", .(variant, schedule, n, window_B0_median, window_median, HL_B0_median, HL_median)]
  }
}
trans <- rbindlist(trans_list)
setcolorder(trans, c("variant", "model_label", "model_label_en", "schedule", "grp", "n", "n_subjects", "pct_of_subjects", "pct_of_B0_reliable"))
trans <- trans[order(match(variant, VARIANTS), sched_levels(schedule), match(grp, c("lost_span", "stayed_reliable", "gained_from_span")))]
fwrite(trans, file.path(out_dir, "reliability_span_transition.csv"))
wins <- rbindlist(win_list); setnames(wins, "N", "n"); setcolorder(wins, c("variant", "schedule"))
fwrite(wins, file.path(out_dir, "reliability_span_transition_windows.csv"))
CV <- rbindlist(conv_list)
fwrite(CV, file.path(out_dir, "reliability_short_window_conversion.csv"))

# ---------------------------------------------------------------------------------------------
# 6. 변형별 B0 두 세트 요약과 본문 두 모델 범위
b0s <- by_sched[schedule == "B0"]
rng <- b0s[variant %in% TWO_MODELS, .(i_min = min(reliable_i_pct), i_max = max(reliable_i_pct), ii_min = min(reliable_ii_pct), ii_max = max(reliable_ii_pct),
                                       loss_min = min(span_only_loss_pct), loss_max = max(span_only_loss_pct))]
p1 <- proj_path("results", "rationale", "pillar1_two_model_range.csv")
if (file.exists(p1)) { pr <- fread(p1); stopifnot(abs(pr$reliable_min - rng$ii_min) < 1e-9, abs(pr$reliable_max - rng$ii_max) < 1e-9) }
eng <- fread(proj_path("results", "nca_engine", "engine_difference_individual_B0.csv"))[grepl("D-039", engine)]
two <- b0s[, .(variant, model_label, model_label_en, n, lambda_ok_pct, reliable_i_pct, reliable_i_lo, reliable_i_hi,
               reliable_ii_pct, reliable_ii_lo, reliable_ii_hi, span_only_loss_pct,
               two_model_i_min = rng$i_min, two_model_i_max = rng$i_max, two_model_ii_min = rng$ii_min, two_model_ii_max = rng$ii_max,
               two_model_range_i_text = sprintf("%s to %s", f1(rng$i_min), f1(rng$i_max)), two_model_range_ii_text = sprintf("%s to %s", f1(rng$ii_min), f1(rng$ii_max)))]
# 참고: scripts/29(엔진 차이)는 같은 대상자에 B0 격자만으로 관측 시각·잔차를 따로 뽑았다 → 별도 표본. 비교용으로만 옆에 둔다.
two <- merge(two, eng[, .(variant = model, nca_engine_B0only_i_pct = reliable_rsq_extrap_only_pct, nca_engine_B0only_ii_pct = reliable_pct)], by = "variant", all.x = TRUE)
two <- two[order(match(variant, VARIANTS))]
fwrite(two, file.path(out_dir, "reliability_two_flag_sets_summary.csv"))

# ---------------------------------------------------------------------------------------------
# 4-3. 기전 판정: 손실(lost_span) 대상자, 기전 변형 4개 × D1–D4, 대상자 50명 이상인 칸
L <- trans[grp == "lost_span" & schedule %in% DENSE & n >= 50]
S <- trans[grp == "stayed_reliable" & schedule %in% DENSE]
G <- trans[grp == "gained_from_span" & schedule %in% DENSE & n >= 50]
m3_b0 <- min3_span(B0_DAYS); m3_dense <- vapply(DENSE, function(s) min3_span(days_of(s)), numeric(1))
# 전제 검사(주장의 방향). 하나라도 어긋나면 문구를 쓰지 않고 중단한다.
stopifnot(
  nrow(L) == length(MECH_VARIANTS) * length(DENSE),     # 모든 칸에 판정할 만큼의 손실 대상자
  m3_b0 == 14, all(m3_dense == 7),                      # 명목 3점 창 최소 폭: B0 14일, 추가 일정 7일
  all(L$window_shorter_pct >= 90), all(L$d_window_median < 0),          # 창이 짧아진다(거의 전원)
  all(L$d_lower_median > 0), all(L$lower_later_pct >= 90),              # 창 시작이 늦어진다
  all(L$npts3_pct > L$npts3_B0_pct), all(L$npts3_pct >= 75),            # 3점 창으로 바뀐다
  all(L$d_HL_median < 0), all(L$HL_longer_pct < 25),                    # 반감기는 길어지지 않는다(대부분 짧아진다)
  all(L$mean_dlog_window < L$mean_dlog_span), all(L$mean_dlog_span < 0), # 창 폭 감소가 span 감소 전부보다 크고 반감기 변화가 일부 상쇄
  all(L$d_rsq_median > 0), all(L$rsq_higher_pct > 50),                  # 새 창의 adj R²가 더 높다(Best Fit가 그 창을 고르는 이유)
  all(L$window_with_added_day_pct >= 90),                               # 새 창에 추가 채혈일이 들어 있다
  all(L$HL_median > m3_dense[1] / 2), all(L$HL_median < m3_b0 / 2))     # 새 창 반감기 중앙값: 7일 창 기준(3.5일)보다 길고 14일 창 기준(7일)보다 짧다
stopifnot(all(L$nominal_map_consistent_pct >= 95))                     # 명목일 대응의 타당성
stopifnot(all(G$upper_later_pct >= 50), all(G$window_median > G$window_B0_median))   # 회복군: 창 상한(Tlast)이 늦어지고 창이 길어진다
# 회복군은 창이 길어지는 것만이 아니다: 반감기도 짧아진다(span 증가의 일부 또는 대부분). 기여 비교는 칸 수로 보고
stopifnot(all(G$HL_shorter_pct > 50), all(G$mean_dlog_HL < 0), all(G$mean_dlog_span > 0), !anyNA(G$window_share_of_dlog_span))
G_hl_major <- sum(G$window_share_of_dlog_span < 0.5)     # 반감기 단축의 기여(Δlog 반감기의 절댓값)가 창 폭 증가보다 큰 칸 수
# 반증 검사(짧은 창 ≠ 플래그): 최단 명목 3점 창으로 바뀐 B0 충족자 중 span 플래그는 소수이고, 그 창의 반감기 중앙값은 한계(창 폭/2)보다 짧고 B0 창보다 짧다
stopifnot(nrow(CV) == length(MECH_VARIANTS) * length(DENSE), all(CV$min3_window_days == m3_dense[1]), all(CV$n_short >= 50),
          all(CV$span_flag_pct < 50), all(CV$HL_median < CV$min3_window_days / 2), all(CV$HL_median < CV$HL_B0_median), all(CV$npts3_pct >= 99))
# 절벽(D-041, 1일 정의) 길이와 추가 일정의 Day 22 이후 최소 명목 간격: 절벽이 간격보다 짧으면 3점 창이 절벽 안에 모두 들어갈 수 없다(명목일 기준)
min_gap_dense <- min(vapply(DENSE, function(s) { d <- sort(days_of(s)); d <- d[d >= 21]; min(diff(d)) }, numeric(1)))
cliff_f <- proj_path("results", "cliff", "cliff_summary.csv")
cliff_p95 <- if (file.exists(cliff_f)) fread(cliff_f)[weight == "base", max(len1_p95)] else NA_real_
if (!is.na(cliff_p95)) stopifnot(cliff_p95 < min_gap_dense)
# 범위 표기: 한국어는 "a–b"(음수가 있으면 "a ~ b"), 영문은 "a to b". 양끝이 같으면 한 값
rtx <- function(x, fmt = f1, lang = "ko") {
  r <- range(x, na.rm = TRUE); a <- fmt(r[1]); b <- fmt(r[2])
  if (a == b) return(a)
  if (lang == "en") return(paste(a, "to", b))
  if (r[1] < 0) paste(a, "~", b) else paste0(a, "–", b)
}
rk <- function(col, fmt = f1, d = L) rtx(d[[col]], fmt, "ko"); re <- function(col, fmt = f1, d = L) rtx(d[[col]], fmt, "en")
Lb3 <- L[variant == "base" & schedule == "D3"]; Ls3 <- L[variant == "struct2020" & schedule == "D3"]
Gb3 <- G[variant == "base" & schedule == "D3"]; Gs3 <- G[variant == "struct2020" & schedule == "D3"]
pb <- paired[variant == "base"]; ps <- paired[variant == "struct2020"]
pv <- function(p, sh, col) p[schedule == sh][[col]]
ci_txt <- function(p, sh, k, sep = ", ") sprintf("%s [%s%s%s]", s2(pv(p, sh, paste0("gain_", k, "_pp"))), f2(pv(p, sh, paste0("gain_", k, "_lo"))), sep, f2(pv(p, sh, paste0("gain_", k, "_hi"))))
bb <- b0s[variant == "base"]; bs <- b0s[variant == "struct2020"]
wci <- function(r, k) sprintf("%s%% [%s, %s]", f2(r[[paste0("reliable_", k, "_pct")]]), f2(r[[paste0("reliable_", k, "_lo")]]), f2(r[[paste0("reliable_", k, "_hi")]]))
wci_en <- function(r, k) sprintf("%s%% (95%% CI %s to %s)", f2(r[[paste0("reliable_", k, "_pct")]]), f2(r[[paste0("reliable_", k, "_lo")]]), f2(r[[paste0("reliable_", k, "_hi")]]))
any_c_i <- paired[schedule %in% DENSE, any(crit_c_i)]; any_c_ii <- paired[schedule %in% DENSE, any(crit_c_ii)]
rec_change <- paired[!is.na(recommend_i) & recommend_i != recommend_ii]
fs_net <- paired[variant %in% MECH_VARIANTS & schedule %in% DENSE, flag_span_change_pp]
# 두 모델 D1–D4: (ii) 증감은 언제나 (i)보다 낮다(= span 단독 손실이 늘어난다)
P2 <- paired[variant %in% TWO_MODELS & schedule %in% DENSE]
P2[, span_only_loss_change_pp := span_only_loss_pct - span_only_loss_B0_pct]
stopifnot(all(P2$gain_ii_pp < P2$gain_i_pp), all(P2$span_only_loss_change_pp > 0))
P124 <- P2[schedule != "D3"]
# 곡률(순간 반감기 < 1일, D-041) 구간 자체를 고른 경우는 드물다: 새 창 반감기 < 1일 비율
stopifnot(all(L$HL_lt1_pct < 5))
# 대조: B−(Day 50 삭제)는 다른 경로. 창 상한(Tlast)이 앞당겨져 가파른 끝 구간을 잃고 반감기가 길어진다(반감기 증가가 창 폭 감소보다 크다)
Lm <- trans[grp == "lost_span" & schedule == "Bminus"]
fs_bm <- paired[variant %in% MECH_VARIANTS & schedule == "Bminus", flag_span_change_pp]
stopifnot(nrow(Lm) == length(MECH_VARIANTS), all(Lm$d_upper_median < 0), all(Lm$HL_longer_pct > 50), all(Lm$mean_dlog_HL > 0),
          all(Lm$mean_dlog_HL > abs(Lm$mean_dlog_window)), all(fs_bm > 0))
Lm[, upper_earlier := -d_upper_median]

tbl_gain <- function(lang) {
  hdr <- if (lang == "ko") c("| 일정 | 2016 (i) | 2016 (ii) | Model 1 (i) | Model 1 (ii) |", "|---|---|---|---|---|") else c("| Schedule | 2016 model (i) | 2016 model (ii) | Model 1 (i) | Model 1 (ii) |", "|---|---|---|---|---|")
  f <- if (lang == "ko") function(p, sh, k) ci_txt(p, sh, k) else function(p, sh, k) ci_txt(p, sh, k, sep = " to ")
  c(hdr, vapply(setdiff(SCHED, "B0"), function(sh) sprintf("| %s | %s | %s | %s | %s |", sh, f(pb, sh, "i"), f(pb, sh, "ii"), f(ps, sh, "i"), f(ps, sh, "ii")), ""))
}
tbl_mech <- function(lang) {
  hdr <- if (lang == "ko") c("| 모델 | 일정 | 손실 n | λz 창 폭 중앙값 B0 → 일정 (일) | 짧아짐 % | 창 하한 이동 중앙값 (일) | 3점 창 % B0 → 일정 | 반감기 중앙값 B0 → 일정 (일) | 길어짐 % | adj R² 변화 중앙값 | 평균 Δlog 창 폭 / Δlog 반감기 / Δlog span |",
                             "|---|---|---|---|---|---|---|---|---|---|---|")
         else c("| Model | Schedule | Lost n | Window median, B0 to schedule (days) | Shorter % | Window start shift, median (days) | 3-point windows %, B0 to schedule | Half-life median, B0 to schedule (days) | Longer % | Adj R2 change, median | Mean change in log window / log half-life / log span |",
                "|---|---|---|---|---|---|---|---|---|---|---|")
  sl <- if (lang == "ko") SHORT_KO else SHORT_EN; arrow <- if (lang == "ko") " → " else " to "
  c(hdr, L[, sprintf("| %s | %s | %d | %s%s%s | %s | %s | %s%s%s | %s%s%s | %s | %s | %s / %s / %s |", sl[variant], schedule, n, f1(window_B0_median), arrow, f1(window_median), f1(window_shorter_pct),
                     s2(d_lower_median), f1(npts3_B0_pct), arrow, f1(npts3_pct), f2(HL_B0_median), arrow, f2(HL_median), f1(HL_longer_pct), sprintf("%+.3f", d_rsq_median),
                     f2(mean_dlog_window), f2(mean_dlog_HL), f2(mean_dlog_span))])
}
top_win <- function(v, sh) wins[variant == v & schedule == sh][1]
tw_b <- top_win("base", "D3"); tw_s <- top_win("struct2020", "D3")
win_ko <- function(w) sprintf("%s(%d점) → %s(%d점)", w$window_B0, w$npts_B0, w$window_sched, w$npts_sched)
win_en <- function(w) sprintf("%s (%d points) to %s (%d points)", w$window_B0, w$npts_B0, w$window_sched, w$npts_sched)
eng_ii <- function(v) f1(eng[model == v, reliable_pct])
ko <- c("# 신뢰 플래그 두 세트와 촘촘한 후기 채혈의 신뢰 충족률 기전 (D-039, D-044)", "",
  sprintf("자료: 저장된 개인 수준 NCA(`results/individual/nca_<변형>_20000.rds`), 모델별 가상 대상자 20,000명, B0 조건(60–90 kg, 채혈 허용창·잔차·BLQ), 일정 %s. 새 모의 없음. 생성: `scripts/39_reliability_flags.R`.", paste(SCHED, collapse = ", ")), "",
  "## 1. 두 플래그 세트의 B0 충족률", "",
  "- (i) = λz 산출 가능 & adj R² ≥ 0.80 & 외삽 ≤ 20% (D-039 이전 정의). (ii) = (i) & span ratio ≥ 2 (현행 reliable 열과 전원 일치).",
  sprintf("- 2016 모델: (i) %s, (ii) %s, span 단독 손실 %s%%p. Model 1: (i) %s, (ii) %s, 손실 %s%%p (Wilson 95%% 구간).",
          wci(bb, "i"), wci(bb, "ii"), f2(bb$span_only_loss_pct), wci(bs, "i"), wci(bs, "ii"), f2(bs$span_only_loss_pct)),
  sprintf("- 본문 두 모델 범위: (i) %s–%s%%, (ii) %s–%s%% ((ii) 범위는 `results/rationale/pillar1_two_model_range.csv`와 일치).", f1(rng$i_min), f1(rng$i_max), f1(rng$ii_min), f1(rng$ii_max)),
  sprintf("- 참고: `results/nca_engine/`의 B0 값((ii) %s%%, %s%%)은 같은 대상자에 B0 격자만으로 관측 시각·잔차를 따로 뽑은 별도 표본이라 위 값과 다르다.", eng_ii("base"), eng_ii("struct2020")), "",
  "## 2. B0 대비 쌍대 증감(%p, 같은 20,000명, 95% 구간)", "", tbl_gain("ko"), "",
  sprintf("- 두 모델의 D1–D4 모두 (ii) 증감이 (i)보다 낮다: span 단독 손실이 %s%%p 늘어난다. D1·D2·D4의 증감은 (i) %s%%p, (ii) %s%%p. D3는 (i) 2016 %s, Model 1 %s%%p, (ii) %s, %s%%p.",
          rtx(P2$span_only_loss_change_pp, f2), rtx(P124$gain_i_pp, f2), rtx(P124$gain_ii_pp, f2), s2(pv(pb, "D3", "gain_i_pp")), s2(pv(ps, "D3", "gain_i_pp")), s2(pv(pb, "D3", "gain_ii_pp")), s2(pv(ps, "D3", "gain_ii_pp"))),
  sprintf("- 기준 (c)(+%g%%p 이상)는 %s. 기준 (d) 값은 플래그와 무관해 두 세트에서 같다. (ii) 증감은 `schedule_decision_<변형>.csv`의 c_reliable_gain_pp와 0.01%%p 이내로 일치.",
          dr$reliability_gain_pp_min, if (!any_c_i && !any_c_ii) "두 세트 모두 어느 변형·추가 일정에서도 충족되지 않는다" else "일부 칸에서 충족된다(표 참조)"),
  sprintf("- 두 세트의 권고 판정이 달라지는 칸: %s.", if (nrow(rec_change)) paste(rec_change[, paste(variant, schedule)], collapse = ", ") else "없음"), "",
  "## 3. 기전: 왜 촘촘한 후기 채혈에서 span 플래그가 생기는가", "",
  "**판정: Best Fit가 더 늦고 짧은 λz 창(대개 추가 채혈일이 든 마지막 3점)을 고르기 때문이다. 반감기가 길어져서가 아니다.** 반감기는 오히려 짧아져 span 감소를 일부 상쇄한다. 새 창은 adj R²가 더 높아 Best Fit가 선택한다.", "",
  sprintf("- 구조: Day 22 이후 명목일로 만들 수 있는 가장 짧은 3점 창은 B0 %g일, 추가 일정(D1–D4) %g일이다. 3점 창에서 span ≥ 2가 되려면 반감기가 B0 %g일, 추가 일정 %g일 이하여야 한다.", m3_b0, m3_dense[1], m3_b0 / 2, m3_dense[1] / 2),
  sprintf("- 손실 대상자(B0 (ii) 충족 → 추가 일정 span 플래그, 기전 변형 4개 × D1–D4 = %d칸, 범위는 칸 사이): 창 폭이 짧아진 비율 %s%%, 창 폭 중앙값 B0 %s일 → %s일, 창 하한이 중앙값 %s일 늦어짐(늦어진 비율 %s%%), 3점 창 %s%% → %s%%.",
          nrow(L), rk("window_shorter_pct"), rk("window_B0_median"), rk("window_median"), rk("d_lower_median"), rk("lower_later_pct"), rk("npts3_B0_pct"), rk("npts3_pct")),
  sprintf("- 같은 대상자의 반감기 중앙값은 %s일 → %s일로 짧아졌고 길어진 비율은 %s%%뿐이다. 새 창에 추가 채혈일이 든 비율 %s%%, adj R²는 중앙값 %s 높아졌다(높아진 비율 %s%%).",
          rk("HL_B0_median", f2), rk("HL_median", f2), rk("HL_longer_pct"), rk("window_with_added_day_pct"), rtx(L$d_rsq_median, function(x) sprintf("%+.3f", x)), rk("rsq_higher_pct")),
  sprintf("- 로그 분해 log span = log 창 폭 − log 반감기: 평균 변화 창 폭 %s, 반감기 %s, span %s. 창 폭 감소가 span 감소의 %s%%에 해당하고(100%% 초과), 반감기 단축이 그 초과분을 상쇄한다.",
          rk("mean_dlog_window", f2), rk("mean_dlog_HL", f2), rk("mean_dlog_span", f2), rtx(100 * L$window_share_of_dlog_span, f1)),
  sprintf("- 예(D3, 가장 흔한 창 변화, 연구일): 2016 %s(손실의 %s%%), Model 1 %s(%s%%). 새 창 반감기 중앙값(%s일)은 %g일 창의 기준 %g일보다 길고 %g일 창의 기준 %g일보다 짧다 → 같은 반감기가 B0의 3점 창에서는 통과, 추가 일정의 3점 창에서는 플래그.",
          win_ko(tw_b), f1(tw_b$pct_of_lost), win_ko(tw_s), f1(tw_s$pct_of_lost), rk("HL_median", f2), m3_dense[1], m3_dense[1] / 2, m3_b0, m3_b0 / 2),
  sprintf("- 새 창의 적합 반감기 < 1일(D-041 절벽 정의의 순간 반감기 기준) 비율은 %s%%다. 고른 창은 절벽 자체가 아니라 그 직전의 가속 감소 구간이다%s.", rk("HL_lt1_pct"),
          if (is.na(cliff_p95)) "" else sprintf("(절벽 길이 95백분위 %s일, 두 모델 60–90 kg < 추가 일정의 Day 22 이후 최소 명목 간격 %g일 → 명목일 기준 3점 창이 절벽 안에 모두 들어갈 수 없다; `results/cliff/cliff_summary.csv`)", f2(cliff_p95), min_gap_dense)),
  sprintf("- 대조: B−(Day 50 삭제)의 span 플래그 증가(%s%%p)는 다른 경로다. 손실 대상자의 창 상한이 중앙값 %s일 앞당겨지고(가파른 끝 구간을 잃음) 반감기가 길어진 비율 %s%%, 평균 Δlog 반감기 %s가 Δlog 창 폭 %s보다 크다.",
          rtx(fs_bm, s2), rtx(Lm$upper_earlier, f1), rtx(Lm$HL_longer_pct, f1), rtx(Lm$mean_dlog_HL, s2), rtx(Lm$mean_dlog_window, f2)), "",
  tbl_mech("ko"), "",
  "## 4. 순효과는 모델에 따라 다르다", "",
  sprintf("- 반대 방향(B0 span 플래그 → 추가 일정 (ii) 충족, 회복)의 경로는 일부 다르다: 창 상한(Tlast)이 늦어진 비율 %s%%, 반감기가 짧아진 비율 %s%%. 평균 Δlog 창 폭은 Δlog span 증가의 %s%%이고, 반감기 단축의 기여가 창 폭 증가보다 큰 칸이 %d/%d칸이다.",
          rk("upper_later_pct", d = G), rk("HL_shorter_pct", d = G), rtx(100 * G$window_share_of_dlog_span, f1), G_hl_major, nrow(G)),
  sprintf("- 창이 짧아진다고 곧 플래그가 되지는 않는다: B0 (ii) 충족자의 %s%%가 추가 일정에서 가장 짧은 3점 창(명목 %g일, Day 22 이후)으로 바뀌었고, 이 가운데 span 플래그는 %s%%뿐이다. 그 창의 적합 반감기 중앙값이 %s일(같은 대상자의 B0 창 %s일)로 한계 %g일보다 짧기 때문이다. 이들이 손실 대상자의 %s%%다(`reliability_short_window_conversion.csv`). span 플래그는 새 창의 길이만이 아니라 그 창의 반감기에 달려 있고, 순변화는 손실과 회복의 차이다.",
          rtx(CV$short_pct_of_B0_reliable, f1), m3_dense[1], rtx(CV$span_flag_pct, f1), rtx(CV$HL_median, f2), rtx(CV$HL_B0_median, f2), m3_dense[1] / 2, rtx(CV$lost_in_short_pct, f1)),
  sprintf("- D3 순변화: 2016 손실 %s명·회복 %s명, span 플래그 %s%%p, (ii) %s%%p. Model 1 손실 %s명·회복 %s명, span 플래그 %s%%p, (ii) %s%%p. D1–D4 × 기전 변형 4개에서 span 플래그 순변화 %s%%p.",
          format(Lb3$n, big.mark = ","), format(Gb3$n, big.mark = ","), s2(pv(pb, "D3", "flag_span_change_pp")), s2(pv(pb, "D3", "gain_ii_pp")),
          format(Ls3$n, big.mark = ","), format(Gs3$n, big.mark = ","), s2(pv(ps, "D3", "flag_span_change_pp")), s2(pv(ps, "D3", "gain_ii_pp")), rtx(fs_net, s2)), "",
  "## 전제 검사(stopifnot)", "",
  sprintf("- 손실 대상자(D1–D4): 창 폭 짧아짐 ≥ 90%%, 창 하한이 늦어짐 ≥ 90%%, 3점 창 비율 증가(≥ 75%%), 반감기 중앙값 감소와 길어짐 < 25%%, 평균 Δlog 창 폭 < 평균 Δlog span < 0, adj R² 중앙값 증가(높아진 비율 > 50%%), 추가 채혈일 포함 ≥ 90%%, 새 창 반감기 중앙값 %g–%g일 사이, 새 창 적합 반감기 < 1일 < 5%%, 명목일 대응 일치 ≥ 95%%. 회복 대상자: 창 상한이 늦어짐 ≥ 50%%, 창 폭 중앙값 증가, 반감기 짧아짐 > 50%%, 평균 Δlog 반감기 < 0. 최단 3점 창으로 바뀐 B0 충족자: span 플래그 < 50%%, 반감기 중앙값 < %g일이고 B0 창보다 짧음, 3점 창 ≥ 99%%. B−: 창 상한 중앙값 감소, 반감기 길어짐 > 50%%, Δlog 반감기 > |Δlog 창 폭|. 두 모델 D1–D4: (ii) 증감 < (i) 증감. 절벽 길이 95백분위 < 추가 일정 최소 명목 간격.", m3_dense[1] / 2, m3_b0 / 2, m3_dense[1] / 2),
  "- (ii) = 저장된 reliable, (i) = λz 산출 & 두 플래그 없음, 개인 수준 표·판정표·두 모델 범위와 수치 일치.")
en <- c("# Two reliability flag sets and why denser late sampling lowers reliability (D-039, D-044)", "",
  sprintf("Data: saved subject-level NCA (`results/individual/nca_<variant>_20000.rds`), 20,000 virtual subjects per model, B0 conditions (60 to 90 kg, visit windows, residual error, BLQ), schedules %s. No new simulation. Generated by `scripts/39_reliability_flags.R`.", paste(SCHED, collapse = ", ")), "",
  "## 1. Reliability at B0 under the two flag sets", "",
  "- (i) = lambda-z estimable, adjusted R2 at least 0.80 and extrapolation at most 20% (definition before D-039). (ii) = (i) plus span ratio at least 2 (identical to the stored reliable flag for every subject).",
  sprintf("- 2016 model: (i) %s, (ii) %s, loss from the span ratio alone %s percentage points. Model 1: (i) %s, (ii) %s, loss %s points (Wilson intervals).",
          wci_en(bb, "i"), wci_en(bb, "ii"), f2(bb$span_only_loss_pct), wci_en(bs, "i"), wci_en(bs, "ii"), f2(bs$span_only_loss_pct)),
  sprintf("- Two-model range used in the texts: (i) %s to %s%%, (ii) %s to %s%% (the (ii) range matches `results/rationale/pillar1_two_model_range.csv`).", f1(rng$i_min), f1(rng$i_max), f1(rng$ii_min), f1(rng$ii_max)),
  sprintf("- Note: the B0 values in `results/nca_engine/` ((ii) %s%% and %s%%) come from a separate draw of observation times and residuals for the same subjects on the B0-only grid, so they differ from the values above.", eng_ii("base"), eng_ii("struct2020")), "",
  "## 2. Paired change versus B0 (percentage points, same 20,000 subjects, 95% CI)", "", tbl_gain("en"), "",
  sprintf("- In both models and all of D1 to D4 the (ii) change is below the (i) change: the loss from the span ratio alone grows by %s points. For D1, D2 and D4 the change is %s points under (i) and %s points under (ii). For D3 it is %s (2016 model) and %s (Model 1) under (i), %s and %s under (ii).",
          rtx(P2$span_only_loss_change_pp, f2, "en"), rtx(P124$gain_i_pp, f2, "en"), rtx(P124$gain_ii_pp, f2, "en"),
          s2(pv(pb, "D3", "gain_i_pp")), s2(pv(ps, "D3", "gain_i_pp")), s2(pv(pb, "D3", "gain_ii_pp")), s2(pv(ps, "D3", "gain_ii_pp"))),
  sprintf("- Criterion (c) (gain of at least %g points) is %s. The criterion (d) quantities do not depend on the flags and are the same under both sets. The (ii) gains match c_reliable_gain_pp in `schedule_decision_<variant>.csv` within 0.01 points.",
          dr$reliability_gain_pp_min, if (!any_c_i && !any_c_ii) "met under neither flag set for any variant or added-sampling schedule" else "met in some cells (see table)"),
  sprintf("- Cells where the recommendation differs between the two flag sets: %s.", if (nrow(rec_change)) paste(rec_change[, paste(variant, schedule)], collapse = ", ") else "none"), "",
  "## 3. Mechanism: why denser late sampling adds span flags", "",
  "**Verdict: Best Fit picks a later and shorter lambda-z window (usually the last 3 points, including added sampling days). The half-life does not get longer.** It gets shorter, which partly offsets the loss of span. The new window has a higher adjusted R2, which is why Best Fit selects it.", "",
  sprintf("- Structure: the shortest 3-point window that the nominal days after Day 22 allow is %g days at B0 and %g days in the added-sampling schedules (D1 to D4). A 3-point window reaches span ratio 2 only if the half-life is at most %g days at B0 and %g days with added sampling.", m3_b0, m3_dense[1], m3_b0 / 2, m3_dense[1] / 2),
  sprintf("- Lost subjects (reliable under (ii) at B0, span-flagged in the added-sampling schedule; 4 model variants x D1 to D4 = %d cells, ranges across cells): the window got shorter in %s%%, median window %s days at B0 versus %s days, window start later by a median %s days (later in %s%%), 3-point windows %s%% versus %s%%.",
          nrow(L), re("window_shorter_pct"), re("window_B0_median"), re("window_median"), re("d_lower_median"), re("lower_later_pct"), re("npts3_B0_pct"), re("npts3_pct")),
  sprintf("- In the same subjects the median half-life was %s days at B0 and %s days with added sampling (shorter, not longer); it got longer in only %s%% of subjects. The new window contains an added sampling day in %s%%, and its adjusted R2 was higher by a median %s (higher in %s%%).",
          re("HL_B0_median", f2), re("HL_median", f2), re("HL_longer_pct"), re("window_with_added_day_pct"), rtx(L$d_rsq_median, function(x) sprintf("%+.3f", x), "en"), re("rsq_higher_pct")),
  sprintf("- Log decomposition, log span = log window minus log half-life: mean change %s for the window, %s for the half-life, %s for the span. The window shortening equals %s%% of the span decrease (more than 100%%); the shorter half-life offsets the excess.",
          re("mean_dlog_window", f2), re("mean_dlog_HL", f2), re("mean_dlog_span", f2), rtx(100 * L$window_share_of_dlog_span, f1, "en")),
  sprintf("- Example (D3, most common window change, study days): 2016 model %s (%s%% of lost subjects), Model 1 %s (%s%%). The median half-life of the new window (%s days) is above the %g-day limit of a %g-day window and below the %g-day limit of a %g-day window, so the same half-life passes with a 3-point window at B0 and is flagged with a 3-point window under added sampling.",
          win_en(tw_b), f1(tw_b$pct_of_lost), win_en(tw_s), f1(tw_s$pct_of_lost), re("HL_median", f2), m3_dense[1] / 2, m3_dense[1], m3_b0 / 2, m3_b0),
  sprintf("- A fitted half-life below 1 day (the instantaneous half-life of the D-041 cliff definition) occurs in %s%% of the new windows: the selected window is the accelerating decline before the cliff, not the cliff itself%s.", re("HL_lt1_pct"),
          if (is.na(cliff_p95)) "" else sprintf(" (the 95th percentile of cliff length, %s days in both models at 60 to 90 kg, is below the %g-day minimum nominal interval after Day 22 in the added-sampling schedules, so a 3-point window cannot lie entirely in the cliff at nominal days; `results/cliff/cliff_summary.csv`)", f2(cliff_p95), min_gap_dense)),
  sprintf("- Contrast: the increase in span flags with B minus (Day 50 removed; %s points) takes a different route. In lost subjects the window end moves a median %s days earlier (the steep final segment is lost), the half-life gets longer in %s%%, and the mean change in log half-life (%s) exceeds the change in log window (%s).",
          rtx(fs_bm, s2, "en"), rtx(Lm$upper_earlier, f1, "en"), rtx(Lm$HL_longer_pct, f1, "en"), rtx(Lm$mean_dlog_HL, s2, "en"), rtx(Lm$mean_dlog_window, f2, "en")), "",
  tbl_mech("en"), "",
  "## 4. The net effect depends on the model", "",
  sprintf("- The reverse flow (span-flagged at B0, reliable under (ii) with added sampling) takes a partly different route: the window end (Tlast) moved later in %s%% and the half-life got shorter in %s%%. The mean change in log window is %s%% of the rise in log span, and the shorter half-life contributes more than the longer window in %d of %d cells.",
          re("upper_later_pct", d = G), re("HL_shorter_pct", d = G), rtx(100 * G$window_share_of_dlog_span, f1, "en"), G_hl_major, nrow(G)),
  sprintf("- A shorter window alone does not flag a subject: %s%% of the subjects reliable under (ii) at B0 moved to a shortest 3-point window (nominal %g days, after Day 22) with added sampling, and only %s%% of these were span-flagged, because the half-life fitted in that window had a median of %s days (%s days in the same subjects' B0 window), below the %g-day limit. They are %s%% of the lost subjects (`reliability_short_window_conversion.csv`). The span flag therefore depends on the half-life in the new window, not only on its length, and the net change is the difference between the lost and regained flows.",
          rtx(CV$short_pct_of_B0_reliable, f1, "en"), m3_dense[1], rtx(CV$span_flag_pct, f1, "en"), rtx(CV$HL_median, f2, "en"), rtx(CV$HL_B0_median, f2, "en"), m3_dense[1] / 2, rtx(CV$lost_in_short_pct, f1, "en")),
  sprintf("- Net change at D3: 2016 model %s lost and %s regained, span flag %s points, (ii) %s points. Model 1 %s lost and %s regained, span flag %s points, (ii) %s points. Across D1 to D4 and the 4 model variants the net change in span flags is %s points.",
          format(Lb3$n, big.mark = ","), format(Gb3$n, big.mark = ","), s2(pv(pb, "D3", "flag_span_change_pp")), s2(pv(pb, "D3", "gain_ii_pp")),
          format(Ls3$n, big.mark = ","), format(Gs3$n, big.mark = ","), s2(pv(ps, "D3", "flag_span_change_pp")), s2(pv(ps, "D3", "gain_ii_pp")), rtx(fs_net, s2, "en")), "",
  "## Premise checks (stopifnot)", "",
  sprintf("- Lost subjects (D1 to D4): window shorter in at least 90%%, window start later in at least 90%%, more 3-point windows (at least 75%%), lower median half-life and longer half-life in under 25%%, mean change in log window below mean change in log span below 0, higher median adjusted R2 (higher in over 50%%), added sampling day in the new window in at least 90%%, median half-life of the new window between %g and %g days, fitted half-life below 1 day in under 5%%, nominal-day mapping consistent in at least 95%%. Regained subjects: window end later in at least 50%%, longer median window, shorter half-life in over 50%% and negative mean change in log half-life. B0-reliable subjects moved to a shortest 3-point window: span-flagged in under 50%%, median half-life below %g days and below their B0 window, 3-point windows in at least 99%%. B minus: lower median window end, longer half-life in over 50%%, change in log half-life larger than the absolute change in log window. Both models, D1 to D4: (ii) change below (i) change. 95th percentile of cliff length below the minimum nominal interval of the added-sampling schedules.", m3_dense[1] / 2, m3_b0 / 2, m3_dense[1] / 2),
  "- (ii) equals the stored reliable flag, (i) equals lambda-z estimable with neither of the other two flags, and the numbers match the subject-level tables, the decision tables and the two-model range.")
# 영문: 한글·em dash·en dash 없음
stopifnot(!any(grepl("[ᄀ-ᇿ㄰-㆏가-힣]", en)), !any(grepl("—|–", en)))
writeLines(ko, file.path(out_dir, "mechanism_summary_ko.md")); writeLines(en, file.path(out_dir, "mechanism_summary_en.md"))

# ---------------------------------------------------------------------------------------------
# 5. 그림: 위 = λz 창 폭·span ratio 누적분포(일정별, 색 + 선 모양), 아래 = 손실 대상자의 창 폭·반감기 중앙값 B0 → 추가 일정
SCHED_COL <- c(B0 = VIZ$ink, Bminus = VIZ$muted, D1 = VIZ$s1, D2 = VIZ$s3, D3 = VIZ$s2, D4 = "#4a3aa7")   # 참조 팔레트 1–3·7번, 전쌍 검증 통과
SCHED_LT <- c(B0 = "solid", Bminus = "dotted", D1 = "dashed", D2 = "longdash", D3 = "twodash", D4 = "dotdash")
SCHED_SH <- c(B0 = 16, Bminus = 1, D1 = 17, D2 = 15, D3 = 18, D4 = 25)
sched_lab <- function(lang) {
  vapply(SCHED, function(sh) {
    if (sh == "B0") return(if (lang == "ko") "B0 (현행)" else "B0 (current)")
    rm_ <- setdiff(B0_DAYS, days_of(sh)) + 1; ad <- added_days(sh) + 1
    if (length(rm_)) return(if (lang == "ko") sprintf("B− (Day %s 삭제)", paste(rm_, collapse = "·")) else sprintf("B minus (Day %s removed)", paste(rm_, collapse = ", ")))
    sprintf("%s (+Day %s)", sh, paste(ad, collapse = if (lang == "ko") "·" else ", "))
  }, "")
}
TX <- list(
  ko = list(model = c(base = "Kovalenko 2016 (주)", struct2020 = "Kovalenko 2020 Model 1"), window = "λz 창 폭 (일)", span = "Span ratio (창 폭 / 반감기)",
            y_top = "누적 비율 (λz 산출 가능 대상자)", thr = "span 2 (플래그 기준)", m3 = c("7일: 추가 일정의\n최단 3점 창", "14일: B0의\n최단 3점 창"),
            share = "span ratio < 2 비율", sched_short = c(B0 = "B0", Bminus = "B−", D1 = "D1", D2 = "D2", D3 = "D3", D4 = "D4"),
            t_top = "A. 일정별 λz 창 폭(Lambda_z_upper − Lambda_z_lower)과 span ratio 분포",
            s_top = "모델별 가상 대상자 20,000명, B0 조건(60–90 kg, 채혈 허용창·잔차·BLQ), λz 산출 가능 대상자, 저장된 NCA.\nSpan ratio는 0–3.5 구간만 표시(비율의 분모는 λz 산출 가능 전원).",
            t_bot = "B. B0에서 신뢰 충족 → 추가 일정에서 span 플래그: 같은 대상자의 중앙값 변화",
            s_bot = "B0 값(원)에서 추가 일정 값(색 삼각형)으로.\n창 폭은 짧아지고 반감기도 짧아진다(길어지지 않음).",
            bot_metric = c(window = "λz 창 폭 중앙값 (일)", HL = "반감기 HL_Lambda_z 중앙값 (일)"), stage = c(B0 = "B0 중앙값", sched = "추가 일정 중앙값"), n = "%s (손실 %s명)",
            title = "그림. 촘촘한 후기 채혈에서 span 플래그가 생기는 기전: Best Fit가 더 짧은 λz 창을 고른다",
            cap = "플래그 세트 (ii) = λz 산출 & adj R² ≥ 0.80 & 외삽 ≤ 20% & span ratio ≥ 2.\n수치: results/reliability/reliability_span_transition.csv, reliability_lz_window_by_schedule.csv"),
  en = list(model = c(base = "Kovalenko 2016 (primary)", struct2020 = "Kovalenko 2020 Model 1"), window = "Lambda-z window length (days)", span = "Span ratio (window length / half-life)",
            y_top = "Cumulative share (subjects with estimable lambda-z)", thr = "span 2 (flag limit)", m3 = c("7 days: shortest 3-point\nwindow, added sampling", "14 days: shortest\n3-point window, B0"),
            share = "Share with span ratio < 2", sched_short = c(B0 = "B0", Bminus = "B minus", D1 = "D1", D2 = "D2", D3 = "D3", D4 = "D4"),
            t_top = "A. Lambda-z window length (Lambda_z_upper minus Lambda_z_lower) and span ratio by schedule",
            s_top = "20,000 virtual subjects per model, B0 conditions (60 to 90 kg, visit windows, residual error, BLQ),\nsubjects with estimable lambda-z, saved NCA. Span ratio shown from 0 to 3.5 only (shares use all subjects with estimable lambda-z).",
            t_bot = "B. Reliable at B0, span-flagged with added sampling: change in medians of the same subjects",
            s_bot = "From the B0 value (circle) to the added-sampling value (coloured triangle).\nThe window gets shorter and so does the half-life (it does not get longer).",
            bot_metric = c(window = "Median lambda-z window length (days)", HL = "Median half-life HL_Lambda_z (days)"), stage = c(B0 = "B0 median", sched = "Added-sampling median"), n = "%s (lost %s)",
            title = "Figure. Why denser late sampling adds span flags: Best Fit picks a shorter lambda-z window",
            cap = "Flag set (ii) = lambda-z estimable, adjusted R2 at least 0.80, extrapolation at most 20%, span ratio at least 2.\nNumbers: results/reliability/reliability_span_transition.csv, reliability_lz_window_by_schedule.csv"))
stopifnot(!any(grepl("[ᄀ-ᇿ㄰-㆏가-힣]|—|–", unlist(TX$en))))
# 누적분포는 격자에서 계산(표시 범위 밖 값도 분모에 포함)
grid_w <- seq(0, 45, by = 0.1); grid_s <- seq(0, 3.5, by = 0.01)
ecdf_d <- rbindlist(lapply(TWO_MODELS, function(v) NCA[[v]][lambda_ok == TRUE, {
  fw <- ecdf(window); fs <- ecdf(Span_ratio)
  rbind(data.table(metric = "window", x = grid_w, y = fw(grid_w)), data.table(metric = "span", x = grid_s, y = fs(grid_s)))
}, by = schedule][, variant := v]))
share_d <- lz[variant %in% TWO_MODELS, .(variant, schedule, share = span_flag_pct_of_lambda_ok)]
bot_d <- rbindlist(fig_pairs)
bot_l <- rbind(bot_d[, .(variant, schedule, n, metric = "window", B0 = window_B0_median, sched = window_median)],
               bot_d[, .(variant, schedule, n, metric = "HL", B0 = HL_B0_median, sched = HL_median)])
for (lg in names(TX)) {
  T_ <- TX[[lg]]; sfx <- if (lg == "en") "_en" else ""; SL <- sched_lab(lg)
  panel_lab <- function(v, m) paste0(T_$model[v], "\n", c(window = T_$window, span = T_$span)[m])
  plev <- as.vector(t(outer(TWO_MODELS, c("window", "span"), function(v, m) panel_lab(v, m))))
  e <- copy(ecdf_d)
  e[, `:=`(panel = factor(panel_lab(variant, metric), levels = plev), sched = factor(SL[schedule], levels = SL))]
  thr <- data.table(panel = factor(panel_lab(TWO_MODELS, "span"), levels = plev), x = rel$span_ratio_min)
  m3d <- CJ(v = TWO_MODELS, k = 1:2)[, `:=`(panel = factor(panel_lab(v, "window"), levels = plev), x = c(m3_dense[1], m3_b0)[k], y = c(0.98, 0.13)[k], lab = T_$m3[k])]
  sh <- share_d[, .(txt = paste0(T_$share, ":\n", paste(sprintf("%s %s%%", T_$sched_short[schedule], f1(share)), collapse = ", "))), by = variant]
  sh[, `:=`(panel = factor(panel_lab(variant, "span"), levels = plev), x = 0, y = Inf)]
  sh[, txt := sub(", D1", "\nD1", txt)]                                          # 두 줄로
  col <- setNames(SCHED_COL[SCHED], SL); lt <- setNames(SCHED_LT[SCHED], SL)
  gA <- ggplot(e, aes(x, y, colour = sched, linetype = sched)) +
    geom_vline(data = m3d, aes(xintercept = x), colour = VIZ$muted, linewidth = 0.35, linetype = "13") +
    geom_text(data = m3d, aes(x = x, y = y, label = lab), inherit.aes = FALSE, hjust = -0.04, vjust = 1, size = 2.5, colour = VIZ$ink2, lineheight = 0.9) +
    geom_vline(data = thr, aes(xintercept = x), colour = VIZ$ink, linewidth = 0.45, linetype = "22") +
    geom_text(data = thr, aes(x = x, y = 0, label = T_$thr), inherit.aes = FALSE, hjust = -0.05, vjust = 0, size = 2.5, colour = VIZ$ink) +
    geom_text(data = sh, aes(x = x, y = y, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 2.5, colour = VIZ$ink2, lineheight = 0.95) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = col, name = NULL) + scale_linetype_manual(values = lt, name = NULL) +
    scale_y_continuous(labels = function(z) paste0(round(100 * z), "%")) +
    facet_wrap(~panel, ncol = 2, scales = "free") + labs(x = NULL, y = T_$y_top, title = T_$t_top, subtitle = T_$s_top) +
    theme_dupi() + theme(legend.key.width = grid::unit(2.2, "lines")) + guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
  b <- copy(bot_l)
  b[, `:=`(model = factor(T_$model[variant], levels = T_$model), metric = factor(T_$bot_metric[metric], levels = T_$bot_metric),
           ylab = sprintf(T_$n, schedule, prettyNum(n, big.mark = ",")))]
  b[, ylab := factor(ylab, levels = rev(unique(ylab[order(match(variant, TWO_MODELS), match(schedule, DENSE))])))]
  b[, schedule_f := factor(SL[schedule], levels = SL)]
  gB <- ggplot(b, aes(y = ylab)) +
    geom_segment(aes(x = B0, xend = sched, yend = ylab), colour = VIZ$muted, linewidth = 0.5, arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")) +
    geom_point(aes(x = B0, shape = T_$stage[["B0"]]), colour = VIZ$ink, fill = VIZ$surface, size = 2.4, stroke = 0.8) +
    geom_point(aes(x = sched, colour = schedule_f, shape = T_$stage[["sched"]]), size = 2.6) +
    geom_text(aes(x = B0, label = f1(B0)), vjust = -1, size = 2.5, colour = VIZ$ink2) +
    geom_text(aes(x = sched, label = f1(sched)), vjust = -1, size = 2.5, colour = VIZ$ink2) +
    scale_colour_manual(values = col, guide = "none", drop = FALSE) +
    scale_shape_manual(values = setNames(c(21, 17), T_$stage), breaks = unname(T_$stage), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    facet_grid(model ~ metric, scales = "free") + labs(x = NULL, y = NULL, title = T_$t_bot, subtitle = T_$s_bot) + theme_dupi()
  p <- gA / gB + plot_layout(heights = c(1.2, 1)) +
    plot_annotation(title = T_$title, caption = T_$cap, theme = theme_dupi() + theme(plot.caption = element_text(colour = VIZ$ink2, size = 8, hjust = 0)))
  ggsave(file.path(out_dir, sprintf("fig_reliability_mechanism%s.png", sfx)), p, width = 12, height = 13, dpi = 120)
}
append_run_log(logfile, "done")
cat(ko, sep = "\n")
