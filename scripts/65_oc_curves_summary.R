#!/usr/bin/env Rscript
# 65_oc_curves_summary.R — S9 재설계 요약: 판정 구성별 통과 확률 곡선, 1종·2종 오류, 쌍대 차이, 2종 오류 분해, 간이 계산 대조, 문구 규칙
# (지시 2026-09-29 "S9 재설계: 1종·2종 오류", config/prereg_20260929_oc.yaml section8, D-064).
# 입력: 경계·동일 제품 = results/oc_models(oc_models_be/crit, 연장 칸은 oc_models_ext_be/ext_crit), 그 밖의 칸 = results/oc_curves(scripts/64).
# 검사(다르면 중단): 경계 P2·F3B·G2_B(M1)를 type1_models.csv, G2_A_iii(M1)를 criteria_g2_type1.csv, 동일 제품 P2·F3B·G2_B(M1)를 power_models.csv와
#   (시험 수, 통과 수) 정확히 대조. 새 칸의 시험 수가 section8 반복 수와 같은지.
# 산출(results/oc_curves): oc_curves_pass.csv, oc_type1_summary.csv, oc_type2_summary.csv, oc_paired.csv, oc_decomposition.csv,
#   oc_analytic_check.csv, oc_wording_rules.csv, oc_identity_checks.csv, oc_bias_relation.csv(별첨 A5d 상자, 사전 등록 뒤 추가한 기술 요약)
source("R/00_setup.R"); source_project()
s8 <- read_cfg("prereg_20260929_oc.yaml")$section8; design <- read_cfg("trial_design.yaml"); oc <- read_cfg("oc_design.yaml")
args <- commandArgs(trailingOnly = TRUE); partial <- "--partial" %in% args          # --partial <dir>: 시험이 덜 끝난 상태의 시험용(시험 수 검사 생략, 결과는 <dir>에, 실행 기록 없음)
out_dir <- if (partial) args[which(args == "--partial") + 1L] else proj_path(s8$output_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
MODELS <- unlist(s8$pk_models); AM <- s8$analysis_model; stopifnot(AM == "M1")
CFG <- lapply(s8$configurations, function(x) unlist(x$endpoints)); CFG_N <- names(CFG)
MAIN <- c("P2", "F3A_iii", "G2A_iii"); SENS <- c("F3B", "G2B")
EPN <- unique(unlist(CFG)); stopifnot(setequal(EPN, c("AUClast", "Cmax", "AUCinf_Aiii", "AUCinf_B", "AUCinf_true")))
lims <- c(0.80, 1.25); alpha_pct <- 100 * (1 - design$be$ci_level) / 2; stopifnot(abs(alpha_pct - 5) < 1e-9)
bnd <- as.numeric(unlist(oc$boundary_targets)); n_near <- as.integer(s8$reps$near_one$trials); n_other <- as.integer(s8$reps$other$trials)
near_t <- as.numeric(unlist(s8$reps$near_one$targets)); other_t <- as.numeric(unlist(s8$reps$other$targets))
logfile <- if (partial) NULL else start_run_log("oc_curves_summary", master_seed = as.integer(oc$trials$master_seed), run_mode = "final", extra = list(models = paste(MODELS, collapse = ",")))
say <- function(msg) { cat(msg, "\n"); if (!is.null(logfile)) append_run_log(logfile, msg) }
on_t <- function(x, tt) vapply(x, function(v) any(abs(v - tt) < 1e-9), logical(1))

# ---- 칸 표 ----------------------------------------------------------------------------------------------------------
inv_all <- fread(proj_path("results", "oc", "inversion_all.csv"))
cells_of <- function(m) {
  inv <- inv_all[model == m & reachable == TRUE & mechanism %in% unlist(s8$cells$mechanisms) & on_t(target, as.numeric(unlist(s8$cells$targets))) & abs(target - 1) > 1e-9]
  inv[, code := sprintf("%s_%s_%03d", mechanism, direction, round(100 * target))]
  if (anyDuplicated(inv[, .(mechanism, target)])) stop("한 기전·목표에 도달 방향이 둘입니다: ", m)
  s0 <- data.table(code = "S00", mechanism = "identity", direction = "none", target = 1, multiplier = 1, auc_ratio = 1, cmax_ratio = 1)
  x <- rbind(s0, inv[, .(code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio)])
  x[, kind := fifelse(code == "S00", "identity", fifelse(on_t(target, bnd), "boundary", "curve"))]
  x[, cmax_in_limits := cmax_ratio >= lims[1] & cmax_ratio <= lims[2]]
  unreach <- inv_all[model == m & reachable == FALSE & mechanism %in% unlist(s8$cells$mechanisms) & on_t(target, as.numeric(unlist(s8$cells$targets))) & abs(target - 1) > 1e-9]
  unreach <- unreach[!paste(mechanism, target) %in% inv[, paste(mechanism, target)], .(mechanism, target)]
  list(cells = x[, pk_model := m][], unreachable = unique(unreach)[, pk_model := m][])
}
CL <- lapply(setNames(MODELS, MODELS), cells_of)

# ---- 시험별 구성 통과 -------------------------------------------------------------------------------------------------
rd <- function(f, scen) { if (!file.exists(f)) stop("입력 없음: ", f); x <- fread(f)[scenario %in% scen & endpoint %in% EPN & model %in% c("M0", "M1")]; x }
per_trial <- function(m) {
  cl <- CL[[m]]$cells; old <- cl[kind != "curve", code]; new <- cl[kind == "curve", code]
  om <- proj_path("results", "oc_models"); ocv <- if ("--input" %in% args) args[which(args == "--input") + 1L] else proj_path("results", "oc_curves")   # --input: 시험용 사본 폴더
  be <- rbind(rd(file.path(om, sprintf("oc_models_be_%s.csv.gz", m)), old), rd(file.path(om, sprintf("oc_models_crit_be_%s.csv.gz", m)), old),
              rd(file.path(ocv, sprintf("oc_curves_be_%s.csv.gz", m)), new), rd(file.path(ocv, sprintf("oc_curves_crit_be_%s.csv.gz", m)), new))
  for (f in c(sprintf("oc_models_ext_be_%s.csv.gz", m), sprintf("oc_models_ext_crit_be_%s.csv.gz", m))) if (file.exists(file.path(om, f))) be <- rbind(be, rd(file.path(om, f), old))
  if (anyDuplicated(be, by = c("trial", "scenario", "endpoint", "model"))) stop("중복 행: ", m)
  w <- dcast(be, trial + scenario + model ~ endpoint, value.var = "pass", fun.aggregate = function(v) as.integer(v %in% TRUE), fill = NA_integer_)
  if (anyNA(w[, ..EPN])) stop("평가변수가 빠진 시험: ", m)
  for (cf in CFG_N) w[, (cf) := as.integer(Reduce(`&`, lapply(CFG[[cf]], function(e) get(e) == 1L)))]
  w[, AUClast_only := AUClast][, AUCinf_Aiii_only := AUCinf_Aiii]
  nn <- dcast(be[endpoint %in% c("AUCinf_Aiii", "AUCinf_B", "AUClast")], trial + scenario + model ~ endpoint, value.var = c("n_R", "n_T"))
  ses <- be[endpoint %in% c("AUClast", "AUCinf_Aiii") & model == "M0", .(trial, scenario, endpoint, se, n_R, n_T)]
  est <- be[model == AM & scenario %in% cl[kind == "boundary", code] & endpoint %in% c("AUClast", "AUCinf_Aiii", "AUCinf_B", "AUCinf_true"), .(trial, scenario, endpoint, est, se, n_R, n_T, pass)]
  list(w = w[, pk_model := m][], n = nn[, pk_model := m][], se = ses[, pk_model := m][], est = est[, pk_model := m][])
}
PT <- lapply(setNames(MODELS, MODELS), per_trial)

# ---- 칸별 통과 확률 ----------------------------------------------------------------------------------------------------
ALLC <- c(CFG_N, "AUClast_only", "AUCinf_Aiii_only")
pass_tab <- rbindlist(lapply(MODELS, function(m) {
  w <- PT[[m]]$w; cl <- CL[[m]]$cells
  lp <- melt(w, id.vars = c("pk_model", "trial", "scenario", "model"), measure.vars = ALLC, variable.name = "config", value.name = "ok")
  s <- lp[, .(n_trials = .N, n_pass = sum(ok)), by = .(pk_model, scenario, analysis_model = model, config)]
  s[, c("pass_pct", "lo", "hi") := wilson_ci(n_pass, n_trials)]
  merge(cl, s, by.x = c("pk_model", "code"), by.y = c("pk_model", "scenario"))
}))
pass_tab[, class := fifelse(lo > alpha_pct, "exceeding", fifelse(hi < alpha_pct, "conservative", "nominal"))]

# 시험 수 검사
chk_n <- pass_tab[analysis_model == AM & config == "P2"]
exp_n <- function(kind, target, pk, code) fifelse(kind == "curve", fifelse(on_t(target, near_t), n_near, n_other), NA_integer_)
chk_n[, expected := exp_n(kind, target)]
if (!partial && any(chk_n[kind == "curve", n_trials != expected])) stop("새 칸의 시험 수가 section8과 다릅니다")
if (any(chk_n[kind == "identity", n_trials != 10000L])) stop("동일 제품 칸의 시험 수가 10,000이 아닙니다")

# ---- 저장 요약과 대조 --------------------------------------------------------------------------------------------------
t1 <- fread(proj_path("results", "oc_models", "type1_models.csv"))[analysis_model == AM & config %in% c("P2", "F3B", "G2_B")]
t1[, config := fifelse(config == "G2_B", "G2B", config)]
cg <- fread(proj_path("results", "criteria", "criteria_g2_type1.csv"))[analysis_model == AM & config == "G2_A_iii"][, config := "G2A_iii"]
cg[, n_pass := round(pass_pct * n_trials / 100)]
pw <- fread(proj_path("results", "oc_models", "power_models.csv"))[analysis_model == AM & scenario == "S00" & config %in% c("P2", "F3B", "G2_B")]
pw[, config := fifelse(config == "G2_B", "G2B", config)]
ref <- rbind(t1[, .(pk_model, scenario, config, n_trials, n_pass, source = "type1_models.csv")], cg[, .(pk_model, scenario, config, n_trials, n_pass, source = "criteria_g2_type1.csv")],
             pw[, .(pk_model, scenario, config, n_trials, n_pass, source = "power_models.csv")])
idc <- merge(ref, pass_tab[analysis_model == AM, .(pk_model, scenario = code, config, n_trials_new = n_trials, n_pass_new = n_pass)], by = c("pk_model", "scenario", "config"), all.x = TRUE)
idc[, equal := n_trials == n_trials_new & n_pass == n_pass_new]
fwrite(idc, file.path(out_dir, "oc_identity_checks.csv"))
if (!all(idc$equal %in% TRUE)) { print(idc[!equal %in% TRUE]); stop("저장된 요약과 다릅니다") }
say(sprintf("identity checks: %d stored summaries (type1_models, criteria_g2_type1, power_models) reproduced exactly", nrow(idc)))

pm <- pass_tab[analysis_model == AM]
fwrite(pass_tab[, .(pk_model, code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio, cmax_in_limits, kind, analysis_model, config, n_trials, n_pass, pass_pct, lo, hi, class)],
       file.path(out_dir, "oc_curves_pass.csv"))

# ---- 1종 오류(경계 16칸) ---------------------------------------------------------------------------------------------
b <- pm[kind == "boundary"]
type1 <- rbind(b[, .(scope = "both", n_cells = .N, max_pct = max(pass_pct), max_cell = paste(pk_model, code)[which.max(pass_pct)], n_gt5 = sum(pass_pct > alpha_pct),
                     n_wilson_gt5 = sum(lo > alpha_pct), min_pct = min(pass_pct)), by = config],
               b[, .(scope = pk_model[1], n_cells = .N, max_pct = max(pass_pct), max_cell = paste(pk_model, code)[which.max(pass_pct)], n_gt5 = sum(pass_pct > alpha_pct),
                     n_wilson_gt5 = sum(lo > alpha_pct), min_pct = min(pass_pct)), by = .(config, pk_model)][, pk_model := NULL])
fwrite(type1, file.path(out_dir, "oc_type1_summary.csv"))

# ---- 2종 오류 ---------------------------------------------------------------------------------------------------------
pm[, type2_pct := 100 - pass_pct][, type2_lo := 100 - hi][, type2_hi := 100 - lo]
grp <- function(target) fifelse(abs(target - 1) < 1e-9, "identity", fifelse(on_t(target, c(0.95, 1.05)), "near_095_105", fifelse(on_t(target, c(0.90, 1.11)), "mid_090_111", "other")))
pm[, t2group := grp(target)]
t2src <- pm[t2group != "other" & (t2group == "identity" | cmax_in_limits)]
type2 <- rbind(t2src[, .(scope = "both", n_cells = .N, min_pct = min(type2_pct), max_pct = max(type2_pct), min_cell = paste(pk_model, code)[which.min(type2_pct)],
                         max_cell = paste(pk_model, code)[which.max(type2_pct)]), by = .(config, t2group)],
               t2src[, .(scope = pk_model[1], n_cells = .N, min_pct = min(type2_pct), max_pct = max(type2_pct), min_cell = code[which.min(type2_pct)],
                         max_cell = code[which.max(type2_pct)]), by = .(config, t2group, pk_model)][, pk_model := NULL])
excl <- pm[t2group %in% c("near_095_105", "mid_090_111") & !cmax_in_limits & config == "P2", .(pk_model, code, cmax_ratio)]
fwrite(type2, file.path(out_dir, "oc_type2_summary.csv"))
say(sprintf("type II cells excluded for Cmax outside the limits: %s", if (nrow(excl)) paste(excl$pk_model, excl$code, collapse = ", ") else "none"))

# ---- 쌍대 차이와 분해 ---------------------------------------------------------------------------------------------------
pd <- function(a, b) paired_prop_diff_ci(a, b)                       # 100 x mean(a - b), 정규 근사 95% 구간
paired <- rbindlist(lapply(MODELS, function(m) {
  w <- PT[[m]]$w[model == AM]; cl <- CL[[m]]$cells
  rbindlist(lapply(c("F3A_iii", "G2A_iii", "F3B", "G2B", "G2_true"), function(X) {
    w[, { d <- pd(get(X), P2); .(comparison = paste(X, "- P2"), diff_pass_pp = d$est, lo = d$lo, hi = d$hi, n = d$n) }, by = scenario]
  }))[, pk_model := m][cl, on = c(scenario = "code"), nomatch = NULL]
}))
fwrite(paired[, .(pk_model, scenario, mechanism, direction, target, auc_ratio, cmax_ratio, cmax_in_limits, kind, comparison, diff_pass_pp, lo, hi, n)], file.path(out_dir, "oc_paired.csv"))

decomp <- rbindlist(lapply(MODELS, function(m) {
  w <- PT[[m]]$w[model == AM]; cl <- CL[[m]]$cells; nn <- PT[[m]]$n[model == AM]
  cells_t2 <- cl[grp(target) != "other" & (code == "S00" | cmax_in_limits), code]
  rbindlist(lapply(c(G2 = "G2", F3 = "F3"), function(X) {
    A <- if (X == "G2") "G2A_iii" else "F3A_iii"; B <- if (X == "G2") "G2B" else "F3B"
    w[scenario %in% cells_t2, {
      tot <- pd(P2, get(A)); est <- pd(P2, get(B)); red <- pd(get(B), get(A))
      .(config = X, total_pp = tot$est, total_lo = tot$lo, total_hi = tot$hi, estimation_pp = est$est, estimation_lo = est$lo, estimation_hi = est$hi,
        reduction_pp = red$est, reduction_lo = red$lo, reduction_hi = red$hi, n_trials = tot$n) }, by = scenario]
  }))[, pk_model := m][
    nn[scenario %in% cells_t2, .(n_arm_Aiii = median((n_R_AUCinf_Aiii + n_T_AUCinf_Aiii) / 2), n_arm_B = median((n_R_AUCinf_B + n_T_AUCinf_B) / 2),
                                 n_arm_AUClast = median((n_R_AUClast + n_T_AUClast) / 2)), by = scenario], on = "scenario"][cl, on = c(scenario = "code"), nomatch = NULL]
}))
decomp[, additive_gap := total_pp - estimation_pp - reduction_pp]
if (any(abs(decomp$additive_gap) > 1e-9)) stop("분해의 두 부분이 합계와 맞지 않습니다")
fwrite(decomp[, .(pk_model, scenario, mechanism, direction, target, config, total_pp, total_lo, total_hi, estimation_pp, estimation_lo, estimation_hi,
                  reduction_pp, reduction_lo, reduction_hi, n_trials, n_arm_AUClast, n_arm_B, n_arm_Aiii)], file.path(out_dir, "oc_decomposition.csv"))

# ---- 간이 계산 대조(동일 제품) ------------------------------------------------------------------------------------------
an_pass <- function(cv, n) { s <- sqrt(log(1 + (cv / 100)^2)); se <- s * sqrt(2 / n); t <- qt(0.95, 2 * n - 2)
  pmax(0, pnorm((log(lims[2]) - t * se) / se) - pnorm((log(lims[1]) + t * se) / se)) }
ana_a <- rbindlist(lapply(MODELS, function(m) {
  se <- PT[[m]]$se[scenario == "S00"]; se[, s_w := se / sqrt(1 / n_R + 1 / n_T)]
  st <- se[, .(s_median = median(s_w), n_arm = median((n_R + n_T) / 2)), by = endpoint][, cv_pct := 100 * sqrt(exp(s_median^2) - 1)]
  a <- rbind(data.table(case = "directive", endpoint = c("AUClast", "AUCinf_Aiii"), cv_pct = 40, n_arm = c(117, 75)), st[, .(case = "study", endpoint, cv_pct, n_arm)])
  a[, type2_analytic_pct := 100 * (1 - an_pass(cv_pct, n_arm))][, pk_model := m]
  sim <- pass_tab[pk_model == m & code == "S00" & config %in% c("AUClast_only", "AUCinf_Aiii_only", "P2", "G2A_iii")]
  sim_w <- dcast(sim, config ~ analysis_model, value.var = "pass_pct")
  a[, sim_single_M1_pct := 100 - sim_w[config == fifelse(endpoint == "AUClast", "AUClast_only", "AUCinf_Aiii_only")[1], M1], by = endpoint]
  a[, sim_single_M0_pct := 100 - sim_w[config == fifelse(endpoint == "AUClast", "AUClast_only", "AUCinf_Aiii_only")[1], M0], by = endpoint]
  a[, sim_config_M1_pct := 100 - sim_w[config == fifelse(endpoint == "AUClast", "P2", "G2A_iii")[1], M1], by = endpoint]
  a[, sim_config := fifelse(endpoint == "AUClast", "P2", "G2A_iii")]
  a
}))
fwrite(ana_a, file.path(out_dir, "oc_analytic_check.csv"))

# ---- 편향과 1종 오류의 관계(경계 16칸, 단일 지표; 별첨 A5d 상자. 사전 등록 뒤 더한 기술 요약이며 문구 규칙에 쓰지 않는다, D-065) ------------------
# 가까운 한계 쪽 1종 오류 ≈ Φ(b/SD - t): b = 로그 기하평균비 평균이 참 AUCinf 비보다 1 쪽으로 치우친 양(로그), SD = 시험 간 표준편차, t = qt(0.95, n - 3)(M1).
# 예측값은 가까운 한계만 쓴 정규 근사: Φ((평균 - log 0.80 - t x SE 중앙값)/SD)(참 비 < 1), Φ((log 1.25 - 평균 - t x SE 중앙값)/SD)(참 비 > 1).
relation <- rbindlist(lapply(MODELS, function(m) {
  e <- PT[[m]]$est; cl <- CL[[m]]$cells[kind == "boundary"]
  s <- e[, .(n_trials = .N, mean_est = mean(est, na.rm = TRUE), sd_est = sd(est, na.rm = TRUE), se_median = median(se, na.rm = TRUE), df_median = median(n_R + n_T - 3, na.rm = TRUE),
             sim_pass_pct = 100 * mean(pass %in% TRUE)), by = .(scenario, endpoint)]
  s <- cl[, .(pk_model, scenario = code, mechanism, direction, target, auc_ratio)][s, on = "scenario"]
  s[, sgn := fifelse(auc_ratio < 1, 1, -1)][, bias_log := sgn * (mean_est - log(auc_ratio))][, bias_pct := 100 * (exp(bias_log) - 1)][, b_over_sd := bias_log / sd_est]
  s[, t_q := qt(1 - alpha_pct / 100, df_median)]
  s[, pred_pass_pct := 100 * fifelse(auc_ratio < 1, pnorm((mean_est - log(lims[1]) - t_q * se_median) / sd_est), pnorm((log(lims[2]) - mean_est - t_q * se_median) / sd_est))]
  s[, sgn := NULL][]
}))
fwrite(relation, file.path(out_dir, "oc_bias_relation.csv"))
say(sprintf("bias relation: %d boundary cell x endpoint rows; |predicted - simulated| single-endpoint pass at most %.2f percentage points",
            nrow(relation), relation[, max(abs(pred_pass_pct - sim_pass_pct))]))

# ---- 문구 규칙(section8 wording) -------------------------------------------------------------------------------------
per_m <- function(cf, m, what) {
  if (what == "t1max") return(pm[config == cf & pk_model == m & kind == "boundary", max(pass_pct)])
  if (what == "t1n") return(pm[config == cf & pk_model == m & kind == "boundary", sum(pass_pct > alpha_pct)])
  if (what == "t2id") return(pm[config == cf & pk_model == m & kind == "identity", type2_pct])
}
t2cells_m <- function(m) pm[pk_model == m & config == "P2" & t2group != "other" & (t2group == "identity" | cmax_in_limits), code]
f3_t2inc <- rbindlist(lapply(MODELS, function(m) paired[pk_model == m & comparison == "F3A_iii - P2" & scenario %in% t2cells_m(m), .(pk_model = m, scenario, inc_pp = -diff_pass_pp)]))
f3_t1red <- paired[comparison == "F3A_iii - P2" & kind == "boundary", .(pk_model, scenario, red_pp = -diff_pass_pp)]
rules <- rbindlist(lapply(MODELS, function(m) data.table(pk_model = m,
  P2_t1_max = per_m("P2", m, "t1max"), P2_t2_id = per_m("P2", m, "t2id"), F3_t2_id = per_m("F3A_iii", m, "t2id"), G2_t2_id = per_m("G2A_iii", m, "t2id"),
  G2_t1_max = per_m("G2A_iii", m, "t1max"), G2B_t1_max = per_m("G2B", m, "t1max"), P2_t1_n = per_m("P2", m, "t1n"), G2_t1_n = per_m("G2A_iii", m, "t1n"), G2B_t1_n = per_m("G2B", m, "t1n"),
  F3_t2_inc_max = f3_t2inc[pk_model == m, max(inc_pp)], F3_t2_inc_max_cell = f3_t2inc[pk_model == m][which.max(inc_pp), scenario],
  F3_t1_red_max = f3_t1red[pk_model == m, max(red_pp)], F3_t1_red_max_cell = f3_t1red[pk_model == m][which.max(red_pp), scenario])))
rules[, `:=`(rule_T = P2_t1_max <= 6.0 & P2_t2_id <= F3_t2_id & P2_t2_id <= G2_t2_id,
             rule_4a_keep_type2 = (G2_t2_id - P2_t2_id) >= 2.0,
             rule_4b_small_cost = F3_t2_inc_max < 2.0,
             rule_protect_almost_none = F3_t1_red_max < 1.0,
             rule_G2_type1 = (G2_t1_max - P2_t1_max) >= 2.0 & G2_t1_n > P2_t1_n,
             rule_G2B_type1 = (G2B_t1_max - P2_t1_max) >= 2.0 & G2B_t1_n > P2_t1_n)]
c4 <- paired[comparison == "G2A_iii - P2" & on_t(target, c(0.90, 0.95)) & lo > 0, .(pk_model, scenario, mechanism, target, diff_pass_pp, lo, hi)]
fwrite(rules, file.path(out_dir, "oc_wording_rules.csv"))
fwrite(c4, file.path(out_dir, "oc_rule4c_cells.csv"))
all_rule <- function(r) all(rules[[r]])
say(sprintf("rules (both models must hold): T %s, 4a keep type II claim %s, 4b small F3 cost %s, protect almost none %s, G2 type I %s, G2B type I %s; 4c cells %d",
            all_rule("rule_T"), all_rule("rule_4a_keep_type2"), all_rule("rule_4b_small_cost"), all_rule("rule_protect_almost_none"), all_rule("rule_G2_type1"), all_rule("rule_G2B_type1"), nrow(c4)))
for (m in MODELS) say(sprintf("%s: type I max P2 %.2f, F3 %.2f, G2 %.2f, G2B %.2f; identity type II P2 %.2f, F3 %.2f, G2 %.2f; F3 type II increase max %.2f (%s); F3 type I reduction max %.2f (%s)",
                              m, rules[pk_model == m, P2_t1_max], pm[config == "F3A_iii" & pk_model == m & kind == "boundary", max(pass_pct)], rules[pk_model == m, G2_t1_max], rules[pk_model == m, G2B_t1_max],
                              rules[pk_model == m, P2_t2_id], rules[pk_model == m, F3_t2_id], rules[pk_model == m, G2_t2_id], rules[pk_model == m, F3_t2_inc_max], rules[pk_model == m, F3_t2_inc_max_cell],
                              rules[pk_model == m, F3_t1_red_max], rules[pk_model == m, F3_t1_red_max_cell]))
fwrite(rbindlist(lapply(MODELS, function(m) CL[[m]]$unreachable)), file.path(out_dir, "oc_unreachable_cells.csv"))
say("done")
