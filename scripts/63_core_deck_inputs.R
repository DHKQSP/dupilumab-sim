#!/usr/bin/env Rscript
# 63_core_deck_inputs.R: 핵심 덱(v1.1, 지시 2026-09-29)의 요약 파일(results/core_deck/). 사전 등록 config/prereg_20260929.yaml section7 (D-063).
# 새 시험 모의 없음. 저장된 대상자 수준 결과(.rds, git 밖)를 요약하고, 개체 수준 파일이 없는 곡선 모양 변형(scripts/21)과
# 대표 대상자의 관측값(scripts/10)은 원래 시드 규칙으로 다시 만든다. 다시 만든 값은 커밋된 요약·저장된 NCA와 같아야 하며(허용 오차 1e-9) 다르면 멈춘다.
# 산출(results/core_deck/):
#   coverage_by_case.csv      7a 창 포착률(참 AUC0-tlast / 참 AUC0-inf) 케이스별 n, 중앙값, 5·95백분위, 최솟값, 80% 미만 비율
#   title_rule.csv            7b S8 제목 규칙 판정(모든 케이스 최솟값 >= 84%인지)과 제목에 쓰는 값
#   rep_subjects.csv          7c 대표 대상자(두 모델: 중앙값·5백분위·최솟값) NCA와 참값, 외삽 면적(NCA 대 참)
#   rep_profiles.csv          7c 대표 대상자 참 농도 곡선(0.05일 간격, 연구일 1~70)
#   rep_obs.csv               7c 대표 대상자 B0 관측값(연구일, 농도, BLQ 여부, 참 농도)
#   reliable_iii_by_schedule.csv 7d 세트 (iii) 충족 비율, 일정별, B0 대비 차(%p, 같은 대상자)
#   provenance.csv            입력 SHA-256과 대조 결과
source("R/00_setup.R"); source_project()
pr <- read_cfg("prereg_20260929.yaml")$section7
design <- read_cfg("trial_design.yaml"); sc <- load_scenarios()
MASTER_SEED <- 20260923L                                                  # scripts/10, scripts/21과 같은 값
out_dir <- proj_path(pr$output_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("core_deck_inputs", master_seed = MASTER_SEED, run_mode = "final", extra = list(prereg = "config/prereg_20260929.yaml section7"))
sha <- function(f) digest::digest(file = f, algo = "sha256")
premise <- function(ok, what) if (!isTRUE(ok)) stop("premise failed: ", what, call. = FALSE)
PROV <- list()
chk <- function(what, a, b, tol = 1e-9) {
  ok <- length(a) == length(b) && isTRUE(all(abs(a - b) <= tol))
  show <- function(x) if (length(x) <= 10L) paste(signif(x, 12), collapse = ";") else sprintf("%d values, sha256 %s", length(x), digest::digest(signif(x, 12), algo = "sha256"))   # 긴 벡터는 개수와 해시로 기록
  PROV[[length(PROV) + 1L]] <<- data.table(check = what, value = show(a), committed = show(b), equal = ok)
  if (!ok) stop("mismatch with committed result: ", what)
}
src <- function(f) { PROV[[length(PROV) + 1L]] <<- data.table(check = paste("input", f), value = sha(proj_path(f)), committed = NA_character_, equal = NA); proj_path(f) }
cstats <- function(v) { v <- v[is.finite(v)]
  list(n = length(v), median = median(v), p05 = quantile(v, 0.05, names = FALSE), p95 = quantile(v, 0.95, names = FALSE), min = min(v), pct_lt80 = 100 * mean(v < 0.80)) }
# 커밋된 요약의 외삽(%) 통계와 대조: 외삽 = 100 x (1 - 창 포착률)이므로 외삽 95백분위 = 100 x (1 - 창 포착률 5백분위)(type 7은 선형 변환에 대칭)
chk_extrap <- function(what, v, ref, cols = c("n", "extrap_true_median", "extrap_true_p95", "coverage_lt80_pct"), max_col = NULL) {
  s <- cstats(v); a <- c(s$n, 100 * (1 - s$median), 100 * (1 - s$p05), s$pct_lt80); b <- unlist(ref[, ..cols])
  if (!is.null(max_col)) { a <- c(a, 100 * (1 - s$min)); b <- c(b, ref[[max_col]]) }
  chk(what, a, b, tol = 1e-8)
}

# ---- 7a 케이스별 창 포착률 -------------------------------------------------------------------------------------------------------------------
IND <- c(k2016 = "results/individual/nca_base_20000.rds", k2020 = "results/individual/nca_struct2020_20000.rds")
ind <- lapply(IND, function(f) readRDS(src(f)))
tcv <- fread(src("results/trialpop/tp_coverage_individual.csv"))[grepl("^window coverage", metric)]
rows <- list()
add_case <- function(case, label, model, v, source, checked) rows[[length(rows) + 1L]] <<- as.data.table(c(list(case = case, label = label, model = model), cstats(v), list(source = source, identity_checked = checked)))
for (m in names(IND)) {
  v <- ind[[m]][schedule == "B0", coverage_true]; s <- cstats(v); r <- tcv[pk_model == m]
  chk(sprintf("%s base B0: n, median, p05, p95, min vs tp_coverage_individual", m), c(s$n, s$median, s$p05, s$p95, s$min), c(r$n, r$median, r$p05, r$p95, r$min))
  add_case(paste0(m, "_base"), sprintf("%s model, base", substr(m, 2, 5)), m, v, paste0(IND[[m]], " [schedule == 'B0']"), TRUE)
}
# 곡선 모양 변형(scripts/21과 같은 호출; 태그 curve_<변형>). base도 함께 다시 만들어 대조(행이 아니라 노트용)
CS <- fread(src("results/curve_shape/curve_shape_B0.csv"))
curve <- list()
for (v in c("base", "vmax080_both", "vmax125_both", "km05_both", "km10_both")) {
  rv <- resolve_variant(v, design, sc); p <- rv$p
  pop <- run_individual_population(20000L, p, design, "B0", MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = paste0("curve_", v))
  x <- pop$nca; curve[[v]] <- x
  chk_extrap(sprintf("curve-shape %s regenerated vs curve_shape_B0.csv (n, extrap median, p95, lt80, max)", v), x$coverage_true, CS[variant == v], max_col = "extrap_true_max")
  cat(sprintf("  curve %s regenerated %s\n", v, format(Sys.time(), "%H:%M:%S")))
}
CL <- c(vmax080_both = "vmax080", vmax125_both = "vmax125", km05_both = "km05", km10_both = "km10")
CLAB <- c(vmax080_both = "Vmax x0.8 (both arms, 2016 model)", vmax125_both = "Vmax x1.25 (both arms, 2016 model)", km05_both = "Km x0.5 (both arms, 2016 model)", km10_both = "Km x10 (both arms, 2016 model)")
for (v in names(CL)) add_case(CL[[v]], CLAB[[v]], "k2016", curve[[v]]$coverage_true, sprintf("regenerated scripts/21 variant %s (tag curve_%s)", v, v), TRUE)
# LLOQ(2016 모델, 추정 잔차)
LI <- readRDS(src("results/lloq/lloq_individual_nca_k2016.rds")); LT <- fread(src("results/lloq/lloq_individual_table.csv"))
for (lq in c(0.02, 0.5)) {
  v <- LI[abs(lloq - lq) < 1e-9 & resid == "fixed", coverage_true]
  chk_extrap(sprintf("LLOQ %s (k2016, fixed) vs lloq_individual_table.csv", lq), v, LT[model == "k2016" & resid == "fixed" & abs(lloq - lq) < 1e-9])
  add_case(if (lq < 0.1) "lloq002" else "lloq05", sprintf("LLOQ %s mg/L (2016 model, estimated residual)", format(lq)), "k2016", v, sprintf("results/lloq/lloq_individual_nca_k2016.rds [lloq == %s & resid == 'fixed']", format(lq)), TRUE)
}
# 비례 잔차 12%
R12 <- readRDS(src("results/individual/nca_resid12_20000.rds"))[schedule == "B0"]; IR <- fread(src("results/individual/individual_resid12.csv"))
chk_extrap("resid12 B0 vs individual_resid12.csv", R12$coverage_true, IR[schedule == "B0"])
add_case("resid12", "proportional residual 12% (2016 model)", "k2016", R12$coverage_true, "results/individual/nca_resid12_20000.rds [schedule == 'B0']", TRUE)
# 체중 층(2016 모델 기본, B0): 분할점은 trial_design.yaml
split_kg <- as.numeric(design$stratification$split_kg$value); b16 <- ind$k2016[schedule == "B0"]
premise(nrow(b16[WT <= split_kg]) + nrow(b16[WT > split_kg]) == nrow(b16), "weight strata cover the base population")
add_case("wt60_75", sprintf("weight stratum 60-%g kg (2016 model, base)", split_kg), "k2016", b16[WT <= split_kg, coverage_true], sprintf("%s [B0 & WT <= %g]", IND[["k2016"]], split_kg), FALSE)
add_case("wt75_90", sprintf("weight stratum >%g-90 kg (2016 model, base)", split_kg), "k2016", b16[WT > split_kg, coverage_true], sprintf("%s [B0 & WT > %g]", IND[["k2016"]], split_kg), FALSE)
CV <- rbindlist(rows)[, role := "row"]
cb <- as.data.table(c(list(case = "curve_base", label = "curve-shape study base draw (2016 model; note only)", model = "k2016"), cstats(curve$base$coverage_true),
                      list(source = "regenerated scripts/21 variant base (tag curve_base)", identity_checked = TRUE, role = "note")))
CV <- rbind(CV, cb)
fwrite(CV, file.path(out_dir, "coverage_by_case.csv"))

# ---- 7b 제목 규칙 ------------------------------------------------------------------------------------------------------------------------
rw <- CV[role == "row"]
TR <- data.table(all_cases_min_ge_84 = all(rw$min >= 0.84), base_min_k2016 = rw[case == "k2016_base", min], base_min_k2020 = rw[case == "k2020_base", min],
                 all_min = min(rw$min), all_min_case = rw[which.min(min), case], max_pct_lt80 = max(rw$pct_lt80), max_pct_lt80_case = rw[which.max(pct_lt80), case],
                 min_p05 = min(rw$p05), min_p05_case = rw[which.min(p05), case], n_cases = nrow(rw))
fwrite(TR, file.path(out_dir, "title_rule.csv"))

# ---- 7c 대표 대상자 -----------------------------------------------------------------------------------------------------------------------
VAR <- c(k2016 = "base", k2020 = "struct2020"); scheds <- design$schedule_analysis; n_ind <- as.integer(design$mc$n_individual)
LLOQ <- study_lloq(); CHK_COLS <- c("tlast", "Clast", "lambda_z", "AUClast", "AUClast_true", "AUCinf_true", "coverage_true")
RS <- list(); RP <- list(); RO <- list()
for (m in names(VAR)) {
  v <- VAR[[m]]; rv <- resolve_variant(v, design, sc); p <- rv$p
  pop <- run_individual_population(n_ind, p, design, scheds, MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = paste0("ind_", v))
  a <- pop$nca[schedule == "B0"][order(id)]; b <- ind[[m]][schedule == "B0"][order(id)]
  for (cc in CHK_COLS) { x1 <- a[[cc]]; x2 <- b[[cc]]; chk(sprintf("%s regenerated B0 NCA equals stored: %s (NA pattern and values)", m, cc), c(sum(is.na(x1)), x1[!is.na(x1)]), c(sum(is.na(x2)), x2[!is.na(x2)])) }
  cv <- b$coverage_true; q <- c(median = median(cv), p05 = quantile(cv, 0.05, names = FALSE))
  pick <- function(target) { d <- abs(cv - target); b$id[which(d == min(d))[1]] }                   # 같은 거리면 작은 id(정렬됨)
  ids <- c(median = pick(q[["median"]]), p05 = pick(q[["p05"]]), min = b$id[which(cv == min(cv))[1]])
  sel <- b[id %in% ids]; sel[, reliable_iii := crit_ok(sel, "iii")]
  premise(nrow(sel) == 3, sprintf("%s: three representative subjects", m))                       # λz 산출 불가 대상자도 그대로 둔다(외삽 면적은 NA)
  RS[[m]] <- data.table(model = m, role = names(ids), id = unname(ids))[sel, on = "id"][, .(model, role, id, WT, coverage_true, pct_extrap_true, tlast, study_day_tlast = tlast + 1,
    Clast, lambda_z, lambda_z_lower = Lambda_z_lower, lambda_z_upper = Lambda_z_upper, adj_r2 = Rsq_adjusted, pct_extrap_nca = `AUC_%Extrap_obs`, lambda_ok = lambda_ok %in% TRUE, AUClast, AUClast_true, AUCinf_true,
    AUCinf_nca = AUCINF_obs, extrap_area_nca = AUCINF_obs - AUClast, extrap_area_true = AUCinf_true - AUClast_true, reliable_iii)]
  RS[[m]][, extrap_ratio_nca_to_true := extrap_area_nca / extrap_area_true]
  ob <- subset_schedule(pop$obs, get_schedule(design, "B0"))[id %in% ids & planned > 0]
  RO[[m]] <- ob[, .(model = m, id, planned_day = planned, time_after_dose = time, study_day = time + 1, conc_obs = conc, blq, conc_true = C, auc_true = auc)]
  # 관측 tlast 시점의 참 누적 AUC가 저장된 AUClast_true와 같은지
  for (i in ids) chk(sprintf("%s subject %d: true AUC at observed tlast", m, i), RO[[m]][id == i & time_after_dose == RS[[m]][id == i, tlast], auc_true], RS[[m]][id == i, AUClast_true])
  ipr <- individual_params(p, pop$subj[id %in% ids])
  # 같은 개인 파라미터로 관측 시각을 다시 풀면 모의 관측의 참 농도와 같아야 한다(적분 격자가 달라 상대 오차 1e-6 허용)
  re <- as.data.table(solve_model(ipr, ob[, .(id, time)], design$dose_mg, p$model_id))[, .(id, time, C2 = C)][ob, on = c("id", "time")]
  chk(sprintf("%s representative subjects: re-solved true concentration at the observed times (max relative difference)", m), max(abs(re$C2 / re$C - 1)), 0, tol = 1e-6)
  pf <- as.data.table(solve_model(ipr, CJ(id = unname(ids), time = seq(0.05, 69, by = 0.05)), design$dose_mg, p$model_id))
  RP[[m]] <- pf[, .(model = m, id, time_after_dose = time, study_day = time + 1, conc = C, auc = auc)]
  cat(sprintf("  representative subjects %s regenerated %s\n", m, format(Sys.time(), "%H:%M:%S")))
}
fwrite(rbindlist(RS), file.path(out_dir, "rep_subjects.csv")); fwrite(rbindlist(RO), file.path(out_dir, "rep_obs.csv")); fwrite(rbindlist(RP), file.path(out_dir, "rep_profiles.csv"))

# ---- 7d 채혈 일정별 세트 (iii) 충족 비율 --------------------------------------------------------------------------------------------------
TF <- fread(src("results/trialpop/tp_failure_by_set.csv"))
RL <- rbindlist(lapply(names(IND), function(m) {
  x <- ind[[m]]; b0 <- x[schedule == "B0"][order(id)]
  chk(sprintf("%s set (iii) B0 failing share vs tp_failure_by_set.csv", m), 100 * mean(!crit_ok(b0, "iii")), TF[pk_model == m & set == "iii", fail_pct])
  rbindlist(lapply(c("B0", "D1", "D2", "D3", "D4"), function(sh) { y <- x[schedule == sh][order(id)]; premise(identical(y$id, b0$id), "same subjects in every schedule")
    data.table(model = m, schedule = sh, n = nrow(y), pct_reliable_iii = 100 * mean(crit_ok(y, "iii")), diff_vs_B0_pp = 100 * (mean(crit_ok(y, "iii")) - mean(crit_ok(b0, "iii")))) }))
}))
fwrite(RL, file.path(out_dir, "reliable_iii_by_schedule.csv"))

fwrite(rbindlist(PROV), file.path(out_dir, "provenance.csv"))
print(CV[, .(case, model, n, median = round(100 * median, 2), p05 = round(100 * p05, 2), min = round(100 * min, 2), pct_lt80, role)]); print(TR); print(rbindlist(RS)[, .(model, role, id, coverage_true, extrap_ratio_nca_to_true)]); print(RL)
append_run_log(logfile, "done")
