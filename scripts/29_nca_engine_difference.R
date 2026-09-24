#!/usr/bin/env Rscript
# 기존 엔진(run_nca_legacy) 대비 새 엔진(run_nca, Phoenix 호환) 차이표 (검토 의견 2026-09-24 §1-2, §1-3).
# 같은 가상 대상자·같은 관측치(B0, 60–90 kg, 채혈 허용창·잔차·BLQ)에 두 엔진을 적용한다. 모델별 20,000명.
# 산출: 개인 수준 지표 비교, 대상자별 변화 유형(λz 창·Tlast·신뢰 판정), 탈락 사유별 비율(중복 포함).
source("R/00_setup.R"); source_project()
design <- read_cfg("trial_design.yaml"); MASTER_SEED <- 20260923L; N <- 20000L
out_dir <- proj_path("results", "nca_engine"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("nca_engine_difference", master_seed = MASTER_SEED, run_mode = "final", extra = list(n = N))
res <- lapply(c(k2016 = "base", k2020 = "struct2020"), function(v) {
  rv <- resolve_variant(v, design)
  pop <- run_individual_population(N, rv$p, design, "B0", MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = rv$p$model_id, tag = paste0("ind_", v))
  ob <- subset_schedule(pop$obs, get_schedule(design, "B0"))
  new <- attach_truth(run_nca(ob), ob, pop$truth); old <- attach_truth(run_nca_legacy(ob), ob, pop$truth)
  old[, reliable_old_rule := reliable]
  list(new = new, old = old, v = v)
})
metric <- function(x, eng) data.table(engine = eng, n = nrow(x), AUClast_geo = geo_mean(x$AUClast), tlast_median = median(x$tlast, na.rm = TRUE),
  lambda_ok_pct = 100 * mean(x$lambda_ok), reliable_pct = 100 * mean(x$reliable),
  reliable_rsq_extrap_only_pct = 100 * mean(x$lambda_ok & x$adj_r2 >= 0.80 & x$pct_extrap <= 20, na.rm = TRUE),
  extrap_median = median(x$pct_extrap, na.rm = TRUE), extrap_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE),
  n_lambda_median = median(x$n_lambda, na.rm = TRUE), AUCinf_rel_geo = geo_mean(x$AUCinf[x$reliable]), AUCinf_all_geo = geo_mean(x$AUCinf[x$lambda_ok]),
  extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE))
cmp <- rbindlist(lapply(res, function(r) rbind(metric(r$old, "이전 엔진(D-010)"), metric(r$new, "새 엔진(Phoenix 호환, D-039)"))[, model := r$v]))
setcolorder(cmp, c("model", "engine"))
# 대상자별 변화 유형
chg <- rbindlist(lapply(res, function(r) {
  m <- merge(r$old[, .(id, tlast_o = tlast, AUClast_o = AUClast, lz_o = lambda_z, n_o = n_lambda, rel_o = reliable)],
             r$new[, .(id, tlast_n = tlast, AUClast_n = AUClast, lz_n = lambda_z, n_n = n_lambda, rel_n = reliable, flag_span, flag_rsq, flag_extrap, lambda_ok)], by = "id")
  data.table(model = r$v, n = nrow(m),
             tlast_changed_pct = 100 * mean(abs(m$tlast_o - m$tlast_n) > 1e-9, na.rm = TRUE),
             AUClast_changed_pct = 100 * mean(abs(m$AUClast_o / m$AUClast_n - 1) > 1e-9, na.rm = TRUE),
             lambda_window_changed_pct = 100 * mean(xor(is.na(m$lz_o), is.na(m$lz_n)) | (!is.na(m$lz_o) & !is.na(m$lz_n) & (m$n_o != m$n_n | abs(m$lz_o - m$lz_n) > 1e-9))),
             reliable_old_to_unreliable_new_pct = 100 * mean(m$rel_o & !m$rel_n), reliable_new_not_old_pct = 100 * mean(!m$rel_o & m$rel_n),
             lost_by_span_only_pct = 100 * mean(m$rel_o & !m$rel_n & m$flag_span & !m$flag_rsq & !m$flag_extrap))
}))
# 탈락 사유(새 엔진, 중복 포함)
reasons <- rbindlist(lapply(res, function(r) {
  x <- r$new[n_quant > 0]
  combo <- x[, .(pct = 100 * .N / nrow(x)), by = .(lambda_fail = !lambda_ok, rsq = flag_rsq, extrap = flag_extrap, span = flag_span)][order(-pct)]
  combo[, model := r$v][]
}))
reasons_marg <- rbindlist(lapply(res, function(r) dropout_reasons(r$new)[, model := r$v]))
fwrite(cmp, file.path(out_dir, "engine_difference_individual_B0.csv")); fwrite(chg, file.path(out_dir, "engine_difference_subject_changes.csv"))
fwrite(reasons, file.path(out_dir, "dropout_reason_combinations_B0.csv")); fwrite(reasons_marg, file.path(out_dir, "dropout_reasons_B0.csv"))
print(cmp, digits = 4); print(chg, digits = 3); print(reasons_marg, digits = 3); print(reasons, digits = 3)
append_run_log(logfile, "done")
