#!/usr/bin/env Rscript
# 50_sap_support.R — 통계분석계획(SAP) 제안 문안(§4)의 보조 수치. 새 모의 없음: 저장된 20,000명 B0 NCA(results/individual/nca_<v>_20000.rds)에서
# AUC0-inf 처리 규칙별 제외 대상자의 비율·체중·노출을 계산한다(제외가 무작위가 아님을 보이는 근거; D-056).
# 산출: results/sap/sap_exclusion_support.csv
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "sap"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
geo <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
res <- rbindlist(lapply(c(k2016 = "base", k2020 = "struct2020"), function(v) {
  x <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", v)))[schedule == "B0"]
  x[, `:=`(ok_B = lambda_ok %in% TRUE, ok_i = (lambda_ok & !(flag_rsq %in% TRUE) & !(flag_extrap %in% TRUE)) %in% TRUE, ok_ii = reliable %in% TRUE)]
  rbindlist(lapply(c(rule_B = "ok_B", criteria_i = "ok_i", criteria_ii = "ok_ii"), function(cn) {
    ok <- x[[cn]]
    t_ <- t.test(x$WT[!ok], x$WT[ok])
    data.table(n = nrow(x), n_excluded = sum(!ok), excluded_pct = 100 * mean(!ok), wt_retained_mean = mean(x$WT[ok]), wt_excluded_mean = mean(x$WT[!ok]),
               wt_diff = mean(x$WT[!ok]) - mean(x$WT[ok]), wt_diff_lo = t_$conf.int[1], wt_diff_hi = t_$conf.int[2],
               auclast_gm_ratio_excluded_to_retained = geo(x$AUClast[!ok]) / geo(x$AUClast[ok]), aucinf_true_gm_ratio_excluded_to_retained = geo(x$AUCinf_true[!ok]) / geo(x$AUCinf_true[ok]))
  }), idcol = "set")
}), idcol = "model")
res[, source := "results/individual/nca_<variant>_20000.rds, schedule B0"]
fwrite(res, file.path(out_dir, "sap_exclusion_support.csv"))
print(res[, .(model, set, excluded_pct, wt_retained_mean, wt_excluded_mean, wt_diff, auclast_gm_ratio_excluded_to_retained)])
