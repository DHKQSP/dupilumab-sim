#!/usr/bin/env Rscript
# §4-2: 500회 결과와 5,000회 결과 비교. 5,000회의 1–500번 시험은 기존 500회와 같은 난수 구조이므로 먼저 동일성을 확인하고,
# 500회 추정치를 독립인 501–5,000번 시험과 비교한다(z = 차이 / sqrt(p(1−p)(1/500 + 1/4500))).
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "trials5000")
b5 <- read_raw(proj_path("results", "trials5000", "products5000_be_raw_base.csv"))
b500 <- read_raw(proj_path("results", "trials", "products_be_raw.csv"))[method == "pooled_t"]
common <- intersect(unique(b5$scenario), unique(b500$scenario)); eps <- c("Cmax", "AUClast", "AUCinf_all", "AUCinf_reliable", "AUCinf_true")
idc <- merge(b500[scenario %in% common & endpoint %in% eps, .(trial, scenario, endpoint, GMR_500 = GMR)], b5[trial <= 500 & endpoint %in% eps, .(trial, scenario, endpoint, GMR_5000 = GMR)], by = c("trial", "scenario", "endpoint"))
ident <- data.table(rows_compared = nrow(idc), max_abs_diff = max(abs(idc$GMR_500 - idc$GMR_5000)))
cmp <- rbindlist(lapply(common, function(s_) rbindlist(lapply(eps, function(ep) {
  a <- b500[scenario == s_ & endpoint == ep, pass]; r <- b5[scenario == s_ & endpoint == ep & trial > 500, pass]; all5 <- b5[scenario == s_ & endpoint == ep, pass]
  p <- mean(all5); se <- sqrt(max(p * (1 - p), 1e-12) * (1 / length(a) + 1 / length(r)))
  data.table(scenario = s_, endpoint = ep, pass_500 = 100 * mean(a), pass_rest = 100 * mean(r), pass_all = 100 * p, n_all = length(all5), z = (mean(a) - mean(r)) / se)
}))))
cmp[, outside_mc := abs(z) > 1.96]
fwrite(ident, file.path(out_dir, "mc_consistency_identity.csv")); fwrite(cmp, file.path(out_dir, "mc_consistency_500_vs_5000.csv"))
print(ident); print(cmp[order(-abs(z))][1:15]); cat("MC 오차 밖(|z| > 1.96):", sum(cmp$outside_mc), "/", nrow(cmp), "(우연 기대치 약", round(0.05 * nrow(cmp), 1), ")\n")
