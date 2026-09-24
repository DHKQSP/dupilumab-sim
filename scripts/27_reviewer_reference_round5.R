#!/usr/bin/env Rscript
# 통합본 §5-1·§6-1의 독립 구현 기준값(8,000명)과 이 구현(20,000명) 비교 — 서술용, 판정 기준 아님.
# 신뢰 충족률은 기준값과 같은 정의(λz 산출·adj R² ≥ 0.80·외삽 ≤ 20%, span 플래그 제외, D-039)로 비교한다.
# 입력: results/curve_shape/curve_shape_B0.csv, results/weight_generalization/weight_bands_B0_abcd.csv, config/crossval_reference.yaml
source("R/00_setup.R"); source_project()
ref <- read_cfg("crossval_reference.yaml"); out_dir <- proj_path("results", "crossval")
cs <- fread(proj_path("results", "curve_shape", "curve_shape_B0.csv"))
wb <- fread(proj_path("results", "weight_generalization", "weight_bands_B0_abcd.csv"))
row <- function(analysis, cond, metric, r, s) data.table(analysis = analysis, condition = cond, metric = metric, ref_8000 = r, sim_20000 = s)
out <- list()
kr <- ref$curve_shape_round5$km_range; k <- cs[variant %in% unlist(kr$variants)]
out[[length(out) + 1]] <- rbind(row("5-1", "Km ×0.5–×10", "extrap_true_median_min", kr$extrap_true_median_min, min(k$extrap_true_median)),
                                row("5-1", "Km ×0.5–×10", "extrap_true_median_max", kr$extrap_true_median_max, max(k$extrap_true_median)),
                                row("5-1", "Km ×0.5–×10", "extrap_true_p95_min", kr$extrap_true_p95_min, min(k$extrap_true_p95)),
                                row("5-1", "Km ×0.5–×10", "extrap_true_p95_max", kr$extrap_true_p95_max, max(k$extrap_true_p95)),
                                row("5-1", "Km ×0.5–×10", "coverage_lt80_pct_max", kr$coverage_lt80_pct, max(k$coverage_lt80_pct)))
for (v in c("vmax080_both", "vmax050_both")) { r <- ref$curve_shape_round5[[v]]; x <- cs[variant == v]
  out[[length(out) + 1]] <- rbindlist(lapply(names(r), function(m) row("5-1", v, m, r[[m]], x[[m]]))) }
for (mk in c("a", "c", "d")) for (b in names(ref$weight_bands_round5[[mk]])) { r <- ref$weight_bands_round5[[mk]][[b]]; x <- wb[model == mk & band == b]
  out[[length(out) + 1]] <- rbindlist(lapply(names(r), function(m) row("6-1", paste0("(", mk, ") ", b, " kg"), m, r[[m]],
    if (m == "reliable_pct" && "reliable_rsq_extrap_pct" %in% names(x)) x[["reliable_rsq_extrap_pct"]] else x[[m]]))) }   # 기준값은 adj R²·외삽만의 신뢰 정의
cmp <- rbindlist(out)
# 채혈 허용창 편차 없이 같은 대상자로 재계산한 결과(24_weight_generalization.R a nojitter ...)가 있으면 병기: 차이의 원인이 편차인지 확인
fnj <- proj_path("results", "weight_generalization", "weight_bands_B0_a_nojitter.csv")
if (file.exists(fnj)) { nj <- fread(fnj)
  cmp[, sim_20000_nojitter := NA_real_]
  for (i in which(cmp$analysis == "6-1" & startsWith(cmp$condition, "(a) "))) { b <- sub("^\\(a\\) (.*) kg$", "\\1", cmp$condition[i]); x <- nj[band == b]
    mc_ <- if (cmp$metric[i] == "reliable_pct" && "reliable_rsq_extrap_pct" %in% names(x)) "reliable_rsq_extrap_pct" else cmp$metric[i]
    if (nrow(x) && mc_ %in% names(x)) set(cmp, i, "sim_20000_nojitter", x[[mc_]]) } }
cmp[, abs_diff := sim_20000 - ref_8000]
cmp[, rel_diff_pct := fifelse(ref_8000 != 0, 100 * (sim_20000 / ref_8000 - 1), NA_real_)]
cmp[, note := fifelse(metric %in% c("coverage_lt80_pct", "coverage_lt80_pct_max", "extrap_true_max"), "희귀 사건·극값: 표본 수 차이로 변동 큼(참고)", "")]
fwrite(cmp, file.path(out_dir, "reviewer_reference_round5.csv"))
print(cmp, digits = 3)
cat(sprintf("|상대 차이| ≤ 10%%: %d / %d (희귀 사건·극값 제외 %d / %d)\n", sum(abs(cmp$rel_diff_pct) <= 10, na.rm = TRUE), sum(!is.na(cmp$rel_diff_pct)),
            cmp[note == "" & !is.na(rel_diff_pct), sum(abs(rel_diff_pct) <= 10)], cmp[note == "" & !is.na(rel_diff_pct), .N]))
