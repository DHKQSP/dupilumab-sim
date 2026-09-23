#!/usr/bin/env Rscript
# 단계 1 미달 원인 탐색(진단 전용, 파라미터 조정 아님 — SPEC §1.3).
# 어떤 가정(체중/Vc, Vmax, F, ke, LLOQ)이 얼마나 바뀌어야 각 용량의 대표 개체 tlast가 목표와 맞는지 스캔한다.
source("R/00_setup.R"); source_project()
p <- load_params("dev"); dz <- read_cfg("design_clot2021.yaml")
days <- as.numeric(dz$sample_days_post_dose); lloq <- dz$lloq_mg_L
targets <- rbindlist(lapply(dz$targets, function(x) data.table(dose = x$dose_mg, target_tlast = x$tlast_median_day)))
out_dir <- proj_path("results", "dev"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tlast_for <- function(mod_theta = list(), WT = NULL, lloq_val = lloq) {
  ip <- typical_subject(p, WT = WT)
  for (nm in names(mod_theta)) set(ip, j = nm, value = ip[[nm]] * mod_theta[[nm]])
  rbindlist(lapply(targets$dose, function(d) {
    tr <- true_auc(ip, d, t_grid = days); sched <- tr$profile[time > 0]
    tl <- true_auc_to_tlast(sched, lloq_val)
    data.table(dose = d, tlast = tl$tlast_true, C_at_target = sched[time == targets[dose == d, target_tlast], C],
               ratio_true_pct = 100 * tl$AUClast_true / tr$inf$AUCinf_true)
  }))
}
# 대표 개체의 LLOQ 교차 시각(연속) — 주간 채혈 tlast의 바로 위 해상도
crossing <- rbindlist(lapply(targets$dose, function(d) {
  ip <- typical_subject(p)
  tr <- true_auc(ip, d, t_grid = seq(0, 84, by = 1 / 24)); prof <- tr$profile[time > 0]
  above <- prof[C >= lloq]; cross <- max(above$time)
  data.table(dose = d, t_cross_lloq_day = cross, target_tlast = targets[dose == d, target_tlast],
             target_window = sprintf("[%d, %d)", targets[dose == d, target_tlast], targets[dose == d, target_tlast] + 7),
             in_window = cross >= targets[dose == d, target_tlast] & cross < targets[dose == d, target_tlast] + 7,
             shortfall_day = pmax(0, targets[dose == d, target_tlast] - cross))
}))
fwrite(crossing, file.path(out_dir, "step1_lloq_crossing.csv"))
cat("\n대표 개체 LLOQ 교차 시각 (연속시간):\n"); print(crossing)

scan <- rbindlist(list(
  rbindlist(lapply(c(1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6), function(m) tlast_for(list(Vmax = m))[, `:=`(factor = "Vmax x", level = m)])),
  rbindlist(lapply(c(1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4), function(m) tlast_for(list(F = m))[, `:=`(factor = "F x", level = m)])),
  rbindlist(lapply(c(1.0, 0.9, 0.8, 0.7), function(m) tlast_for(list(ke = m))[, `:=`(factor = "ke x", level = m)])),
  rbindlist(lapply(c(75, 70, 65, 60, 55, 50), function(w) tlast_for(WT = w)[, `:=`(factor = "WT kg (Vc, theta_WT=1 DEV)", level = w)])),
  rbindlist(lapply(c(0.078, 0.05, 0.02, 0.01), function(l) tlast_for(lloq_val = l)[, `:=`(factor = "LLOQ mg/L", level = l)]))
))
scan <- merge(scan, targets, by = "dose")[, match := tlast == target_tlast]
setcolorder(scan, c("factor", "level", "dose", "tlast", "target_tlast", "match", "C_at_target", "ratio_true_pct"))
setorder(scan, factor, -level, dose)
fwrite(scan, file.path(out_dir, "step1_sensitivity_scan.csv"))
w <- dcast(scan, factor + level ~ dose, value.var = "tlast")
setnames(w, c("200", "300", "600"), c("tlast_200", "tlast_300", "tlast_600"))
w[, all_match := tlast_200 == 28 & tlast_300 == 42 & tlast_600 == 56]
print(w)
