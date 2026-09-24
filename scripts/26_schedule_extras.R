#!/usr/bin/env Rscript
# §2 보고서 보강: (1) B− 명시 표, (2) AUClast CI 폭 변화의 원인 분해(잔차 없는 참 농도로 같은 일정 재계산)
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "trials")
dec <- fread(proj_path("results", "trials", "schedule_decision_base.csv"))
bm <- dec[schedule == "Bminus", .(schedule, a_width_B0_pp, a_width_pp, a_mean_width_rel_decrease, a_paired_lo, a_paired_hi, a_label,
                                  c_reliable_B0, c_reliable, c_reliable_gain_pp, c_gain_boot_lo, c_gain_boot_hi, d_extrap20_B0, d_extrap20, d_extrap20_ratio, d_ratio_boot_lo, d_ratio_boot_hi)]
bm[, recommendation := "Day 50 삭제 권고하지 않음 (AUC0-inf 이차 평가변수와 fallback용 말단 점 확보)"]
fwrite(bm, file.path(out_dir, "bminus_explicit.csv"))
pt_base <- fread(proj_path("results", "trials", "paired_trial_vs_B0_base.csv"))[scenario == "S00"]
pt_nr <- fread(proj_path("results", "trials", "paired_trial_vs_B0_noresid.csv"))[scenario == "S00"]
dcmp <- merge(pt_base[, .(schedule, base = width_rel_decrease, base_lo = width_rel_decrease_lo, base_hi = width_rel_decrease_hi)],
              pt_nr[, .(schedule, noresid = width_rel_decrease, noresid_lo = width_rel_decrease_lo, noresid_hi = width_rel_decrease_hi)], by = "schedule")
dcmp[, base_widened := base_hi < 0]
dcmp[, noresid_widened := noresid_hi < 0]
dcmp[, cause := fifelse(!base_widened, "해당 없음(확대 아님)", fifelse(noresid_widened, "tlast 연장에 따른 꼬리 면적 이질성", "추가된 저농도 점의 측정오차"))]
fwrite(dcmp, file.path(out_dir, "width_decomposition.csv")); print(bm); print(dcmp)
