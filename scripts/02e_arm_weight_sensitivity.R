#!/usr/bin/env Rscript
# 완전 외부 gate 보조: Li 2020 arm(체중 미보고)의 모의/관측 비가 가정 체중에 얼마나 의존하는지. 판정 변경에 쓰지 않는다(서술용).
# 사용법: Rscript scripts/02e_arm_weight_sensitivity.R [k2016|k2020]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); MODEL <- if (length(args)) args[1] else "k2016"
p <- load_params(MODEL); dz <- read_cfg("design_clot2021.yaml"); sr <- dz$cohort_sim$sex_ratio_male$value
out_dir <- proj_path("results", if (MODEL == "k2016") "step1" else paste0("step1_", MODEL))
arms <- rbindlist(lapply(dz$arm_checks_300mg$arms, function(a) data.table(study = a$study, arm = a$arm, schedule = a$schedule, obs = a$auclast_mean)))
res <- rbindlist(lapply(c(70, 74, 78, 82, 86, 90), function(W) rbindlist(lapply(unique(arms$schedule), function(sch) {
  sim <- simulate_dataset(p, 300, W, 10, schedule_days(dz, sch), 10000, sprintf("armwt_%s_%g", sch, W), sr, c(50, 110), c(40, 100))
  m <- mean(sim$nca$AUClast, na.rm = TRUE)
  arms[schedule == sch, .(model = MODEL, weight_mean = W, study, arm, obs, sim = m, ratio = m / obs, within_15 = abs(m / obs - 1) <= 0.15)]
}))))
fwrite(res, file.path(out_dir, "step1i_arm_weight_sensitivity.csv"))
print(dcast(res[, .(study, arm, weight_mean, r = round(ratio, 3))], study + arm ~ weight_mean, value.var = "r"))
