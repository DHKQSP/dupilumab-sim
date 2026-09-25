#!/usr/bin/env Rscript
# §2 절벽 채혈 무의미성 분석 (config/oc_design.yaml cliff, D-041).
# 절벽 = tmax 이후 순간 반감기 ln2/(−d ln C/dt)가 처음 1일(민감도 2일) 미만이 되는 시점부터 참 농도가 연구 LLOQ(config/assay.yaml, 현재 0.078 mg/L)에 닿는 시점까지.
# d ln C/dt는 ODE 상태값으로 정확히 계산: dC/dt = (ka·흡수구획 − (ke + k12)·central + k21·periph − Vmax·C/(Km + C)·Vc)/Vc. 교차 시점은 0.05일 격자 사이 선형 보간.
# 조건: 두 모델 × 체중(60–90 kg 기본, 90–110, 110–130, 130–150 kg 균등) × 20,000명 × 채혈 시각(명목, 방문 허용창).
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff; design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", if (length(commandArgs(trailingOnly = TRUE))) "cliff_test" else "cliff"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args)) as.integer(args[1]) else as.integer(cf$n_subjects); LLOQ <- study_lloq(); STEP <- as.numeric(cf$grid_step_day); SEED <- as.integer(cf$seed)
if (!isTRUE(all.equal(LLOQ, as.numeric(cf$lloq)))) message(sprintf("연구 LLOQ %g mg/L(config/assay.yaml)가 사전 고정 절벽 설계의 %g mg/L(oc_design.yaml, 779e068)와 다릅니다: 단일 출처 값을 씁니다(D-054)", LLOQ, as.numeric(cf$lloq)))
defs <- as.numeric(unlist(cf$definition_days))
logfile <- start_run_log("cliff_analysis", master_seed = SEED, run_mode = "final", extra = list(n = N, step = STEP))
b0 <- get_schedule(design, "B0")
SCHED <- list(current = b0, plus_39_46 = sort(c(b0, 38, 45)), plus_39_46_53 = sort(c(b0, 38, 45, 52)),
              plus_32_39_46_53 = sort(c(b0, 31, 38, 45, 52)), plus_40_47 = sort(c(b0, 39, 46)), daily_29_57 = sort(unique(c(b0, 28:56))))
SCHED_LABEL <- c(current = "현행", plus_39_46 = "+Day 39·46", plus_39_46_53 = "+Day 39·46·53", plus_32_39_46_53 = "+Day 32·39·46·53",
                 plus_40_47 = "+Day 40·47", daily_29_57 = "Day 29–57 매일")
min_interval <- vapply(SCHED, function(s) { x <- s[s >= 21]; min(diff(x)) }, numeric(1))    # 경과일 21(Day 22) 이후 최소 방문 간격

# 계산은 R/cliff.R(cliff_subjects, cliff_count_points; LLOQ 벡터 지원, D-054). 여기서는 연구 LLOQ 하나로 부른다.
cliff_subjects_one <- function(model, wname) {
  cs <- cliff_subjects(model, wname, N, SEED, STEP, LLOQ, defs, cf$weights[[wname]], design)
  cs[, lloq := NULL][]
}
count_points <- function(cs, timing, model, wname) cliff_count_points(cs, timing, model, wname, SCHED, SCHED_LABEL, min_interval, defs, SEED, design)

all_cs <- list(); pts <- list()
for (model in c("k2016", "k2020")) for (wname in names(cf$weights)) {
  t0 <- Sys.time(); cs <- cliff_subjects_one(model, wname)
  all_cs[[paste(model, wname)]] <- cs
  for (timing in unlist(cf$timing)) pts[[paste(model, wname, timing)]] <- count_points(cs, timing, model, wname)
  cat(sprintf("[%s] %s %s done in %s\n", format(Sys.time(), "%H:%M:%S"), model, wname, format(Sys.time() - t0)))
}
cs_all <- rbindlist(all_cs); pts <- rbindlist(pts)
saveRDS(cs_all, file.path(out_dir, "cliff_subjects.rds"))
qq <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- cs_all[, .(n = .N, lloq_never_pct = 100 * mean(is.na(t_lloq)),
                   lloq_studyday_median = median(t_lloq + 1, na.rm = TRUE), lloq_studyday_p05 = qq(t_lloq + 1, 0.05), lloq_studyday_p95 = qq(t_lloq + 1, 0.95),
                   lloq_after_day58_pct = 100 * mean(t_lloq + 1 > 58, na.rm = TRUE),
                   c_start1_median = median(c_start1, na.rm = TRUE), c_start2_median = median(c_start2, na.rm = TRUE),
                   len1_median = median(len1, na.rm = TRUE), len1_p05 = qq(len1, 0.05), len1_p95 = qq(len1, 0.95),
                   len2_median = median(len2, na.rm = TRUE), len2_p05 = qq(len2, 0.05), len2_p95 = qq(len2, 0.95),
                   tmax_median = median(tmax)), by = .(model, weight)]
fwrite(summ, file.path(out_dir, "cliff_summary.csv")); fwrite(pts, file.path(out_dir, "cliff_points.csv"))
fwrite(data.table(schedule = names(SCHED), schedule_label = unname(SCHED_LABEL[names(SCHED)]), min_interval_day_from_day22 = unname(min_interval),
                  added_study_days = vapply(SCHED, function(s) paste(setdiff(s, b0) + 1, collapse = " "), character(1))), file.path(out_dir, "cliff_schedules.csv"))
print(summ, digits = 3); print(pts[weight == "base" & timing == "nominal"], digits = 3)

source(proj_path("scripts", "34b_cliff_figures.R"))   # 그림 2-1 ~ 2-4 (저장된 결과에서 다시 그릴 수 있게 분리)
append_run_log(logfile, "done")
