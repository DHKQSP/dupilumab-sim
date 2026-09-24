#!/usr/bin/env Rscript
# §2 보고 문구를 절벽 분석 결과로 생성(한국어·영문). 검토자 문구의 전제(구간 길이 < 최소 방문 간격, 고정 일정 2점 이상 0%, 매일 채혈 3점 0%)를 검사하고,
# 어긋나는 부분은 결과대로 바꿔 쓴다(방문 허용창 결과 포함). 기준값(독립 구현) 비교표도 만든다.
source("R/00_setup.R"); source_project()
d <- proj_path("results", "cliff")
s <- fread(file.path(d, "cliff_summary.csv")); pt <- fread(file.path(d, "cliff_points.csv")); sch <- fread(file.path(d, "cliff_schedules.csv"))
b <- s[model == "k2016" & weight == "base"]; b20 <- s[model == "k2020" & weight == "base"]
fixed <- c("current", "plus_39_46", "plus_39_46_53", "plus_32_39_46_53", "plus_40_47")
min_fixed <- min(sch[schedule %in% fixed & schedule != "current", min_interval_day_from_day22])
nom <- pt[weight == "base" & timing == "nominal"]; win <- pt[weight == "base" & timing == "windowed"]
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)
len_ok <- max(s[weight == "base", len2_p95]) < min_fixed                  # 2일 정의 95백분위도 최소 간격보다 짧은가
cs <- readRDS(file.path(d, "cliff_subjects.rds"))[weight == "base"]
pct_len2_short <- cs[, 100 * mean(len2 < min_fixed, na.rm = TRUE), by = model]; pct_len1_short <- cs[, 100 * mean(len1 < min_fixed, na.rm = TRUE), by = model]
ge2_nom <- max(nom[schedule %in% fixed, pct_ge2]); ge2_nom_1d <- max(nom[schedule %in% fixed & definition_day == 1, pct_ge2])
ge2_win_1d <- max(win[schedule %in% fixed & definition_day == 1, pct_ge2]); ge2_win_2d <- max(win[schedule %in% fixed & definition_day == 2, pct_ge2])
ge3_daily_1d_nom <- max(nom[schedule == "daily_29_57" & definition_day == 1, pct_ge3]); ge3_daily_1d_win <- max(win[schedule == "daily_29_57" & definition_day == 1, pct_ge3])
ge2_nom_2d <- max(nom[schedule %in% fixed & definition_day == 2, pct_ge2])
which_2d <- nom[schedule %in% fixed & definition_day == 2][which.max(pct_ge2), paste(c(k2016 = "2016", k2020 = "Model 1")[model], schedule_label)]
dl <- nom[model == "k2016" & schedule == "daily_29_57" & definition_day == 1]
out_before <- cs[model == "k2016", 100 * mean(start1 < 28, na.rm = TRUE)]; out_after <- cs[model == "k2016", 100 * mean(t_lloq > 56, na.rm = TRUE)]
ko <- c(sprintf("절벽 구간 길이는 1일 정의 중앙값 %s일(5–95백분위 %s–%s), 2일 정의 %s일(%s–%s)이다(2016 모델, 60–90 kg, 20,000명; Model 1 %s일, %s일). 추가 채혈 일정의 Day 22 이후 최소 방문 간격 %g일보다 짧은 대상자는 1일 정의 %s%%, 2일 정의 %s%%(2016).",
              f2(b$len1_median), f2(b$len1_p05), f2(b$len1_p95), f2(b$len2_median), f2(b$len2_p05), f2(b$len2_p95), f2(b20$len1_median), f2(b20$len2_median), min_fixed,
              f1(pct_len1_short[model == "k2016", V1]), f1(pct_len2_short[model == "k2016", V1])),
        sprintf("명목 채혈일에서 고정 방문 일정이 한 대상자의 절벽 구간에 2점 이상 찍는 비율은 1일 정의 %s%%%s, 2일 정의 최대 %s%%%s다. 방문 허용창(±1일)을 적용하면 1일 정의 최대 %s%%, 2일 정의 최대 %s%%(허용창 때문에 두 방문이 절벽 길이 안으로 붙는 경우).",
                f1(ge2_nom_1d), if (ge2_nom_1d == 0) "(구조적으로 불가능)" else "", f1(ge2_nom_2d), if (ge2_nom_2d > 0) sprintf("(최대: %s)", which_2d) else "",
                f1(ge2_win_1d), f1(ge2_win_2d)),
        sprintf("λz 산출에 필요한 3점은 Day 29–57 매일 채혈에서도 1일 정의 기준 명목일 %s%%(방문 허용창 적용 시 %s%%)다. 설령 절벽에서 λz를 구해도 순간 반감기가 0에 가까운 기울기라 외삽 면적은 0에 수렴하고 AUCinf는 AUClast로 수렴한다.",
                f1(ge3_daily_1d_nom), f1(ge3_daily_1d_win)),
        sprintf("Day 29–57 매일 채혈에서 절벽에 1점 이상은 전체 대상자 기준 %s%%, 절벽이 Day 29–57 안에 있는 대상자 기준 %s%%다(2016, 1일 정의). 차이는 절벽이 Day 29 이전에 시작하는 %s%%와 LLOQ 도달이 Day 57 이후인 %s%% 때문이다.",
                f1(dl$pct_ge1), f1(dl$pct_ge1_cliff_in_day29_57), f1(out_before), f1(out_after)),
        sprintf("정량한계 도달일 중앙값 Day %s(5–95백분위 Day %s–%s), Day 58 이후 %s%%; 절벽 시작 농도 중앙값 %s mg/L(2016 모델, 60–90 kg). 비만 밴드(90–150 kg)에서도 절벽 길이 중앙값 %s–%s일.",
                f1(b$lloq_studyday_median), f1(b$lloq_studyday_p05), f1(b$lloq_studyday_p95), f1(b$lloq_after_day58_pct), f2(b$c_start1_median),
                f2(min(s[weight != "base", len1_median])), f2(max(s[weight != "base", len1_median]))))
en <- c(sprintf("The cliff (instantaneous half-life below 1 day until the true concentration reaches the LLOQ) lasts a median %s days (5th to 95th percentile %s to %s) with the 1-day definition and %s days (%s to %s) with the 2-day definition (2016 model, 60 to 90 kg, 20,000 subjects; Model 1 %s and %s days). It is shorter than the minimum visit interval after Day 22 of any added-sampling schedule (%g days) in %s%% (1-day) and %s%% (2-day) of subjects.",
              f2(b$len1_median), f2(b$len1_p05), f2(b$len1_p95), f2(b$len2_median), f2(b$len2_p05), f2(b$len2_p95), f2(b20$len1_median), f2(b20$len2_median), min_fixed,
              f1(pct_len1_short[model == "k2016", V1]), f1(pct_len2_short[model == "k2016", V1])),
        sprintf("At nominal sampling days, the share of subjects with two or more samples in the cliff under a fixed schedule is %s%% with the 1-day definition%s and at most %s%% with the 2-day definition. With visit windows of plus or minus 1 day it is at most %s%% (1-day) and %s%% (2-day).",
                f1(ge2_nom_1d), if (ge2_nom_1d == 0) " (structurally impossible)" else "", f1(ge2_nom_2d), f1(ge2_win_1d), f1(ge2_win_2d)),
        sprintf("Three points, the minimum for lambda-z, occur in %s%% of subjects even with daily sampling from Day 29 to Day 57 (1-day definition, nominal days; %s%% with visit windows). Even if lambda-z were estimated on the cliff, the instantaneous half-life is close to zero, so the extrapolated area tends to zero and AUCinf tends to AUClast.",
                f1(ge3_daily_1d_nom), f1(ge3_daily_1d_win)),
        sprintf("With daily sampling from Day 29 to Day 57, %s%% of all subjects and %s%% of subjects whose cliff lies within Day 29 to Day 57 have at least one sample in the cliff (2016 model, 1-day definition); the difference comes from cliffs starting before Day 29 (%s%%) or ending after Day 57 (%s%%).",
                f1(dl$pct_ge1), f1(dl$pct_ge1_cliff_in_day29_57), f1(out_before), f1(out_after)))
# 전제 검사: 문구에 "구조적으로 불가능"을 쓰면 명목일 고정 일정 2점 이상이 모두 0이어야 한다
if (any(grepl("구조적으로 불가능", ko))) stopifnot(ge2_nom_1d == 0)
writeLines(ko, file.path(d, "cliff_conclusion_ko.md")); writeLines(en, file.path(d, "cliff_conclusion_en.md"))
# 독립 구현 기준값 비교(2016, 60–90 kg, 명목일)
ref <- read_cfg("crossval_reference.yaml")$cliff_round6
if (!is.null(ref)) {
  rr <- rbindlist(list(
    data.table(metric = c("LLOQ 도달일 중앙값", "5백분위", "95백분위", "Day 58 이후 %", "절벽 시작 농도 중앙값", "구간 길이 1일 정의", "구간 길이 2일 정의"),
               ref = c(ref$lloq_day_median, ref$lloq_day_p05, ref$lloq_day_p95, ref$after_day58_pct, ref$c_start_median, ref$len1_median, ref$len2_median),
               sim = c(b$lloq_studyday_median, b$lloq_studyday_p05, b$lloq_studyday_p95, b$lloq_after_day58_pct, b$c_start1_median, b$len1_median, b$len2_median)),
    rbindlist(lapply(names(ref$points), function(k) { r <- ref$points[[k]]; x <- nom[model == "k2016" & schedule == r$schedule & definition_day == r$def]
      data.table(metric = sprintf("%s, %d일 정의: 1점/2점/3점 이상 %%", r$schedule, r$def), ref = NA_real_, sim = NA_real_,
                 ref_txt = paste(unlist(r$pct), collapse = " / "), sim_txt = paste(f1(c(x$pct_ge1, x$pct_ge2, x$pct_ge3)), collapse = " / "),
                 sim_in_window_txt = paste(f1(c(x$pct_ge1_cliff_in_day29_57, x$pct_ge2_cliff_in_day29_57, x$pct_ge3_cliff_in_day29_57)), collapse = " / ")) }))), fill = TRUE)
  fwrite(rr, file.path(d, "cliff_reference_comparison.csv")); print(rr)
}
cat(ko, sep = "\n\n")
