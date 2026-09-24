#!/usr/bin/env Rscript
# §8 문헌 정합성 표: 출처 | 조건 | 문헌 값 | 같은 조건 시뮬레이션 값(비구획·참값, 두 모델) | 판정. 정성 행 포함.
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "literature"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cv <- rbind(fread(proj_path("results", "step1", "step1g_literature_coverage.csv"))[, model := "k2016"], fread(proj_path("results", "step1_k2020", "step1g_literature_coverage.csv"))[, model := "k2020"], fill = TRUE)
f <- function(x) sprintf("%.1f%%", 100 * x)
lit <- data.table(
  item = c("clot200", "clot300", "clot600", "pkm12350"),
  source = c("Clot 2021 Table 3", "Clot 2021 Table 3", "Clot 2021 Table 3", "FDA BLA 761055 Clin Pharm Review Table 4.2.c (PKM12350 대조군)"),
  condition = c("200 mg, Day 2–57 채혈, LLOQ 0.078 mg/L, 체중 62.2 kg (제형 차이: 외부 점검 성격)", "300 mg, 같은 채혈, 체중 59.9 kg", "600 mg, 같은 채혈, 체중 58.3 kg", "300 mg, PKM12350 채혈, 체중 78 kg(확정 전 가정)"),
  dose_mg = c(200, 300, 600, 300),
  lit_mean_ratio = c(402 / 425, 792 / 807, 2110 / 2150, 500 / 521), lit_geo_ratio = c(391 / 418, 710 / 726, 2060 / 2100, NA),
  lit_tlast = c("28 (21–35)", "42 (28–49)", "56 (56–56)", NA),
  note = c("", "", "", "시험군 488 대 554(88.1%)는 AUC0-inf 산출 대상자 수가 다를 가능성(공개 자료로 확인 불가)"))
rows <- rbindlist(lapply(seq_len(nrow(lit)), function(i) { L <- lit[i]
  s <- cv[if (L$item == "pkm12350") grepl("PKM12350", source) else (grepl("Clot", source) & dose_mg == L$dose_mg)]
  g <- function(m, col) { v <- s[model == m][[col]]; if (length(v)) v else NA_real_ }
  tl <- function(m) { x <- s[model == m]; if (!nrow(x) || is.na(x$tlast_cohort_median)) sprintf("대상자 %g (5–95%% %g–%g)", x$tlast_median_subj, x$tlast_p05_subj, x$tlast_p95_subj) else sprintf("코호트 중앙값 %g (5–95%% %g–%g); 대상자 %g", x$tlast_cohort_median, x$tlast_cohort_p05, x$tlast_cohort_p95, x$tlast_median_subj) }
  data.table(source = L$source, condition = L$condition, dose_mg = L$dose_mg,
             literature = paste0("평균비 ", f(L$lit_mean_ratio), if (!is.na(L$lit_geo_ratio)) paste0(", 기하 ", f(L$lit_geo_ratio)) else "", if (!is.na(L$lit_tlast)) paste0("; tlast ", L$lit_tlast, "일") else ""),
             sim_k2016 = sprintf("비구획 평균비 %s (신뢰군 %s), 기하 %s; 참값 평균비 %s; tlast %s", f(g("k2016", "mean_ratio_nca_all")), f(g("k2016", "mean_ratio_nca_reliable")), f(g("k2016", "geo_ratio_nca_all")), f(g("k2016", "mean_ratio_true")), tl("k2016")),
             sim_k2020 = sprintf("비구획 평균비 %s (신뢰군 %s), 기하 %s; 참값 평균비 %s; tlast %s", f(g("k2020", "mean_ratio_nca_all")), f(g("k2020", "mean_ratio_nca_reliable")), f(g("k2020", "geo_ratio_nca_all")), f(g("k2020", "mean_ratio_true")), tl("k2020")),
             nca_mean_ratio_k2016 = g("k2016", "mean_ratio_nca_all"), true_mean_ratio_k2016 = g("k2016", "mean_ratio_true"), nca_mean_ratio_k2020 = g("k2020", "mean_ratio_nca_all"), true_mean_ratio_k2020 = g("k2020", "mean_ratio_true"),
             lit_mean_ratio = L$lit_mean_ratio, note = L$note)
}))
rows[, judgment := fifelse(abs(nca_mean_ratio_k2016 - lit_mean_ratio) <= 0.03 & abs(nca_mean_ratio_k2020 - lit_mean_ratio) <= 0.03, "일치(비구획 평균비 ±3%p)",
                   fifelse(nca_mean_ratio_k2016 > lit_mean_ratio | nca_mean_ratio_k2020 > lit_mean_ratio, "모의가 높음", "모의가 낮음(AUClast 커버리지 과소평가, 보수적 방향)"))]
rows[grepl("200 mg", condition), judgment := paste0(judgment, " — 제형 차이로 외부 점검 성격")]
qual <- data.table(source = c("Kovalenko 2020 Results", "Kovalenko 2020 Discussion, Supplementary Figure 1, Figure 4", "Kovalenko 2021 Discussion (p.1353–1354)",
                              "Kovalenko 2016 Methods", "Li 2020 p.750", "Cohen 2022 p.677–678", "Clot 2021 Discussion"),
                   statement = c("LLOQ에 가까워질수록 말단 기울기가 음의 무한대로 가며, 표적 매개 구간은 거의 수직선이다.",
                                 "의미 있는 말단 반감기는 계산할 수 없고 규제기관에 제출하지 않았다. 순간 반감기는 0으로 감소한다. 최대 총 청소율이 선형 청소율보다 훨씬 커서 가파른 구간이 생긴다.",
                                 "Km 0.01 mg/L에서 0.09 mg/L 농도는 순환 표적의 90%를 제거한다. 순간 반감기는 베타기 약 25일에서 표적 매개기 0 근처로 변한다.",
                                 "가파른 표적 매개 구간을 기술할 저농도 정량값이 적어 BLQ 값을 M3로 사용했다.",
                                 "농도가 낮을수록 소실 기울기가 가팔라지고 AUClast가 용량 비례 이상으로 증가한다.",
                                 "말단 비선형 때문에 반감기와 AUCinf를 산출하지 않았다.",
                                 "저농도에서 표적 매개로 소실이 빨라지는 다지수 감소."))
interp <- "문헌의 커버리지 값은 비구획 AUCinf 기준이라 외삽이 부풀려져 있으므로 실제 커버리지의 하한이다. 정량한계 아래 곡선 모양은 관측으로 확인할 수 없으므로 곡선 모양 민감도(Km·Vmax)로 대응한다."
x6 <- rows[dose_mg == 600]
if (nrow(x6) && x6$true_mean_ratio_k2016 < x6$lit_mean_ratio)
  interp <- paste(interp, sprintf("600 mg에서는 모의 참값 커버리지(%s, %s)가 문헌 비구획 값(%s)보다 낮다. 문헌 값이 하한이므로 두 모델은 이 용량에서 Day 56 이후 꼬리를 관측보다 크게 예측하며, AUClast 커버리지를 과소평가하는 보수적 방향이다. 연구 용량 300 mg에서는 비구획 평균비 차이가 %.1f%%p로 ±3%%p 이내다.",
    f(x6$true_mean_ratio_k2016), f(x6$true_mean_ratio_k2020), f(x6$lit_mean_ratio), 100 * abs(rows[dose_mg == 300 & grepl("Clot", source), nca_mean_ratio_k2016] - rows[dose_mg == 300 & grepl("Clot", source), lit_mean_ratio])))
fwrite(rows, file.path(out_dir, "literature_numeric.csv")); fwrite(qual, file.path(out_dir, "literature_qualitative.csv")); writeLines(interp, file.path(out_dir, "literature_interpretation.txt"))
print(rows[, .(source, lit_mean_ratio, nca_mean_ratio_k2016, true_mean_ratio_k2016, nca_mean_ratio_k2020, true_mean_ratio_k2020, judgment)], digits = 3)
