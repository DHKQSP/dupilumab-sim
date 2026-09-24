#!/usr/bin/env Rscript
# §2 그림 2-1 ~ 2-4: results/cliff/의 저장 결과에서 그린다(34_cliff_analysis.R 끝에서 호출, 단독 실행 가능).
if (!exists("cs_all")) {
  source("R/00_setup.R"); source_project(); suppressPackageStartupMessages(library(ggplot2))
  oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff; design <- read_cfg("trial_design.yaml")
  out_dir <- proj_path("results", "cliff"); N <- as.integer(cf$n_subjects); LLOQ <- as.numeric(cf$lloq); SEED <- as.integer(cf$seed)
  cs_all <- readRDS(file.path(out_dir, "cliff_subjects.rds")); pts <- fread(file.path(out_dir, "cliff_points.csv"))
  b0 <- get_schedule(design, "B0")
  SCHED_LABEL <- c(current = "현행", plus_39_46 = "+Day 39·46", plus_39_46_53 = "+Day 39·46·53", plus_32_39_46_53 = "+Day 32·39·46·53", plus_40_47 = "+Day 40·47", daily_29_57 = "Day 29–57 매일")
  sch <- fread(file.path(out_dir, "cliff_schedules.csv")); min_interval <- setNames(sch$min_interval_day_from_day22, sch$schedule)
  wt_spec_for <- function(w) if (identical(w, "base")) weight_spec_from_design(design, "base") else list(dist = w$dist, trunc = as.numeric(unlist(w$trunc)))
}
# ----- 그림 2-1 ~ 2-4 (2016 주 모델, 60–90 kg, 명목일; 캡션에 모델·조건·대상자 수) ------------------------------------
cap <- function(model, w, extra = "") sprintf("%s, %s, 가상 대상자 %s명%s", if (model == "k2016") "Kovalenko 2016 (주)" else "Kovalenko 2020 Model 1",
                                             if (w == "base") "60–90 kg" else w, format(N, big.mark = ","), extra)
base_cs <- cs_all[model == "k2016" & weight == "base" & !is.na(t_lloq)]
cur_days <- c(22, 29, 36, 43, 50, 57)
h <- base_cs[, .(day = t_lloq + 1)]
h[, bin := floor(day / 2) * 2]
hb <- h[, .N, by = bin][, pct := 100 * N / nrow(h)][order(bin)]
hb[, has_visit := vapply(bin, function(b) any(cur_days >= b & cur_days < b + 2), logical(1))]
f21 <- ggplot(hb, aes(bin + 1, pct, fill = has_visit)) + geom_col(width = 1.8) +
  scale_fill_manual(values = c(`TRUE` = VIZ$s2, `FALSE` = VIZ$s1), labels = c(`TRUE` = "현행 채혈일(Day 22·29·36·43·50·57)이 든 구간", `FALSE` = "그 외"), name = NULL) +
  labs(x = "참 농도가 LLOQ(0.078 mg/L)에 닿는 연구일 (2일 간격)", y = "대상자 비율 (%)", title = "그림 2-1. 정량한계 도달일 분포",
       subtitle = cap("k2016", "base", ", 명목 채혈일")) + theme_dupi()
ggsave(file.path(out_dir, "fig2_1_lloq_day.png"), f21, width = 8, height = 4.2, dpi = 120)

pp <- pts[model == "k2016" & weight == "base"]
pl <- melt(pp, id.vars = c("timing", "definition_day", "schedule_label"), measure.vars = c("pct_ge1", "pct_ge2", "pct_ge3"), variable.name = "k", value.name = "pct")
pl[, k := factor(k, levels = c("pct_ge1", "pct_ge2", "pct_ge3"), labels = c("1점 이상", "2점 이상", "3점 이상"))]
pl[, schedule_label := factor(schedule_label, levels = unname(SCHED_LABEL))]
pl[, panel := sprintf("%s 정의 · %s", paste0(definition_day, "일"), fifelse(timing == "nominal", "명목일", "방문 허용창"))]
f22 <- ggplot(pl, aes(schedule_label, pct, fill = k)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.4, colour = VIZ$ink2) +
  scale_fill_manual(values = c(VIZ$s1, VIZ$s2, VIZ$s3), name = NULL) + facet_wrap(~panel, ncol = 2) +
  labs(x = NULL, y = "대상자 비율 (%)", title = "그림 2-2. 절벽 구간에 찍힌 채혈점 수", subtitle = cap("k2016", "base")) +
  theme_dupi() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "fig2_2_points_in_cliff.png"), f22, width = 10, height = 6.5, dpi = 120)

ld <- rbind(base_cs[, .(len = len1, def = "1일 정의")], base_cs[, .(len = len2, def = "2일 정의")])
mi <- data.table(schedule_label = unname(SCHED_LABEL), mi = unname(min_interval))[, lab := sprintf("%s: %g일", schedule_label, mi)]
mi_u <- mi[, .(lab = paste(schedule_label, collapse = ", ")), by = mi][order(mi)]
ymax <- max(hist(ld$len, breaks = seq(0, max(ld$len, na.rm = TRUE) + 0.05, by = 0.05), plot = FALSE)$counts)
mi_u[, y := ymax * c(1.02, 0.80, 1.02)[seq_len(.N)]]
f23 <- ggplot(ld, aes(len, fill = def)) + geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.85) +
  geom_vline(data = mi_u, aes(xintercept = mi), colour = VIZ$ink, linetype = "22", linewidth = 0.5) +
  geom_label(data = mi_u, aes(x = mi, y = y, label = sprintf("최소 방문 간격 %g일\n(%s)", mi, lab)), inherit.aes = FALSE, vjust = 1, hjust = -0.03, size = 2.5, colour = VIZ$ink, fill = VIZ$surface, label.size = 0) +
  scale_fill_manual(values = c(VIZ$s1, VIZ$s2), name = NULL) + coord_cartesian(xlim = c(0, 8)) +
  labs(x = "절벽 구간 길이 (일)", y = "대상자 수", title = "그림 2-3. 절벽 구간 길이와 일정별 최소 방문 간격(Day 22 이후)", subtitle = cap("k2016", "base")) + theme_dupi()
ggsave(file.path(out_dir, "fig2_3_cliff_length.png"), f23, width = 9, height = 4.4, dpi = 120)

# 대표 대상자 3명(LLOQ 도달일 25, 50, 75백분위에 가장 가까운 대상자)
qs <- quantile(base_cs$t_lloq, c(0.25, 0.5, 0.75))
rep_ids <- vapply(qs, function(q) base_cs$id[which.min(abs(base_cs$t_lloq - q))], numeric(1))
p16 <- load_params("k2016")
subj16 <- with_seed(derive_seed(SEED, "k2016", "base", "subj"), make_subjects(N, p16, wt_spec_for(cf$weights$base), 0.5, 0))
iprep <- individual_params(p16, subj16[id %in% rep_ids])
prof <- solve_model(iprep, CJ(id = rep_ids, time = seq(0.05, 70, by = 0.05)), design$dose_mg, p16$model_id)
lab_rep <- setNames(sprintf("%s백분위 대상자 (LLOQ 도달 Day %.1f)", c("25", "50", "75"), base_cs[match(rep_ids, id), t_lloq] + 1), rep_ids)
prof[, who := factor(lab_rep[as.character(id)], levels = unname(lab_rep))]
shade <- base_cs[id %in% rep_ids, .(id, start1, t_lloq)][, who := factor(lab_rep[as.character(id)], levels = unname(lab_rep))]
sp <- rbind(CJ(id = rep_ids, time = b0)[, set := "현행 채혈점"], CJ(id = rep_ids, time = c(38, 45, 52))[, set := "+Day 39·46·53"])
sp[, C := { pr <- prof[id == .BY$id]; exp(approx(pr$time, log(pr$C), time)$y) }, by = id]   # 각 대상자 곡선에서 로그 보간
sp[, who := factor(lab_rep[as.character(id)], levels = unname(lab_rep))]
f24 <- ggplot(prof[C > 0.01], aes(time, C)) +
  geom_rect(data = shade, aes(xmin = start1, xmax = t_lloq, ymin = 0.01, ymax = Inf), inherit.aes = FALSE, fill = VIZ$s2, alpha = 0.35) +
  geom_vline(data = shade, aes(xintercept = start1), colour = VIZ$s2, linewidth = 0.4) +
  geom_line(colour = VIZ$ink2, linewidth = 0.5) + geom_hline(yintercept = LLOQ, colour = VIZ$muted, linetype = "22") +
  geom_point(data = sp[C > 0.01], aes(shape = set, colour = set), size = 2.6) +
  scale_shape_manual(values = c(16, 17), name = NULL) + scale_colour_manual(values = c(VIZ$s1, VIZ$s3), name = NULL) +
  scale_y_log10() + coord_cartesian(xlim = c(0, 60)) + facet_wrap(~who, ncol = 3) +
  labs(x = "투여 후 경과일", y = "참 농도 (mg/L, 로그)", title = "그림 2-4. 대표 대상자의 농도 곡선과 절벽 구간(음영, 1일 정의)",
       subtitle = cap("k2016", "base", ", 점선 = LLOQ 0.078 mg/L")) + theme_dupi()
ggsave(file.path(out_dir, "fig2_4_representative.png"), f24, width = 11, height = 4.2, dpi = 120)
