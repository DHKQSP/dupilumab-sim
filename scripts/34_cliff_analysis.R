#!/usr/bin/env Rscript
# §2 절벽 채혈 무의미성 분석 (config/oc_design.yaml cliff, D-041).
# 절벽 = tmax 이후 순간 반감기 ln2/(−d ln C/dt)가 처음 1일(민감도 2일) 미만이 되는 시점부터 참 농도가 LLOQ(0.078 mg/L)에 닿는 시점까지.
# d ln C/dt는 ODE 상태값으로 정확히 계산: dC/dt = (ka·흡수구획 − (ke + k12)·central + k21·periph − Vmax·C/(Km + C)·Vc)/Vc. 교차 시점은 0.05일 격자 사이 선형 보간.
# 조건: 두 모델 × 체중(60–90 kg 기본, 90–110, 110–130, 130–150 kg 균등) × 20,000명 × 채혈 시각(명목, 방문 허용창).
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff; design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", if (length(commandArgs(trailingOnly = TRUE))) "cliff_test" else "cliff"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args)) as.integer(args[1]) else as.integer(cf$n_subjects); LLOQ <- as.numeric(cf$lloq); STEP <- as.numeric(cf$grid_step_day); SEED <- as.integer(cf$seed)
defs <- as.numeric(unlist(cf$definition_days))
logfile <- start_run_log("cliff_analysis", master_seed = SEED, run_mode = "final", extra = list(n = N, step = STEP))
b0 <- get_schedule(design, "B0")
SCHED <- list(current = b0, plus_39_46 = sort(c(b0, 38, 45)), plus_39_46_53 = sort(c(b0, 38, 45, 52)),
              plus_32_39_46_53 = sort(c(b0, 31, 38, 45, 52)), plus_40_47 = sort(c(b0, 39, 46)), daily_29_57 = sort(unique(c(b0, 28:56))))
SCHED_LABEL <- c(current = "현행", plus_39_46 = "+Day 39·46", plus_39_46_53 = "+Day 39·46·53", plus_32_39_46_53 = "+Day 32·39·46·53",
                 plus_40_47 = "+Day 40·47", daily_29_57 = "Day 29–57 매일")
min_interval <- vapply(SCHED, function(s) { x <- s[s >= 21]; min(diff(x)) }, numeric(1))    # 경과일 21(Day 22) 이후 최소 방문 간격

wt_spec_for <- function(w) if (identical(w, "base")) weight_spec_from_design(design, "base") else list(dist = w$dist, trunc = as.numeric(unlist(w$trunc)))

# 한 모델·체중군: 20,000명 → 대상자별 tmax, 절벽 시작(정의별), LLOQ 도달 시각, 절벽 시작 농도
cliff_subjects <- function(model, wname) {
  p <- load_params(model); mod <- get_model(p$model_id)
  subj <- with_seed(derive_seed(SEED, model, wname, "subj"), make_subjects(N, p, wt_spec_for(cf$weights[[wname]]), 0.5, 0))
  ip <- individual_params(p, subj)
  times <- seq(STEP, 100, by = STEP)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  res <- list(); chunk <- 1000L
  for (b in seq_len(ceiling(N / chunk))) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, N)
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    ev <- rxode2::add.sampling(rxode2::et(amt = design$dose_mg, cmt = "depot", time = 0), times)
    s <- rxode2::rxSolve(mod, params = prm, events = ev, atol = 1e-10, rtol = 1e-8, maxsteps = 2000000L, cores = 1L, returnType = "data.table", addDosing = FALSE)
    s[, id := prm$id[sim.id]]
    s <- merge(s, ip[rows, .(id, Vc_i, ke_i = ke, k12_i = k12, k21_i = k21, Vmax_i = Vmax, Km_i = Km, ka_i = ka)], by = "id")
    absorb <- if ("absc" %in% names(s)) s$absc else s$depot
    s[, dCdt := (ka_i * absorb - (ke_i + k12_i) * central + k21_i * periph - Vmax_i * C / (Km_i + C) * Vc_i) / Vc_i]
    s[, r := -dCdt / C]                                                   # −d ln C/dt
    setorder(s, id, time)
    res[[b]] <- s[, {
      i <- which.max(C); tt <- time; cc <- C; rr <- r; n <- .N
      post <- (i + 1):n
      kl <- post[which(cc[post] < LLOQ)[1]]
      t_lloq <- if (is.na(kl)) NA_real_ else tt[kl - 1] + (log(LLOQ) - log(cc[kl - 1])) / (log(cc[kl]) - log(cc[kl - 1])) * (tt[kl] - tt[kl - 1])
      st <- vapply(defs, function(dd) {
        thr <- log(2) / dd; k <- post[which(rr[post] >= thr)[1]]
        if (is.na(k)) NA_real_ else if (k == i + 1) tt[k] else tt[k - 1] + (thr - rr[k - 1]) / (rr[k] - rr[k - 1]) * (tt[k] - tt[k - 1]) }, numeric(1))
      cst <- vapply(st, function(x) if (is.na(x)) NA_real_ else exp(approx(tt, log(cc), x)$y), numeric(1))
      .(tmax = tt[i], Cmax = cc[i], t_lloq = t_lloq, start1 = st[1], start2 = st[2], c_start1 = cst[1], c_start2 = cst[2])
    }, by = id]
  }
  out <- rbindlist(res)
  out[, `:=`(len1 = pmax(t_lloq - start1, 0), len2 = pmax(t_lloq - start2, 0))]
  out[, `:=`(model = model, weight = wname, WT = subj$WT[match(id, subj$id)])][]
}

count_points <- function(cs, timing, model, wname) {
  rbindlist(lapply(names(SCHED), function(sn) {
    days <- SCHED[[sn]]
    ob <- with_seed(derive_seed(SEED, model, wname, sn, "jitter"), make_obs_times(cs$id, days, design, jitter = timing == "windowed"))[planned > 0]
    m <- merge(ob, cs[, .(id, start1, start2, t_lloq)], by = "id")
    rbindlist(lapply(seq_along(defs), function(k) {
      stc <- if (k == 1) "start1" else "start2"
      cnt <- m[, .(n = sum(!is.na(t_lloq) & time >= get(stc) & time < t_lloq), st = get(stc)[1], tl = t_lloq[1]), by = id]
      w1 <- wilson_ci(sum(cnt$n >= 1), nrow(cnt)); w2 <- wilson_ci(sum(cnt$n >= 2), nrow(cnt)); w3 <- wilson_ci(sum(cnt$n >= 3), nrow(cnt))
      inw <- cnt[!is.na(tl) & st >= 28 & tl <= 56]                     # 절벽이 Day 29–57(경과일 28–56) 안에 있는 대상자(매일 채혈 창과 비교용)
      data.table(model = model, weight = wname, timing = timing, definition_day = defs[k], schedule = sn, schedule_label = SCHED_LABEL[[sn]], n_subjects = nrow(cnt),
                 min_interval_day = min_interval[[sn]], pct_ge1 = w1$est, ge1_lo = w1$lo, ge1_hi = w1$hi, pct_ge2 = w2$est, ge2_lo = w2$lo, ge2_hi = w2$hi,
                 pct_ge3 = w3$est, ge3_lo = w3$lo, ge3_hi = w3$hi, n_cliff_in_day29_57 = nrow(inw),
                 pct_ge1_cliff_in_day29_57 = 100 * mean(inw$n >= 1), pct_ge2_cliff_in_day29_57 = 100 * mean(inw$n >= 2), pct_ge3_cliff_in_day29_57 = 100 * mean(inw$n >= 3))
    }))
  }))
}

all_cs <- list(); pts <- list()
for (model in c("k2016", "k2020")) for (wname in names(cf$weights)) {
  t0 <- Sys.time(); cs <- cliff_subjects(model, wname)
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
mi_u <- mi[, .(lab = paste(schedule_label, collapse = ", ")), by = mi]
f23 <- ggplot(ld, aes(len, fill = def)) + geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.85) +
  geom_vline(data = mi_u, aes(xintercept = mi), colour = VIZ$ink, linetype = "22", linewidth = 0.5) +
  geom_text(data = mi_u, aes(x = mi, y = Inf, label = sprintf("최소 방문 간격 %g일\n(%s)", mi, lab)), inherit.aes = FALSE, vjust = 1.2, hjust = -0.03, size = 2.6, colour = VIZ$ink) +
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
sp <- merge(sp, prof[, .(id, time, C)], by = c("id", "time"), all.x = TRUE)
sp[is.na(C), C := vapply(seq_len(.N), function(k) approx(prof[id == sp$id[k]]$time, prof[id == sp$id[k]]$C, sp$time[k])$y, numeric(1))]
sp[, who := factor(lab_rep[as.character(id)], levels = unname(lab_rep))]
f24 <- ggplot(prof[C > 0.01], aes(time, C)) +
  geom_rect(data = shade, aes(xmin = start1, xmax = t_lloq, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = VIZ$s2, alpha = 0.18) +
  geom_line(colour = VIZ$ink2, linewidth = 0.5) + geom_hline(yintercept = LLOQ, colour = VIZ$muted, linetype = "22") +
  geom_point(data = sp[C > 0.01], aes(shape = set, colour = set), size = 2.6) +
  scale_shape_manual(values = c(16, 17), name = NULL) + scale_colour_manual(values = c(VIZ$s1, VIZ$s3), name = NULL) +
  scale_y_log10() + facet_wrap(~who, ncol = 3) +
  labs(x = "투여 후 경과일", y = "참 농도 (mg/L, 로그)", title = "그림 2-4. 대표 대상자의 농도 곡선과 절벽 구간(음영, 1일 정의)",
       subtitle = cap("k2016", "base", ", 점선 = LLOQ 0.078 mg/L")) + theme_dupi()
ggsave(file.path(out_dir, "fig2_4_representative.png"), f24, width = 11, height = 4.2, dpi = 120)
append_run_log(logfile, "done")
