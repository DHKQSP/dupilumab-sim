# 보강 A-①(v1.2: 본문 S6에서 별첨으로 옮김; 표시 이름 "보강 A-①"): 왼쪽 = 절벽 폭 대 방문 간격 개념도(현행 7일, 추가 후보의 최소 간격 3일;
# S4 대표 대상자 = 2020 모델 절벽 대표 대상자, 사전 등록 config/prereg_20260929.yaml 7g, results/core_deck/cliff_rep_k2020.csv)와
# 채혈 추가 시 세트 (iii) 충족 비율 변화(D1~D4 - B0, 같은 대상자, 두 모델; results/core_deck/reliable_iii_by_schedule.csv),
# 오른쪽 = 세트 (iii) 탈락자 대 유지자의 참 AUCinf 기하평균비(95% CI, 두 모델; tp_characteristics.csv). 본문 2줄.
slide_B1 <- function() {
  RL <- "core_deck/reliable_iii_by_schedule.csv"; TCH <- "trialpop/tp_characteristics.csv"; CP <- "cliff/cliff_points.csv"; CS <- "cliff/cliff_summary.csv"; RS <- "core_deck/cliff_rep_k2020.csv"   # S4와 같은 대상자(2020 모델)
  deck_slide("B1", tag = "sim")
  d <- copy(rows(RL, "schedule!='B0'")); premise(nrow(d) == 8, "D1-D4 x two models")
  rng <- dderived("change in the share meeting set (iii), D1-D4 vs B0, range over schedules and models (percentage points)", RL, "schedule!='B0' :: range(diff_vs_B0_pp)",
                  range(d$diff_vs_B0_pp), sprintf("%s~%s%%p", fnum(min(d$diff_vs_B0_pp), 1), sprintf("%+.1f", round(max(d$diff_vs_B0_pp) + 1e-9, 1))))
  premise(min(d$diff_vs_B0_pp) < 0 && max(d$diff_vs_B0_pp) > 0, "the change has both signs (title)")
  premise(all(rows(TCH, "set=='iii'")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUCinf than retained ones, both models (title)")
  f <- list(rng = rng)
  y0 <- core_title(tx("B1.title", f), tx("B1.kicker"))
  L <- DK$txt$B1$fig; ML <- DK$txt$common$models_short

  # ---- 왼쪽 위: 개념도(Day 29~57) ----
  sched <- .read("config/trial_design.yaml")$schedules
  b0 <- unlist(sched$B0$days) + 1; d3 <- unlist(sched$D3$days) + 1
  late <- function(x) x[x >= 29 - 1e-9]
  vb <- data.table(day = late(b0), y = 1); va <- data.table(day = setdiff(late(d3), late(b0)), y = 1)
  subj <- rows(RS, "model=='k2020' & percentile==50"); premise(nrow(subj) == 1, "one 2020 model cliff representative subject (S4)"); cx0 <- subj$start1 + 1; cx1 <- subj$t_lloq + 1
  i7 <- diff(late(b0)); premise(length(unique(i7)) == 1, "equal late B0 interval"); iv <- unique(i7)
  gmin <- min(diff(sort(late(d3)))); premise(gmin == 3, "smallest late interval with D3 is 3 days (concept)")
  k36 <- vb$day[vb$day < cx0 & vb$day + iv > cx1]; premise(length(k36) == 1, "the representative cliff lies inside one B0 interval")
  ki <- k36 + iv                                                     # 간격 화살표는 절벽 다음 간격에(글자가 겹치지 않게)
  ka <- max(va$day); kb <- max(vb$day[vb$day < ka]); premise(ka - kb == gmin || any(diff(sort(late(d3))) == gmin), "an added visit gives the smallest interval")
  kb2 <- sort(late(d3))[max(which(diff(sort(late(d3))) == gmin))]; ka2 <- kb2 + gmin
  pa <- ggplot() +
    annotate("rect", xmin = cx0, xmax = cx1, ymin = 0.6, ymax = 1.4, fill = PAL$orange, alpha = 0.35) +
    annotate("segment", x = min(vb$day) - 0.5, xend = max(vb$day) + 0.5, y = 1, yend = 1, colour = PAL$grid, linewidth = 1.2) +
    geom_point(data = vb, aes(day, y), shape = 21, fill = PAL$ink, colour = "white", size = 4.2, stroke = 0.6) +
    geom_point(data = va, aes(day, y), shape = 21, fill = "white", colour = PAL$ink, size = 3.8, stroke = 1.1) +
    annotate("text", x = (cx0 + cx1) / 2, y = 1.75, label = fill(L$cliff, list(len = fnum(subj$len1, 1))), size = PT(15), family = FONT, colour = PAL$orange, fontface = "bold") +
    annotate("segment", x = ki, xend = ki + iv, y = 1.75, yend = 1.75, colour = PAL$ink2, linewidth = 0.6, arrow = grid::arrow(ends = "both", length = grid::unit(0.08, "in"))) +
    annotate("text", x = ki + iv / 2, y = 2.2, label = fill(L$interval, list(d = fnum(iv, 0))), size = PT(15), family = FONT, colour = PAL$ink2) +
    annotate("segment", x = kb2, xend = ka2, y = 0.45, yend = 0.45, colour = PAL$ink2, linewidth = 0.6, arrow = grid::arrow(ends = "both", length = grid::unit(0.08, "in"))) +
    annotate("text", x = (kb2 + ka2) / 2, y = 0.2, vjust = 1, label = fill(L$added, list(d = fnum(gmin, 0))), size = PT(15), family = FONT, colour = PAL$ink2) +
    scale_x_continuous(breaks = vb$day, limits = c(min(vb$day) - 0.8, max(vb$day) + 0.8), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(-0.35, 2.4), breaks = NULL) +
    labs(x = L$xlab, y = NULL, subtitle = fill(L$concept, list(d = fnum(min(vb$day), 0)))) + theme_core(16) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  # ---- 왼쪽 아래: 세트 (iii) 충족 비율 변화 ----
  d[, sch := factor(schedule, levels = c("D4", "D3", "D2", "D1"))][, mod := factor(unlist(ML[model]), levels = unlist(ML))]
  pb <- ggplot(d, aes(diff_vs_B0_pp, sch, shape = mod, colour = mod)) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_point(size = 3.4) +
    scale_shape_manual(values = unname(CORE_MODEL_SHAPE), name = NULL) + scale_colour_manual(values = unname(CORE_MODEL_COL), name = NULL) +   # 범례: 주 모델(2020) 먼저
    scale_x_continuous(labels = function(v) ifelse(abs(v) < 1e-9, "0", sprintf("%+g", v)), expand = expansion(add = 0.3)) +
    labs(x = L$dx, y = NULL, subtitle = L$change) + theme_core(16) + theme(legend.position = "right", legend.justification = c(0, 0.5))
  # ---- 오른쪽: 탈락자 대 유지자의 참 AUCinf 비 ----
  r <- copy(rows(TCH, "set=='iii'"))[, mod := factor(unlist(ML[pk_model]), levels = rev(unlist(ML)))]
  r[, lab := sprintf("%s (%s~%s)", fnum(true_aucinf_gmr, 2), fnum(true_aucinf_gmr_lo, 2), fnum(true_aucinf_gmr_hi, 2))]
  pc <- ggplot(r, aes(true_aucinf_gmr, mod, shape = mod)) +
    geom_vline(xintercept = 1, colour = PAL$ink2, linetype = "22", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = true_aucinf_gmr_lo, xmax = true_aucinf_gmr_hi), height = 0.18, colour = PAL$orange, linewidth = 0.9) +
    geom_point(size = 4.2, colour = PAL$orange) + scale_shape_manual(values = setNames(unname(CORE_MODEL_SHAPE), unlist(ML)), guide = "none") +
    geom_text(aes(label = lab), vjust = -1.2, size = PT(16), family = FONT, colour = PAL$ink) +
    annotate("text", x = 1, y = 0.45, label = L$equal, hjust = 1.15, size = PT(14), family = FONT, colour = PAL$ink2) +
    scale_x_continuous(limits = c(0.8, 1.04), breaks = c(0.8, 0.9, 1.0)) + scale_y_discrete(expand = expansion(add = c(0.7, 0.7))) +
    labs(x = L$rx, y = NULL, subtitle = L$ratio) + theme_core(16)
  p <- patchwork::wrap_plots(patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1.15, 1)), pc, widths = c(1.35, 1))

  ge2 <- { w <- "weight=='base' & definition_day==1 & timing=='windowed' & grepl('^plus_', schedule)"; x <- max(rows(CP, w)$pct_ge2)
    dderived("share with two or more samples on the cliff, candidate schedules with visit windows, largest over schedules and models", CP, sprintf("%s :: max(pct_ge2)", w), x, paste0(fnum(cl(x, 2), 2), "%")) }
  body <- tx("B1.body", list(min = fnum(gmin, 0) |> (\(x) dderived("smallest late sampling interval with candidate schedule D3 (days)", "config/trial_design.yaml", "schedules.D3.days :: min(diff(days >= 28))", gmin, x))(),
                             len = drange(CS, "weight=='base'", "len1_median", 1, "", "median cliff length, two models"),
                             ge2 = ge2,
                             gmr = drange(TCH, "set=='iii'", "true_aucinf_gmr", 2, "", "true AUCinf ratio, failing to retained subjects, set (iii), two models")))
  by <- core_body(body, GEO$BODY_BOTTOM)
  deck_figure(p, "b1_schedule_and_bias", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(RL, TCH, RS))
  deck_notes(tx("B1.notes", list(ge2 = ge2, iv = dderived("late B0 sampling interval (days)", "config/trial_design.yaml", "schedules.B0.days :: diff(days >= 28)", iv, fnum(iv, 0)),
    b16 = dv(RL, "model=='k2016' & schedule=='B0'", "pct_reliable_iii", 1, "%", "share meeting set (iii), B0, 2016"), b20 = dv(RL, "model=='k2020' & schedule=='B0'", "pct_reliable_iii", 1, "%", "share meeting set (iii), B0, 2020"),
    d3 = drange(RL, "schedule=='D3'", "diff_vs_B0_pp", 1, "%p", "change with D3, two models"),
    w = drange(TCH, "set=='iii'", "wt_diff_kg", 1, "", "weight difference failing minus retained (kg), set (iii), two models"),
    al = drange(TCH, "set=='iii'", "auclast_gmr", 2, "", "AUClast ratio failing to retained, set (iii), two models"),
    gmin = dderived("smallest late sampling interval with candidate schedule D3 (days)", "config/trial_design.yaml", "schedules.D3.days :: min(diff(days >= 28))", gmin, fnum(gmin, 0)),
    nsub = dint(RL, "model=='k2016' & schedule=='B0'", "n", "virtual subjects per model"),
    ad = drange("trialpop/tp_arm_difference.csv", "set=='iii' & grepl('_(080|125)$', scenario) & !grepl('^ka_', scenario)", "diff_mean", 1, "%p",
                "arm difference in the set (iii) failing share (test - reference), boundary cells except absorption rate, two models"))))
  deck_end()
}
