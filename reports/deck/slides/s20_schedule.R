# S20 제안의 실행 · 채혈 일정 판단: 현행 B0 유지. 사전 규정 결정 규칙(보고서 3.9절, D-026)의 기준 (a)~(d)를 후보 D1~D4와 B-(Day 50 제거)에
# 두 모델로 적용한 결과표(보고서 Table 5-14와 같은 열)와, Day 57 이후 채혈의 가치(참 농도가 LLOQ 위에 남는 대상자 비율).
# 시험 모집단(건강인 60~90 kg, 두 구조 모델 base·struct2020)만. 체중 50~90 kg 분포(weight_alt)는 시험 모집단 밖이라 쓰지 않는다(부록 A5).
# 자료: results/trials/schedule_decision_<base|struct2020>.csv(대상자 20,000명, 시험 500회), results/reliability/reliability_paired_vs_B0.csv(세트 (i)),
#       results/individual200k/criterion_d_200k_<variant>.csv(200,000명 재평가), results/cliff/cliff_summary.csv(LLOQ 도달일, Day 58 이후 비율),
#       results/deck_inputs/cliff_lloq_day_hist.csv(LLOQ 도달 연구일 1일 간격 분포; 커밋된 요약과 대조: deck_inputs/checks.csv),
#       results/individual/individual_<base|struct2020>.csv(마지막 채혈 시료 정량 비율).
# 주의: 기준 (a)·(b)는 합동 t 검정(M0, trial_design.yaml be.method)으로만 계산했다. M1 일정 판정 파일은 없다(미산출).
#       Day 57 이후 채혈 일정은 어떤 결과 파일에서도 모의하지 않았다(Day 85 연장안은 DECISIONS D-016에서 철회). 그래서 "Day 85 불필요"를 직접 보이는 결과는 없고,
#       Day 58 이후에도 참 농도가 LLOQ 위인 대상자 비율로만 말한다.
S20_V <- c(k2016 = "base", k2020 = "struct2020")
S20_SCHED <- c("D1", "D2", "D3", "D4", "Bminus")
s20_sd <- function(v) sprintf("trials/schedule_decision_%s.csv", v)
s20_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?%?$", p)) p else paste0("+", p)
# 표시층: 음수 부호에서 줄이 바뀌지 않도록 줄바꿈 없는 붙임표(U+2011; check_deck이 숫자 검사 전에 "-"로 되돌린다)
s20_nb <- function(p) gsub("(^|~|\\s)-(?=[0-9])", "\\1\u2011", p, perl = TRUE)
# config 채혈일(투여 후 일) -> 연구일, 일정 간 차이(추가 또는 제거된 연구일)
s20_days <- function(s) { d <- unlist(.read("config/trial_design.yaml")$schedules[[s]]$days); premise(length(d) > 0, paste("schedule", s)); round(d + 1, 6) }
s20_diff <- function(s, what = c("added", "removed")) {
  what <- match.arg(what); a <- if (what == "added") setdiff(s20_days(s), s20_days("B0")) else setdiff(s20_days("B0"), s20_days(s))
  a <- sort(a); premise(length(a) >= 1 && all(a == round(a)), paste("whole study days", what, s))
  dderived(sprintf("schedule %s: study days %s relative to B0", s, what), "config/trial_design.yaml",
           sprintf("%s :: schedules.%s.days + 1", if (what == "added") sprintf("setdiff(%s, B0)", s) else sprintf("setdiff(B0, %s)", s), s), a, paste(fnum(a, 0), collapse = ", "))
}
# 두 모델 값 "2016 / 2020"
s20_pair <- function(s, col, d, item, scale = 1, signed = FALSE) {
  # 행 조건은 "schedule == 'X'"(공백 포함): 보고서 bm_d가 같은 파일의 "schedule=='Bminus' :: d_extrap20_ratio"에 여러 변이의 범위를 기록해 두어 대조(검사 8)가 엇갈리지 않게 한다
  v <- vapply(names(S20_V), function(m) { p <- dv(s20_sd(S20_V[[m]]), sprintf("schedule == '%s'", s), col, d, "", sprintf("%s, schedule %s, %s", item, s, m), scale)
    if (signed) s20_signed(p) else p }, "")
  s20_nb(paste(v, collapse = " / "))
}
# 두 결과 파일에 걸친 범위(부호 붙임, 파생값; 두 파일을 locator에 적는다)
s20_rng2 <- function(where, col, d, item, scale = 1) {
  x <- unlist(lapply(S20_V, function(v) rows(s20_sd(v), where)[[col]] * scale))
  dsrc(item, s20_sd(S20_V[["k2020"]]), "(table)")
  p <- paste0(s20_signed(fnum(min(x), d)), "~", s20_signed(fnum(max(x), d)))
  s20_nb(dderived(item, s20_sd(S20_V[["k2016"]]), sprintf("range over trials/schedule_decision_base.csv and trials/schedule_decision_struct2020.csv [%s] :: %s x %s", where, col, scale), range(x), p))
}
# 두 결과 파일에 걸친 최댓값(파생값; 두 파일을 locator에 적는다)
s20_max2 <- function(where, col, d, item, scale = 1) {
  x <- vapply(S20_V, function(v) max(rows(s20_sd(v), where)[[col]] * scale), 0)
  dsrc(item, s20_sd(S20_V[["k2020"]]), "(table)")
  dderived(item, s20_sd(S20_V[["k2016"]]), sprintf("max over trials/schedule_decision_base.csv and trials/schedule_decision_struct2020.csv [%s] :: %s x %s", where, col, scale), max(x), fnum(max(x), d))
}

slide_S20 <- function() {
  RP <- "reliability/reliability_paired_vs_B0.csv"; CS <- "cliff/cliff_summary.csv"; CH <- "deck_inputs/cliff_lloq_day_hist.csv"
  D200 <- function(v) sprintf("individual200k/criterion_d_200k_%s.csv", v); IND <- function(v) sprintf("individual/individual_%s.csv", v)
  RW <- "variant %in% c('base','struct2020') & schedule %in% c('D1','D2','D3','D4')"
  deck_slide("S20", tag = "sim")
  L <- DK$txt$S20

  # ---- 전제: 문장이 기대는 사실 ----
  for (v in S20_V) {
    r <- rows(s20_sd(v)); premise(setequal(r$schedule, S20_SCHED), paste("candidate schedules D1 to D4 and Bminus in", v))
    premise(!any(r$crit_a %in% TRUE | r$crit_b %in% TRUE | r$crit_c %in% TRUE | r$crit_d %in% TRUE | r$recommend %in% TRUE), paste("no criterion met and no recommendation,", v))
    premise(all(r[schedule == "Bminus", a_mean_width_rel_decrease] < .read("config/trial_design.yaml")$decision_rule$ci_width_mean_rel_decrease_min) && all(r[schedule == "Bminus", d_extrap20_ratio] > 1),
            paste("removing Day 50: width decrease below criterion (a) and more subjects above 20% extrapolation,", v))
  }
  premise(identical(rows(s20_sd("base"))[order(schedule), added_visits_total], rows(s20_sd("struct2020"))[order(schedule), added_visits_total]), "added visits equal in both models")
  premise(.read("config/schedule_decision.yaml")$final_schedule == "B0", "final schedule decision is B0")
  premise(.read("config/trial_design.yaml")$be$method == "pooled_t", "criteria (a) and (b) come from trials analysed with the pooled t-test (M0)")
  sch <- .read("config/trial_design.yaml")$schedules
  premise(max(unlist(lapply(sch, function(z) unlist(z$days)))) + 1 == max(s20_days("B0")), "no configured schedule samples after the last B0 study day (Day 57)")
  rp <- rows(RP, RW); premise(nrow(rp) == 8 && !any(rp$crit_c_i %in% TRUE) && !any(rp$crit_c_ii %in% TRUE), "criterion (c) met under neither criteria set")
  premise(rp[which.max(gain_i_pp), variant == "base" & schedule == "D3"], "largest set (i) gain: 2016 model, D3")
  for (v in S20_V) premise(!rows(D200(v))$d_meets_point, paste("criterion (d) not met at 200,000 subjects,", v))
  h <- rows(CH); cs <- rows(CS, "weight=='base'")
  for (m in names(S20_V)) premise(abs(sum(h[pk_model == m & day_bin >= 58, pct]) - cs[model == m, lloq_after_day58_pct]) < 1e-6 && abs(sum(h[pk_model == m, pct]) - 100) < 1e-6,
                                  paste("histogram tail from Day 58 equals the committed share above the LLOQ after Day 58,", m))

  # ---- 제목 ----
  last <- f_study_days("B0", "last")
  f <- list(n = f_study_days("B0", "n"), last = last, k = dcount(s20_sd("base"), "TRUE", "candidate schedules (added or removed samples)"))
  deck_kicker(tx("S20.kicker")); deck_title(tx("S20.title", f))

  # ---- 표: 후보 일정별 기준 (a)~(d), 두 모델 ----
  thr <- list(a = dcfg("trial_design.yaml", c("decision_rule", "ci_width_mean_rel_decrease_min"), "criterion (a) threshold (%)", function(x) fnum(100 * x, 0)),
              b = dcfg("trial_design.yaml", c("decision_rule", "ke110_pass_gain_pp_min"), "criterion (b) threshold (points)", num_fmt(0)),
              c = dcfg("trial_design.yaml", c("decision_rule", "reliability_gain_pp_min"), "criterion (c) threshold (points)", num_fmt(0)),
              d = dcfg("trial_design.yaml", c("decision_rule", "extrap_gt20_ratio_max"), "criterion (d) threshold (ratio to B0)", num_fmt(1)),
              x = dcfg("nca_rules.yaml", c("standard", "reliability", "extrap_max_pct"), "NCA extrapolation flag threshold (%)", num_fmt(0)),
              ci = f_ci_level(), ke = dcfg("scenarios.yaml", c("scenarios", "KE110", "T_multipliers", "ke"), "ke multiplier of the scenario for criterion (b)", num_fmt(2)),
              cost = dcfg("trial_design.yaml", c("decision_rule", "cost_per_added_point", "subjects"), "added visits per added sample point", function(x) fnum(as.numeric(x), 0, big = TRUE)))
  slab <- function(s) if (s == "Bminus") fill(L$table$removed, list(d = s20_diff(s, "removed"))) else fill(L$table$added, list(s = s, d = s20_diff(s, "added")))
  df <- data.frame(a = vapply(S20_SCHED, slab, ""),
                   v = vapply(S20_SCHED, function(s) s20_nb(dint(s20_sd("base"), sprintf("schedule=='%s'", s), "added_visits_total", sprintf("added visits, %s", s))), ""),
                   ca = vapply(S20_SCHED, s20_pair, "", col = "a_mean_width_rel_decrease", d = 2, item = "(a) relative decrease of mean AUC0-last CI width (%)", scale = 100, signed = TRUE),
                   cb = vapply(S20_SCHED, s20_pair, "", col = "b_ke110_gain_pp", d = 1, item = "(b) AUC0-last pass change at ke x1.10 (points)", signed = TRUE),
                   cc = vapply(S20_SCHED, s20_pair, "", col = "c_reliable_gain_pp", d = 2, item = "(c) reliability change, set (ii) (points)", signed = TRUE),
                   cd = vapply(S20_SCHED, s20_pair, "", col = "d_extrap20_ratio", d = 2, item = "(d) share with NCA extrapolation above 20%, ratio to B0"),
                   met = rep(L$table$none, length(S20_SCHED)), stringsAsFactors = FALSE, check.names = FALSE)
  # 첫 행: 사전 기준(머리글과 같은 색). 기준값은 config에서
  crit <- vapply(unlist(L$table$crit), fill, "", facts = thr, USE.NAMES = FALSE)
  df <- rbind(setNames(as.data.frame(as.list(crit), stringsAsFactors = FALSE), names(df)), df)
  names(df) <- tx("S20.table.head", thr)
  TH <- 2.46   # 렌더링 실측(빈 칸이 있는 기준 행은 조금 높다)
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, GEO$CW, TH), widths = c(2.4, 1.12, 2.0, 2.0, 2.0, 2.0, 0.71), size = 12, highlight = 1, highlight_fill = PAL$tint_blue, label = "table_schedule")
  FY <- GEO$BODY_TOP + TH + 0.02; FTH <- 0.4
  deck_text(tx("S20.foot"), c(GEO$ML, FY, GEO$CW, FTH), size = 16, color = PAL$ink2, label = "text_table_note", gap_pt = 0)

  # ---- 왼쪽 아래: 요점 ----
  BY <- FY + FTH + 0.02; LW <- 7.4
  ab <- lapply(names(S20_V), function(m) dv(CS, sprintf("model=='%s' & weight=='base'", m), "lloq_after_day58_pct", 1, "%", sprintf("share above the LLOQ after study Day 58, %s", m)))
  ADD <- "schedule %in% c('D1','D2','D3','D4')"
  b <- list(arng = s20_rng2(ADD, "a_mean_width_rel_decrease", 2, "relative decrease of the mean AUC0-last CI width with added samples (%), range over D1 to D4, two models", scale = 100),
            amax = s20_signed(s20_max2(ADD, "a_mean_width_rel_decrease", 2, "largest relative decrease of the mean AUC0-last CI width with added samples (%), two models", scale = 100)),
            ci = s20_signed(dext(RP, RW, "gain_i_pp", max, 2, "", "largest reliability gain versus B0, set (i), trial population, D1 to D4")),
            cii = s20_signed(dext(RP, RW, "gain_ii_pp", max, 2, "", "largest reliability gain versus B0, set (ii), trial population, D1 to D4")),
            c = thr$c, a = thr$a, x = thr$x,
            bma = s20_pair("Bminus", "a_mean_width_rel_decrease", 2, "(a) CI width decrease when Day 50 is removed (%)", scale = 100),
            bmd = s20_pair("Bminus", "d_extrap20_ratio", 2, "(d) ratio of the share above 20% extrapolation when Day 50 is removed"),
            bmc = s20_pair("Bminus", "c_reliable_gain_pp", 2, "(c) reliability change, set (ii), when Day 50 is removed (points)"),
            d50 = s20_diff("Bminus", "removed"), last = last,
            d58 = dderived("study-day threshold in column name lloq_after_day58_pct", CS, "column name lloq_after_day58_pct (100 x mean(t_lloq + 1 > 58))", 58, "58"),
            a16 = ab[[1]], a20 = ab[[2]])
  premise(as.numeric(b$d58) == as.numeric(last) + 1, "Day 58 is the day after the last B0 sample")
  premise(rp[which.max(gain_ii_pp), variant == "base" & schedule == "D3"] && rp[which.max(gain_i_pp), variant == "base" & schedule == "D3"], "largest set (i) and set (ii) gains both at D3 (2016 model)")
  premise(grepl("AUC0-inf", .read("config/schedule_decision.yaml")$day50_removal) && grepl("fallback", .read("config/schedule_decision.yaml")$day50_removal), "recorded reason for keeping Day 50: AUC0-inf secondary endpoint and fallback terminal points")
  for (v in S20_V) premise(all(rows(s20_sd(v), "schedule=='Bminus'")$c_reliable_gain_pp < 0), paste("removing Day 50 lowers reliability under set (ii),", v))
  deck_bullets(tx("S20.bullets", b), box = c(GEO$ML, BY, LW, GEO$BODY_BOTTOM - BY), size = 16, gap_pt = 5)

  # ---- 오른쪽 아래: 참 농도가 LLOQ 위인 대상자 비율(연구일별), B0 채혈일, Day 57 이후(확대 그림 포함) ----
  XR <- GEO$ML + LW + 0.2; WR <- GEO$W - GEO$MR - XR
  dmax <- max(h$day_bin) + 1; grid_ <- seq(min(h$day_bin), dmax)
  sv <- rbindlist(lapply(names(S20_V), function(m) { hh <- h[pk_model == m]; data.table(model = m, day = grid_, pct = vapply(grid_, function(d) sum(hh[day_bin >= d, pct]), 0)) }))
  sv[, mod := factor(model_lab()[model], levels = model_lab())]
  b0 <- s20_days("B0"); b0t <- b0[b0 >= min(grid_) & b0 == round(b0)]; lastn <- max(b0)
  tail_ <- sv[day == lastn + 1][, lab := paste0(fnum(pct, 1), "%")]
  premise(identical(tail_[order(model), lab], vapply(ab, identity, "")), "figure labels equal the committed shares above the LLOQ after Day 58")
  FL <- L$fig
  FH <- GEO$BODY_BOTTOM - BY
  scm <- list(scale_colour_manual(values = unname(MODEL_COL)), scale_linetype_manual(values = unname(MODEL_LT)), scale_shape_manual(values = unname(MODEL_SHAPE)))
  # 확대 그림: Day 57 이후, 세로축 0~(최댓값 x 1.3)
  si <- sv[day >= lastn]; yin <- max(si$pct) * 1.25
  pin <- ggplot(si, aes(day, pct, colour = mod, linetype = mod)) + geom_step(linewidth = 0.8, direction = "hv") +
    geom_point(data = tail_, aes(day, pct, colour = mod, shape = mod), size = 1.8, inherit.aes = FALSE) +
    geom_label(data = tail_, aes(day + 1.3, pct, label = lab, colour = mod), hjust = 0, vjust = 0.5, size = 3.5, family = FONT, fontface = "bold", fill = "white",
               label.size = 0, label.padding = grid::unit(0.08, "lines"), inherit.aes = FALSE) +
    scm + scale_x_continuous(breaks = seq(lastn, dmax, by = 14), expand = expansion(add = c(0.3, 0.3))) +
    scale_y_continuous(limits = c(0, yin), breaks = scales::breaks_pretty(3), expand = expansion(add = c(0.1, 0))) +
    labs(x = NULL, y = NULL, title = fill(FL$inset, list(d = fnum(lastn, 0), d1 = fnum(lastn + 1, 0)))) + theme_deck(10) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.border = element_rect(colour = PAL$muted, fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 9.5, colour = PAL$ink2, margin = margin(0, 0, 2, 0)), plot.title.position = "plot",
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(2, 4, 2, 2))
  p <- ggplot(sv, aes(day, pct, colour = mod, linetype = mod)) +
    annotate("rect", xmin = lastn, xmax = Inf, ymin = -Inf, ymax = Inf, fill = PAL$tint_orange, alpha = 0.8) +
    geom_vline(xintercept = setdiff(b0t, lastn), colour = PAL$grid, linewidth = 0.5) +
    geom_vline(xintercept = lastn, colour = PAL$ink2, linewidth = 0.6, linetype = "22") +
    geom_step(linewidth = 1.0, direction = "hv") +
    annotate("text", x = lastn + 0.8, y = 99, hjust = 0, vjust = 1, size = 3.5, family = FONT, colour = PAL$ink2, label = fill(FL$last, list(d = fnum(lastn, 0)))) +
    annotation_custom(ggplotGrob(pin), xmin = lastn + 1, xmax = dmax + 0.5, ymin = 6, ymax = 89) +
    scm + scale_x_continuous(breaks = b0t, expand = expansion(add = c(0.5, 0.5))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(add = c(1, 1))) +
    labs(x = FL$xlab, y = NULL, subtitle = FL$ylab) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.major.x = element_blank(), plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))
  deck_figure(p, "s20_lloq_tail", c(XR, BY, WR, FH), src = c(CH, CS))

  # ---- 노트 ----
  d2 <- function(v, item) dci(D200(v), "TRUE", "extrap_gt20_ratio", "d_ratio_boot_lo", "d_ratio_boot_hi", 3, "", item)
  vm <- "trials/schedule_decision_vmax080_both.csv"
  nts <- list(
    a = thr$a, b = thr$b, c = thr$c, d = thr$d, x = thr$x, ci = thr$ci, ke = thr$ke,
    cost = dcfg("trial_design.yaml", c("decision_rule", "cost_per_added_point", "subjects"), "added visits per added sample point", num_fmt(0)),
    nsub = dcfg("trial_design.yaml", c("mc", "n_individual"), "subjects per model and variant", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    ntr = dcfg("trial_design.yaml", c("mc", "n_trials"), "trials per model and variant", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    date = dcfg("schedule_decision.yaml", "decision_date", "schedule decision date", function(x) as.character(x)),
    d3i16 = dci(RP, "variant=='base' & schedule=='D3'", "gain_i_pp", "gain_i_lo", "gain_i_hi", 2, "%p", "reliability change set (i), D3, k2016 (points)"),
    d3i20 = dv(RP, "variant=='struct2020' & schedule=='D3'", "gain_i_pp", 2, "%p", "reliability change set (i), D3, k2020 (points)"),
    irng = drange(RP, RW, "gain_i_pp", 2, "", "reliability change set (i), D1 to D4, two models (points)"),
    iirng = drange(RP, RW, "gain_ii_pp", 2, "", "reliability change set (ii), D1 to D4, two models (points)"),
    dd16 = s20_pair("D3", "d_extrap20_ratio", 2, "(d) D3 ratio"), dboot = s20_nb(paste(vapply(S20_V, function(v) dci(s20_sd(v), "schedule=='D3'", "d_extrap20_ratio", "d_ratio_boot_lo", "d_ratio_boot_hi", 2, "", sprintf("(d) D3 ratio with bootstrap interval, %s", v)), ""), collapse = " / ")),
    n200 = dint(D200("base"), "TRUE", "n_subjects", "subjects in the criterion (d) re-evaluation"),
    e16 = d2("base", "criterion (d) at 200,000 subjects, D3, k2016"), e20 = d2("struct2020", "criterion (d) at 200,000 subjects, D3, k2020"),
    vmx = dcfg("scenarios.yaml", c("sensitivity_variants", "vmax080_both", "theta_multipliers", "Vmax"), "Vmax multiplier of the curvature variant (both arms)", num_fmt(1)),
    vmd = dci(vm, "schedule=='D3'", "d_extrap20_ratio", "d_ratio_boot_lo", "d_ratio_boot_hi", 3, "", "criterion (d), Vmax x0.8 both arms, D3, 20,000 subjects"),
    vm200 = d2("vmax080_both", "criterion (d) at 200,000 subjects, Vmax x0.8 both arms, D3"),
    vmabs = dv(D200("vmax080_both"), "TRUE", "d_abs_change_per_arm", 2, "", "absolute change per arm at 200,000 subjects, Vmax x0.8, D3"),
    vmv = dint(vm, "schedule=='D3'", "added_visits_total", "added visits, D3"),
    m16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, k2016"),
    m20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_median", 1, "", "median study day at LLOQ, k2020"),
    p16 = dv(CS, "model=='k2016' & weight=='base'", "lloq_studyday_p95", 1, "", "95th percentile study day at LLOQ, k2016"),
    p20 = dv(CS, "model=='k2020' & weight=='base'", "lloq_studyday_p95", 1, "", "95th percentile study day at LLOQ, k2020"),
    q16 = dv(IND("base"), "schedule=='B0'", "quant_at_last_pct", 1, "%", "quantifiable at the last planned sample, k2016"),
    q20 = dv(IND("struct2020"), "schedule=='B0'", "quant_at_last_pct", 1, "%", "quantifiable at the last planned sample, k2020"),
    wr = dcfg("trial_design.yaml", c("weight", "sensitivity", "trunc"), "alternative weight distribution, truncation (kg), outside the trial population", function(x) rng_fmt(x[1], x[2], 0)),
    a16 = b$a16, a20 = b$a20, d58 = b$d58, last = last, arng = b$arng, amax = b$amax, bma = b$bma, bmd = b$bmd, bmc = b$bmc, d50 = b$d50)
  premise(nrow(rows(vm, "recommend==TRUE")) == 1 && rows(vm, "recommend==TRUE")$schedule == "D3" && isTRUE(rows(vm, "recommend==TRUE")$crit_d), "Vmax x0.8 (both arms): D3 recommended on criterion (d) alone at 20,000 subjects (notes)")
  deck_notes(tx("S20.notes", nts))
  deck_end()
}
