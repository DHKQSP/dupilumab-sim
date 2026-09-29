# 별첨 A7 채혈 일정(지시 2026-09-29): 현행 B0(Syneos 기본 일정) 유지. 결과보고 덱 S20(slides/s20_schedule.R)의 자료 논리를 따른다.
#  표: 후보 일정(D1~D4 = 후반 채혈 추가, B- = Day 50 제거) x 결정 기준 (a)~(d), 두 모델(trials/schedule_decision_<base|struct2020>.csv;
#      기준값 config/trial_design.yaml decision_rule; 최종 판정 config/schedule_decision.yaml). 첫 본문 행 = 추가 권고 기준.
#      마지막 열(참고, 기준 아님): 신뢰할 수 있는 AUCinf(세트 (iii)) 비율의 B0 대비 변화(results/core_deck/reliable_iii_by_schedule.csv; 같은 대상자, B-는 미산출).
#  결정 규칙의 시점: 첫 일정 모의 뒤, 결과 보고 전(regulatory/tables/prespecification_register.csv). 그래서 '사전 규정'이라 부르지 않는다.
#  (a)·(b)는 일정 판정 모의의 합동 t 검정(M0) 결과다(M1 재계산 없음). (c)·(d)와 세트 (iii)은 대상자 수준이라 분석 모형과 무관하다.
#  시험 모집단(건강인, 두 구조 모델 base·struct2020)만. 채혈일: 연구일 = 투여 후 일 + 1.
A7_V <- c(k2016 = "base", k2020 = "struct2020")
A7_SCHED <- c("D1", "D2", "D3", "D4", "Bminus")
a7_sd <- function(v) sprintf("trials/schedule_decision_%s.csv", v)
a7_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?%?$", p)) p else paste0("+", p)
# config 채혈일(투여 후 일) -> 연구일, B0와의 차이(추가 또는 제거된 연구일)
a7_days <- function(s) { d <- unlist(.read("config/trial_design.yaml")$schedules[[s]]$days); premise(length(d) > 0, paste("schedule", s)); round(d + 1, 6) }
a7_diff <- function(s, what = c("added", "removed")) {
  what <- match.arg(what); a <- if (what == "added") setdiff(a7_days(s), a7_days("B0")) else setdiff(a7_days("B0"), a7_days(s))
  a <- sort(a); premise(length(a) >= 1 && all(a == round(a)), paste("whole study days", what, s))
  dderived(sprintf("schedule %s: study days %s relative to B0", s, what), "config/trial_design.yaml",
           sprintf("%s :: schedules.%s.days + 1", if (what == "added") sprintf("setdiff(%s, B0)", s) else sprintf("setdiff(B0, %s)", s), s), a, paste(fnum(a, 0), collapse = ", "))
}
# 두 모델 값 "2016 / 2020"(행 조건은 결과보고 덱 S20과 같은 "schedule == 'X'")
a7_pair <- function(s, col, d, item, scale = 1, signed = FALSE) {
  v <- vapply(names(A7_V), function(m) { p <- dv(a7_sd(A7_V[[m]]), sprintf("schedule == '%s'", s), col, d, "", sprintf("%s, schedule %s, %s", item, s, m), scale)
    if (signed) a7_signed(p) else p }, "")
  paste(v, collapse = " / ")
}
a7_same <- function(s, col, d, item, scale = 1) {
  v <- vapply(names(A7_V), function(m) dv(a7_sd(A7_V[[m]]), sprintf("schedule == '%s'", s), col, d, "", sprintf("%s, schedule %s, %s", item, s, m), scale), "")
  premise(length(unique(v)) == 1, paste("same printed value in both models:", item)); v[[1]]
}
# 두 결과 파일에 걸친 최댓값(파생값; 두 파일을 locator에 적는다)
a7_max2 <- function(where, col, d, item, scale = 1) {
  x <- vapply(A7_V, function(v) max(rows(a7_sd(v), where)[[col]] * scale), 0)
  dsrc(item, a7_sd(A7_V[["k2020"]]), "(table)")
  dderived(item, a7_sd(A7_V[["k2016"]]), sprintf("max over trials/schedule_decision_base.csv and trials/schedule_decision_struct2020.csv [%s] :: %s x %s", where, col, scale), max(x), fnum(max(x), d))
}


slide_A7 <- function() {
  RL <- "core_deck/reliable_iii_by_schedule.csv"; RP <- "reliability/reliability_paired_vs_B0.csv"; PSR <- "regulatory/tables/prespecification_register.csv"
  RW <- "variant %in% c('base','struct2020') & schedule %in% c('D1','D2','D3','D4')"
  deck_slide("A7", tag = "sim")
  L <- DK$txt$A7
  td <- .read("config/trial_design.yaml")

  # ---- 전제: 문장이 기대는 사실 ----
  for (v in A7_V) {
    r <- rows(a7_sd(v)); premise(setequal(r$schedule, A7_SCHED), paste("candidate schedules D1 to D4 and Bminus in", v))
    premise(!any(r$crit_a %in% TRUE | r$crit_b %in% TRUE | r$crit_c %in% TRUE | r$crit_d %in% TRUE | r$recommend %in% TRUE), paste("no criterion met and no recommendation,", v))
    premise(all(r$a_mean_width_rel_decrease < td$decision_rule$ci_width_mean_rel_decrease_min) && all(r$b_ke110_gain_pp < td$decision_rule$ke110_pass_gain_pp_min) &&
              all(r$c_reliable_gain_pp < td$decision_rule$reliability_gain_pp_min) && all(r$d_extrap20_ratio > td$decision_rule$extrap_gt20_ratio_max),
            paste("every value misses its threshold (table),", v))
    premise(all(r[schedule == "Bminus", d_extrap20_ratio] > 1) && all(r[schedule == "Bminus", c_reliable_gain_pp] < 0), paste("removing Day 50: more subjects above 20% extrapolation, lower set (ii) reliability,", v))
  }
  premise(identical(rows(a7_sd("base"))[order(schedule), added_visits_total], rows(a7_sd("struct2020"))[order(schedule), added_visits_total]), "added visits equal in both models")
  SDY <- .read("config/schedule_decision.yaml")
  premise(SDY$final_schedule == "B0", "final schedule decision is B0")
  # 민감도 변형: 규칙상 추가 채혈이 권고된 곳은 양 군 Vmax x0.8의 D3(기준 d 단독) 하나뿐이고, 판정 기록은 판단으로 B0를 유지했다(D-031, D-037)
  sdv <- sub("^schedule_decision_(.*)\\.csv$", "\\1", list.files(proj_path("results", "trials"), pattern = "^schedule_decision_.*\\.csv$"))
  rec <- rbindlist(lapply(sdv, function(v) { r <- rows(a7_sd(v), "recommend == TRUE"); if (nrow(r)) r[, .(v = v, schedule, crit_a, crit_b, crit_c, crit_d)] else NULL }))
  premise(nrow(rec) == 1 && rec$v == "vmax080_both" && rec$schedule == "D3" && rec$crit_d && !rec$crit_a && !rec$crit_b && !rec$crit_c, "the only rule recommendation over all variants: D3 in Vmax x0.8 (both arms), criterion (d) alone")
  VM <- a7_sd("vmax080_both"); vd <- row1(VM, "schedule == 'D3'"); v2 <- row1("individual200k/criterion_d_200k_vmax080_both.csv", "TRUE")
  premise(vd$d_ratio_boot_lo < td$decision_rule$extrap_gt20_ratio_max && vd$d_ratio_boot_hi > td$decision_rule$extrap_gt20_ratio_max, "Vmax x0.8 D3: bootstrap interval of criterion (d) includes the threshold")
  premise(abs(vd$d_abs_change_per_arm) < 1 && !v2$d_meets_point && v2$extrap_gt20_ratio > td$decision_rule$extrap_gt20_ratio_max, "Vmax x0.8 D3: absolute change below one subject per arm; not met at 200,000 subjects")
  nt_ <- SDY$sensitivity_d3_rule_triggers$note
  premise(grepl("0.5 포함", nt_, fixed = TRUE) && grepl("arm당 1명 미만", nt_, fixed = TRUE) && grepl(sprintf("추가 방문 %s회", format(vd$added_visits_total, big.mark = ",")), nt_, fixed = TRUE) && grepl("절대 하한", SDY$lesson, fixed = TRUE),
          "decision record: reasons for keeping B0 (interval includes 0.5, below one subject per arm, added visits, no absolute floor in the rule)")
  premise(grepl("AUC0-inf", SDY$day50_removal, fixed = TRUE) && grepl("fallback", SDY$day50_removal, fixed = TRUE) && grepl("말단 점", SDY$day50_removal, fixed = TRUE), "decision record: Day 50 kept for the AUCinf secondary endpoint and the fallback")
  premise(as.Date(SDY$decision_date) < as.Date(.read("config/prereg_20260926.yaml")$registered_on), "the schedule simulations and decision precede the v1.0.1 pre-registration (v1.0 results)")
  premise(td$be$method == "pooled_t", "criteria (a) and (b) come from trials analysed with the pooled t-test (M0)")
  premise(max(unlist(lapply(td$schedules, function(z) unlist(z$days)))) + 1 == max(a7_days("B0")), "no configured schedule samples after the last B0 study day")
  ps <- rows(PSR, "grepl('^Sampling-schedule decision rule', Item)")
  premise(nrow(ps) == 1 && grepl("after the first schedule simulations", ps$Status) && grepl("before any result was reported", ps$Status),
          "decision rule fixed after the first schedule simulations, before any result was reported (not pre-specified)")
  rl <- rows(RL); premise(setequal(rl$schedule, c("B0", "D1", "D2", "D3", "D4")) && all(rl[schedule == "B0", diff_vs_B0_pp] == 0), "set (iii) by schedule: B0 and D1-D4 only (Bminus not computed)")
  premise(all(abs(rl[schedule != "B0", diff_vs_B0_pp]) < td$decision_rule$reliability_gain_pp_min), "set (iii) changes are all smaller than the criterion (c) threshold (text)")
  premise(rl[model == "k2016"][which.max(diff_vs_B0_pp), schedule] == "D3" && rl[model == "k2020"][which.max(diff_vs_B0_pp), schedule] == "D3", "set (iii) rises only with D3 (notes)")
  premise(all(rl[schedule %in% c("D1", "D2", "D4"), diff_vs_B0_pp] < 0) && all(rl[schedule == "D3", diff_vs_B0_pp] > 0), "D1, D2, D4 lower and D3 raises set (iii) (notes)")
  aw <- rbindlist(lapply(names(A7_V), function(m) copy(rows(a7_sd(A7_V[[m]]), "schedule %in% c('D1','D2','D3','D4')"))[, pk := m]))
  premise(nrow(aw) == 8 && sum(aw$a_mean_width_rel_decrease > 0) == 1 && max(aw$a_mean_width_rel_decrease) > 0 && aw[a_mean_width_rel_decrease > 0, schedule == "D3" & pk == "k2020"],
          "added samples: the mean CI width narrows only in D3 of the 2020 model (notes)")

  # ---- 제목 ----
  f <- list(n = f_study_days("B0", "n"), last = f_study_days("B0", "last"))
  y0 <- core_title(tx("A7.title", f), tx("A7.kicker"))

  # ---- 표: 후보 일정별 기준 (a)~(d), 두 모델 + 참고 열(세트 (iii)) ----
  thr <- list(a = dcfg("trial_design.yaml", c("decision_rule", "ci_width_mean_rel_decrease_min"), "criterion (a) threshold (%)", function(x) fnum(100 * x, 0)),
              b = dcfg("trial_design.yaml", c("decision_rule", "ke110_pass_gain_pp_min"), "criterion (b) threshold (points)", num_fmt(0)),
              c = dcfg("trial_design.yaml", c("decision_rule", "reliability_gain_pp_min"), "criterion (c) threshold (points)", num_fmt(0)),
              d = dcfg("trial_design.yaml", c("decision_rule", "extrap_gt20_ratio_max"), "criterion (d) threshold (ratio to B0)", num_fmt(1)),
              x = dcfg("nca_rules.yaml", c("standard", "reliability", "extrap_max_pct"), "NCA extrapolation flag threshold (%)", num_fmt(0)),
              ci = f_ci_level(), ke = dcfg("scenarios.yaml", c("scenarios", "KE110", "T_multipliers", "ke"), "ke multiplier of the scenario for criterion (b)", num_fmt(2)))
  slab <- function(s) if (s == "Bminus") fill(L$table$removed, list(d = a7_diff(s, "removed"))) else fill(L$table$added, list(s = s, d = a7_diff(s, "added")))
  r3 <- function(s) if (s == "Bminus") L$table$na else paste(vapply(names(A7_V), function(m)
    a7_signed(dv(RL, sprintf("model=='%s' & schedule=='%s'", m, s), "diff_vs_B0_pp", 1, "", sprintf("change in the share meeting set (iii) vs B0, schedule %s, %s (points)", s, m))), ""), collapse = " / ")
  df <- data.frame(a = vapply(A7_SCHED, slab, ""),
                   ca = vapply(A7_SCHED, a7_pair, "", col = "a_mean_width_rel_decrease", d = 2, item = "(a) relative decrease of mean AUC0-last CI width (%)", scale = 100, signed = TRUE),
                   cb = vapply(A7_SCHED, a7_pair, "", col = "b_ke110_gain_pp", d = 1, item = "(b) AUC0-last pass change at ke x1.10 (points)", signed = TRUE),
                   cc = vapply(A7_SCHED, a7_pair, "", col = "c_reliable_gain_pp", d = 2, item = "(c) reliability change, set (ii) (points)", signed = TRUE),
                   cd = vapply(A7_SCHED, a7_pair, "", col = "d_extrap20_ratio", d = 2, item = "(d) share with NCA extrapolation above 20%, ratio to B0"),
                   c3 = vapply(A7_SCHED, r3, ""), stringsAsFactors = FALSE, check.names = FALSE)
  crit <- vapply(unlist(L$table$crit), fill, "", facts = thr, USE.NAMES = FALSE)
  df <- rbind(setNames(as.data.frame(as.list(crit), stringsAsFactors = FALSE), names(df)), df)
  names(df) <- tx("A7.table.head", thr)
  WD <- c(2.45, 1.95, 1.95, 1.9, 1.9, 2.05)
  TH <- deck_table_h(df, GEO$CW, WD, 14, pad = 2) + 0.02
  deck_table(df, box = c(GEO$ML, y0, GEO$CW, TH), widths = WD, size = 14, highlight = 1, highlight_fill = PAL$tint_blue, label = "table_schedule", pad = 2)   # 셀 위아래 여백 2pt: 표 아래 요점·캡션 자리

  # ---- 요점과 캡션 ----
  ADD <- "schedule %in% c('D1','D2','D3','D4')"
  b <- list(amax = (a7_max2(ADD, "a_mean_width_rel_decrease", 2, "largest relative decrease of the mean AUC0-last CI width with added samples (%), two models", scale = 100)),
            a = thr$a, c = thr$c, x = thr$x,
            cii = a7_signed(dext(RP, RW, "gain_ii_pp", max, 2, "", "largest reliability gain versus B0, set (ii), trial population, D1 to D4")),
            r3 = { d <- rl[schedule != "B0"]; dderived("change in the share meeting set (iii), D1-D4 vs B0, range over schedules and models (percentage points)", RL, "schedule!='B0' :: range(diff_vs_B0_pp)",
                   range(d$diff_vs_B0_pp), sprintf("%s~%s", fnum(min(d$diff_vs_B0_pp), 1), a7_signed(fnum(max(d$diff_vs_B0_pp), 1)))) },
            d50 = a7_diff("Bminus", "removed"),
            bma = a7_pair("Bminus", "a_mean_width_rel_decrease", 2, "(a) CI width decrease when Day 50 is removed (%)", scale = 100),
            bmd1 = a7_same("Bminus", "d_extrap20_ratio", 2, "(d) ratio of the share above 20% extrapolation when Day 50 is removed, same in both models"))
  vs <- list(vm = dcfg("scenarios.yaml", c("sensitivity_variants", "vmax080_both", "theta_multipliers", "Vmax"), "Vmax multiplier of the sensitivity variant (both arms)", num_fmt(1)),
             d1 = dv(VM, "schedule == 'D3'", "d_extrap20_ratio", 3, "", "(d) D3 ratio, Vmax x0.8 both arms (20,000 subjects)"),
             n1 = dcfg("trial_design.yaml", c("mc", "n_individual"), "subjects per model and variant", function(x) fnum(as.numeric(x), 0, big = TRUE)),
             n2 = dint("individual200k/criterion_d_200k_vmax080_both.csv", "TRUE", "n_subjects", "subjects in the criterion (d) re-evaluation, Vmax x0.8"),
             d2 = dv("individual200k/criterion_d_200k_vmax080_both.csv", "TRUE", "extrap_gt20_ratio", 3, "", "(d) D3 ratio at 200,000 subjects, Vmax x0.8 both arms"),
             ci = dspan(VM, "schedule == 'D3'", "d_ratio_boot_lo", "d_ratio_boot_hi", 2, "", "(d) D3 bootstrap 95% interval, Vmax x0.8 both arms"),
             d = thr$d, ab = dv(VM, "schedule == 'D3'", "d_abs_change_per_arm", 2, "", "(d) D3 absolute change in subjects above 20% extrapolation per arm, Vmax x0.8", scale = -1),
             vis = dint(VM, "schedule == 'D3'", "added_visits_total", "added visits, D3"))
  bl <- tx("A7.bullets", c(b, vs))
  cs_ <- .read("config/prereg_20260926.yaml")$section4$criteria_sets
  premise(all(c(cs_$ii$extrap_max_pct, cs_$iii$extrap_max_pct) == .read("config/nca_rules.yaml")$standard$reliability$extrap_max_pct), "sets (ii) and (iii) use the same extrapolation limit as the NCA flag (caption)")
  cap <- tx("A7.caption", c(vs, list(x = thr$x, r2i = f_set("ii", "r2"), sp2 = f_set("ii", "span"), r2 = f_set("iii", "r2"))))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  BY <- y0 + TH + 0.02
  deck_bullets(bl, box = c(GEO$ML, BY, GEO$CW, capy - 0.02 - BY), size = 18, gap_pt = 3)

  # ---- 노트 ----
  D200 <- function(v) sprintf("individual200k/criterion_d_200k_%s.csv", v)
  d2 <- function(v, item) dci(D200(v), "TRUE", "extrap_gt20_ratio", "d_ratio_boot_lo", "d_ratio_boot_hi", 3, "", item)
  for (v in A7_V) premise(!rows(D200(v))$d_meets_point, paste("criterion (d) not met at 200,000 subjects,", v))
  premise(grepl("AUC0-inf", .read("config/schedule_decision.yaml")$day50_removal) && grepl("fallback", .read("config/schedule_decision.yaml")$day50_removal), "recorded reason for keeping Day 50")
  deck_notes(tx("A7.notes", list(
    a = thr$a, b = thr$b, c = thr$c, d = thr$d, x = thr$x, ci = thr$ci, ke = thr$ke,
    cost = dcfg("trial_design.yaml", c("decision_rule", "cost_per_added_point", "subjects"), "added visits per added sample point", num_fmt(0)),
    nsub = dcfg("trial_design.yaml", c("mc", "n_individual"), "subjects per model and variant", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    ntr = dcfg("trial_design.yaml", c("mc", "n_trials"), "trials per model and variant", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    date = dcfg("schedule_decision.yaml", "decision_date", "schedule decision date", function(x) as.character(x)),
    b0 = f_study_days("B0", "all"), n = f$n, last = f$last,
    h = { d0 <- min(unlist(td$schedules$B0$days)); premise(d0 < 1 && sum(a7_days("B0") != round(a7_days("B0"))) == 1, "one B0 sample within the first day (the others on whole study days)")
      dderived("first B0 sample, hours after dose", "config/trial_design.yaml", "schedules.B0.days :: min(days) x 24", d0 * 24, fnum(d0 * 24, 0)) },
    vis = paste(vapply(c("D1", "D2", "D3", "D4", "Bminus"), function(s) dint(a7_sd("base"), sprintf("schedule=='%s'", s), "added_visits_total", sprintf("added visits, %s", s)), ""), collapse = ", "),
    amax = b$amax, cii = b$cii, r3 = b$r3, d50 = b$d50, bma = b$bma, bmd = a7_pair("Bminus", "d_extrap20_ratio", 2, "(d) ratio of the share above 20% extrapolation when Day 50 is removed"),
    bmc = a7_pair("Bminus", "c_reliable_gain_pp", 2, "(c) reliability change, set (ii), when Day 50 is removed (points)"),
    irng = drange(RP, RW, "gain_i_pp", 2, "", "reliability change set (i), D1 to D4, two models (points)"),
    iirng = drange(RP, RW, "gain_ii_pp", 2, "", "reliability change set (ii), D1 to D4, two models (points)"),
    b16 = dv(RL, "model=='k2016' & schedule=='B0'", "pct_reliable_iii", 1, "%", "share meeting set (iii), B0, 2016"),
    b20 = dv(RL, "model=='k2020' & schedule=='B0'", "pct_reliable_iii", 1, "%", "share meeting set (iii), B0, 2020"),
    d3 = drange(RL, "schedule=='D3'", "diff_vs_B0_pp", 1, "%p", "change with D3, two models"),
    n3 = dint(RL, "model=='k2016' & schedule=='B0'", "n", "virtual subjects per model"),
    dboot = paste(vapply(A7_V, function(v) dci(a7_sd(v), "schedule=='D3'", "d_extrap20_ratio", "d_ratio_boot_lo", "d_ratio_boot_hi", 2, "", sprintf("(d) D3 ratio with bootstrap interval, %s", v)), ""), collapse = " / "),
    n200 = dint(D200("base"), "TRUE", "n_subjects", "subjects in the criterion (d) re-evaluation"),
    e16 = d2("base", "criterion (d) at 200,000 subjects, D3, k2016"), e20 = d2("struct2020", "criterion (d) at 200,000 subjects, D3, k2020"),
    st = dderived("sampling-schedule decision rule: timing status in the prespecification register", PSR, "grepl('^Sampling-schedule decision rule', Item) :: Status", ps$Status, ps$Status),
    sd = dderived("sampling-schedule decision rule: date (evidence) in the prespecification register", PSR, "grepl('^Sampling-schedule decision rule', Item) :: Date (evidence)", ps[["Date (evidence)"]], ps[["Date (evidence)"]]),
    vm = vs$vm, d1 = vs$d1, n1 = vs$n1, n2 = vs$n2, d2 = vs$d2, bci = vs$ci, ab = vs$ab, v3 = vs$vis,
    d2ci = dspan("individual200k/criterion_d_200k_vmax080_both.csv", "TRUE", "d_ratio_boot_lo", "d_ratio_boot_hi", 3, "", "(d) D3 bootstrap 95% interval at 200,000 subjects, Vmax x0.8"),
    pd = dcfg("prereg_20260926.yaml", "registered_on", "registration date of the v1.0.1 pre-registration", function(x) as.character(x)))))
  deck_end()
}
