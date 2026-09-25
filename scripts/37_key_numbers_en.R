#!/usr/bin/env Rscript
# results/key_numbers_en.md: key numbers in English with their source files (review 2026-09-24 section 6).
# Rules as summary_en.md: English only, no em-dash; proportions with 95% intervals where they come from trials.
source("R/00_setup.R"); source_project()
R <- function(...) { f <- proj_path("results", ...); if (file.exists(f) || file.exists(paste0(f, ".gz"))) read_raw(f) else NULL }
# 반올림은 fmt_num과 같은 round_half_away(0에서 먼 쪽): summary_en.md·mechanism_summary와 같은 값이 같은 표기로 찍히도록
fd <- function(x, d) formatC(round_half_away(x, d), format = "f", digits = d)
f1 <- function(x) fd(x, 1); f2 <- function(x) fd(x, 2); f3 <- function(x) fd(x, 3)
ci <- function(e, lo, hi, d = 1) sprintf("%s%% (95%% CI %s to %s)", fd(e, d), fd(lo, d), fd(hi, d))
out <- c("# Key numbers (dupilumab biosimilar Phase 1 PK simulation)", "", sprintf("Generated %s. Each line names its source file under results/.", format(Sys.Date())), "")
add <- function(section, ...) out <<- c(out, if (!is.null(section)) c(paste0("## ", section), ""), ...)
ML <- c(k2016 = "2016 model", k2020 = "Model 1", base = "2016 model", struct2020 = "Model 1")
# 신뢰 충족률·탈락률과 그 증감은 reliability_facts()에서(전제 검사 포함, 어긋나면 중단). 플래그 세트 (i) 주값, (ii) 대괄호 병기
RF <- reliability_facts()
FPe <- function(i, ii, d = 1, unit = "%", signed = FALSE, each = FALSE) fmt_flag_pair(i, ii, d, unit, " to ", signed, each)
REL_CONV_EN <- "span ratio at least 2 is a convention used only by some statistical analysis plans (not a Phoenix feature)"

ev <- R("nca_engine", "engine_validation_summary.csv"); dr <- R("nca_engine", "dropout_reasons_B0.csv"); ed <- R("nca_engine", "engine_difference_individual_B0.csv")
if (!is.null(ev)) add("NCA engine (Phoenix WinNonlin-compatible)",
  sprintf("- Own engine vs NonCompart 0.8.4 vs PKNCA 0.12.1: lambda-z points identical in %s of %s profile comparisons; largest relative parameter difference %s (nca_engine/engine_validation_summary.csv).",
          format(sum(ev$lz_points_identical), big.mark = ","), format(sum(ev$n_profiles), big.mark = ","), formatC(max(ev$max_rel_diff), format = "e", digits = 1)),
  if (!is.null(RF)) sprintf("- AUCinf reliability failure at B0 (lambda-z not estimable or flagged; flag set (i) = adjusted R-squared at least 0.80 and extrapolation at most 20%%, (ii) adds span ratio at least 2): %s (2016 model) and %s (Model 1); span ratio below 2 flags %s%% and %s%% (reliability/dropout_reasons_by_schedule.csv, reliability/reliability_two_flag_sets_summary.csv).",
          FPe(RF$b0[model == "k2016", fail_i], RF$b0[model == "k2016", fail_ii]), FPe(RF$b0[model == "k2020", fail_i], RF$b0[model == "k2020", fail_ii]),
          fmt_num(RF$b0[model == "k2016", flag_span]), fmt_num(RF$b0[model == "k2020", flag_span]))
  else if (!is.null(dr)) sprintf("- AUCinf reliability failure at B0, flag set (ii) (any flag or lambda-z not estimable): %s%% (2016 model) and %s%% (Model 1); span ratio below 2 flags %s%% and %s%% (nca_engine/dropout_reasons_B0.csv).",
          f1(dr[model == "base", any_flag_or_fail_pct]), f1(dr[model == "struct2020", any_flag_or_fail_pct]), f1(dr[model == "base", flag_span_pct]), f1(dr[model == "struct2020", flag_span_pct])),
  if (!is.null(RF) && !is.null(RF$eng)) sprintf("- Reliability rate, previous vs new engine on the same observations (separate B0-only draw): %s; the span ratio flag alone accounts for %s of the drop (nca_engine/engine_difference_individual_B0.csv).",
          paste(sprintf("%s%% to %s (%s)", fmt_num(RF$eng$old), FPe(RF$eng$new_i, RF$eng$new_ii, each = TRUE), ML[RF$eng$model]), collapse = ", "), paste(sprintf("%s%%", fmt_num(RF$eng$span_share_pct)), collapse = " and "))
  else if (!is.null(ed)) sprintf("- Reliability rate, previous vs new engine, flag set (ii): %s%% to %s%% (2016 model), %s%% to %s%% (Model 1) (nca_engine/engine_difference_individual_B0.csv).",
          f1(ed[model == "base"][1, reliable_pct]), f1(ed[model == "base"][2, reliable_pct]), f1(ed[model == "struct2020"][1, reliable_pct]), f1(ed[model == "struct2020"][2, reliable_pct])), "")

# 신뢰 플래그 두 세트와 촘촘한 후기 채혈(검토 의견 W2 §1)
st <- R("reliability", "reliability_span_transition.csv")
if (!is.null(RF)) { b <- RF$b0; d3 <- RF$dense$d3
  lost <- if (!is.null(st)) st[grp == "lost_span" & schedule %in% c("D1", "D2", "D3", "D4") & variant %in% c("base", "struct2020", "vmax080_both", "vmax125_both")] else NULL
  rg <- function(x, d = 1) sprintf("%s to %s", fmt_num(min(x), d), fmt_num(max(x), d))
  add("AUCinf reliability: two flag sets (20,000 subjects per model, B0 conditions)",
      sprintf("- Flag set (i): lambda-z estimable, adjusted R-squared at least 0.80, extrapolation at most 20%%; flag set (ii): (i) plus span ratio at least 2; %s. Cited as (i) [(ii)].", REL_CONV_EN),
      sprintf("- Reliable at B0: %s (2016 model, 95%% CI %s to %s [%s to %s]) and %s (Model 1, 95%% CI %s to %s [%s to %s]); loss from the span ratio alone %s and %s percentage points (reliability/reliability_two_flag_sets_summary.csv).",
              FPe(b[model == "k2016", rel_i], b[model == "k2016", rel_ii]), fmt_num(b[model == "k2016", rel_i_lo]), fmt_num(b[model == "k2016", rel_i_hi]), fmt_num(b[model == "k2016", rel_ii_lo]), fmt_num(b[model == "k2016", rel_ii_hi]),
              FPe(b[model == "k2020", rel_i], b[model == "k2020", rel_ii]), fmt_num(b[model == "k2020", rel_i_lo]), fmt_num(b[model == "k2020", rel_i_hi]), fmt_num(b[model == "k2020", rel_ii_lo]), fmt_num(b[model == "k2020", rel_ii_hi]),
              fmt_num(b[model == "k2016", span_loss], 2), fmt_num(b[model == "k2020", span_loss], 2)),
      sprintf("- Two-model range at B0: reliable %s; failing %s.", FPe(RF$rng$rel_i, RF$rng$rel_ii), FPe(RF$rng$fail_i, RF$rng$fail_ii)),
      sprintf("- D3 versus B0 (paired, same subjects): %s (reliability/reliability_paired_vs_B0.csv; (ii) equals c_reliable_gain_pp of trials/schedule_decision_<variant>.csv).",
              paste(sprintf("%s %s", ML[d3$model], sprintf("(i) %s [(ii) %s] percentage points", sprintf("%s (95%% CI %s to %s)", fmt_num(d3$gain_i_pp, 2, TRUE), fmt_num(d3$gain_i_lo, 2), fmt_num(d3$gain_i_hi, 2)),
                                                            sprintf("%s (95%% CI %s to %s)", fmt_num(d3$gain_ii_pp, 2, TRUE), fmt_num(d3$gain_ii_lo, 2), fmt_num(d3$gain_ii_hi, 2)))), collapse = "; ")),
      sprintf("- D1, D2 and D4 versus B0, both models: %s; in all of D1 to D4 the loss from the span ratio alone rises by %s percentage points; sign reversed by the span ratio in %s.",
              FPe(RF$dense$other_i, RF$dense$other_ii, 2, " percentage points", TRUE), rg(RF$dense$rise, 2), if (nrow(RF$flips)) paste(sprintf("%s %s", ML[RF$flips$model], RF$flips$schedule), collapse = ", ") else "no cell"),
      if (!is.null(RF$lz_min)) sprintf("- Shortest 3-point lambda-z window from nominal days after Day 22: %s days at B0 and %s days in D1 to D4, i.e. half-life limits of %s and %s days for span ratio 2 (reliability/reliability_lz_window_by_schedule.csv).",
              RF$lz_min[["B0"]], RF$lz_min[["dense"]], fmt_num(RF$lz_min[["B0"]] / 2), fmt_num(RF$lz_min[["dense"]] / 2)),
      if (!is.null(lost) && nrow(lost)) sprintf("- Subjects reliable at B0 and span-flagged with added sampling (4 model variants x D1 to D4): window shorter in %s%%, median window %s days at B0 versus %s days, median half-life %s days versus %s days (longer in only %s%%), added sampling day in the new window in %s%% (reliability/reliability_span_transition.csv).",
              rg(lost$window_shorter_pct), rg(lost$window_B0_median), rg(lost$window_median), rg(lost$HL_B0_median, 2), rg(lost$HL_median, 2), rg(lost$HL_longer_pct), rg(lost$window_with_added_day_pct)),
      sprintf("- Criterion (c) (gain of at least %s percentage points) is met under neither flag set for any variant or schedule; the recommendation is the same under both sets (reliability/reliability_paired_vs_B0.csv).",
              { th <- unique(RF$all$crit_c_threshold_pp); stopifnot(length(th) == 1, is.finite(th)); format(th) }),
      "- Numbers that change with the span ratio criterion: reliability and dropout rates, the change in reliability with added sampling (criterion (c)), the rule A analysis set and the pass rates of G2, F3-A and F3-C (rules A and C). Unchanged: criterion (d), AUClast, Cmax, rule B and the schedule recommendation.", "") }

g2f <- R("oc", "g2_rules_flags.csv")
if (!is.null(g2f)) { x <- g2f[region == "boundary" & label %in% c("P2", "G2-A(ii)", "G2-B", "G2-C(ii)", "G2-A(i)", "G2-C(i)")]
  y <- x[!is.na(pass_pct), .SD[which.max(pass_pct)], by = .(model, label)][, .(model, label, pass_pct, lo, hi, scenario, n_trials)]
  cnt <- x[!is.na(pass_pct), .(n_gt5 = sum(pass_pct > 5), n_sc = .N), by = .(model, label)]; y <- merge(y, cnt, by = c("model", "label"))
  na_sets <- unique(x[is.na(pass_pct), label])
  add("G2 by AUCinf handling rule and flag set (boundary scenarios, true AUC0-inf ratio 0.80 or 1.25)",
      sprintf("- Largest boundary pass rate: %s (oc/g2_rules_flags.csv).", paste(sprintf("%s %s %s in %s (%s trials; above 5%% in %d of %d)", ML[y$model], y$label, ci(y$pass_pct, y$lo, y$hi, 2), y$scenario, formatC(y$n_trials, format = "d", big.mark = ","), y$n_gt5, y$n_sc), collapse = "; ")),
      if (length(na_sets)) sprintf("- Not yet computed (rejudge file absent): %s.", paste(na_sets, collapse = ", ")), "") }

cs <- R("cliff", "cliff_summary.csv"); cp <- R("cliff", "cliff_points.csv")
if (!is.null(cs)) { b <- cs[model == "k2016" & weight == "base"]; n1 <- cp[model == "k2016" & weight == "base" & timing == "nominal"]
  add("Sampling cliff (2016 model, 60 to 90 kg, 20,000 subjects)",
      sprintf("- True LLOQ reached at study Day %s (5th to 95th percentile %s to %s); after Day 58 in %s%% (cliff/cliff_summary.csv).", f1(b$lloq_studyday_median), f1(b$lloq_studyday_p05), f1(b$lloq_studyday_p95), f1(b$lloq_after_day58_pct)),
      sprintf("- Cliff length: %s days (1-day definition), %s days (2-day definition); cliff starts at %s mg/L (median).", f2(b$len1_median), f2(b$len2_median), f2(b$c_start1_median)),
      sprintf("- Two or more samples in the cliff at nominal days, fixed schedules, 1-day definition: at most %s%%; three or more with daily Day 29 to 57 sampling: %s%% (cliff/cliff_points.csv).",
              f1(max(n1[definition_day == 1 & schedule != "daily_29_57", pct_ge2])), f1(n1[definition_day == 1 & schedule == "daily_29_57", pct_ge3])), "") }

bt <- R("oc", "boundary_type1.csv"); pw <- R("oc", "power.csv"); rk <- R("oc", "random_space_risks.csv"); pr <- R("oc", "prereg.csv")
# 보고 우선순위(지시 2026-09-25 §4): 1차 지표 = 경계 1종 오류, 무작위 제품 공간 = 보조 지표. 최대 칸·연장·Km·VM150·KE120은 key_facts(전제 검사 포함)
oc_cfg <- read_cfg("oc_design.yaml"); KF <- key_facts(oc = oc_cfg)
if (!is.null(bt)) {
  x <- bt[config == "P2"][which.max(pass_pct)]; g <- bt[config == "G2"][which.max(pass_pct)]
  stopifnot(!is.null(KF$p2max), identical(KF$p2max$row$scenario, x$scenario), identical(KF$p2max$row$model, x$model))
  add("Operating characteristics (pre-specified design; primary metric: boundary type I error)",
      if (!is.null(pr)) sprintf("- Design pre-registered in commit %s (oc/prereg.csv).", substr(pr$prereg_commit, 1, 7)),
      sprintf("- Largest boundary type I error of P2 (AUClast + Cmax): %s, %s, %s %s, true AUC0-inf ratio %s; %s (oc/boundary_type1.csv%s).", ci(x$pass_pct, x$lo, x$hi, 2), ML[[x$model]], x$mechanism, x$direction, f2(x$target),
              key_trials_text(KF$p2max$row, KF$p2max$ext, "en"), if (!is.null(KF$p2max$ext)) ", oc/extension_decision.csv" else ""),
      sprintf("- Largest boundary type I error of G2 (AUCinf + Cmax): %s, %s, %s %s, true ratio %s; %s.", ci(g$pass_pct, g$lo, g$hi, 2), ML[[g$model]], g$mechanism, g$direction, f2(g$target), key_trials_text(KF$g2max$row, KF$g2max$ext, "en")),
      sprintf("- Boundary scenarios with P2 above 5%%: %d of %d.", nrow(bt[config == "P2" & pass_pct > 5]), nrow(bt[config == "P2"])))
  if (!is.null(pw)) { s0 <- pw[scenario == "S00" & config %in% c("P2", "F3A", "G2")]
    out <<- c(out, sprintf("- Power for identical products (%s trials): %s (oc/power.csv).", if ("n_trials" %in% names(s0)) paste(unique(formatC(s0$n_trials, format = "d", big.mark = ",")), collapse = " or ") else "10,000", paste(sprintf("%s %s %s", ML[s0$model], c(P2 = "P2", F3A = "F3-A", G2 = "G2")[s0$config], ci(s0$pass_pct, s0$lo, s0$hi)), collapse = "; "))) }
  if (!is.null(rk)) { r <- rk[truth == "AUC0-inf" & scope == "전체" & config %in% c("P2", "G2")]
    out <<- c(out, sprintf("- Random product space (secondary metric; depends on the assumed virtual product distribution, log-uniform multipliers %s), consumer risk (truth outside, %s products): %s (oc/random_space_risks.csv).",
                          random_space_ranges_text(oc_cfg, " to "), format(oc_cfg$random_space$n_products, big.mark = ","),
                          paste(sprintf("%s %s %s", ML[r[startsWith(metric, "소비자")]$model], r[startsWith(metric, "소비자")]$config, ci(r[startsWith(metric, "소비자")]$pct, r[startsWith(metric, "소비자")]$lo, r[startsWith(metric, "소비자")]$hi, 2)), collapse = "; "))) }
  out <- c(out, "")
}

# Pillar 3 첫 문장(결합 상수 Km, 역산)과 예비 제품 시나리오(VM150, KE120): key_facts의 수치와 전제 검사
if (!is.null(KF$km)) add("Binding constant (Pillar 3 lead sentence)",
  sprintf("- %s Sources: oc/inversion_all.csv, oc/inversion_scan_k2016_Km.csv, oc/inversion_scan_k2020_Km.csv.", key_km_sentence(KF$km, "en")), "")
if (!is.null(KF$vm150) || !is.null(KF$ke120)) { v <- KF$vm150; k2 <- KF$ke120; d_ <- function(e) if (e < 1) 3 else 2
  add("Preliminary product scenarios (arbitrary multipliers, 2016 model; AUCinf reliability under flag set (ii))",
      if (!is.null(v)) sprintf("- VM150 (Vmax x1.50), %s trials: AUCinf alone (reliable subjects) passes in %s, above the nominal 5%% although the true AUC0-inf ratio %s is outside the limits; AUClast alone passes in %s (trials5000/products5000_props_base.csv; true ratio = mean trial GMR of the individual model AUC0-inf, fallback/consumer_risk.csv).",
                               formatC(v$n_trials, format = "d", big.mark = ","), ci(v$rel$est, v$rel$lo, v$rel$hi, d_(v$rel$est)), f3(v$true_ratio), ci(v$last$est, v$last$lo, v$last$hi, d_(v$last$est))),
      if (!is.null(k2)) sprintf("- KE120 (ke x1.20), %s trials: AUClast pass with AUCinf (reliable subjects) fail in %s; the true AUC0-inf ratio %s is inside the limits, so these are false negatives of AUCinf (fallback/discordance_classification.csv).",
                                formatC(k2$n_trials, format = "d", big.mark = ","), ci(k2$est, k2$lo, k2$hi, d_(k2$est)), f3(k2$true_ratio)), "") }

dec <- R("trials", "schedule_decision_base.csv"); k <- R("individual200k", "criterion_d_200k_base.csv")
if (!is.null(dec)) { d3 <- dec[schedule == "D3"]
  add("Sampling density decision (2016 model)",
      sprintf("- Final schedule B0. Candidates meeting any pre-specified criterion: %s (trials/schedule_decision_base.csv).", if (any(dec$recommend)) paste(dec[recommend == TRUE, schedule], collapse = ", ") else "none"),
      sprintf("- D3 versus B0: AUClast CI width change %s%% (positive = narrower), reliability gain %s percentage points, NCA extrapolation above 20%% ratio %s%s.", f2(100 * d3$a_mean_width_rel_decrease),
              { pv <- if (!is.null(RF)) RF$dense$d3[variant == "base"] else NULL; if (!is.null(pv) && nrow(pv) == 1) FPe(pv$gain_i_pp, pv$gain_ii_pp, 2, "", TRUE) else sprintf("[(ii) %s]", f2(d3$c_reliable_gain_pp)) }, f3(d3$d_extrap20_ratio),
              if (!is.null(k)) sprintf("; 200,000 subjects %s (95%% CI %s to %s)", f3(k[schedule == "D3", extrap_gt20_ratio]), f3(k[schedule == "D3", d_ratio_boot_lo]), f3(k[schedule == "D3", d_ratio_boot_hi])) else ""), "") }

p1 <- R("rationale", "pillar1_coverage_B0.csv")
if (!is.null(p1)) { x <- p1[group == "전체"]
  add("Coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 subjects per model)",
      sprintf("- True extrapolated share median %s%% (2016) and %s%% (Model 1); 95th percentile %s%% and %s%%; true coverage below 80%%: %s and %s (rationale/pillar1_coverage_B0.csv).",
              f2(x[model == "k2016", extrap_true_median]), f2(x[model == "k2020", extrap_true_median]), f2(x[model == "k2016", extrap_true_p95]), f2(x[model == "k2020", extrap_true_p95]),
              x[model == "k2016", coverage_lt80_pct_ci], x[model == "k2020", coverage_lt80_pct_ci]), "") }

txt <- paste(out, collapse = "\n")
if (grepl("—|–|\u2212", txt)) stop("key_numbers_en.md contains an em-dash, en-dash or minus sign (U+2212)")
if (grepl("[가-힣]", txt)) { bad <- regmatches(txt, gregexpr("[^\n]*[가-힣][^\n]*", txt))[[1]]; stop("key_numbers_en.md contains Korean text: ", paste(head(bad, 3), collapse = " || ")) }
writeLines(txt, proj_path("results", "key_numbers_en.md"), useBytes = TRUE)
cat(txt, "\n")
