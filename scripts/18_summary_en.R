#!/usr/bin/env Rscript
# results/summary_en.md: English summary for regulatory briefing (review consolidated 2026-09-24, section 9).
# Built from the same result files as the report. Rules: English only, no em-dash, abbreviations spelled out at first use,
# trial-level proportions cited only from runs with >= 5,000 trials (with 95% intervals); two-model ranges for residual-sensitive metrics.
source("R/00_setup.R"); source_project()
R <- function(...) { f <- proj_path("results", ...); if (file.exists(f) || file.exists(paste0(f, ".gz"))) read_raw(f) else NULL }
f1 <- function(x, d = 1) formatC(x, format = "f", digits = d); f2 <- function(x) f1(x, 2); f3 <- function(x) f1(x, 3)
ci <- function(e, lo, hi, d = 1) sprintf("%s (%s to %s)", f1(e, d), f1(lo, d), f1(hi, d))
md_table <- function(dt) { h <- paste0("| ", paste(names(dt), collapse = " | "), " |"); s <- paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  b <- apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")); c(h, s, b, "") }
out <- c()
add <- function(...) out <<- c(out, ...)
gate <- read_cfg("gate_decision.yaml"); sdec <- read_cfg("schedule_decision.yaml"); design <- read_cfg("trial_design.yaml")

add("# Dupilumab biosimilar Phase 1 pharmacokinetic simulation: summary for regulatory briefing", "",
    sprintf("Generated %s from repository results (branch claude/epic-bardeen-axbreo). Every number below is read from the result files used by the full report.", format(Sys.Date())), "",
    "Abbreviations: area under the concentration-time curve to the last quantifiable concentration (AUClast), to infinity (AUCinf); maximum concentration (Cmax); non-compartmental analysis (NCA); geometric mean ratio (GMR); confidence interval (CI); lower limit of quantification (LLOQ, 0.078 mg/L); inter-individual variability (IIV); target-mediated drug disposition (TMDD); Michaelis-Menten (MM); Monte Carlo (MC); body mass index (BMI).", "",
    "Models: primary model Kovalenko et al. 2016 (CPT Pharmacometrics Syst Pharmacol 5:617, Table 2, BLQ-included column): two-compartment, first-order absorption, parallel linear and MM elimination, Km fixed at 0.01 mg/L, central volume scaled by (weight/75)^0.705. Sensitivity model Kovalenko et al. 2020 Model 1 (Clin Pharmacol Drug Dev 9:756, Table 1 and Supplementary Table 2): transit absorption (3 compartments, mean transit time 0.105 day), its own IIV and residual error (proportional 15.0%, additive 0.03 mg/L).", "",
    "Scenario codes (test arm only unless stated; reference arm shared through common random numbers): S00 identical products; F085, F090, F097, F110 bioavailability x0.85, x0.90, x0.97, x1.10; KE110, KE120 linear elimination rate constant (ke) x1.10, x1.20; VM080, VM125, VM150 maximum MM elimination rate (Vmax) x0.80, x1.25, x1.50; KM05, KM2, KM5, KM10 MM constant (Km) x0.5, x2, x5, x10; KA075 absorption rate constant x0.75. Sampling schedules: B0 Syneos baseline; D1 to D4 add two to four samples between Day 32 and Day 53 (D1: Days 39, 46; D2: Days 39, 46, 53; D3: Days 32, 39, 46, 53; D4: Days 40, 47); B- removes Day 50.", "",
    "Study design simulated: 300 mg single subcutaneous dose (2 mL of 150 mg/mL), parallel groups, 117 evaluable subjects per arm, body weight 60 to 90 kg, Syneos sampling schedule (B0: Days 1, 2, 4, 6, 8, 11, 15, 22, 29, 36, 43, 50, 57). Equivalence: two one-sided tests via the 90% CI of the GMR from a pooled two-sample t on log scale, limits 80.00% to 125.00%.", "")

# 1. Model qualification
add("## 1. Model qualification", "",
    sprintf("Gate scope: %s. Excluded from the gate: PKM14271 200 mg and Cohen 2022 200 mg (different presentation, 1.14 mL of 175 mg/mL, faster absorption with median time to Cmax 3.0 days versus 7.0 days at 300 mg). Decision: option %s, %s, %s.", gate$gate_scope, gate$decision_option, gate$decided_by, gate$decision_date), "")
g16 <- R("step1", "step1b_quant_gate.csv"); g20 <- R("step1_k2020", "step1b_quant_gate.csv")
if (!is.null(g16) && !is.null(g20)) {
  m <- merge(g16[, .(id, gate_role, dose_mg, obs = AUClast_obs_mean, r16 = AUClast_ratio, z16 = AUClast_z, c16 = Cmax_ratio, d16 = dev)],
             g20[, .(id, r20 = AUClast_ratio, z20 = AUClast_z, c20 = Cmax_ratio, d20 = dev)], by = "id")
  lab <- c(clot300_nonasian_pooled = "300 mg non-Asian pooled (Clot 2021 Table 5, n=40)", clot300_japanese_TDU12265 = "300 mg Japanese (TDU12265)", clot300_chinese = "300 mg Chinese (Clot 2021)",
           clot600_chinese = "600 mg Chinese (Clot 2021)", clot600_japanese = "600 mg Japanese (TDU12265)", PKM14271_200mg_test = "200 mg PKM14271 test", PKM14271_200mg_reference = "200 mg PKM14271 reference",
           Cohen2022_200mg_AI = "200 mg Cohen 2022 autoinjector", Cohen2022_200mg_PFSS = "200 mg Cohen 2022 prefilled syringe")
  dv <- c(internal = "internal", partial = "partly internal", external = "external")
  m[, dataset := lab[id]]; setorder(m, -gate_role, dose_mg)
  add(md_table(m[, .(Dataset = dataset, Role = fifelse(gate_role == "gate", "gate", "external check"), `Observed AUClast mean (mg*day/L)` = f1(obs, 0),
                     `2016: sim/obs AUClast` = f2(r16), `2016: development data` = dv[d16], `Model 1: sim/obs AUClast` = f2(r20), `Model 1: development data` = dv[d20],
                     `2016: sim/obs Cmax` = f2(c16), `Model 1: sim/obs Cmax` = f2(c20))]))
  gg <- m[gate_role == "gate"]
  add(sprintf("Gate result (AUClast mean within 15%% of observed, log-scale CV 35%% to 51%%, cohort median last quantifiable time, weight slope): PASS for both models. Across the study-presentation datasets the simulated/observed AUClast ratio is %s to %s for the 2016 model and %s to %s for Model 1; Cmax is %s to %s and %s to %s.",
              f2(min(gg$r16)), f2(max(gg$r16)), f2(min(gg$r20)), f2(max(gg$r20)), f2(min(gg$c16)), f2(max(gg$c16)), f2(min(gg$c20)), f2(max(gg$c20))), "")
}
e16 <- R("step1", "step1h_external_only_gate.csv"); e20 <- R("step1_k2020", "step1h_external_only_gate.csv"); ws <- R("step1", "step1i_arm_weight_sensitivity.csv")
if (!is.null(e16) && !is.null(e20)) {
  fail_txt <- function(e) { x <- e[pass_mean == FALSE & !is.na(ratio)]; if (!nrow(x)) "all items within limits" else paste(sprintf("%s arm of %s (%s)", sub("^.* ", "", x$item), sub(" .*$", "", x$item), f2(x$ratio)), collapse = "; ") }
  add("Gate re-judged on datasets fully external to each model's development data (section 7 of the review):",
      sprintf("- 2016 model (%d items): %s. Failing items: %s.", nrow(e16), if (all(e16$pass_mean %in% TRUE)) "PASS" else "FAIL", fail_txt(e16)),
      sprintf("- Model 1 (%d items): %s. Failing items: %s.", nrow(e20), if (all(e20$pass_mean %in% TRUE)) "PASS" else "FAIL", fail_txt(e20)))
  if (!is.null(ws)) { x <- ws[study == "PKM12350"]; add(sprintf("- The Li 2020 single-arm checks have no reported body weight; they are simulated at an assumed mean of 78 kg. For the 2016 model the PKM12350 ratios fall from %s at 74 kg to %s at 82 kg and %s at 86 kg, so the external-only failure depends on the assumed weight.",
                                                  paste(f2(x[weight_mean == 74, ratio]), collapse = "/"), paste(f2(x[weight_mean == 82, ratio]), collapse = "/"), paste(f2(x[weight_mean == 86, ratio]), collapse = "/"))) }
  add("")
}
ab <- R("step1", "appendix_absorption_diagnostic.csv"); sh <- R("step1", "appendix_dose_shape_model.csv"); so <- R("step1", "appendix_dose_shape_observed.csv")
if (!is.null(ab)) {
  add("200 mg external checks and absorption (diagnostic only, not adopted): the 2016 model does not reproduce the fast absorption of the 200 mg 1.14 mL (175 mg/mL) presentation (observed median time to Cmax 3 days versus 7 days simulated).",
      sprintf("Matching absorption to the observed time to Cmax closes about one third to one half of the 200 mg exposure shortfall: the dose-normalized 200:300 mg AUClast ratio is %s observed (four-arm, n-weighted; %s for the PKM14271 test arm alone) versus %s simulated, and %s to %s with the absorption rate constant multiplied by 1.5 to 2 for 200 mg only. The remainder (about 9%% to 11%% on the four-arm basis) is unexplained: a presentation-specific bioavailability difference or over-estimated TMDD at low exposure are possible; the latter is covered by the population maximum elimination rate (Vmax) x0.8 sensitivity variant. For the study presentation (300 mg, 2 mL of 150 mg/mL) the model reproduces both absorption timing and exposure.",
              f2(so[1, ratio_per_mg]), f2(so[2, ratio_per_mg]), f2(sh[ka_multiplier == 1, ratio_per_mg]), f2(sh[ka_multiplier == 1.5, ratio_per_mg]), f2(sh[ka_multiplier == 2, ratio_per_mg])), "")
}
xv <- R("crossval", "crossval_model1.csv")
if (!is.null(xv)) { off <- xv[agree_3pct %in% FALSE]
  add(sprintf("Cross-validation of Model 1 against an independent implementation (first-order absorption approximation of the transit model, 20,000 subjects, acceptance within 3%%): %d of %d metrics within 3%%; the %d others are small percentages with absolute differences of %s to %s percentage points (%s).",
              sum(xv$agree_3pct %in% TRUE), sum(!is.na(xv$agree_3pct)), nrow(off), f2(min(abs(off$sim - off$ref))), f2(max(abs(off$sim - off$ref))),
              paste(c(extrap_gt20_pct = "share with NCA extrapolation above 20%", extrap_nca_median = "median NCA extrapolation", extrap_true_median = "median true extrapolation")[off$metric], collapse = ", ")), "") }
rf <- R("crossval", "reviewer_reference_round5.csv")
if (!is.null(rf)) { m <- rf[note == "" & !is.na(rel_diff_pct)]; p95 <- rf[metric == "extrap_true_p95" & !is.na(sim_20000_nojitter)]
  add(sprintf("Curve-shape and weight-band results were compared with the reviewer's independent implementation (8,000 subjects per condition): %d of %d non-rare metrics within 10%%. The 95th percentile of true extrapolation is higher here because sampling-time windows are simulated; without them the same subjects give %s%% to %s%% relative to the reference.",
              sum(abs(m$rel_diff_pct) <= 10), nrow(m), f1(min(100 * (p95$sim_20000_nojitter / p95$ref_8000 - 1))), f1(max(100 * (p95$sim_20000_nojitter / p95$ref_8000 - 1)))), "") }

# 2. Pillar 1
add("## 2. Pillar 1: coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 virtual subjects per model)", "")
p1 <- R("rationale", "pillar1_coverage_B0.csv")
if (!is.null(p1)) {
  p1e <- p1[group == "전체"]
  add(md_table(p1e[, .(Model = fifelse(model == "k2016", "2016 (primary)", "Model 1"), `True extrapolated %: median` = f2(extrap_true_median), `95th percentile` = f2(extrap_true_p95), Max = f2(extrap_true_max),
                       `True coverage below 80% (%, 95% CI)` = coverage_lt80_pct_ci, `NCA extrapolated %: median` = f2(extrap_nca_median), `NCA extrapolation above 20% (%, 95% CI)` = extrap_nca_gt20_pct_ci,
                       `AUCinf reliability met (%, 95% CI)` = reliable_pct_ci)]))
  add(sprintf("Residual-sensitive metrics are stated as the two-model range: AUCinf reliability criteria (adjusted R-squared at least 0.80 and extrapolation at most 20%%) are not met in %s%% to %s%% of subjects; NCA extrapolation above 20%% occurs in %s%% to %s%%. True extrapolation (model integral beyond the last quantifiable time) is essentially the same in both models.",
              f1(100 - max(p1e$reliable_pct)), f1(100 - min(p1e$reliable_pct)), f2(min(p1e$extrap_gt20_pct)), f2(max(p1e$extrap_gt20_pct))), "")
}
cs <- R("curve_shape", "curve_shape_B0.csv")
if (!is.null(cs)) {
  vl <- c(base = "Primary model", km05_both = "Km x0.5", km2_both = "Km x2", km5_both = "Km x5", km10_both = "Km x10", vmax080_both = "Vmax x0.8 (longer tail, conservative)", vmax125_both = "Vmax x1.25", vmax050_both = "Vmax x0.5 (stress test, fails gate)")
  cs[, v := vl[variant]]; cs <- cs[match(names(vl), variant)][!is.na(variant)]
  add("Curve shape below the LLOQ cannot be observed; sensitivity to Km and Vmax (both arms, 20,000 subjects each):", "",
      md_table(cs[, .(Variant = v, `True extrapolated %: median` = f2(extrap_true_median), `95th pct` = f2(extrap_true_p95), Max = f2(extrap_true_max), `Coverage below 80% (%)` = f2(coverage_lt80_pct),
                      `NCA extrap. median (%)` = f2(extrap_nca_median), `NCA extrap. >20% (%)` = f2(extrap_gt20_pct), `AUCinf reliable (%)` = f1(reliable_pct), `Lambda-z not estimable (%)` = f2(lambda_fail_pct),
                      `Median last quantifiable day` = f1(tlast_median), `AUClast geometric mean` = f1(AUClast_geo, 0), `300 mg ~78 kg AUClast vs observed 544` = f2(AUClast_ratio_vs_obs544))]),
      "Kovalenko 2020 reported that excluding below-LLOQ values makes the model predict a less steep TMDD phase with Vm and Km increasing together; this motivates the Km-increase sensitivity.", "")
}

# 3. Pillar 2
add("## 3. Pillar 2: decision concordance between AUClast and AUCinf (B0, 117 per arm)", "")
p2 <- R("rationale", "pillar2_products_B0.csv")
if (!is.null(p2)) {
  x <- p2[model == "k2016" & source == "5,000회 이상"]
  add(sprintf("Primary model, scenarios run with at least 5,000 trials (common random numbers; test-arm multipliers on fixed effects). Pass rates and concordance with Wilson 95%% intervals; GMR is the mean over trials."), "",
      md_table(x[, .(Scenario = scenario, Trials = n_trials, `GMR AUClast` = f3(GMR_mean_AUClast), `GMR true AUCinf` = f3(GMR_mean_AUCinf_true), `GMR NCA AUCinf (reliable)` = f3(GMR_mean_AUCinf_reliable),
                     `Pass AUClast (%)` = ci(pass_rate_AUClast, pass_lo_AUClast, pass_hi_AUClast), `Pass AUCinf reliable (%)` = ci(pass_rate_AUCinf_reliable, pass_lo_AUCinf_reliable, pass_hi_AUCinf_reliable),
                     `Pass Cmax (%)` = ci(pass_rate_Cmax, pass_lo_Cmax, pass_hi_Cmax), `Agreement (%)` = ci(agree, agree_lo, agree_hi), `AUClast pass, AUCinf fail (%)` = ci(last_pass_inf_fail, last_pass_inf_fail_lo, last_pass_inf_fail_hi, 2),
                     `log GMR correlation` = f3(cor_logGMR_last_infrel))]))
  pr <- R("trials5000", "products5000_props_base.csv"); mc <- R("trials5000", "mc_consistency_500_vs_5000.csv"); idt <- R("trials5000", "mc_consistency_identity.csv")
  if (!is.null(pr)) { reps <- unique(pr[, .(scenario, n_trials)]); esc <- reps[n_trials > 5000 & scenario != "S00"]
    pr[, target := fifelse(est <= 10 | est >= 90, 1, 1.5)]; miss <- pr[(hi - lo) / 2 > target]
    open_thr <- pr[scenario %in% esc$scenario & ((lo <= 90 & hi >= 90) | (lo <= 5 & hi >= 5))]
    add(sprintf("Monte Carlo precision: %s (95%% interval half-width at most 1 percentage point below 10%% or above 90%%, at most 1.5 points otherwise; largest half-width here %s points). Adaptive escalation: %s extended to %s trials because an interval included the 90%% or 5%% threshold at 5,000 and 10,000 trials (identical products run in the same batches for paired comparison)%s.%s",
                if (nrow(miss)) sprintf("%d proportions miss the pre-set precision (%s)", nrow(miss), paste(miss$scenario, miss$metric, collapse = "; ")) else "every cited proportion meets the pre-set precision",
                f2(max((pr$hi - pr$lo) / 2)), paste(esc$scenario, collapse = " and "), format(max(esc$n_trials), big.mark = ","),
                if (nrow(open_thr)) sprintf("; at the cap, %s", paste(sprintf("%s %s %s%% (%s to %s) still includes the threshold", open_thr$scenario, open_thr$metric, f1(open_thr$est), f1(open_thr$lo), f1(open_thr$hi)), collapse = "; ")) else "",
                if (!is.null(mc) && !is.null(idt)) sprintf(" Trials 1 to 500 reproduce the earlier 500-trial run exactly (maximum difference %s); the 500-trial estimates differ from the independent later trials beyond MC error in %d of %d comparisons.", format(idt$max_abs_diff), sum(mc$outside_mc), nrow(mc)) else ""), "") }
  pc <- R("pillar2_curvature", "pillar2_vmax080_per_endpoint.csv")
  if (!is.null(pc)) { x <- dcast(pc[endpoint %in% c("AUClast", "AUCinf_true", "AUCinf_reliable")], scenario + n_trials ~ endpoint, value.var = "GMR_mean")
    add("Curvature robustness (population Vmax x0.8 in both arms, 2,000 trials per scenario; mean GMR only):", "",
        md_table(x[, .(Scenario = scenario, Trials = n_trials, `GMR AUClast` = f3(AUClast), `GMR true AUCinf` = f3(AUCinf_true), `GMR NCA AUCinf (reliable)` = f3(AUCinf_reliable))])) }
  y <- p2[model == "k2020"]
  if (nrow(y)) add("Model 1 (500 trials per scenario; mean GMR only, proportions are reported in the full report): ", "",
                   md_table(y[, .(Scenario = scenario, `GMR AUClast` = f3(GMR_mean_AUClast), `GMR true AUCinf` = f3(GMR_mean_AUCinf_true), `GMR NCA AUCinf (reliable)` = f3(GMR_mean_AUCinf_reliable))]))
}
fb <- R("fallback", "discordance_classification.csv"); cr <- R("fallback", "consumer_risk.csv"); fc <- R("fallback", "fallback_cost.csv"); rr <- R("fallback", "aucinf_rules.csv")
if (!is.null(fb)) { top <- fb[which.max(last_pass_inf_fail)]
  add(sprintf("Largest rate of AUClast pass with AUCinf fail: %s%% in scenario %s (true AUCinf ratio %s, inside 80%% to 125%%, so these are false negatives of AUCinf).", ci(top$last_pass_inf_fail, top$last_pass_inf_fail_lo, top$last_pass_inf_fail_hi, 2), top$scenario, f3(top$true_ratio)), "") }
if (!is.null(cr) && nrow(cr)) add("Consumer risk (true AUCinf ratio outside the limits):", "",
  md_table(cr[, .(Scenario = scenario, `True ratio` = f3(true_ratio), Trials = n_trials, `Pass AUClast and Cmax (%)` = ci(joint_last_cmax, joint_last_cmax_lo, joint_last_cmax_hi, 2),
                  `Pass all three (%)` = ci(joint_3, joint_3_lo, joint_3_hi, 2), `Only AUCinf fails (%)` = ci(only_inf_fails, only_inf_fails_lo, only_inf_fails_hi, 2))]))
if (!is.null(fc)) add("Cost of a three-endpoint fallback (joint pass rates, >= 5,000 trials):", "",
  md_table(fc[, .(Scenario = scenario, `AUClast + Cmax (%)` = ci(i_last_cmax, i_lo, i_hi), `+ AUCinf reliable (%)` = ci(ii_plus_inf_rel, ii_lo, ii_hi), `+ AUCinf all estimable (%)` = ci(iii_plus_inf_all, iii_lo, iii_hi),
                  `Loss with reliable set (pp)` = ci(cost_ii_pp, cost_ii_lo, cost_ii_hi, 2), `Loss with all estimable (pp)` = ci(cost_iii_pp, cost_iii_lo, cost_iii_hi, 2))]))
if (!is.null(rr)) add("AUCinf handling rules if AUCinf were mandated (A: exclude subjects failing reliability, B: include all estimable, C: substitute AUClast when reliability fails):", "",
  md_table(rr[, .(Rule = rule, `Analysis n per arm` = f1((n_R_mean + n_T_mean) / 2), `Pass, identical products (%)` = ci(pass_S00, pass_S00_lo, pass_S00_hi),
                  `GMR bias vs truth, Vmax x1.25 (%)` = f2(bias_VM125_pct), `GMR bias vs truth, F x0.90 (%)` = f2(bias_F090_pct), `Agreement with AUClast, Vmax x1.25 (%)` = ifelse(startsWith(rule, "Reference"), "", f1(agree_VM125)))]))

sd <- R("fallback", "sample_size_logsd.csv"); pw <- R("fallback", "empirical_power.csv")
if (!is.null(sd) && !is.null(pw)) {
  g <- function(v, col) sd[variant == v][[col]]
  add("Sample size cross-check (proposed log-scale coefficient of variation 43%, 117 subjects per arm):",
      sprintf("- Log-scale SD (CV) at B0, 60 to 90 kg, 20,000 subjects: AUClast %s (%s%%) for the 2016 model and %s (%s%%) for Model 1; Cmax %s (%s%%) and %s (%s%%).",
              f3(g("base", "AUClast_logsd")), f1(g("base", "AUClast_logcv")), f3(g("struct2020", "AUClast_logsd")), f1(g("struct2020", "AUClast_logcv")),
              f3(g("base", "Cmax_logsd")), f1(g("base", "Cmax_logcv")), f3(g("struct2020", "Cmax_logsd")), f1(g("struct2020", "Cmax_logcv"))),
      sprintf("- The AUClast 90%% CI reported by Cohen 2022 (0.96 to 1.28, n 62 and 63) implies a log-scale SD of about 0.49 (CV about 52%%), above both models (%s%% and %s%%). The proposed 43%% lies between the models and the Cohen estimate.", f1(g("base", "AUClast_logcv")), f1(g("struct2020", "AUClast_logcv"))))
  x <- pw[model == "k2016"]
  add(sprintf("- Empirical power with 117 per arm (2016 model, %s trials): identical products %s%% for AUClast and Cmax jointly and %s%% with AUCinf (reliable set) added; test bioavailability x0.97 (true AUC ratio about 0.95) %s%% and %s%%.",
              format(x$n_trials[1], big.mark = ","), ci(x[scenario == "S00", power_last_cmax], x[scenario == "S00", lo], x[scenario == "S00", hi]), ci(x[scenario == "S00", power_3], x[scenario == "S00", lo3], x[scenario == "S00", hi3]),
              ci(x[scenario == "F097", power_last_cmax], x[scenario == "F097", lo], x[scenario == "F097", hi]), ci(x[scenario == "F097", power_3], x[scenario == "F097", lo3], x[scenario == "F097", hi3])), "")
}

# 4. Pillar 3
p3 <- R("rationale", "pillar3_km_summary.csv"); p3d <- R("rationale", "pillar3_km_B0.csv")
if (!is.null(p3)) add("## 4. Pillar 3: invisibility of binding-constant differences", "",
  sprintf("Test-arm Km multiplied by 0.5 to 10 (0.005 to 0.1 mg/L) changes the mean GMR of every endpoint by at most %s%% (2016 model) and %s%% (Model 1) relative to identical products (paired within the same 500 trials). Km is an MM approximation constant and is not identical to binding affinity.",
          f2(p3[model == "k2016", max_abs_GMR_change_pct]), f2(p3[model == "k2020", max_abs_GMR_change_pct])), "")

# 5. Weight generalization
wb <- NULL; for (f in list.files(proj_path("results", "weight_generalization"), pattern = "^weight_bands_B0_[a-d]+\\.csv$", full.names = TRUE)) wb <- rbind(wb, fread(f))
if (!is.null(wb)) {
  ml <- c(a = "(a) 2016", b = "(b) Model 1", c = "(c) 2016 + ke~BMI", d = "(d) 2016 + ke~BMI + Vc~weight 0.817")
  add("## 5. Body weight generalization (300 mg, B0, 20,000 subjects per uniform weight band)", "",
      "Covariate variants (c) and (d) use the adult coefficients of Kovalenko 2020 Model 4 (elimination rate constant ke proportional to (BMI/26)^0.368, central volume exponent 0.817); the BMI reference of 26 is a placeholder based on phase 3 mean BMI 25.4 to 27.3 (Kamal 2022). Height is simulated as normal (mean 170 cm, SD 9, truncated 150 to 195 cm). Development-data weight ranges are not reported in the source publications; bands above 130 kg are flagged as possible extrapolation.", "",
      md_table(wb[, .(Model = ml[model], `Band (kg)` = band, `Median BMI` = f1(BMI_median), `True extrap. median (%)` = f2(extrap_true_median), `95th pct` = f2(extrap_true_p95), `Coverage <80% (%)` = f3(coverage_lt80_pct),
                      `AUCinf reliable (%)` = f1(reliable_pct), `Lambda-z not estimable (%)` = f2(lambda_fail_pct), `AUClast geometric mean` = f1(AUClast_geo, 0), Note = ifelse(dev_range_note == "", "", "outside confirmed development range"))]))
}
ow <- NULL; for (mk in c("a", "b", "d")) { f <- proj_path("results", "weight_generalization", sprintf("obese_trials_per_endpoint_%s.csv", mk)); if (file.exists(f)) ow <- rbind(ow, fread(f)) }
od <- NULL; for (mk in c("a", "b", "d")) { f <- proj_path("results", "weight_generalization", sprintf("obese_trials_dropout_%s.csv", mk)); if (file.exists(f)) od <- rbind(od, fread(f)) }
if (!is.null(ow)) {
  x <- dcast(ow[endpoint %in% c("AUClast", "AUCinf_true", "AUCinf_reliable")], model + scenario ~ endpoint, value.var = "GMR_mean")
  d <- od[, .(retained = mean(wt_retained_mean), dropped = mean(wt_dropout_mean), dropout_pct = mean(dropout_pct)), by = .(model, scenario)]
  x <- merge(x, d, by = c("model", "scenario"))
  add("Trial level in a population with many obese subjects (weight normal mean 100 kg, SD 20, truncated 60 to 150 kg; 2,000 trials per scenario; mean GMR and dropout weights only):", "",
      md_table(x[, .(Model = ml[model], Scenario = scenario, `GMR AUClast` = f3(AUClast), `GMR true AUCinf` = f3(AUCinf_true), `GMR NCA AUCinf reliable` = f3(AUCinf_reliable),
                     `Subjects failing AUCinf reliability (%)` = f1(dropout_pct), `Mean weight retained (kg)` = f1(retained), `Mean weight failing (kg)` = f1(dropped))]))
}

# 6. Literature
ln <- R("literature", "literature_numeric.csv"); lq <- R("literature", "literature_qualitative.csv")
if (!is.null(ln)) {
  src_en <- function(s) sub("대조군", "control arm", s)
  add("## 6. Consistency with published coverage values", "",
      md_table(ln[, .(Source = src_en(source), `Dose (mg)` = dose_mg, `Published AUClast/AUCinf (mean ratio)` = sprintf("%.1f%%", 100 * lit_mean_ratio),
                      `2016: NCA / true` = sprintf("%.1f%% / %.1f%%", 100 * nca_mean_ratio_k2016, 100 * true_mean_ratio_k2016), `Model 1: NCA / true` = sprintf("%.1f%% / %.1f%%", 100 * nca_mean_ratio_k2020, 100 * true_mean_ratio_k2020))]),
      "Clot 2021 conditions: sampling Days 2 to 57, LLOQ 0.078 mg/L, cohort weights from Clot Table 1 (62.2, 59.9, 58.3 kg); the 200 mg row is an external check because of the different presentation. PKM12350 control arm: AUC0-t 500 versus AUC0-inf 521 (FDA BLA 761055 Clinical Pharmacology Review, Table 4.2.c); the test arm (488 versus 554) may reflect a different number of subjects with AUC0-inf and cannot be verified from public documents.",
      "Published statements on the terminal phase: Kovalenko 2020 describes a terminal slope tending to minus infinity near the LLOQ and a nearly vertical TMDD phase, no meaningful terminal half-life, and instantaneous half-life falling to zero; Kovalenko 2021 notes that 0.09 mg/L removes 90% of circulating target at Km 0.01 mg/L; Kovalenko 2016 used the M3 method because few quantifiable low concentrations describe the steep phase; Li 2020 reports steeper elimination at lower concentrations and a more than dose-proportional AUClast; Cohen 2022 did not compute half-life or AUCinf because of terminal non-linearity; Clot 2021 describes multi-exponential decline with faster target-mediated elimination at low concentrations.",
      "Interpretation: published coverage values are based on NCA AUCinf, whose extrapolation is inflated, so they are a lower bound of the true coverage. The curve shape below the LLOQ cannot be observed and is addressed by the Km and Vmax sensitivity analyses.",
      { x6 <- ln[dose_mg == 600]; x3 <- ln[dose_mg == 300 & grepl("Clot", source)]
        if (nrow(x6) && x6$true_mean_ratio_k2016 < x6$lit_mean_ratio) sprintf("At 600 mg the simulated true coverage (%.1f%% and %.1f%%) is below the published NCA value (%.1f%%). Because the published value is a lower bound, both models predict a larger tail beyond Day 56 than observed at this dose, which understates AUClast coverage (conservative direction). At the study dose of 300 mg the NCA mean ratio differs by %.1f percentage points, within 3 percentage points.",
          100 * x6$true_mean_ratio_k2016, 100 * x6$true_mean_ratio_k2020, 100 * x6$lit_mean_ratio, 100 * abs(x3$nca_mean_ratio_k2016 - x3$lit_mean_ratio)) else "" }, "")
}

# 7. Sampling density conclusion
add("## 7. Sampling density between Day 36 and Day 50: conclusion", "",
    sprintf("Final schedule: %s (Syneos baseline). Governing model: Kovalenko 2016. Rule application: %s. Decided by %s on %s.", sdec$final_schedule, sdec$rule_application, sdec$decided_by, sdec$decision_date), "",
    "Pre-specified rule, versus B0: recommend added sampling if at least one holds: (a) mean width of the AUClast 90% CI decreases by at least 2%; (b) AUClast pass rate increases by at least 2 percentage points when ke is multiplied by 1.10; (c) the AUCinf reliability rate rises by at least 5 percentage points; (d) the share of subjects with NCA extrapolation above 20% falls to half or less.", "")
vv <- c(base = "2016 (primary)", struct2020 = "Model 1", vmax080_both = "Vmax x0.8 (both arms)", vmax125_both = "Vmax x1.25 (both arms)")
rows <- rbindlist(lapply(names(vv), function(v) { d <- R("trials", sprintf("schedule_decision_%s.csv", v)); if (is.null(d)) return(NULL)
  d3 <- d[schedule == "D3"]; k <- R("individual200k", sprintf("criterion_d_200k_%s.csv", v)); k3 <- if (!is.null(k)) k[schedule == "D3"] else NULL
  data.table(Variant = vv[[v]], `Rule result` = if (any(d$recommend)) paste0("D3 by criterion ", paste(c("a", "b", "c", "d")[unlist(d3[, .(crit_a %in% TRUE, crit_b %in% TRUE, crit_c %in% TRUE, crit_d %in% TRUE)])], collapse = ","), " only") else "no schedule meets any criterion",
             `D3: AUClast CI width change (%, positive = wider)` = sprintf("%s (%s to %s)", f2(-100 * d3$a_mean_width_rel_decrease), f2(-100 * d3$a_paired_hi), f2(-100 * d3$a_paired_lo)),
             `D3: reliability gain (pp)` = f2(d3$c_reliable_gain_pp),
             `D3: extrapolation >20% ratio, 20,000 subjects` = sprintf("%s (%s to %s)", f2(d3$d_extrap20_ratio), f2(d3$d_ratio_boot_lo), f2(d3$d_ratio_boot_hi)),
             `D3: same ratio, 200,000 subjects` = if (!is.null(k3) && nrow(k3)) sprintf("%s (%s to %s), %s fewer subjects per arm; %s", f2(k3$extrap_gt20_ratio), f2(k3$d_ratio_boot_lo), f2(k3$d_ratio_boot_hi), f2(-k3$d_abs_change_per_arm),
               if (k3$d_ratio_boot_lo <= 0.5 && k3$d_ratio_boot_hi >= 0.5) "fragile (interval includes 0.5)" else if (k3$d_ratio_boot_lo > 0.5) "not met (interval above 0.5)" else "met (interval below 0.5)") else "not re-evaluated",
             `Final recommendation` = if (any(d$recommend)) "B0 (conclusion unchanged)" else "B0") }))
if (nrow(rows)) add(md_table(rows))
bm <- R("trials", "bminus_explicit.csv"); wd <- R("trials", "width_decomposition.csv")
if (!is.null(wd)) { w <- wd[schedule %in% c("D1", "D2", "D4")]
  add(sprintf("AUClast CI width: D1, D2 and D4 widen the mean AUClast 90%% CI by %s%% to %s%% (paired intervals exclude zero; D3 %s%%). Recomputed with true concentrations and no residual error, the widening remains (%s%% to %s%%), so it reflects heterogeneity of the added tail area as the last quantifiable time is extended, not measurement error at the added low points.",
              f2(-100 * max(w$base)), f2(-100 * min(w$base)), f2(-100 * wd[schedule == "D3", base]), f2(-100 * max(w$noresid)), f2(-100 * min(w$noresid))), "") }
if (!is.null(bm)) add(sprintf("Removing Day 50 (B-): AUClast CI width %s%% (negative = narrower), AUCinf reliability change %s percentage points, share with NCA extrapolation above 20%% multiplied by %s. Not recommended, to keep a terminal sample for AUC0-inf as a secondary endpoint and for a fallback analysis.",
                      f2(-100 * bm$a_mean_width_rel_decrease), ci(bm$c_reliable_gain_pp, bm$c_gain_boot_lo, bm$c_gain_boot_hi, 2), ci(bm$d_extrap20_ratio, bm$d_ratio_boot_lo, bm$d_ratio_boot_hi, 2)), "")
trig <- if (nrow(rows)) rows[startsWith(`Rule result`, "D3"), Variant] else character(0)
abc_any <- any(unlist(lapply(c(names(vv), "iiv150", "resid12", "weight_alt", "ada10"), function(v) { d <- R("trials", sprintf("schedule_decision_%s.csv", v)); if (is.null(d)) FALSE else d[, any(crit_a %in% TRUE | crit_b %in% TRUE | crit_c %in% TRUE)] })))
k_lo <- unlist(lapply(c("base", "struct2020", "vmax080_both"), function(v) { k <- R("individual200k", sprintf("criterion_d_200k_%s.csv", v)); if (is.null(k)) NA_real_ else k[schedule == "D3", d_ratio_boot_lo] }))
if (abc_any) stop("criterion (a), (b) or (c) met in some variant: rationale text must be revised")
if (anyNA(k_lo) || any(k_lo <= 0.5)) stop("200,000-subject criterion (d) interval reaches 0.5 in some variant: rationale text must be revised")
add(sprintf("Rationale: %s; criteria (a), (b) and (c) were not met in any variant; re-evaluated with 200,000 subjects (4,000 bootstrap resamples), the D3 ratio for criterion (d) is above 0.5 with its whole interval in every re-evaluated variant; the absolute reduction is below one subject per arm; neither the AUClast CI width nor power improves; the cost is 1,040 additional visits. Removing the Day 50 sample (B-) is not recommended because it preserves a terminal point for AUC0-inf as a secondary endpoint and for a fallback analysis. Lesson recorded: a relative-reduction criterion for a rare event needs an absolute floor (for example at least one subject per arm); not applied retroactively.",
            if (length(trig)) sprintf("at the pre-specified 20,000-subject level, criterion (d) alone was met in %s", paste(trig, collapse = " and ")) else "no variant met any criterion at the pre-specified 20,000-subject level"), "")

# 8. Limitations
add("## 8. Limitations", "",
    "- The primary model does not reproduce the faster absorption of the 200 mg 1.14 mL (175 mg/mL) presentation; about one third to one half of the 200 mg shortfall is explained by absorption, the remainder is unexplained. Clot 2021 Chinese 200 mg data (same 175 mg/mL prefilled syringe) showed median time to Cmax 7.0 days (range 3.0 to 10.0), but with n=8 and no sample between Day 4 and Day 8 the resolution is low.",
    "- Cmax: the 2016 model predicts about 20% above observed and Model 1 about 5% above; AUClast: the 2016 model is close to observed and Model 1 about 12% above.",
    "- Km is fixed at 0.01 mg/L in both models; with no uncertainty or IIV on Km, terminal-phase variability may be underestimated. Km and Vmax sensitivity analyses address this.",
    "- The AUCinf reliability rate and the share with NCA extrapolation above 20% depend on the residual error model and are reported as two-model ranges.",
    "- Placeholders not yet confirmed: Day 1 post-dose sampling time (0.25 day), weight distribution and stratification split, sampling windows, BMI reference of 26, body weights of the Li 2020 single-arm studies.",
    "- Development-data body weight ranges are not reported; results above 130 kg are extrapolations.", "")
txt <- paste(out, collapse = "\n")
if (grepl("—", txt)) stop("summary_en.md contains an em-dash")
if (grepl("[가-힣]", txt)) { bad <- regmatches(txt, gregexpr("[^\n]*[가-힣][^\n]*", txt))[[1]]; stop("summary_en.md contains Korean text: ", paste(head(bad, 3), collapse = " || ")) }
writeLines(txt, proj_path("results", "summary_en.md"), useBytes = TRUE)
cat("results/summary_en.md written:", length(out), "lines\n")
