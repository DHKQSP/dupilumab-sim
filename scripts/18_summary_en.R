#!/usr/bin/env Rscript
# results/summary_en.md: English summary for regulatory briefing (review consolidated 2026-09-24, section 9).
# Built from the same result files as the report. Rules: English only, no em-dash, abbreviations spelled out at first use,
# trial-level proportions cited only from runs with >= 5,000 trials (with 95% intervals); two-model ranges for residual-sensitive metrics.
source("R/00_setup.R"); source_project()
R <- function(...) { f <- proj_path("results", ...); if (file.exists(f) || file.exists(paste0(f, ".gz"))) read_raw(f) else NULL }
# 반올림은 fmt_num과 같은 round_half_away(0에서 먼 쪽): 20,000명 비율(0.005의 배수)이 표와 본문(reliability_facts, mechanism_summary)에서 같게 찍히도록
f1 <- function(x, d = 1) formatC(round_half_away(x, d), format = "f", digits = d); f2 <- function(x) f1(x, 2); f3 <- function(x) f1(x, 3)
ci <- function(e, lo, hi, d = 1) sprintf("%s (%s to %s)", f1(e, d), f1(lo, d), f1(hi, d))
md_table <- function(dt) { h <- paste0("| ", paste(names(dt), collapse = " | "), " |"); s <- paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|")
  b <- apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")); c(h, s, b, "") }
out <- c()
add <- function(...) out <<- c(out, ...)
# 그림: results/ 기준 상대 경로의 영문판(_en.png). 캡션에 모델·조건·반복 수
fig <- function(rel, caption) { if (file.exists(proj_path("results", rel))) add(sprintf("![%s](%s)", caption, rel), "", sprintf("*%s*", caption), "") }
oc_cfg <- read_cfg("oc_design.yaml"); nb_ <- format(oc_cfg$trials$reps_boundary, big.mark = ","); no_ <- format(oc_cfg$trials$reps_other, big.mark = ",")
# OC 그림·표의 시험 수: 결과 파일의 실제 수(적응적 연장 포함, scripts/42와 R/oc.R oc_n_trials_text). 결과 파일이 없으면 설계값
oc_cv_ <- R("oc", "oc_curves.csv"); oc_dec_ <- R("oc", "extension_decision.csv"); bnd_ <- as.numeric(unlist(oc_cfg$boundary_targets))
oc_ntxt <- function(m_, part) {                                     # part: "bnd" 경계 + 동일 제품, "bnd_only" 경계만, "other" 나머지
  if (is.null(oc_cv_)) return(sprintf("%s each", if (part == "other") no_ else nb_))
  x <- oc_cv_[oc_cv_$model == m_ & oc_cv_$config == "P2"]; isb <- abs(x$target - bnd_[1]) < 1e-9 | abs(x$target - bnd_[2]) < 1e-9
  x <- switch(part, bnd = x[isb | x$scenario == "S00"], bnd_only = x[isb], other = x[!isb & x$scenario != "S00"])
  if (nrow(x)) oc_n_trials_text(x, "en", oc_dec_) else "none"
}
MLF <- c(k2016 = "Kovalenko 2016 (primary)", k2020 = "Kovalenko 2020 Model 1")
gate <- read_cfg("gate_decision.yaml"); sdec <- read_cfg("schedule_decision.yaml"); design <- read_cfg("trial_design.yaml")
# 신뢰 충족률·탈락률과 그 증감: 결과 파일에서 만들고 전제를 검사한다(reliability_facts, 어긋나면 중단). 플래그 세트 (i) 주값, (ii) 대괄호 병기
RF <- reliability_facts()
FPe <- function(i, ii, d = 1, unit = "%", signed = FALSE, each = FALSE) fmt_flag_pair(i, ii, d, unit, " to ", signed, each)
cie <- function(e, lo, hi, d = 2) sprintf("%s (%s to %s)", fmt_num(e, d), fmt_num(lo, d), fmt_num(hi, d))   # 신뢰 충족 증감: 39 문구와 같은 반올림
premise <- function(ok, msg) if (!isTRUE(all(ok))) stop("summary_en.md premise not met: ", msg, call. = FALSE)
REL_DEF_EN <- "flag set (i): lambda-z estimable, adjusted R-squared at least 0.80 and extrapolated AUC at most 20%; flag set (ii): (i) plus span ratio at least 2"
REL_CONV_EN <- "A span ratio of at least 2 is a convention used only by some statistical analysis plans (not a Phoenix feature, D-039)"
ML2e <- c(k2016 = "2016 model", k2020 = "Model 1")
rel_changes_en <- if (!is.null(RF)) sprintf("Numbers that change because of the span ratio criterion: reliability at B0, %s (loss from the span ratio alone %s to %s percentage points); the change in reliability with added sampling (the criterion (c) quantity; D3: %s; sign reversed in %s); the analysis set and pass rates of configurations that use rules A or C (G2, F3-A, F3-C); dropout rates and the weight of dropouts. Unchanged: criterion (d) (extrapolation above 20%%), AUClast, Cmax, rule B and the schedule recommendation (no cell where the two flag sets recommend differently).",
    FPe(RF$rng$rel_i, RF$rng$rel_ii), fmt_num(RF$rng$span_loss[1]), fmt_num(RF$rng$span_loss[2]),
    paste(sprintf("%s %s", ML2e[RF$dense$d3$model], FPe(RF$dense$d3$gain_i_pp, RF$dense$d3$gain_ii_pp, 2, " points", TRUE, TRUE)), collapse = ", "),
    if (nrow(RF$flips)) paste(sprintf("%s %s", ML2e[RF$flips$model], RF$flips$schedule), collapse = ", ") else "no cell") else ""
# 결과 .md를 이 문서에 넣는다: 첫 줄 제목(# ...)을 빼고 제목 수준을 shift만큼 낮춘다
md_body <- function(path, shift = 1) { if (!file.exists(path)) return(character(0)); x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(x) && startsWith(x[1], "# ")) x <- x[-1]; c(sub("^(#+) ", paste0("\\1", strrep("#", shift), " "), x), "") }

add("# Dupilumab biosimilar Phase 1 pharmacokinetic simulation: summary for regulatory briefing", "",
    sprintf("Generated %s from repository results (branch claude/epic-bardeen-axbreo). Every number below is read from the result files used by the full report.", format(Sys.Date())), "",
    paste0("Abbreviations: area under the concentration-time curve to the last quantifiable concentration (AUClast), to infinity (AUCinf); maximum concentration (Cmax); non-compartmental analysis (NCA); geometric mean ratio (GMR); confidence interval (CI); lower limit of quantification (LLOQ, " , format(study_lloq()), " mg/L); inter-individual variability (IIV); target-mediated drug disposition (TMDD); Michaelis-Menten (MM); Monte Carlo (MC); body mass index (BMI); operating characteristic (OC); terminal elimination rate constant (lambda-z); common random numbers (CRN).", ""),
    "Models: primary model Kovalenko et al. 2016 (CPT Pharmacometrics Syst Pharmacol 5:617, Table 2, BLQ-included column): two-compartment, first-order absorption, parallel linear and MM elimination, Km fixed at 0.01 mg/L, central volume scaled by (weight/75)^0.705. Sensitivity model Kovalenko et al. 2020 Model 1 (Clin Pharmacol Drug Dev 9:756, Table 1 and Supplementary Table 2): transit absorption (3 compartments, mean transit time 0.105 day), its own IIV and residual error (proportional 15.0%, additive 0.03 mg/L).", "",
    "Scenario codes (test arm only unless stated; reference arm shared through common random numbers): S00 identical products; F085, F090, F097, F110 bioavailability x0.85, x0.90, x0.97, x1.10; KE110, KE120 linear elimination rate constant (ke) x1.10, x1.20; VM080, VM125, VM150 maximum MM elimination rate (Vmax) x0.80, x1.25, x1.50; KM05, KM2, KM5, KM10 MM constant (Km) x0.5, x2, x5, x10; KA075 absorption rate constant x0.75. Sampling schedules: B0 Syneos baseline; D1 to D4 add two to four samples between Day 32 and Day 53 (D1: Days 39, 46; D2: Days 39, 46, 53; D3: Days 32, 39, 46, 53; D4: Days 40, 47); B- removes Day 50.", "",
    "Study design simulated: 300 mg single subcutaneous dose (2 mL of 150 mg/mL), parallel groups, 117 evaluable subjects per arm, body weight 60 to 90 kg, Syneos sampling schedule (B0: Days 1, 2, 4, 6, 8, 11, 15, 22, 29, 36, 43, 50, 57). Equivalence: two one-sided tests via the 90% CI of the GMR from a pooled two-sample t on log scale, limits 80.00% to 125.00%.", "")

# 0. Key conclusions: operating characteristics (pre-specified design) and the sampling cliff
ocen <- proj_path("results", "oc", "oc_conclusion_en.md"); clen <- proj_path("results", "cliff", "cliff_conclusion_en.md"); pr <- R("oc", "prereg.csv")
if (file.exists(ocen)) add("## Key conclusion: operating characteristics of the co-primary endpoint configurations", "",
  sprintf("Design fixed before any result in config/oc_design.yaml (commit %s%s). Truth is the population GMR of model-integrated AUC0-inf (no residual error; the same 200,000 virtual subjects receive both products, CRN), not AUClast; the true Cmax ratio is reported alongside. Configurations: P2 = AUClast + Cmax (proposed); F3-A, F3-B, F3-C = P2 + AUCinf under handling rule A, B or C; G2 = AUCinf (rule A) + Cmax (guideline default). Mechanisms F, ka, ke, Vmax, Km and peripheral volume V2 in both directions, with multipliers inverted by bisection to target true AUC0-inf ratios 0.70 to 1.43.",
          if (!is.null(pr)) substr(pr$prereg_commit, 1, 7) else "not found", if (!is.null(pr) && isTRUE(pr$changed_since)) ", changed since" else ", unchanged since"), "",
  c(rbind(readLines(ocen, warn = FALSE), "")))   # 결론 문구 한 줄 = 한 문단(첫 문단 = P2 경계 1종 오류)
# 핵심 수치 요약표(지시 2026-09-25 §4): 1차 지표(경계 1종 오류 최대), 예비 제품 시나리오 VM150·KE120, Pillar 3(Km). 수치·전제 검사는 key_facts(R/summarize.R)
KF <- key_facts(oc = oc_cfg); ktab_ <- key_summary_table(KF, "en")
if (nrow(ktab_)) add("## Key numbers", "",
  sprintf("The primary metric is the boundary type I error (pass rate at a true AUC0-inf ratio of %s; mechanism by direction by configuration by model). The consumer and producer risks of the random product space are a secondary metric that depends on the assumed virtual product distribution. VM150 and KE120 are preliminary scenarios with arbitrary multipliers (primary model, at least 5,000 trials); their AUCinf reliability uses flag set (ii), and their true ratio is the mean trial GMR of the individual model AUC0-inf of the same trial subjects, a different definition from the 200,000-subject integral of the inverted scenarios. Every value below is read from the named result file (proportions with Wilson 95%% intervals).",
          paste(fmt_num(KF$bnd, 2), collapse = " or ")), "",
  md_table(ktab_))
if (file.exists(clen)) { add("## Key conclusion: sampling on the terminal cliff", "", readLines(clen, warn = FALSE), "")
  ccf <- oc_cfg$cliff; nsub <- format(ccf$n_subjects, big.mark = ",")
  fig("cliff/fig2_1_lloq_day_en.png", sprintf("Figure 2-1. Study day when the true concentration reaches the LLOQ (Kovalenko 2016 primary model, 60 to 90 kg, 300 mg, %s virtual subjects, nominal days; bins containing a current sampling day highlighted).", nsub))
  fig("cliff/fig2_2_points_in_cliff_en.png", sprintf("Figure 2-2. Share of subjects with 1, 2 or 3 or more samples on the cliff, by schedule, cliff definition (instantaneous half-life below 1 or 2 days) and timing (nominal days or visit windows of plus or minus 1 day); Kovalenko 2016, 60 to 90 kg, %s virtual subjects.", nsub))
  fig("cliff/fig2_3_cliff_length_en.png", sprintf("Figure 2-3. Cliff length under the 1-day and 2-day definitions against the minimum visit interval after Day 22 of each schedule (Kovalenko 2016, 60 to 90 kg, %s virtual subjects).", nsub))
  fig("cliff/fig2_4_representative_en.png", "Figure 2-4. True concentration of three representative subjects (25th, 50th and 75th percentile of the day of reaching the LLOQ) with the cliff shaded (1-day definition), current samples and the added Day 39, 46, 53 samples (Kovalenko 2016, 60 to 90 kg).") }

# 0b. NCA engine
ev <- R("nca_engine", "engine_validation_summary.csv"); dr <- R("nca_engine", "dropout_reasons_B0.csv"); edif <- R("nca_engine", "engine_difference_individual_B0.csv")
if (!is.null(ev)) {
  add("## NCA engine: Phoenix WinNonlin-compatible rules", "",
      "Rules: BLQ handled explicitly before NCA (zero before the first quantifiable value, missing between quantifiable values, excluded after the last one, and quantifiable values after two consecutive BLQ set to missing); linear-up log-down AUC; lambda-z Best Fit (last 3, 4, 5, ... positive concentrations after Cmax, windows with positive slope excluded, largest adjusted R-squared, ties within 0.0001 resolved to more points). Reliability flags (statistical analysis plan conventions, not Phoenix features): adjusted R-squared below 0.80, extrapolated AUC above 20% and, in flag set (ii) only, span ratio below 2. AUCinf handling rules: A excludes flagged subjects, B includes all subjects with lambda-z, C substitutes AUClast for flagged subjects. Reliability and dropout rates are cited for flag set (i) with flag set (ii) in brackets.",
      sprintf("Validation against the reference implementation NonCompart 0.8.4 and PKNCA 0.12.1 on Theoph (12 profiles), Indometh (6) and %s simulated dupilumab profiles: lambda-z points identical in every profile and every pair of implementations; largest relative difference in any parameter %s (criterion 1e-6).",
              format(max(ev$n_profiles), big.mark = ","), formatC(max(ev$max_rel_diff), format = "e", digits = 1)), "")
  if (!is.null(RF)) { x <- RF$b0; rg <- function(v) sprintf("%s%% to %s%%", fmt_num(min(v)), fmt_num(max(v)))
    add(sprintf("Subjects failing the reliability criteria at B0 (20,000 per model, two-model ranges): %s. By reason (flags overlap): lambda-z not estimable %s, adjusted R-squared below 0.80 %s, extrapolation above 20%% %s, span ratio below 2 %s (counted in flag set (ii) only; the span ratio alone removes %s). %s.",
                FPe(RF$rng$fail_i, RF$rng$fail_ii), rg(x$lambda_fail), rg(x$flag_rsq), rg(x$flag_extrap), rg(x$flag_span), rg(x$span_loss), REL_CONV_EN), "")
    if (!is.null(RF$eng)) add(sprintf("Previous engine (D-010) versus the new engine on the same observations (a separate B0-only draw, so the values differ slightly from the main analysis sample): reliability %s; the span ratio flag alone accounts for %s of the drop.",
                                      paste(sprintf("%s%% to %s (%s)", fmt_num(RF$eng$old), FPe(RF$eng$new_i, RF$eng$new_ii, each = TRUE), ML2e[RF$eng$model]), collapse = " and "), paste(sprintf("%s%%", fmt_num(RF$eng$span_share_pct)), collapse = " and ")), "")
  } else if (!is.null(dr)) { x <- dr[model %in% c("base", "struct2020")]
    add(sprintf("Subjects failing the reliability criteria at B0, flag set (ii) (20,000 per model, separate B0-only draw of the engine comparison, flags overlap): %s%% to %s%% in total; span ratio below 2 %s%% to %s%%. %s.",
                f1(min(x$any_flag_or_fail_pct)), f1(max(x$any_flag_or_fail_pct)), f1(min(x$flag_span_pct)), f1(max(x$flag_span_pct)), REL_CONV_EN), "") }
}

# 0c. 신뢰 플래그 두 세트와 촘촘한 후기 채혈(검토 의견 W2 §1): scripts/39 결과
rbs <- R("reliability", "reliability_by_schedule.csv")
if (!is.null(RF) && !is.null(rbs)) {
  add("## AUCinf reliability: two flag sets and denser late sampling (review W2, section 1)", "",
      sprintf("Definitions: %s. %s. %s", REL_DEF_EN, REL_CONV_EN, rel_changes_en), "",
      sprintf("Structure: after Day 22 the shortest 3-point lambda-z window allowed by the nominal days is %s days at B0 and %s days with added sampling (D1 to D4), so a 3-point window reaches span ratio 2 only if the half-life is at most %s and %s days. Saved subject-level NCA of 20,000 virtual subjects per model under B0 conditions; no new simulation (scripts/39_reliability_flags.R).",
              RF$lz_min[["B0"]], RF$lz_min[["dense"]], fmt_num(RF$lz_min[["B0"]] / 2), fmt_num(RF$lz_min[["dense"]] / 2)), "")
  x <- merge(rbs[variant %in% c("base", "struct2020")], RF$paired[, .(variant, schedule, gain_i_pp, gain_ii_pp)], by = c("variant", "schedule"), all.x = TRUE)
  x <- x[order(match(variant, c("base", "struct2020")), match(schedule, c("B0", "Bminus", "D1", "D2", "D3", "D4")))]
  add(md_table(x[, .(Model = ML2e[c(base = "k2016", struct2020 = "k2020")[variant]], Schedule = schedule, Samples = n_points,
                     `Reliable, flag set (i) (%, 95% CI)` = ci(reliable_i_pct, reliable_i_lo, reliable_i_hi, 2), `Reliable, flag set (ii) (%, 95% CI)` = ci(reliable_ii_pct, reliable_ii_lo, reliable_ii_hi, 2),
                     `Change vs B0, (i) [(ii)] (pp)` = ifelse(is.na(gain_i_pp), "", FPe(gain_i_pp, gain_ii_pp, 2, "", TRUE, TRUE)), `Loss from span ratio alone (pp)` = f2(span_only_loss_pct))]))
  fig("reliability/fig_reliability_mechanism_en.png", "Figure R-1. Why denser late sampling adds span flags: lambda-z window length and span ratio by schedule (top) and, for subjects reliable at B0 and span-flagged with added sampling, the change in median window length and half-life (bottom). Kovalenko 2016 and Model 1, 20,000 virtual subjects per model, B0 conditions (60 to 90 kg, visit windows, residual error, BLQ).")
  add(md_body(proj_path("results", "reliability", "mechanism_summary_en.md"), 1))
}

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
  s1_ <- sh[ka_multiplier == 1, ratio_per_mg]; sk_ <- sh[ka_multiplier %in% c(1.5, 2), ratio_per_mg]; closed_ <- 100 * (sk_ - s1_) / (so[1, ratio_per_mg] - s1_); remain_ <- 100 * (1 - sk_ / so[1, ratio_per_mg])
  premise(length(sk_) == 2 && all(closed_ >= 25 & closed_ <= 60), "absorption matching closes about one third to one half of the 200 mg shortfall")
  add("200 mg external checks and absorption (diagnostic only, not adopted): the 2016 model does not reproduce the fast absorption of the 200 mg 1.14 mL (175 mg/mL) presentation (observed median time to Cmax 3 days versus 7 days simulated).",
      sprintf("Matching absorption to the observed time to Cmax closes about one third to one half (%s%% to %s%%) of the 200 mg exposure shortfall: the dose-normalized 200:300 mg AUClast ratio is %s observed (four-arm, n-weighted; %s for the PKM14271 test arm alone) versus %s simulated, and %s to %s with the absorption rate constant multiplied by 1.5 to 2 for 200 mg only. The remainder (%s%% to %s%% on the four-arm basis) is unexplained: a presentation-specific bioavailability difference or over-estimated TMDD at low exposure are possible; the latter is covered by the population maximum elimination rate (Vmax) x0.8 sensitivity variant. For the study presentation (300 mg, 2 mL of 150 mg/mL) the model reproduces both absorption timing and exposure.",
              fmt_num(min(closed_), 0), fmt_num(max(closed_), 0), f2(so[1, ratio_per_mg]), f2(so[2, ratio_per_mg]), f2(sh[ka_multiplier == 1, ratio_per_mg]), f2(sh[ka_multiplier == 1.5, ratio_per_mg]), f2(sh[ka_multiplier == 2, ratio_per_mg]),
              fmt_num(min(remain_), 0), fmt_num(max(remain_), 0)), "")
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
  # (i)·(ii) 모두 reliability_facts 값에서 같은 반올림(f2)으로: pillar1 CSV의 미리 만든 문자열(%.2f)은 80.735 → 80.73으로 본문(80.74)과 달라진다.
  # (ii)는 pillar1과 같은 표본임을 reliability_facts가 0.01 이내로 검사한다
  if (!is.null(RF)) {
    rb <- RF$b0[match(p1e$model, RF$b0$model)]
    p1e[, `:=`(rel_i_ci = sprintf("%s [%s, %s]", f2(rb$rel_i), f2(rb$rel_i_lo), f2(rb$rel_i_hi)), reliable_pct_ci = sprintf("%s [%s, %s]", f2(rb$rel_ii), f2(rb$rel_ii_lo), f2(rb$rel_ii_hi)))]
  } else p1e[, rel_i_ci := "not available"]
  add(md_table(p1e[, .(Model = fifelse(model == "k2016", "2016 (primary)", "Model 1"), `True extrapolated %: median` = f2(extrap_true_median), `95th percentile` = f2(extrap_true_p95), Max = f2(extrap_true_max),
                       `Window coverage below 80% (%, 95% CI)` = coverage_lt80_pct_ci, `NCA extrapolated %: median` = f2(extrap_nca_median), `NCA extrapolation above 20% (%, 95% CI)` = extrap_nca_gt20_pct_ci,
                       `AUCinf reliability met, flag set (i) (%, 95% CI)` = rel_i_ci, `AUCinf reliability met, flag set (ii) (%, 95% CI)` = reliable_pct_ci)]))
  add(sprintf("Residual-sensitive metrics are stated as the two-model range: AUCinf reliability criteria are not met in %s of subjects (%s); NCA extrapolation above 20%% occurs in %s%% to %s%%. True extrapolation (model integral beyond the last quantifiable time) is essentially the same in both models.",
              if (!is.null(RF)) FPe(RF$rng$fail_i, RF$rng$fail_ii) else sprintf("[(ii) %s%% to %s%%]", f1(100 - max(p1e$reliable_pct)), f1(100 - min(p1e$reliable_pct))), REL_DEF_EN,
              f2(min(p1e$extrap_gt20_pct)), f2(max(p1e$extrap_gt20_pct))), "")
}
cs <- R("curve_shape", "curve_shape_B0.csv")
if (!is.null(cs)) {
  vl <- c(base = "Primary model", km05_both = "Km x0.5", km2_both = "Km x2", km5_both = "Km x5", km10_both = "Km x10", vmax080_both = "Vmax x0.8 (longer tail, conservative)", vmax125_both = "Vmax x1.25", vmax050_both = "Vmax x0.5 (stress test, fails gate)")
  cs[, v := vl[variant]]; cs <- cs[match(names(vl), variant)][!is.na(variant)]
  add("Curve shape below the LLOQ cannot be observed; sensitivity to Km and Vmax (both arms, 20,000 subjects each):", "",
      md_table(cs[, .(Variant = v, `True extrapolated %: median` = f2(extrap_true_median), `95th pct` = f2(extrap_true_p95), Max = f2(extrap_true_max), `Window coverage below 80% (%)` = f2(coverage_lt80_pct),
                      `NCA extrap. median (%)` = f2(extrap_nca_median), `NCA extrap. >20% (%)` = f2(extrap_gt20_pct), `AUCinf reliable, (i) [(ii)] (%)` = sprintf("%s [%s]", f1(reliable_rsq_extrap_pct), f1(reliable_pct)), `Lambda-z not estimable (%)` = f2(lambda_fail_pct),
                      `Median last quantifiable day` = f1(tlast_median), `AUClast geometric mean` = f1(AUClast_geo, 0), `300 mg ~78 kg AUClast vs observed 544` = f2(AUClast_ratio_vs_obs544))]),
      "Kovalenko 2020 reported that excluding below-LLOQ values makes the model predict a less steep TMDD phase with Vm and Km increasing together; this motivates the Km-increase sensitivity.", "")
}

# 3. Pillar 2
tpc <- proj_path("results", "trialpop", "tp_conclusion_en.md")
if (file.exists(tpc)) { x <- readLines(tpc, encoding = "UTF-8"); x <- x[-1]; x <- sub("^## ", "### ", x)
  add("## 2A. Reliability of AUCinf in the trial population (healthy adults, 60 to 90 kg, weight-stratified randomization)", "",
      "Window coverage is the true AUC from dosing to the last quantifiable sample divided by the true AUCinf; the observed-to-true ratio is the observed AUClast (or NCA AUCinf) divided by the true AUCinf and also contains the residual error and the trapezoidal approximation.", "", x, "",
      "Source: results/trialpop/ (scripts/54_trial_population_drop.R, scripts/55_trial_population_summary.R; registered in config/prereg_20260926.yaml section 6).", "") }
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
if (!is.null(cr) && nrow(cr)) add("Preliminary exploration only (arbitrary multipliers; the out-of-range judgement is superseded by the pre-specified operating-characteristic analysis below):", "",
  md_table(cr[, .(Scenario = scenario, `True ratio` = f3(true_ratio), Trials = n_trials, `Pass AUClast and Cmax (%)` = ci(joint_last_cmax, joint_last_cmax_lo, joint_last_cmax_hi, 2),
                  `Pass all three (%)` = ci(joint_3, joint_3_lo, joint_3_hi, 2), `Only AUCinf fails (%)` = ci(only_inf_fails, only_inf_fails_lo, only_inf_fails_hi, 2))]))
if (!is.null(fc)) add("Cost of a three-endpoint fallback (joint pass rates, >= 5,000 trials):", "",
  md_table(fc[, .(Scenario = scenario, `AUClast + Cmax (%)` = ci(i_last_cmax, i_lo, i_hi), `+ AUCinf reliable (%)` = ci(ii_plus_inf_rel, ii_lo, ii_hi), `+ AUCinf all estimable (%)` = ci(iii_plus_inf_all, iii_lo, iii_hi),
                  `Loss with reliable set (pp)` = ci(cost_ii_pp, cost_ii_lo, cost_ii_hi, 2), `Loss with all estimable (pp)` = ci(cost_iii_pp, cost_iii_lo, cost_iii_hi, 2))]))
if (!is.null(rr)) add("AUCinf handling rules if AUCinf were mandated (A: exclude subjects failing reliability, B: include all estimable, C: substitute AUClast when reliability fails; reliability under flag set (ii)):", "",
  md_table(rr[, .(Rule = rule, `Analysis n per arm` = f1((n_R_mean + n_T_mean) / 2), `Pass, identical products (%)` = ci(pass_S00, pass_S00_lo, pass_S00_hi),
                  `GMR bias vs truth, Vmax x1.25 (%)` = f2(bias_VM125_pct), `GMR bias vs truth, F x0.90 (%)` = f2(bias_F090_pct), `Agreement with AUClast, Vmax x1.25 (%)` = ifelse(startsWith(rule, "Reference"), "", f1(agree_VM125)))]))

sd <- R("fallback", "sample_size_logsd.csv"); pw <- R("fallback", "empirical_power.csv")
if (!is.null(sd) && !is.null(pw)) {
  g <- function(v, col) sd[variant == v][[col]]
  add("Sample size cross-check (proposed log-scale coefficient of variation 43%, 117 subjects per arm):",
      sprintf("- Log-scale SD (CV) at B0, 60 to 90 kg, 20,000 subjects: AUClast %s (%s%%) for the 2016 model and %s (%s%%) for Model 1; Cmax %s (%s%%) and %s (%s%%).",
              f3(g("base", "AUClast_logsd")), f1(g("base", "AUClast_logcv")), f3(g("struct2020", "AUClast_logsd")), f1(g("struct2020", "AUClast_logcv")),
              f3(g("base", "Cmax_logsd")), f1(g("base", "Cmax_logcv")), f3(g("struct2020", "Cmax_logsd")), f1(g("struct2020", "Cmax_logcv"))),
      { cv_m <- c(g("base", "AUClast_logcv"), g("struct2020", "AUClast_logcv")); cv_c <- sd[startsWith(variant, "Cohen"), AUClast_logcv]; cv_s <- sd[!startsWith(variant, "Cohen") & !variant %in% c("base", "struct2020"), AUClast_logcv]
        premise(length(cv_c) == 1 && length(cv_s) == 1 && max(cv_m) < cv_s && cv_s < cv_c, "the proposed CV lies between the models and the Cohen 2022 estimate")
        sprintf("- The AUClast 90%% CI reported by Cohen 2022 (0.96 to 1.28, n 62 and 63) implies a log-scale SD of about 0.49 (CV about %s%%), above both models (%s%% and %s%%). The proposed %s%% lies between the models and the Cohen estimate.", fmt_num(cv_c, 0), f1(cv_m[1]), f1(cv_m[2]), fmt_num(cv_s, 0)) })
  x <- pw[model == "k2016"]
  add(sprintf("- Empirical power with 117 per arm (2016 model, %s trials): identical products %s%% for AUClast and Cmax jointly and %s%% with AUCinf (reliable set) added; test bioavailability x0.97 (true AUC ratio about 0.95) %s%% and %s%%.",
              format(x$n_trials[1], big.mark = ","), ci(x[scenario == "S00", power_last_cmax], x[scenario == "S00", lo], x[scenario == "S00", hi]), ci(x[scenario == "S00", power_3], x[scenario == "S00", lo3], x[scenario == "S00", hi3]),
              ci(x[scenario == "F097", power_last_cmax], x[scenario == "F097", lo], x[scenario == "F097", hi]), ci(x[scenario == "F097", power_3], x[scenario == "F097", lo3], x[scenario == "F097", hi3])), "")
}

# OC detail(지시 2026-09-25 §4의 보고 우선순위): 1차 지표 경계 1종 오류(역산 경계표, 표, 그림 3-B, 치우침, P2 해석, G2 규칙 × 플래그 세트) →
# 운용특성 곡선(3-A) → 검정력 → 구성 비교 → 역산 전체(3-D) → 보조 지표 무작위 제품 공간(3-C, 가정한 가상 제품 분포에 의존)
bt <- R("oc", "boundary_type1.csv"); rk <- R("oc", "random_space_risks.csv"); pw_oc <- R("oc", "power.csv"); cc_oc <- R("oc", "config_comparison.csv")
ml <- c(k2016 = "2016", k2020 = "Model 1"); isb_ <- function(t) abs(t - bnd_[1]) < 1e-9 | abs(t - bnd_[2]) < 1e-9
bnd_txt <- paste(fmt_num(bnd_, 2), collapse = " or ")
if (!is.null(bt)) {
  ntr_ <- if ("n_trials" %in% names(bt)) oc_n_trials_text(bt[config == "P2"], "en", oc_dec_, ml) else paste(nb_, "each")   # 반복 수는 결과 파일에서(적응적 연장·부분 실행 반영)
  add(sprintf("## Operating characteristics, primary metric: boundary type I error (true AUC0-inf ratio %s)", bnd_txt), "",
      "The primary metric is the boundary type I error: the pass rate when the true AUC0-inf ratio is exactly at an equivalence limit, by mechanism, direction, configuration and model. Operating characteristic curves, power and the configuration comparison follow; the random product space is a secondary metric (last subsection).", "")
  if (!is.null(KF$inv_bnd)) {
    ivn_ <- unique(R("oc", "inversion_all.csv")$n_subjects)
    add(sprintf("Boundary scenarios: test-arm multiplier giving a true AUC0-inf ratio of %s and the true Cmax ratio at that multiplier (%s CRN subjects, 60 to 90 kg, %s mg, log-multiplier bisection to within %s%%). Unreachable: the range-end multiplier on the side of the target and its true AUC0-inf ratio.",
                paste(fmt_num(bnd_, 2), collapse = " and "), format(ivn_, big.mark = ","), oc_cfg$estimand$dose_mg, format(100 * oc_cfg$inversion$tolerance_rel)), "",
        md_table(key_inv_bnd_table(KF$inv_bnd, "en"))) }
  add(md_table(dcast(bt[config %in% c("P2", "F3A", "F3C", "G2")], model + mechanism + direction + target + auc_ratio + cmax_ratio ~ config, value.var = "pass_pct")[
        , .(Model = ml[model], Mechanism = mechanism, Direction = direction, `Target` = f2(target), `True AUC0-inf ratio` = f3(auc_ratio), `True Cmax ratio` = f3(cmax_ratio),
            `P2 (%)` = f2(P2), `F3-A (%)` = f2(F3A), `F3-C (%)` = f2(F3C), `G2 (%)` = f2(G2))]),
      sprintf("Wilson 95%% intervals are in the full report (Appendix B). Trials per scenario: %s.", ntr_), "")
  for (m_ in c("k2016", "k2020")) fig(sprintf("oc/fig3B_boundary_type1_%s_en.png", m_), sprintf("Figure 3-B (%s). Boundary type I error by mechanism and configuration, true AUC0-inf ratio %s (trials per scenario %s; Wilson 95%% intervals; dashed line 5%%).", MLF[[m_]], bnd_txt, oc_ntxt(m_, "bnd_only")))
}
gb <- R("oc", "gmr_by_endpoint.csv")
if (!is.null(gb)) {
  x <- dcast(gb[endpoint %in% c("AUClast", "AUCinf_A", "AUCinf_B", "AUCinf_C", "Cmax") & isb_(target)], model + mechanism + direction + target + auc_ratio + cmax_ratio ~ endpoint, value.var = "rel_bias_pct")
  add("Relative bias (%) of the geometric mean of trial GMRs against the truth at the boundaries (AUC endpoints against the true AUC0-inf ratio, Cmax against the true Cmax ratio). Rule A combines estimation and selection (flagged subjects excluded), rule B is estimation only, rule C substitutes AUClast for flagged subjects:", "",
      md_table(x[order(model, mechanism, target), .(Model = ml[model], Mechanism = mechanism, Direction = direction, Target = f2(target), `AUClast` = f1(AUClast), `AUCinf rule A` = f1(AUCinf_A),
                                                    `AUCinf rule B` = f1(AUCinf_B), `AUCinf rule C` = f1(AUCinf_C), `Cmax` = f1(Cmax))]))
}
# P2 해석(scripts/43 결과): 불편 추정량 기준, 분류, 원인
p2en <- proj_path("results", "oc", "p2_interpretation_en.md")
if (file.exists(p2en)) add("## P2 boundary type I error: interpretation (directive 2026-09-25, section 3)", "",
  "The P2 boundary type I error split into the baseline of an unbiased estimator (theory and the individual model AUC0-inf of the trial subjects), a classification of each cell (Wilson 95% interval against 5%) and its causes (AUClast bias, Cmax failures), from the saved trial-level results (scripts/43_p2_interpretation.R; no new simulation; numbers in oc/p2_interpretation.csv).", "",
  md_body(p2en, 1))
# G2: AUCinf 처리 규칙 × 플래그 세트(검토 의견 W2 §2, scripts/41 결과)
g2en <- proj_path("results", "oc", "g2_rules_conclusion_en.md"); pbz <- R("oc", "p2_bias_boundary.csv")
if (file.exists(g2en)) {
  add("## G2 boundary type I error by AUCinf handling rule and flag set (review W2, section 2)", "",
      sprintf("G2 (AUCinf + Cmax) re-judged from the saved trial-level results under rules A, B and C and flag sets (i) and (ii) (scripts/41_oc_rules_flags.R; no new simulation). %s; rule B does not depend on the flags. Full tables: oc/g2_rules_flags.csv and oc/g2_decomposition.csv.", REL_CONV_EN), "",
      md_body(g2en, 1))
  if (!is.null(pbz)) { x <- dcast(pbz, model + scenario + target ~ endpoint, value.var = "bias_pct")
    add("Relative bias (%) of the geometric mean of trial GMRs against the truth in the boundary scenarios (AUC endpoints against the true AUC0-inf ratio, Cmax against the true Cmax ratio; AUCinf rules A and C under flag set (ii); AUCinf true = individual model AUC0-inf of the trial subjects, reference; oc/p2_bias_boundary.csv):", "",
        md_table(x[order(model, target, scenario), .(Model = ml[model], Scenario = scenario, Target = f2(target), AUClast = f2(AUClast), `AUCinf rule A` = f2(AUCinf_A), `AUCinf rule B` = f2(AUCinf_B),
                        `AUCinf rule C` = f2(AUCinf_C), `AUCinf rule A, set (i)` = ifelse(is.na(AUCinf_Ai), "NA", f2(AUCinf_Ai)), `AUCinf rule C, set (i)` = ifelse(is.na(AUCinf_Ci), "NA", f2(AUCinf_Ci)), `AUCinf true` = f2(AUCinf_true), Cmax = f2(Cmax))])) }
}
# 운용특성 곡선, 검정력, 구성 비교. 시험 수준 비율은 5,000회 이상 실행에서만 구간과 함께 인용(이 문서 규칙); 2,000회 시나리오는 보고서
if (!is.null(bt) || !is.null(pw_oc) || !is.null(cc_oc)) {
  add("## Operating characteristics: OC curves, power and configuration comparison", "")
  for (m_ in c("k2016", "k2020")) fig(sprintf("oc/fig3A_oc_curves_%s_en.png", m_), sprintf("Figure 3-A (%s). Pass probability of P2, F3-A, F3-C and G2 against the true AUC0-inf ratio, by mechanism (B0, 117 per arm, 60 to 90 kg; trials per scenario: boundaries and identical products %s, elsewhere %s).", MLF[[m_]], oc_ntxt(m_, "bnd"), oc_ntxt(m_, "other")))
  CFGe <- c(P2 = "P2", F3A = "F3-A", F3C = "F3-C", G2 = "G2")
  if (!is.null(pw_oc)) { big <- pw_oc[config %in% names(CFGe) & n_trials >= 5000]; small <- unique(pw_oc[n_trials < 5000, n_trials])
    if (nrow(big)) add(sprintf("Power (pass rate, %%, Wilson 95%% CI) in scenarios run with at least 5,000 trials%s:", if (length(small)) sprintf("; power at true ratios %s (%s trials per mechanism) is tabulated in the full report", paste(fmt_num(sort(unique(pw_oc[n_trials < 5000, target])), 2), collapse = " and "), paste(format(sort(small), big.mark = ","), collapse = " or ")) else ""), "",
      md_table(dcast(big[, .(model, scenario, target, auc_ratio, n_trials, config, v = ci(pass_pct, lo, hi, 1))], model + scenario + target + auc_ratio + n_trials ~ config, value.var = "v")[
        , .(Model = ml[model], Scenario = scenario, Target = f2(target), `True AUC0-inf ratio` = f3(auc_ratio), Trials = format(n_trials, big.mark = ","), `P2 (%)` = P2, `F3-A (%)` = F3A, `F3-C (%)` = F3C, `G2 (%)` = G2)])) }
  if (!is.null(cc_oc)) { x <- cc_oc[isb_(target) & n_trials >= 5000]
    if (nrow(x)) add("Configuration comparison at the boundaries (paired within the same trials, percentage points, 95% CI): P2 minus F3 is the additional protection of F3 (P2 passes, F3 fails), G2 minus P2 the difference in pass rate; inside and outside the limits the ranges are in the key conclusion and the full table in the full report:", "",
      md_table(x[order(model, mechanism, target), .(Model = ml[model], Mechanism = mechanism, Direction = direction, Target = f2(target), Trials = format(n_trials, big.mark = ","),
                                                    `P2 minus F3-A` = ci(P2_minus_F3A, P2_minus_F3A_lo, P2_minus_F3A_hi, 2), `P2 minus F3-C` = ci(P2_minus_F3C, P2_minus_F3C_lo, P2_minus_F3C_hi, 2),
                                                    `G2 minus P2` = ci(G2_minus_P2, G2_minus_P2_lo, G2_minus_P2_hi, 2))])) }
  fig("oc/fig3D_inversion_multipliers_en.png", "Figure 3-D. Test-arm multiplier needed for each target true AUC0-inf ratio, by mechanism and model (200,000 CRN subjects, 60 to 90 kg, 300 mg; bisection to within 0.1%; x marks targets not reachable within the search range).")
}
# 보조 지표: 무작위 제품 공간(가정한 가상 제품 분포에 의존)
if (!is.null(rk)) {
  rs_rng_en <- random_space_ranges_text(oc_cfg, " to "); nprod_ <- format(oc_cfg$random_space$n_products, big.mark = ",")
  x <- rk[config %in% c("P2", "F3A", "F3C", "G2") & truth == "AUC0-inf"]
  x[, metric_en := fifelse(startsWith(metric, "소비자"), "Consumer risk (pass when truth outside)", "Producer risk (fail when truth inside)")]
  x[, scope_en := fifelse(scope == "전체", "all", "near boundary")]
  add("## Secondary metric: random product space (depends on the assumed virtual product distribution)", "",
      sprintf("Secondary metric. The consumer and producer risks below depend on the assumed distribution of virtual products (Latin hypercube of %s products, log-uniform multipliers %s; config/oc_design.yaml). They are averages over that distribution and change if the distribution changes; configurations are judged on the primary metric, the boundary type I error.", nprod_, rs_rng_en), "",
      sprintf("Random product space (one trial per product, B0, 117 per arm; truth from %s CRN subjects per product):", format(oc_cfg$random_space$truth_subjects, big.mark = ",")), "",
      md_table(x[, .(Model = ml[model], Configuration = c(P2 = "P2", F3A = "F3-A", F3C = "F3-C", G2 = "G2")[config], Metric = metric_en, Scope = scope_en,
                     `Rate (%)` = ci(pct, lo, hi, 2), Products = format(n_products, big.mark = ",", trim = TRUE))]))
  for (m_ in c("k2016", "k2020")) fig(sprintf("oc/fig3C_random_space_%s_en.png", m_), sprintf("Figure 3-C (%s). Random product space, secondary metric that depends on the assumed virtual product distribution: %s Latin hypercube products with log-uniform multipliers (%s), one trial each (B0, 117 per arm); top, distribution of the true AUC0-inf ratio; bottom, pass rate per bin and logistic spline smooth.", MLF[[m_]], nprod_, rs_rng_en))
}

# 4. Pillar 3: 첫 문장은 역산(Km 0.01–100배)에서(key_facts, 전제 검사), 500회 쌍대 표는 보조 근거
p3 <- R("rationale", "pillar3_km_summary.csv"); p3d <- R("rationale", "pillar3_km_B0.csv")
if (!is.null(p3) || !is.null(KF$km)) add("## 4. Pillar 3: invisibility of binding-constant differences", "",
  if (!is.null(KF$km)) c(sprintf("**%s**", key_km_sentence(KF$km, "en")), "", "Source: oc/inversion_all.csv (range ends and reachable rows) and oc/inversion_scan_<model>_Km.csv (screening scan); boundary multipliers of all mechanisms are in the boundary scenario table above.", ""),
  if (!is.null(p3)) c(sprintf("Supporting evidence: test-arm Km multiplied by 0.5 to 10 (0.005 to 0.1 mg/L) changes the mean GMR of every endpoint by at most %s%% (2016 model) and %s%% (Model 1) relative to identical products (paired within the same 500 trials). Km is an MM approximation constant and is not identical to binding affinity.",
          f2(p3[model == "k2016", max_abs_GMR_change_pct]), f2(p3[model == "k2020", max_abs_GMR_change_pct])), ""))

# 5. Weight generalization
wb <- NULL; for (f in list.files(proj_path("results", "weight_generalization"), pattern = "^weight_bands_B0_[a-d]+\\.csv$", full.names = TRUE)) wb <- rbind(wb, fread(f))
if (!is.null(wb)) {
  ml <- c(a = "(a) 2016", b = "(b) Model 1", c = "(c) 2016 + ke~BMI", d = "(d) 2016 + ke~BMI + Vc~weight 0.817")
  add("## 5. Robustness across body weight (300 mg, B0, 20,000 subjects per uniform weight band)", "",
      "These results are robustness checks, not evidence for the endpoint proposal, which rests on the trial population (healthy adults of 60 to 90 kg, randomization stratified by body weight; section 2A). The adult atopic dermatitis body-weight distribution is reported only in Appendix I of the Modeling and Simulation Report.", "",
      "Covariate variants (c) and (d) use the adult coefficients of Kovalenko 2020 Model 4 (elimination rate constant ke proportional to (BMI/26)^0.368, central volume exponent 0.817); the BMI reference of 26 is a placeholder based on phase 3 mean BMI 25.4 to 27.3 (Kamal 2022). Height is simulated as normal (mean 170 cm, SD 9, truncated 150 to 195 cm). Development-data weight ranges are not reported in the source publications; bands above 130 kg are flagged as possible extrapolation.", "",
      md_table(wb[, .(Model = ml[model], `Band (kg)` = band, `Median BMI` = f1(BMI_median), `True extrap. median (%)` = f2(extrap_true_median), `95th pct` = f2(extrap_true_p95), `Window coverage <80% (%)` = f3(coverage_lt80_pct),
                      `AUCinf reliable, (i) [(ii)] (%)` = sprintf("%s [%s]", f1(reliable_rsq_extrap_pct), f1(reliable_pct)), `Lambda-z not estimable (%)` = f2(lambda_fail_pct), `AUClast geometric mean` = f1(AUClast_geo, 0), Note = ifelse(dev_range_note == "", "", "outside confirmed development range"))]))
}
ow <- NULL; for (mk in c("a", "b", "d")) { f <- proj_path("results", "weight_generalization", sprintf("obese_trials_per_endpoint_%s.csv", mk)); if (file.exists(f)) ow <- rbind(ow, fread(f)) }
od <- NULL; for (mk in c("a", "b", "d")) { f <- proj_path("results", "weight_generalization", sprintf("obese_trials_dropout_%s.csv", mk)); if (file.exists(f)) od <- rbind(od, fread(f)) }
if (!is.null(ow)) {
  x <- dcast(ow[endpoint %in% c("AUClast", "AUCinf_true", "AUCinf_reliable")], model + scenario ~ endpoint, value.var = "GMR_mean")
  d <- od[, .(retained = mean(wt_retained_mean), dropped = mean(wt_dropout_mean), dropout_pct = mean(dropout_pct)), by = .(model, scenario)]
  x <- merge(x, d, by = c("model", "scenario"))
  add("Trial level in a heavier population (stress test, not a population estimate; weight normal mean 100 kg, SD 20, truncated 60 to 150 kg; 2,000 trials per scenario; mean GMR and dropout weights only):", "",
      md_table(x[, .(Model = ml[model], Scenario = scenario, `GMR AUClast` = f3(AUClast), `GMR true AUCinf` = f3(AUCinf_true), `GMR NCA AUCinf reliable` = f3(AUCinf_reliable),
                     `Subjects failing AUCinf reliability, flag set (ii) (%)` = f1(dropout_pct), `Mean weight retained (kg)` = f1(retained), `Mean weight failing (kg)` = f1(dropped))]))
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
             `D3: reliability gain, (i) [(ii)] (pp)` = { pv <- if (!is.null(RF)) RF$all[RF$all$variant == v & RF$all$schedule == "D3"] else NULL
               if (!is.null(pv) && nrow(pv) == 1) FPe(pv$gain_i_pp, pv$gain_ii_pp, 2, "", TRUE) else sprintf("[(ii) %s]", f2(d3$c_reliable_gain_pp)) },
             `D3: extrapolation >20% ratio, 20,000 subjects` = sprintf("%s (%s to %s)", f2(d3$d_extrap20_ratio), f2(d3$d_ratio_boot_lo), f2(d3$d_ratio_boot_hi)),
             `D3: same ratio, 200,000 subjects` = if (!is.null(k3) && nrow(k3)) sprintf("%s (%s to %s), %s fewer subjects per arm; %s", f2(k3$extrap_gt20_ratio), f2(k3$d_ratio_boot_lo), f2(k3$d_ratio_boot_hi), f2(-k3$d_abs_change_per_arm),
               if (k3$d_ratio_boot_lo <= 0.5 && k3$d_ratio_boot_hi >= 0.5) "fragile (interval includes 0.5)" else if (k3$d_ratio_boot_lo > 0.5) "not met (interval above 0.5)" else "met (interval below 0.5)") else "not re-evaluated",
             `Final recommendation` = if (any(d$recommend)) "B0 (conclusion unchanged)" else "B0") }))
if (nrow(rows)) add(md_table(rows))
bm <- R("trials", "bminus_explicit.csv"); wd <- R("trials", "width_decomposition.csv")
if (!is.null(wd)) { w <- wd[schedule %in% c("D1", "D2", "D4")]
  add(sprintf("AUClast CI width: D1, D2 and D4 widen the mean AUClast 90%% CI by %s%% to %s%% (paired intervals exclude zero; D3 %s%%). Recomputed with true concentrations and no residual error, the widening remains (%s%% to %s%%), so it reflects heterogeneity of the added tail area as the last quantifiable time is extended, not measurement error at the added low points.",
              f2(-100 * max(w$base)), f2(-100 * min(w$base)), f2(-100 * wd[schedule == "D3", base]), f2(-100 * max(w$noresid)), f2(-100 * min(w$noresid))), "") }
if (!is.null(bm)) { bmr <- if (!is.null(RF)) RF$dense$bminus[model == "k2016"] else NULL
  add(sprintf("Removing Day 50 (B-): AUClast CI width %s%% (negative = narrower), AUCinf reliability change %s percentage points (2016 model, paired 95%% CI), share with NCA extrapolation above 20%% multiplied by %s. Not recommended, to keep a terminal sample for AUC0-inf as a secondary endpoint and for a fallback analysis.",
              f2(-100 * bm$a_mean_width_rel_decrease),
              if (!is.null(bmr) && nrow(bmr) == 1) sprintf("(i) %s [(ii) %s]", cie(bmr$gain_i_pp, bmr$gain_i_lo, bmr$gain_i_hi), cie(bmr$gain_ii_pp, bmr$gain_ii_lo, bmr$gain_ii_hi)) else sprintf("[(ii) %s]", ci(bm$c_reliable_gain_pp, bm$c_gain_boot_lo, bm$c_gain_boot_hi, 2)),
              ci(bm$d_extrap20_ratio, bm$d_ratio_boot_lo, bm$d_ratio_boot_hi, 2)), "") }
trig <- if (nrow(rows)) rows[startsWith(`Rule result`, "D3"), Variant] else character(0)
abc_any <- any(unlist(lapply(c(names(vv), "iiv150", "resid12", "weight_alt", "ada10"), function(v) { d <- R("trials", sprintf("schedule_decision_%s.csv", v)); if (is.null(d)) FALSE else d[, any(crit_a %in% TRUE | crit_b %in% TRUE | crit_c %in% TRUE)] })))
k_lo <- unlist(lapply(c("base", "struct2020", "vmax080_both"), function(v) { k <- R("individual200k", sprintf("criterion_d_200k_%s.csv", v)); if (is.null(k)) NA_real_ else k[schedule == "D3", d_ratio_boot_lo] }))
if (abc_any) stop("criterion (a), (b) or (c) met in some variant: rationale text must be revised")
if (anyNA(k_lo) || any(k_lo <= 0.5)) stop("200,000-subject criterion (d) interval reaches 0.5 in some variant: rationale text must be revised")
# 전제: 절대 감소 arm당 1명 미만(규칙상 권고 행과 20만 명 재평가), 주 모델 D1–D4에서 AUClast CI 폭·ke ×1.10 통과율 개선 없음. 추가 방문 수는 판정표에서
dd_ <- rbindlist(lapply(names(vv), function(v) { d <- R("trials", sprintf("schedule_decision_%s.csv", v)); if (is.null(d)) NULL else d[recommend %in% TRUE] }), fill = TRUE)
k_abs <- unlist(lapply(c("base", "struct2020", "vmax080_both"), function(v) { k <- R("individual200k", sprintf("criterion_d_200k_%s.csv", v)); if (is.null(k)) NA_real_ else k[schedule == "D3", d_abs_change_per_arm] }))
premise(abs(c(dd_$d_abs_change_per_arm, k_abs)) < 1, "absolute reduction below one subject per arm")
db_ <- R("trials", "schedule_decision_base.csv")
premise(db_[schedule %in% c("D1", "D2", "D3", "D4"), all(a_mean_width_rel_decrease <= 0 & b_ke110_gain_pp <= 0)], "neither the AUClast CI width nor the ke x1.10 pass rate improves in the primary model")
if (!is.null(RF)) premise(!(RF$all$crit_c_i %in% TRUE) & !(RF$all$crit_c_ii %in% TRUE), "criterion (c) not met under either flag set")
add(sprintf("Rationale: %s; criteria (a), (b) and (c) were not met in any variant (criterion (c) under neither reliability flag set); re-evaluated with 200,000 subjects (4,000 bootstrap resamples), the D3 ratio for criterion (d) is above 0.5 with its whole interval in every re-evaluated variant; the absolute reduction is below one subject per arm; in the primary model neither the AUClast CI width nor power improves with any added schedule; the cost of D3 is %s additional visits. Removing the Day 50 sample (B-) is not recommended because it preserves a terminal point for AUC0-inf as a secondary endpoint and for a fallback analysis. Lesson recorded: a relative-reduction criterion for a rare event needs an absolute floor (for example at least one subject per arm); not applied retroactively.",
            if (length(trig)) sprintf("at the pre-specified 20,000-subject level, criterion (d) alone was met in %s", paste(trig, collapse = " and ")) else "no variant met any criterion at the pre-specified 20,000-subject level",
            format(db_[schedule == "D3", added_visits_total], big.mark = ",")), "")

# 8. Limitations
# 예측 편향: 연구 제형 gate 데이터셋의 모의/관측 비(전제: AUClast는 2016 모델이 Model 1보다 관측에 가깝다)
pb_txt <- if (!is.null(g16) && !is.null(g20)) { a16 <- g16[gate_role == "gate"]; a20 <- g20[gate_role == "gate"]; r2 <- function(x) sprintf("%s to %s", f2(min(x)), f2(max(x)))
  premise(mean(abs(a16$AUClast_ratio - 1)) < mean(abs(a20$AUClast_ratio - 1)), "AUClast of the 2016 model is closer to observed than Model 1")
  sprintf("- Simulated/observed ratios over the %d study-presentation gate datasets: Cmax %s (2016 model) and %s (Model 1); AUClast %s (2016 model, close to observed) and %s (Model 1).", nrow(a16), r2(a16$Cmax_ratio), r2(a20$Cmax_ratio), r2(a16$AUClast_ratio), r2(a20$AUClast_ratio)) } else character(0)
resid_en <- if (!is.null(RF) && !is.null(RF$resid)) sprintf(" (proportional residual 12%%: reliability %s, extrapolation above 20%% %s%%; 2016 model %s and %s%%); two-model ranges: reliability %s, extrapolation above 20%% %s%% to %s%%",
  FPe(RF$resid$rel_i, RF$resid$rel_ii), f2(RF$resid$gt20), FPe(RF$b0[model == "k2016", rel_i], RF$b0[model == "k2016", rel_ii]), f2(RF$resid$gt20_base), FPe(RF$rng$rel_i, RF$rng$rel_ii),
  f2(RF$rng$gt20[1]), f2(RF$rng$gt20[2])) else ""
add("## 8. Limitations", "",
    "- The primary model does not reproduce the faster absorption of the 200 mg 1.14 mL (175 mg/mL) presentation; about one third to one half of the 200 mg shortfall is explained by absorption, the remainder is unexplained. Clot 2021 Chinese 200 mg data (same 175 mg/mL prefilled syringe) showed median time to Cmax 7.0 days (range 3.0 to 10.0), but with n=8 and no sample between Day 4 and Day 8 the resolution is low.",
    pb_txt,
    "- Km is fixed at 0.01 mg/L in both models; with no uncertainty or IIV on Km, terminal-phase variability may be underestimated. Km and Vmax sensitivity analyses address this.",
    sprintf("- The AUCinf reliability rate and the share with NCA extrapolation above 20%% depend on the residual error model and are reported as two-model ranges%s.", resid_en),
    if (!is.null(RF)) sprintf("- %s. %s", REL_CONV_EN, rel_changes_en),
    "- Placeholders not yet confirmed: Day 1 post-dose sampling time (0.25 day), weight distribution and stratification split, sampling windows, BMI reference of 26, body weights of the Li 2020 single-arm studies.",
    "- Development-data body weight ranges are not reported; results above 130 kg are extrapolations.",
    "- Only the automatic lambda-z Best Fit is simulated; in a real study a pharmacokineticist may review and adjust the lambda-z points.",
    sprintf("- The random product space is a secondary metric: its consumer and producer risks depend on the assumed distribution of virtual products (log-uniform multipliers %s). The truth of each product is computed from %s common virtual subjects (not 200,000) for computational reasons, as pre-specified; the Monte Carlo standard error of each product's truth is reported.",
            random_space_ranges_text(oc_cfg, " to "), format(oc_cfg$random_space$truth_subjects, big.mark = ",")), "")
# 9. 재현성 점검(scripts/38_repro_check.R 로컬, 38b GitHub 기록): 결과 파일이 있을 때만
rl_ <- R("repro", "repro_check.csv"); rgr_ <- R("repro", "repro_github_run.csv"); rgi_ <- R("repro", "repro_github_items.csv"); rlm_ <- R("repro", "repro_check_meta.csv")
if (!is.null(rl_)) {
  kv_ <- function(d, k) if (is.null(d)) NA_character_ else d$value[d$field == k]
  lab_ <- c(a1 = "Exact regeneration of committed operating-characteristic trials (2016 model, trials 1 to 100, identical products and Vmax down 1.25): largest relative difference in GMR and 90% CI limits",
            a2 = "Same regeneration: rows whose pass flag, missing status or arm sizes differ (of 1,200)",
            a3_P2_pass_rate_S00_pct = "Same regeneration: P2 pass rate (%), identical products",
            a3_P2_pass_rate_Vmax_down_125_pct = "Same regeneration: P2 pass rate (%), Vmax down 1.25 (boundary)",
            b1 = "True AUC0-inf ratio at the Vmax inversion multiplier: first 20,000 of the 200,000 common virtual subjects",
            b2 = "Inversion screening consistency on the same 20,000 subjects",
            c1 = "Median true extrapolated share at B0 (%): 2,000 new subjects versus the committed 20,000")
  g8 <- function(v) format(v, digits = 8)
  lab_of <- function(it) ifelse(it %in% names(lab_), lab_[it], lab_[substr(it, 1, 2)])   # 전체 항목 이름 우선, 없으면 접두(a1, b1, ...)
  tb_ <- data.table(Item = unname(lab_of(rl_$item)), Committed = vapply(rl_$committed, g8, ""), `Reproduced (local)` = vapply(rl_$reproduced, g8, ""),
                    Tolerance = paste(vapply(rl_$tolerance, function(v) format(v, digits = 3), ""), rl_$tolerance_type), `Pass (local)` = ifelse(rl_$pass, "yes", "NO"))
  if (anyNA(tb_$Item)) stop("repro_check.csv has an item without an English label: ", paste(rl_$item[is.na(tb_$Item)], collapse = ", "))
  if (!is.null(rgi_)) { m_ <- match(rl_$item, rgi_$item); if (anyNA(m_)) stop("repro_github_items.csv does not cover every local item")
    tb_[, `Reproduced (GitHub)` := vapply(rgi_$reproduced_github[m_], g8, "")][, `Pass (GitHub)` := ifelse(rgi_$pass_github[m_], "yes", "NO")] }
  gh_txt <- if (!is.null(rgr_)) sprintf(" GitHub Actions run %s (commit %s, %s, %s, rxode2 %s, clean renv restore, %s cores; job %s s): %s of %s items pass; all %s reproduced values are %s to the local run at the 8 significant digits printed in the job log (%s).",
      kv_(rgr_, "run_number"), substr(kv_(rgr_, "commit"), 1, 7), kv_(rgr_, "runner"), kv_(rgr_, "r_version"), kv_(rgr_, "rxode2_version"), kv_(rgr_, "cores_used"), kv_(rgr_, "job_duration_s"),
      kv_(rgr_, "n_pass_github"), kv_(rgr_, "n_items"), kv_(rgr_, "n_items"), if (identical(kv_(rgr_, "all_identical_to_local_8sig"), "TRUE")) "identical" else "NOT all identical", kv_(rgr_, "url")) else " The GitHub Actions run has not been recorded yet."
  add("## 9. Reproducibility check", "",
      sprintf("Pre-specified tolerances (config/repro_check.yaml), script scripts/38_repro_check.R. Local run: commit %s, %s, rxode2 %s, %s of %s items pass.%s",
              substr(kv_(rlm_, "commit"), 1, 7), kv_(rlm_, "r_version"), kv_(rlm_, "rxode2_version"), sum(rl_$pass), nrow(rl_), gh_txt), "",
      md_table(tb_))
}
# 10. v1.0.1 추가 분석(config/prereg_20260926.yaml section1-4): 분석 모형, 기준 관행, LLOQ, 표본 수. 각 요약 스크립트가 만든 영문 결론을 그대로 넣는다
v101 <- list(c("10.1 Analysis model: pooled t-test (M0), weight-stratum ANOVA (M1), log-weight ANCOVA (M2)", proj_path("results", "oc_models", "oc_models_conclusion_en.md")),
             c("10.2 Lambda-z reliability criteria convention (sets (i) to (iv))", proj_path("results", "criteria", "criteria_conclusion_en.md")),
             c("10.3 Study-assay LLOQ sensitivity", proj_path("results", "lloq", "lloq_conclusion_en.md")),
             c("10.4 Sample size", proj_path("results", "sample_size", "ss_conclusion_en.md")))
if (any(vapply(v101, function(z) file.exists(z[2]), TRUE))) {
  add("## 10. Version 1.0.1 analyses (pre-registered in config/prereg_20260926.yaml)", "",
      "The trial population is the evidence base (section 2A); the adult atopic dermatitis results in 10.2 are a robustness check only.", "")
  for (z in v101) if (file.exists(z[2])) add(paste0("### ", z[1]), "", md_body(z[2], shift = 2))
}
txt <- paste(out, collapse = "\n")
if (grepl("—|–|\u2212", txt)) { bad <- regmatches(txt, gregexpr("[^\n]*(—|–|\u2212)[^\n]*", txt))[[1]]; stop("summary_en.md contains an em-dash, en-dash or minus sign (U+2212): ", paste(substr(head(bad, 3), 1, 160), collapse = " || ")) }
if (grepl("[가-힣]", txt)) { bad <- regmatches(txt, gregexpr("[^\n]*[가-힣][^\n]*", txt))[[1]]; stop("summary_en.md contains Korean text: ", paste(head(bad, 3), collapse = " || ")) }
writeLines(txt, proj_path("results", "summary_en.md"), useBytes = TRUE)
cat("results/summary_en.md written:", length(out), "lines\n")
