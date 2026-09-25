# Builds the appendix tables of the regulatory package (regulatory/tables/*.csv), in English.
# Parameter provenance comes from the model configuration files; commit dates from the version-control history;
# verification results from the result files. Curated text (assumptions, pre-specification status) is written here and
# reviewed with the sponsor.
# Usage: Rscript regulatory/src/build_tables.R   (called by scripts/60_regulatory_package.R)
source("R/00_setup.R"); source("regulatory/src/reg_helpers.R")
out <- proj_path("regulatory", "tables"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
rd <- function(p) fread(proj_path(p), encoding = "UTF-8")
git1 <- function(...) { r <- suppressWarnings(system2("git", shQuote(c("-C", PROJ_ROOT, ...)), stdout = TRUE, stderr = FALSE)); if (length(r)) r[1] else NA_character_ }
cdate <- function(h) { d <- git1("show", "-s", "--format=%cd", "--date=format-local:%Y-%m-%d %H:%M UTC", h); premise(!is.na(d) && nzchar(d), paste("git date for commit", h)); sprintf("%s (commit %s)", d, h) }
Sys.setenv(TZ = "UTC")

# ---- A. parameter provenance (from config) ----
pt <- yaml::read_yaml(proj_path("config", "params_typical.yaml")); pv <- yaml::read_yaml(proj_path("config", "params_variability.yaml")); p20 <- yaml::read_yaml(proj_path("config", "params_k2020_model1.yaml"))
src_en <- function(x) {
  x <- gsub("\\(고정\\)", "(fixed)", x); x <- gsub("\\(고정, profiling으로 결정\\)", "(fixed; determined by likelihood profiling)", x)
  x <- gsub("kpc, Mpc 0.686에서 유도", "kpc, derived from Mpc 0.686", x); x <- gsub("\\(지시서 §2[^)]*\\)", "", x)
  x <- gsub("\\(기준 75 kg, 검토 의견 통합본 §0-2\\)", "(reference weight 75 kg)", x); x <- gsub("\\(비례 ([0-9.]+)% CV\\)", "(proportional \\1% CV)", x)
  x <- gsub("\\(비례 ([0-9.]+)%\\)", "(proportional \\1%)", x); x <- gsub("\\(가산 0.03 mg/L 고정\\)", "(additive 0.03 mg/L, fixed)", x)
  x <- gsub("지시서 §2: k23 IIV 없음", "Kovalenko 2016 Table 2 (no IIV estimated)", x); x <- gsub("지시서 §2: k32 IIV 없음", "Kovalenko 2016 Table 2 (no IIV estimated)", x)
  x <- gsub("지시서 §2: F IIV 없음", "Kovalenko 2016 Table 2 (no IIV estimated)", x); x <- gsub("^Km 고정$", "Km fixed; no IIV", x)
  x <- gsub("Clot 2021 §2.4; 지시서 §2", "Clot 2021 section 2.4 (assay of the originator; sponsor assay to be confirmed)", x)
  x <- gsub("Clot 2021 §2.4", "Clot 2021 section 2.4 (assay of the originator; sponsor assay to be confirmed)", x)
  trimws(gsub("\\s+", " ", x))
}
th16 <- rbindlist(lapply(names(pt$theta), function(k) { e <- pt$theta[[k]]; data.table(Model = "Kovalenko 2016 (primary)", Parameter = k, Value = e$value, Unit = e$unit, Status = e$status, Source = src_en(e$source)) }))
iiv16 <- rbindlist(lapply(names(pv$iiv_omega2), function(k) { e <- pv$iiv_omega2[[k]]; data.table(Model = "Kovalenko 2016 (primary)", Parameter = paste0("IIV variance ", k), Value = e$omega2, Unit = "log-scale variance", Status = e$status, Source = src_en(e$source)) }))
res16 <- data.table(Model = "Kovalenko 2016 (primary)", Parameter = c("Residual proportional", "Residual additive"), Value = c(pv$residual$sigma_prop$value, pv$residual$sigma_add$value),
                    Unit = c("fraction (SD)", "mg/L (SD)"), Status = "confirmed", Source = src_en(c(pv$residual$sigma_prop$source, pv$residual$sigma_add$source)))
cov16 <- data.table(Model = "Kovalenko 2016 (primary)", Parameter = c("Weight exponent on Vc", "Reference weight"), Value = c(pt$covariates$WT_on_Vc$theta_WT$value, pt$covariates$WT_on_Vc$WT_ref$value),
                    Unit = c("", "kg"), Status = "confirmed", Source = "Kovalenko 2016 Table 2")
th20 <- rbindlist(lapply(names(p20$theta), function(k) { e <- p20$theta[[k]]; data.table(Model = "Kovalenko 2020 Model 1 (sensitivity)", Parameter = k, Value = e$value, Unit = if (is.null(e$unit)) "" else e$unit, Status = e$status, Source = src_en(e$source)) }))
iiv20 <- rbindlist(lapply(names(p20$iiv$sd), function(k) { e <- p20$iiv$sd[[k]]; data.table(Model = "Kovalenko 2020 Model 1 (sensitivity)", Parameter = paste0("IIV variance ", k), Value = e$omega2, Unit = sprintf("log-scale variance (SD %s)", e$sd), Status = "confirmed",
                                                                                          Source = if (e$omega2 > 0) "Kovalenko 2020 Supplementary Table 2 (Model 1; reported as SD, squared; covariances not reported, diagonal)" else "not estimated in Kovalenko 2020 Model 1") }))
res20 <- data.table(Model = "Kovalenko 2020 Model 1 (sensitivity)", Parameter = c("Residual proportional", "Residual additive"), Value = c(p20$residual$sigma_prop$value, p20$residual$sigma_add$value),
                    Unit = c("fraction (SD)", "mg/L (SD)"), Status = "confirmed", Source = src_en(c(p20$residual$sigma_prop$source, p20$residual$sigma_add$source)))
cov20 <- data.table(Model = "Kovalenko 2020 Model 1 (sensitivity)", Parameter = c("Weight exponent on Vc", "Reference weight"), Value = c(p20$covariates$WT_on_Vc$theta_WT$value, p20$covariates$WT_on_Vc$WT_ref$value),
                    Unit = c("", "kg"), Status = "confirmed", Source = c("Kovalenko 2020 Table 1 Model 1", "Kovalenko 2020 Model 1 (reference weight 75 kg)"))
asy <- yaml::read_yaml(proj_path("config", "assay.yaml"), fileEncoding = "UTF-8")
lloq <- data.table(Model = "both", Parameter = "LLOQ (study assay, config/assay.yaml)", Value = asy$lloq_mg_L$value, Unit = "mg/L", Status = "assumption for the sponsor's assay", Source = asy$lloq_mg_L$source)
pa_ <- yaml::read_yaml(proj_path("config", "population_atopic.yaml")); vb_ <- yaml::read_yaml(proj_path("config", "scenarios.yaml"))$sensitivity_variants$k2016_bmi_vc0817
premise(!is.null(vb_$ke_bmi$exponent) && !is.null(vb_$theta_WT_override), "phase 3 covariate variant in config/scenarios.yaml")
AT_SRC <- c(primary = "Kamal 2022 (phase 2b, 379 adults with atopic dermatitis: mean body weight 74.0 to 80.6 kg by dose group); Kovalenko 2021 (adult phase 3 population PK data: BMI SD 5.47 kg/m2); NCT03389893 (71 adults: 80.2 kg, SD 19.0); used only in the Appendix I robustness check",
            sens_nct03389893 = "NCT03389893 (71 adults with atopic dermatitis: 80.2 kg, SD 19.0)", sens_75_18 = "sensitivity with a lower mean (analyst)")
premise(identical(sort(names(pa_$distributions)), sort(names(AT_SRC))), "atopic distributions in config/population_atopic.yaml")
atw <- rbindlist(lapply(names(pa_$distributions), function(d) { x <- pa_$distributions[[d]]
  data.table(Model = "Appendix I robustness check (adult atopic dermatitis body weight)", Parameter = sprintf("Body weight, %s", if (d == "primary") "primary distribution" else "sensitivity distribution"), Value = x$label, Unit = "kg",
             Status = "placeholder (assumption)", Source = AT_SRC[[d]]) }))
atw <- rbind(atw, data.table(Model = "Appendix I robustness check (adult atopic dermatitis body weight)", Parameter = c("Height (for BMI)", "Male fraction"),
                             Value = c(sprintf("normal, mean %s cm, SD %s cm, truncated %s to %s cm, independent of weight", pa_$height$mean, pa_$height$sd, pa_$height$trunc[[1]], pa_$height$trunc[[2]]), as.character(pa_$sex_ratio_male)),
                             Unit = c("cm", ""), Status = "placeholder (assumption)", Source = c("NCT03389893 (172.5 cm, SD 10.4)", "assumption")),
             data.table(Model = "2016 model with phase 3 covariates (weight sensitivity)", Parameter = c("Exponent of BMI on linear elimination (ke)", "Reference BMI", "Weight exponent on Vc (replaces the 2016 value)"),
                        Value = c(as.character(vb_$ke_bmi$exponent), as.character(vb_$ke_bmi$bmi_ref), as.character(vb_$theta_WT_override)), Unit = c("", "kg/m2", ""), Status = c("confirmed", "assumption", "confirmed"),
                        Source = c("Kovalenko 2020 Table 3 and Supplementary Table 3 (Model 4)", "phase 3 mean BMI 25.4 to 27.3 (Kamal 2022)", "Kovalenko 2020 Table 3 and Supplementary Table 3 (Model 4)")))
pp <- rbind(th16, cov16, iiv16, res16, th20, cov20, iiv20, res20, lloq, atw)
pp[, Parameter := gsub("^Vc$", "Vc (central volume, L at 75 kg; V2 in Kovalenko 2016)", Parameter)]
pp[, Parameter := gsub("^k12$", "k12 (k23 in Kovalenko 2016; kcp in 2020)", Parameter)]; pp[, Parameter := gsub("^k21$", "k21 (k32 in Kovalenko 2016; kpc in 2020)", Parameter)]
pp[, Value := as.character(Value)]
fwrite(pp, file.path(out, "parameter_provenance.csv"))

# ---- version 1.0.1 facts (analysis model, LLOQ, sample size) used in the registers below ----
lt_ <- rd("results/lloq/lloq_individual_table.csv")[model == "k2016" & resid == "fixed"]; ltt_ <- rd("results/lloq/lloq_trial_type1.csv")[resid == "fixed" & config == "P2"]
t1_ <- rd("results/oc_models/type1_models.csv"); ss_ <- rd("results/oc_models/sd_se_models.csv")[endpoint == "AUCinf_true" & !scenario %in% c("S00", "F097")]
pp_ <- rd("results/sample_size/ss_table_power.csv")[input_model == "k2016" & gmr == 0.95 & n == 117]
f1_ <- function(x) formatC(round(x, 1) + 0, format = "f", digits = 1); f2_ <- function(x) formatC(round(x, 2) + 0, format = "f", digits = 2); f3_ <- function(x) formatC(round(x, 3) + 0, format = "f", digits = 3)
rg_ <- function(x, f = f2_) sprintf("%s to %s", f(min(x)), f(max(x)))
cls_ <- function(a) { x <- t1_[analysis_model == a & config == "P2"]; sprintf("%d conservative, %d nominal, %d exceeding", sum(x$class == "conservative"), sum(x$class == "nominal"), sum(x$class == "exceeding")) }
lloq_impact <- sprintf("Simulated at 0.02 to 0.5 mg/L with the same subjects (Section 5.10): true coverage of AUC0-last below 80%% in %s%% of subjects; reliability under criteria (i) %s%% (0.02 mg/L) to %s%% (0.5 mg/L); boundary type I error of AUC0-last + Cmax %s%% (M0) and %s%% (M1) in three boundary scenarios.",
                       rg_(lt_$coverage_lt80_pct), f1_(lt_[abs(lloq - 0.02) < 1e-9, reliable_no_span_pct]), f1_(lt_[abs(lloq - 0.5) < 1e-9, reliable_no_span_pct]), rg_(ltt_[model == "M0", pass_pct]), rg_(ltt_[model == "M1", pass_pct]))
strat_impact <- sprintf("The pooled t-test (M0) does not model the stratum: between-trial SD / within-trial SE of the unbiased reference %s (M0) against %s with the stratum in the model (M1). AUC0-last + Cmax: M0 %s; M1 %s (Section 5.9).",
                        rg_(ss_[analysis_model == "M0", sd_se_ratio], f3_), rg_(ss_[analysis_model == "M1", sd_se_ratio], f3_), cls_("M0"), cls_("M1"))
ac_ <- rd("results/atopic/atopic_criteria_by_band.csv")[distribution == "primary" & band == "all"]; ax_ <- rd("results/atopic/atopic_exploratory_height_corr.csv"); aw_ <- rd("results/atopic/atopic_weight_table.csv")
at_impact <- sprintf("%s%% of simulated patients outside 60 to 90 kg; without a reliable AUC0-inf: %s%% (criteria set (i)) and %s%% (set (iii)) across three model variants; sensitivity distributions (80/19 and 75/18 kg): %s%% and %s%% (Appendix I; robustness check, not evidence).",
                     f1_(aw_[distribution == "primary" & grepl("^simulated", source), outside_60_90]), rg_(ac_[set == "i", fail_pct], f1_), rg_(ac_[set == "iii", fail_pct], f1_),
                     rg_(rd("results/atopic/atopic_criteria_by_band.csv")[distribution != "primary" & band == "all" & set == "i", fail_pct], f1_), rg_(rd("results/atopic/atopic_criteria_by_band.csv")[distribution != "primary" & band == "all" & set == "iii", fail_pct], f1_))
ht_impact <- sprintf("Overstates the BMI SD (%s against 5.45 kg/m2). With the correlation implied by NCT03389893 (%s), failure of the BMI-covariate variant changes from %s%% to %s%% (set (i)); exploratory, not pre-registered (Appendix I).",
                     f2_(ax_[height_model == "independent", bmi_sd]), f2_(ax_[height_model == "correlated", rho_target]), f2_(ax_[height_model == "independent", fail_pct_i]), f2_(ax_[height_model == "correlated", fail_pct_i]))
n_impact <- sprintf("P2 power with 117 evaluable per arm at a true GMR of 0.95: %s%% (M0) and %s%% (M1) at CV 43%%, %s%% and %s%% at CV 50%% (Section 5.11).",
                    f1_(pp_[cv == 43 & analysis_model == "M0", analytic_pct]), f1_(pp_[cv == 43 & analysis_model == "M1", analytic_pct]), f1_(pp_[cv == 50 & analysis_model == "M0", analytic_pct]), f1_(pp_[cv == 50 & analysis_model == "M1", analytic_pct]))

# ---- B. assumptions register (curated) ----
td <- yaml::read_yaml(proj_path("config", "trial_design.yaml"))
ab <- data.table(
  ID = sprintf("B%02d", 1:17),
  Assumption = c("Bioanalytical LLOQ of the sponsor's assay", "Day 1 post-dose sampling time", "Visit windows", "Body weight distribution and sex ratio of enrolled subjects",
                 "Randomization stratification", "Evaluable subjects per arm; dropout", "Km fixed, without between-subject variability", "No parameter uncertainty layer",
                 "Residual error model", "Immunogenicity (ADA)", "Automated lambda-z selection without manual review", "Presentation of test and reference products", "Population: healthy subjects",
                 "Body weight of single-arm literature studies (validation only)", "BMI reference for the ke-BMI covariate (weight generalization only)",
                 "Body weight of adults with atopic dermatitis (Appendix I robustness check)", "Height independent of body weight (Appendix I robustness check)"),
  `Value used` = c(sprintf("%s mg/L", asy$lloq_mg_L$value), sprintf("%s day after dosing", td$day1_postdose_time$value), "plus or minus 2 h to Day 1, 6 h to Day 14, 1 day thereafter",
                   sprintf("normal, mean %s kg, SD %s kg, truncated %s to %s kg; male fraction %s", td$weight$base$mean, td$weight$base$sd, td$weight$base$trunc[1], td$weight$base$trunc[2], td$weight$sex_ratio_male$value),
                   sprintf("two strata split at %s kg, 1:1 within strata", td$stratification$split_kg$value), sprintf("%s evaluable per arm (%s randomized); dropout not simulated", td$n_per_arm, td$n_randomized_per_arm),
                   "0.01 mg/L in both models", "single-layer Monte Carlo (parameters fixed at published estimates)", "as published: proportional 24.2% (2016) and 15.0% (2020), additive 0.03 mg/L",
                   "not in the operating characteristics; ADA-like subgroup sensitivity (10% of subjects, ke x2 after Day 14)", "Phoenix Best Fit rules", "300 mg as 2 mL of 150 mg/mL for both", "healthy adults, BMI 18 to 32 (BMI not modelled)",
                   "normal, mean 78 kg, SD 10 kg", "26 kg/m2", pa_$distributions$primary$label, sprintf("height normal, mean %s cm, SD %s cm, independent of weight", pa_$height$mean, pa_$height$sd)),
  Basis = c("Clot 2021 (originator's assay)", "placeholder; protocol not final", "placeholder; typical protocol windows", "placeholder based on the inclusion range 60 to 90 kg",
            "placeholder", "sample-size assumption of the protocol", "as published (Km determined by likelihood profiling; the objective function was insensitive below 0.01 mg/L)", "parameter estimation uncertainty was not propagated",
            "as published", "no ADA model in the published PK models", "standard NCA practice", "study design", "study design", "weights not reported in Li 2020", "phase 3 mean BMI 25.4 to 27.3 (Kamal 2022)", "published summaries (Kamal 2022, Kovalenko 2021, NCT03389893)", "specified in the request"),
  `Impact evidence` = c(lloq_impact,
                        "Affects Cmax sampling only; terminal-phase results unaffected.", "Simulated in all trial-level analyses; cliff analysis with and without windows.",
                        "Alternative distribution (mean 72, SD 10, 50 to 90 kg) and weight bands 40 to 150 kg: coverage preserved; schedule recommendation unchanged.",
                        strat_impact, n_impact,
                        "Km x0.5 to x10 in both arms and Km x0.01 to x100 in the test arm: coverage and conclusions unchanged.", "Two structural models, variance x1.5 and curve-shape sensitivity; conclusions unchanged.",
                        "Reliability rates depend on the residual error (reported as two-model ranges; proportional 12% sensitivity); operating-characteristic conclusions identical in both models.",
                        "In the ADA-like subgroup coverage is preserved and the schedule recommendation unchanged.", "Manual review could change lambda-z windows and AUC0-inf; AUC0-last and Cmax unaffected.",
                        "The 200 mg 175 mg/mL presentation (faster absorption) is not reproduced; not relevant if both products use the 300 mg 150 mg/mL presentation.", "Results apply to healthy adults 60 to 90 kg.",
                        "External-only validation result depends on it (Section 4.2); not used in any analysis of the proposal.", "Used only for the weight generalization sensitivity.", at_impact, ht_impact),
  `Action before submission` = c("Set the LLOQ of the validated assay in config/assay.yaml and regenerate all results (scripts/run_all.sh).", "Set to the protocol time.", "Align with the protocol.",
                                 "Confirm expected enrolment; rerun only if materially different.", "Confirm the randomization plan; select the primary analysis model (M0 or M1) in the SAP.", "Select the target power (90% or 85%) and the sample size (Section 5.11).",
                                 "None; state as limitation.", "None; state as limitation.", "None; report residual-sensitive metrics as two-model ranges.", "Consider the ADA incidence reported for the reference product in healthy subjects; describe ADA handling in the SAP.",
                                 "Pre-specify lambda-z rules in the SAP and document any manual changes.", "Confirm presentations of both products.", "None.", "None.", "None.",
                                 "None; robustness check only, not evidence for the proposal (principle of section 6).", "None; robustness check only."))
fwrite(ab, file.path(out, "assumptions_register.csv"))

# ---- C. verification and QC record ----
ev <- .translate(rd("results/nca_engine/engine_validation_summary.csv")); cv1 <- rd("results/crossval/crossval_model1.csv"); rr5 <- rd("results/crossval/reviewer_reference_round5.csv")
idn <- rd("results/trials5000/mc_consistency_identity.csv"); mc5 <- rd("results/trials5000/mc_consistency_500_vs_5000.csv"); rg <- rd("results/repro/repro_github_run.csv"); rl <- rd("results/repro/repro_check.csv")
ci <- rd("results/ci/ci_runs_after_fix.csv"); rj <- c(readLines(proj_path("logs", "oc_rejudge_k2016.out")), readLines(proj_path("logs", "oc_rejudge_k2020.out")))
rj_n <- as.numeric(sub("^done: ([0-9]+) stored rows matched.*$", "\\1", grep("^done: [0-9]+ stored rows matched", rj, value = TRUE)))
premise(length(rj_n) == 2, "two rejudge completion lines")
rr5n <- rr5[is.na(note) | note == ""]
gv <- function(k) rg$value[rg$field == k]
qc <- data.table(
  Activity = c("Model code against the analytical solution", "Mass balance, Michaelis-Menten zero-order limit, transit mean transit time", "NCA rules on hand-calculated examples",
               "NCA engine against reference implementations", "Equivalence statistics", "Scenario multipliers act on the test arm only", "Configuration schema and coding rules",
               "Deterministic seeds", "Independent re-implementation (Model 1)", "Independent re-implementation (curve shape and weight bands)", "Trial regeneration identity",
               "Monte Carlo consistency", "Regeneration for post hoc rules (same seeds)", "Reproducibility in a clean environment", "Continuous integration",
               "Premise checks of generated conclusions", "Analysis models M1 and M2 (closed form)", "Regeneration for the analysis-model re-judgement (same seeds)",
               "LLOQ re-censoring at the current LLOQ", "Sample-size analytic approximation", "Restart of the regenerated trials (criteria sets (iii) and (iv) added)",
               "AUC0-inf + Cmax under the new criteria file versus version 1.0", "Adult atopic dermatitis population summaries", "Independent human QC of this report"),
  Method = c("Automated test: linear limit (Vmax = 0) against the closed-form two-compartment solution", "Automated tests", "Automated tests (linear-up log-down AUC, lambda-z windows, BLQ rules, ties)",
             "Theoph, Indometh and 1,000 simulated profiles: this engine, NonCompart 0.8.4, PKNCA 0.12.1", "Automated tests against t.test (var.equal) and lm/confint", "Automated tests over all product scenarios",
             "Automated tests (YAML boolean keys, types, data.table scoping lint)", "Automated tests", "Reviewer's Python implementation, 20,000 subjects", "Reviewer's Python implementation, 8,000 subjects per condition",
             "5,000-trial set versus the earlier 500-trial run (trials 1 to 500)", "500-trial estimates versus the independent later trials", "Regenerated boundary trials versus stored rows (6 original endpoints)",
             "Pre-specified tolerances (config/repro_check.yaml); GitHub-hosted runner with a fresh renv restore", "Fast-scope tests on every commit that changes code, configuration or dependencies",
             "Generated conclusion texts stop with an error if a stated premise is contradicted by the results", "Automated tests against lm/confint, including unequal arms and a singular design",
             "Regenerated M0 versus stored rows of the operating-characteristic and product runs", "Individual, cliff and trial results at 0.078 mg/L versus the stored results; rescaled residual at scale 1 versus the stored residual",
             "Simulation of every grid cell (5,000 trials) and PK-model trials at n = 117", "Rows of the trials completed before the restart versus the same rows after it",
             "M0 pass rates of rules A and C (criteria (i), (ii)) and rule B versus results/oc/g2_rules_flags.csv where the trial counts are equal",
             "Analytic weight shares versus the independent reference; simulated versus analytic shares; stored flags recomputed; stored run reproduced before the exploratory check",
             "Line-by-line check of every number against regulatory/traceability.csv and of the text against the results"),
  `Acceptance criterion` = c("Relative difference below tolerance", "As specified in the tests", "As specified in the tests", "Identical lambda-z windows; relative difference at most 1e-6", "Identical to base R", "Reference arm unchanged; test arm shifted by the multiplier",
                             "No violation", "Same inputs give same seeds", "Within 3% (small percentages: absolute difference)", "Within 10% (non-rare metrics)", "Maximum difference 0", "Differences within Monte Carlo error",
                             "Relative 1e-6 and identical pass flags", "All items within tolerance", "No failure", "No premise violated", "Relative difference at most 1e-10",
                             "Relative 1e-6 (1e-12 for product runs) and identical pass flags", "Identical (tolerance 0) or relative 1e-9", "Analytic value inside the Wilson 95% interval", "Identical", "Identical (tolerance 1e-9)", "Within 0.06 points; binomial z at most 4; identical", "No discrepancy"),
  Result = c("Pass", "Pass", "Pass",
             sprintf("Pass: %s of %s window comparisons identical; largest relative difference %s", format(sum(ev$lz_points_identical), big.mark = ","), format(sum(ev$lz_points_identical + ev$lz_points_mismatch), big.mark = ","), format(signif(max(ev$max_rel_diff), 2))),
             "Pass", "Pass", "Pass", "Pass", sprintf("%d of %d metrics within 3%%; the others are small percentages that differ by fractions of a percentage point", sum(cv1$agree_3pct, na.rm = TRUE), sum(!is.na(cv1$agree_3pct))),
             sprintf("%d of %d within 10%%", sum(abs(rr5n$rel_diff_pct) <= 10, na.rm = TRUE), nrow(rr5n)), sprintf("%s rows, maximum difference %s", format(idn$rows_compared, big.mark = ","), idn$max_abs_diff),
             sprintf("%d of %d comparisons outside Monte Carlo error", sum(mc5$outside_mc), nrow(mc5)), sprintf("%s stored rows matched (%s per model)", format(sum(rj_n), big.mark = ","), format(rj_n[1], big.mark = ",")),
             sprintf("Local %d of %d; clean runner %s of %s, identical to the local values to 8 significant digits", sum(rl$pass), nrow(rl), gv("n_pass_github"), gv("n_items")),
             sprintf("%d successful runs after the fix; latest code commit %s: %s", sum(ci$workflow == "tests" & ci$conclusion == "success"), ci[workflow == "tests"][.N, commit], ci[workflow == "tests"][.N, conclusion]),
             "Pass (all generated texts produced)", "Pass",
             { rv <- rd("results/oc_models/m0_reverification.csv"); premise(all(rv$ok), "M0 re-verification"); sprintf("Pass: %s rows", format(sum(rv$n_rows), big.mark = ",")) },
             { ck <- rbind(rd("results/lloq/lloq_individual_check_k2016.csv"), rd("results/lloq/lloq_individual_check_k2020.csv"), rd("results/lloq/lloq_cliff_check.csv"), fill = TRUE); premise(all(ck$pass), "LLOQ checks"); sprintf("Pass: %d checks", nrow(ck)) },
             { mc <- rd("results/sample_size/ss_power_mc.csv"); pk <- rd("results/sample_size/ss_pk_check.csv")[analysis_model != "M2"]
               sprintf("Analytic inside the interval in %d of %d grid cells (largest difference %s points); PK-model trials: %d of %d", sum(mc$analytic_in_ci), nrow(mc), f2_(max(abs(mc$diff_pp))), sum(pk$analytic_in_ci), nrow(pk)) },
             { rc <- rd("results/criteria/restart_identity_check.csv"); premise(all(rc$ok), "restart identity"); sprintf("Pass: %s rows (%s trials)", format(sum(rc$rows_matched), big.mark = ","), paste(format(rc$trials_compared, big.mark = ","), collapse = " and ")) },
             { vc <- rd("results/criteria/criteria_check_vs_v10.csv"); premise(all(vc[same_n == TRUE, identical]), "criteria vs v1.0"); sprintf("Pass: %d of %d cells compared (equal trial counts), all identical", sum(vc$same_n), nrow(vc)) },
             "Pass (the summary stops on any failed check)",
             "PENDING (sponsor)"),
  Evidence = c("tests/testthat/test-model-structure.R", "tests/testthat/test-model-structure.R", "tests/testthat/test-nca.R, test-nca-wnl.R", "results/nca_engine/engine_validation_summary.csv",
               "tests/testthat/test-be-stats.R", "tests/testthat/test-scenario-propagation.R", "tests/testthat/test-config-schema.R, test-lint-datatable-scope.R", "tests/testthat/test-seeds.R",
               "results/crossval/crossval_model1.csv", "results/crossval/reviewer_reference_round5.csv", "results/trials5000/mc_consistency_identity.csv", "results/trials5000/mc_consistency_500_vs_5000.csv",
               "logs/oc_rejudge_k2016.out, logs/oc_rejudge_k2020.out", "results/repro/repro_check.csv, results/repro/repro_github_items.csv, results/repro/repro_github_run.csv",
               "results/ci/ci_runs_after_fix.csv", "scripts/33_oc_summary.R, 36_cliff_conclusion.R, 39_reliability_flags.R, 43_p2_interpretation.R", "tests/testthat/test-oc-models.R",
               "results/oc_models/m0_reverification.csv, logs/oc_models_k2016.out, logs/oc_models_k2020.out", "results/lloq/lloq_individual_check_k2016.csv, lloq_individual_check_k2020.csv, lloq_cliff_check.csv, logs/lloq_trials.out",
               "results/sample_size/ss_power_mc.csv, results/sample_size/ss_pk_check.csv", "results/criteria/restart_identity_check.csv", "results/criteria/criteria_check_vs_v10.csv",
               "scripts/53_atopic_summary.R; results/atopic/atopic_weight_table.csv, atopic_exploratory_height_corr.csv", "signature page of the report"))
fwrite(qc, file.path(out, "verification_qc.csv"))

# ---- D. pre-specification and post hoc register (dates from git) ----
fcommit <- function(path) { h <- suppressWarnings(system2("git", shQuote(c("-C", PROJ_ROOT, "log", "--diff-filter=A", "--format=%h", "--", path)), stdout = TRUE, stderr = FALSE))
  premise(length(h) >= 1 && nzchar(h[length(h)]), paste("result file committed:", path)); cdate(h[length(h)]) }
dd <- data.table(
  Item = c("Models and parameter values", "Model validation criteria", "Scope of model validation restricted to the study presentation", "Sampling-schedule decision rule (criteria a to d)",
           "NCA engine replaced by a Phoenix-compatible engine", "Operating-characteristic design (mechanisms, targets, trials, configurations, seeds, truth definition)", "Cliff analysis design",
           "Reliability criteria set (i) (without the span ratio)", "AUC0-inf + Cmax under rules B and C, and rules A and C under criteria set (i)", "Extension to 20,000 trials of a boundary scenario whose Wilson interval included 5%",
           "Interpretation of AUC0-last + Cmax at the boundaries (unbiased baseline, classification, decomposition)", "Random product space reported as a secondary metric",
           "Analysis-model re-judgement: M0, M1, M2 on regenerated trials; extension rule per model; expectations and reporting rule (version 1.0.1)",
           "Numeric criteria for the expectation checks", "Single study-LLOQ source and LLOQ sensitivity, including the scaled additive residual variant", "Sample-size table (analytic, simulation, PK-model check)",
           "Criteria sets (iii) and (iv): individual level from stored NCA, AUC0-inf + Cmax by rule and set, bias, decision instability", "Regenerated trials restarted from trial 1 to add the criteria-set endpoints",
           "Adult atopic dermatitis body-weight analyses (three model variants, three distributions; trials)", "Coverage measure in the summary sentence of the atopic analysis (sampling-window coverage of the true AUC0-inf)", "Weight-height correlation check (BMI covariate variant)",
           "Principle: the case against AUC0-inf rests on the trial population; atopic and out-of-range results reclassified as robustness checks", "Trial-population analyses: failure by criteria set, strata, treatment-dependent failure, failing-subject characteristics, coverage versus observed-to-true ratio"),
  Status = c("pre-specified", "pre-specified", "post hoc (after the first validation results)", "specified after the first schedule simulations were run, before any result was reported",
             "change after initial results; all dependent outputs regenerated", "pre-specified", "pre-specified", "post hoc", "post hoc", "post hoc (data-dependent)", "post hoc", "post hoc reporting priority; the analysis itself was pre-specified",
             "pre-registered (config/prereg_20260926.yaml section1)", "registered amendment: after the runs started, before any result was read", "pre-registered (section2); the scaled residual variant was added by the analyst, not requested", "pre-registered (section3)",
             "pre-registered (section4)", "registered change before any result of the restarted runs; identity of the earlier rows checked", "pre-registered (section5)",
             "decided after the first individual-level summary; the pre-registered observed-AUC0-last ratio is reported alongside", "post hoc, exploratory (after the simulated BMI SD was seen)",
             "sponsor decision after the atopic results were seen; no analysis changed, only where the results are reported", "pre-registered (section6); analyses seen before registration are listed in the registration"),
  `Date (evidence)` = c(cdate("8d69d82"), cdate("8d69d82"), cdate("c018729"), cdate("c018729"), cdate("d26f169"), cdate("779e068"), cdate("779e068"), cdate("0f6bf72"), cdate("0f6bf72"), cdate("68be707"), cdate("68be707"), cdate("68be707"),
                        cdate("521645a"), cdate("c5b304b"), cdate("211e8ed"), cdate("211e8ed"), cdate("325e63b"), cdate("325e63b"), cdate("325e63b"), cdate("167e2f5"), cdate("167e2f5"), cdate("bae3574"), cdate("bae3574")),
  `First results` = c(cdate("d2592d7"), cdate("d2592d7"), "same commit as the decision", "first schedule simulations: see Date of the models row; results withheld until validation was accepted",
                      "outputs regenerated before the operating characteristics were run", cdate("b8a5351"), cdate("5f972bd"), cdate("0f6bf72"), cdate("bf9a5b1"), cdate("bf9a5b1"), cdate("bf9a5b1"), cdate("b8a5351"),
                      fcommit("results/oc_models/type1_models.csv"), fcommit("results/oc_models/expectations_check.csv"), fcommit("results/lloq/lloq_trial_type1.csv"), fcommit("results/sample_size/ss_table_n_needed.csv"),
                      fcommit("results/criteria/criteria_g2_type1.csv"), fcommit("results/criteria/restart_identity_check.csv"), fcommit("results/atopic/atopic_criteria_by_band.csv"), fcommit("results/atopic/atopic_coverage_failing.csv"), fcommit("results/atopic/atopic_exploratory_height_corr.csv"), "not applicable (reporting decision)", fcommit("results/trialpop/tp_failure_by_set.csv")),
  `How reported` = c("Appendix A", "Section 4.2", "200 mg data sets reported as external checks with their ratios; re-judgement on fully external data reported", "Section 5.7; recommendation unchanged in every variant",
                     "Previous versus new engine differences tabulated (results/nca_engine/); conclusions unchanged", "Primary metric; every mechanism and configuration reported (Appendix F)", "Section 5.1",
                     "Reported first, with set (ii) alongside (set (ii) is the pre-specified definition)", "Labelled post hoc in Table 5-6; same trials regenerated with the same seeds", "Both the pre-specified 10,000-trial value and the 20,000-trial value are reported",
                     "Section 5.4", "Section 5.5", "Section 5.9; M0 and M1 side by side; primary model left to the sponsor", "Table 5-15", "Section 5.10 (residual as estimated first, scaled variant alongside)", "Section 5.11",
                     "Sections 5.3 and 5.4", "Section 3.13; results/criteria/restart_identity_check.csv", "Appendix I (robustness check)", "Appendix I", "Appendix I, labelled exploratory", "Sections 5.3, 5.8 and Appendix I", "Sections 5.2 and 5.3"))
fwrite(dd, file.path(out, "prespecification_register.csv"))

# ---- E. software environment (renv.lock) ----
lk <- jsonlite::fromJSON(proj_path("renv.lock"), simplifyVector = FALSE)
pv_ <- function(p) { x <- lk$Packages[[p]]$Version; premise(!is.null(x), paste("renv.lock:", p)); x }
se <- data.table(Software = c("R", "rxode2", "data.table", "NonCompart", "PKNCA", "digest", "yaml", "ggplot2", "patchwork", "knitr", "rmarkdown", "testthat", "renv"),
                 Version = c(lk$R$Version, vapply(c("rxode2", "data.table", "NonCompart", "PKNCA", "digest", "yaml", "ggplot2", "patchwork", "knitr", "rmarkdown", "testthat", "renv"), pv_, "")),
                 Role = c("Computing environment", "ODE solution of the PK models", "Data handling", "Reference NCA implementation (verification)", "Second reference NCA implementation (verification)", "Seed derivation and SHA-256",
                          "Configuration files", "Figures", "Figure layout", "Report generation", "Report generation", "Automated tests", "Package version locking"))
se <- rbind(se, data.table(Software = "Own code", Version = "this repository", Role = "Virtual population, sampling, NCA (Phoenix WinNonlin-compatible rules), equivalence statistics, operating characteristics (R/, scripts/)"))
fwrite(se, file.path(out, "software_environment.csv"))

# ---- label translation (Korean category labels in some analysis outputs) ----
fwrite(data.table(label_in_result_file = names(LABEL_EN), english_code = unname(LABEL_EN)), file.path(out, "label_translation.csv"))

for (f in c("parameter_provenance.csv", "assumptions_register.csv", "verification_qc.csv", "prespecification_register.csv", "software_environment.csv")) check_english(file.path(out, f))
cat("regulatory tables written to", out, "\n")
