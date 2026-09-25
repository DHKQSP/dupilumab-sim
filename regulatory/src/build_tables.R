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
asy <- read_cfg("assay.yaml")
lloq <- data.table(Model = "both", Parameter = "LLOQ (study assay, config/assay.yaml)", Value = asy$lloq_mg_L$value, Unit = "mg/L", Status = "assumption for the sponsor's assay", Source = asy$lloq_mg_L$source)
pp <- rbind(th16, cov16, iiv16, res16, th20, cov20, iiv20, res20, lloq)
pp[, Parameter := gsub("^Vc$", "Vc (central volume, L at 75 kg; V2 in Kovalenko 2016)", Parameter)]
pp[, Parameter := gsub("^k12$", "k12 (k23 in Kovalenko 2016; kcp in 2020)", Parameter)]; pp[, Parameter := gsub("^k21$", "k21 (k32 in Kovalenko 2016; kpc in 2020)", Parameter)]
pp[, Value := as.character(Value)]
fwrite(pp, file.path(out, "parameter_provenance.csv"))

# ---- B. assumptions register (curated) ----
td <- yaml::read_yaml(proj_path("config", "trial_design.yaml"))
ab <- data.table(
  ID = sprintf("B%02d", 1:15),
  Assumption = c("Bioanalytical LLOQ of the sponsor's assay", "Day 1 post-dose sampling time", "Visit windows", "Body weight distribution and sex ratio of enrolled subjects",
                 "Randomization stratification", "Evaluable subjects per arm; dropout", "Km fixed, without between-subject variability", "No parameter uncertainty layer",
                 "Residual error model", "Immunogenicity (ADA)", "Automated lambda-z selection without manual review", "Presentation of test and reference products", "Population: healthy subjects",
                 "Body weight of single-arm literature studies (validation only)", "BMI reference for the ke-BMI covariate (weight generalization only)"),
  `Value used` = c(sprintf("%s mg/L", asy$lloq_mg_L$value), sprintf("%s day after dosing", td$day1_postdose_time$value), "plus or minus 2 h to Day 1, 6 h to Day 14, 1 day thereafter",
                   sprintf("normal, mean %s kg, SD %s kg, truncated %s to %s kg; male fraction %s", td$weight$base$mean, td$weight$base$sd, td$weight$base$trunc[1], td$weight$base$trunc[2], td$weight$sex_ratio_male$value),
                   sprintf("two strata split at %s kg, 1:1 within strata", td$stratification$split_kg$value), sprintf("%s evaluable per arm (%s randomized); dropout not simulated", td$n_per_arm, td$n_randomized_per_arm),
                   "0.01 mg/L in both models", "single-layer Monte Carlo (parameters fixed at published estimates)", "as published: proportional 24.2% (2016) and 15.0% (2020), additive 0.03 mg/L",
                   "not in the operating characteristics; ADA-like subgroup sensitivity (10% of subjects, ke x2 after Day 14)", "Phoenix Best Fit rules", "300 mg as 2 mL of 150 mg/mL for both", "healthy adults, BMI 18 to 32 (BMI not modelled)",
                   "normal, mean 78 kg, SD 10 kg", "26 kg/m2"),
  Basis = c("Clot 2021 (originator's assay)", "placeholder; protocol not final", "placeholder; typical protocol windows", "placeholder based on the inclusion range 60 to 90 kg",
            "placeholder", "sample-size assumption of the protocol", "as published (Km determined by likelihood profiling; the objective function was insensitive below 0.01 mg/L)", "parameter estimation uncertainty was not propagated",
            "as published", "no ADA model in the published PK models", "standard NCA practice", "study design", "study design", "weights not reported in Li 2020", "phase 3 mean BMI 25.4 to 27.3 (Kamal 2022)"),
  `Impact evidence` = c("The cliff starts at a median of about 0.9 mg/L and lasts about 1.4 days, so a different LLOQ below that level is expected to move tlast by at most about the cliff length and to change coverage little; not simulated.",
                        "Affects Cmax sampling only; terminal-phase results unaffected.", "Simulated in all trial-level analyses; cliff analysis with and without windows.",
                        "Alternative distribution (mean 72, SD 10, 50 to 90 kg) and weight bands 40 to 150 kg: coverage preserved; schedule recommendation unchanged.",
                        "Candidate cause of the between-trial SD / within-trial SE ratio below 1 (conservative).", "Operating characteristics computed at the evaluable number; power about 99% for identical products.",
                        "Km x0.5 to x10 in both arms and Km x0.01 to x100 in the test arm: coverage and conclusions unchanged.", "Two structural models, variance x1.5 and curve-shape sensitivity; conclusions unchanged.",
                        "Reliability rates depend on the residual error (reported as two-model ranges; proportional 12% sensitivity); operating-characteristic conclusions identical in both models.",
                        "In the ADA-like subgroup coverage is preserved and the schedule recommendation unchanged.", "Manual review could change lambda-z windows and AUC0-inf; AUC0-last and Cmax unaffected.",
                        "The 200 mg 175 mg/mL presentation (faster absorption) is not reproduced; not relevant if both products use the 300 mg 150 mg/mL presentation.", "Results apply to healthy adults 60 to 90 kg.",
                        "External-only validation result depends on it (Section 4.2); not used in any analysis of the proposal.", "Used only for the weight generalization sensitivity."),
  `Action before submission` = c("Confirm the LLOQ of the validated assay; if it differs by more than about 2-fold, rerun coverage and reliability at the actual LLOQ.", "Set to the protocol time.", "Align with the protocol.",
                                 "Confirm expected enrolment; rerun only if materially different.", "Confirm the randomization plan and whether the analysis model includes the stratum.", "Confirm the sample-size justification (log-scale CV assumption).",
                                 "None; state as limitation.", "None; state as limitation.", "None; report residual-sensitive metrics as two-model ranges.", "Consider the ADA incidence reported for the reference product in healthy subjects; describe ADA handling in the SAP.",
                                 "Pre-specify lambda-z rules in the SAP and document any manual changes.", "Confirm presentations of both products.", "None.", "None.", "None."))
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
               "Premise checks of generated conclusions", "Independent human QC of this report"),
  Method = c("Automated test: linear limit (Vmax = 0) against the closed-form two-compartment solution", "Automated tests", "Automated tests (linear-up log-down AUC, lambda-z windows, BLQ rules, ties)",
             "Theoph, Indometh and 1,000 simulated profiles: this engine, NonCompart 0.8.4, PKNCA 0.12.1", "Automated tests against t.test (var.equal) and lm/confint", "Automated tests over all product scenarios",
             "Automated tests (YAML boolean keys, types, data.table scoping lint)", "Automated tests", "Reviewer's Python implementation, 20,000 subjects", "Reviewer's Python implementation, 8,000 subjects per condition",
             "5,000-trial set versus the earlier 500-trial run (trials 1 to 500)", "500-trial estimates versus the independent later trials", "Regenerated boundary trials versus stored rows (6 original endpoints)",
             "Pre-specified tolerances (config/repro_check.yaml); GitHub-hosted runner with a fresh renv restore", "Fast-scope tests on every commit that changes code, configuration or dependencies",
             "Generated conclusion texts stop with an error if a stated premise is contradicted by the results", "Line-by-line check of every number against regulatory/traceability.csv and of the text against the results"),
  `Acceptance criterion` = c("Relative difference below tolerance", "As specified in the tests", "As specified in the tests", "Identical lambda-z windows; relative difference at most 1e-6", "Identical to base R", "Reference arm unchanged; test arm shifted by the multiplier",
                             "No violation", "Same inputs give same seeds", "Within 3% (small percentages: absolute difference)", "Within 10% (non-rare metrics)", "Maximum difference 0", "Differences within Monte Carlo error",
                             "Relative 1e-6 and identical pass flags", "All items within tolerance", "No failure", "No premise violated", "No discrepancy"),
  Result = c("Pass", "Pass", "Pass",
             sprintf("Pass: %s of %s window comparisons identical; largest relative difference %s", format(sum(ev$lz_points_identical), big.mark = ","), format(sum(ev$lz_points_identical + ev$lz_points_mismatch), big.mark = ","), format(signif(max(ev$max_rel_diff), 2))),
             "Pass", "Pass", "Pass", "Pass", sprintf("%d of %d metrics within 3%%; the others are small percentages that differ by fractions of a percentage point", sum(cv1$agree_3pct, na.rm = TRUE), sum(!is.na(cv1$agree_3pct))),
             sprintf("%d of %d within 10%%", sum(abs(rr5n$rel_diff_pct) <= 10, na.rm = TRUE), nrow(rr5n)), sprintf("%s rows, maximum difference %s", format(idn$rows_compared, big.mark = ","), idn$max_abs_diff),
             sprintf("%d of %d comparisons outside Monte Carlo error", sum(mc5$outside_mc), nrow(mc5)), sprintf("%s stored rows matched (%s per model)", format(sum(rj_n), big.mark = ","), format(rj_n[1], big.mark = ",")),
             sprintf("Local %d of %d; clean runner %s of %s, identical to the local values to 8 significant digits", sum(rl$pass), nrow(rl), gv("n_pass_github"), gv("n_items")),
             sprintf("%d successful runs after the fix; latest code commit %s: %s", sum(ci$workflow == "tests" & ci$conclusion == "success"), ci[workflow == "tests"][.N, commit], ci[workflow == "tests"][.N, conclusion]),
             "Pass (all generated texts produced)", "PENDING (sponsor)"),
  Evidence = c("tests/testthat/test-model-structure.R", "tests/testthat/test-model-structure.R", "tests/testthat/test-nca.R, test-nca-wnl.R", "results/nca_engine/engine_validation_summary.csv",
               "tests/testthat/test-be-stats.R", "tests/testthat/test-scenario-propagation.R", "tests/testthat/test-config-schema.R, test-lint-datatable-scope.R", "tests/testthat/test-seeds.R",
               "results/crossval/crossval_model1.csv", "results/crossval/reviewer_reference_round5.csv", "results/trials5000/mc_consistency_identity.csv", "results/trials5000/mc_consistency_500_vs_5000.csv",
               "logs/oc_rejudge_k2016.out, logs/oc_rejudge_k2020.out", "results/repro/repro_check.csv, results/repro/repro_github_items.csv, results/repro/repro_github_run.csv",
               "results/ci/ci_runs_after_fix.csv", "scripts/33_oc_summary.R, 36_cliff_conclusion.R, 39_reliability_flags.R, 43_p2_interpretation.R", "signature page of the report"))
fwrite(qc, file.path(out, "verification_qc.csv"))

# ---- D. pre-specification and post hoc register (dates from git) ----
dd <- data.table(
  Item = c("Models and parameter values", "Model validation criteria", "Scope of model validation restricted to the study presentation", "Sampling-schedule decision rule (criteria a to d)",
           "NCA engine replaced by a Phoenix-compatible engine", "Operating-characteristic design (mechanisms, targets, trials, configurations, seeds, truth definition)", "Cliff analysis design",
           "Reliability criteria set (i) (without the span ratio)", "AUC0-inf + Cmax under rules B and C, and rules A and C under criteria set (i)", "Extension to 20,000 trials of a boundary scenario whose Wilson interval included 5%",
           "Interpretation of AUC0-last + Cmax at the boundaries (unbiased baseline, classification, decomposition)", "Random product space reported as a secondary metric"),
  Status = c("pre-specified", "pre-specified", "post hoc (after the first validation results)", "specified after the first schedule simulations were run, before any result was reported",
             "change after initial results; all dependent outputs regenerated", "pre-specified", "pre-specified", "post hoc", "post hoc", "post hoc (data-dependent)", "post hoc", "post hoc reporting priority; the analysis itself was pre-specified"),
  `Date (evidence)` = c(cdate("8d69d82"), cdate("8d69d82"), cdate("c018729"), cdate("c018729"), cdate("d26f169"), cdate("779e068"), cdate("779e068"), cdate("0f6bf72"), cdate("0f6bf72"), cdate("68be707"), cdate("68be707"), cdate("68be707")),
  `First results` = c(cdate("d2592d7"), cdate("d2592d7"), "same commit as the decision", "first schedule simulations: see Date of the models row; results withheld until validation was accepted",
                      "outputs regenerated before the operating characteristics were run", cdate("b8a5351"), cdate("5f972bd"), cdate("0f6bf72"), cdate("bf9a5b1"), cdate("bf9a5b1"), cdate("bf9a5b1"), cdate("b8a5351")),
  `How reported` = c("Appendix A", "Section 4.2", "200 mg data sets reported as external checks with their ratios; re-judgement on fully external data reported", "Section 5.7; recommendation unchanged in every variant",
                     "Previous versus new engine differences tabulated (results/nca_engine/); conclusions unchanged", "Primary metric; every mechanism and configuration reported (Appendix F)", "Section 5.1",
                     "Reported first, with set (ii) alongside (set (ii) is the pre-specified definition)", "Labelled post hoc in Table 5-5; same trials regenerated with the same seeds", "Both the pre-specified 10,000-trial value and the 20,000-trial value are reported",
                     "Section 5.4", "Section 5.5"))
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
