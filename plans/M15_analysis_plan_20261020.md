# Model analysis plan: model-informed design of the dupilumab biosimilar phase 1 pharmacokinetic study

Prepared under the principles of the International Council for Harmonisation of Technical Requirements for Pharmaceuticals for Human Use (ICH) guideline M15, general principles for model-informed drug development. Registered before the full runs of the analyses it describes (section 0, development checks).

## 0. Document control

| Item | Entry |
|---|---|
| Plan identifier | Model analysis plan 2026-10-20, version 1 |
| Version history | 1, 2026-10-20, first registration |
| Status | Registered by the git commit that first adds this file together with its configuration file. Any later change is a deviation (section 13). |
| Author | Drafted with an artificial intelligence (AI) coding assistant for the sponsor's modelling team |
| Independent reviewer | Not yet assigned (pending before any regulatory use) |
| Approver | Not yet assigned (pending before any regulatory use) |
| Audience | Sponsor-internal pre-registration, written so that it can be shared with the US Food and Drug Administration (FDA) as the model analysis plan. Codes used in the configuration and result files are listed in Appendix B; internal traceability is in Appendix C. |
| Machine-readable registration | `config/prereg_20261020_m15.yaml`. The scripts read every value specific to this analysis from that file; design constants come from `config/trial_design.yaml`, `config/design_clot2021.yaml` and `config/assay.yaml`. |
| Decision record | Decision log entry D-068 and specification section 19 (both in Korean, internal) |
| Development checks before registration | Disclosed in section 13.2: the code was written and unit-tested; the external-validation code was run at 20 to 60 replicates and its verdicts were seen; the failure rule was revised between those runs; the population step and small test runs of the other scripts were executed. |
| Language and style | English; ranges are written with the word "to"; abbreviations are defined at first use (Appendix A) |

## 1. Introduction and background

The planned study is a randomized, parallel-group, single-dose study in healthy adults that compares a proposed dupilumab biosimilar (test) with the reference product after one subcutaneous (SC) dose. Similarity is concluded from the non-compartmental analysis (NCA) of the observed serum concentrations: the 90% confidence interval (CI) of the geometric mean ratio (GMR, test over reference) must lie within 0.80 to 1.25 for both the area under the concentration-time curve from time zero to the last quantifiable concentration (AUC0-last) and the maximum concentration (Cmax).

Dupilumab shows target-mediated pharmacokinetics (PK): parallel linear and saturable Michaelis-Menten (MM) elimination. Concentrations fall steeply once the saturable pathway dominates, so the terminal phase, the time of the last quantifiable concentration (Tlast) and the share of the area extrapolated beyond Tlast depend on body weight, dose and the sampling schedule. Body weight is the main covariate of the published population PK models.

This plan uses published population PK models already implemented in the sponsor's repository, without re-estimation or structural change. It evaluates how the design factors of the study affect its ability to detect PK differences between test and reference while still passing products that are similar. The design factors are the weight criteria and stratification, the sampling schedule, the sample size, the study population and the dose.

## 2. Question of interest, context of use and model risk

| Element | Statement |
|---|---|
| Question of interest | Does the phase 1 design (weight inclusion criteria and stratification, sampling schedule, sample size, study population, dose) detect PK differences between test and reference products sufficiently, while concluding similarity for products that are similar? |
| Context of use | Design rationale only. The equivalence conclusion of the study is made by NCA of the study's observed data; model results are not used in that conclusion. The model-based results support the choice of design factors and are reported as operating characteristics of the design. |
| Model influence | Medium. The model results are one of several inputs to design choices (with statistical sample-size methods, regulatory precedent and feasibility). They do not replace the observed-data analysis. |
| Consequence of a wrong decision | Medium. A model-informed design that is wrong could give a study with insufficient power (a similar product fails) or insufficient sensitivity (a design that is less able to detect a true PK difference). The study's two one-sided tests keep the nominal 5% consumer risk for the GMR of the measured endpoints (AUC0-last and Cmax). They do not guarantee it for total exposure when a product difference changes the terminal phase, because AUC0-last then under-represents the difference in the area to infinity; section 8.2 quantifies this. The consequence is rated medium because the study design can be revised before dosing and the observed-data NCA remains the basis of the conclusion. |
| Model risk | Medium (medium influence and medium consequence). |
| Model impact | Low. The model evidence informs design choices only and replaces no clinical study or observed-data analysis. |
| Regulatory interaction | Shared with FDA, if requested, as the model analysis plan for the phase 1 design; results are reported in a model analysis report. |
| Credibility activities proportionate to the risk | Code verification (section 6.1), external validation by posterior predictive checks against the original healthy-volunteer studies with a pre-specified failure rule and its false-alarm and detection probabilities (section 6.2), an applicability assessment (section 6.3), technical criteria with actions (section 6.4), a structural sensitivity model (Kovalenko 2016) for every main result, and an explicit statement of the uncertainty that is not propagated (section 9). |

### 2.1 Pre-specified reading of the results (descriptive; does not affect the study's NCA decision)

1. A design option is called adequate when, under the primary model (Kovalenko 2020 Model 1, 10,000 simulated studies at the margins), the simulated power at a true GMR of 0.95 is at least 90% and the upper Wilson 95% bound of the joint consumer risk at true GMRs of 0.80 and 1.25 is at most 5%. Kovalenko 2016 results are reported descriptively.
2. A design factor is reported as changing the operating characteristics when the paired Tango 95% CI of the difference in pass probability excludes 0 and the difference is at least 2 percentage points.
3. A model that passes the external-validation failure rule is described as "reproduced" only where the detection probability of the combined rule at a 25% bias (factors 0.75 and 1.25; the smaller of the two) is at least 50%; otherwise it is described as "not contradicted".

These readings are fixed at registration; a later change is a deviation.

### 2.2 Departures from the requested analysis

| Requested | Planned | Reason | Effect on conclusions | Confirmation needed |
|---|---|---|---|---|
| Primary model: Kovalenko 2020 adult covariate model (phase 3 Model 4) | Kovalenko 2020 Model 1 | The Model 4 parameter table is not recorded in the repository and the publication cannot be opened from the analysis environment; structure and estimates may not be changed | Weight is the only covariate; no race, albumin or disease term (section 4.1) | No (confirmed by the requester on 2026-10-21: Model 1 is Kovalenko 2020 Table 1 Model 1; section 16.14) |
| Parameter uncertainty propagated by the existing Monte Carlo and bootstrap | Not propagated; method registered as pending | No relative standard error (RSE), covariance matrix or parameter bootstrap is recorded (`config/param_uncertainty.yaml`) | All results are conditional on the published point estimates (section 9) | Yes (Model 1 bootstrap intervals, FDA 2017 review Table 4.6.6); amendment 1 (section 16.8) propagates parameter uncertainty from the Kovalenko 2016 bootstrap percentiles supplied on 2026-10-21 and, for Model 1, from the Kovalenko 2016 relative widths (borrowed uncertainty) until Table 4.6.6 is supplied (section 16.14) |
| Korean population from the Korea National Health and Nutrition Examination Survey (KNHANES) raw data | Weight range of the Japanese TDU12265 study (52.1 to 78.3 kg), as the request allows when raw data are absent | KNHANES raw data are not available in the analysis environment | The Korean population differs from the global one by weight only (section 7.1) | Yes (supply Korean sex-specific weight and height summaries) |
| Li 2020 Table 2 and Table 3 values and Table S1 weights for every requested arm | Only values already recorded in the repository | The publication cannot be opened | 4 of the 14 requested arms have no observed value; several studies lack Cmax or weights (section 5) | Yes (data request, section 6.2) |
| TDU12265 300 and 600 mg values from Li 2020 Table 2 | Values of Clot 2021 Table 5 for the same study | Recorded source in the repository; equality with Li 2020 not checked | Possible mismatch with the requested source | Yes |
| TDU12265 600 mg number of subjects | 6, from the request | The repository records no number for this arm | Prediction interval width for this arm | Yes |

### 2.3 Interpretations of the request

| Item | Interpretation used |
|---|---|
| Sample size 80 to 160 per arm | Randomized subjects per arm; evaluable subjects are about 0.9 times this after 10% withdrawal |
| Protocol assumption "CV 43%, 117 per arm" | Coefficient of variation (CV) of AUC0-last 43% (Syneos proposal) with 117 evaluable subjects per arm (130 randomized) |
| Stratification boundary "median of the enrolled population" | The median of the enrolled-population weight distribution, fixed a priori from 200,000 draws, because a per-study sample median is unknown at sequential randomization |
| Partial AUC "28 to 57 days" | 28 to 56 days after the dose, that is study Day 29 to study Day 57 (study Day 1 is the dosing day) |
| "Vm" in the mechanism analysis | The maximum MM elimination rate Vmax of the repository models |

## 3. Objectives

1. Confirm the model inputs used for healthy volunteers and check that the models reproduce the original healthy-volunteer studies of Li 2020 (section 6.2).
2. Estimate the operating characteristic (OC) curves of the similarity conclusion over a grid of true GMR values for each design scenario (section 8.1).
3. Estimate mechanism-specific OC curves when the test product differs from the reference in one PK parameter: bioavailability (F), absorption rate constant (ka), linear elimination rate constant (ke) or Vmax (section 8.2).
4. Quantify the secondary discriminating ability of endpoints that have no equivalence margin by the area under the receiver operating characteristic curve (AUROC). The endpoints are the area under the curve to infinity (AUC0-inf) and partial AUCs (section 8.3).
5. Evaluate the sampling schedules by the bias and precision of the NCA parameters against the model truth (section 8.4).
6. Give power and assurance against sample size, and compare with the protocol assumption (section 8.5).
7. Present the body-weight effect on exposure directly from the model covariate (section 8.6).

## 4. Models

### 4.1 Models used and their roles

| Role | Model | Repository source | Covariates |
|---|---|---|---|
| Primary | Kovalenko 2020 Model 1: two-compartment model with linear and MM elimination, absorption through three transit compartments (transit rate (n + 1) / MTT, where MTT is the mean transit time) and first-order absorption | `config/params_k2020_model1.yaml` (Kovalenko 2020 Table 1 Model 1 and Supplementary Table 2) | Body weight (WT, kg) on the central volume of distribution (Vc): Vc = 2.48 L x (WT / 75 kg)^0.711 |
| Sensitivity | Kovalenko 2016: two-compartment model with linear and MM elimination and first-order absorption | `config/params_typical.yaml` and `config/params_variability.yaml` (Kovalenko 2016 Table 2, estimates with data below the limit of quantification (BLQ) included) | Vc = 2.74 L x (WT / 75 kg)^0.705 |

The typical values, inter-individual variability (IIV) and residual error of both models are those in the configuration files and are not re-estimated or changed. The concentration-based Vmax enters the amount equation as Vmax x C / (Km + C) x Vc, with the Michaelis-Menten constant (Km) fixed at 0.01 mg/L. Bioavailability has no IIV in either model, so a multiplier on F changes every subject's F by exactly that factor.

The request names the Kovalenko 2020 adult covariate model as primary. The repository implements Kovalenko 2020 Model 1 only. The adult phase 3 covariate model (Model 4) has body weight, body mass index (BMI), albumin, race and baseline Eczema Area and Severity Index (EASI) terms. It is not implemented because its full parameter table is not recorded in the repository and the publication cannot be opened from the analysis environment. The repository had earlier judged Models 3 and 4 unsuitable for single-dose healthy-volunteer simulation (IIV on Vc and ke only, additive residual error of about 6 mg/L; decision log D-029). Two Model 4 coefficients (ke on BMI with exponent 0.368, Vc weight exponent 0.817) exist only as hybrid variants on the 2016 model; they are not used, because they would change the published estimates. Under the rule that the model structure and estimates are not changed, Model 1 is the primary model of this plan. Under Model 1 the differences between the global and Korean populations come from weight alone; no ethnic effect on PK parameters is represented.

Earlier records disagree with this choice of primary model, and they are disclosed here. The repository's earlier qualification gate (`config/gate_decision.yaml`) and the sampling-schedule decision (`config/schedule_decision.yaml`) were made with Kovalenko 2016 as primary; the schedule decision was later confirmed under Model 1. In the fully external re-judgement of that qualification, Model 1 met the plus or minus 15% AUC0-last criterion in 5 of 8 items and Kovalenko 2016 in 9 of 10 (results/step1_k2020/step1h_external_only_gate.csv). The 3 Model 1 misses were Li 2020 single arms without reported weights (PKM12350 test 1.27, reference 1.24; R668-HV-1108 reference 1.17), and they depend on the assumed 78 kg mean weight. For this context of use Model 1 is primary because the request names the Kovalenko 2020 model; section 6.2 tests it again. The earlier modelling and simulation report, which keeps Kovalenko 2016 as primary for the endpoint question, is not changed, and any FDA-facing summary of this plan states both roles.

### 4.2 Covariate inputs for healthy volunteers

| Covariate | In Model 1 | In Kovalenko 2016 | Input used for healthy volunteers | Flag |
|---|---|---|---|---|
| Body weight | Yes (Vc) | Yes (Vc) | Simulated from the weight distribution of each study or population (sections 5 and 7.1) | None |
| EASI | No | No | None (not an input of either model) | No patient EASI value is used; nothing to flag |
| Albumin | No | No | None | None |
| Race | No | No | None; ethnicity enters only through the weight distribution | None |
| BMI | No | No | Generated for reporting only (section 7.1); no effect on PK | None |
| Sex | No | No | Sex ratio 0.5 (assumption); used only for sex-specific weight bounds | None |

Kovalenko 2020 Model 1 was developed from rich-sampling studies in healthy subjects and in patients with atopic dermatitis (R668-AS-0907, TDU12265, PKM14161, R668-AD-1117) without a disease covariate. This is part of the applicability assessment (section 6.3).

### 4.3 Simulation settings

Individual parameters are computed in R from the published typical values, the subject's weight and the subject's random effects, and passed to the ordinary differential equation (ODE) solver of rxode2 5.1.7. The solver uses LSODA (Livermore Solver for Ordinary Differential Equations with automatic method switching), with absolute tolerance 1e-10 and relative tolerance 1e-8. No random numbers are drawn inside the solver.

## 5. Data for model evaluation

The external validation uses the original healthy-volunteer studies summarized in Li 2020. Only values already recorded in the repository are used, because the publication cannot be opened from the analysis environment.

| Study (dose) | Subjects per arm | Observed values in the repository | Weight distribution used |
|---|---|---|---|
| R668-AS-0907 SC (150 mg, 300 mg) | 6 | None for the single arms. The 300 mg arm is a component of the Clot 2021 non-Asian pooled 300 mg group (n 40: Cmax 32.0 +/- 10.5 mg/L, AUC0-last 544 +/- 218 mg x day/L) | Range 58.0 to 95.1 kg (range relayed by the sponsor from Li 2020; not checked against the publication) |
| TDU12265, Japanese (75, 150, 300, 600 mg) | 6 | 300 mg and 600 mg Cmax and AUC0-last mean and standard deviation (SD) (Clot 2021 Table 5); none at 75 and 150 mg | Range 52.1 to 78.3 kg (relayed by the sponsor) |
| PKM12350 (300 mg) | 15 test, 13 reference | AUC0-last mean and SD per arm (Li 2020 Table 3); GMR of the area to the last measurable time (AUC0-t, equivalent to AUC0-last) 0.90 (90% CI 0.71 to 1.16) from the FDA review | Fixed: normal 78 kg, SD 10 kg. Second judgement: uncertain study mean |
| PKM14161 (300 mg) | 19, 19 | AUC0-last mean and SD per arm | As PKM12350 |
| PKM14271 (200 mg) | 19, 19 | Cmax and AUC0-last mean and SD per arm | Fixed: normal 80 kg, SD 10 kg. Second judgement: uncertain study mean |
| R668-HV-1108 (300 mg) | 18, 18 | AUC0-last mean and SD per arm | As PKM12350 |

Where only a range is known, the weight distribution is a normal distribution with mean at the midpoint and SD equal to one quarter of the range, truncated to the range, for both sexes. Where nothing is known, the registered assumption B14 is used: normal, mean 78 kg (80 kg for PKM14271), SD 10 kg, with sex-specific bounds of 50 to 100 kg (men) and 40 to 90 kg (women). For the second judgement, each replicate study draws one study-level mean weight from a uniform distribution from the assumed mean minus 8 kg to the assumed mean plus 8 kg, shared by all arms of the study. The 8 kg half-width was proposed by the statistical review and has no source independent of the earlier results, so it does not decide the outcome. The lower limit of quantification (LLOQ) of these literature data is 0.078 mg/L (`config/design_clot2021.yaml`, the originator's assay). It is kept separate from the LLOQ of the planned study (`config/assay.yaml`, a placeholder; sections 7 and 8). Sampling schedules are the study-specific schedules in `config/design_clot2021.yaml`; the R668-AS-0907 schedule is used for both of its doses.

Items that were requested but are not in the repository are reported as not evaluable:

- the Li 2020 Table 2 values for R668-AS-0907 and for TDU12265 at 75 and 150 mg;
- Cmax for PKM12350, PKM14161 and R668-HV-1108;
- the Table 3 GMRs other than the PKM12350 value;
- the Table S1 weight distributions.

Model predictions are still produced for these arms, so that observed values can be compared when they are supplied.

## 6. Model credibility plan

### 6.1 Verification

1. Unit tests (testthat, `tests/testthat/test-m15.R`) for every new function, with tolerances in section 6.4:
   - the closed-form analysis models (with and without stratum) against `lm` and `confint`;
   - the cumulative-sum computation for nested sample sizes against direct fits;
   - the AUC0-last, Cmax and Tlast helper against the project NCA engine;
   - the partial AUC function against PKNCA 0.12.1 `pk.calc.aucint`;
   - the AUROC estimator against a direct pairwise computation, and the bootstrap index sharing;
   - the Tango score interval against a numerically maximized restricted likelihood;
   - the population sampler (sex-specific bounds, marginal distribution, weight-height correlation, common random numbers across populations);
   - the dense-sampling AUC0-last truth against a 0.01-day reference grid;
   - the concentration-scaling truth;
   - the extrapolated partial AUC against numerical integration.
2. Partial AUC verification on simulated data: the partial AUC function is compared with PKNCA 0.12.1 `pk.calc.aucint` (method lin up/log down, auc.type AUClast). The comparison uses 1,000 simulated profiles per schedule with visit-window jitter and residual error, after the standard BLQ rules. Agreement must be within a relative difference of 1e-8 where PKNCA returns a value. Where PKNCA returns no value (an interval that starts at or after Tlast), the function must return 0, and the count is reported (`results/m15/individual/pauc_pknca_check.csv`).
3. Legacy check of the individual-level engine. With the legacy trial population (normal 75 kg, SD 9 kg, truncated to 60 to 90 kg) and the confirmed schedule, three quantities must agree with the stored results within 3 standard errors (SE) of the difference: the share of subjects with a computable terminal elimination rate constant (lambda_z), and the median and 95th percentile of the extrapolated percentage. The stored results are `results/individual/individual_base.csv` (Kovalenko 2016) and `individual_struct2020.csv` (Model 1). The SE of the difference is sqrt(SE_new^2 + SE_stored^2), with SE_stored taken equal to SE_new because the stored runs used the same number of subjects, and SE_new from a 200-resample subject bootstrap. The new engine uses a different random-number layout, so bitwise agreement is not expected.
4. Identity check: for every simulated study, the pass result of section 8.1 at true GMR 1.00 with the 75 kg boundary must equal the pass result of the identity version of section 8.2. This is checked for both analysis models, all schedules and the pass counts at every sample size.
5. Analytic cross-check: the simulated pass probability at a true GMR of 1.00 and 0.95 (confirmed schedule, 75 kg boundary, 130 randomized per arm) is compared with the analytic joint power of `R/samplesize.R` at the CV of section 8.4.

### 6.2 Posterior predictive check (Li 2020 healthy-volunteer studies)

This section was headed "External validation: posterior predictive check" at registration; amendment 1 (section 16.5) renames it and the failure rule becomes the misfit rule. The rules, thresholds, units and judgements below are unchanged.

For each arm in section 5 and each PK model, 5,000 replicate studies are simulated under the fixed weight assumption. For studies without weight information, 5,000 more are simulated under the uncertain-mean weight assumption. Each replicate uses the study's number of subjects, dose and schedule, sex ratio 0.5, nominal sampling times, residual error and BLQ censoring at the LLOQ. Replicates 1 to 1,000 are also simulated with fixed mean weights of 70 and 86 kg (descriptive sensitivity). The replicates are analysed with the project NCA rules. The prediction interval (PI) is the 2.5th to 97.5th percentile over replicates (95% PI), except for the external-unit rule, which uses the 0.5th to 99.5th percentile (99% PI).

Rule units are the five independent studies with an observed AUC0-last at 300 mg or 600 mg: TDU12265 300 mg, TDU12265 600 mg, PKM12350, PKM14161 and R668-HV-1108. In a parallel-group study both arms received dupilumab, so the unit statistic pools both arms. The unit mean is the arithmetic mean over all subjects of the unit. The unit variability is the pooled within-arm CV, sqrt(sum (n_i - 1) SD_i^2 / (N - k)) / mean, computed the same way for observed and simulated units, so that a difference between the products of an observed study does not enter it. In this formula n_i and SD_i are the size and SD of arm i, N is the total number of subjects and k is the number of arms.

Pre-specified failure rule (per model):

- Mean rule (R1): the observed AUC0-last mean lies outside its 95% PI in 2 or more of the 5 units.
- Variability rule (R2): the observed pooled within-arm AUC0-last CV lies outside its 95% PI in 2 or more of the 5 units.
- External-unit rule (R4): the observed AUC0-last mean of any unit that is external to the model lies outside its 99% PI. The external units are PKM12350 and R668-HV-1108 for Model 1, and TDU12265 300 mg, TDU12265 600 mg, PKM12350 and PKM14161 for Kovalenko 2016.
- Cmax rule (R5): the observed Cmax mean lies outside its 95% PI in both TDU12265 units.

Reproduction failure for a model is R1 or R2 or R4 or R5 under the fixed weight assumption. The same rule is also judged under the uncertain-mean weight assumption. The outcome is reported as one of four categories: fails under both assumptions, fails under the fixed assumption only, fails under the uncertain-mean assumption only, or passes under both.

Under independent units and a correct model the false-alarm probabilities are about 2.3% (R1), 2.3% (R2), 2.0% for Model 1 and 3.9% for Kovalenko 2016 (R4), and 0.25% (R5). The false-alarm probability of each rule and of the combined rule is also estimated by treating each simulated replicate in turn as the observed data. The detection probability is estimated in the same way after multiplying every replicate's AUC0-last and Cmax means by 0.75, 0.85, 1.15 and 1.25 (CVs unchanged). Both are reported, and no threshold is changed after they are seen. Each observed percentile rank is reported with a Wilson 95% CI. A comparison is labelled borderline when that CI contains the 2.5th or 97.5th percentile.

The following are reported but are not part of the rule:

- the single arms;
- the non-Asian pooled 300 mg group. It contains the PKM12350 and PKM14161 test arms (34 of its 40 subjects), and the assumed split of 6, 15 and 19 subjects does not fit the stored SDs: the pooled total sum of squares, 1,853,436, is below the within-arm sums of squares of the two Li 2020 components, 2,199,434. Its replicate is assembled from the same replicate draws as its component arms;
- the parallel-group GMRs. The predictive distribution uses the same estimator as the observed value: each replicate's arm arithmetic mean and SD are converted to a geometric mean (GM) by GM = mean / sqrt(1 + CV^2), and GMR = GM_test / GM_reference. The pooled two-sample t estimate on subject logarithms is shown beside it and is compared with the reported PKM12350 GMR;
- the within-study pooled log variance, sum((n_i - 1) log(1 + CV_i^2)) / sum(n_i - 1);
- the 200 mg PKM14271 study, a different presentation with faster absorption;
- the arms without observed values.

The table reports, for each comparison: the observed value, the predicted median, the PI, the percentile rank and whether the study is internal or external for each model. Being internal means being part of the model's development data. R668-AS-0907 is internal to both models; TDU12265 and PKM14161 are internal to Model 1; R668-HV-1108 is internal to Kovalenko 2016; PKM12350 and PKM14271 are external to both.

Disclosure: the rule was written after the arm-level comparisons of the earlier qualification were known (results/step1/step1f_300mg_arm_check.csv and results/step1_k2020/step1f_300mg_arm_check.csv); those results are shown beside the new comparison. The rule was also revised after development runs of this analysis were seen (section 13.2).

Action: the section 6.2 table, the list of items that cannot be evaluated, a data request and the question on Model 1 as primary are sent to the requester before the population step and sections 8.1 to 8.6 start, whatever the outcome. If the primary model meets the failure rule, sections 8.1 to 8.6 are still run, and every later result is labelled as conditional on a model that did not pass the pre-specified external check. The requester may stop or redirect the work at any time, and the reply is recorded in the decision log. If only the sensitivity model passes, the model roles are not swapped; both models are reported as registered.

Data requested from the requester:

- per study and arm: n, Cmax and AUC0-last mean and SD, and the GMR with its 90% CI (Li 2020 Tables 2 and 3);
- the Table S1 weight mean, SD and range;
- the RSE or covariance matrix of both models;
- Korean sex-specific weight and height summaries;
- the Kovalenko 2020 Model 4 parameter table.

Values supplied after registration are added in a dated addendum. The failure rule set is not changed, and new arms are reported descriptively.

### 6.3 Applicability

The plan reports the weight range of the development data against the global and Korean population ranges. It notes the presence of Japanese subjects (TDU12265) in the Model 1 development data, which is the only ethnic information relevant to the Korean population, and the absence of a disease covariate. It also reports the known mis-predictions for presentations other than the study presentation (results/step1/step1b_quant_gate.csv; results/step1_k2020/step1b_quant_gate.csv):

- 200 mg presentation (1.14 mL of 175 mg/mL, with faster absorption than the model), simulated over observed AUC0-last: 0.82 to 0.91 for Kovalenko 2016 and 0.90 to 1.00 for Model 1 (under-prediction in 3 of 4 datasets);
- 200 mg presentation, simulated over observed Cmax: 0.92 to 1.01 and 0.86 to 0.94;
- study presentation (300 mg), Cmax over-predicted: 1.08 to 1.27 for Kovalenko 2016 and 1.00 to 1.17 for Model 1.

The 200 mg results are labelled "limited by the known AUC0-last under-prediction for the 200 mg presentation". Any enrolled weight outside the development-data weight range is reported as extrapolation.

### 6.4 Technical criteria

| Check | Metric | Threshold | Action if not met |
|---|---|---|---|
| Analysis models against lm and confint | absolute difference of estimate and CI | 1e-10 | Stop; fix before any run |
| Partial AUC against PKNCA | relative difference where PKNCA returns a value; 0 where PKNCA returns no value | 1e-8; all zeros | Stop; fix before any run |
| AUROC against direct pairwise computation | absolute difference | exact (1e-12) | Stop; fix before any run |
| Tango interval | restricted-likelihood score at the bounds | equal to plus or minus 1.96 within 1e-5 | Stop; fix before any run |
| Dense AUC0-last truth against a 0.01-day grid | relative difference of the AUC | 1e-4 | Stop; fix before any run |
| Concentration-scaling truth | log ratio against log g for AUC0-last | below 0.002 | Stop; fix before any run |
| Population sampler | Kolmogorov-Smirnov test against the truncated normal; correlation | p above 0.001; within 0.02 of 0.513 | Stop; fix before any run |
| Legacy check | difference against stored results | within 3 SE of the difference | Full trial runs do not start; cause reported |
| Identity check | per-study pass agreement | no mismatch | Full trial runs do not start; cause reported |
| Analytic cross-check (pilot) | simulated against analytic power | within max(3 Monte Carlo (MC) SE, 2.5 percentage points) | Full trial runs do not start; cause reported; any correction is a deviation |
| Monotonicity (pilot) | decrease between adjacent g values moving toward 1.00 | no more than 2 MC SE | As above |
| Applicability | enrolled weight outside the development-data range | reported | Labelled as extrapolation |

MC SE is the Monte Carlo standard error.

## 7. Simulation design

### 7.1 Populations

| Population | Weight distribution before inclusion | Height | Source |
|---|---|---|---|
| Global healthy-volunteer population | Normal, mean 76.55 kg, SD 9.275 kg, truncated to 58.0 to 95.1 kg | Normal, mean 172.5 cm, SD 10.4 cm | Weight: the R668-AS-0907 range (relayed by the sponsor). Height: the NCT03389893 summary of 71 adults with atopic dermatitis recorded in the repository, used for healthy volunteers (assumption) |
| Korean healthy-volunteer population | Normal, mean 65.2 kg, SD 6.55 kg, truncated to 52.1 to 78.3 kg | Normal, mean 167.0 cm, SD 8.68 cm (both sexes pooled) | Weight: the TDU12265 range, the request's fallback because KNHANES raw data are not available. Height: Size Korea 7th survey (2015), excerpt retrieved by web search on 2026-10-20 (mean height of men 173.9 cm in their 20s and 173.1 cm in their 30s; women 160.9 and 160.2 cm). The excerpt reports height means only, with no SD and no weight, so it cannot replace the weight range. The sex means 173.5 and 160.5 cm with assumed SDs of 6.0 and 5.5 cm are pooled at a sex ratio of 0.5, to match the sex-common weight distribution |

Weight and height are generated jointly by a Gaussian copula. The weight is drawn by inverse cumulative distribution function (CDF) sampling from its truncated distribution. The height (HT, cm) is drawn from its normal distribution conditional on the weight's normal score, with correlation 0.513 (the value derived in `scripts/53_atopic_summary.R` from the NCT03389893 summary statistics). BMI = WT / (HT / 100)^2 is reported but changes no result. Neither model has a BMI covariate, and the protocol's BMI 18 to 32 criterion is not applied as an inclusion filter; the share of enrolled draws outside 18 to 32 is reported. The sex ratio is 0.5 (assumption), and both sexes share the weight distribution before the sex-specific bounds. Inclusion criteria are men 60 to 90 kg and women 50 to 90 kg. The enrolled population is the population distribution truncated by these bounds: for the global population, men 60 to 90 kg and women 58.0 to 90 kg; for the Korean population, men 60 to 78.3 kg and women 52.1 to 78.3 kg. A weight-shape sensitivity uses uniform distributions over the same ranges (sections 8.1 and 8.6).

### 7.2 Design factors

| Factor | Levels |
|---|---|
| Population | Global, Korean |
| Stratification boundary | 75 kg; median of the enrolled population (fixed a priori from 200,000 draws, rounded to 0.1 kg) |
| Dose | 300 mg (primary), 200 mg (comparison) |
| Sampling schedule (days after the dose; study Day d is day d - 1 after the dose) | Confirmed schedule: 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56 (study Day 1, 2, 4, 6, 8, 11, 15, 22, 29, 36, 43, 50, 57). Reduced schedule: the confirmed schedule without study Days 36, 43 and 50. Extended schedule: the confirmed schedule with study Days 64 and 71 added |
| Randomized subjects per arm | 80 to 160 in steps of 10 |
| Withdrawal | 10%, missing completely at random (MCAR), independently per randomized subject; a withdrawn subject is not evaluable for any endpoint |
| Analysis model | Stratified analysis (treatment and weight stratum); unstratified analysis (pooled two-sample t) |
| PK model | Model 1 (primary), Kovalenko 2016 (sensitivity) |

Randomization uses stratified permuted blocks of two within each weight stratum in enrollment order. A study with n randomized subjects per arm consists of the first n subjects allocated to each arm, so every smaller study is nested in the larger one (common random numbers across n). Each simulated study enrolls 330 subjects, enough for 160 per arm. The base design for the mechanism and AUROC analyses is 130 randomized per arm (117 evaluable on average), the repository design. Expected share above 75 kg (pre-registration population check): about 55% in the global population and about 5% in the Korean population. The boundary comparison is therefore informative mainly in the Korean population. The number of evaluable subjects above the boundary is reported per study.

### 7.3 Data generation

For each subject the true concentration is solved once on the union of the three schedules. Sampling times are drawn uniformly within the visit windows of `config/trial_design.yaml`: plus or minus 2 hours up to day 1, 6 hours up to day 14 and 1 day afterwards (assumption). Time 0 is fixed, times are not negative, and the planned order is kept with a minimum gap of 1 hour. Observed concentrations are y = C x (1 + e_prop) + e_add with the model's residual SDs, censored as BLQ when y < LLOQ (0.078 mg/L, the planned-study placeholder) or at the pre-dose time. All schedules of a subject share the same jitter and residual draws, so schedule comparisons are paired.

### 7.4 Non-compartmental analysis

The project NCA engine (`R/nca.R`, Phoenix WinNonlin compatible, validated against NonCompart 0.8.4 and PKNCA 0.12.1) is used without change, with the rules of `config/nca_rules.yaml`:

- linear-up/log-down AUC on actual times;
- BLQ values before the first quantifiable value set to zero, intermediate BLQ values treated as missing, BLQ values after Tlast excluded, and values after two consecutive BLQ samples truncated;
- lambda_z by best fit over the last three or more points after Cmax (Cmax excluded), maximizing the adjusted R-squared with tie tolerance 1e-4 and the more-points rule, with positive slopes excluded;
- AUC0-inf from the observed last concentration (AUCINF_obs) = AUC0-last + Clast / lambda_z, where Clast is the last quantifiable concentration.

The percentage of AUC0-inf extrapolated beyond Tlast is AUC_%Extrap_obs. For the operating-characteristic analysis, which needs only AUC0-last, Cmax and Tlast, a vectorized helper with the same BLQ handling and integration is used; a unit test shows that it gives identical values. Partial AUCs are computed by a new function (section 8.3).

### 7.5 Analysis models and decision rule

Endpoints are log-transformed. N is the number of evaluable subjects. The stratified analysis is log y ~ treatment + weight stratum (two levels, the stratum used for allocation), with N - 3 degrees of freedom (df). The unstratified analysis is a pooled two-sample t test on log y, with N - 2 df. The 90% CI of the GMR is exp(estimate +/- t(0.95, df) x SE). Similarity is concluded when the 90% CIs for both AUC0-last and Cmax are within 0.80 to 1.25 (inclusive, unrounded). An endpoint that cannot be analysed counts as not concluded. If a stratum has no evaluable subject in a study, the stratified analysis drops the stratum term, and the number of such studies is reported.

### 7.6 Random numbers and common random numbers

The master seed is 20261020. Seeds are derived with the project rule (`R/seeds.R`, rule version 1):

- Studies in sections 8.1 to 8.3 use derive_seed(20261020, "trial", j, stream), with stream subj, alloc, drop, jitter or eps. The seed table is built for studies 1 to 10,000. A collision in the 28-bit seed space is resolved by appending "r1" (then "r2") to the later (study, stream) pair, and resolved pairs are logged.
- The external validation uses derive_seed(20261020, "val", arm, stream) for subjects and residuals, and derive_seed(20261020, "val", study, "wmean") for the study-level mean weight.
- The other analyses use the tags listed in the configuration file (seeds block).

No tag contains the PK model, population, boundary, dose, schedule, sample size, analysis model, true GMR or mechanism version. Within a study, every random number is indexed by enrollment position and planned time, never by arm: sex, weight uniform, height score, random-effect scores, allocation uniform, withdrawal uniform, and the jitter and residual draws on the union grid. Random effects are drawn as standard normal scores per parameter name and scaled by each model's IIV, so the two PK models also share subjects. These factors are therefore compared on the same subjects and residuals, and their differences are paired.

## 8. Analyses

### 8.1 Operating characteristic curves

Construction of the true GMR: for every test-arm subject the model concentration-time profile is multiplied by g before residual error and BLQ censoring; the additive residual is not scaled. The true GMR is then exactly g for AUC0-inf, Cmax and fixed-interval partial AUCs. For the dense-sampling AUC0-last truth (section 8.2(a)) it is g up to the shift of the LLOQ crossing (absolute difference of the log ratio and log g below 0.002, unit test); decisions are judged against g. This construction is mechanism-free by design; mechanism-specific differences are analysed in section 8.2.

Grid: 19 values of g equally spaced on the log scale from 0.80 to 1.25 (including 1.00), plus 0.95 and 1/0.95, so 21 values. For each PK model, population and dose, 1,000 studies are simulated. Every study is analysed for both boundaries, all three schedules, all nine sample sizes, both analysis models and all 21 values of g. Output: the probability of concluding similarity per cell with its Wilson 95% CI. The scenario table gives the power at g = 0.95 and 1.00. It gives the consumer risk in three ways at g = 0.80 and 1.25: the joint pass with both endpoints at the margin, the AUC0-last-only pass and the Cmax-only pass. The larger single-endpoint value approximates the largest consumer risk of the joint rule. Paired differences on the same studies are reported with the Tango score 95% CI and the discordant counts. The differences are: with versus without stratum, 75 kg versus median boundary, each schedule versus the confirmed schedule, global versus Korean population, and 200 mg versus 300 mg.

Consumer-risk precision: for the primary model, studies 1,001 to 10,000 are added at g = 0.80 and 1.25 for both populations and both doses, with all other factors. The consumer risk is then based on 10,000 studies (Wilson 95% half-width about 0.4 percentage points at 5%).

Weight-shape sensitivity: for the primary model at 300 mg and the confirmed schedule, the curves at g = 0.80, 0.95, 1.00 and 1.25 are repeated with uniform weight distributions over the global and Korean ranges (1,000 studies).

### 8.2 Mechanism-specific operating characteristic curves

One of F, ka, ke or Vmax is multiplied by 0.80, 0.85, 0.90, 0.95, 1.05, 1.10, 1.15 or 1.20 (32 versions plus the identity). The multiplier is applied to the test-arm subjects of each study. Every version uses the same test-arm subjects, jitter and residual draws as the identity version; the reference arm is solved once and is unchanged. Under Model 1 the ka multiplier leaves the mean transit time unchanged. For each PK model, population and dose, 1,000 studies are simulated with the 75 kg boundary, and the test arm is re-solved for every version.

(a) Induced true GMR: for each version the paired population GMR, exp(mean(log(X_test,i / X_reference,i))), is computed over 20,000 subjects of the enrolled population without residual error, where X is one of the following metrics:

- AUC0-inf: the model AUC to 400 days, extended to 4,000 days when the remaining amount exceeds 1e-7 of the dose;
- Cmax: the maximum of the true concentration on a 0.25-day grid to 60 days with three-point parabola refinement;
- the dense-sampling AUC0-last truth for each schedule: the model AUC from 0 to the earlier of two times. The first is the continuous true Tlast, the last time at which the true concentration is at or above the LLOQ, interpolated log-linearly on a 0.25-day grid to day 70. The second is the schedule's last planned time;
- the fixed-interval partial AUC truths, the model AUC between the interval boundaries.

(b) Probability of concluding similarity per version, schedule, sample size and analysis model, with Wilson 95% CIs; the main presentation uses 130 randomized subjects per arm and the stratified analysis.

Summary metric: the 50% crossing is the first point, moving outward from multiplier 1.00 (the identity is a node), at which the conclusion probability falls below 50%. It is interpolated linearly on the log multiplier. Its 95% CI comes from a bootstrap over studies (2,000 resamples). It is reported as a percent difference, 100 x (multiplier - 1), with the induced true AUC0-last, AUC0-inf and Cmax GMRs interpolated at that point. If the probability stays at or above 50% within the grid, the result is "not reached within the grid", together with the probability at the end of the grid. Comparisons are global versus Korean population, 300 mg versus 200 mg, and the three schedules.

Supplementary analysis beyond the requested range: for the primary model at 300 mg with the confirmed schedule, ke, Vmax and ka are also run at multipliers 0.60, 0.70, 1.30, 1.40 and 1.50 in both populations (130 randomized per arm, stratified analysis, 1,000 studies), so that the 50% crossing can be located.

### 8.3 Secondary discrimination: AUROC

This analysis is restricted to endpoints without an equivalence margin:

- AUC0-inf (AUCINF_obs). The primary set is the reliable-terminal-phase set, which excludes subjects with a reliability flag or without lambda_z. All subjects with lambda_z are used as a sensitivity.
- Partial AUCs from 0 to 14, 14 to 28 and 28 to 56 days after the dose (study Day 1 to 15, 15 to 29 and 29 to 57).

Partial AUC rule: linear-up/log-down integration on actual times after the standard BLQ rules. The concentration at an interval boundary is interpolated within its segment, linearly if rising and log-linearly if falling. There is no contribution after Tlast. A partial AUC of 0 is a valid result and is kept; subjects with no quantifiable sample are excluded. The share of zero values by arm, schedule and endpoint is reported. A cell where the zero share in the reference arm exceeds 50% is reported as not interpretable.

Statistic, at 130 randomized subjects per arm:

- AUC0-inf: the absolute estimated log GMR from the stratified analysis;
- partial AUCs: |log(arithmetic mean partial AUC (pAUC)_test / arithmetic mean pAUC_reference)|, with zeros kept.

Sensitivities for the partial AUCs:

- the absolute estimated log GMR from the stratified analysis fitted to log(pAUC + c), with c = 0.078 mg x day/L (the LLOQ times 1 day);
- lambda_z extrapolation for the 28 to 56 day interval: Clast x exp(-lambda_z (t - Tlast)) is added from the later of 28 days and Tlast to 56 days when lambda_z is computable, and both statistics are computed.

A study with fewer than 2 analysable subjects in an arm counts as missing for that endpoint. AUROC uses the studies in which both the positive and the negative are estimable, and the count is reported per cell.

Positives are the studies of section 8.2 with an injected difference of 5%, 10% or 15% (multipliers 0.95 and 1.05, 0.90 and 1.10, 0.85 and 1.15) for each mechanism. They are reported per direction and with both directions pooled per size. Negatives are the identity studies. AUROC is the Mann-Whitney estimate (ties count one half). The 95% CI comes from a percentile cluster bootstrap over studies (2,000 resamples; all versions of a resampled study are kept together). The CIs of all AUROC differences between schedules, endpoints, populations and doses come from the same bootstrap, because one set of resampled study indices is shared by every cell. The DeLong interval is shown for reference only, because positives and negatives share studies. Output: AUROC against difference size per mechanism, endpoint, population, dose and schedule.

### 8.4 Sampling evaluation

For each PK model, population and dose, 20,000 subjects are simulated on the three schedules (paired). The truths are:

- AUC0-last, primary: the dense-sampling truth of section 8.2(a), the same definition for both analyses;
- AUC0-last, secondary: the model AUC to the observed Tlast, labelled "integration error only";
- AUC0-inf: the model AUC to infinity, as in section 8.2;
- Cmax: the true maximum.

Bias is the median percent error, 100 x (estimate / truth - 1). Precision is the root mean square error (RMSE) of the percent error, together with the RMSE of the log ratio and the median absolute percent error, because the percent-error RMSE is unstable when lambda_z rests on few points. AUC0-inf is reported for all subjects with lambda_z and for the reliable-terminal-phase set.

Also reported:

- the distribution of AUC_%Extrap_obs (5th, 25th, 50th, 75th, 95th and 99th percentiles, and the share above 20%);
- the lambda_z success rate, defined as lambda_z computable and AUC_%Extrap_obs below 20%. It is given over all subjects and over subjects with a computable lambda_z, with Wilson 95% CIs;
- for the extended schedule, the shares of observed and true quantifiable samples at study Days 64 and 71;
- the share of subjects whose Tlast sample has a true concentration below the LLOQ.

The effect of the reduced and extended schedules on sections 8.1 to 8.3 is reported as paired differences against the confirmed schedule.

### 8.5 Sample size

1. Simulation-based power against randomized subjects per arm (80 to 160) at true GMRs of 0.95 and 1.00, for every scenario, from section 8.1.
2. Analytic joint power from `R/samplesize.R` (unchanged). The inputs are re-estimated from section 8.4 per PK model, population, dose and schedule: the CVs, the AUC0-last and Cmax correlation and the log-weight slopes. The log-weight variance components are computed from the simulated enrolled weights per boundary, and 90% of the randomized subjects are evaluable. The number of subjects needed comes from the analytic search (integer evaluable n; randomized = ceiling(n / 0.9)); the simulated power is interpolated between grid sizes only for display.
3. Assurance over CV uncertainty. The source is the AUC0-last CVs of the six Li 2020 300 mg arms: PKM12350 with 15 and 13 subjects, PKM14161 with 19 and 19, and R668-HV-1108 with 18 and 18 (df = 96).
   - Pooled log variance: s^2 = sum((n_i - 1) log(1 + CV_i^2)) / df.
   - Draws: sigma^2 = c x df x s^2 / chi-square(df). The calibration factor c is the true log variance divided by the mean of the arithmetic-CV-based pooled estimator, from 20,000 log-normal replicate sets at the six arms' sizes at the observed pooled value.
   - Transport to the enrolled population: sigma^2_trial = sigma^2 - beta_AUC^2 x (var_logWT_source - var_logWT_enrolled). Here beta_AUC is the model's log-weight slope. The source variance var_logWT_source is the mean over 2,000 simulated studies (study mean weight uniform from 70 to 86 kg, SD 10 kg, default bounds, 100 subjects each) of the within-study variance of log weight.
   - Cmax SD: the AUC0-last SD times the model's ratio of Cmax SD to AUC0-last SD, taken from section 8.4 for the same model, population, dose and schedule.
   - Power is evaluated on a CV grid from 20% to 80% in steps of 0.25 points and linearly interpolated at each draw. Draws outside the grid are clamped, and their share is reported. Assurance is the mean power over 10,000 draws.
   - It is reported for every randomized n per arm from 80 to 160, at true GMRs of 0.95 and 1.00, for both analysis models, each PK model, population, dose, schedule and boundary. The difference between the Li 2020 weight distributions and the enrolled population is a limitation, partly addressed by the transport.
4. Assurance over parameter uncertainty. The method is registered: log-normal sampling from RSEs (independent), 1,000 parameter sets of 2,000 subjects each, with the CVs and slopes re-estimated per set and power averaged. Its status is pending, because no RSE or covariance is recorded in `config/param_uncertainty.yaml`. Amendment 1 (section 16.8) implements it with the Kovalenko 2016 bootstrap percentiles supplied on 2026-10-21 (Model 1: borrowed Kovalenko 2016 relative widths until FDA 2017 review Table 4.6.6 is supplied), with the departures stated there.
5. Protocol assumption. The comparison is made at a CV of 43% (Syneos proposal, `config/prereg_20260926.yaml`) and 117 evaluable (130 randomized) subjects per arm (`config/trial_design.yaml`). It reports the power at true GMRs of 0.95 and 1.00 and the subjects per arm needed for 85% and 90% power. These are set side by side with the same quantities at several other CVs:
   - the observed pooled CV of the Li 2020 arms, raw and calibrated and transported;
   - the model-predicted CV of those arms (section 6.2);
   - the model CV of the enrolled population.

   Both Cmax rules are shown: the ratio rule, and the fixed rule, which uses the model Cmax SD and is the rule of the stored protocol figures.

### 8.6 Body-weight effect

The weight effect is simulated directly through the model covariate. No derived concept (such as an effective dose) is used. For each population (global and Korean, and their uniform-shape variants as a sensitivity), the 5th, 25th, 50th, 75th and 95th percentiles of the enrolled weight distribution are taken from 200,000 draws. For each PK model, dose (300 mg and 200 mg) and percentile weight, the metrics X are AUC0-inf, Cmax and the dense-sampling AUC0-last truth for the confirmed schedule. The following are reported:

1. Point estimate: the metric of the typical subject (all random effects zero) at the percentile weight, divided by the typical subject at 75 kg.
2. 90% interval of the ratio: the 5th to 95th percentile over 20,000 subjects of the paired individual ratio X_i(w) / X_i(75 kg), with the same random effects at both weights. This is the between-subject spread of the weight effect itself, caused by the nonlinear elimination. The geometric mean of the paired ratios is also given.
3. A separate, labelled column: the 5th to 95th percentile of X_i(w) / typical X(75 kg), which includes between-subject variability. The parameter-uncertainty interval is pending (section 9). Amendment 1 (section 16.8) computes the parameter-uncertainty interval.
4. Normalized sensitivity coefficients with respect to F, ke and Vmax for the typical subject at each percentile weight: (ln X(1.01 theta) - ln X(0.99 theta)) / (ln 1.01 - ln 0.99), with Cmax from a 0.01-day grid.

The forest plot shows the ratios with their intervals on a log scale.

## 9. Uncertainty

| Source | Handling |
|---|---|
| Between-subject variability and residual error | Simulated (published IIV and residual SDs) |
| Monte Carlo error | Wilson 95% CIs for proportions; Tango score CIs for paired differences on common random numbers; cluster bootstrap CIs for AUROC and crossings; replicate numbers in section 8 |
| Parameter (estimation) uncertainty | Not propagated. The request asks for propagation with the existing Monte Carlo and bootstrap, but the repository has no RSE, covariance matrix or bootstrap parameter sets (`config/param_uncertainty.yaml` status pending). The only bootstrap in the repository resamples subjects for a schedule criterion, not parameters. The registered method (section 8.5 item 4, and the same parameter sets for the interval in section 8.6) runs once the RSEs are entered and committed. All results are conditional on the published point estimates. Amendment 1 (section 16.8) propagates it from the Kovalenko 2016 bootstrap percentiles supplied on 2026-10-21 (Model 1: borrowed Kovalenko 2016 relative widths until FDA 2017 review Table 4.6.6 is supplied). |
| Structural uncertainty | Every main result is repeated with Kovalenko 2016 |
| Population assumptions | Shapes derived from ranges (normal with SD of one quarter of the range); uniform-shape sensitivity in sections 8.1 and 8.6 |
| Validation weight assumptions | Fixed assumption (decides), uncertain-mean assumption (second judgement), fixed 70 and 86 kg (descriptive) |
| Assay and design placeholders | LLOQ 0.078 mg/L, visit windows, Day 1 post-dose time 0.25 day, sex ratio 0.5 (section 12) |

## 10. Outputs

Figures, each paired with a table of the same values and written in English (with Korean companions for internal use):

1. predicted against observed for section 6.2;
2. OC curves for section 8.1;
3. mechanism OC curves for section 8.2;
4. AUROC against difference size for section 8.3;
5. bias and precision box plots for section 8.4;
6. power against sample size for section 8.5;
7. a forest plot for section 8.6.

Deliverables:

- this plan;
- one results summary table. Its rows are:
  - the external validation outcome;
  - for the primary scenario (global population, 300 mg, confirmed schedule, 75 kg boundary, 130 randomized per arm, stratified analysis): the power at true GMRs of 0.95 and 1.00, and the joint and single-endpoint consumer risk at 0.80 and 1.25;
  - the effects of the stratum term, the boundary, the population, the dose and the schedule changes;
  - the 50% crossing per mechanism;
  - AUROC at 10% differences;
  - the sample size for 85% and 90% power and the assurance;
  - the comparison with the protocol assumption;
  - the weight effect at the 5th and 95th percentiles.

  Its columns are the primary model, the sensitivity model and the source file;
- at most five lines of key conclusions, generated from the result files with premise checks, using the readings of section 2.1;
- the list of all assumptions in the format of `regulatory/tables/assumptions_register.csv` (section 12).

Every number is read from a result file.

## 11. Computation, pilot and records

Hardware: four cores, 15 GB memory. Order of work:

1. section 6.2 and its report to the requester;
2. population quantiles and boundaries;
3. truths;
4. the legacy check;
5. section 8.4;
6. the partial AUC verification;
7. the pilot: studies 1 to 500 of sections 8.1 and 8.2 for the primary model, global population, 300 mg, and studies 1 to 50 for the primary model, Korean population, 200 mg;
8. the pilot checks (section 6.4);
9. the full runs of sections 8.1 and 8.2;
10. the supplementary mechanism run;
11. the consumer-risk extension;
12. the weight-shape sensitivity;
13. section 8.6;
14. section 8.5;
15. the summary.

The pilot studies are part of the final run, so the pilot changes no design element; any change after the pilot is a deviation (section 13). The pilot also reports the run time, the number of studies in which the stratum term was dropped, and the per-endpoint counts of missing values.

Estimated cost from single-core timings of the analysis scripts is about 28 core-hours, or 7 to 9 hours of wall time on four cores. If the pilot projects more than 10 hours, every Model 1 cell runs first, and the Kovalenko 2016 mechanism runs at 200 mg run last.

Every run writes a log (`R/logging.R`) with the seed rule, the git commit, the list of uncommitted changes and the package versions. Logs, configuration and result files are committed, and the software environment is recorded in `renv.lock`.

## 12. Assumptions

| ID | Assumption | Value used | Basis | Impact evidence | Action before submission |
|---|---|---|---|---|---|
| B01 | Bioanalytical LLOQ of the planned study | 0.078 mg/L | Originator's assay (placeholder, `config/assay.yaml`) | Earlier LLOQ sensitivity (assumptions register B01) | Set the validated assay LLOQ and regenerate |
| B02 | Day 1 post-dose sampling time | 0.25 day | Placeholder | Affects Cmax sampling only | Set to the protocol time |
| B03 | Visit windows | plus or minus 2 h to Day 1, 6 h to Day 14, 1 day thereafter | Placeholder | Simulated in all trial-level analyses | Align with the protocol |
| B04 | Body weight distribution of enrolled subjects (replaced for this analysis) | Global: normal 76.55 kg, SD 9.275, truncated 58.0 to 95.1 kg; Korean: normal 65.2 kg, SD 6.55, truncated 52.1 to 78.3 kg; inclusion men 60 to 90 kg, women 50 to 90 kg | Ranges relayed by the sponsor from Li 2020; shape rule of section 5 | Uniform-shape sensitivity (sections 8.1, 8.6) | Confirm ranges against the publication; supply Korean data |
| B05 | Randomization stratification (replaced) | Two strata at 75 kg or at the enrolled median; permuted blocks of two | Request | Both boundaries simulated | Align with the protocol |
| B06 | Evaluable subjects (replaced) | 80 to 160 randomized per arm; 10% MCAR withdrawal; withdrawn subjects not evaluable for any endpoint | Request | Sample-size grid | Align with the protocol |
| B07 | Km fixed, without between-subject variability | 0.01 mg/L in both models | As published | Earlier sensitivity (register B07) | None |
| B08 | Parameter uncertainty (replaced by amendment 1, corrected on 2026-10-21) | Kovalenko 2016: standard errors from the Table 2 bootstrap 5th and 95th percentiles (BLQ-included column; log scale (ln P95 - ln P5)/3.29 for positive parameters, (P95 - P5)/3.29 for the exponent and the variances); Model 1: borrowed uncertainty, the Kovalenko 2016 relative widths applied to Vc, ke, kcp, kpc, ka, Vm, F and the weight exponent until its own bootstrap intervals (FDA 2017 review Table 4.6.6) are supplied; independence | Values supplied on 2026-10-21; Model 1 intervals to be supplied | Assurance and weight-ratio intervals over 1,000 parameter sets (section 16.8) | Supply FDA 2017 review Table 4.6.6 |
| B09 | Residual error model, down to the LLOQ | As published (proportional 15.0% and 24.2%, additive 0.03 mg/L) | As published | Not evaluated here | None |
| B10 | Immunogenicity | Not modelled | No anti-drug antibody (ADA) model in the published PK models | Not evaluated here | None |
| B11 | Automated lambda_z selection without manual review | Phoenix best-fit rules | Standard NCA practice | AUC0-last and Cmax unaffected | Pre-specify in the statistical analysis plan |
| B12 | Same presentation for both products; 200 mg presentation modelled as the 300 mg model | Model parameters unchanged | Study design; no presentation-specific absorption | Section 6.3 limitation | Label 200 mg results as limited |
| B13 | Population: healthy adults (replaced) | BMI 18 to 32 criterion not applied as a filter | BMI changes no result (section 7.1) | Share outside 18 to 32 reported | None |
| B14 | Body weight of Li 2020 arms without weight information | Normal 78 kg (80 kg for PKM14271), SD 10 kg, bounds 50 to 100 kg (men), 40 to 90 kg (women) | Weights not recorded | Second judgement with uncertain study mean; fixed 70 and 86 kg; pooled ratio reported with and without PKM12350 | Supply Table S1 |
| B18 | Primary model | Kovalenko 2020 Table 1 Model 1 (initial base structural model of FDA 2017 review section 4.6) | Section 4.1; confirmed by the requester on 2026-10-21 (section 16.14) | Kovalenko 2016 repeats every main result | None |
| B19 | Covariates | Weight on Vc only; no EASI, albumin, race or BMI effect | Model structure of Model 1 and Kovalenko 2016 | Model 3 and Model 4 only in the optional P3 comparison (section 16.10) | None |
| B20 | Height and weight-height relation | Global: normal 172.5 cm, SD 10.4 (atopic dermatitis summary used for healthy adults); Korean: normal 167.0 cm, SD 8.68 (Size Korea means, assumed SDs); Gaussian copula correlation 0.513 | Section 7.1 | BMI only; no effect on PK | None |
| B21 | Sex ratio and sex-common weight shape | 0.5; same weight distribution for both sexes before sex-specific bounds | Repository assumption | Not evaluated | Confirm expected enrolment |
| B22 | Median boundary | Fixed a priori from 200,000 draws | Section 2.3 | Both boundaries simulated | None |
| B23 | Jitter | Uniform within each visit window | Repository rule | Simulated | Align with the protocol |
| B24 | True GMR in section 8.1 | Concentration scaling by g; additive residual not scaled | Section 8.1 | Unit test of the scaling truth | None |
| B25 | Representation of a product difference | One parameter multiplied in the test arm; Model 1 ka multiplier with MTT unchanged; Vm read as Vmax | Section 8.2 | Four mechanisms reported | None |
| B26 | Variability of test and reference | Identical IIV and residual error; diagonal IIV | Published models | Not evaluated | None |
| B27 | Li 2020 arm sizes and TDU12265 values | Sizes from the request; TDU12265 300 and 600 mg values from Clot 2021 Table 5, assumed equal to Li 2020 Table 2; TDU12265 600 mg n = 6 | Section 2.2 | Not evaluated | Confirm against the publication |
| B28 | Range-only weight distributions | Normal with midpoint mean and SD of one quarter of the range, both sexes within the range | Section 5 | Not evaluated | Supply Table S1 |
| B29 | Non-Asian pooled 300 mg composition | 6, 15 and 19 subjects (inconsistent with the stored SDs; descriptive only) | `config/design_clot2021.yaml` | Excluded from the rule | Confirm composition |
| B30 | R668-AS-0907 schedule for 150 mg | Same as 300 mg | No dose-specific schedule recorded | Not evaluated | Confirm |
| B31 | LLOQ of the literature data | 0.078 mg/L | `config/design_clot2021.yaml` | Not evaluated | None |
| B32 | Truth definitions | Dense-sampling AUC0-last truth (section 8.2(a)); AUC0-inf to 400 or 4,000 days | Section 8 | Secondary integration-only truth reported | None |
| B33 | Partial AUC handling | No contribution after Tlast; zeros kept; arithmetic-mean ratio statistic; offset 0.078 mg x day/L in the sensitivity | Section 8.3 | Two sensitivities | None |
| B34 | CV source for assurance | Li 2020 arm CVs, calibrated and transported; Cmax SD by the model ratio | Section 8.5 | Raw and transported values reported | Supply study-specific variability |
| B35 | Observed GMR of the Li 2020 parallel studies | Log-normal conversion of arm means and SDs (descriptive) | Section 6.2 | Pooled t shown beside it | Supply reported GMRs |
| B36 | Interpretations of the request | Section 2.3 | Request wording | Not applicable | Confirm with the requester |
| B37 | Race effect form in the Model 4 configuration (optional P3 comparison only) | ke x (1 + 0.123 x white), non-white reference | Coefficient supplied; form not recorded | Exponential form differs by less than 1% | Supply the equation if P3 is requested |
| B38 | EASI, albumin and BMI in a Model 4 healthy-volunteer mapping (optional P3 comparison only) | Factor 1 (EASI effect neutral; albumin and BMI at the reference) | Equations and reference values not recorded | Not run | Supply the equations if P3 is requested |
| B39 | Reference weight of the Model 3 and Model 4 configurations (optional P3 comparison only) | 75 kg | Not supplied | Not run | Supply if P3 is requested |
| B40 | Variability of the Model 3 and Model 4 configurations (optional P3 comparison only) | Model 1 between-subject variability and residual error | Phase 3 variability not supplied (D-029: variability on Vc and ke only) | Not run | Supply if P3 is requested |
| B41 | BMI reference for the Model 4 weight effect (optional P3 comparison only) | 26 kg/m2 (register B15) | Phase 3 mean BMI 25.4 to 27.3 | At a fixed height the ratio does not depend on it | None |
| B42 | Visit windows of anti-drug-antibody samples | As the PK windows (B03); elapsed day 84 at least 1 hour after day 70 | No protocol windows | Not evaluated | Align with the protocol |
| B43 | Between-study exposure variation in the posterior predictive check | Log-normal factor per study, CV 10%, median 1; used only when between-arm variation is found (section 16.7 item 3) | Request of 2026-10-21 and its correction | False-positive and detection probabilities with and without the factor, if run | None |
| B44 | NCA convention of the Li 2020 values | Unknown; the registered project rule decides, the request's rule (post-dose BLQ = 0) is reported beside it | Publications not accessible | Both rules reported | Confirm the convention |
| B45 | Product of each Li 2020 arm | Recorded only for the PKM14161 test arm (product #7, C2P1 prefilled syringe); assay AV-09095 for PKM12350 and R668-HV-1108, AV-13074 for PKM14161 and PKM14271 | Requester, from FDA 2017 review Tables 4.1.1.a, 4.2.a and 4.2.b (review not accessible) | Descriptive columns only | Supply the product of each arm |

Register IDs B01 to B14 reuse the register in `regulatory/tables/assumptions_register.csv`; B04, B05, B06 and B13 are replaced for this analysis only; B15 to B17 of the register are not used here (B41 cites B15); new rows start at B18; B37 to B45 were added by amendment 1 (section 16.10); B37 to B41 apply to the optional P3 comparison only (correction of 2026-10-21).

## 13. Deviations, post hoc analyses and disclosures

### 13.1 Rule for changes

Any change after registration is recorded in the decision log with the reason and in a deviations table of the results summary. Analyses not described here are labelled post hoc.

### 13.2 Known before registration

1. The arm-level comparisons of the earlier qualification of both models were known when the failure rule of section 6.2 was written, including the Model 1 misses on PKM12350 and R668-HV-1108 (section 4.1).
2. Development runs of the section 6.2 code, whose verdicts were seen:
   - 20 replicates at 12:25 and 12:38 Coordinated Universal Time (UTC) on 2026-10-20, under two earlier rule versions. The first version had nine arms and a GMR rule. The second, revised between these two runs, had eight arms, external-arm and Cmax rules, and the pooled group excluded.
   - 60 replicates at 12:46 UTC under the five-unit rule as proposed by the statistical review, with the uncertain-mean weight assumption deciding.
   - 40 replicates in a scratch copy of the repository at 13:00 UTC, under the same rule, while testing the code; its verdict lines were seen.
   - 60 replicates in the scratch copy at 13:10 UTC under the registered rule, to test the revised code; only the structure of the output tables was inspected.
3. Changes to the rule after the 12:46 run, made in response to the consistency review:
   - the external-unit rule was changed from "2 or more external units outside the 95% interval", which is implied by the mean rule, to "any external unit outside its 99% interval";
   - the variability rule was changed to the pooled within-arm CV;
   - the fixed weight assumption B14 was made the deciding assumption, because the uncertain-mean assumption had no independent source.
4. The population step was run (enrolled medians 76.0 kg for the global and 66.1 kg for the Korean population; shares above 75 kg about 55% and 5%).
5. The truth, sampling-evaluation, weight-effect, trial and sample-size scripts were run with 300 to 600 subjects and 4 to 8 studies in the scratch copy, to test the code. These outputs are not results, are not stored, and no design element was changed in response to them. The truth grid was refined from 0.5 to 0.25 day in response to a unit test against a fine reference grid, not to an outcome.

### 13.3 Deviations after registration

1. Seed lookup error in the trial script (found by the registered pilot check). In `scripts/73_m15_trials.R`, the per-study seed was looked up in a data.table whose study-index column had the same name as the function argument, so every simulated study used the same seed and the 500 pilot studies were identical. The analytic cross-check of section 6.4 failed (simulated power 100% against an analytic 94% at a true GMR of 0.95), and the full run did not start. The lookup was replaced by a matrix lookup with a unit test, a static check for this pattern was added to the test suite, and the defective pilot outputs were deleted before they were committed (the pilot run logs and the failed check lines in `results/m15/run_progress.log` are kept as the record, and the outputs can be regenerated from the code at commit 3acf813) and the pilot was rerun from the start. No design element, threshold or seed rule was changed. The same review found an argument and column name clash in a filter of the summary script, which had not yet produced any output; it was fixed before use.

2. Pre-run code audit. A multi-agent code audit before the full run (independent recomputation of one OC study and one mechanism study: no mismatch) led to changes in output handling only:
   - the trial script now stops when a parallel worker returns no result and checks that every batch file holds exactly its studies. The change was made by atomic file replacement after the first three full-run stages had started (the Model 1 operating characteristic runs for the global population at 300 and 200 mg, which were complete, and for the Korean population at 300 mg, which was running; all three used the previous version of the script); computed values are unchanged, and the summary script checks every cell for missing or duplicated studies;
   - the summary script censors bootstrap resamples that do not reach the 50% crossing beyond the grid end instead of dropping them;
   - it matches bootstrap indices to study numbers;
   - it reads empty fields as missing;
   - it reports the test-arm zero share and the interpretability flag on AUROC differences;
   - it labels the summary table comma-separated values (CSV) file as conditional;
   - it reports the required sample size for every cell.
3. Run order. The full run follows the order in `scripts/run_m15.sh` (Model 1 operating characteristic and mechanism runs, the supplementary runs, Kovalenko 2016 operating characteristic runs at 300 and 200 mg and mechanism runs at 300 mg, the consumer-risk extension, the weight-shape sensitivity, Kovalenko 2016 mechanism runs at 200 mg), which differs from the order listed in the configuration file. Results do not depend on the order, because every study has its own seeds.
4. Post-run correction. The Wilson 95% CI of the lambda_z success rate among subjects with a computable lambda_z (section 8.4) was missing from the sampling-evaluation table; it was added from the stored subject-level results without re-simulation.
5. Figures. Korean-labelled companion figures are not produced; the Korean summary text is the internal companion.
6. Post-run verification (before the post stage). A multi-agent verification of the preview outputs recomputed 157 reported numbers independently from the raw files (operating characteristics, Tango intervals, crossings, AUROC, calibration, assurance, weight ratios, bias): all matched. It led to changes in reporting and output handling only; no computed value changed:
   - the completeness check now requires every cell to hold every study;
   - the paired-difference table carries the randomized n per arm. Paired differences are reported at 130 randomized per arm and true GMRs of 0.80, 0.95, 1.00 and 1.25 only, which section 8.1 did not state;
   - in differences between 50% crossings, a lower-side resample that does not reach 50% is censored at minus infinity (it was 0), so a bound that falls on censored resamples is reported as not determinable instead of as a finite value;
   - AUROC rows report the test-arm zero share of the positive versions, and the identity-version column is named as such;
   - the pilot report counts studies, not study-version rows, with a missing estimate;
   - the figure 5 subtitle gives the subject counts from the data (the AUC0-inf boxes use the reliable-terminal-phase subset);
   - every result table written by the summary, sample-size and weight steps carries a conditionality column, and a notice file covers the other result folders;
   - the summary table shows the protocol comparison under both Cmax rules with the sample sizes, plain labels instead of internal codes, and an abbreviation line; the key conclusions have an abbreviation line, and the last clause of conclusion 5 also requires the registered design-difference reading;
   - the English text check covers every .csv and .md file under `results/m15` except the Korean summary; the Hangul check also covers Jamo;
   - in `results/m15/individual/bias_precision.csv`, the endpoint codes AUCinf_ruleA and AUCinf_computable denote the sets coded AUCinf_A and AUCinf_B in Appendix B.
7. Editorial. Abbreviations that this plan used without expansion are now expanded at first use and listed in Appendix A.
8. Reporting choices after the Model 1 preview (before the post stage; no computed value changed). Quality checks of a preview made from the Model 1 results changed the layouts of figures 1 to 6 and gave conclusion 3 readable design names and the direction and range of the qualifying paired differences. They also moved the parameter-uncertainty caveat to the header of the key conclusions and of the summary table, and replaced conclusion 5 by a contrast of the reduced schedule's power gain with its measurement quality (model CV of AUC0-last, error against the dense truth, share of reliable AUC0-inf, AUROC). This choice was made after seeing in the preview that the gain comes from a shorter AUC0-last; its numbers and the condition of its last clause are generated from the result files.
9. Summary script and run driver. `scripts/75_m15_summary.R` and `scripts/run_m15.sh` were not in the registration commit (section 14 lists them as analysis code). They were added at 13:32 UTC on 2026-10-20, before the section 6.2 results and before any run of sections 7 and 8, and were later revised as described in items 1, 2, 6 and 8 and section 13.4.
10. Final verification (after the post stage). A second multi-agent verification of the final outputs (290 checks, including the Kovalenko 2016 results, the 10,000-study consumer-risk extension, every number of the summary documents, the figures, the documents and reproducibility) found no wrong number. It led to these changes, none of which changes a computed value:
   - an applicability table (section 6.3): the full weight range of the development data is not recorded in the repository, so the enrolled weights are compared with the weight ranges of the internal studies whose ranges are known (R668-AS-0907, internal to both models; TDU12265, internal to Model 1), and the share outside is reported as extrapolation;
   - the label for 200 mg results required by section 6.3, in a limitation column of every result table with 200 mg rows, in the captions of figures 2, 3, 5 and 6, in the summary table and in conclusion 3;
   - a deviations and post hoc table in the results summary (section 13.1), generated from this section;
   - the AUROC labels state that both directions are pooled; figures 3 and 6 name the primary model; the figure 7 axis title is wrapped; the Korean summary states its sources;
   - corrections to items 1, 2 and 3 of this section (pilot outputs, number of stages, run order);
   - the run log records an unknown code state when git is not available;
   - the post stage was rerun through the driver at a clean commit, and the outputs were compared with the previous ones.

11. Computations before the amendment (amendment 1, D-069). On 2026-10-21, before amendment 1 was committed, read-only mapping agents computed values from stored files and scratch copies, contrary to the request's rule that the amendment is committed before any new computation. The values and the design elements chosen after seeing them are listed in section 16.11; the registered runs regenerate every reported value.
12. Readings and labels (amendment 1, D-069). Section 2.1 item 3 (the detection-probability wording) is suspended until section 16.7 is complete and is then replaced by the fixed template of section 16.7 item 7; "External validation" is renamed "posterior predictive check (Li 2020 healthy-volunteer studies)", the failure rule the misfit rule and rule R4 the rule on out-of-sample units, with no change to any rule, unit or judgement; enrolled weights outside the documented ranges are labelled "outside the documented range" instead of "extrapolation" (sections 6.3, 6.4 and 16.9); sampling conclusions follow the rule of section 16.3 item 4 (coverage on the Kovalenko 2016 basis with Model 1 in parentheses, as corrected on 2026-10-21; other levels from the less favourable model) instead of section 2.1 item 1, as requested.
13. Parameter-uncertainty method (amendment 1, D-069). Section 16.8 departs from section 8.5 item 4: standard errors come from bootstrap percentiles or are borrowed; the weight exponent and the variances are drawn on the linear scale; draws with a non-positive variance or with F of at least 1 are redrawn; the same 2,000 subjects are used for every set, with a control-variate correction.
14. Model 3 and Model 4 configurations and code (amendment 1, D-069). The Model 3 and Model 4 configurations (optional P3 comparison only, section 16.10) borrow the Model 1 variability, not the phase 3 variability the request named (register B40). `scripts/77_m15_amend.R` is written after the amendment commit and committed before its first run (section 16.12).

15. Correction of the amendment (amendment 1, D-070). On 2026-10-21 the requester corrected the first request: Model 1 is Kovalenko 2020 Table 1 Model 1 as registered, the supplied standard errors belong to Model 4 and the "base model" values to Model 3. Section 16 was revised before any further computation (section 16.14): in-sample labels of Model 1 follow the FDA 2017 review; the healthy-volunteer mappings are deleted and the Model 3 and Model 4 comparison becomes an optional P3 item; the between-study sensitivity becomes conditional on a between-arm variation check; Model 1 parameter uncertainty borrows the Kovalenko 2016 relative widths instead of the Model 4 standard errors; the applicability range of 42.3 to 175.4 kg and the EASI dimension are withdrawn; coverage statements quote Kovalenko 2016. The between-arm criterion, the systematic-bias reading and the Kovalenko 2016 coverage basis were registered after the stored-replicate arm table and the six-arm pooled ratio had been seen (section 16.11); mode `ppcstored` had run once before the correction, during the `mechsuppx` run (section 16.12). The relabelling in `scripts/75_m15_summary.R` (figure 1 out-of-sample labels, decision-log column of the deviations table) is part of the correction commit. The weight check met its registered tolerance of 1e-6 in 232 of 252 rows; 20 rows of Model 1 200 mg Cmax differ by at most 2.0e-6 because the stored values were solved with the registered truth tolerances (D-042, subject-level differences up to 5e-6); the criterion is not changed and the rows are reported as not meeting it, as agreed by the requester. During the runs, script 77 was corrected before use of any affected output: the identity check counted the 18-row versions file whole (the data matched exactly), the weight check first used a 4,000-day solver step that underflows for transit-structure models and a finer Cmax grid than registered, and the seven-arm table first required weight-mode replicates for the range-only R668-AS-0907 arm. The first run of mode `p1report` stopped on a formatting error (a vector condition in the robustness text); it and the following findings of the review of the correction were fixed before any P1 output was used (D-070): conclusion 3 named Model 1 as the higher anti-drug-antibody concentration on every day and printed Day 57 and Day 85 medians as 0.00 (now the higher model per day from the result file, small values with three significant digits); the appendix coverage sentence gave a range over both models (now Kovalenko 2016 with Model 1 in parentheses); the weight exponents and the weight elasticity were fixed numbers (now from the result files); the conclusion 1 template slot carried the words "over-prediction"; abbreviations of the P1 report were not expanded; the P1 table had no assay column; the parameter-set redraw stopped after 20 instead of the registered 10 redraws; the Model 1 rows of the final table lacked the "borrowed uncertainty" label; the posterior predictive check still wrote a diagnostic-mapping rule row.

### 13.4 Post hoc analyses

1. Assurance with the CV shifted to each cell (registered on 2026-10-20 during the full run, before its computation in the post stage). The assurance of section 8.5 item 3 draws its CV from the Li 2020 300 mg arms for every dose and schedule, so it carries no model effect of dose or schedule on variability. This was noticed when a preview of figure 6, made from the Model 1 results available at the time, showed the 200 mg assurance above the power at the model CV. The post hoc variant adds the model difference in log variance between the cell and a reference cell (same PK model and population, 300 mg, confirmed schedule): sigma^2_cell = sigma^2 - beta_AUC,ref^2 x (var_logWT_source - var_logWT_enrolled) + [ln(1 + CV_model,cell^2) - ln(1 + CV_model,ref^2)]. It uses the same draws, grid and interpolation, and equals the registered assurance in the reference cells (checked in the script). The registered assurance is reported unchanged, and the summary table, which uses the reference cell, is not affected.

2. Amendment 1 analyses (amendment 1, D-069). Every analysis of section 16 is post hoc relative to the registration of 2026-10-20; each was registered in section 16 and the configuration block amendment_20261021 before its computation, except as stated in section 13.3 items 11 and 15 (section 16.11) and section 13.4 item 3.
3. Calibration of the between-arm variation criterion (amendment 1, D-070). The share of stored replicates with a chi-square p-value of Cochran's Q below 0.10 is reported beside the decision of section 16.7 item 3. A reviewer of the correction computed it from the stored replicates before it was registered (about 16% for each model) and after the observed Q (p = 0.36) had been computed; it does not change the decision, which stays as registered.

## 14. Software and quality control

Software: R 4.3.3; rxode2 5.1.7; data.table 1.14.10; ggplot2 3.4.4 and patchwork 1.2.0 (figures); the project NCA engine; PKNCA 0.12.1 (verification of the partial AUC only). No new package is added; AUROC and the Tango interval are computed in base R.

New code:

- `R/m15.R`, the functions of this analysis, sourced by the scripts;
- `scripts/70_m15_validation.R` to `scripts/75_m15_summary.R`;
- tests in `tests/testthat/test-m15.R`.

An independent human quality control (QC) of `R/m15.R` against sections 7 and 8 has not been performed before registration; it is pending before any regulatory use, and the unit tests of section 6.1 are the verification in the meantime. English outputs, including this plan, are checked for Hangul, en dash, em dash and the minus sign U+2212 with the project check (`regulatory/src/reg_helpers.R` check_english), and the plan is also checked to contain only American Standard Code for Information Interchange (ASCII) characters (unit test).

## 15. References

1. ICH. M15 General Principles for Model-Informed Drug Development. Step 4, 2026 (not checked against the current text).
2. Kovalenko P, et al. Exploratory population PK analysis of dupilumab. CPT Pharmacometrics Syst Pharmacol 2016;5(11):617-624.
3. Kovalenko P, et al. Base and covariate population pharmacokinetic analyses of dupilumab using phase 3 data. Clin Pharmacol Drug Dev 2020;9(6):756-767.
4. Li Z, et al. Pharmacokinetics, pharmacodynamics, safety, and tolerability of dupilumab in healthy adult subjects. Clin Pharmacol Drug Dev 2020;9(6):742-755.
5. Clot PA, et al. Pharmacokinetics of dupilumab in Chinese and Japanese healthy subjects. Int Immunopharmacol 2021;99:107985.
6. FDA. Biologics License Application (BLA) 761055 (dupilumab) Clinical Pharmacology Review, 2016 (Table 4.2.c, p.27). Cited by the requester, and throughout this plan, as the "FDA 2017 review" (the review supporting the March 2017 approval); its section 4.6 and Tables 4.1.1.a, 4.2.a, 4.2.b and 4.6.6 are cited from the requester's statements of 2026-10-21.
7. Wilson EB. Probable inference, the law of succession, and statistical inference. J Am Stat Assoc 1927;22:209-212.
8. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics 1988;44:837-845.
9. Tango T. Equivalence test and confidence interval for the difference in proportions for the paired-sample design. Stat Med 1998;17:891-908.
10. Size Korea. 7th Korean anthropometric survey, 2015 (excerpt retrieved by web search on 2026-10-20).
11. ClinicalTrials.gov NCT03389893 (summary statistics recorded in `config/population_atopic.yaml`).

## 16. Amendment 1 (2026-10-21): review response and support for the European Medicines Agency (EMA) scientific advice briefing

### 16.1 Scope, priorities and rules

The request of 2026-10-21 asks for a review response in three priorities: P1, needed before the briefing document deadline of 2026-10-26 and reported first; P2, diagnostics and parameter uncertainty; P3, a conditional re-run. Rules:

- This amendment, its configuration block (`config/prereg_20261020_m15.yaml`, key `amendment_20261021`) and the input file `config/param_uncertainty_m15.yaml` are committed before any new computation. Every change is recorded in the deviations table (sections 13.3 and 13.4), the decision log (D-069; the correction of section 16.14, D-070) and the assumptions register.
- The structure and the estimates of the registered models are not changed. An implementation error (defined in section 16.10) is corrected and disclosed as a deviation; it is not a change of model. Model 1 stays the primary model and Kovalenko 2016 the sensitivity model.
- No derived concept such as an effective dose is used. The body-weight effect is simulated through the model covariate.
- Cells already computed are extracted, not recomputed. Where a new run necessarily recomputes a stored cell, the stored value is used and the recomputation is an identity check: every numeric column is compared with tolerance 0 (1e-12 relative for floating-point values read from CSV files). On any mismatch the run stops, the mismatch is recorded as a deviation, and neither value is used until the cause is found.
- New outputs are written under `results/m15/amend1/`, `results/m15/trials/mechsuppx/` and `results/m15/truth/truth_gmr_ext.csv`. Registered outputs are not overwritten, except for the relabelling of section 16.5, which changes text only. The label changes (the conditionality text in `R/m15.R`; the labels, the outcome wording, conclusion 1 and the applicability label in `scripts/75_m15_summary.R`) are part of the amendment commit, before any amendment run, or of the correction commit 9b61f0b (the out-of-sample labels of figure 1 and the decision-log column of the deviations table; section 13.3 item 15), before the post-stage re-run; one post-stage re-run at a clean commit propagates them. `results/m15/validation/` is not regenerated and keeps its registered wording; the new names apply from the summary stage onward. After the re-run, every CSV file under `results/m15/summary`, `results/m15/samplesize`, `results/m15/weight` and `results/m15/figures` is compared with the file before the re-run, column by column; only the conditionality column and changed label text may differ. The deviations table of the summary takes the decision log entry of each item from its text (D-069 for the items of this amendment).

### 16.2 Model identity

- Model 1, the registered primary model, is Kovalenko 2020 Table 1 Model 1 as recorded in `config/params_k2020_model1.yaml` (decision D-029). The requester confirmed on 2026-10-21 (correction, section 16.14) that this is the "initial base structural model (Model 1)" of section 4.6 of the FDA 2017 clinical pharmacology review of BLA 761055 (reference 6 of this plan; "FDA 2017 review" in this section), developed from four rich-sampling studies (R668-AS-0907, TDU12265, PKM14161 and R668-AD-1117). Model 1 stays primary as registered: Vc 2.48 L, ke 0.0534 1/day, kcp 0.213 1/day, kpc 0.310 1/day, ka 0.256 1/day, MTT 0.105 day with three transit compartments, Vm 1.07 mg/L/day, Km 0.01 mg/L, F 0.643, weight exponent on Vc 0.711 (reference 75 kg). It has no covariate other than weight.
- Rationale stated in the appendix: Model 1 is the first stage of the stepwise models reviewed by FDA; it was developed from rich data including PKM14161, which used the commercial prefilled syringe (product #7); the later models fixed F, ka, Vm and the distribution parameters at the values of this model.
- The values supplied in the first request of 2026-10-21 belong to other models: the standard errors (Vc 2.74 L and the covariate coefficients) are those of Model 4, the phase 3 covariate model, and the "base model" values (Vc 2.76 L, ke 0.0448 1/day, weight exponent 0.919) are those of Model 3, the phase 3 base model with only V2 and ke estimated. They are not used for Model 1. Their configuration files (`config/params_k2020cov.yaml` for Model 4, `config/params_k2020base.yaml` for Model 3; status assumption) are kept for the optional P3 comparison of section 16.10 only.
- The original publications and the FDA review could not be opened from the analysis environment on 2026-10-21; the repository records and the requester's statements are the sources.

### 16.3 P1-1 Sampling

Definitions used in this section and in the appendix:
- lambda_z success: lambda_z computable and AUC_%Extrap_obs below 20%, over all subjects, with its Wilson 95% CI (section 8.4; column success_pct of `results/m15/individual/extrapolation_lambda_z.csv`), labelled "AUC0-inf estimable with extrapolation below 20%". This is the request's "reliable AUC0-inf estimation (extrapolated AUC below 20%)".
- reliable-terminal-phase share: the share of subjects without a reliability flag (adjusted R-squared at least 0.80, extrapolation at most 20%, span ratio at least 2; section 8.3, set AUCinf_A; column reliable_pct). The AUROC uses set AUCinf_A as registered. The AUC0-inf error is against the model AUC to infinity, for set AUCinf_A and for all subjects with a computable lambda_z.
- The words "reliable AUC0-inf" are not used.

1. Extended schedule (Days 64 and 71 added) against the confirmed schedule, 300 mg, global and Korean populations, both models: power at a true GMR of 0.95 (130 randomized per arm, stratified analysis, 75 kg boundary; section 8.1) with the paired Tango 95% CI of the difference; lambda_z success and reliable-terminal-phase share; the error of AUC0-last and AUC0-inf against the dense-sampling truth (median percent error and root mean square error (RMSE) of the log ratio; section 8.4); the AUROC of AUC0-inf for a 10% bioavailability difference (both directions pooled; section 8.3) with its shared cluster bootstrap CI; and the 50% crossings for F, ke, Vmax and ka (section 8.2). All are extracted from stored results except the crossings for ke, Vmax and ka beyond plus or minus 20%, which need item 5. Paired schedule differences of the shares and errors use a paired subject bootstrap (2,000 resamples of the 20,000 subjects of section 8.4, the same indices for all schedules and both models, seed derive_seed(20261020, 'amend_boot')).
2. Reduced schedule (Days 36 to 50 removed): the same error, lambda_z success, reliable-terminal-phase share and AUROC quantities for Kovalenko 2016 (stored), beside Model 1.
3. Coverage, both models side by side: the share of subjects whose Day 57 sample is quantifiable under three definitions (raw: observed value at or above the LLOQ; true: true concentration at the actual time at or above the LLOQ; retained: kept as a quantifiable value in the NCA after the BLQ rules), for Days 57, 64 and 71 where the schedule has them, with lambda_z success and the reliable-terminal-phase share.
4. Rule for sampling conclusions. Coverage statements (the Day 57, 64 and 71 quantifiable shares, lambda_z success and the reliable-terminal-phase share) quote the Kovalenko 2016 value with the Model 1 value in parentheses, because Model 1 over-predicts the healthy-volunteer 300 mg exposure (pooled predicted/observed AUC0-last ratio 1.115 against 1.027 for Kovalenko 2016, section 16.7 item 2) and a higher predicted exposure overstates late quantifiable concentrations (correction of 2026-10-21, section 16.14). For other levels (power, AUROC, crossings, errors) the less favourable model's value is quoted, with the favourable directions power higher, median percent error closer to 0, RMSE lower, AUROC higher and absolute crossing smaller. For a schedule contrast, "X is better than Y" is stated only for class 1 of section 16.10 and is quoted with the smaller of the two models' differences; "X is not worse than Y" is stated only when neither model shows a material difference against X (class 2, or class 3 in favour of X), and is quoted with the larger of the two differences against X. Anti-drug-antibody concentrations are reported as the range over both models (the higher value is the conservative input for a drug-tolerance criterion). Conclusions use the primary scenario cell (global population, 300 mg) and are repeated for the Korean population. This rule deviates from section 2.1 item 1 for sampling conclusions only, as requested; the stored values were seen before it was written (section 16.11).
5. Anti-drug-antibody (ADA) candidate days 15, 29, 57, 64, 71 and 85 (elapsed days 14, 28, 56, 63, 70 and 84) after a single 300 mg dose, global and Korean populations, both models, the 20,000 subjects of section 8.4: distribution of the true individual serum concentration at the nominal time (primary: 1st, 5th, 50th, 95th and 99th percentiles, the 5th to 95th percentiles being the 90% interval, and the share below the LLOQ) and of the observed concentration with residual error at the actual time (5th, 50th and 95th percentiles with BLQ values entered as 0, and the BLQ share). Visit windows are those of the PK samples (assumption B03). Elapsed day 84 is not on the section 8.4 grid: its jitter and residuals use derive_seed(20261020, 'ada', 'jitter') and derive_seed(20261020, 'ada', 'eps'), one draw per subject in subject id order, and its actual time is at least 1 hour after the actual elapsed-day-70 time. No drug-tolerance threshold is applied; the distributions are the input for setting it.
6. New simulation, mode `mechsuppx` of `scripts/73_m15_trials.R`: both models, 300 mg, global and Korean populations, all three schedules, ke, Vmax and ka at multipliers 0.60, 0.70, 1.30, 1.40 and 1.50, and ka also at 1.75 and 2.00, 1,000 studies with the registered trial seeds (common random numbers), 75 kg boundary, 130 randomized per arm, stratified analysis. The ka multipliers 1.75 and 2.00 are added because the stored Model 1 results did not reach the upper ka crossing at 1.50 (60.7% and 65.1% at ka x 1.50; section 16.11). The version order is that of the stored supplementary run (versions 1 to 16), with the two ka versions appended (17 and 18). Identity checks (section 16.1): (i) Model 1, confirmed schedule: rows with vi 1 to 16, keyed by (trial, vi), equal the stored rows of `results/m15/trials/mechsupp/k2020_<population>_300` in every column, the count rows keyed by (vi, n_rand, model) equal the stored counts, and rows 1 to 16 of `versions.csv` are identical; the stored values are used; (ii) both models and all three schedules: the identity rows (vi 1), keyed by (schedule, trial), equal the vi 1 rows of `results/m15/trials/mech/<model>_<population>_300` in every column. The pilot is trials 1 to 50 of Model 1, global population; the full run starts only if both checks pass on it and its run time projects the full run; any mismatch stops the run. The truth extension (mode `truthx` of `scripts/71_m15_populations.R`) computes the induced true GMRs of every version not in `results/m15/truth/truth_gmr.csv` on the section 8.2 truth subjects (derive_seed(20261020, 'truth')): Kovalenko 2016, ke, Vmax and ka at 0.60, 0.70, 1.30, 1.40 and 1.50 and ka at 1.75 and 2.00 (17 versions per population); Model 1, ka at 1.75 and 2.00 (2 versions per population). The recomputed identity metrics must equal `results/m15/truth/truth_reference.csv` within a relative 1e-10. Crossings follow section 8.2 with the registered bootstrap index m15_boot_index(1000, 2000, derive_seed(20261020, 'crossing_boot')), so differences stay paired with the stored results, and with the grid ends at 0.60 and 1.50 (2.00 for ka). They are computed by `scripts/77_m15_amend.R` and written to `results/m15/amend1/` (mech_pass_ext.csv, mech_crossing_ext.csv, mech_crossing_diff_ext.csv); `scripts/75_m15_summary.R` is not changed for them. Check: for every stored cell reached within its stored grid with no unreached bootstrap resample, the recomputed crossing and bounds equal `results/m15/summary/mech_crossing.csv`.

### 16.4 P1-2 Body weight and populations

1. Weight-effect check. The review noted that the 5th and 95th percentile AUC0-inf ratios are 1.207 and 0.841 in both models. Each model's ratios in `results/m15/weight/weight_ratio.csv` are recomputed from that model's own parameters for every typical_ratio row (both doses, both populations, all percentiles; AUC0-inf and Cmax from an independent solver call, AUC0-last from the registered function) with a relative tolerance of 1e-6, and reported with six decimals (`results/m15/amend1/weight_check.csv`). Unit tests: (a) with Vmax multiplied by 1e-9, the typical AUC0-inf ratio at 63.2 kg equals (63.2/75)^(-theta_WT) within 1e-6 for each model; (b) under (a), the two models' ratios differ by more than 5e-4; (c) for each model, the full-model ratio at 63.2 kg changes by more than 1e-3 when theta_WT increases by 0.01; (d) the stored ratios at the 5th and 95th percentile weights equal the independent recomputation within 1e-6. The near-equality is a coincidence, not an expectation: weight acts only through Vc in both models, and each metric's weight elasticity is about theta_WT times its normalized sensitivity to F (Model 1: 0.711 x 1.561 = 1.110; Kovalenko 2016: 0.705 x 1.571 = 1.107). With Vmax near zero the 63.2 kg ratios differ more (1.129428 and 1.128269) than with the full models (1.206272 and 1.206050 at 300 mg). The weight effect is not a defect (correction of 2026-10-21): the two models carry different exponents (0.711 and 0.705). The unit tests are kept; the wording that suspected a defect is removed. The weight ratios of the Model 3 and Model 4 configurations (`results/m15/amend1/weight_ratio_diagnostic.csv`, computed before the correction) are not part of P1 and belong to the optional P3 comparison.
2. Cell detail: both models x global and Korean populations x boundary (75 kg, enrolled median) x analysis model (stratified, unstratified), 300 mg, confirmed schedule, 130 randomized per arm: power at a true GMR of 0.95 with its Wilson 95% interval (Monte Carlo interval: simulation error of 1,000 studies only; it does not include parameter or model uncertainty); the total CV of AUC0-last and Cmax (section 8.4) and the residual CV of the analysis model, sqrt(exp(s^2) - 1), with s^2 = SD^2 - beta^2 x var_between for the stratified analysis and s^2 = SD^2 for the unstratified analysis (inputs of section 8.5 item 2); the share of subjects in the upper stratum and its expected number among 117 evaluable subjects per arm.
3. Variance decomposition per cell: the share of the variance of log AUC0-last explained by log weight (continuous; R-squared of a linear regression) and by the stratum indicator (R-squared), and their ratio (the part of the weight effect absorbed by the stratum term). Empirical values come from the 20,000 subjects of section 8.4 (confirmed schedule, 300 mg); analytic values from the sample-size inputs (beta^2 x variance / SD^2).
4. Global against Korean population, 300 mg, confirmed schedule, both models: 50% crossings for F, ke and Vmax on both sides with bootstrap intervals; the AUROC of AUC0-inf for a 10% F difference; power at 0.95; Korean minus global differences with paired intervals (shared bootstrap indices; Tango interval for power). The Kovalenko 2016 crossings for ke and Vmax beyond plus or minus 20% use section 16.3 item 6.

### 16.5 P1-3 Labels and conclusion 1

- "External validation" is renamed "posterior predictive check (Li 2020 healthy-volunteer studies)" and the "failure rule" is renamed the "misfit rule". Rule R4 is renamed the "rule on out-of-sample units". The rules, thresholds, units and judgements of section 6.2 do not change.
- In-sample status is a reporting column (configuration key amendment_20261021.labels.in_sample); `config/design_clot2021.yaml` model_development_data does not change, and the status agrees with it:
  - Kovalenko 2016: R668-AS-0907 and R668-HV-1108 in-sample (Kovalenko 2016 Table 1); PKM12350, PKM14161, PKM14271 and TDU12265 out-of-sample.
  - Model 1: R668-AS-0907, TDU12265 and PKM14161 in-sample; PKM12350, PKM14271 and R668-HV-1108 out-of-sample (FDA 2017 review, section 4.6).
- Section 2.1 item 3 (the "reproduced" or "not contradicted" wording based on the detection probability at a 25% bias) is suspended until section 16.7 is complete and is then replaced by the template of section 16.7 item 7. This is a deviation from a fixed reading (section 13.3). Neither model is described as "reproduced" or "not contradicted", and the detection probability is reported only as a property of the rule.
- Interim conclusion 1: "Posterior predictive check (Li 2020 healthy-volunteer studies): Model 1 met the misfit rule under the fixed weight assumption only (rule R4: the PKM12350 AUC0-last mean is outside its 99% prediction interval) and did not meet it under the uncertain-mean weight assumption; Kovalenko 2016 did not meet it under either assumption. The cause is under diagnosis (amendment section 16.7)." The numbers and the categories are generated from the result files.
- The conditionality label of the result files changes accordingly (section 16.1 for the propagation and the check).

### 16.6 P1-4 Briefing appendix

- An English modelling and simulation appendix of 5 to 8 pages, tables first, with no em dash or en dash, structured as in ICH M15: question of interest; context of use (design support only; the equivalence conclusion rests on the non-compartmental analysis of the trial); model influence and model risk, with the consequence of a wrong decision; models and sources (with the rationale for Model 1 of section 16.2); posterior predictive check results, stating the in-sample status, the product and assay version of each arm, the between-arm variation check, and the Model 1 over-prediction of the healthy-volunteer 300 mg exposure as found, with and without PKM12350; simulation design; results (operating characteristics, sampling evaluation, sample size, body weight and population); assumptions and limitations. The appendix states that the design conclusions are the same under Kovalenko 2016 where the classification of section 16.10 shows it.
- Body conclusions use only rows of the robustness classification (section 16.10) in class 1 (same direction) or class 2 (no material effect in either model), and sample-size rows classified "both meet" or "neither meets". Rows in class 3 appear in the tables labelled "effect shown in one model only" and support no conclusion; class 4 rows are reported as disagreement in the limitations. Coverage levels follow section 16.3 item 4. For the P1 appendix the classification uses Model 1 as implemented and Kovalenko 2016.
- Excluded from the appendix: 200 mg results; the Model 3 and Model 4 configurations; Korean-population results not contained in the P1 report, including any race or ethnic sensitivity. Korean-population results of sections 16.3 and 16.4 in the P1 report are used under the classes of section 16.10.
- The appendix is generated from the result files (markdown, and Word through pandoc with input format markdown-smart, so that no dash is converted) by `scripts/77_m15_amend.R` mode `appendix`. The markdown source and the text of the Word file (extracted with officer::docx_summary) are checked with the project English check and an ASCII check. Model role labels come from `scripts/75_m15_summary.R` (Model 1 primary, Kovalenko 2016 sensitivity). The posterior predictive check section of the P1 appendix uses section 16.7 items 1 to 3, which need stored replicates only.

### 16.7 P2-1 Posterior predictive diagnostics

1. Full table, from stored replicates without new simulation: every simulated arm x (Cmax, AUC0-last) x both models, fixed and uncertain-mean weight assumptions: observed mean where available, predicted median, 95% and 99% prediction intervals, predicted/observed ratio, percentile rank, whether the observation is outside the 99% interval (verdict), in-sample status, and the product and assay version of the arm (FDA 2017 review Tables 4.1.1.a, 4.2.a and 4.2.b, as supplied by the requester: products #2 (C1P2 vial), #6 (C2P1 vial) and #7 (C2P1 prefilled syringe, the phase 3 and commercial product); assay AV-09095 for PKM12350 and R668-HV-1108, AV-13074 for PKM14161 and PKM14271). Only the product of the PKM14161 test arm (#7) is stated arm by arm in the material available; the other arms read "product not recorded (FDA 2017 review Table 4.1.1.a)" until supplied, and the assay version of TDU12265 and R668-AS-0907 reads "not recorded". The PKM14161 test arm (#7) is flagged as the check closest to the European Union (EU) Dupixent reference product.
2. The seven 300 mg arms that the request describes as about 80 kg (R668-AS-0907; PKM12350, PKM14161 and R668-HV-1108 with two arms each; assumed means 76.55 kg for R668-AS-0907 and 78 kg under B14, realized mean weight of the simulated subjects reported). R668-AS-0907 300 mg has no observed AUC0-last in the repository, so six arms (n 102) carry observed values. Reported per model: the n-weighted observed mean AUC0-last; the predicted median of the n-weighted replicate mean divided by the observed mean, with the 95% and 99% prediction intervals of that ratio (not confidence intervals) and the percentile rank of the observed value; for the six arms and for the four arms without PKM12350; under the fixed, uncertain-mean, 70 kg and 86 kg assumptions. Reading, fixed before the split is computed: Model 1 over-prediction is called systematic across these arms when the observed pooled mean lies below the 2.5th percentile of the predicted pooled mean both with and without PKM12350 (fixed weight assumption); driven by PKM12350 when it does so only with PKM12350; not shown otherwise. The same reading is applied to Kovalenko 2016.
3. Between-arm variation check (correction of 2026-10-21): for the six arms with observed values, y_i = ln(observed mean AUC0-last) with within-arm sampling variance v_i = (SD_i / mean_i)^2 / n_i (delta method). Cochran's Q = sum w_i (y_i - y_w)^2 with w_i = 1 / v_i, 5 degrees of freedom; also I^2 = max(0, (Q - 5) / Q), the DerSimonian-Laird between-arm variance tau^2 and the ratio of the observed variance of y_i to the mean v_i. A simulation companion computes the same Q in each stored replicate of each model (replicate arm means and SDs, fixed weight assumption) and reports the share of replicates with Q at least the observed Q. Between-arm variation is considered present when the Cochran Q p-value is below 0.10 (fixed here, before Q is computed). Beside the decision, the share of stored replicates of each model with a chi-square p-value below 0.10 is reported as a calibration of the criterion (descriptive; added after review, section 13.4 item 3). The two arms of a study are different products, so this variation includes product and study differences.
4. Covariate inputs. Model 1 and Kovalenko 2016 have weight as their only covariate; for healthy volunteers no EASI, albumin, race or BMI input exists or is needed (section 4.2). No healthy-volunteer mapping alternative is computed for them (correction of 2026-10-21). A comparison with Models 3 and 4 is an optional P3 item (section 16.10).
5. False-positive rate when the model is correct, for Model 1 and Kovalenko 2016, from the replicates (each replicate in turn treated as the observed data, limits from all replicates; 5,000 replicates, 1,000 for the 70 and 86 kg assumptions). This replaces the requested 1,000 regenerations of the six studies (same design, n, schedule and NCA rules). Rules:
   - the registered misfit rule (R1, R2, R4 and R5), recomputed from the replicates and checked against `validation_rule.csv`;
   - a per-study rule (descriptive; its outcome on the observed data is already known for PKM12350 and it is not a decision rule): at least one of the following lies outside its two-sided 99% prediction interval (0.5th to 99.5th percentiles of the replicate statistic): the unit AUC0-last means (both arms pooled, section 6.2) of TDU12265 300 mg, TDU12265 600 mg, PKM12350, PKM14161 and R668-HV-1108, and the unit Cmax means of TDU12265 300 mg and 600 mg (7 comparisons). A second version adds the PKM14271 200 mg unit AUC0-last and Cmax means (9 comparisons) and is labelled as limited by the 200 mg presentation. R668-AS-0907 has no observed value and is not included, so the rule covers five of the six requested studies (six with PKM14271);
   - the fixed weight assumption decides; the uncertain-mean assumption is a sensitivity.
   NCA rules: the registered project rules of section 7.4 (intermediate BLQ values missing, BLQ values after Tlast excluded, truncation after two consecutive BLQ values) and the request's rule (every post-dose BLQ value set to 0; AUC0-last to the last positive concentration, linear-up log-down with linear segments ending at 0). The registered rule decides, because the observed data were judged under it; the request's rule is reported beside it, with observed values judged against limits from replicates of the same rule, and with the share of subject profiles whose AUC0-last differs between the rules per arm. The NCA convention of the Li 2020 values is unknown (assumption). Replicates under the request's rule come from `scripts/70_m15_validation.R` mode `alt` with the labels a_model1 and a_k2016 (the registered models, the registered seeds per arm, 5,000 replicates, fixed and uncertain-mean assumptions, output `results/m15/amend1/validation_alt/`); their replicates under the registered rule must equal `results/m15/validation/validation_reps.csv.gz` in every column, or the run stops (about 10 minutes on 4 cores).
6. Between-study exposure variation (sensitivity, only when item 3 finds between-arm variation; correction of 2026-10-21): one log-normal factor per study and replicate, CV 10% (median 1, log-scale SD sqrt(log 1.01)), multiplying the judged replicate's unit AUC0-last and Cmax means, with limits from the replicates without the factor (post-processing of the stored replicates). The factor for study k and replicate r is the r-th draw of rnorm from derive_seed(20261020, 'val', study, 'bsf'); it is shared by both models, all weight assumptions (the 70 and 86 kg assumptions use draws 1 to 1,000), all arms and both TDU12265 dose units of the study, and both AUC0-last and Cmax. CVs are unchanged, so rule R2 is unaffected. Reported with the factor: false-positive rates and detection probabilities (factors 0.75, 0.85, 1.15 and 1.25 applied on top), and the observed data judged against limits computed from the replicates multiplied by the factor (descriptive). If item 3 finds no between-arm variation, the sensitivity is not computed and this is stated.
7. Conclusion 1 is rewritten with this fixed template; only the bracketed values and one category sentence are filled from the result files: "Posterior predictive check (Li 2020 healthy-volunteer studies; [k] of 5 rule units in the development data of [model]): [model] [met / did not meet] the registered misfit rule ([registered outcome category]) through [rules and units]. When [model] is correct, the rule is met by chance with probability [FP0]%[ without and [FP10]% with 10% between-study exposure variation]; the probability of detecting a 25% bias is [D0]%[ and [D10]%] (the smaller of factors 0.75 and 1.25). Between-arm variation among the 300 mg arms: [Cochran Q p-value; present or not detected]. Pooled 300 mg AUC0-last predicted/observed: [r6] with and [r4] without PKM12350 ([systematic / driven by PKM12350 / not shown]). [Category sentence]" Category sentence, judged under the fixed weight assumption and the registered NCA rule, with FP = FP10 when item 6 was computed and FP = FP0 otherwise: (A) rule met and FP at most 10%: "The misfit is unlikely to be explained by chance[ or by 10% between-study variation]."; (B) rule met and FP above 10%: "The misfit is compatible with chance[ and 10% between-study variation] and is not taken as evidence of model error."; (C) rule not met: "The check does not contradict the model." The 10% threshold is fixed here. The same template is written for Kovalenko 2016.

### 16.8 P2-2 Parameter uncertainty

- Input: `config/param_uncertainty_m15.yaml` (each estimate is checked against the model configuration, relative 1e-9, and against its 5th to 95th percentile range). The existing `config/param_uncertainty.yaml` and `scripts/67_day57_param_uncertainty.R` are not changed.
- Kovalenko 2016: the "BLQ data included" column of Kovalenko 2016 Table 2, which is the repository implementation (`config/params_typical.yaml`, `config/params_variability.yaml`), estimates with bootstrap 5th and 95th percentiles as supplied. Standard errors: positive parameters on the log scale, (ln P95 - ln P5)/3.29; the weight exponent and the variances on the linear scale, (P95 - P5)/3.29; independence assumed. Full set: V2, ke, k23, k32, ka, Vm, F, weight exponent, omega2 of V2, ke, ka and Vm, and the proportional residual (log scale). Fixed: Km and the additive residual. Matched-scope set: the eight structural parameters only (V2, ke, k23, k32, ka, Vm, F, weight exponent), for comparison with Model 1.
- Model 1 (correction of 2026-10-21): the Model 4 standard errors are not borrowed, because they leave out the main uncertainties of Model 1 (F, ka, Vm, k23, k32, MTT). Model 1's own bootstrap confidence intervals (Kovalenko 2020 Table 1, FDA 2017 review Table 4.6.6) will be supplied by the requester and then used in a dated addendum, with the same method. Until then, the Kovalenko 2016 bootstrap relative widths are applied to the corresponding Model 1 parameters (Vc for V2, ke, kcp for k23, kpc for k32, ka, Vm, F and the weight exponent), labelled "borrowed uncertainty": positive parameters take the Kovalenko 2016 log-scale SD; the weight exponent takes the Kovalenko 2016 relative standard error, SE = (P95 - P5)/3.29 / 0.705 x 0.711. MTT, Km, the between-subject variances and the residual of Model 1 are not varied. The covariance is not published; independence is assumed and disclosed.
- Scores: z(s) = rnorm(13) from derive_seed(20261020, 'pset', s), assigned in this order: Vc, ke, k12, k21, ka, Vmax, F, theta_WT, omega2 of Vc, ke, ka and Vmax, sigma prop (Kovalenko 2016 V2, k23 and k32 map to Vc, k12 and k21). Each set uses the scores of the parameters it varies, so matched parameters share scores across models and sets. Set 0 is the estimates and equals the model configuration exactly. Redraw: if any variance or residual is not positive or F is at least 1, the whole vector of set s is redrawn from derive_seed(20261020, 'pset', s, 'r1'), then 'r2', up to 'r10'; redraws are counted, and if more than 10 of 1,000 sets are redrawn the result is labelled as from a truncated distribution. Normal draws for the exponent and the variances, and the redraw rule, depart from the log-normal method of section 8.5 item 4 and are recorded as a deviation.
- Method: 1,000 parameter sets per set definition. The same 2,000 subjects (seeds derive_seed(20261020, 'pu', stream), stream in subj, jitter, eps; the section 8.4 pipeline on the union grid, NCA on the confirmed schedule) are used for every set of a model and population, so differences between sets reflect parameters only. Per set and population: the CVs, the correlation and the log-weight slopes are re-estimated and the analytic power (section 8.5; stratified analysis, 75 kg boundary weight moments of `results/m15/samplesize/weight_moments.csv`, Cmax SD from the model; with the model's own CV both Cmax rules agree) is computed at 107 and 130 randomized per arm (96.3 and 117 evaluable), true GMR 0.95, global and Korean populations, 300 mg, confirmed schedule. Each set's power is corrected by a control variate: power_corr(s) = power(s) - power(0; same 2,000 subjects) + power(0; 20,000 subjects stored in `results/m15/individual/ss_inputs_m15.csv`). Assurance is the mean of power_corr over sets 1 to 1,000; the median, the 2.5th and 97.5th percentiles and the share of sets with power_corr of at least 90% are labelled "parameter-uncertainty distribution (not a confidence interval)"; uncorrected values are reported beside them.
- Checks: the per-set pipeline at set 0, run with the section 8.4 subjects and draws (20,000 subjects, derive_seed(20261020, 'ind', stream), jitter and residuals on the union grid, the basic NCA function, which equals the full NCA for AUC0-last and Cmax), reproduces the 300 mg confirmed-schedule rows of `results/m15/individual/ss_inputs_m15.csv` within a relative 1e-10 in every column, or the run stops; the 2,000-subject sets use the 'pu' streams; the subject-sampling error of power(0; 2,000) is the SD over 30 replicate 2,000-subject samples (seeds derive_seed(20261020, 'pu', 'mc', r, stream)).
- The same parameter sets give the parameter-uncertainty interval of section 8.6 item 3 (typical AUC0-inf and Cmax ratios at the 5th and 95th percentile weights against 75 kg, 300 mg, both populations; 2.5th to 97.5th percentiles over sets 1 to 1,000; set 0 checked against `weight_ratio.csv`, relative 1e-6), as registered in section 9.
- Code: `R/pu.R` and `scripts/76_m15_param_assurance.R` (modes check0, mcnoise, sets, summary); outputs under `results/m15/amend1/param_assurance/`.

### 16.9 P2-3 Applicability

- Model 1 is judged against the documented weight ranges of its four development studies: R668-AS-0907 58.0 to 95.1 kg and TDU12265 52.1 to 78.3 kg; the ranges of PKM14161 and R668-AD-1117 are not recorded. The range of 42.3 to 175.4 kg given in the first request belongs to the phase 3 data of Models 3 and 4 and is withdrawn (correction of 2026-10-21).
- Kovalenko 2016: the only documented healthy-volunteer range is 58.0 to 95.1 kg (R668-AS-0907).
- For both models, enrolled weights outside the documented ranges are labelled "outside the documented range" instead of the registered "extrapolation" (sections 6.3 and 6.4; deviation), because the full ranges are unknown. No disease (EASI) dimension is assessed: neither model has a disease covariate (correction of 2026-10-21).

### 16.10 P3 and final outputs

- P3 trigger, fixed before any section 16.7 output is seen. An implementation error of a registered model is declared only when (i) a source excerpt (a publication table or equation, supplied by the requester or opened from the analysis environment, and committed to the repository) shows that a typical value, between-subject variance, residual or covariate coefficient in `config/params_k2020_model1.yaml`, `config/params_typical.yaml` or `config/params_variability.yaml` differs from the published value of that registered model by more than rounding, that an equation in `R/model.R` differs from the published equation of that model, or that a covariate input that the published registered model needs for healthy volunteers is missing or set wrongly; or (ii) a code defect is shown by a unit test or an independent recomputation. Results of section 16.7 (misfit, bias, between-arm variation, false-positive rates) are never by themselves evidence of an implementation error, and Model 3 or Model 4 is never a corrected mapping of Model 1. When P3 is triggered, sections 8.1 to 8.6 are re-run for the affected model with the registered seeds, under a new model identifier with its own parameter file, writing under `results/m15/amend1/rerun/` (never into the registered trial folders, where existing batch files would be skipped); its modes, populations, doses and studies are registered before it runs; the stored results are reported beside the new ones, and the P1 tables and the appendix are updated. If no error is found by the date of the final outputs, P3 is not triggered and this is stated. Replacing Model 1 by Model 3 or Model 4 would be a change of model, not a correction, and is not done without a decision by the requester.
- Optional P3 comparison with Models 3 and 4 (correction of 2026-10-21): run only if the requester asks for it. The configuration files (`config/params_k2020base.yaml` for Model 3, `config/params_k2020cov.yaml` for Model 4) and the mapping labels of the first registration (c_cov, bc_cov_hv80, bc_cov_hv0, d_base and the reference-weight sensitivities; configuration key amendment_20261021.p3_optional) are kept for it; its scope and reporting are registered in a dated addendum before it runs. The Model 3 and Model 4 weight ratios already computed (`results/m15/amend1/weight_ratio_diagnostic.csv`) belong to it.
- Robustness table. Each row is one paired contrast in the primary scenario cell (global population, 300 mg, confirmed schedule, 75 kg boundary, stratified analysis, 130 randomized per arm, true GMR 0.95), except for the factor being contrasted; the number of other cells (Korean population; other boundary, schedule or analysis model) in the same class is shown beside it. Rows:
  - schedule (reduced minus confirmed; extended minus confirmed): power (Tango 95% CI); AUROC of AUC0-inf for a 10% F difference, directions pooled (shared cluster bootstrap 95% CI); F 50% crossing on each side (shared bootstrap 95% CI); lambda_z success and reliable-terminal-phase share; median percent error and RMSE of the log ratio of AUC0-last against the dense truth (paired subject bootstrap of section 16.3 item 1);
  - stratification (stratified minus unstratified; median minus 75 kg boundary): power (Tango);
  - dose (200 minus 300 mg): power (Tango); limited by the 200 mg presentation (section 6.3) and not used in the appendix;
  - Korean minus global: power (Tango), AUROC (bootstrap), F, ke and Vmax crossings (bootstrap);
  - sample size: per model, whether simulated power at 130 randomized per arm meets section 2.1 item 1 for power (lower Wilson 95% bound at least 90%), with the analytic n for 90% power.
  Sign per model: +1 or -1 when the 95% CI excludes 0 and the absolute difference is at least 2 percentage points (power, shares, percent errors, crossings) or at least 0.02 (AUROC, RMSE of the log ratio); otherwise 0. These thresholds mirror section 2.1 item 2; they were set after the mapping agents had seen stored values (section 16.11) and are applied to every row without exception. Classes: (1) same direction: equal non-zero signs; (2) no material effect in either model: both signs 0; (3) material in one model only, with the other model's point estimate stated as same or opposite sign; (4) opposite directions: +1 and -1. The sample-size row is classified as both meet, neither meets, or one meets. A metric without a CI is shown as "no interval, descriptive" and is not classified. Columns: Model 1 as implemented; Model 1 with a corrected mapping, which reads "not applicable (no implementation error found)" unless P3 is triggered (no values are copied into it); Kovalenko 2016; class.
- Assumptions register (plan section 12, exported to `results/m15/summary/m15_assumptions_register.csv`). B08 is replaced: Kovalenko 2016 standard errors from Table 2 bootstrap percentiles (BLQ-included column, divisor 3.29); Model 1 borrowed uncertainty (Kovalenko 2016 relative widths on Vc, ke, kcp, kpc, ka, Vm, F and the weight exponent) until its own bootstrap intervals (FDA 2017 review Table 4.6.6) are supplied; independence. B14 keeps its value; its impact evidence adds the pooled ratio with and without PKM12350 (section 16.7 item 2). B18 records the model identity confirmed on 2026-10-21. New rows: product per arm recorded only for the PKM14161 test arm; ADA visit windows equal to the PK windows; between-study exposure factor log-normal, CV 10%, median 1, used only when between-arm variation is found; NCA convention of the Li 2020 values unknown. Rows on the Model 3 and Model 4 configurations (race factor form, EASI, albumin and BMI factors, reference weight, borrowed variability, BMI reference) apply to the optional P3 comparison only.
- Deviations (section 13.3): one row per item of sections 16.3 to 16.10 that departs from the registered plan, and the items of section 16.11.
- Final outputs: the robustness table, one summary table, at most five conclusion lines, the updated assumptions register (including B08 and B14) and the deviations table.

### 16.11 Known before this amendment

On 2026-10-21, before this amendment was committed, read-only mapping agents computed from stored files or in scratch copies:
- the typical-subject weight ratios of both models to six decimals and the explanation of their near-equality, and prototype ratios for the covariate model (about 1.246 at 63.2 kg with the BMI factor fixed, about 1.289 with BMI varying) and the base model (about 1.284);
- the Day 57 quantifiable shares, the lambda_z success and reliable-terminal-phase shares, the errors against the dense truth and the variance decomposition for the stored 300 mg cells;
- a preview of the Model 1 global-population ADA concentrations;
- the pooled AUC0-last ratio of the six about-80 kg arms with observed values (Model 1 predicted/observed about 1.12, Kovalenko 2016 about 1.03) and a CV-0 reproduction check of the between-study post-processing (the CV 10% result was not computed);
- set-0 power at 107 and 130 per arm, Monte Carlo noise and a one-at-a-time sensitivity of power with invented standard errors;
- extracts of stored power, AUROC and crossing results by schedule and population.
These computations preceded the commit of this amendment, contrary to the request's rule, and this is recorded as a deviation in section 13.3. Elements chosen after these values were seen: the robustness classification of section 16.10 (an earlier four-class version was prototyped on the stored power, AUROC and crossing contrasts); the ka multipliers 1.75 and 2.00; the rule for sampling conclusions of section 16.3 item 4 (registered as the conservative-side rule); the per-study rule of section 16.7 item 5 (the PKM12350 99% outcome was known from section 6.2); the 2,000-subject set size and the control variate (Monte Carlo noise check); and the post-processing form of the between-study factor (CV-0 reproduction). A scratch file with per-arm 99% flags (`map/validation/out_cv0/validation_arms_pi99.csv`) was written by the reproduction check; the mapping report states it was not opened, and it is not used. No other element was chosen with knowledge of a result. The registered runs regenerate every reported value.
The correction of 2026-10-21 (section 16.14) arrived after the stored-replicate arm table and the six-arm pooled ratio of section 16.7 items 1 and 2 had been computed and seen (Model 1 1.115, Kovalenko 2016 1.027; no single arm outside its 99% interval). The between-arm variation criterion (Cochran Q p-value below 0.10), the systematic-bias reading and the coverage basis were written with that knowledge; the four-arm pooled ratio without PKM12350 and Cochran's Q had not been computed.

### 16.12 Computation order, code and the P1 report

- Order: (1) commit of this amendment, its configuration block, `config/param_uncertainty_m15.yaml` and the code changes listed below (commit 323adc9); (2) section 16.4 item 1 checks and unit tests; (3) script 72 mode `amend` (section 16.3 items 3 and 5, section 16.4 item 3); (4) the `mechsuppx` pilot with its identity checks, the full run, `truthx` and the identity check of all cells; (5) commit of the correction of section 16.14 with its code changes; (6) the post-stage re-run at a clean commit with its check (about 15 minutes); (7) the crossing tables, the P1 extractions (section 16.3 items 1, 2 and 4; section 16.4 items 1, 2 and 4), section 16.7 items 1 to 3 from stored replicates, the robustness classification with the P1 results table and conclusions (mode `p1report`), and the appendix draft (mode `appendix`, which reads the classification); (8) the P1 report sent to the requester; (9) the section 16.7 run (about 10 minutes) and post-processing; (10) section 16.8 (about 45 minutes) and section 16.9; (11) P3 and the final outputs of section 16.10. Steps 2 to 4 ran before the correction (section 16.14). Out of this order, mode `ppcstored` (section 16.7 items 1 and 2) first ran at 03:18 UTC on 2026-10-21, during step 4 and before the correction and the post-stage re-run; it was regenerated in step 7 with the corrected labels and the new columns (section 13.3 item 15).
- The P1 report is sent to the requester as soon as step 7 is complete and before step 9 starts: one P1 results table (sampling, coverage, ADA concentrations, weight-effect check, cell detail, global against Korean population, posterior predictive check with products, assays and the between-arm variation check; both models side by side); at most five key conclusion lines written under section 16.3 item 4 and section 16.6; the appendix draft of section 16.6; and the departures table of section 16.13.
- Changes to registered code, committed with this amendment: `R/params.R` (model-to-file map, transit structure detected from the model identifier, borrowed variability), `R/model.R` (transit rate by model identifier, per-subject ke factor `ke_cov`), `R/m15.R` (ke_cov pass-through, `m15_nca_blq0`, the conditionality text), `scripts/70_m15_validation.R` (mode `alt`), `scripts/71_m15_populations.R` (mode `truthx`), `scripts/72_m15_individual.R` (mode `amend`, observation return), `scripts/73_m15_trials.R` (mode `mechsuppx`) and `scripts/75_m15_summary.R` (labels, conclusion 1, decision log entry). New files: `config/params_k2020cov.yaml` (Model 4) and `config/params_k2020base.yaml` (Model 3), for the optional P3 comparison only; `config/param_uncertainty_m15.yaml` (`config/param_uncertainty.yaml` of scripts/67 does not change); `R/pu.R`; `scripts/76_m15_param_assurance.R`; `scripts/77_m15_amend.R`. A unit test pins the registered models: the digests of the parameters and of individual_params for 1,000 global-population subjects (seed 20261021) equal their values at commit 4a948df for Kovalenko 2016 and Model 1. Further unit tests cover the weight effect (section 16.4 item 1), the Model 3 and Model 4 configurations, the request's NCA rule and the standard-error rules.
- `scripts/77_m15_amend.R` implements the extraction, diagnostics, crossing, appendix and robustness steps as specified in this section. It was written after commit 323adc9, while the simulation runs proceeded, and each mode was committed before its first run (section 13.3 item 14).
- Seeds, run logs and the run progress log follow section 11.

### 16.13 Departures from the request

| Requested | Planned | Reason | Effect on conclusions | Confirmation needed |
|---|---|---|---|---|
| 1,000 regenerations of the six studies, post-dose BLQ = 0 | 5,000 stored replicates (registered NCA rule) and 5,000 re-simulated replicates under both NCA rules | More replicates; the observed data were judged under the registered rule | Rates reported under both rules | No |
| Per-study false-positive rate over six studies | Five studies (six with PKM14271); R668-AS-0907 has no observed value | No observation to compare | None | No |
| Seven about-80 kg arms | Six arms with observed AUC0-last (R668-AS-0907 300 mg has none); the seventh is shown as prediction only | No observation to compare | None | No |
| Product and assay per arm (FDA 2017 review Tables 4.1.1.a, 4.2.a, 4.2.b) | Assay by study as supplied; product only for the PKM14161 test arm (#7) | The review cannot be opened from the analysis environment | Product column incomplete | Yes (product per arm) |
| Model 1 parameter uncertainty from Model 1's own bootstrap intervals (FDA 2017 review Table 4.6.6) | Kovalenko 2016 relative widths borrowed until the table is supplied | Table not yet supplied and not accessible | Assurance spread may differ | Yes (Table 4.6.6) |
| Conclusion 1 "under diagnosis" until P2 | Interim wording of section 16.5; fixed template of section 16.7 item 7 | Pre-specification before the rates are seen | None | No |

### 16.14 Correction of 2026-10-21

The requester corrected the first request of 2026-10-21 the same day, after the amendment had been registered (commit 323adc9) and while the P1 runs were in progress. The correction states that the first request described Model 1 wrongly as the phase 3 covariate model and that the repository record is right. Changes to this section (deviation, section 13.3 item 15; decision log D-070):
1. Model identity (section 16.2): Model 1 is Kovalenko 2020 Table 1 Model 1, the initial base structural model of section 4.6 of the FDA 2017 review; the question to the requester in the P1 report is dropped; the supplied standard errors belong to Model 4 and the "base model" values to Model 3.
2. Weight effect (section 16.4 item 1): not a defect; the unit tests stay; the wording that suspected a defect is removed.
3. In-sample status (section 16.5): Model 1 PKM12350, PKM14271 and R668-HV-1108 are out-of-sample (FDA 2017 review, section 4.6); the "202 healthy volunteers" argument is removed. Kovalenko 2016 is unchanged.
4. Product and assay version columns, and the flag on the PKM14161 test arm (section 16.7 item 1).
5. Between-arm variation check and the pooled ratio with and without PKM12350 (section 16.7 items 2 and 3).
6. Healthy-volunteer mappings (b) race and (c) EASI neutral are deleted, because Model 1 and Kovalenko 2016 have no such covariate; the Model 3 and Model 4 comparison becomes an optional P3 item (section 16.10). The mapping labels of the first registration other than a_model1 and a_k2016 are not run.
7. The false-positive rate stays for Model 1 and Kovalenko 2016; the 10% between-study variation sensitivity runs only when item 3 of section 16.7 finds between-arm variation (section 16.7 items 6 and 7).
8. Parameter uncertainty (section 16.8): the Model 4 standard errors are not borrowed; Model 1's own bootstrap intervals will be supplied; until then the Kovalenko 2016 relative widths are applied to the corresponding Model 1 parameters as borrowed uncertainty.
9. Applicability (section 16.9): the 42.3 to 175.4 kg range and the EASI dimension are withdrawn.
10. Coverage statements quote Kovalenko 2016 with Model 1 in parentheses (section 16.3 item 4).
11. The tolerance shortfall of the weight check (20 rows of Model 1 200 mg Cmax, at most 2e-6 against 1e-6, from the registered solver tolerance) is kept as a deviation without changing the criterion, as agreed by the requester (section 13.3 item 15).
12. The appendix states the rationale for Model 1 of section 16.2, the healthy-volunteer 300 mg over-prediction of Model 1, and that the design conclusions are the same under Kovalenko 2016 where the classification shows it (section 16.6).
Results computed before the correction and kept: the script 72 mode `amend` outputs, the `mechsuppx` run and its identity checks, the truth extension, the weight check and the stored-replicate tables of section 16.7 items 1 and 2 (arm table and six-arm pooled ratio). None of them depends on the corrected items, except that the stored-replicate table is regenerated with the corrected in-sample labels and the new columns, and the Model 3 and Model 4 weight ratios move to the optional P3 comparison.

## Appendix A. Abbreviations

ADA, anti-drug antibodies; AI, artificial intelligence; ASCII, American Standard Code for Information Interchange; AUC, area under the concentration-time curve; AUC0-inf, from time zero to infinity; AUC0-last, from time zero to the last quantifiable concentration; AUC0-t, to the last measurable time, as reported in the FDA review; AUC_%Extrap_obs, percentage of AUC0-inf extrapolated beyond Tlast; AUCINF_obs, AUC0-inf from the observed last concentration; AUROC, area under the receiver operating characteristic curve; BLA, Biologics License Application; BLQ, below the limit of quantification; BMI, body mass index; CDF, cumulative distribution function; CI, confidence interval; Clast, last quantifiable concentration; Cmax, maximum concentration; CSV, comma-separated values; CV, coefficient of variation; df, degrees of freedom; EASI, Eczema Area and Severity Index; EMA, European Medicines Agency; F, bioavailability; FP0 and FP10, false-positive rates of the misfit rule without and with 10% between-study exposure variation (D0 and D10, the matching detection probabilities); FDA, US Food and Drug Administration; EU, European Union; GM, geometric mean; GMR, geometric mean ratio; HT, height; I^2, I-squared, the share of the variation between arm means beyond within-arm sampling error; ICH, International Council for Harmonisation of Technical Requirements for Pharmaceuticals for Human Use; IIV, inter-individual variability; ka, absorption rate constant; kcp and kpc, central-to-peripheral and peripheral-to-central rate constants (k12 and k21; k23 and k32 in Kovalenko 2016); ke, linear elimination rate constant; Km, Michaelis-Menten constant; KNHANES, Korea National Health and Nutrition Examination Survey; lambda_z, terminal elimination rate constant; LLOQ, lower limit of quantification; LSODA, Livermore Solver for Ordinary Differential Equations with automatic method switching; MC, Monte Carlo; MCAR, missing completely at random; MM, Michaelis-Menten; MTT, mean transit time; N, number of evaluable subjects; NCA, non-compartmental analysis; OC, operating characteristic; ODE, ordinary differential equation; pAUC, partial AUC; PI, prediction interval; PK, pharmacokinetic(s); QC, quality control; RMSE, root mean square error; RSE, relative standard error; SC, subcutaneous; SD, standard deviation; SE, standard error; Tlast, time of the last quantifiable concentration; UTC, Coordinated Universal Time; V2, central volume of distribution in Kovalenko 2016; Vc, central volume of distribution; Vm, maximum elimination rate as named in the publications (Vmax); Vmax, maximum MM elimination rate; WT, body weight.

## Appendix B. Code key (configuration and result files)

| Name in this plan | Code |
|---|---|
| Global healthy-volunteer population | S1 (uniform shape: S1_uniform) |
| Korean healthy-volunteer population | S2 (uniform shape: S2_uniform) |
| Confirmed schedule | B0 |
| Reduced schedule (study Days 36, 43 and 50 removed) | RM36_50 |
| Extended schedule (study Days 64 and 71 added) | ADD64_71 |
| Unstratified analysis (pooled two-sample t) | M0 |
| Stratified analysis (treatment and weight stratum) | M1 |
| 75 kg boundary; median boundary | fixed75; median |
| Kovalenko 2020 Model 1; Kovalenko 2016 | k2020; k2016 |
| Reliable-terminal-phase set; all subjects with lambda_z | AUCinf_A; AUCinf_B |
| Mean, variability, out-of-sample-unit and Cmax rules | R1, R2, R4, R5 |
| Fixed weight assumption; uncertain-mean assumption; fixed 70 and 86 kg | primary; wunc; w70, w86 |
| Dense-sampling AUC0-last truth | AUClast_dense_<schedule> |

## Appendix C. Internal traceability

The request items map to sections as follows: [0] this plan; [1] sections 4.2 and 6.2; [2] section 7; [3] section 8.1; [4] section 8.2; [5] section 8.3; [6] section 8.4; [7] section 8.5; [8] section 8.6; [9] and [10] section 10. Decision log entries: D-023 (qualification gate), D-029 (Models 3 and 4 judged unsuitable), D-030 (fully external re-judgement), D-031 and D-037 (schedule decision), D-034 (covariates not used for healthy subjects), D-039 (NCA engine), D-055 (stored protocol power figures), D-065 and D-066 (Model 1 as primary in later presentations), D-068 (this plan; D-067 is reserved for a separate pending analysis).
