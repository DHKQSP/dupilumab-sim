# Modelling and simulation appendix: design of the phase 1 pharmacokinetic similarity study

Draft for the European Medicines Agency (EMA) scientific advice briefing document. Structured as in the International Council for Harmonisation (ICH) M15 guideline. Generated on 2026-10-21 from the result files at commit 3a6b66a (plan `plans/M15_analysis_plan_20261020.md`, amendment 1, section 16). Every number is read from a result file.

Model 1 is Kovalenko 2020 Table 1 Model 1, the initial base structural model reviewed by the US Food and Drug Administration (FDA) in its 2017 clinical pharmacology review (section 4.6); Kovalenko 2016 is the sensitivity model.

Conditionality: Conditional on a primary model that met the pre-specified misfit rule of the posterior predictive check (plan section 6.2; cause under diagnosis, amendment section 16.7).

Abbreviations: AUC0-inf, area under the concentration-time curve from time zero to infinity; AUC0-last, area to the last quantifiable concentration; AUROC, area under the receiver operating characteristic curve; CI, confidence interval; Cmax, maximum concentration; CV, coefficient of variation; EU, European Union; F, bioavailability; FDA, US Food and Drug Administration; GMR, geometric mean ratio; I2, I-squared, the share of the variation between arm means beyond within-arm sampling error; ka, absorption rate constant; ke, linear elimination rate constant; lambda_z, terminal elimination rate constant; LLOQ, lower limit of quantification; NCA, non-compartmental analysis; PK, pharmacokinetic; RMSE, root mean square error; SD, standard deviation; Vmax, maximum Michaelis-Menten elimination rate.

## 1. Question of interest, context of use and model risk

| Element | Statement |
|---|---|
| Question of interest | Does the planned phase 1 design (single 300 mg subcutaneous dose in healthy adults, parallel groups, weight-stratified randomization, sampling schedule, sample size) detect pharmacokinetic (PK) differences between the proposed biosimilar and the reference product while concluding similarity for similar products? |
| Context of use | Design support only. The similarity conclusion of the study rests on the non-compartmental analysis (NCA) of the observed data: the 90% confidence interval (CI) of the geometric mean ratio (GMR) within 0.80 to 1.25 for AUC0-last and Cmax. Model results are not used in that conclusion. |
| Model influence | Medium: one of several inputs to the design choices, with statistical sample-size methods, precedent and feasibility. |
| Consequence of a wrong decision | Medium: an underpowered study (a similar product fails) or a design less able to detect a true PK difference; the design can be revised before dosing, and the observed-data NCA remains the basis of the conclusion. |
| Model risk | Medium (medium influence, medium consequence). Credibility activities: code verification, posterior predictive check against the original healthy-volunteer studies, applicability assessment, a structural sensitivity model for every main result, and the uncertainty that is not propagated stated explicitly. |

## 2. Models and sources

| Role | Model | Structure | Covariates | Recorded development data |
|---|---|---|---|---|
| Primary | Kovalenko 2020, Table 1 Model 1 | Two compartments, parallel linear and Michaelis-Menten elimination, three transit compartments and first-order absorption | Body weight on the central volume (exponent 0.711, reference 75 kg); no other covariate | R668-AS-0907, TDU12265, PKM14161, R668-AD-1117 (repository record) |
| Sensitivity | Kovalenko 2016, Table 2 (data below the LLOQ included) | Two compartments, parallel linear and Michaelis-Menten elimination, first-order absorption | Body weight on the central volume (exponent 0.705, reference 75 kg) | R668-AS-0907, R668-HV-1108 and patient studies (Kovalenko 2016 Table 1) |

Why Model 1: it is the first stage of the stepwise population PK models reviewed by FDA; it was developed from rich-sampling data that include PKM14161, which used the commercial prefilled syringe (product #7); and the later models fixed bioavailability, the absorption rate constant, the maximum Michaelis-Menten elimination rate and the distribution parameters at the values of this model. For this context of use: Model 1 over-predicts the pooled 300 mg AUC0-last of the healthy-volunteer studies (predicted/observed 1.115, 95% prediction interval 1.026 to 1.209; section 3).

Estimates and structure of both models are used as published; nothing is re-estimated. The publications and the FDA review could not be opened from the analysis environment, so the repository records and the sponsor's records are the source.

## 3. Posterior predictive check (Li 2020 healthy-volunteer studies)

Arm means simulated 5,000 times with each study's design, dose, schedule and LLOQ (0.078 mg/L) under the registered fixed weight assumption (normal 78 kg, SD 10 kg, where weights are not reported). Each cell: predicted median (99% prediction interval); predicted/observed ratio; in-sample status.

| Study arm | Product; assay | Endpoint | Observed mean | Model 1 | Kovalenko 2016 |
|---|---|---|---|---|---|
| PKM12350 reference arm (300 mg, n 13) | product not recorded; assay AV-09095 | AUC0-last | 500.0 | 618.7 (451.6 to 817.9); 1.24; out-of-sample | 569.3 (426.4 to 747.5); 1.14; out-of-sample |
| PKM12350 test arm (300 mg, n 15) | product not recorded; assay AV-09095 | AUC0-last | 487.5 | 618.5 (470.0 to 815.8); 1.27; out-of-sample | 571.0 (437.0 to 746.0); 1.17; out-of-sample |
| PKM14161 reference arm (300 mg, n 19) | product not recorded; assay AV-13074 | AUC0-last | 575.0 | 616.6 (481.1 to 792.2); 1.07; in-sample | 568.5 (446.8 to 719.8); 0.99; out-of-sample |
| PKM14161 test arm (300 mg, n 19); closest to the European Union (EU) reference product | #7, C2P1 prefilled syringe (phase 3 and commercial product); assay AV-13074 | AUC0-last | 587.0 | 617.9 (484.7 to 790.9); 1.05; in-sample | 569.2 (452.7 to 719.9); 0.97; out-of-sample |
| R668-HV-1108 reference arm (300 mg, n 18) | product not recorded; assay AV-09095 | AUC0-last | 530.0 | 618.4 (477.8 to 785.9); 1.17; out-of-sample | 569.7 (447.5 to 710.4); 1.07; in-sample |
| R668-HV-1108 test arm (300 mg, n 18) | product not recorded; assay AV-09095 | AUC0-last | 630.0 | 617.0 (476.4 to 793.5); 0.98; out-of-sample | 568.4 (447.1 to 724.7); 0.90; in-sample |
| TDU12265 300 mg (300 mg, n 6) | product not recorded; assay not recorded | AUC0-last | 700.0 | 734.4 (474.7 to 1,114.1); 1.05; in-sample | 680.2 (451.7 to 1,001.9); 0.97; out-of-sample |
| TDU12265 300 mg (300 mg, n 6) | product not recorded; assay not recorded | Cmax | 38.3 | 37.3 (26.5 to 51.4); 0.97; in-sample | 40.4 (28.3 to 55.6); 1.06; out-of-sample |
| TDU12265 600 mg (600 mg, n 6) | product not recorded; assay not recorded | AUC0-last | 1,780.0 | 1,967.6 (1,361.0 to 2,774.3); 1.11; in-sample | 1,859.7 (1,274.3 to 2,587.2); 1.04; out-of-sample |
| TDU12265 600 mg (600 mg, n 6) | product not recorded; assay not recorded | Cmax | 70.1 | 80.2 (59.2 to 107.3); 1.14; in-sample | 87.5 (64.3 to 116.7); 1.25; out-of-sample |

Products as recorded in the FDA 2017 review Tables 4.1.1.a, 4.2.a and 4.2.b: #2 (C1P2 vial), #6 (C2P1 vial), #7 (C2P1 prefilled syringe, the phase 3 and commercial product). The product of each arm other than the PKM14161 test arm is not yet recorded in the analysis repository.

Pooled 300 mg arms with an observed AUC0-last (n 102, observed n-weighted mean 557 mg x day/L): predicted/observed ratio 1.115 (95% prediction interval 1.026 to 1.209) for Model 1 and 1.027 (0.950 to 1.110) for Kovalenko 2016; the simulated mean weight was 76.8 kg. Without PKM12350 (four arms, n 74): 1.068 (0.968 to 1.175) and 0.983 (0.896 to 1.078). Reading for Model 1: over-prediction driven by PKM12350; for Kovalenko 2016: over-prediction not shown.

Between-arm variation of the six observed arm means: Cochran Q 5.48 on 5 degrees of freedom (p 0.360), I2 9%; variation beyond the within-arm sampling error is not detected. Under each model, the stored replicates give a chi-square p-value below 0.10 in 16.3% (Model 1) and 16.1% (Kovalenko 2016) of replicates (calibration of the criterion). The two arms of each study are different products, so this check covers product and study differences together.

Posterior predictive check (Li 2020 healthy-volunteer studies): Model 1 met the misfit rule under the fixed weight assumption only (rule R4: the PKM12350 AUC0-last mean is outside its 99% prediction interval) and did not meet it under the uncertain-mean weight assumption; Kovalenko 2016 did not meet it under either weight assumption. The cause is under diagnosis (amendment section 16.7).

The false-positive rate of the misfit rule when the model is correct is computed in the next stage and is not reported here.

## 4. Simulation design

| Factor | Setting |
|---|---|
| Populations | Global: normal 76.55 kg (SD 9.275), truncated 58 to 95.1 kg; Korean: normal 65.20 kg (SD 6.550), truncated 52.1 to 78.3 kg; inclusion men 60 to 90 kg, women 50 to 90 kg |
| Dose | 300 mg single subcutaneous dose (200 mg results are not used in this appendix) |
| Sampling schedules | confirmed schedule: elapsed days 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56; reduced schedule: elapsed days 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 56; extended schedule: elapsed days 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56, 63, 70 |
| Sample size | 80 to 160 randomized per arm (base 130); 10% withdrawal completely at random; permuted blocks of two within weight strata |
| Strata and analysis | Two weight strata at 75 kg or at the enrolled median; stratified (treatment + stratum) and unstratified analysis of log AUC0-last and log Cmax |
| Product differences | True GMR grid 0.80 to 1.25 (21 points) applied to test-arm concentrations; mechanism versions with test-arm F, ka, ke or Vmax multiplied by 0.60 to 2.00 |
| Replicates and seeds | 1,000 simulated studies per cell (10,000 at true GMRs 0.80 and 1.25 for Model 1); common random numbers across designs; 20,000 subjects for the sampling evaluation |
| Software | R 4.3.3, rxode2 5.1.7, project NCA engine (linear-up log-down; Phoenix best-fit lambda_z rules) |

## 5. Results

### 5.1 Operating characteristics (300 mg, confirmed schedule, 75 kg boundary, stratified analysis, 130 randomized per arm)

| Population | Model | Power at GMR 0.95 (%) | Power at GMR 1.00 (%) | Consumer risk at 0.80 (%): joint; largest single endpoint | Consumer risk at 1.25 (%): joint; largest single endpoint | Studies at the margins |
|---|---|---|---|---|---|---|
| global | Model 1 | 95.2 (93.7 to 96.4) | 99.2 (98.4 to 99.6) | 2.6 (2.3 to 3.0); 5.3 | 2.3 (2.0 to 2.6); 4.9 | 10,000 |
| global | Kovalenko 2016 | 96.2 (94.8 to 97.2) | 99.7 (99.1 to 99.9) | 1.9 (1.2 to 2.9); 4.7 | 2.4 (1.6 to 3.5); 5.5 | 1,000 |
| Korean | Model 1 | 94.7 (93.1 to 95.9) | 99.2 (98.4 to 99.6) | 2.6 (2.3 to 2.9); 5.2 | 2.6 (2.3 to 2.9); 5.3 | 10,000 |
| Korean | Kovalenko 2016 | 97.1 (95.9 to 98.0) | 99.9 (99.4 to 100.0) | 2.4 (1.6 to 3.5); 5.5 | 2.1 (1.4 to 3.2); 5.2 | 1,000 |

Intervals are Wilson 95% intervals of the simulation (Monte Carlo) error only; they do not include parameter or model uncertainty.

### 5.2 Sampling evaluation (300 mg; Model 1 / Kovalenko 2016)

| Population | Schedule | Power at GMR 0.95 (%) | AUC0-inf estimable, extrapolation below 20% (%) | Reliable terminal phase (%) | AUC0-last median error vs dense truth (%) | AUROC, 10% F difference |
|---|---|---|---|---|---|---|
| global | confirmed schedule | 95.2 / 96.2 | 98.5 / 97.4 | 83.2 / 80.0 | -3.3 / -3.8 | 0.960 / 0.966 |
| global | extended schedule | 95.0 / 96.2 | 98.6 / 97.5 | 83.6 / 80.2 | -3.4 / -3.9 | 0.960 / 0.967 |
| global | reduced schedule | 97.2 / 98.5 | 60.5 / 70.1 | 36.2 / 45.5 | -9.8 / -8.9 | 0.856 / 0.903 |
| Korean | confirmed schedule | 94.7 / 97.1 | 99.1 / 98.5 | 85.3 / 82.2 | -3.1 / -3.5 | 0.952 / 0.963 |
| Korean | extended schedule | 94.5 / 97.1 | 99.2 / 98.5 | 86.0 / 82.5 | -3.1 / -3.6 | 0.955 / 0.965 |
| Korean | reduced schedule | 97.8 / 98.2 | 56.5 / 66.7 | 33.4 / 42.5 | -10.6 / -9.7 | 0.862 / 0.887 |

Day 57 sample on the confirmed schedule, quantifiable and retained in the NCA (Kovalenko 2016, with Model 1 in parentheses, as registered for coverage statements in plan section 16.3 item 4): 3.1% (5.4%) of global and 5.5% (8.8%) of Korean subjects. Rows with a confidence interval are compared between models in section 6.

Predicted true serum concentrations at candidate anti-drug antibody sampling days (single 300 mg dose; range over both models and both populations; 20,000 subjects each):

| Study day | Median (mg/L) | 5th to 95th percentile (mg/L) | Below the LLOQ (%) |
|---|---|---|---|
| Day 15 | 20.87 to 24.87 | 11.70 to 39.44 | 0.0 to 0.0 |
| Day 29 | 6.30 to 10.23 | 4.26e-03 to 22.68 | 3.0 to 9.5 |
| Day 57 | 3.80e-05 to 5.01e-04 | 9.24e-07 to 1.58 | 91.5 to 97.1 |
| Day 64 | 5.32e-06 to 1.98e-04 | 1.11e-07 to 1.95e-03 | 97.5 to 99.3 |
| Day 71 | 7.39e-07 to 7.95e-05 | 1.31e-08 to 7.32e-04 | 99.3 to 99.9 |
| Day 85 | 1.50e-08 to 1.31e-05 | 1.85e-10 to 1.23e-04 | 100.0 to 100.0 |

### 5.3 Sample size (true GMR 0.95, stratified analysis, 10% withdrawal)

| Population | Model | Model CV of AUC0-last (%) | Randomized per arm for 85% power | Randomized per arm for 90% power | Simulated power at 130 per arm (%) | Assurance over CV uncertainty at 130 (%) |
|---|---|---|---|---|---|---|
| global | Model 1 | 42.9 | 92 | 107 | 95.2 (93.7 to 96.4) | 94.5 |
| global | Kovalenko 2016 | 40.0 | 83 | 96 | 96.2 (94.8 to 97.2) | 94.0 |
| Korean | Model 1 | 41.2 | 88 | 104 | 94.7 (93.1 to 95.9) | 94.1 |
| Korean | Kovalenko 2016 | 38.6 | 82 | 95 | 97.1 (95.9 to 98.0) | 93.6 |

### 5.4 Body weight and population

| Quantity | Model 1 | Kovalenko 2016 |
|---|---|---|
| Typical AUC0-inf ratio to 75 kg, global 5th / 95th percentile weight | 1.207 / 0.841 | 1.207 / 0.841 |
| Typical AUC0-inf ratio to 75 kg, Korean 5th / 95th percentile weight | 1.337 / 0.999 | 1.337 / 0.999 |
| Share of log AUC0-last variance explained by log weight (global; Korean) | 0.076; 0.051 | 0.087; 0.058 |
| Part of it absorbed by a 75 kg stratum term (global; Korean) | 73%; 14% | 74%; 14% |
| Share of enrolled subjects above 75 kg (global; Korean) | 54.9%; 5.2% | same populations |

The two models' weight ratios are near-equal by coincidence: weight acts only on the central volume in both models, with exponents 0.711 (Model 1) and 0.705 (Kovalenko 2016), and the implied weight elasticity of AUC0-inf at the global 5th percentile weight is 1.095 and 1.094. Each model was confirmed to run with its own exponent (independent recomputation of every stored ratio).

## 6. Comparison between the two models

Each design contrast is classified by its sign in each model: a sign counts when the 95% interval excludes 0 and the difference is at least 2 percentage points (0.02 for AUROC and for the root mean square error (RMSE) of the log ratio). Class 1: same direction in both models; class 2: no material effect in either model; class 3: material in one model only; class 4: opposite directions.

Rows are the primary scenario cell (global population; Korean minus global for the population contrast); the last column gives the class of the same contrast in the Korean population.

| Contrast | Metric | Model 1 | Kovalenko 2016 | Class | Class, Korean population |
|---|---|---|---|---|---|
| schedule: extended schedule minus confirmed | power at true GMR 0.95 | -0.2 pp, 95% CI -0.7 to 0.2 | 0.0 pp, 95% CI -0.4 to 0.4 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | AUROC (10% F difference) | -0.001, 95% CI -0.003 to 0.002 | 0.000, 95% CI -0.001 to 0.002 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | F lower 50% crossing | 0.0 pp, 95% CI 0.0 to 0.0 | 0.0 pp, 95% CI 0.0 to 0.0 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | F upper 50% crossing | 0.0 pp, 95% CI -0.1 to 0.0 | 0.0 pp, 95% CI 0.0 to 0.0 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | lambda_z success | 0.0 pp, 95% CI 0.0 to 0.1 | 0.0 pp, 95% CI 0.0 to 0.1 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | reliable-terminal-phase share | 0.4 pp, 95% CI 0.3 to 0.5 | 0.2 pp, 95% CI 0.1 to 0.2 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | median AUC0-last error | 0.0 pp, 95% CI 0.0 to 0.0 | 0.0 pp, 95% CI 0.0 to 0.0 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: extended schedule minus confirmed | RMSE of the AUC0-last log ratio | 0.000, 95% CI 0.000 to 0.000 | 0.000, 95% CI 0.000 to 0.000 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: reduced schedule minus confirmed | power at true GMR 0.95 | 2.0 pp, 95% CI 1.3 to 3.1 | 2.3 pp, 95% CI 1.5 to 3.4 | 1, same direction | 3, effect shown in one model only |
| schedule: reduced schedule minus confirmed | AUROC (10% F difference) | -0.105, 95% CI -0.120 to -0.089 | -0.064, 95% CI -0.075 to -0.052 | 1, same direction | 1, same direction |
| schedule: reduced schedule minus confirmed | F lower 50% crossing | -1.0 pp, 95% CI -1.1 to -0.9 | -0.9 pp, 95% CI -1.1 to -0.8 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: reduced schedule minus confirmed | F upper 50% crossing | 1.6 pp, 95% CI 1.4 to 1.9 | 1.8 pp, 95% CI 1.6 to 2.0 | 2, no material effect in either model | 2, no material effect in either model |
| schedule: reduced schedule minus confirmed | lambda_z success | -38.1 pp, 95% CI -38.7 to -37.4 | -27.4 pp, 95% CI -28.0 to -26.7 | 1, same direction | 1, same direction |
| schedule: reduced schedule minus confirmed | reliable-terminal-phase share | -47.0 pp, 95% CI -47.7 to -46.2 | -34.5 pp, 95% CI -35.3 to -33.8 | 1, same direction | 1, same direction |
| schedule: reduced schedule minus confirmed | median AUC0-last error | -6.4 pp, 95% CI -6.5 to -6.3 | -5.0 pp, 95% CI -5.2 to -4.9 | 1, same direction | 1, same direction |
| schedule: reduced schedule minus confirmed | RMSE of the AUC0-last log ratio | 0.072, 95% CI 0.071 to 0.073 | 0.047, 95% CI 0.046 to 0.048 | 1, same direction | 1, same direction |
| stratification: stratified minus unstratified analysis | power at true GMR 0.95 | 0.7 pp, 95% CI 0.3 to 1.4 | 0.2 pp, 95% CI -0.3 to 0.8 | 2, no material effect in either model | 2, no material effect in either model |
| stratification: median minus 75 kg boundary | power at true GMR 0.95 | -0.5 pp, 95% CI -2.3 to 1.3 | 0.0 pp, 95% CI -1.6 to 1.6 | 2, no material effect in either model | 2, no material effect in either model |
| population: Korean minus global | power at true GMR 0.95 | -0.5 pp, 95% CI -2.3 to 1.3 | 0.9 pp, 95% CI -0.6 to 2.4 | 2, no material effect in either model | - |
| population: Korean minus global | AUROC (10% F difference) | -0.009, 95% CI -0.020 to 0.003 | -0.003, 95% CI -0.013 to 0.008 | 2, no material effect in either model | - |
| population: Korean minus global | F lower 50% crossing | -0.2 pp, 95% CI -0.3 to 0.0 | -0.2 pp, 95% CI -0.3 to 0.0 | 2, no material effect in either model | - |
| population: Korean minus global | F upper 50% crossing | 0.5 pp, 95% CI 0.2 to 0.8 | 0.5 pp, 95% CI 0.1 to 0.8 | 2, no material effect in either model | - |
| population: Korean minus global | ke lower 50% crossing | 0.2 pp, 95% CI -0.4 to 0.9 | 0.4 pp, 95% CI -0.1 to 0.9 | 2, no material effect in either model | - |
| population: Korean minus global | ke upper 50% crossing | -0.9 pp, 95% CI -1.7 to -0.1 | -1.0 pp, 95% CI -1.8 to -0.2 | 2, no material effect in either model | - |
| population: Korean minus global | Vmax lower 50% crossing | -2.1 pp, 95% CI -2.8 to -1.4 | -1.7 pp, 95% CI -2.3 to -1.1 | 3, effect shown in one model only | - |
| population: Korean minus global | Vmax upper 50% crossing | 1.2 pp, 95% CI 0.5 to 1.9 | 0.9 pp, 95% CI 0.3 to 1.6 | 2, no material effect in either model | - |

Of the 44 design contrasts above, 42 point the same way in both models or show no material effect in either (classes 1 and 2), so the design conclusions drawn from them are the same under Kovalenko 2016 (Kovalenko 2016 is consistent with the pooled 300 mg AUC0-last of the healthy-volunteer studies; predicted/observed 1.027); 2 show an effect in one model only and 0 disagree in direction.

Conclusions supported by both models (classes 1 and 2 only):

- schedule: extended schedule minus confirmed (global): no material change in either model: power at true GMR 0.95 (largest difference -0.2 pp, 95% CI -0.7 to 0.2), AUROC (10% F difference) (largest difference -0.001, 95% CI -0.003 to 0.002), F lower 50% crossing (largest difference 0.0 pp, 95% CI 0.0 to 0.0), F upper 50% crossing (largest difference 0.0 pp, 95% CI -0.1 to 0.0), lambda_z success (largest difference 0.0 pp, 95% CI 0.0 to 0.1), reliable-terminal-phase share (largest difference 0.4 pp, 95% CI 0.3 to 0.5), median AUC0-last error (largest difference 0.0 pp, 95% CI 0.0 to 0.0), RMSE of the AUC0-last log ratio (largest difference 0.000, 95% CI 0.000 to 0.000).
- schedule: reduced schedule minus confirmed (global): higher in both models: power at true GMR 0.95 (smaller difference 2.0 pp, 95% CI 1.3 to 3.1), RMSE of the AUC0-last log ratio (smaller difference 0.047, 95% CI 0.046 to 0.048); lower in both models: AUROC (10% F difference) (smaller difference -0.064, 95% CI -0.075 to -0.052), lambda_z success (smaller difference -27.4 pp, 95% CI -28.0 to -26.7), reliable-terminal-phase share (smaller difference -34.5 pp, 95% CI -35.3 to -33.8), median AUC0-last error (smaller difference -5.0 pp, 95% CI -5.2 to -4.9); no material change in either model: F lower 50% crossing (largest difference -1.0 pp, 95% CI -1.1 to -0.9), F upper 50% crossing (largest difference 1.8 pp, 95% CI 1.6 to 2.0).
- stratification: stratified minus unstratified analysis (global): no material change in either model: power at true GMR 0.95 (largest difference 0.7 pp, 95% CI 0.3 to 1.4).
- stratification: median minus 75 kg boundary (global): no material change in either model: power at true GMR 0.95 (largest difference -0.5 pp, 95% CI -2.3 to 1.3).
- schedule: extended schedule minus confirmed (Korean): no material change in either model: power at true GMR 0.95 (largest difference -0.2 pp, 95% CI -0.7 to 0.2), AUROC (10% F difference) (largest difference 0.004, 95% CI 0.000 to 0.007), F lower 50% crossing (largest difference 0.0 pp, 95% CI 0.0 to 0.0), F upper 50% crossing (largest difference -0.1 pp, 95% CI -0.1 to 0.0), lambda_z success (largest difference 0.1 pp, 95% CI 0.0 to 0.1), reliable-terminal-phase share (largest difference 0.7 pp, 95% CI 0.6 to 0.9), median AUC0-last error (largest difference 0.0 pp, 95% CI 0.0 to 0.0), RMSE of the AUC0-last log ratio (largest difference 0.000, 95% CI 0.000 to 0.000).
- schedule: reduced schedule minus confirmed (Korean): higher in both models: RMSE of the AUC0-last log ratio (smaller difference 0.057, 95% CI 0.056 to 0.059); lower in both models: AUROC (10% F difference) (smaller difference -0.077, 95% CI -0.089 to -0.064), lambda_z success (smaller difference -31.8 pp, 95% CI -32.4 to -31.2), reliable-terminal-phase share (smaller difference -39.6 pp, 95% CI -40.4 to -38.9), median AUC0-last error (smaller difference -6.2 pp, 95% CI -6.3 to -6.0); no material change in either model: F lower 50% crossing (largest difference -1.3 pp, 95% CI -1.5 to -1.1), F upper 50% crossing (largest difference 1.5 pp, 95% CI 1.3 to 1.7).
- stratification: stratified minus unstratified analysis (Korean): no material change in either model: power at true GMR 0.95 (largest difference 0.1 pp, 95% CI -0.5 to 0.8).
- stratification: median minus 75 kg boundary (Korean): no material change in either model: power at true GMR 0.95 (largest difference 0.6 pp, 95% CI -1.2 to 2.4).
- population: Korean minus global: no material change in either model: power at true GMR 0.95 (largest difference 0.9 pp, 95% CI -0.6 to 2.4), AUROC (10% F difference) (largest difference -0.009, 95% CI -0.020 to 0.003), F lower 50% crossing (largest difference -0.2 pp, 95% CI -0.3 to 0.0), F upper 50% crossing (largest difference 0.5 pp, 95% CI 0.2 to 0.8), ke lower 50% crossing (largest difference 0.4 pp, 95% CI -0.1 to 0.9), ke upper 50% crossing (largest difference -1.0 pp, 95% CI -1.8 to -0.2), Vmax upper 50% crossing (largest difference 1.2 pp, 95% CI 0.5 to 1.9).
- Sample size: simulated power at 130 randomized per arm has a lower 95% bound of at least 90%: global population, both meet; Korean population, both meet.

## 7. Assumptions and limitations

| Item | Statement |
|---|---|
| Healthy-volunteer exposure | Model 1 over-predicts the pooled 300 mg AUC0-last of the healthy-volunteer studies (section 3); under the fixed weight assumption Model 1 met and Kovalenko 2016 did not meet the pre-specified misfit rule. Coverage statements (Day 57 quantifiable share, AUC0-inf estimability) use Kovalenko 2016, with Model 1 alongside (section 16.3 item 4 of the plan). |
| In-sample status | Model 1: R668-AS-0907, TDU12265, PKM14161 in its development data, PKM12350, R668-HV-1108, PKM14271 out-of-sample (FDA 2017 review section 4.6). Kovalenko 2016: R668-AS-0907, R668-HV-1108 in-sample, PKM12350, PKM14161, PKM14271, TDU12265 out-of-sample. |
| Parameter uncertainty | Not propagated in these results; published point estimates are used. Propagation is in progress, with Model 1 uncertainty borrowed from the Kovalenko 2016 bootstrap until Model 1's own bootstrap intervals are entered. |
| Populations | The Korean population differs from the global one by body weight only (Japanese TDU12265 weight range in place of national survey data); no race effect is in the registered models. |
| Weights of literature studies | Studies without reported weights assume a normal 78 kg (SD 10 kg) mean weight; results depend on it. |
| Assay and windows | LLOQ 0.078 mg/L and visit windows are placeholders until the assay and the protocol are final. |
| Immunogenicity | Anti-drug antibodies are not modelled; the concentrations in section 5.2 are inputs for setting the assay drug tolerance. |
| 200 mg presentation | Excluded from this appendix (known AUC0-last under-prediction for that presentation). |
| Disagreement between the models | No design contrast changes in opposite directions in the two models (class 4). |

