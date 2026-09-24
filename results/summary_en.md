# Dupilumab biosimilar Phase 1 pharmacokinetic simulation: summary for regulatory briefing

Generated 2026-09-24 from repository results (branch claude/epic-bardeen-axbreo). Every number below is read from the result files used by the full report.

Abbreviations: area under the concentration-time curve to the last quantifiable concentration (AUClast), to infinity (AUCinf); maximum concentration (Cmax); non-compartmental analysis (NCA); geometric mean ratio (GMR); confidence interval (CI); lower limit of quantification (LLOQ, 0.078 mg/L); inter-individual variability (IIV); target-mediated drug disposition (TMDD); Michaelis-Menten (MM); Monte Carlo (MC); body mass index (BMI).

Models: primary model Kovalenko et al. 2016 (CPT Pharmacometrics Syst Pharmacol 5:617, Table 2, BLQ-included column): two-compartment, first-order absorption, parallel linear and MM elimination, Km fixed at 0.01 mg/L, central volume scaled by (weight/75)^0.705. Sensitivity model Kovalenko et al. 2020 Model 1 (Clin Pharmacol Drug Dev 9:756, Table 1 and Supplementary Table 2): transit absorption (3 compartments, mean transit time 0.105 day), its own IIV and residual error (proportional 15.0%, additive 0.03 mg/L).

Scenario codes (test arm only unless stated; reference arm shared through common random numbers): S00 identical products; F085, F090, F097, F110 bioavailability x0.85, x0.90, x0.97, x1.10; KE110, KE120 linear elimination rate constant (ke) x1.10, x1.20; VM080, VM125, VM150 maximum MM elimination rate (Vmax) x0.80, x1.25, x1.50; KM05, KM2, KM5, KM10 MM constant (Km) x0.5, x2, x5, x10; KA075 absorption rate constant x0.75. Sampling schedules: B0 Syneos baseline; D1 to D4 add two to four samples between Day 32 and Day 53 (D1: Days 39, 46; D2: Days 39, 46, 53; D3: Days 32, 39, 46, 53; D4: Days 40, 47); B- removes Day 50.

Study design simulated: 300 mg single subcutaneous dose (2 mL of 150 mg/mL), parallel groups, 117 evaluable subjects per arm, body weight 60 to 90 kg, Syneos sampling schedule (B0: Days 1, 2, 4, 6, 8, 11, 15, 22, 29, 36, 43, 50, 57). Equivalence: two one-sided tests via the 90% CI of the GMR from a pooled two-sample t on log scale, limits 80.00% to 125.00%.

## 1. Model qualification

Gate scope: study presentation only (300 mg as 2 mL of 150 mg/mL; 600 mg as 2 x 300 mg). Excluded from the gate: PKM14271 200 mg and Cohen 2022 200 mg (different presentation, 1.14 mL of 175 mg/mL, faster absorption with median time to Cmax 3.0 days versus 7.0 days at 300 mg). Decision: option 1, reviewer, sponsor approved (Donghyun Kim), 2026-09-23.

| Dataset | Role | Observed AUClast mean (mg*day/L) | 2016: sim/obs AUClast | 2016: development data | Model 1: sim/obs AUClast | Model 1: development data | 2016: sim/obs Cmax | Model 1: sim/obs Cmax |
|---|---|---|---|---|---|---|---|---|
| 300 mg Chinese (Clot 2021) | gate | 792 | 0.94 | external | 1.02 | external | 1.16 | 1.07 |
| 300 mg Japanese (TDU12265) | gate | 700 | 1.00 | external | 1.08 | internal | 1.08 | 1.00 |
| 300 mg non-Asian pooled (Clot 2021 Table 5, n=40) | gate | 544 | 1.04 | partly internal | 1.13 | partly internal | 1.12 | 1.03 |
| 600 mg Chinese (Clot 2021) | gate | 2110 | 0.98 | external | 1.03 | external | 1.22 | 1.13 |
| 600 mg Japanese (TDU12265) | gate | 1780 | 1.07 | external | 1.13 | internal | 1.27 | 1.17 |
| 200 mg Cohen 2022 autoinjector | external check | 311 | 0.84 | external | 0.92 | external | 0.94 | 0.88 |
| 200 mg Cohen 2022 prefilled syringe | external check | 284 | 0.91 | external | 1.00 | external | 1.01 | 0.94 |
| 200 mg PKM14271 reference | external check | 323 | 0.86 | external | 0.94 | external | 0.94 | 0.87 |
| 200 mg PKM14271 test | external check | 339 | 0.82 | external | 0.90 | external | 0.92 | 0.86 |

Gate result (AUClast mean within 15% of observed, log-scale CV 35% to 51%, cohort median last quantifiable time, weight slope): PASS for both models. Across the study-presentation datasets the simulated/observed AUClast ratio is 0.94 to 1.07 for the 2016 model and 1.02 to 1.13 for Model 1; Cmax is 1.08 to 1.27 and 1.00 to 1.17.

Gate re-judged on datasets fully external to each model's development data (section 7 of the review):
- 2016 model (10 items): FAIL. Failing items: test arm of PKM12350 (1.17).
- Model 1 (8 items): FAIL. Failing items: test arm of PKM12350 (1.27); reference arm of PKM12350 (1.24); reference arm of HV-1108 (1.17).
- The Li 2020 single-arm checks have no reported body weight; they are simulated at an assumed mean of 78 kg. For the 2016 model the PKM12350 ratios fall from 1.23/1.20 at 74 kg to 1.10/1.07 at 82 kg and 1.04/1.02 at 86 kg, so the external-only failure depends on the assumed weight.

200 mg external checks and absorption (diagnostic only, not adopted): the 2016 model does not reproduce the fast absorption of the 200 mg 1.14 mL (175 mg/mL) presentation (observed median time to Cmax 3 days versus 7 days simulated).
Matching absorption to the observed time to Cmax closes about one third to one half of the 200 mg exposure shortfall: the dose-normalized 200:300 mg AUClast ratio is 0.84 observed (four-arm, n-weighted; 0.93 for the PKM14271 test arm alone) versus 0.70 simulated, and 0.75 to 0.77 with the absorption rate constant multiplied by 1.5 to 2 for 200 mg only. The remainder (about 9% to 11% on the four-arm basis) is unexplained: a presentation-specific bioavailability difference or over-estimated TMDD at low exposure are possible; the latter is covered by the population maximum elimination rate (Vmax) x0.8 sensitivity variant. For the study presentation (300 mg, 2 mL of 150 mg/mL) the model reproduces both absorption timing and exposure.

Cross-validation of Model 1 against an independent implementation (first-order absorption approximation of the transit model, 20,000 subjects, acceptance within 3%): 16 of 19 metrics within 3%; the 3 others are small percentages with absolute differences of 0.03 to 0.12 percentage points (share with NCA extrapolation above 20%, median NCA extrapolation, median true extrapolation).

Curve-shape and weight-band results were compared with the reviewer's independent implementation (8,000 subjects per condition): 48 of 51 non-rare metrics within 10%. The 95th percentile of true extrapolation is higher here because sampling-time windows are simulated; without them the same subjects give -1.5% to 2.7% relative to the reference.

## 2. Pillar 1: coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 virtual subjects per model)

| Model | True extrapolated %: median | 95th percentile | Max | True coverage below 80% (%, 95% CI) | NCA extrapolated %: median | NCA extrapolation above 20% (%, 95% CI) | AUCinf reliability met (%, 95% CI) |
|---|---|---|---|---|---|---|---|
| 2016 (primary) | 0.64 | 3.61 | 15.22 | 0.00 [0.00, 0.02] | 2.71 | 1.09 [0.95, 1.24] | 86.63 [86.15, 87.09] |
| Model 1 | 0.63 | 3.52 | 15.68 | 0.00 [0.00, 0.02] | 2.92 | 0.78 [0.67, 0.92] | 91.88 [91.49, 92.25] |

Residual-sensitive metrics are stated as the two-model range: AUCinf reliability criteria (adjusted R-squared at least 0.80 and extrapolation at most 20%) are not met in 8.1% to 13.4% of subjects; NCA extrapolation above 20% occurs in 0.78% to 1.09%. True extrapolation (model integral beyond the last quantifiable time) is essentially the same in both models.

Curve shape below the LLOQ cannot be observed; sensitivity to Km and Vmax (both arms, 20,000 subjects each):

| Variant | True extrapolated %: median | 95th pct | Max | Coverage below 80% (%) | NCA extrap. median (%) | NCA extrap. >20% (%) | AUCinf reliable (%) | Lambda-z not estimable (%) | Median last quantifiable day | AUClast geometric mean | 300 mg ~78 kg AUClast vs observed 544 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Primary model | 0.64 | 3.57 | 16.59 | 0.00 | 2.70 | 1.13 | 86.1 | 1.35 | 34.8 | 546 | 1.04 |
| Km x0.5 | 0.66 | 3.65 | 17.61 | 0.00 | 2.77 | 1.09 | 86.4 | 1.43 | 34.7 | 542 | 1.03 |
| Km x2 | 0.64 | 3.55 | 17.71 | 0.00 | 2.67 | 1.23 | 86.4 | 1.30 | 34.8 | 544 | 1.04 |
| Km x5 | 0.54 | 3.24 | 16.37 | 0.00 | 2.25 | 0.79 | 86.2 | 1.07 | 34.9 | 545 | 1.04 |
| Km x10 | 0.39 | 2.77 | 14.04 | 0.00 | 1.54 | 0.48 | 86.2 | 0.71 | 35.3 | 551 | 1.06 |
| Vmax x0.8 (longer tail, conservative) | 0.53 | 2.79 | 14.21 | 0.00 | 2.24 | 0.44 | 90.5 | 0.39 | 41.3 | 622 | 1.18 |
| Vmax x1.25 | 0.91 | 5.27 | 19.89 | 0.00 | 3.58 | 3.21 | 79.6 | 3.74 | 28.2 | 467 | 0.89 |
| Vmax x0.5 (stress test, fails gate) | 0.63 | 6.41 | 23.08 | 0.04 | 2.33 | 2.20 | 92.9 | 0.01 | 55.0 | 775 | 1.48 |

Kovalenko 2020 reported that excluding below-LLOQ values makes the model predict a less steep TMDD phase with Vm and Km increasing together; this motivates the Km-increase sensitivity.

## 3. Pillar 2: decision concordance between AUClast and AUCinf (B0, 117 per arm)

Primary model, scenarios run with at least 5,000 trials (common random numbers; test-arm multipliers on fixed effects). Pass rates and concordance with Wilson 95% intervals; GMR is the mean over trials.

| Scenario | Trials | GMR AUClast | GMR true AUCinf | GMR NCA AUCinf (reliable) | Pass AUClast (%) | Pass AUCinf reliable (%) | Pass Cmax (%) | Agreement (%) | AUClast pass, AUCinf fail (%) | log GMR correlation |
|---|---|---|---|---|---|---|---|---|---|---|
| F085 |  5000 | 0.762 | 0.767 | 0.782 | 0.3 (0.2 to 0.5) | 1.9 (1.5 to 2.3) | 16.3 (15.3 to 17.4) | 98.1 (97.7 to 98.5) | 0.18 (0.09 to 0.34) | 0.851 |
| F090 |  5000 | 0.840 | 0.843 | 0.853 | 23.4 (22.3 to 24.6) | 33.4 (32.1 to 34.7) | 73.2 (72.0 to 74.4) | 83.8 (82.8 to 84.8) | 3.14 (2.69 to 3.66) | 0.859 |
| F097 |  5000 | 0.952 | 0.953 | 0.956 | 97.0 (96.5 to 97.4) | 96.8 (96.2 to 97.2) | 99.7 (99.5 to 99.8) | 97.0 (96.5 to 97.5) | 1.60 (1.29 to 1.99) | 0.872 |
| F110 |  5000 | 1.169 | 1.165 | 1.155 | 39.3 (38.0 to 40.7) | 49.2 (47.8 to 50.6) | 83.6 (82.5 to 84.6) | 83.5 (82.5 to 84.5) | 3.30 (2.84 to 3.83) | 0.886 |
| KE110 |  5000 | 0.949 | 0.950 | 0.950 | 96.5 (95.9 to 96.9) | 96.1 (95.5 to 96.6) | 100.0 (99.9 to 100.0) | 96.4 (95.8 to 96.9) | 1.98 (1.63 to 2.40) | 0.873 |
| KE120 |  5000 | 0.903 | 0.905 | 0.905 | 76.7 (75.5 to 77.9) | 77.2 (76.0 to 78.4) | 99.7 (99.5 to 99.8) | 87.2 (86.2 to 88.1) | 6.14 (5.51 to 6.84) | 0.872 |
| S00 | 20000 | 1.001 | 1.001 | 1.001 | 99.6 (99.5 to 99.6) | 99.5 (99.4 to 99.6) | 100.0 (99.9 to 100.0) | 99.4 (99.3 to 99.5) | 0.30 (0.24 to 0.39) | 0.874 |
| VM125 |  5000 | 0.858 | 0.866 | 0.889 | 36.9 (35.6 to 38.3) | 62.9 (61.5 to 64.2) | 99.2 (98.9 to 99.4) | 71.5 (70.3 to 72.8) | 1.26 (0.99 to 1.61) | 0.839 |
| VM150 | 20000 | 0.743 | 0.758 | 0.801 | 0.1 (0.0 to 0.1) | 4.7 (4.4 to 5.0) | 89.7 (89.3 to 90.1) | 95.4 (95.1 to 95.7) | 0.00 (0.00 to 0.02) | 0.797 |

Monte Carlo precision: every cited proportion meets the pre-set precision (95% interval half-width at most 1 percentage point below 10% or above 90%, at most 1.5 points otherwise; largest half-width here 1.39 points). Adaptive escalation: VM150 extended to 20,000 trials because an interval included the 90% or 5% threshold at 5,000 and 10,000 trials (identical products run in the same batches for paired comparison); at the cap, VM150 Cmax 89.7% (89.3 to 90.1) still includes the threshold. Trials 1 to 500 reproduce the earlier 500-trial run exactly (maximum difference 0); the 500-trial estimates differ from the independent later trials beyond MC error in 0 of 40 comparisons.

Curvature robustness (population Vmax x0.8 in both arms, 2,000 trials per scenario; mean GMR only):

| Scenario | Trials | GMR AUClast | GMR true AUCinf | GMR NCA AUCinf (reliable) |
|---|---|---|---|---|
| F090 | 2000 | 0.848 | 0.850 | 0.857 |
| KE120 | 2000 | 0.895 | 0.896 | 0.897 |
| KM10 | 2000 | 1.008 | 1.006 | 0.994 |
| S00 | 2000 | 1.000 | 1.000 | 1.001 |
| VM125 | 2000 | 0.876 | 0.881 | 0.897 |

Model 1 (500 trials per scenario; mean GMR only, proportions are reported in the full report): 

| Scenario | GMR AUClast | GMR true AUCinf | GMR NCA AUCinf (reliable) |
|---|---|---|---|
| F085 | 0.764 | 0.768 | 0.782 |
| F090 | 0.841 | 0.844 | 0.853 |
| F097 | 0.952 | 0.953 | 0.955 |
| F110 | 1.167 | 1.163 | 1.155 |
| KA075 | 0.954 | 0.954 | 0.970 |
| KE110 | 0.947 | 0.948 | 0.948 |
| KE120 | 0.900 | 0.901 | 0.902 |
| KM05 | 1.000 | 1.000 | 1.001 |
| KM10 | 1.009 | 1.007 | 0.994 |
| KM2 | 1.002 | 1.001 | 1.000 |
| KM5 | 1.005 | 1.003 | 0.998 |
| S00 | 1.001 | 1.000 | 1.001 |
| VM080 | 1.139 | 1.133 | 1.118 |
| VM125 | 0.860 | 0.867 | 0.888 |
| VM150 | 0.746 | 0.758 | 0.798 |

Largest rate of AUClast pass with AUCinf fail: 6.14 (5.51 to 6.84)% in scenario KE120 (true AUCinf ratio 0.905, inside 80% to 125%, so these are false negatives of AUCinf).

Consumer risk (true AUCinf ratio outside the limits):

| Scenario | True ratio | Trials | Pass AUClast and Cmax (%) | Pass all three (%) | Only AUCinf fails (%) |
|---|---|---|---|---|---|
| F085 | 0.767 |  5000 | 0.32 (0.20 to 0.52) | 0.16 (0.08 to 0.32) | 0.16 (0.05 to 0.27) |
| VM150 | 0.758 | 20000 | 0.08 (0.05 to 0.13) | 0.08 (0.05 to 0.13) | 0.00 (0.00 to 0.00) |

Cost of a three-endpoint fallback (joint pass rates, >= 5,000 trials):

| Scenario | AUClast + Cmax (%) | + AUCinf reliable (%) | + AUCinf all estimable (%) | Loss with reliable set (pp) | Loss with all estimable (pp) |
|---|---|---|---|---|---|
| S00 | 99.5 (99.4 to 99.6) | 99.2 (99.1 to 99.3) | 99.5 (99.4 to 99.6) | 0.30 (0.23 to 0.38) | 0.01 (-0.00 to 0.03) |
| KE110 | 96.5 (95.9 to 96.9) | 94.5 (93.8 to 95.1) | 96.1 (95.6 to 96.6) | 1.98 (1.59 to 2.37) | 0.32 (0.16 to 0.48) |
| F097 | 96.9 (96.4 to 97.4) | 95.4 (94.8 to 95.9) | 96.8 (96.3 to 97.3) | 1.56 (1.22 to 1.90) | 0.14 (0.04 to 0.24) |

AUCinf handling rules if AUCinf were mandated (A: exclude subjects failing reliability, B: include all estimable, C: substitute AUClast when reliability fails):

| Rule | Analysis n per arm | Pass, identical products (%) | GMR bias vs truth, Vmax x1.25 (%) | GMR bias vs truth, F x0.90 (%) | Agreement with AUClast, Vmax x1.25 (%) |
|---|---|---|---|---|---|
| A | 101.1 | 99.5 (99.4 to 99.6) | 2.65 | 1.17 | 71.5 |
| B | 115.4 | 99.8 (99.7 to 99.8) | 1.87 | 0.78 | 74.4 |
| C | 117.0 | 99.6 (99.5 to 99.7) | -0.47 | -0.14 | 95.4 |
| Reference: AUClast | 117.0 | 99.6 (99.5 to 99.6) | -0.90 | -0.38 |  |

Sample size cross-check (proposed log-scale coefficient of variation 43%, 117 subjects per arm):
- Log-scale SD (CV) at B0, 60 to 90 kg, 20,000 subjects: AUClast 0.388 (40.3%) for the 2016 model and 0.410 (42.7%) for Model 1; Cmax 0.334 (34.3%) and 0.330 (33.9%).
- The AUClast 90% CI reported by Cohen 2022 (0.96 to 1.28, n 62 and 63) implies a log-scale SD of about 0.49 (CV about 52%), above both models (40.3% and 42.7%). The proposed 43% lies between the models and the Cohen estimate.
- Empirical power with 117 per arm (2016 model, 20,000 trials): identical products 99.5 (99.4 to 99.6)% for AUClast and Cmax jointly and 99.2 (99.1 to 99.3)% with AUCinf (reliable set) added; test bioavailability x0.97 (true AUC ratio about 0.95) 96.9 (96.4 to 97.4)% and 95.4 (94.8 to 95.9)%.

## 4. Pillar 3: invisibility of binding-constant differences

Test-arm Km multiplied by 0.5 to 10 (0.005 to 0.1 mg/L) changes the mean GMR of every endpoint by at most 1.18% (2016 model) and 0.85% (Model 1) relative to identical products (paired within the same 500 trials). Km is an MM approximation constant and is not identical to binding affinity.

## 5. Body weight generalization (300 mg, B0, 20,000 subjects per uniform weight band)

Covariate variants (c) and (d) use the adult coefficients of Kovalenko 2020 Model 4 (elimination rate constant ke proportional to (BMI/26)^0.368, central volume exponent 0.817); the BMI reference of 26 is a placeholder based on phase 3 mean BMI 25.4 to 27.3 (Kamal 2022). Height is simulated as normal (mean 170 cm, SD 9, truncated 150 to 195 cm). Development-data weight ranges are not reported in the source publications; bands above 130 kg are flagged as possible extrapolation.

| Model | Band (kg) | Median BMI | True extrap. median (%) | 95th pct | Coverage <80% (%) | AUCinf reliable (%) | Lambda-z not estimable (%) | AUClast geometric mean | Note |
|---|---|---|---|---|---|---|---|---|---|
| (a) 2016 | 40-60 | 17.2 | 0.50 | 2.69 | 0.005 | 90.8 | 0.23 | 860 |  |
| (a) 2016 | 60-75 | 23.2 | 0.60 | 3.23 | 0.000 | 88.1 | 0.74 | 615 |  |
| (a) 2016 | 75-90 | 28.5 | 0.72 | 4.01 | 0.000 | 84.5 | 1.88 | 488 |  |
| (a) 2016 | 90-110 | 34.5 | 0.88 | 5.18 | 0.000 | 80.5 | 3.48 | 386 |  |
| (a) 2016 | 110-130 | 41.4 | 1.10 | 6.61 | 0.000 | 75.6 | 5.76 | 307 |  |
| (a) 2016 | 130-150 | 48.3 | 1.37 | 8.42 | 0.010 | 71.2 | 8.43 | 251 | outside confirmed development range |
| (b) Model 1 | 40-60 | 17.3 | 0.50 | 2.88 | 0.000 | 94.1 | 0.14 | 916 |  |
| (b) Model 1 | 60-75 | 23.2 | 0.57 | 3.26 | 0.005 | 92.5 | 0.41 | 658 |  |
| (b) Model 1 | 75-90 | 28.5 | 0.70 | 3.91 | 0.000 | 90.7 | 0.73 | 522 |  |
| (b) Model 1 | 90-110 | 34.4 | 0.82 | 4.99 | 0.005 | 88.1 | 1.76 | 414 |  |
| (b) Model 1 | 110-130 | 41.4 | 1.05 | 6.14 | 0.005 | 85.2 | 2.90 | 328 |  |
| (b) Model 1 | 130-150 | 48.4 | 1.21 | 7.33 | 0.005 | 82.2 | 4.08 | 272 | outside confirmed development range |
| (c) 2016 + ke~BMI | 40-60 | 17.2 | 0.53 | 3.20 | 0.010 | 90.5 | 0.21 | 938 |  |
| (c) 2016 + ke~BMI | 60-75 | 23.3 | 0.59 | 3.17 | 0.000 | 88.0 | 0.73 | 626 |  |
| (c) 2016 + ke~BMI | 75-90 | 28.5 | 0.75 | 4.13 | 0.000 | 84.7 | 1.83 | 480 |  |
| (c) 2016 + ke~BMI | 90-110 | 34.5 | 0.93 | 5.42 | 0.000 | 80.5 | 3.57 | 365 |  |
| (c) 2016 + ke~BMI | 110-130 | 41.4 | 1.23 | 7.22 | 0.005 | 75.2 | 6.53 | 282 |  |
| (c) 2016 + ke~BMI | 130-150 | 48.4 | 1.52 | 9.15 | 0.000 | 70.1 | 9.28 | 227 | outside confirmed development range |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 40-60 | 17.2 | 0.53 | 3.57 | 0.020 | 90.4 | 0.14 | 1007 |  |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 60-75 | 23.3 | 0.60 | 3.17 | 0.000 | 88.3 | 0.90 | 640 |  |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 75-90 | 28.5 | 0.74 | 4.17 | 0.000 | 84.4 | 1.90 | 472 |  |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 90-110 | 34.5 | 1.00 | 5.62 | 0.005 | 79.4 | 4.30 | 346 |  |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 110-130 | 41.4 | 1.33 | 8.06 | 0.015 | 73.0 | 7.75 | 256 |  |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | 130-150 | 48.3 | 1.74 | 9.85 | 0.010 | 66.0 | 11.43 | 199 | outside confirmed development range |

Trial level in a population with many obese subjects (weight normal mean 100 kg, SD 20, truncated 60 to 150 kg; 2,000 trials per scenario; mean GMR and dropout weights only):

| Model | Scenario | GMR AUClast | GMR true AUCinf | GMR NCA AUCinf reliable | Subjects failing AUCinf reliability (%) | Mean weight retained (kg) | Mean weight failing (kg) |
|---|---|---|---|---|---|---|---|
| (a) 2016 | F090 | 0.831 | 0.835 | 0.851 | 21.6 | 99.6 | 104.8 |
| (a) 2016 | KE120 | 0.908 | 0.911 | 0.912 | 19.9 | 99.7 | 104.9 |
| (a) 2016 | KM10 | 1.013 | 1.009 | 0.979 | 18.7 | 99.9 | 104.3 |
| (a) 2016 | S00 | 1.000 | 1.000 | 1.001 | 19.6 | 99.7 | 104.9 |
| (a) 2016 | VM125 | 0.839 | 0.849 | 0.885 | 24.2 | 99.5 | 104.7 |
| (b) Model 1 | F090 | 0.832 | 0.836 | 0.850 | 13.4 | 100.1 | 105.0 |
| (b) Model 1 | KE120 | 0.906 | 0.908 | 0.908 | 12.1 | 100.2 | 105.1 |
| (b) Model 1 | KM10 | 1.013 | 1.009 | 0.984 | 12.5 | 100.3 | 104.0 |
| (b) Model 1 | S00 | 1.001 | 1.001 | 1.001 | 12.1 | 100.2 | 105.1 |
| (b) Model 1 | VM125 | 0.840 | 0.849 | 0.882 | 15.3 | 100.0 | 105.0 |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | F090 | 0.832 | 0.836 | 0.855 | 23.3 | 99.3 | 105.5 |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | KE120 | 0.906 | 0.908 | 0.910 | 21.5 | 99.4 | 105.7 |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | KM10 | 1.014 | 1.010 | 0.973 | 19.8 | 99.7 | 104.9 |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | S00 | 1.001 | 1.000 | 1.001 | 21.2 | 99.5 | 105.6 |
| (d) 2016 + ke~BMI + Vc~weight 0.817 | VM125 | 0.840 | 0.850 | 0.890 | 26.0 | 99.1 | 105.4 |

## 6. Consistency with published coverage values

| Source | Dose (mg) | Published AUClast/AUCinf (mean ratio) | 2016: NCA / true | Model 1: NCA / true |
|---|---|---|---|---|
| Clot 2021 Table 3 | 200 | 94.6% | 94.7% / 98.6% | 95.1% / 98.7% |
| Clot 2021 Table 3 | 300 | 98.1% | 96.7% / 99.2% | 96.7% / 99.1% |
| Clot 2021 Table 3 | 600 | 98.1% | 92.3% / 96.8% | 92.8% / 96.9% |
| FDA BLA 761055 Clin Pharm Review Table 4.2.c (PKM12350 control arm) | 300 | 96.0% | 96.3% / 99.0% | 96.4% / 99.1% |

Clot 2021 conditions: sampling Days 2 to 57, LLOQ 0.078 mg/L, cohort weights from Clot Table 1 (62.2, 59.9, 58.3 kg); the 200 mg row is an external check because of the different presentation. PKM12350 control arm: AUC0-t 500 versus AUC0-inf 521 (FDA BLA 761055 Clinical Pharmacology Review, Table 4.2.c); the test arm (488 versus 554) may reflect a different number of subjects with AUC0-inf and cannot be verified from public documents.
Published statements on the terminal phase: Kovalenko 2020 describes a terminal slope tending to minus infinity near the LLOQ and a nearly vertical TMDD phase, no meaningful terminal half-life, and instantaneous half-life falling to zero; Kovalenko 2021 notes that 0.09 mg/L removes 90% of circulating target at Km 0.01 mg/L; Kovalenko 2016 used the M3 method because few quantifiable low concentrations describe the steep phase; Li 2020 reports steeper elimination at lower concentrations and a more than dose-proportional AUClast; Cohen 2022 did not compute half-life or AUCinf because of terminal non-linearity; Clot 2021 describes multi-exponential decline with faster target-mediated elimination at low concentrations.
Interpretation: published coverage values are based on NCA AUCinf, whose extrapolation is inflated, so they are a lower bound of the true coverage. The curve shape below the LLOQ cannot be observed and is addressed by the Km and Vmax sensitivity analyses.
At 600 mg the simulated true coverage (96.8% and 96.9%) is below the published NCA value (98.1%). Because the published value is a lower bound, both models predict a larger tail beyond Day 56 than observed at this dose, which understates AUClast coverage (conservative direction). At the study dose of 300 mg the NCA mean ratio differs by 1.4 percentage points, within 3 percentage points.

## 7. Sampling density between Day 36 and Day 50: conclusion

Final schedule: B0 (Syneos baseline). Governing model: Kovalenko 2016. Rule application: pre-specified rule applied to the primary model; sensitivity variants assess robustness only. Decided by reviewer, sponsor approved (Donghyun Kim) on 2026-09-24.

Pre-specified rule, versus B0: recommend added sampling if at least one holds: (a) mean width of the AUClast 90% CI decreases by at least 2%; (b) AUClast pass rate increases by at least 2 percentage points when ke is multiplied by 1.10; (c) the AUCinf reliability rate rises by at least 5 percentage points; (d) the share of subjects with NCA extrapolation above 20% falls to half or less.

| Variant | Rule result | D3: AUClast CI width change (%, positive = wider) | D3: reliability gain (pp) | D3: extrapolation >20% ratio, 20,000 subjects | D3: same ratio, 200,000 subjects | Final recommendation |
|---|---|---|---|---|---|---|
| 2016 (primary) | no schedule meets any criterion | 0.06 (0.03 to 0.09) | 2.89 | 0.59 (0.52 to 0.66) | 0.59 (0.57 to 0.62), 0.53 fewer subjects per arm; not met (interval above 0.5) | B0 |
| Model 1 | no schedule meets any criterion | -0.11 (-0.13 to -0.08) | 1.30 | 0.67 (0.60 to 0.75) | 0.71 (0.68 to 0.73), 0.27 fewer subjects per arm; not met (interval above 0.5) | B0 |
| Vmax x0.8 (both arms) | D3 by criterion d only | -0.40 (-0.43 to -0.36) | 2.69 | 0.48 (0.36 to 0.61) | 0.59 (0.55 to 0.62), 0.23 fewer subjects per arm; not met (interval above 0.5) | B0 (conclusion unchanged) |
| Vmax x1.25 (both arms) | no schedule meets any criterion | 0.49 (0.47 to 0.51) | 2.77 | 0.77 (0.74 to 0.81) | not re-evaluated | B0 |

AUClast CI width: D1, D2 and D4 widen the mean AUClast 90% CI by 0.36% to 0.41% (paired intervals exclude zero; D3 0.06%). Recomputed with true concentrations and no residual error, the widening remains (0.43% to 0.52%), so it reflects heterogeneity of the added tail area as the last quantifiable time is extended, not measurement error at the added low points.

Removing Day 50 (B-): AUClast CI width -0.66% (negative = narrower), AUCinf reliability change -0.07 (-0.27 to 0.12) percentage points, share with NCA extrapolation above 20% multiplied by 1.29 (1.22 to 1.39). Not recommended, to keep a terminal sample for AUC0-inf as a secondary endpoint and for a fallback analysis.

Rationale: at the pre-specified 20,000-subject level, criterion (d) alone was met in Vmax x0.8 (both arms); criteria (a), (b) and (c) were not met in any variant; re-evaluated with 200,000 subjects (4,000 bootstrap resamples), the D3 ratio for criterion (d) is above 0.5 with its whole interval in every re-evaluated variant; the absolute reduction is below one subject per arm; neither the AUClast CI width nor power improves; the cost is 1,040 additional visits. Removing the Day 50 sample (B-) is not recommended because it preserves a terminal point for AUC0-inf as a secondary endpoint and for a fallback analysis. Lesson recorded: a relative-reduction criterion for a rare event needs an absolute floor (for example at least one subject per arm); not applied retroactively.

## 8. Limitations

- The primary model does not reproduce the faster absorption of the 200 mg 1.14 mL (175 mg/mL) presentation; about one third to one half of the 200 mg shortfall is explained by absorption, the remainder is unexplained. Clot 2021 Chinese 200 mg data (same 175 mg/mL prefilled syringe) showed median time to Cmax 7.0 days (range 3.0 to 10.0), but with n=8 and no sample between Day 4 and Day 8 the resolution is low.
- Cmax: the 2016 model predicts about 20% above observed and Model 1 about 5% above; AUClast: the 2016 model is close to observed and Model 1 about 12% above.
- Km is fixed at 0.01 mg/L in both models; with no uncertainty or IIV on Km, terminal-phase variability may be underestimated. Km and Vmax sensitivity analyses address this.
- The AUCinf reliability rate and the share with NCA extrapolation above 20% depend on the residual error model and are reported as two-model ranges.
- Placeholders not yet confirmed: Day 1 post-dose sampling time (0.25 day), weight distribution and stratification split, sampling windows, BMI reference of 26, body weights of the Li 2020 single-arm studies.
- Development-data body weight ranges are not reported; results above 130 kg are extrapolations.

