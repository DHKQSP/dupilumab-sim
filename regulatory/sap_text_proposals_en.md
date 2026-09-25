- [1. Primary pharmacokinetic endpoints and primary analysis
  model](#primary-pharmacokinetic-endpoints-and-primary-analysis-model)
- [2. AUC0-inf as a secondary
  endpoint](#auc0-inf-as-a-secondary-endpoint)
- [3. Fallback if FDA requires AUC0-inf as a co-primary
  endpoint](#fallback-if-fda-requires-auc0-inf-as-a-co-primary-endpoint)
- [4. Anti-drug antibodies](#anti-drug-antibodies)
- [5. Concentrations below the lower limit of quantification and
  calculation of PK
  parameters](#concentrations-below-the-lower-limit-of-quantification-and-calculation-of-pk-parameters)
- [6. Sample size](#sample-size)
- [7. Summary of decisions](#summary-of-decisions)

- **Status:** 1.0.1 (draft for sponsor review). Proposed text for review
  by the sponsor’s clinical pharmacology and biostatistics functions;
  not approved.
- **Source commit:** 8d0ebc6b88e7f378b65f3fdc6c74c8ea98708e21. Every
  number below is read from a committed result or configuration file
  when this document is generated, and is listed with its source in
  regulatory/traceability.csv.
- **Use:** each section gives the proposed statistical analysis plan
  (SAP) text, the simulation evidence that supports it, and the
  decisions left to the sponsor. Text in square brackets is to be
  completed by the sponsor.

Abbreviations are defined at first use. The simulations are described in
the Modeling and Simulation Report (regulatory/MS_report).

# 1. Primary pharmacokinetic endpoints and primary analysis model

**Proposed text.**

> The co-primary pharmacokinetic (PK) endpoints are the area under the
> serum concentration-time curve from time zero to the last quantifiable
> concentration (AUC0-last, also written AUC0-t) and the maximum
> observed serum concentration (Cmax). PK similarity is concluded if,
> for both endpoints, the 90% confidence interval (CI) of the ratio of
> geometric means (test/reference; geometric mean ratio, GMR) lies
> within 80.00% to 125.00%.
>
> Option M0: the natural-log-transformed endpoint is compared between
> treatments with a two-sample t-test with pooled variance
> (equivalently, an analysis of variance (ANOVA) with treatment as the
> only factor).
>
> Option M1: the natural-log-transformed endpoint is analysed with an
> ANOVA with treatment and body-weight stratum (60 to 75 kg, above 75 to
> 90 kg; the stratification factor of the randomization) as fixed
> effects. The difference in least-squares means between treatments and
> its 90% CI are back-transformed to the GMR and its 90% CI.
>
> Primary analysis model: \[M0 or M1, to be selected by the sponsor\].
> The other option is a sensitivity analysis. An analysis of covariance
> (ANCOVA) with treatment and log body weight is a further sensitivity
> analysis.

**Evidence.** In 16 boundary scenarios (true AUC0-inf ratio of exactly
0.80 or 1.25, 8 mechanisms and directions in each of two PK models), the
probability of concluding similarity with AUC0-last and Cmax (P2),
i.e. the boundary type I error, was:

- M0: 0.00% to 5.18%; 15 cells conservative (Wilson 95% upper bound
  below 5%), 1 nominal and 0 exceeding (lower bound above 5%). Highest:
  2020 model, peripheral volume up (multiplier 2.28), 5.18% (95% CI 4.88
  to 5.50).
- M1: 0.00% to 5.63%; 15 conservative, 0 nominal and 1 exceeding.
  Highest: 2020 model, peripheral volume up (multiplier 2.28), 5.63%
  (95% CI 5.31 to 5.95). Cells with a point estimate above 5%: 1 under
  M1 and 1 under M0.
- Why M1 differs: the randomization is stratified by body weight, which
  explains part of the between-subject variance, and M0 does not model
  the stratum. With M0 the between-trial standard deviation (SD) of the
  log GMR is smaller than the within-trial standard error (SE) (ratio
  0.961 to 0.993 for the unbiased reference, the true AUC0-inf of each
  subject); with M1 the ratio is 0.991 to 1.020. The unbiased reference
  passes in 4.12% to 5.43% of trials under M0 and 4.53% to 5.83% under
  M1, against 5% expected for an unbiased estimator with a correctly
  sized SE.
- Power of P2 for identical products: 99.6% (95% CI 99.4 to 99.7) (M0)
  and 99.7% (95% CI 99.5 to 99.7) (M1), 2016 model; for a test product
  with bioavailability x0.97: 96.8% (95% CI 96.4 to 97.1) and 97.1% (95%
  CI 96.7 to 97.4).
- Basis: ICH E9 (section 5.7, covariates and subgroups) states that
  factors on which the randomization has been stratified should be
  accounted for in the analysis, and the FDA guidance Adjusting for
  Covariates in Randomized Clinical Trials for Drugs and Biological
  Products (2023) states that a covariate adjustment model should
  generally include the strata variables. Under M0 the stratification is
  not reflected in the analysis.

**Decision left to the sponsor:** the primary analysis model (M0 or M1).

# 2. AUC0-inf as a secondary endpoint

**Proposed text.**

> AUC from time zero to infinity (AUC0-inf) is a secondary endpoint.
> AUC0-inf = AUC0-last + Clast/lambda-z, where Clast is the last
> observed quantifiable concentration and lambda-z the terminal
> elimination rate constant (Phoenix WinNonlin AUCINF_obs). Lambda-z is
> estimated by unweighted log-linear regression using the Best Fit
> method: candidate windows are the last 3, 4, 5 or more positive
> concentrations after Cmax; windows with a positive slope are excluded;
> the window with the largest adjusted coefficient of determination
> (adjusted R-squared) is selected, ties within 0.0001 going to the
> window with more points. Lambda-z is not estimable if fewer than 3
> positive concentrations follow Cmax.
>
> AUC0-inf is considered reliable if the adjusted R-squared of the
> lambda-z regression is at least 0.80 and the extrapolated area
> (Clast/lambda-z) is at most 20% of AUC0-inf. No criterion on the span
> of the lambda-z window is applied. These criteria are entered as the
> Lambda Z Acceptance Criteria of Phoenix WinNonlin (Rules tab), which
> applies no criteria by default and flags, rather than excludes,
> profiles that do not meet them. As a supplementary summary, the number
> of subjects meeting an adjusted R-squared of at least 0.90 is also
> reported for each treatment.
>
> The GMR and 90% CI of AUC0-inf are reported for two analysis sets: (1)
> all subjects with an estimable lambda-z; (2) subjects meeting the
> reliability criteria. For each treatment, the number of subjects with
> lambda-z not estimable and the number not meeting each criterion are
> reported, with the body weight and AUC0-last of the excluded and
> retained subjects. Individual lambda-z windows are listed.
>
> Sensitivity analysis: AUC0-last is substituted for AUC0-inf in
> subjects whose AUC0-inf is not reliable. AUC0-inf does not enter the
> similarity decision.

**Evidence.**

- Lambda-z is estimable in 98.5% to 99.4% of simulated subjects (planned
  schedule, 60 to 90 kg, 20,000 subjects per model). The reliability
  criteria above are met by 86.6% to 91.8%; adding a span-ratio
  criterion (window of at least 2 half-lives) lowers this to 80.7% to
  83.9%.
- The threshold is a convention: 0.80 is the lower end of the
  conventional range, and 0.90 is common in public statistical analysis
  plans (NCT04117607: adjusted R-squared at least 0.90, span at least
  3.0 half-lives, at least 3 points after tmax; NCT04441905: adjusted
  R-squared at least 0.90, at least 3 points; NCT04700163: adjusted
  R-squared above 0.90, otherwise AUC0-inf and the other terminal
  parameters are excluded from summaries and statistical analysis,
  verified from a web-search excerpt). In the trial population (healthy
  adults, 60 to 90 kg, two models) the share of subjects without a
  reliable AUC0-inf is 8.2% to 13.4% with 0.80 and 33.2% to 38.8% with
  0.90 (report Section 5.3). Reporting both counts shows the reader how
  much the secondary analysis depends on the convention.
- A span criterion penalizes denser sampling rather than poor
  estimation: adding late samples (4 alternative schedules, two models)
  increased the share excluded by the span criterion alone by 0.50 to
  3.17 percentage points, because Best Fit then selects a later and
  shorter window.
- Exclusion is not random. Subjects failing the criteria are heavier (by
  1.1 kg in the 2016 model and 0.8 kg in the 2020 model) and have lower
  exposure (ratio of geometric mean true AUC0-inf, excluded to retained:
  0.81 and 0.79); 13.4% and 8.2% of subjects are excluded. Reporting
  both analysis sets, with the characteristics of the excluded subjects,
  lets the reader judge the effect of the exclusion.
- The sensitivity rule (AUC0-last substituted) keeps the boundary type I
  error close to that of the primary endpoints: with AUC0-inf and Cmax
  its highest boundary value is 4.63% (95% CI 4.34 to 4.92) (M0) and
  5.06% (95% CI 4.76 to 5.37) (M1).

# 3. Fallback if FDA requires AUC0-inf as a co-primary endpoint

**Proposed text (to be used only if AUC0-inf is required as a co-primary
endpoint).**

> The co-primary PK endpoints are AUC0-last, AUC0-inf and Cmax. PK
> similarity is concluded if the 90% CI of the GMR lies within 80.00% to
> 125.00% for all three endpoints. For AUC0-inf, all subjects with an
> estimable lambda-z are included; subjects whose AUC0-inf does not meet
> the reliability criteria of Section 2 are flagged and listed but not
> excluded.
>
> Sensitivity analyses of AUC0-inf: (A) excluding the flagged subjects;
> (C) substituting AUC0-last for AUC0-inf in the flagged subjects.

**Evidence.**

- Including all subjects with an estimable lambda-z excludes only 1.5%
  (2016 model) and 0.6% (2020 model) of subjects, and keeps the analysis
  set closest to the randomized set.
- Keeping AUC0-last as a co-primary endpoint bounds the boundary type I
  error of the three-endpoint decision by that of P2: with AUC0-inf
  under this rule it is 0.00% to 3.57% under M0 (16 conservative, 0
  nominal, 0 exceeding; highest 3.57% (95% CI 3.32 to 3.84)) and 0.00%
  to 3.95% under M1 (16, 0, 0; highest 3.95% (95% CI 3.59 to 4.35)).
- Without AUC0-last (AUC0-inf under this rule with Cmax only), 8 of 16
  cells exceed 5% under M0 (highest 15.06% (95% CI 14.37 to 15.77)) and
  10 under M1 (highest 16.03% (95% CI 15.32 to 16.76)): the NCA AUC0-inf
  of this product is biased toward a ratio of 1 at the boundaries.
- Under M1 the three-endpoint decision with sensitivity rule A gives
  0.00% to 3.39% and with rule C 0.00% to 5.16%.
- Power for identical products: 99.5% (95% CI 99.4 to 99.7) (M0) and
  99.6% (95% CI 99.5 to 99.7) (M1), against 99.7% (95% CI 99.5 to 99.7)
  for P2 under M1 (2016 model).

# 4. Anti-drug antibodies

**Proposed text.**

> Serum samples for anti-drug antibodies (ADA) are collected before
> dosing on Day 1 and on Days 15, 29 and 57 (or at early termination).
> ADA status is classified as negative, or positive with
> treatment-emergent or treatment-boosted status and titer; neutralizing
> activity is reported if assessed.
>
> PK parameters are summarized descriptively by treatment and ADA status
> (negative; positive at any post-dose time point; by titer category).
> The primary PK analysis includes all subjects irrespective of ADA
> status. A sensitivity analysis repeats the primary analysis after
> excluding, in both treatments, subjects who are ADA-positive at any
> post-dose time point.

**Evidence.** In a simulated ADA-like subgroup (10% of subjects whose
linear clearance doubles after Day 14; an assumption, not estimated from
data), the median time of the last quantifiable concentration was 28.4
days after dosing, against 34.7 days in the other subjects, and the
median true extrapolated share of AUC0-inf was 0.65% against 0.68%. ADA
mainly shortens the terminal phase, which is where AUC0-inf is least
reliable; AUC0-last remains close to total exposure. Sampling on Day 15
detects early responses and Days 29 and 57 cover the period over which
clearance could change.

# 5. Concentrations below the lower limit of quantification and calculation of PK parameters

**Proposed text.**

> Serum dupilumab is measured with a validated assay with a lower limit
> of quantification (LLOQ) of \[LLOQ of the validated assay, mg/L\]. PK
> parameters are calculated by non-compartmental analysis (NCA) with
> actual sampling times.
>
> Concentrations below the LLOQ (BLQ) are handled as follows: BLQ values
> before the first quantifiable concentration are set to zero (the
> pre-dose value is set to zero); a BLQ value between two quantifiable
> concentrations is set to missing; after two consecutive BLQ values,
> all later values are set to missing; BLQ values after the last
> quantifiable concentration are excluded. \[Handling of a quantifiable
> pre-dose concentration: to be specified by the sponsor.\]
>
> Cmax and the time of Cmax are the observed values (first occurrence).
> AUC0-last is calculated with the linear-up log-down trapezoidal rule
> from time zero to the time of the last quantifiable concentration.
> Lambda-z and AUC0-inf are as in Section 2. No concentrations are
> imputed.

**Evidence.** The simulations used an LLOQ of 0.078 mg/L (the
originator’s assay, Clot 2021), set in a single configuration file so
that all analyses can be regenerated with the sponsor’s assay value.
Across LLOQs of 0.02 to 0.5 mg/L (same subjects and residual draws,
residual error as estimated), in the 2016 model the share of subjects
meeting the reliability criteria of Section 2 ranged from 84.2% (LLOQ
0.02) to 88.0% (LLOQ 0.5), the share with window coverage (true
AUC0-tlast / true AUC0-inf) below 80% from 0.00% to 0.01%, and the
boundary type I error of P2 at the three boundary scenarios examined
(target-mediated elimination capacity in both directions,
bioavailability down) ranged 2.22% to 4.22% under M0 and 2.60% to 4.62%
under M1.

# 6. Sample size

**Proposed text.**

> Assuming a between-subject coefficient of variation (CV) of 43% for
> AUC0-last, the CV and correlation of Cmax from the PK model, and a
> true GMR of 0.95 for both endpoints, \[n\] evaluable subjects per arm
> give \[target\]% power to conclude similarity on both co-primary
> endpoints with the primary analysis model. \[n randomized\] subjects
> per arm are randomized to allow for \[10\]% non-evaluable subjects.

**Evidence** (P2 power; evaluable n per arm for the target power;
randomized n assumes 117 evaluable of 130 randomized):

- With 117 evaluable subjects per arm and a CV of 43%, power is 94.0%
  (M0) and 94.5% (M1); with a CV of 50%, 87.2% and 88.0%.
- For 90% power at CV 43%: 100 (M0) or 97 (M1) evaluable subjects per
  arm; at CV 50%: 129 or 126.
- For 85% power at CV 43%: 85 (M0) or 83 (M1); at CV 50%: 110 or 107.
- CV sources: Li 2020 Table 3 (300 mg arms, SD/mean 35% to 51%); Cohen
  2022 (about 52%, 200 mg, different presentation, device comparison);
  PK model 40% (2016 model) and 43% (2020 model). The full table is in
  results/sample_size/ss_table_n_needed.csv.

**Decisions left to the sponsor:** target power (90% or 85%) and the
sample size.

# 7. Summary of decisions

| Item                                                                                                                      | Status                                           |
|---------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| AUC0-inf as a secondary endpoint (Section 2)                                                                              | confirmed by the sponsor’s clinical pharmacology |
| Reliability criteria: adjusted R-squared at least 0.80 as the definition, count at 0.90 reported alongside (Section 2)    | proposed (sponsor to confirm)                    |
| Disclosure of the boundary scenario whose type I error is above 5%                                                        | confirmed                                        |
| Fallback if AUC0-inf is required as co-primary: all subjects with estimable lambda-z, flagged subjects listed (Section 3) | confirmed                                        |
| CV 43% as the base assumption, 50% as sensitivity                                                                         | confirmed                                        |
| Primary analysis model (M0 or M1)                                                                                         | pending (sponsor)                                |
| LLOQ of the validated assay                                                                                               | pending (sponsor)                                |
| Target power and sample size                                                                                              | pending (sponsor)                                |
