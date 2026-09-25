# v1.0.1

Regulatory package version 1.0.1, generated 2026-09-25 15:57 UTC from commit `8d0ebc6b88e7f378b65f3fdc6c74c8ea98708e21`. Status: draft for sponsor review; not approved. Previous package: document version 0.9 (commit ed9e8ba, not tagged).

## Pre-registration

- Analyses of this version were registered in `config/prereg_20260926.yaml` and committed before their results: analysis model (commit 521645a), numeric criteria for the expectation checks (c5b304b, before any result was read), LLOQ sensitivity and sample size (211e8ed), lambda-z criteria convention and adult atopic dermatitis body weight (325e63b), trial-population analyses of the correction directive (bae3574).

## Results added

- Analysis model (report Section 5.9). Boundary type I error of AUC0-last + Cmax in 16 cells: pooled t-test (M0) 15 conservative, 1 nominal, 0 exceeding (0.00% to 5.18%); ANOVA with the randomization weight stratum (M1) 15 conservative, 0 nominal, 1 exceeding (0.00% to 5.63%). Cells above 5% (point estimate): 2020 model V2_up_080 under M0: 5.18% [4.88, 5.50], nominal, 20,000 trials; 2020 model V2_up_080 under M1: 5.63% [5.31, 5.95], exceeding, 20,000 trials; 2020 model V2_up_080 under M2: 5.75% [5.43, 6.08], exceeding, 20,000 trials.
- Expectations recorded before the results: M1: unbiased reference (AUCinf_true only) close to 5% (consistent); M1: between-trial SD / within-trial SE close to 1 (consistent); P2 rises slightly under M1 (consistent); Model 1 V2 x2.28 cell about 5.6% under M1 (exceeding) (consistent).
- LLOQ (report Section 5.10). The study LLOQ is set in `config/assay.yaml` (single source). Between 0.02 and 0.5 mg/L: window coverage of AUC0-last below 80% in 0.00% to 0.01% of subjects; AUC0-inf reliability, criteria (i), 84.2% (0.02 mg/L) to 88.0% (0.5 mg/L); boundary type I error of AUC0-last + Cmax 2.22% to 4.22% (M0) and 2.60% to 4.62% (M1) in three boundary scenarios.
- Sample size (report Section 5.11). At a true GMR of 0.95, 117 evaluable subjects per arm give P2 power 94.0% (M0) and 94.5% (M1) at CV 43%, 87.2% and 88.0% at CV 50%. Evaluable per arm for 90% power: 100 (M0) and 97 (M1) at CV 43%, 129 and 126 at CV 50%.
- Lambda-z criteria convention (report Sections 3.4, 5.3, 5.4). Phoenix WinNonlin applies no reliability criteria unless the user enters them; adjusted R-squared 0.80 is the lower end of the conventional range. Study population without a reliable AUC0-inf: 8.2% to 13.4% (0.80), 33.2% to 38.8% (0.90), 65.9% to 72.5% (0.90 and span 3). AUC0-inf + Cmax with rule A exceeds 5% in 8, 12, 9 and 12 of 16 boundary scenarios under criteria sets (i) to (iv) (M0).
- Trial population (report Sections 5.2 and 5.3; healthy adults, 60 to 90 kg, weight-stratified randomization, two models). Without a reliable AUC0-inf: 8.2% to 13.4% (criteria set (i)) and 33.2% to 38.8% (set (iii)); retained per arm of 117 under rule A: median 101 to 107 (set (i)) and 72 to 78 (set (iii)). The heavier stratum fails more often (set (i) 1.48 to 2.69 percentage points), and the share failing differs between arms when the products differ (Vmax x1.25: test minus reference 4.27 to 6.64 percentage points, set (i)). The geometric mean of trial AUC0-last GMRs is within 0.012 of the true AUC0-inf ratio in every scenario examined.
- Robustness across body weight (report Appendix I): uniform weight bands, a population with a mean of 100 kg and an adult atopic dermatitis body-weight distribution are robustness checks, not evidence for the proposal (sponsor principle, correction directive).
- Proposed statistical analysis plan text for the PK analyses: `sap_text_proposals_en.md` (criteria convention added to Section 2).
- Anticipated FDA questions: Q1, Q7, Q11, Q12 and Q16 updated; Q17 (analysis model and stratification), Q18 (assay LLOQ), Q19 (is 0.80 a Phoenix WinNonlin criterion) and Q20 (patients outside the study weight range: validity judged in the trial population, wider weight range a robustness check) added.
- Terms: window coverage (true AUC0-tlast / true AUC0-inf) and the observed-to-true ratio (observed AUC0-last or NCA AUC0-inf / true AUC0-inf) are defined and used consistently.
- Report tables and figures of Section 5 renumbered in order (a duplicated table number in version 1.0 corrected).

## Decisions

- Confirmed: AUC0-inf as a secondary endpoint (Best Fit lambda-z; adjusted R-squared at least 0.80 and extrapolation at most 20%, no span criterion; two analysis sets with excluded subjects and their body weight; substitution rule as sensitivity); disclosure of the boundary cell above 5%; fallback if AUC0-inf is required as co-primary (all subjects with an estimable lambda-z, flagged subjects listed; exclusion and substitution as sensitivity analyses); CV 43% as the base assumption and 50% as sensitivity.
- Proposed (sponsor to confirm): adjusted R-squared at least 0.80 as the reliability definition for the secondary AUC0-inf, with the count at 0.90 reported alongside.
- Pending (sponsor): primary analysis model (M0 or M1); LLOQ of the validated assay; target power and sample size.

## Verification

- The regenerated pooled t-test results equal the stored results in 1,553,000 rows (results/oc_models/m0_reverification.csv). LLOQ re-censoring at 0.078 mg/L reproduces the stored individual, cliff and trial results. Automated tests pass.
- Regeneration restarted to add the criteria-set endpoints: 2,040,000 rows of the trials completed before the restart are identical after it (results/criteria/restart_identity_check.csv). AUC0-inf + Cmax under rules A, B and C from the new files equals the version 1.0 values in all 80 cells with equal trial counts (results/criteria/criteria_check_vs_v10.csv).
- Trial-population regeneration: per-arm counts equal the section1 n_R and n_T in all 264,000 compared rows (results/trialpop/tp_identity_check.csv).

## Release assets

- `dupilumab-sim-v1.0.1-source.tar.gz`: repository at the tag (git archive).
- `dupilumab-sim-v1.0.1-regulatory.zip`: the regulatory/ folder.
- `dupilumab-sim-v1.0.1-manifest_sha256.csv`: SHA-256 of every program, configuration, cited result and document.
- `SHA256SUMS.txt`: SHA-256 of the release assets.
