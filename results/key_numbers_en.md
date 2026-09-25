# Key numbers (dupilumab biosimilar Phase 1 PK simulation)

Generated 2026-09-25. Each line names its source file under results/.

## NCA engine (Phoenix WinNonlin-compatible)

- Own engine vs NonCompart 0.8.4 vs PKNCA 0.12.1: lambda-z points identical in 3,054 of 3,054 profile comparisons; largest relative parameter difference 4.9e-13 (nca_engine/engine_validation_summary.csv).
- AUCinf reliability failure at B0 (lambda-z not estimable or flagged; flag set (i) = adjusted R-squared at least 0.80 and extrapolation at most 20%, (ii) adds span ratio at least 2): (i) 13.4% [(ii) 19.3%] (2016 model) and (i) 8.2% [(ii) 16.1%] (Model 1); span ratio below 2 flags 8.3% and 9.2% (reliability/dropout_reasons_by_schedule.csv, reliability/reliability_two_flag_sets_summary.csv).
- Reliability rate, previous vs new engine on the same observations (separate B0-only draw): 86.4% to (i) 86.3% [(ii) 80.3%] (2016 model), 91.4% to (i) 91.4% [(ii) 83.6%] (Model 1); the span ratio flag alone accounts for 98.6% and 99.3% of the drop (nca_engine/engine_difference_individual_B0.csv).

## AUCinf reliability: two flag sets (20,000 subjects per model, B0 conditions)

- Flag set (i): lambda-z estimable, adjusted R-squared at least 0.80, extrapolation at most 20%; flag set (ii): (i) plus span ratio at least 2; span ratio at least 2 is a convention used only by some statistical analysis plans (not a Phoenix feature). Cited as (i) [(ii)].
- Reliable at B0: (i) 86.6% [(ii) 80.7%] (2016 model, 95% CI 86.1 to 87.0 [80.2 to 81.3]) and (i) 91.8% [(ii) 83.9%] (Model 1, 95% CI 91.4 to 92.2 [83.4 to 84.4]); loss from the span ratio alone 5.83 and 7.92 percentage points (reliability/reliability_two_flag_sets_summary.csv).
- Two-model range at B0: reliable (i) 86.6 to 91.8% [(ii) 80.7 to 83.9%]; failing (i) 8.2 to 13.4% [(ii) 16.1 to 19.3%].
- D3 versus B0 (paired, same subjects): 2016 model (i) +2.76 (95% CI 2.25 to 3.26) [(ii) +1.69 (95% CI 1.07 to 2.31)] percentage points; Model 1 (i) +1.24 (95% CI 0.80 to 1.67) [(ii) -1.94 (95% CI -2.58 to -1.29)] percentage points (reliability/reliability_paired_vs_B0.csv; (ii) equals c_reliable_gain_pp of trials/schedule_decision_<variant>.csv).
- D1, D2 and D4 versus B0, both models: (i) -1.02 to +0.38 percentage points [(ii) -3.50 to -0.27 percentage points]; in all of D1 to D4 the loss from the span ratio alone rises by 0.50 to 3.17 percentage points; sign reversed by the span ratio in 2016 model D1, 2016 model D2, Model 1 D3.
- Shortest 3-point lambda-z window from nominal days after Day 22: 14 days at B0 and 7 days in D1 to D4, i.e. half-life limits of 7.0 and 3.5 days for span ratio 2 (reliability/reliability_lz_window_by_schedule.csv).
- Subjects reliable at B0 and span-flagged with added sampling (4 model variants x D1 to D4): window shorter in 98.0 to 100.0%, median window 18.5 to 21.5 days at B0 versus 6.7 to 7.3 days, median half-life 6.27 to 7.45 days versus 3.94 to 4.83 days (longer in only 3.2 to 18.0%), added sampling day in the new window in 97.1 to 99.9% (reliability/reliability_span_transition.csv).
- Criterion (c) (gain of at least 5 percentage points) is met under neither flag set for any variant or schedule; the recommendation is the same under both sets (reliability/reliability_paired_vs_B0.csv).
- Numbers that change with the span ratio criterion: reliability and dropout rates, the change in reliability with added sampling (criterion (c)), the rule A analysis set and the pass rates of G2, F3-A and F3-C (rules A and C). Unchanged: criterion (d), AUClast, Cmax, rule B and the schedule recommendation.

## G2 by AUCinf handling rule and flag set (boundary scenarios, true AUC0-inf ratio 0.80 or 1.25)

- Largest boundary pass rate: 2016 model G2-A(i) 19.91% (95% CI 19.14 to 20.70) in Vmax_up_080 (10,000 trials; above 5% in 4 of 8); 2016 model G2-A(ii) 24.58% (95% CI 23.75 to 25.43) in Vmax_down_125 (10,000 trials; above 5% in 6 of 8); 2016 model G2-B 15.06% (95% CI 14.37 to 15.77) in Vmax_up_080 (10,000 trials; above 5% in 4 of 8); 2016 model G2-C(i) 4.13% (95% CI 3.76 to 4.54) in F_down_080 (10,000 trials; above 5% in 0 of 8); 2016 model G2-C(ii) 3.85% (95% CI 3.49 to 4.25) in ke_up_080 (10,000 trials; above 5% in 0 of 8); 2016 model P2 3.97% (95% CI 3.60 to 4.37) in V2_up_080 (10,000 trials; above 5% in 0 of 8); Model 1 G2-A(i) 17.37% (95% CI 16.64 to 18.13) in Vmax_up_080 (10,000 trials; above 5% in 4 of 8); Model 1 G2-A(ii) 27.60% (95% CI 26.73 to 28.48) in Vmax_down_125 (10,000 trials; above 5% in 6 of 8); Model 1 G2-B 11.87% (95% CI 11.25 to 12.52) in Vmax_up_080 (10,000 trials; above 5% in 4 of 8); Model 1 G2-C(i) 4.63% (95% CI 4.34 to 4.92) in V2_up_080 (20,000 trials; above 5% in 0 of 8); Model 1 G2-C(ii) 5.01% (95% CI 4.71 to 5.32) in V2_up_080 (20,000 trials; above 5% in 1 of 8); Model 1 P2 5.18% (95% CI 4.88 to 5.50) in V2_up_080 (20,000 trials; above 5% in 1 of 8) (oc/g2_rules_flags.csv).

## Sampling cliff (2016 model, 60 to 90 kg, 20,000 subjects)

- True LLOQ reached at study Day 38.5 (5th to 95th percentile 27.2 to 54.2); after Day 58 in 2.5% (cliff/cliff_summary.csv).
- Cliff length: 1.38 days (1-day definition), 2.94 days (2-day definition); cliff starts at 0.88 mg/L (median).
- Two or more samples in the cliff at nominal days, fixed schedules, 1-day definition: at most 0.0%; three or more with daily Day 29 to 57 sampling: 0.0% (cliff/cliff_points.csv).

## Operating characteristics (pre-specified design; primary metric: boundary type I error)

- Design pre-registered in commit 779e068 (oc/prereg.csv).
- Largest boundary type I error of P2 (AUClast + Cmax): 5.18% (95% CI 4.88 to 5.50), Model 1, V2 up, true AUC0-inf ratio 0.80; 20,000 trials (adaptive extension; pre-registered 10,000 trials: 5.32% (95% CI 4.90 to 5.78)) (oc/boundary_type1.csv, oc/extension_decision.csv).
- Largest boundary type I error of G2 (AUCinf + Cmax): 27.60% (95% CI 26.73 to 28.48), Model 1, Vmax down, true ratio 1.25; 10,000 trials.
- Boundary scenarios with P2 above 5%: 1 of 16.
- Power for identical products (10,000 trials): 2016 model P2 99.6% (95% CI 99.4 to 99.7); 2016 model F3-A 99.1% (95% CI 98.9 to 99.2); 2016 model G2 99.3% (95% CI 99.2 to 99.5); Model 1 P2 99.0% (95% CI 98.8 to 99.2); Model 1 F3-A 98.3% (95% CI 98.0 to 98.5); Model 1 G2 98.9% (95% CI 98.7 to 99.1) (oc/power.csv).
- Random product space (secondary metric; depends on the assumed virtual product distribution, log-uniform multipliers F 0.80 to 1.25, ka 0.67 to 1.50, ke 0.80 to 1.25, Vmax 0.67 to 1.50, Km 0.20 to 5, V2 0.90 to 1.11), consumer risk (truth outside, 20,000 products): 2016 model P2 0.14% (95% CI 0.08 to 0.25); 2016 model G2 1.08% (95% CI 0.88 to 1.32); Model 1 P2 0.18% (95% CI 0.11 to 0.29); Model 1 G2 1.39% (95% CI 1.17 to 1.67) (oc/random_space_risks.csv).

## Binding constant (Pillar 3 lead sentence)

- Changing the binding constant (Km) from 0.01 to 100 times keeps the true AUC0-inf ratio within 0.999 to 1.059 (both models; range ends and reachable rows with 200,000 common-random-number subjects, screening scan in between with 20,000). Of the 12 pre-specified targets (0.70 to 1.43) the only reachable one is 1.05 (Km about x84 in the 2016 model and x92 in Model 1). Sources: oc/inversion_all.csv, oc/inversion_scan_k2016_Km.csv, oc/inversion_scan_k2020_Km.csv.

## Preliminary product scenarios (arbitrary multipliers, 2016 model; AUCinf reliability under flag set (ii))

- VM150 (Vmax x1.50), 20,000 trials: AUCinf alone (reliable subjects) passes in 6.50% (95% CI 6.16 to 6.84), above the nominal 5% although the true AUC0-inf ratio 0.758 is outside the limits; AUClast alone passes in 0.075% (95% CI 0.045 to 0.124) (trials5000/products5000_props_base.csv; true ratio = mean trial GMR of the individual model AUC0-inf, fallback/consumer_risk.csv).
- KE120 (ke x1.20), 5,000 trials: AUClast pass with AUCinf (reliable subjects) fail in 7.30% (95% CI 6.61 to 8.05); the true AUC0-inf ratio 0.905 is inside the limits, so these are false negatives of AUCinf (fallback/discordance_classification.csv).

## Sampling density decision (2016 model)

- Final schedule B0. Candidates meeting any pre-specified criterion: none (trials/schedule_decision_base.csv).
- D3 versus B0: AUClast CI width change -0.13% (positive = narrower), reliability gain (i) +2.76 [(ii) +1.69] percentage points, NCA extrapolation above 20% ratio 0.599; 200,000 subjects 0.603 (95% CI 0.583 to 0.624).

## Coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 subjects per model)

- True extrapolated share median 0.65% (2016) and 0.64% (Model 1); 95th percentile 3.63% and 3.54%; true coverage below 80%: 0.00 [0.00, 0.02] and 0.00 [0.00, 0.02] (rationale/pillar1_coverage_B0.csv).

