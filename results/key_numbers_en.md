# Key numbers (dupilumab biosimilar Phase 1 PK simulation)

Generated 2026-09-24. Each line names its source file under results/.

## NCA engine (Phoenix WinNonlin-compatible)

- Own engine vs NonCompart 0.8.4 vs PKNCA 0.12.1: lambda-z points identical in 3,054 of 3,054 profile comparisons; largest relative parameter difference 4.9e-13 (nca_engine/engine_validation_summary.csv).
- AUCinf reliability failure at B0 (any flag or lambda-z not estimable): 19.7% (2016 model) and 16.4% (Model 1); span ratio below 2 flags 8.5% and 9.2% (nca_engine/dropout_reasons_B0.csv).
- Reliability rate, previous vs new engine: 86.4% to 80.3% (2016 model), 91.4% to 83.6% (Model 1) (nca_engine/engine_difference_individual_B0.csv).

## Sampling cliff (2016 model, 60 to 90 kg, 20,000 subjects)

- True LLOQ reached at study Day 38.5 (5th to 95th percentile 27.2 to 54.2); after Day 58 in 2.5% (cliff/cliff_summary.csv).
- Cliff length: 1.38 days (1-day definition), 2.94 days (2-day definition); cliff starts at 0.88 mg/L (median).
- Two or more samples in the cliff at nominal days, fixed schedules, 1-day definition: at most 0.0%; three or more with daily Day 29 to 57 sampling: 0.0% (cliff/cliff_points.csv).

## Sampling density decision (2016 model)

- Final schedule B0. Candidates meeting any pre-specified criterion: none (trials/schedule_decision_base.csv).
- D3 versus B0: AUClast CI width change -0.13% (positive = narrower), reliability gain 2.89 percentage points, NCA extrapolation above 20% ratio 0.590; 200,000 subjects 0.595 (95% CI 0.574 to 0.616).

## Coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 subjects per model)

- True extrapolated share median 0.64% (2016) and 0.63% (Model 1); 95th percentile 3.61% and 3.52%; true coverage below 80%: 0.00 [0.00, 0.02] and 0.00 [0.00, 0.02] (rationale/pillar1_coverage_B0.csv).

