#!/bin/bash
# 검토 의견 2026-09-24(WinNonlin 호환 NCA) §1-2: NCA에 의존하는 모든 산출물을 새 엔진으로 재실행. 두 대기열(개체 수준 1코어, 시험 수준 3코어).
# 이전 엔진 결과의 요약 CSV는 results/nca_engine/legacy_snapshot/ 에 보관(차이표용).
cd "$(dirname "$0")/.."
L=logs/round6_rerun_driver_$(date -u +%Y%m%dT%H%M%S).log; echo "start $(date -u)" > $L
R="nice -n 10 Rscript"
qa() {   # 개체 수준·단계 1 (1코어)
  $R scripts/02_validate_step1.R > logs/r6_step1_k2016.out 2>&1; echo "A step1 k2016 $?" >> $L
  $R scripts/02_validate_step1.R k2020 > logs/r6_step1_k2020.out 2>&1; echo "A step1 k2020 $?" >> $L
  for s in 02b_step1_diagnostics 02c_gate_sampling_error 02d_absorption_diagnostic; do $R scripts/$s.R > logs/r6_$s.out 2>&1; echo "A $s $?" >> $L; done
  for m in k2016 k2020; do $R scripts/02e_arm_weight_sensitivity.R $m > logs/r6_02e_$m.out 2>&1; echo "A 02e $m $?" >> $L; done
  for v in base struct2020 iiv150 resid12 weight_alt ada10 vmax080_both vmax125_both; do $R scripts/10_individual_schedules.R $v > logs/r6_ind_$v.out 2>&1; echo "A individual $v $?" >> $L; done
  $R scripts/13b_crossval_model1.R > logs/r6_crossval_model1.out 2>&1; echo "A crossval_model1 $?" >> $L
  $R scripts/21_curve_shape.R > logs/r6_curve_shape.out 2>&1; echo "A curve_shape $?" >> $L
  $R scripts/24_weight_generalization.R > logs/r6_weight_bands.out 2>&1; echo "A weight_bands $?" >> $L
  $R scripts/24_weight_generalization.R a nojitter 40-60 60-75 75-90 130-150 > logs/r6_weight_bands_nojitter.out 2>&1; echo "A weight_bands_nojitter $?" >> $L
  for v in base struct2020 vmax080_both; do $R scripts/16_criterion_d_200k.R $v > logs/r6_criterion_d_$v.out 2>&1; echo "A criterion_d $v $?" >> $L; done
  echo "A DONE $(date -u)" >> $L
}
qb() {   # 시험 수준 (3코어)
  for v in base struct2020 vmax080_both vmax125_both iiv150 resid12 weight_alt ada10 noresid; do $R scripts/11_trial_schedules.R $v 500 3 > logs/r6_sched_$v.out 2>&1; echo "B sched $v $?" >> $L; done
  for v in base struct2020; do $R scripts/12_trial_products.R 500 3 $v > logs/r6_products_$v.out 2>&1; echo "B products $v $?" >> $L; done
  $R scripts/20_products_5000.R 3 base > logs/r6_products5000.out 2>&1; echo "B products5000 $?" >> $L
  for m in a b d; do $R scripts/25_weight_trials.R $m 3 > logs/r6_weight_trials_$m.out 2>&1; echo "B weight_trials $m $?" >> $L; done
  $R scripts/22_pillar2_curvature.R 3 > logs/r6_pillar2_curvature.out 2>&1; echo "B pillar2_curvature $?" >> $L
  $R scripts/13_crossval.R 500 3 > logs/r6_crossval.out 2>&1; echo "B crossval $?" >> $L
  echo "B DONE $(date -u)" >> $L
}
qa & qb & wait
echo "ALL DONE $(date -u)" >> $L
