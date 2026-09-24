#!/bin/bash
# 검토 의견 통합본(2026-09-24) 계산 드라이버. 두 대기열:
#  Q1(3코어, 시험 수준): 20_products_5000 → 25_weight_trials(a,b,d) → 22_pillar2_curvature → 12_trial_products(struct2020, 500) → 11 struct2020 → 11 noresid
#  Q2(1코어, 개체 수준): 02 step1(k2016, k2020) → 10 individual struct2020 → 13b crossval → 21 curve shape → 24 weight bands → 16 criterion d 200k ×3
cd "$(dirname "$0")/.."
L=logs/round5_driver_$(date -u +%Y%m%dT%H%M%S).log; echo "start $(date -u)" > $L
q1() {
  while pgrep -f "scripts/20_products_5000.R" > /dev/null; do sleep 30; done      # 이미 실행 중이면 기다린다(재시작 시 완료 묶음은 건너뜀)
  [ -f results/trials5000/products5000_props_base.csv ] || { Rscript scripts/20_products_5000.R 3 base >> logs/products5000_base.out 2>&1; echo "Q1 products5000 $?" >> $L; }
  echo "Q1 products5000 ready" >> $L
  for m in a b d; do Rscript scripts/25_weight_trials.R $m 3 > logs/weight_trials_$m.out 2>&1; echo "Q1 weight_trials $m $?" >> $L; done
  Rscript scripts/22_pillar2_curvature.R 3 > logs/pillar2_curvature.out 2>&1; echo "Q1 pillar2_curvature $?" >> $L
  Rscript scripts/12_trial_products.R 500 3 struct2020 > logs/trial_products_struct2020.out 2>&1; echo "Q1 products struct2020 $?" >> $L
  while ! grep -q "Q2 individual struct2020" $L; do sleep 20; done
  Rscript scripts/11_trial_schedules.R struct2020 500 3 > logs/trial_sched_struct2020.out 2>&1; echo "Q1 sched struct2020 $?" >> $L
  Rscript scripts/11_trial_schedules.R noresid 500 3 > logs/trial_sched_noresid.out 2>&1; echo "Q1 sched noresid $?" >> $L
  echo "Q1 DONE $(date -u)" >> $L
}
q2() {
  while ! grep -q EXIT logs/step1_k2020_round5.out 2>/dev/null; do sleep 20; done; echo "Q2 step1 done" >> $L
  Rscript scripts/10_individual_schedules.R struct2020 > logs/ind_struct2020.out 2>&1; echo "Q2 individual struct2020 $?" >> $L
  Rscript scripts/13b_crossval_model1.R > logs/crossval_model1.out 2>&1; echo "Q2 crossval_model1 $?" >> $L
  Rscript scripts/21_curve_shape.R > logs/curve_shape.out 2>&1; echo "Q2 curve_shape $?" >> $L
  Rscript scripts/24_weight_generalization.R > logs/weight_bands.out 2>&1; echo "Q2 weight_bands $?" >> $L
  for v in base struct2020 vmax080_both; do Rscript scripts/16_criterion_d_200k.R $v > logs/criterion_d_200k_$v.out 2>&1; echo "Q2 criterion_d_200k $v $?" >> $L; done
  echo "Q2 DONE $(date -u)" >> $L
}
q1 & q2 & wait
echo "ALL DONE $(date -u)" >> $L
