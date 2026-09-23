#!/bin/bash
# 단계 2–4 전체 계산(두 후보 모델). 결과는 results/individual, results/trials. 보고는 gate 판정 후(SPEC §1).
cd "$(dirname "$0")/.."
L=logs/step2_driver_$(date -u +%Y%m%dT%H%M%S).log
echo "start $(date -u)" > $L
# 1) 개인 수준: 변형을 4개 프로세스로 병렬
( Rscript scripts/10_individual_schedules.R base       > logs/ind_base.out 2>&1; echo "ind base $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R struct2020 > logs/ind_struct2020.out 2>&1; echo "ind struct2020 $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R iiv150 resid12 > logs/ind_iiv_resid.out 2>&1; echo "ind iiv150 resid12 $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R weight_alt ada10 > logs/ind_wt_ada.out 2>&1; echo "ind weight_alt ada10 $?" >> $L ) &
wait
Rscript scripts/10_individual_schedules.R __none__ > /dev/null 2>&1 || true
# 2) 시험 수준 일정 비교(4코어)
for v in base struct2020 iiv150 resid12 weight_alt ada10; do
  Rscript scripts/11_trial_schedules.R $v 500 4 > logs/trial_sched_$v.out 2>&1; echo "trial sched $v $?" >> $L
done
# 3) 제품 차이 시나리오
for v in base struct2020; do
  Rscript scripts/12_trial_products.R 500 4 $v > logs/trial_products_$v.out 2>&1; echo "trial products $v $?" >> $L
done
echo "ALL DONE $(date -u)" >> $L
