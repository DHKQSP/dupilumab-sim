#!/bin/bash
# 단계 2–4 전체 계산(검토 의견 3차 기준). 결과: results/individual, results/trials. 로그: logs/
cd "$(dirname "$0")/.."
L=logs/step2_driver_$(date -u +%Y%m%dT%H%M%S).log
echo "start $(date -u)" > $L
# 1) 개인 수준(20,000명 × 6일정): 변형을 4개 프로세스로 병렬
( Rscript scripts/10_individual_schedules.R base iiv150        > logs/ind_a.out 2>&1; echo "ind base iiv150 $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R struct2020 resid12 > logs/ind_b.out 2>&1; echo "ind struct2020 resid12 $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R weight_alt ada10   > logs/ind_c.out 2>&1; echo "ind weight_alt ada10 $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R vmax080_both vmax125_both > logs/ind_d.out 2>&1; echo "ind vmax080 vmax125 $?" >> $L ) &
wait
Rscript scripts/10_individual_schedules.R __none__ > /dev/null 2>&1 || true
# 2) 시험 수준 일정 비교(500회, 4코어)
for v in base struct2020 vmax080_both vmax125_both iiv150 resid12 weight_alt ada10; do
  Rscript scripts/11_trial_schedules.R $v 500 4 > logs/trial_sched_$v.out 2>&1; echo "trial sched $v $?" >> $L
done
# 3) 제품 차이 시나리오(B0, 14개 × 500회)
for v in base struct2020; do
  Rscript scripts/12_trial_products.R 500 4 $v > logs/trial_products_$v.out 2>&1; echo "trial products $v $?" >> $L
done
echo "ALL DONE $(date -u)" >> $L
