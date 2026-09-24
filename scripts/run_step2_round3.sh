#!/bin/bash
# 검토 의견 3차 추가 계산만 실행(기존 결과 재사용). 전체 재현은 run_step2_all.sh.
cd "$(dirname "$0")/.."
L=logs/step2_round3_$(date -u +%Y%m%dT%H%M%S).log
echo "start $(date -u)" > $L
( Rscript scripts/10_individual_schedules.R vmax080_both > logs/ind_vmax080.out 2>&1; echo "ind vmax080_both $?" >> $L ) &
( Rscript scripts/10_individual_schedules.R vmax125_both > logs/ind_vmax125.out 2>&1; echo "ind vmax125_both $?" >> $L ) &
wait
Rscript scripts/10_individual_schedules.R __none__ > /dev/null 2>&1 || true
for v in base struct2020 vmax080_both vmax125_both iiv150 resid12 weight_alt ada10; do
  Rscript scripts/11_trial_schedules.R $v 500 4 > logs/trial_sched_$v.out 2>&1; echo "trial sched $v $?" >> $L
done
echo "ALL DONE $(date -u)" >> $L
