#!/bin/bash
# §3 후속 계산 드라이버: 모델별 역산 6건 완료 → 31(시험 반복) → 32(무작위 공간 참값·시험) → 33(요약).
cd "$(dirname "$0")/.."
L=logs/oc_followup_driver.log; echo "start $(date -u)" >> $L
wait_inv() { while [ $(ls results/oc/inversion_$1_*.csv 2>/dev/null | grep -v scan | wc -l) -lt 6 ]; do sleep 60; done; }
for m in k2016 k2020; do
  wait_inv $m; echo "inversion $m complete $(date -u)" >> $L
  Rscript scripts/31_oc_trials.R $m 3 > logs/oc_trials_$m.out 2>&1; echo "trials $m exit $? $(date -u)" >> $L
done
for m in k2016 k2020; do
  Rscript scripts/32_oc_random_space.R $m truth 3 > logs/oc_random_truth_$m.out 2>&1; echo "random truth $m exit $? $(date -u)" >> $L
  Rscript scripts/32_oc_random_space.R $m trials 3 > logs/oc_random_trials_$m.out 2>&1; echo "random trials $m exit $? $(date -u)" >> $L
done
Rscript scripts/33_oc_summary.R > logs/oc_summary.out 2>&1; echo "summary exit $? $(date -u)" >> $L
echo "ALL DONE $(date -u)" >> $L
