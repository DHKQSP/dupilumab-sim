#!/bin/bash
# 6회차 조립: 재실행·운용특성 계산이 끝난 뒤 요약·차이표·문구·보고서를 순서대로 다시 만든다(시뮬레이션 없음).
cd "$(dirname "$0")/.."
set -e
for s in 14_postprocess 26_schedule_extras 19_mc_consistency 15_rationale_summary 17_literature_table 23_fallback_analyses 27_reviewer_reference_round5 \
         35_engine_difference_trials 36_cliff_conclusion 33_oc_summary 37_key_numbers_en 18_summary_en 05_report; do
  echo "== $s $(date -u)"; Rscript scripts/$s.R > logs/r6_assembly_$s.out 2>&1 || { echo "FAILED $s (logs/r6_assembly_$s.out)"; exit 1; }
done
echo "assembly done $(date -u)"
