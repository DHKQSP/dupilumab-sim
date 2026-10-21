#!/usr/bin/env bash
# run_m15.sh — 계획서 11절 순서대로 [2]~[8] 실행(D-068). 단계마다 진행 줄을 results/m15/run_progress.log에 남긴다.
# 사용법: bash scripts/run_m15.sh <stage> ; stage = prep | pilot | full | post
#   prep  : 집단 → 참값 → 기존 엔진 대조(중단 규칙) → [6] → 부분 AUC 검증(중단 규칙)
#   pilot : 예비 실행(k2020 S1 300 시험 1~500, k2020 S2 200 시험 1~50; [3]·[4]) → 점검(실패 시 종료 코드 3)
#   full  : 본 실행(Model 1 칸 먼저, Kovalenko 2016 [4] 200 mg 마지막; 이미 쓴 묶음은 건너뜀) → 보충 → 확장 → 체중 모양
#   post  : [8] → [7] → 요약 표·그림·요약. 그림 5는 results/m15/individual/individual_subjects.rds(약 17 MB, .gitignore 제외)를 읽는다:
#           새로 받은 저장소에서는 먼저 nice -n 5 Rscript scripts/72_m15_individual.R run 4(약 2분, 결정적)를 실행한다
set -euo pipefail
cd "$(dirname "$0")/.."
P=results/m15/run_progress.log; mkdir -p results/m15
R() {
  local t0 rc; t0=$(date +%s); echo "$(date -u +%FT%TZ) start $*" >> "$P"
  rc=0; nice -n 5 Rscript "$@" > results/m15/last_stage.out 2>&1 || rc=$?
  if [ "$rc" -ne 0 ]; then echo "$(date -u +%FT%TZ) FAILED $* (exit $rc)" >> "$P"; tail -20 results/m15/last_stage.out >> "$P"; exit 3; fi
  echo "$(date -u +%FT%TZ) done $* in $(( $(date +%s) - t0 ))s" >> "$P"
}
case "${1:-}" in
  prep)
    R scripts/71_m15_populations.R populations
    R scripts/71_m15_populations.R truth 4
    R scripts/72_m15_individual.R legacy 2
    Rscript -e 'x <- read.csv("results/m15/individual/legacy_check.csv"); if (!all(x$agree)) quit(status = 3)' || { echo "$(date -u +%FT%TZ) STOP legacy check failed" >> "$P"; exit 3; }
    R scripts/72_m15_individual.R run 4
    R scripts/72_m15_individual.R verify_pauc
    Rscript -e 'x <- read.csv("results/m15/individual/pauc_pknca_check.csv"); if (!all(x$pass)) quit(status = 3)' || { echo "$(date -u +%FT%TZ) STOP partial AUC verification failed" >> "$P"; exit 3; }
    ;;
  pilot)
    R scripts/73_m15_trials.R oc k2020 S1 300 1 500 4
    R scripts/73_m15_trials.R mech k2020 S1 300 1 500 4
    R scripts/73_m15_trials.R oc k2020 S2 200 1 50 4
    R scripts/73_m15_trials.R mech k2020 S2 200 1 50 4
    R scripts/73_m15_trials.R pilotcheck k2020 S1 300 1 500
    R scripts/73_m15_trials.R pilotcheck k2020 S2 200 1 50
    ;;
  full)
    for m in oc mech; do for pop in S1 S2; do for d in 300 200; do R scripts/73_m15_trials.R $m k2020 $pop $d 1 1000 4; done; done; done
    for pop in S1 S2; do R scripts/73_m15_trials.R mechsupp k2020 $pop 300 1 1000 4; done
    for pop in S1 S2; do for d in 300 200; do R scripts/73_m15_trials.R oc k2016 $pop $d 1 1000 4; done; done
    for pop in S1 S2; do R scripts/73_m15_trials.R mech k2016 $pop 300 1 1000 4; done
    for pop in S1 S2; do for d in 300 200; do R scripts/73_m15_trials.R ocext k2020 $pop $d 1001 10000 4; done; done
    for pop in S1_uniform S2_uniform; do R scripts/73_m15_trials.R ocshape k2020 $pop 300 1 1000 4; done
    for pop in S1 S2; do R scripts/73_m15_trials.R mech k2016 $pop 200 1 1000 4; done
    ;;
  post)
    [ -f results/m15/individual/individual_subjects.rds ] || { echo "$(date -u +%FT%TZ) STOP individual_subjects.rds missing: run scripts/72_m15_individual.R run 4 first" >> "$P"; exit 3; }
    R scripts/71_m15_populations.R weight 4
    R scripts/74_m15_samplesize.R 4
    R scripts/75_m15_summary.R tables 4
    R scripts/75_m15_summary.R figures
    R scripts/75_m15_summary.R summary
    ;;
  *) echo "usage: bash scripts/run_m15.sh <prep|pilot|full|post>"; exit 1 ;;
esac
echo "$(date -u +%FT%TZ) stage ${1} complete" >> "$P"
