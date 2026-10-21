#!/usr/bin/env bash
# run_m15.sh — 계획서 11절 순서대로 [2]~[8] 실행(D-068). 단계마다 진행 줄을 results/m15/run_progress.log에 남긴다.
# 사용법: bash scripts/run_m15.sh <stage> ; stage = prep | pilot | full | post | a1prep | a1pilot | a1full | a1post | a1p1 | a1p2a | a1p2pilot | a1p2full | a1final(개정 1, 계획서 16.12절) | tm(추록 1, 계획서 17절)
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
  a1prep)    # 개정 1(계획서 16.12절 3단계): Day 57·64·71 정량 몫, 항약물항체 후보일 농도, 분산 분해(저장 결과와 일치 점검)
    R scripts/72_m15_individual.R amend 4
    ;;
  a1pilot)   # 16.12절 4단계 예비: mechsuppx Model 1 S1 시험 1~50 → 일치 점검(불일치면 종료 코드 3)
    R scripts/73_m15_trials.R mechsuppx k2020 S1 300 1 50 4
    R scripts/77_m15_amend.R mechcheck k2020 S1
    ;;
  a1full)    # 16.12절 4단계 본 실행: 네 칸 1~1,000(이미 쓴 묶음은 건너뜀) → 일치 점검 → 참값 확장
    for m in k2020 k2016; do for pop in S1 S2; do R scripts/73_m15_trials.R mechsuppx $m $pop 300 1 1000 4; done; done
    R scripts/77_m15_amend.R mechcheck
    R scripts/71_m15_populations.R truthx 4
    ;;
  a1post)    # 16.12절 6단계: 표기 반영 사후 단계 재실행(깨끗한 커밋) → 재실행 전 사본과 비교(숫자 열 동일)
    B=results/m15/amend1/pre_rerun_copy; rm -rf "$B"; mkdir -p "$B"
    for d in summary samplesize weight figures; do mkdir -p "$B/$d"; cp results/m15/$d/*.csv "$B/$d/"; done
    [ -f results/m15/individual/individual_subjects.rds ] || { echo "$(date -u +%FT%TZ) STOP individual_subjects.rds missing" >> "$P"; exit 3; }
    R scripts/71_m15_populations.R weight 4
    R scripts/74_m15_samplesize.R 4
    R scripts/75_m15_summary.R tables 4
    R scripts/75_m15_summary.R figures
    R scripts/75_m15_summary.R summary
    R scripts/77_m15_amend.R postcheck "$B"
    ;;
  a1p1)      # 16.12절 7~8단계: 교차 확장 → P1 추출 → 저장 반복 사후 예측 점검 표 → P1 보고 → 부록 초안
    R scripts/77_m15_amend.R mechx
    R scripts/77_m15_amend.R extract 4
    R scripts/77_m15_amend.R ppcstored
    R scripts/77_m15_amend.R p1report
    R scripts/77_m15_amend.R appendix
    ;;
  a1p2a)     # 16.12절 9단계: 요청자 NCA 규칙 반복(등록 두 모델, 약 10분) → 사후 예측 점검 진단(위양성률, 연구별 규칙, 결론 1 고정 서식)
    R scripts/70_m15_validation.R 4 alt
    R scripts/77_m15_amend.R ppc
    ;;
  a1p2pilot) # 16.12절 10단계 예비: 세트 0 재현 → Monte Carlo 잡음 → Model 1 차용 세트 0~49(첫 묶음, 소요 시간 추정; 본 실행이 같은 묶음 파일을 건너뜀)
    R scripts/76_m15_param_assurance.R check0 4
    R scripts/76_m15_param_assurance.R mcnoise 4
    R scripts/76_m15_param_assurance.R sets k2020_b16 0 49 4
    ;;
  a1p2full)  # 16.12절 10단계 본 실행: 세 세트 0~1,000(세트 0 포함, 이미 쓴 묶음은 건너뜀) → 요약
    for sid in k2020_b16 k2016 k2016_thetas; do R scripts/76_m15_param_assurance.R sets $sid 0 1000 4; done
    R scripts/76_m15_param_assurance.R summary
    ;;
  a1final)   # 16.12절 11단계: robustness 표, 적용 가능성, 요약 표, 결론 → 요약 단계(가정·이탈 표)
    R scripts/77_m15_amend.R final
    R scripts/75_m15_summary.R summary
    ;;
  tm)        # 추록 1(계획서 17절, D-071): f_TM 계산(약 20~30분) → 표 2개·그림 2개 → 가정·이탈 표 갱신(요약 단계 재실행, 숫자 열 동일 확인)
    R scripts/78_m15_tm_fraction.R run 4
    R scripts/78_m15_tm_fraction.R report
    B=results/m15/tm/pre_summary_copy; rm -rf "$B"; mkdir -p "$B/summary"; cp results/m15/summary/*.csv "$B/summary/"
    R scripts/75_m15_summary.R summary
    R scripts/77_m15_amend.R postcheck "$B" results/m15/tm/summary_rerun_check.csv
    ;;
  *) echo "usage: bash scripts/run_m15.sh <prep|pilot|full|post|a1prep|a1pilot|a1full|a1post|a1p1|a1p2a|a1p2pilot|a1p2full|a1final|tm>"; exit 1 ;;
esac
echo "$(date -u +%FT%TZ) stage ${1} complete" >> "$P"
