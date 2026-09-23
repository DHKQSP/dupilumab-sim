# dupilumab-sim

듀필루맙 바이오시밀러 Phase 1 comparative PK study에서 **AUC0-last를 1차 평가변수로 두는 것이
AUC0-inf 대비 정보 손실이 없음**을 보이는 시뮬레이션. 기준 문서는 [`SPEC.md`](SPEC.md), 결정 기록은 [`DECISIONS.md`](DECISIONS.md).

## 구성
```
SPEC.md / DECISIONS.md        기준 문서, 결정 기록 (모든 세션은 SPEC.md를 먼저 읽는다)
config/                       출처 주석이 달린 파라미터·설계·시나리오 (PENDING 항목은 사용자 입력 대기)
R/                            모델(rxode2), 가상 집단, 채혈, 관측 생성, NCA·BE 통계(직접 구현), 시험 엔진, 요약
tests/testthat/               단계 1 자동 테스트 등
scripts/                      run_tests.R, 02_validate_step1.R, 02b/02c 진단, 10–13 시뮬레이션, 14 후처리, 05_report.R, run_all.sh
report/report.Rmd             보고서
results/ , logs/              결과(개발 규모는 CSV 커밋), 실행 로그(시드·세션 정보·config 해시)
renv.lock                     패키지 버전 잠금
```

## 재현
```sh
# 1) R 4.3.3 + renv 로 패키지 복원 (CRAN 접근 가능한 환경)
Rscript -e 'install.packages("renv"); renv::restore()'
# 2) 전체 순서 (테스트 → 단계 1 → 진단 → 개인 수준 → 시험 수준 → 제품 시나리오 → 교차검증 → 보고서)
bash scripts/run_all.sh
```
- 개별 실행: `Rscript scripts/run_tests.R`, `02_validate_step1.R [k2016|k2020]`, `02b_step1_diagnostics.R`, `02c_gate_sampling_error.R`, `10_individual_schedules.R [변형...]`, `11_trial_schedules.R [변형] [시험 수] [코어]`, `12_trial_products.R [시험 수] [코어] [변형]`, `13_crossval.R [시험 수] [코어]`, `14_postprocess.R`, `05_report.R`. 단계 2–4 일괄: `bash scripts/run_step2_all.sh`.
- `renv::restore()` 없이 시스템 라이브러리로 실행하려면(CRAN 차단 환경) `RENV_CONFIG_EXTERNAL_LIBRARIES="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"` 를 설정한다. 최종본 전 정리 대상(D-017).
- 단계 1 gate 판정은 `config/gate_decision.yaml`에 기록한다. 보고서는 판정이 `passed` 또는 `accepted_by_reviewer`일 때만 단계 2 이후 결과를 렌더링한다.
- `scripts/01_calibrate_ka.R`은 보관용이다(실행 시 중단, D-015).

## 상태
SPEC.md §4.4(단계 1 결과)와 §11(미결 항목) 참조.
