# dupilumab-sim

듀필루맙 바이오시밀러 Phase 1 comparative PK study에서 **AUC0-last를 1차 평가변수로 두는 것이
AUC0-inf 대비 정보 손실이 없음**을 보이는 시뮬레이션. 공동 1차 평가변수 구성(P2 = AUClast + Cmax, fallback F3, 가이드라인 기본 G2)의
운용특성은 결과 생성 전에 고정한 설계(`config/oc_design.yaml`)로 평가한다. 기준 문서는 [`SPEC.md`](SPEC.md), 결정 기록은 [`DECISIONS.md`](DECISIONS.md).

## 구성
```
SPEC.md / DECISIONS.md        기준 문서, 결정 기록 (모든 세션은 SPEC.md를 먼저 읽는다)
config/                       출처 주석이 달린 파라미터·설계·시나리오 (PENDING 항목은 사용자 입력 대기)
R/                            모델(rxode2), 가상 집단, 채혈, 관측 생성, NCA(Phoenix WinNonlin 호환, NonCompart 대조 검증)·BE 통계, 시험 엔진, 운용특성(oc.R), 요약
tests/testthat/               단계 1 자동 테스트 등
scripts/                      run_tests.R, 02* 단계 1, 10–13 시뮬레이션, 14 후처리, 15–27 요약·추가 분석, 28–29 NCA 엔진 검증,
                              30–33 운용특성, 34–36 절벽 분석, 37 핵심 수치(영문), 05_report.R, run_all.sh, 회차별 드라이버
report/report.Rmd             보고서
results/ , logs/              결과(개발 규모는 CSV 커밋), 실행 로그(시드·세션 정보·config 해시)
renv.lock                     패키지 버전 잠금
```

## 재현
```sh
# 1) R 4.3.3 + renv 로 패키지 복원 (CRAN 접근 가능한 환경)
Rscript -e 'install.packages("renv"); renv::restore()'
# 2) 전체 순서 (테스트 → NCA 엔진 검증 → 단계 1 → 개인·시험 수준 → 교차검증 → 5,000회·20만 명 → 추가 분석 → 절벽 → 운용특성 → 요약 → 보고서; 4코어 약 30시간)
bash scripts/run_all.sh
```
- 개별 실행: `Rscript scripts/run_tests.R`, `02_validate_step1.R [k2016|k2020]`, `02b_step1_diagnostics.R`, `02c_gate_sampling_error.R`, `02d_absorption_diagnostic.R`(부록, 채택 아님), `10_individual_schedules.R [변형...]`, `11_trial_schedules.R [변형] [시험 수] [코어]`, `12_trial_products.R [시험 수] [코어] [변형]`, `13_crossval.R [시험 수] [코어]`, `14_postprocess.R`, `15_rationale_summary.R`, `05_report.R`. 단계 2–4 일괄: `bash scripts/run_step2_all.sh`.
- 통합본(2026-09-24) 추가: `02e_arm_weight_sensitivity.R [k2016|k2020]`(Li 2020 arm 가정 체중), `13b_crossval_model1.R`, `16_criterion_d_200k.R <변형>`(기준 d 20만 명), `17_literature_table.R`, `18_summary_en.R`(영문 요약 `results/summary_en.md`), `19_mc_consistency.R`, `20_products_5000.R [코어] [변형]`(5,000회 + 적응적 상향, 500회 묶음 체크포인트·재시작), `21_curve_shape.R`, `22_pillar2_curvature.R [코어]`, `23_fallback_analyses.R`, `24_weight_generalization.R [모델...] [nojitter] [밴드...]`, `25_weight_trials.R <a|b|d> [코어]`, `26_schedule_extras.R`(B−, 폭 분해; `11_trial_schedules.R noresid` 선행), `27_reviewer_reference_round5.R`. 병렬 드라이버 `bash scripts/run_round5.sh`.
- 검토 의견 통합(2026-09-24, 6회차) 추가: `28_nca_engine_validation.R`(자체 엔진 대 NonCompart 0.8.4·PKNCA 0.12.1, Theoph·Indometh·듀필루맙 1,000명), `29_nca_engine_difference.R`·`35_engine_difference_trials.R`(이전 엔진 대비 차이표; 이전 값은 `results/nca_engine/legacy_snapshot/`), `30_oc_inversion.R <k2016|k2020> prep|<기전...>`(목표 참값 비 역산, 200,000명 공통 난수), `31_oc_trials.R <모델> [코어]`(경계 10,000회·기타 2,000회), `32_oc_random_space.R <모델> <truth|trials> [코어]`(라틴 하이퍼큐브 20,000개 제품), `33_oc_summary.R`(경계 1종 오류, 검정력, 위험, 구성 비교, 그림 3-A~3-D, 자동 결론 문구), `34_cliff_analysis.R`·`34b_cliff_figures.R`·`36_cliff_conclusion.R`(절벽 채혈 분석, 그림 2-1~2-4), `37_key_numbers_en.R`(`results/key_numbers_en.md`). 드라이버 `run_round6_rerun.sh`(새 엔진으로 NCA 의존 산출물 재실행), `run_oc_inversion.sh`, `run_oc_followup.sh`. 그림은 한국어판(보고서)과 영문판(`*_en.png`, `summary_en.md`) 두 벌을 만든다.
- `renv::restore()` 없이 시스템 라이브러리로 실행하려면(CRAN 차단 환경) `RENV_CONFIG_EXTERNAL_LIBRARIES="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"` 를 설정한다. 최종본 전 정리 대상(D-017).
- 단계 1 gate 판정은 `config/gate_decision.yaml`에 기록한다(2026-09-23 옵션 1). 보고서는 `decision_option`이 1 또는 2일 때 단계 2 이후 결과를 렌더링한다.
- `scripts/01_calibrate_ka.R`은 보관용이다(실행 시 중단, D-015).

## 상태
SPEC.md §4.4(단계 1 결과)와 §11(미결 항목) 참조.
