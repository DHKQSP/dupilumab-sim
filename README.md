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
.github/workflows/            CI: tests.yml(fast 범위 테스트), repro.yml(수동 재현성 점검)
scripts/                      run_tests.R, 02* 단계 1, 10–13 시뮬레이션, 14 후처리, 15–27 요약·추가 분석, 28–29 NCA 엔진 검증,
                              30–33 운용특성, 34–36 절벽 분석, 37 핵심 수치(영문), 38 재현성 점검, 05_report.R, run_all.sh, 회차별 드라이버
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

## CI
- **tests** (`.github/workflows/tests.yml`): push·pull_request 중 `R/**`, `tests/**`, `config/**`, `scripts/**`, `renv.lock`, `renv/**`(activate.R·settings.json), `.Rprofile`, `.github/**`가 바뀐 경우에만 실행한다. `results/**`, `logs/**`, `*.md`만 바꾼 커밋은 실행하지 않는다. 수동 실행(workflow_dispatch)도 된다. 같은 ref에 새 push가 오면 진행 중인 실행을 취소한다(concurrency). renv 라이브러리는 `r-lib/actions/setup-renv@v2`가 renv.lock 해시를 키로 캐시한다(성공한 실행에서 저장). 제한 시간 30분, 권한 `contents: read`. 테스트 로그(`tests_<UTC>.log`, `tests_<UTC>_timing.csv`)는 artifact `test-logs`로 올린다.
- **테스트 범위** `DUPI_TEST_SCOPE`: `fast`(CI) = 모델 구조(해석해 일치), NCA 엔진(손계산, Theoph·Indometh·모의 프로필 NonCompart 대조), config 스키마, 판정 규칙, BE 통계, 소규모 스모크 모의, OC 구성 요소. `full`(로컬 기본값) = fast + 느린 테스트(`skip_if_fast_scope`) + 판정성 테스트. 실행: `DUPI_TEST_SCOPE=fast Rscript scripts/run_tests.R`. 파일별 소요 시간은 `results/ci/test_timing_local.csv`. CI는 `DUPI_TEST_STRICT_SKIP=true`로 실행해, 위 두 종류 외의 건너뜀(예: rxode2·NonCompart를 불러오지 못해 `skip_if_not_installed`가 건너뜀)을 실패로 센다.
- **판정성 테스트**: 통과 여부가 모의 결과에 달린 과학적 판정 항목이다. 현재 `tests/testthat/test-step1-clot2021.R`의 단계 1 gate (a)(b)(c)(f)이며 이름이 `[판정]`으로 시작한다. CI에서는 실행하지 않는다(`skip_if_judgment`). full 범위에서는 실행하되, 기대값 불충족은 `scripts/run_tests.R`가 "judgment items"로 따로 보고하고 종료 코드에 넣지 않는다(코드 오류는 실패로 센다). 정식 판정 기록은 `scripts/02_validate_step1.R`가 만드는 `results/step1/` 결과 파일이다.
- **재현성 점검** (`.github/workflows/repro.yml`, 수동 전용): GitHub 저장소의 Actions 탭 > repro > Run workflow, 또는 `gh workflow run repro.yml --ref <브랜치>`. 깨끗한 renv 복원 → fast 테스트 → `scripts/38_repro_check.R` 순서로 실행한다. 비교 항목은 (a) 커밋된 OC 시험(k2016 시험 1–100, S00·Vmax_down_125)의 정확 재생성, (b) Vmax 역산 배율의 참 AUC0-inf 비(공통 난수 앞 20,000명 대 200,000명), (c) B0 개인 참 외삽 비율 중앙값(2,000명 대 20,000명)이다. 허용 오차는 `config/repro_check.yaml`에 사전 고정했다. (a)의 기준은 `results/oc/oc_trials_be_k2016.csv.gz`이며, 이 파일이 저장소에 없으면 발췌본 `results/repro/reference/`(`--extract-reference`로 생성, 둘 다 있으면 일치 확인)을 쓴다. 결과 `results/repro/repro_check.csv`와 `repro_check_meta.csv`(커밋 해시, R·rxode2 버전, 플랫폼, 시각)는 artifact `repro-results`와 작업 요약에 남는다. 하나라도 불통과면 실행이 실패한다. 로컬 실행: `Rscript scripts/38_repro_check.R [코어]`.
- CI 실패 이력과 Actions 사용 시간(분) 추정: `results/ci/`(`ci_failures_before.csv`, `ci_failure_summary.csv`, `ci_minutes_estimate.csv`).

## 상태
SPEC.md §4.4(단계 1 결과)와 §11(미결 항목) 참조.
