# dupilumab-sim

듀필루맙 바이오시밀러 Phase 1 comparative PK study에서 **AUC0-last를 1차 평가변수로 두는 것이
AUC0-inf 대비 정보 손실이 없음**을 보이는 시뮬레이션. 기준 문서는 [`SPEC.md`](SPEC.md), 결정 기록은 [`DECISIONS.md`](DECISIONS.md).

## 구성
```
SPEC.md / DECISIONS.md        기준 문서, 결정 기록 (모든 세션은 SPEC.md를 먼저 읽는다)
config/                       출처 주석이 달린 파라미터·설계·시나리오 (PENDING 항목은 사용자 입력 대기)
R/                            모델(rxode2), 가상 집단, 채혈, 관측 생성, BEmaster 어댑터, 시험 루프, 2층 MC, 요약
tests/testthat/               단계 1 자동 테스트 등
scripts/                      run_tests.R, 01_calibrate_ka.R, 02_validate_step1.R, 03_run_dev.R, 04_run_scenarios.R, 05_report.R
report/report.Rmd             보고서
results/ , logs/              결과(개발 규모는 CSV 커밋), 실행 로그(시드·세션 정보·config 해시)
renv.lock                     패키지 버전 잠금
```

## 재현
```sh
# 1) R 4.3.3 + renv 로 패키지 복원 (CRAN 접근 가능한 환경)
Rscript -e 'install.packages("renv"); renv::restore()'
# 2) 자동 테스트 (단계 1 포함)
Rscript scripts/run_tests.R
# 3) 단계 1 요약, ka 보정, 개발 규모 루프
Rscript scripts/02_validate_step1.R dev
Rscript scripts/02b_step1_sensitivity.R      # 단계 1 미달 원인 스캔(진단)
Rscript scripts/01_calibrate_ka.R dev
Rscript scripts/03_run_dev.R dev S0 d57_base
Rscript scripts/04_run_scenarios.R dev
Rscript scripts/05_report.R
```
- `renv::restore()` 없이 시스템 라이브러리로 바로 실행하려면(예: CRAN 차단 환경) `RENV_CONFIG_EXTERNAL_LIBRARIES="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"` 를 설정한다.
- `run_mode`는 `dev`(개발용 자리표시자 허용) 또는 `final`(출처 미확정 값이 있으면 중단). 결과 보고는 `final`만 쓴다.
- BEmaster: 패키지로 설치하거나 `BEMASTER_PATH`로 위치를 지정한 뒤 `R/nca_bemaster.R`, `R/be_bemaster.R`의 바인딩 함수를 BEmaster 시그니처에 맞춰 채운다(규칙 재구현 금지, SPEC §6).

## 상태
`SPEC.md` §8 미결 항목(Q1–Q9)이 채워지기 전까지 단계 2 이후 결과는 보고하지 않는다.
