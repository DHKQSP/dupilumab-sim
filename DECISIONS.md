# DECISIONS.md — 결정 기록

형식: `D-번호 | 날짜 | 결정 | 근거 | 영향 파일`. 파라미터·가정 변경은 SPEC.md와 같은 커밋에서 갱신.

## D-001 | 2026-09-23 | 실행 환경: R 4.3.3(apt) + rxode2 5.1.7(CRAN 릴리스, GitHub `cran/*` 미러에서 소스 설치)
- 근거: 이 세션의 클라우드 환경은 CRAN·Posit Package Manager·conda 접근이 정책상 차단(403)되어 있고, Ubuntu apt와 GitHub git 읽기만 가능했다. `cran/<pkg>` 미러의 태그는 CRAN 릴리스와 동일 소스이므로 버전 고정 재현성에 문제 없음. 사용자 환경에서는 `renv::restore()`로 CRAN에서 동일 버전을 받으면 된다.
- 영향: `renv.lock`, `README.md`.
- 보충(같은 날): apt의 Rcpp 1.0.12는 rxode2 5.1.7이 쓰는 가변 인자 `NumericVector::create`를 지원하지 않고, Debian 재패키징 BH는 Boost odeint 헤더가 없어 configure가 중단된다. Rcpp 1.1.2, BH 1.90.0-1, RcppParallel 6.2.1을 CRAN 릴리스에서 설치. RcppParallel은 시스템 TBB를 쓰되 `TBB_INC`에 `/usr/include`를 직접 주면 `#include_next` 체인이 깨지므로 심볼릭 링크 디렉터리(`/opt/tbb-inc`)를 경유.

## D-002 | 2026-09-23 | 모델 식에서 Vmax(mg/L/day)에 Vc를 곱해 양(mg) 기준 소실 속도로 변환
- 근거: Kovalenko 2016은 Vmax를 농도 단위(mg/L/day)로 보고. 상태변수를 양으로 두므로 단위 정합을 위해 `Vmax·C/(Km+C)·Vc`.
- 영향: `R/model.R`, SPEC §3.1.

## D-003 | 2026-09-23 | 검증 목표의 "마지막 정량 시점 42/28/56일"은 투여 후 경과일로 해석(스터디 Day 43/29/57)
- 근거: 채혈표(Day 2…57)에 56일이 없어 스터디 Day 해석이 불가능. 투여 Day 1 기준 경과일이면 세 값 모두 채혈 시점과 일치. 사용자 확인 요청(SPEC Q5).
- 영향: `config/design_clot2021.yaml`, `tests/testthat/test-step1-clot2021.R`.

## D-004 | 2026-09-23 | NCA·BE 판정은 BEmaster 어댑터(`R/nca_bemaster.R`, `R/be_bemaster.R`)만 통해 수행. 자체 NCA 미구현
- 근거: 과제 원칙("비구획분석 규칙을 BEmaster와 다르게 구현" 금지). BEmaster 미입수 상태에서는 해당 테스트를 skip하고 파이프라인을 NCA 경계에서 중단. 모델 기반 진적분 AUC는 "truth" 참값으로만 사용하며 NCA 대체가 아님.
- 영향: SPEC §6.

## D-005 | 2026-09-23 | 출처 미확정 값은 `[PENDING]`으로 두고, 개발 실행에만 `config/dev_assumptions.yaml`의 `[DEV-가정]` 값을 허용. `run_mode = "final"`에서는 PENDING 값이 남아 있으면 즉시 중단
- 근거: "출처 없는 값 사용 금지"와 "파이프라인 개발 진행"의 양립. 최종 결과에 임시값이 섞이는 것을 코드로 차단.
- 영향: `R/params.R`, `config/dev_assumptions.yaml`.

## D-006 | 2026-09-23 | 난수 시드 파생 규칙: master seed → 바깥층 i → 안쪽층 j를 `digest` 기반 결정적 해시로 파생
- 근거: 병렬 실행·부분 재실행에서도 동일 시드 재현. 실행 로그에 master seed와 파생 규칙 버전을 기록.
- 영향: `R/seeds.R`.
