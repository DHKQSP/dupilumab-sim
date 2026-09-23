# results/

| 폴더 | 내용 | 보고 상태 |
|---|---|---|
| `step1/` | 단계 1 gate(기본 모델), 진단(`diag_*`) | 보고 가능 |
| `step1_k2020/` | 단계 1 gate(Kovalenko 2020 Model 1) | 보고 가능 |
| `crossval/` | 검토자 독립 구현 대비 교차검증(코드 검증) | 보고 가능 |
| `individual/` | 단계 2 개인 수준(20,000명 × 일정 × 변형) | **보고 보류** — `config/gate_decision.yaml` 판정 전 |
| `trials/` | 단계 2–3 시험 수준(500회), 원자료 `*_be_raw*.csv.gz` | **보고 보류** — 같음 |

- 보고 보류 폴더는 계산 재현과 검토를 위해 저장한 것이며, 단계 1 gate 판정(SPEC §1, §4.4) 전에는 결론으로 인용하지 않는다.
- `*.rds`(개체 수준 NCA 원자료)는 용량 때문에 git에서 제외한다. 스크립트와 시드로 재생성된다(`scripts/10_individual_schedules.R`).
- 시험 수준 요약은 원자료에서 `scripts/14_postprocess.R`로 다시 만들 수 있다.
