#!/bin/bash
# 전체 재현 순서(검토 의견 통합 2026-09-24: Phoenix 호환 NCA 엔진, 절벽 분석, 운용특성 기준). CRAN 차단 환경에서는 RENV_CONFIG_EXTERNAL_LIBRARIES 설정 후 실행(README).
# 4코어 기준 약 30시간(운용특성 역산·시험 반복·무작위 제품 공간이 대부분). 이번 회차의 병렬 드라이버:
#   scripts/run_round6_rerun.sh(새 엔진으로 NCA 의존 산출물 재실행), scripts/run_oc_inversion.sh(역산 12건), scripts/run_oc_followup.sh(31→32→33).
set -e
cd "$(dirname "$0")/.."
Rscript scripts/run_tests.R || echo "[경고] 테스트 실패 — logs/tests_*.log 확인. 단계 1 gate 판정성 항목은 종료 코드에 넣지 않으므로 로그의 judgment items 줄과 config/gate_decision.yaml 참조"
# NCA 엔진 검증(§1): 자체 엔진 대 NonCompart·PKNCA — Theoph, Indometh, 듀필루맙 1,000명. 교체 전 엔진 산출물 스냅숏은 results/nca_engine/legacy_snapshot/
Rscript scripts/28_nca_engine_validation.R
Rscript scripts/29_nca_engine_difference.R
# 단계 1: gate(연구별 채혈 일정, 내부/외부, 완전 외부 재판정), 진단, 흡수 부록, Li 2020 arm 가정 체중 민감도
Rscript scripts/02_validate_step1.R
Rscript scripts/02_validate_step1.R k2020
Rscript scripts/02b_step1_diagnostics.R
Rscript scripts/02c_gate_sampling_error.R
Rscript scripts/02d_absorption_diagnostic.R
Rscript scripts/02e_arm_weight_sensitivity.R
Rscript scripts/02e_arm_weight_sensitivity.R k2020
# 단계 2–4: 개인 수준(20,000명), 일정 비교(500회), 제품 시나리오(500회) — 모든 변형
bash scripts/run_step2_all.sh
Rscript scripts/11_trial_schedules.R noresid 500 4              # CI 폭 확대 원인 분해(잔차 없음)
# 교차검증
Rscript scripts/13_crossval.R
Rscript scripts/13b_crossval_model1.R
# Monte Carlo 정밀도(§4): 주 모델 제품 시나리오 5,000회(+적응적 상향), 기준 (d) 20만 명
Rscript scripts/20_products_5000.R 4 base
for v in base struct2020 vmax080_both; do Rscript scripts/16_criterion_d_200k.R $v; done
# 추가 분석(§5–§6)
Rscript scripts/21_curve_shape.R
Rscript scripts/22_pillar2_curvature.R 4
Rscript scripts/24_weight_generalization.R
Rscript scripts/24_weight_generalization.R a nojitter 40-60 60-75 75-90 130-150   # 독립 구현 조건(채혈 편차 없음) 비교용
for m in a b d; do Rscript scripts/25_weight_trials.R $m 4; done
# 절벽 채혈 분석(§2): 순간 반감기 < 1일(및 < 2일) 구간, 두 모델 × 체중 밴드, 20,000명, 명목일·허용창
Rscript scripts/34_cliff_analysis.R
Rscript scripts/34b_cliff_figures.R
Rscript scripts/36_cliff_conclusion.R
# 운용특성(§3–4): 설계는 config/oc_design.yaml(결과 전 커밋). 역산 → 시험 반복 → 무작위 제품 공간 → 요약
for m in k2016 k2020; do Rscript scripts/30_oc_inversion.R $m prep; done   # 공통 난수 200,000명 대조 참값(results/oc/truth_ref_<model>.rds, 저장소 제외)
bash scripts/run_oc_inversion.sh
bash scripts/run_oc_followup.sh
# 통합 지시 2026-09-25: AUCinf 규칙 A·C의 플래그 세트 (i) 재판정(같은 시드 재생성, 저장 행 대조), P2 경계 적응적 연장, 신뢰 기준 두 세트
for m in k2016 k2020; do Rscript scripts/40_oc_rejudge.R $m 3; done
Rscript scripts/33_oc_summary.R              # 연장 판정에 쓰는 10,000회 boundary_type1.csv
Rscript scripts/42_oc_extend.R 3 20000       # Wilson 구간이 5%를 포함한 P2 경계 시나리오만 20,000회로
Rscript scripts/39_reliability_flags.R
# 추가 지시 2026-09-26(config/prereg_20260926.yaml): 분석 모형 M0·M1·M2 재판정(12·20·40·42 산출과 대조하므로 그 뒤), LLOQ 민감도, 표본 수
for m in k2016 k2020; do Rscript scripts/44_oc_models.R $m 3 main; done
for m in k2016 k2020; do Rscript scripts/44_oc_models.R $m 3 ext; done   # 사전 등록 연장 규칙(M0 또는 M1)으로 고른 경계 칸만
Rscript scripts/45_oc_models_summary.R
for m in k2016 k2020; do Rscript scripts/46_lloq_sensitivity.R individual $m; done
Rscript scripts/46_lloq_sensitivity.R cliff
Rscript scripts/46_lloq_sensitivity.R trial 3
Rscript scripts/48_lloq_summary.R
Rscript scripts/47_sample_size.R stat 3
for m in k2016 k2020; do Rscript scripts/47_sample_size.R pk $m 3; done
Rscript scripts/49_sample_size_summary.R
Rscript scripts/50_sap_support.R
# 후처리·요약·보고서
Rscript scripts/14_postprocess.R
Rscript scripts/26_schedule_extras.R
Rscript scripts/19_mc_consistency.R
Rscript scripts/15_rationale_summary.R
Rscript scripts/17_literature_table.R
Rscript scripts/23_fallback_analyses.R
Rscript scripts/27_reviewer_reference_round5.R
Rscript scripts/35_engine_difference_trials.R
Rscript scripts/33_oc_summary.R
Rscript scripts/41_oc_rules_flags.R
Rscript scripts/43_p2_interpretation.R
Rscript scripts/37_key_numbers_en.R
Rscript scripts/18_summary_en.R
Rscript scripts/05_report.R
Rscript scripts/60_regulatory_package.R        # 규제 패키지(regulatory/): M&S 보고서, 예상 질의응답, SAP 제안 문안, 추적표, 해시 목록
