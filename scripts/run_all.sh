#!/bin/bash
# 전체 재현 순서(검토 의견 통합본 2026-09-24 기준). CRAN 차단 환경에서는 RENV_CONFIG_EXTERNAL_LIBRARIES 설정 후 실행(README).
# 4코어 기준 약 5–6시간. 이번 회차의 병렬 드라이버는 scripts/run_round5.sh(같은 스크립트를 두 대기열로 실행).
set -e
cd "$(dirname "$0")/.."
Rscript scripts/run_tests.R || echo "[경고] 테스트 실패 — logs/tests_*.log 확인. 단계 1 gate 미달이면 config/gate_decision.yaml 참조"
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
# 후처리·요약·보고서
Rscript scripts/14_postprocess.R
Rscript scripts/26_schedule_extras.R
Rscript scripts/19_mc_consistency.R
Rscript scripts/15_rationale_summary.R
Rscript scripts/17_literature_table.R
Rscript scripts/23_fallback_analyses.R
Rscript scripts/27_reviewer_reference_round5.R
Rscript scripts/18_summary_en.R
Rscript scripts/05_report.R
