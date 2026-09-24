#!/bin/bash
# 전체 재현 순서. CRAN 차단 환경에서는 RENV_CONFIG_EXTERNAL_LIBRARIES 설정 후 실행(README).
set -e
cd "$(dirname "$0")/.."
Rscript scripts/run_tests.R || echo "[경고] 테스트 실패 — logs/tests_*.log 확인. 단계 1 gate 미달이면 config/gate_decision.yaml 참조"
Rscript scripts/02_validate_step1.R
Rscript scripts/02_validate_step1.R k2020
Rscript scripts/02b_step1_diagnostics.R
Rscript scripts/02c_gate_sampling_error.R
Rscript scripts/02d_absorption_diagnostic.R
bash scripts/run_step2_all.sh
Rscript scripts/13_crossval.R
Rscript scripts/15_rationale_summary.R
Rscript scripts/05_report.R
