#!/bin/bash
# 전체 재현 순서. CRAN 차단 환경에서는 RENV_CONFIG_EXTERNAL_LIBRARIES 설정 후 실행(README).
set -e
cd "$(dirname "$0")/.."
Rscript scripts/run_tests.R
Rscript scripts/02_validate_step1.R
Rscript scripts/02b_step1_diagnostics.R
Rscript scripts/10_individual_schedules.R
Rscript scripts/11_trial_schedules.R base
for v in struct2020 iiv150 resid12 weight_alt ada10; do Rscript scripts/11_trial_schedules.R $v; done
Rscript scripts/12_trial_products.R
Rscript scripts/13_crossval.R
Rscript scripts/05_report.R
