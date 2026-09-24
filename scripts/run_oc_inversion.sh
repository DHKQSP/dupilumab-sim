#!/bin/bash
# §3-2 역산: 모델 × 기전 12건을 3개 동시 실행(xargs -P 3). 로그: logs/oc_inv_<model>_<mech>.out
cd "$(dirname "$0")/.."
for m in k2016 k2020; do for x in F ke Vmax ka V2 Km; do echo "$m $x"; done; done | \
  xargs -P 3 -L 1 bash -c 'Rscript scripts/30_oc_inversion.R $0 $1 > logs/oc_inv_$0_$1.out 2>&1; echo "$0 $1 exit $?" >> logs/oc_inversion_driver.log'
echo "ALL DONE $(date -u)" >> logs/oc_inversion_driver.log
