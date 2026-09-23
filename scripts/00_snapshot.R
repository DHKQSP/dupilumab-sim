#!/usr/bin/env Rscript
# renv 잠금 파일 생성/갱신. 프로젝트 코드(R/, scripts/, tests/, report/)가 참조하는 패키지와 그 의존성을 기록한다.
# 사용법: Rscript scripts/00_snapshot.R
if (!file.exists("renv/activate.R")) renv::scaffold(project = ".", settings = list(snapshot.type = "implicit"))
renv::snapshot(project = ".", type = "implicit", library = .libPaths(), prompt = FALSE, force = TRUE)
lk <- jsonlite::fromJSON("renv.lock")
cat(sprintf("renv.lock: R %s, %d packages\n", lk$R$Version, length(lk$Packages)))
