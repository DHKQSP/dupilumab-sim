# 00_setup.R — 공통 초기화. 모든 스크립트/테스트는 이 파일을 먼저 source 한다.
# 로케일: config의 한글(UTF-8)을 읽기 위해 UTF-8 CTYPE를 보장한다 (컨테이너는 기본 C 로케일)
if (!grepl("UTF-8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  ok <- FALSE
  for (loc in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "ko_KR.UTF-8")) {
    r <- suppressWarnings(try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE))
    if (!inherits(r, "try-error") && nzchar(r)) { ok <- TRUE; break }
  }
  if (!ok) warning("UTF-8 로케일을 설정하지 못했습니다. config 한글 파싱이 실패할 수 있습니다.")
}
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(digest)
})

PROJ_ROOT <- local({
  # 프로젝트 루트: SPEC.md가 있는 디렉터리를 위로 탐색
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "SPEC.md"))) return(d)
    d <- dirname(d)
  }
  stop("프로젝트 루트(SPEC.md)를 찾지 못했습니다: getwd() = ", getwd())
})

proj_path <- function(...) file.path(PROJ_ROOT, ...)

source_project <- function() {
  files <- c("seeds.R", "params.R", "model.R", "truth.R", "population.R", "sampling.R",
             "nca.R", "be_stats.R", "simulate.R", "scenarios.R", "summarize.R",
             "trial.R", "mc.R", "step1.R", "figures.R", "postprocess.R", "logging.R")
  for (f in files) {
    p <- proj_path("R", f)
    if (file.exists(p)) sys.source(p, envir = globalenv())
  }
  invisible(TRUE)
}
