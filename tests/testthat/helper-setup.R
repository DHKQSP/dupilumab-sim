# helper-setup.R — testthat가 자동으로 먼저 읽는 파일. 프로젝트 루트를 찾아 R/ 모듈을 로드한다.
local({
  d <- normalizePath(getwd(), winslash = "/")
  for (i in 1:6) { if (file.exists(file.path(d, "SPEC.md"))) break; d <- dirname(d) }
  if (!file.exists(file.path(d, "SPEC.md"))) stop("프로젝트 루트(SPEC.md)를 찾지 못했습니다")
  source(file.path(d, "R", "00_setup.R"), local = FALSE)
})
source_project()
skip_if_no_rxode2 <- function() testthat::skip_if_not_installed("rxode2")
skip_if_no_bemaster <- function() if (!bemaster_available()) testthat::skip("BEmaster 미입수 (SPEC Q1) — NCA/BE 기반 검증은 보류")
# IIV·잔차가 출처 확정(confirmed) 상태일 때만 변동성 검증을 실행한다 (DEV-가정 값으로 검증 통과 판단 금지)
skip_if_variability_pending <- function(p) {
  st <- p$status[grepl("^(omega|sigma)\\.", item)]
  if (any(st$status != "confirmed")) testthat::skip("IIV/잔차 값이 출처 미확정 (SPEC Q2) — 변동성 검증 보류")
}
skip_if_ka_pending <- function(p) {
  st <- p$status[item == "theta.ka"]
  if (st$status != "confirmed") testthat::skip("ka 미보정 (SPEC Q3) — ka 확정값 기반 검증 보류")
}
