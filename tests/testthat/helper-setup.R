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

# 테스트 범위 (README "CI"): 환경변수 DUPI_TEST_SCOPE = "fast"(CI) | "full"(로컬 기본, 미설정 시)
# fast: 코드 정확성 테스트만(모델 구조, NCA 엔진, config 스키마, 판정 규칙, BE 통계, 소규모 스모크 모의, OC 구성 요소).
# full: fast + 느린 테스트(로컬 약 30초 초과) + 판정성 테스트.
test_scope <- function() {
  s <- tolower(trimws(Sys.getenv("DUPI_TEST_SCOPE", "full")))
  if (!nzchar(s)) s <- "full"
  if (!s %in% c("fast", "full")) stop("DUPI_TEST_SCOPE는 fast 또는 full이어야 합니다: ", s)
  s
}
# 느린 테스트: fast 범위에서 건너뛴다(통과 여부가 모의 결과에 달려 있지 않은 코드 테스트)
skip_if_fast_scope <- function(reason) if (identical(test_scope(), "fast")) testthat::skip(paste0("fast 범위 제외(느린 테스트): ", reason))
# 판정성 테스트: 통과 여부가 모의 결과(과학적 판정)에 달린 항목. CI(fast)에서는 건너뛰고, full에서는 실행하되
# 테스트 이름을 "[판정]"으로 시작해 scripts/run_tests.R가 기대값 불충족을 스위트 실패가 아닌 보고 항목으로 분류하게 한다
# (코드 오류(error)는 판정성 테스트라도 실패로 센다). 정식 판정 기록은 scripts/02_validate_step1.R 결과(results/step1/).
JUDGMENT_PREFIX <- "[판정]"
skip_if_judgment <- function(reason) if (identical(test_scope(), "fast")) testthat::skip(reason)
