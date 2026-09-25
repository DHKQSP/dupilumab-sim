#!/usr/bin/env Rscript
# 전체 자동 테스트 실행. 결과를 logs/tests_<UTC>.log(요약), logs/tests_<UTC>_timing.csv(테스트별 시간)에 남긴다(DUPI_TEST_LOG_DIR로 변경 가능).
# 범위: 환경변수 DUPI_TEST_SCOPE = fast(CI) | full(기본). fast는 느린 테스트와 판정성 테스트를 건너뛴다(tests/testthat/helper-setup.R).
# 판정성 테스트(이름이 "[판정]"으로 시작): 기대값 불충족은 보고 항목으로 따로 적고 종료 코드에 넣지 않는다. 코드 오류(error)는 실패로 센다.
# JUnit XML(logs/tests_<UTC>.xml)은 xml2가 설치된 경우에만 쓴다(xml2는 testthat의 Suggests이며 renv.lock에 없다: CI 실패 원인, results/ci/).
# DUPI_TEST_STRICT_SKIP=true(CI): 예상하지 않은 건너뜀(fast 범위의 판정성 테스트, skip_if_fast_scope 외. 예: rxode2·NonCompart 미설치로
# skip_if_not_installed가 건너뜀)을 실패로 센다. 패키지 설치 문제로 핵심 테스트가 조용히 빠진 채 통과하는 것을 막는다.
source("R/00_setup.R"); source_project()
scope <- tolower(trimws(Sys.getenv("DUPI_TEST_SCOPE", "full"))); if (!nzchar(scope)) scope <- "full"
if (!scope %in% c("fast", "full")) stop("DUPI_TEST_SCOPE는 fast 또는 full이어야 합니다: ", scope)
Sys.setenv(DUPI_TEST_SCOPE = scope)
log_dir <- Sys.getenv("DUPI_TEST_LOG_DIR", proj_path("logs"))       # CI는 추적 파일과 섞이지 않는 별도 디렉터리를 준다(.github/workflows/tests.yml)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
log_dir <- normalizePath(log_dir, winslash = "/", mustWork = TRUE)  # 절대 경로: test_dir가 작업 디렉터리를 tests/testthat로 바꾼 동안 JunitReporter가 파일을 쓴다
strict_skip <- tolower(Sys.getenv("DUPI_TEST_STRICT_SKIP", "false")) %in% c("true", "1", "yes")
stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
logfile <- file.path(log_dir, sprintf("tests_%s.log", stamp))
has_xml2 <- requireNamespace("xml2", quietly = TRUE)
# max_failures = Inf: 기본값(10)이면 실패 기대값이 10개를 넘는 순간 그 파일 끝에서 나머지 파일을 실행하지 않고 끝난다
# (판정성 테스트의 불충족 10개가 뒤 파일의 코드 실패를 가려 종료 코드 0이 될 수 있다)
reporters <- list(testthat::ProgressReporter$new(show_praise = FALSE, max_failures = Inf))
if (has_xml2) reporters <- c(reporters, testthat::JunitReporter$new(file = sub("\\.log$", ".xml", logfile)))
t0 <- proc.time()
res <- testthat::test_dir(proj_path("tests", "testthat"), reporter = testthat::MultiReporter$new(reporters), stop_on_failure = FALSE)
el <- proc.time() - t0
df <- as.data.frame(res)
# 건너뛴 이유(첫 skip 메시지)
skip_msg <- vapply(res, function(t) {
  s <- Filter(function(e) inherits(e, "expectation_skip"), t$results)
  if (length(s)) sub("^Reason: ", "", conditionMessage(s[[1]])) else ""
}, "")
df$test[is.na(df$test)] <- "(code outside test_that)"                # 파일 최상위 코드의 오류·건너뜀
df$judgment <- startsWith(df$test, "[판정]")
# 예상한 건너뜀: fast 범위의 판정성 테스트(skip_if_judgment), 느린 테스트(skip_if_fast_scope)
df$unexpected_skip <- df$skipped & !((scope == "fast" & df$judgment) | startsWith(skip_msg, "fast 범위 제외"))
df$code_fail <- (df$failed > 0 & !df$judgment) | df$error | (strict_skip & df$unexpected_skip)
df$judgment_unmet <- df$failed > 0 & df$judgment & !df$error
timing <- data.frame(scope = scope, file = df$file, test = df$test, judgment = df$judgment, skipped = df$skipped, failed = df$failed > 0, error = df$error,
                     real_s = round(df$real, 3), cpu_s = round(df$user + df$system, 3), skip_reason = skip_msg, unexpected_skip = df$unexpected_skip)
fwrite(timing, sub("\\.log$", "_timing.csv", logfile))
summ <- data.frame(file = df$file, test = substr(df$test, 1, 70), failed = df$failed, skipped = df$skipped, error = df$error, warning = df$warning, real_s = round(df$real, 2))
jud <- df[df$judgment, ]; usk <- df[df$unexpected_skip, ]
lines <- c(sprintf("tests run at %s (UTC)", format(Sys.time(), tz = "UTC")),
           sprintf("scope: %s (DUPI_TEST_SCOPE)", scope),
           sprintf("rxode2: %s; testthat: %s; JUnit XML: %s", tryCatch(as.character(packageVersion("rxode2")), error = function(e) "NA"),
                   as.character(packageVersion("testthat")), if (has_xml2) "yes" else "no (xml2 not installed)"),
           capture.output(print(summ, row.names = FALSE)),
           if (nrow(jud)) c("judgment items (report only, not counted as suite failure):",
                            sprintf("  %-9s %s", ifelse(jud$skipped, "skipped", ifelse(jud$error, "ERROR", ifelse(jud$failed > 0, "unmet", "met"))), jud$test)),
           if (nrow(usk)) c(sprintf("unexpected skips (%s):", if (strict_skip) "counted as failures, DUPI_TEST_STRICT_SKIP" else "not counted; set DUPI_TEST_STRICT_SKIP=true to fail"),
                            sprintf("  %s | %s | %s", usk$file, usk$test, skip_msg[df$unexpected_skip])),
           sprintf("TOTAL: %d tests, %d failed (%d code, %d judgment unmet), %d errors, %d skipped (%d unexpected); wall %.1f s, cpu %.1f s",
                   nrow(df), sum(df$failed > 0), sum(df$failed > 0 & !df$judgment), sum(df$judgment_unmet), sum(df$error), sum(df$skipped), sum(df$unexpected_skip),
                   el[["elapsed"]], sum(el[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)))
writeLines(lines, logfile)
cat(lines, sep = "\n")
# GitHub Actions 작업 요약
gs <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(gs)) {
  bad <- df[df$code_fail | df$judgment_unmet, ]
  cat(c(sprintf("## Tests (scope: %s)", scope), "",
        sprintf("%d tests, %d code failures, %d errors, %d skipped (%d unexpected), %d judgment items unmet (report only). Wall %.0f s.",
                nrow(df), sum(df$failed > 0 & !df$judgment), sum(df$error), sum(df$skipped), sum(df$unexpected_skip), sum(df$judgment_unmet), el[["elapsed"]]), "",
        if (nrow(bad)) c("| file | test | status |", "|---|---|---|",
                         sprintf("| %s | %s | %s |", bad$file, gsub("|", "/", bad$test, fixed = TRUE), ifelse(bad$error | (bad$failed > 0 & !bad$judgment), "FAIL", ifelse(bad$code_fail, "unexpected skip", "judgment unmet")))), ""),
      file = gs, sep = "\n", append = TRUE)
}
if (any(df$code_fail)) quit(status = 1)
