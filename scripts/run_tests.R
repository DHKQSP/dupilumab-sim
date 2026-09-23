#!/usr/bin/env Rscript
# 전체 자동 테스트 실행. 결과를 logs/tests_<UTC>.log 에 남긴다.
source("R/00_setup.R"); source_project()
dir.create(proj_path("logs"), showWarnings = FALSE)
logfile <- proj_path("logs", sprintf("tests_%s.log", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")))
res <- testthat::test_dir(proj_path("tests", "testthat"), reporter = testthat::MultiReporter$new(list(
  testthat::ProgressReporter$new(), testthat::JunitReporter$new(file = sub("\\.log$", ".xml", logfile)))),
  stop_on_failure = FALSE)
df <- as.data.frame(res)
summ <- data.frame(file = df$file, test = df$test, failed = df$failed, skipped = df$skipped, error = df$error, warning = df$warning)
writeLines(c(sprintf("tests run at %s (UTC)", format(Sys.time(), tz = "UTC")),
             sprintf("rxode2: %s", tryCatch(as.character(packageVersion("rxode2")), error = function(e) "NA")),
             capture.output(print(summ, row.names = FALSE)),
             sprintf("TOTAL: %d tests, %d failed, %d errors, %d skipped", nrow(df), sum(df$failed > 0), sum(df$error), sum(df$skipped))), logfile)
cat(readLines(logfile), sep = "\n")
if (any(df$failed > 0) || any(df$error)) quit(status = 1)
