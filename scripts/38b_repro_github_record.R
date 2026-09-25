#!/usr/bin/env Rscript
# GitHub Actions 재현성 실행 기록 (통합 지시 2026-09-25 §9, D-045·D-050)
# 입력: results/repro/repro_github_job_log_excerpt.txt(GitHub API 메타데이터 + 잡 로그 원문 줄), results/repro/repro_check.csv(로컬 실행)
# 출력: results/repro/repro_github_run.csv(field, value), results/repro/repro_github_items.csv(항목별 GitHub 재계산값·통과, 로컬 재계산값과의 일치)
# GitHub 쪽 재계산값은 로그에 찍힌 8자리 유효숫자이므로, 로컬 재계산값을 같은 형식(format(digits = 8))으로 찍은 문자열과 비교한다.
# 사용법: Rscript scripts/38b_repro_github_record.R
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "repro")
ex_f <- file.path(out_dir, "repro_github_job_log_excerpt.txt"); loc_f <- file.path(out_dir, "repro_check.csv")
if (!file.exists(ex_f) || !file.exists(loc_f)) stop("입력 파일이 없습니다: ", ex_f, ", ", loc_f)
L <- readLines(ex_f, encoding = "UTF-8")

# 메타데이터(key: value) 줄
kv <- grep("^[a-z_0-9]+: ", L, value = TRUE)
md <- setNames(sub("^[a-z_0-9]+: ", "", kv), sub(":.*$", "", kv))
need <- c("run_id", "run_number", "run_attempt", "workflow", "event", "head_sha", "head_branch", "conclusion", "html_url", "runner",
          "job_started_utc", "job_completed_utc", "step_tests", "step_repro", "renv_cache", "artifact")
if (length(miss <- setdiff(need, names(md)))) stop("발췌본에 메타데이터 없음: ", paste(miss, collapse = ", "))
ts <- function(x) as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
dur <- function(a, b) as.numeric(difftime(ts(b), ts(a), units = "secs"))
step_s <- function(k) { p <- strsplit(md[[k]], " ")[[1]]; dur(p[1], p[2]) }

# 점검 항목 줄: [시각] item committed X reproduced Y tol Z PASS|FAIL (scripts/38_repro_check.R의 출력 형식)
rx <- "^\\[[0-9:]+\\] (\\S+)\\s+committed (\\S+)\\s+reproduced (\\S+)\\s+tol (\\S+)\\s+(PASS|FAIL)$"
it <- grep(rx, L, value = TRUE)
gh <- data.table(item = sub(rx, "\\1", it), committed_printed = sub(rx, "\\2", it), reproduced_github_printed = sub(rx, "\\3", it),
                 tolerance_printed = sub(rx, "\\4", it), pass_github = sub(rx, "\\5", it) == "PASS")
loc <- fread(loc_f)
if (!setequal(gh$item, loc$item) || anyDuplicated(gh$item)) stop("GitHub 로그 항목과 로컬 항목이 다릅니다")
x <- merge(loc[, .(item, committed, tolerance, tolerance_type, reproduced_local = reproduced, pass_local = pass)], gh, by = "item", sort = FALSE)
f8 <- function(v) format(v, digits = 8)
x[, reproduced_local_printed := vapply(reproduced_local, f8, "")]
x[, committed_same := vapply(committed, f8, "") == committed_printed]
x[, reproduced_github := as.numeric(reproduced_github_printed)]
x[, identical_to_local_8sig := reproduced_local_printed == reproduced_github_printed]
if (!all(x$committed_same)) stop("GitHub 실행의 committed 값이 로컬 결과 파일과 다릅니다(참조 파일 불일치): ", paste(x[!committed_same, item], collapse = ", "))
x <- x[match(loc$item, item), .(item, committed, tolerance, tolerance_type, reproduced_github, reproduced_github_printed, pass_github,
                                 reproduced_local, reproduced_local_printed, pass_local, identical_to_local_8sig)]
fwrite(x, file.path(out_dir, "repro_github_items.csv"))

tests_line <- grep("^TOTAL: ", L, value = TRUE); ver_line <- grep("^rxode2: ", L, value = TRUE)
if (length(tests_line) != 1 || length(ver_line) != 1) stop("발췌본에 테스트 요약 줄 또는 버전 줄이 한 개가 아닙니다")
# 코어 수: 워크플로 파일의 실행 줄(Rscript scripts/38_repro_check.R <코어>)
wf <- grep("Rscript scripts/38_repro_check.R [0-9]+", readLines(proj_path(".github", "workflows", "repro.yml")), value = TRUE)
if (length(wf) != 1) stop("repro.yml에서 재현성 점검 실행 줄을 한 개 찾지 못했습니다")
cores_wf <- as.integer(sub("^.*38_repro_check.R ([0-9]+).*$", "\\1", wf))
r_ver <- regmatches(md[["renv_cache"]], regexpr("R version [0-9.]+ \\([0-9-]+\\)", md[["renv_cache"]]))
run <- list(workflow = md[["workflow"]], run_id = md[["run_id"]], run_number = md[["run_number"]], run_attempt = md[["run_attempt"]],
            event = md[["event"]], commit = md[["head_sha"]], branch = md[["head_branch"]], conclusion = md[["conclusion"]], url = md[["html_url"]],
            runner = md[["runner"]], r_version = r_ver, rxode2_version = sub("^rxode2: ([^;]+);.*$", "\\1", ver_line),
            testthat_version = sub("^.*testthat: ([^;]+);.*$", "\\1", ver_line), cores_used = cores_wf,
            job_started_utc = md[["job_started_utc"]], job_completed_utc = md[["job_completed_utc"]],
            job_duration_s = dur(md[["job_started_utc"]], md[["job_completed_utc"]]),
            step_tests_s = step_s("step_tests"), step_repro_check_s = step_s("step_repro"),
            tests = sub("^TOTAL: ", "", tests_line), renv_cache = if (grepl("^Cache hit", md[["renv_cache"]])) "hit (primary key)" else "miss",
            renv_status_note = if (any(grepl("out-of-sync", L))) "renv reported the project out-of-sync: patchwork (scripts/39_reliability_flags.R) and xml2 (optional JUnit output in scripts/run_tests.R) were used but not in renv.lock; neither is loaded by the tests or by scripts/38_repro_check.R. Fixed afterwards: patchwork 1.2.0 recorded in renv.lock, xml2 set as an ignored package (renv/settings.json)" else "in sync",
            artifact = md[["artifact"]], values_source = "job log printed values (8 significant digits); artifact download blocked by this session's egress policy",
            n_items = nrow(x), n_pass_github = sum(x$pass_github), all_pass_github = all(x$pass_github),
            n_identical_to_local_8sig = sum(x$identical_to_local_8sig), all_identical_to_local_8sig = all(x$identical_to_local_8sig))
fwrite(data.table(field = names(run), value = unname(vapply(run, function(v) paste(format(v), collapse = " "), ""))), file.path(out_dir, "repro_github_run.csv"))
print(x[, .(item, reproduced_github_printed, reproduced_local_printed, pass_github, identical_to_local_8sig)])
cat(sprintf("GitHub run %s (%s): %d/%d pass, %d/%d identical to local at 8 significant digits\n", run$run_number, substr(run$commit, 1, 7),
            run$n_pass_github, run$n_items, run$n_identical_to_local_8sig, run$n_items))
