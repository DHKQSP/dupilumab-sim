# logging.R — 실행 로그: 시드, 세션 정보, 설정 해시, 시각
start_run_log <- function(run_name, master_seed, run_mode, extra = list()) {
  dir.create(proj_path("logs"), showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  logfile <- proj_path("logs", sprintf("%s_%s.log", run_name, ts))
  cfg_files <- list.files(proj_path("config"), full.names = TRUE)
  cfg_hash <- vapply(cfg_files, function(f) digest::digest(file = f, algo = "sha256"), "")
  names(cfg_hash) <- basename(cfg_files)
  git_sha <- tryCatch(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE), error = function(e) NA_character_)
  # 코드 출처: 추적 파일의 커밋 전 변경(R/, scripts/, config/, tests/)이 있으면 커밋 해시만으로 실행 코드를 특정할 수 없으므로 목록을 남긴다(D-054)
  dirty <- tryCatch(suppressWarnings(system("git status --porcelain --untracked-files=no -- R scripts config tests", intern = TRUE, ignore.stderr = TRUE)), error = function(e) NA_character_)
  lines <- c(
    sprintf("run_name: %s", run_name),
    sprintf("started_utc: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
    sprintf("run_mode: %s", run_mode),
    sprintf("master_seed: %s", master_seed),
    sprintf("seed_rule_version: %s", SEED_RULE_VERSION),
    sprintf("git_sha: %s", if (length(git_sha)) git_sha else NA),
    sprintf("git_code_clean: %s", if (length(dirty) == 1 && is.na(dirty)) NA else length(dirty) == 0),
    if (length(dirty) && !all(is.na(dirty))) sprintf("  uncommitted: %s", trimws(dirty)) else NULL,
    sprintf("R_version: %s", R.version.string),
    sprintf("rxode2: %s", tryCatch(as.character(packageVersion("rxode2")), error = function(e) "not installed")),
    sprintf("BEmaster: %s", tryCatch(as.character(packageVersion("BEmaster")), error = function(e) "not installed")),
    "config_sha256:",
    sprintf("  %s: %s", names(cfg_hash), cfg_hash),
    if (length(extra)) c("extra:", sprintf("  %s: %s", names(extra), vapply(extra, function(x) paste(format(x), collapse = ","), ""))) else NULL,
    "session_info:",
    paste0("  ", capture.output(print(sessionInfo())))
  )
  writeLines(lines, logfile)
  message("run log: ", logfile)
  invisible(logfile)
}

append_run_log <- function(logfile, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S", tz = "UTC"), paste0(..., collapse = "")),
      file = logfile, append = TRUE)
}
