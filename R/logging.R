# logging.R — 실행 로그: 시드, 세션 정보, 설정 해시, 시각
start_run_log <- function(run_name, master_seed, run_mode, extra = list()) {
  dir.create(proj_path("logs"), showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  logfile <- proj_path("logs", sprintf("%s_%s.log", run_name, ts))
  cfg_files <- list.files(proj_path("config"), full.names = TRUE)
  cfg_hash <- vapply(cfg_files, function(f) digest::digest(file = f, algo = "sha256"), "")
  names(cfg_hash) <- basename(cfg_files)
  git_sha <- tryCatch(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE), error = function(e) NA_character_)
  lines <- c(
    sprintf("run_name: %s", run_name),
    sprintf("started_utc: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
    sprintf("run_mode: %s", run_mode),
    sprintf("master_seed: %s", master_seed),
    sprintf("seed_rule_version: %s", SEED_RULE_VERSION),
    sprintf("git_sha: %s", if (length(git_sha)) git_sha else NA),
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
