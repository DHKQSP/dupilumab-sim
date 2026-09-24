# postprocess.R — 시험 수준 원자료(be_raw CSV, gzip 저장)에서 요약 표를 다시 만든다. 11/12 스크립트와 14_postprocess.R 공용.
read_raw <- function(path) {
  if (file.exists(path)) return(fread(path))
  if (!file.exists(paste0(path, ".gz"))) stop("원자료 없음: ", path)
  if (!requireNamespace("R.utils", quietly = TRUE)) stop("gzip 원자료를 읽으려면 R.utils 패키지가 필요합니다")   # fread의 .gz 읽기 의존성
  fread(paste0(path, ".gz"))
}

postprocess_schedules <- function(variant, design = read_cfg("trial_design.yaml"), verbose = TRUE) {
  out_dir <- proj_path("results", "trials")
  be <- read_raw(file.path(out_dir, sprintf("schedules_be_raw_%s.csv", variant)))
  st <- summarize_trials(be, method = "pooled_t"); st_anc <- summarize_trials(be, method = "ancova_weight")
  fwrite(st$per_endpoint, file.path(out_dir, sprintf("schedules_per_endpoint_%s.csv", variant)))
  fwrite(st$concordance, file.path(out_dir, sprintf("schedules_concordance_%s.csv", variant)))
  fwrite(st_anc$per_endpoint, file.path(out_dir, sprintf("schedules_per_endpoint_ancova_%s.csv", variant)))
  scens <- unique(be$scenario); scheds <- unique(be$schedule)
  pt <- rbindlist(lapply(scens, function(s_) paired_trial_width_vs_ref(be, "B0", s_)))
  fwrite(pt, file.path(out_dir, sprintf("paired_trial_vs_B0_%s.csv", variant)))
  dec <- NULL
  indf <- proj_path("results", "individual", sprintf("individual_%s.csv", variant))
  pindf <- proj_path("results", "individual", sprintf("paired_vs_B0_%s.csv", variant))
  if (file.exists(indf) && file.exists(pindf)) {
    ind <- fread(indf)[schedule %in% scheds]
    npts <- data.table(schedule = scheds, n_points = vapply(scheds, function(sh) length(get_schedule(design, sh)), numeric(1)))
    npts[, added_points := n_points - npts[schedule == "B0", n_points]]
    dec <- schedule_decision(ind, fread(pindf)[schedule %in% scheds], st$per_endpoint, pt, design, npts)
    fwrite(dec, file.path(out_dir, sprintf("schedule_decision_%s.csv", variant)))
  }
  if (verbose) cat(sprintf("postprocess schedules %s: %d trials, %d scenarios, %d schedules\n", variant, uniqueN(be$trial), length(scens), length(scheds)))
  invisible(list(per_endpoint = st$per_endpoint, concordance = st$concordance, paired = pt, decision = dec))
}

postprocess_products <- function(variant = "base", verbose = TRUE) {
  out_dir <- proj_path("results", "trials"); sfx <- if (variant == "base") "" else paste0("_", variant)
  be <- read_raw(file.path(out_dir, paste0("products_be_raw", sfx, ".csv")))
  st <- summarize_trials(be, method = "pooled_t"); st_anc <- summarize_trials(be, method = "ancova_weight")
  fwrite(st$per_endpoint, file.path(out_dir, paste0("products_per_endpoint", sfx, ".csv")))
  fwrite(st$concordance, file.path(out_dir, paste0("products_concordance", sfx, ".csv")))
  fwrite(st_anc$per_endpoint, file.path(out_dir, paste0("products_per_endpoint_ancova", sfx, ".csv")))
  if (verbose) cat(sprintf("postprocess products %s: %d trials, %d scenarios\n", variant, uniqueN(be$trial), uniqueN(be$scenario)))
  invisible(list(per_endpoint = st$per_endpoint, concordance = st$concordance))
}

# 개인 수준 대응 비교를 저장된 NCA 원자료(rds)에서 다시 만든다(시뮬레이션 재실행 없음). 부트스트랩 구간 포함.
postprocess_individual_paired <- function(variant, verbose = TRUE) {
  f <- proj_path("results", "individual", sprintf("nca_%s_20000.rds", variant))
  if (!file.exists(f)) { if (verbose) message("NCA 원자료 없음: ", f); return(invisible(NULL)) }
  nca <- readRDS(f)
  pi <- merge(paired_individual_vs_ref(nca, "B0"), paired_bootstrap_cd(nca, "B0"), by = "schedule")
  fwrite(pi, proj_path("results", "individual", sprintf("paired_vs_B0_%s.csv", variant)))
  if (verbose) cat(sprintf("postprocess individual paired %s: %d schedules\n", variant, nrow(pi)))
  invisible(pi)
}
