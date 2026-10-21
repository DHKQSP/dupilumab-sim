#!/usr/bin/env Rscript
# 77_m15_amend.R — 개정 1(계획서 16절, D-069): 일치 점검, 추출, 진단, 교차, 부록, robustness. 등록 산출물은 읽기만 하고 results/m15/amend1/에 쓴다.
# 모드
#   mechcheck [model] [pop]   mechsuppx 일치 점검(16.3절 6): (i) Model 1 B0, 판 1~16 = 저장 mechsupp(모든 열), 묶음별 개수 행, versions.csv 1~16행;
#                             (ii) 두 모델, 세 일정, 항등 판(vi 1) = 저장 mech의 vi 1 행(모든 열). 있는 시험만 비교. 불일치면 종료 코드 3
# 사용법: nice -n 5 Rscript scripts/77_m15_amend.R <mode> ...
source("R/00_setup.R"); source_project(); source(proj_path("R", "m15.R"))
args <- commandArgs(trailingOnly = TRUE)
MODES <- c("mechcheck")
if (!length(args) || !args[1] %in% MODES) stop("사용법: Rscript scripts/77_m15_amend.R <", paste(MODES, collapse = "|"), "> ...")
mode <- args[1]
CFG <- m15_cfg(); MS <- as.integer(CFG$master_seed); stopifnot(MS == 20261020L)
AM <- CFG$amendment_20261021; stopifnot(!is.null(AM))
RD <- proj_path(CFG$output_dir); AD <- file.path(RD, "amend1"); dir.create(AD, showWarnings = FALSE, recursive = TRUE)
TD <- file.path(RD, "trials")

# 두 표가 열·행 모두 같은지(키로 정렬; 숫자는 상대 1e-12, 나머지는 동일). 반환: 최대 상대 차이(행·열 불일치면 Inf)
same_tab <- function(a, b, keys) {
  if (nrow(a) != nrow(b) || !setequal(names(a), names(b))) return(Inf)
  setcolorder(b, names(a)); setkeyv(a, keys); setkeyv(b, keys)
  if (!isTRUE(all.equal(a[, ..keys], b[, ..keys], check.attributes = FALSE))) return(Inf)
  max(vapply(setdiff(names(a), keys), function(cc) { x <- a[[cc]]; y <- b[[cc]]
    if (!identical(is.na(x), is.na(y))) return(Inf)
    z <- !is.na(x); if (!any(z)) return(0)
    if (is.numeric(x) && is.numeric(y)) max(abs(x[z] - y[z]) / pmax(abs(y[z]), 1e-300)) else if (identical(as.character(x[z]), as.character(y[z]))) 0 else Inf }, 1))
}
rd_files <- function(d, pat) { f <- sort(list.files(d, pattern = pat, full.names = TRUE)); if (!length(f)) return(NULL); setNames(lapply(f, fread), basename(f)) }

if (mode == "mechcheck") {
  X <- AM$mechsuppx; ms <- if (length(args) >= 2) args[2] else unlist(X$models); ps <- if (length(args) >= 3) args[3] else unlist(X$populations)
  lf <- start_run_log("m15_amend_mechcheck", master_seed = MS, run_mode = "final", extra = list(models = ms, populations = ps))
  R <- list()
  for (m in ms) for (pv in ps) {
    cell <- sprintf("%s_%s_%d", m, pv, as.integer(X$dose_mg)); dx <- file.path(TD, "mechsuppx", cell); if (!dir.exists(dx)) next
    bx <- rd_files(dx, "^mechsuppx_\\d+_\\d+\\.csv\\.gz$"); cx <- rd_files(dx, "^cnt_mechsuppx_\\d+_\\d+\\.csv\\.gz$"); if (is.null(bx)) next
    BX <- rbindlist(bx); trials <- sort(unique(BX$trial))
    # (ii) 항등 판: 세 일정, 저장 mech의 vi 1 행
    dm <- file.path(TD, "mech", cell); BM <- rbindlist(rd_files(dm, "^mech_\\d+_\\d+\\.csv\\.gz$"))
    a <- BX[vi == 1]; b <- BM[vi == 1 & trial %in% trials]
    R[[length(R) + 1]] <- data.table(cell = cell, check = "identity rows (vi 1) vs mech, all schedules", n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), copy(b), c("schedule", "trial")))
    # (i) Model 1 B0: 판 1~16 = 저장 mechsupp
    if (m == CFG$models$primary) {
      ds <- file.path(TD, "mechsupp", cell); BS <- rbindlist(rd_files(ds, "^mechsupp_\\d+_\\d+\\.csv\\.gz$"))
      a <- BX[schedule == "B0" & vi <= 16]; b <- BS[trial %in% trials]
      R[[length(R) + 1]] <- data.table(cell = cell, check = "B0 versions 1 to 16 vs mechsupp", n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), copy(b), c("schedule", "trial", "vi")))
      cs <- rd_files(ds, "^cnt_mechsupp_\\d+_\\d+\\.csv\\.gz$")
      for (nm in names(cx)) { s_ <- sub("^cnt_mechsuppx_", "cnt_mechsupp_", nm); if (is.null(cs[[s_]])) { R[[length(R) + 1]] <- data.table(cell = cell, check = paste("counts", nm), n_rows = NA_integer_, n_stored = 0L, max_rel_diff = Inf); next }
        a <- cx[[nm]][schedule == "B0" & vi <= 16]; b <- copy(cs[[s_]])
        R[[length(R) + 1]] <- data.table(cell = cell, check = paste("counts", nm), n_rows = nrow(a), n_stored = nrow(b), max_rel_diff = same_tab(copy(a), b, c("schedule", "vi", "n_rand", "model"))) }
      vx <- fread(file.path(dx, "versions.csv")); vs <- fread(file.path(ds, "versions.csv"))
      R[[length(R) + 1]] <- data.table(cell = cell, check = "versions.csv rows 1 to 16", n_rows = nrow(vx[vi <= 16]), n_stored = nrow(vs), max_rel_diff = same_tab(copy(vx[vi <= 16]), copy(vs), "vi"))
    }
  }
  CK <- rbindlist(R); stopifnot(nrow(CK) > 0); CK[, pass := is.finite(max_rel_diff) & max_rel_diff <= as.numeric(AM$identity_check$csv_relative) & n_rows == n_stored]
  fwrite(CK, file.path(AD, "mechsuppx_identity_check.csv")); print(CK)
  append_run_log(lf, sprintf("mechcheck: %d checks, all pass %s", nrow(CK), all(CK$pass)))
  if (!all(CK$pass)) quit(status = 3)
}
