#!/usr/bin/env Rscript
# §3-1·3-2 역산: 기전·방향별로 참 AUC0-inf 비가 목표값이 되는 시험군 배율 (config/oc_design.yaml, D-040).
# 사용법: Rscript scripts/30_oc_inversion.R <k2016|k2020> prep            — 공통 난수 200,000명과 대조 참값 저장
#         Rscript scripts/30_oc_inversion.R <k2016|k2020> <기전> [기전 ...] — 기전별 역산(F, ka, ke, Vmax, Km, V2)
# 절차: (1) 방향별 12점 로그 격자 선별(20,000명, AUC0-inf만) → 목표를 사이에 둔 첫 구간(배율 1에서 바깥쪽으로)
#       (2) 20,000명 로그 배율 이분법(허용 0.02%) (3) 200,000명으로 AUC0-inf·Cmax 참값 → 목표 ±0.1% 밖이면 200,000명으로 이분법 계속
#       도달 불가: 구간 끝 배율의 200,000명 참값(AUC0-inf, Cmax) 보고.
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); model <- args[1]; what <- args[-1]
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml")
p <- load_params(model); dose <- as.numeric(oc$estimand$dose_mg)
out_dir <- proj_path("results", "oc"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
N <- as.integer(oc$estimand$population$n_subjects); NS <- as.integer(oc$inversion$screening_subjects)
seed <- derive_seed(as.integer(oc$truth_seed), "truth200k", model)
subj <- truth_subjects(N, p, design, seed)                    # 결정적: 같은 시드 → 같은 200,000명
ref_f <- file.path(out_dir, sprintf("truth_ref_%s.rds", model))
logfile <- start_run_log(sprintf("oc_inversion_%s_%s", model, paste(what, collapse = "_")), master_seed = seed, run_mode = "final",
                         extra = list(model = model, what = what, n = N, n_screen = NS, atol = TRUTH_ATOL, rtol = TRUTH_RTOL))

if (identical(what, "prep")) {
  t0 <- Sys.time()
  ref <- truth_metrics(individual_params(p, subj), dose, p$model_id, cmax = TRUE)
  saveRDS(ref, ref_f)
  append_run_log(logfile, sprintf("prep done: %d subjects, %s", nrow(ref), format(Sys.time() - t0)))
  cat("prep done\n"); quit(save = "no")
}
ref <- readRDS(ref_f); ref_s <- ref[id <= NS]
subj_s <- subj[id <= NS]
targets <- as.numeric(unlist(oc$targets)); targets <- targets[abs(targets - 1) > 1e-12]
tol <- as.numeric(oc$inversion$tolerance_rel); tol_s <- as.numeric(oc$inversion$screening_tolerance_rel); maxit <- as.integer(oc$inversion$max_iter)

auc_ratio_s <- function(mech, m) truth_ratio(ref_s[, .(id, AUCinf)], truth_metrics(ip_for_mult(p, subj_s, mech_multiplier(mech, m)), dose, p$model_id, cmax = FALSE))$auc_ratio
full_eval <- function(mech, m, cmax = TRUE) { te <- truth_metrics(ip_for_mult(p, subj, mech_multiplier(mech, m)), dose, p$model_id, cmax = cmax)
  truth_ratio(if (cmax) ref else ref[, .(id, AUCinf)], te) }

res <- list(); scans <- list()
for (mech in what) {
  rg <- mech_range(oc, mech, p)
  for (dirn in c("down", "up")) {
    end_m <- if (dirn == "down") rg[1] else rg[2]
    grid_m <- exp(seq(0, log(end_m), length.out = 13))[-1]     # 배율 1에서 바깥쪽으로 12점
    t0 <- Sys.time()
    r_grid <- vapply(grid_m, function(m) auc_ratio_s(mech, m), numeric(1))
    scans[[length(scans) + 1]] <- data.table(model = model, mechanism = mech, direction = dirn, multiplier = c(1, grid_m), auc_ratio_screen = c(1, r_grid))
    cat(sprintf("[%s] %s %s scan %s: ratio %.4f .. %.4f\n", format(Sys.time(), "%H:%M:%S"), mech, dirn, format(Sys.time() - t0), min(r_grid), max(r_grid)))
    mm <- c(1, grid_m); rr <- c(1, r_grid)
    end_full <- NULL
    for (tg in targets) {
      k <- which((rr[-length(rr)] - tg) * (rr[-1] - tg) <= 0)       # 목표를 사이에 둔 구간(들)
      n_cross <- length(k)
      if (!n_cross) {
        if (is.null(end_full)) end_full <- full_eval(mech, end_m)
        res[[length(res) + 1]] <- data.table(model = model, mechanism = mech, direction = dirn, target = tg, reachable = FALSE, n_crossings = 0L,
          multiplier = NA_real_, auc_ratio = NA_real_, auc_se_log = NA_real_, cmax_ratio = NA_real_, cmax_se_log = NA_real_, iter_screen = NA_integer_, iter_full = NA_integer_,
          end_multiplier = end_m, end_auc_ratio = end_full$auc_ratio, end_cmax_ratio = end_full$cmax_ratio, n_subjects = N)
        next
      }
      k <- k[1]
      br <- c(mm[k], mm[k + 1]); rb <- c(rr[k], rr[k + 1])
      if (br[1] > br[2]) { br <- rev(br); rb <- rev(rb) }
      bs <- bisect_log(function(m) auc_ratio_s(mech, m), tg, br, rb, tol_s, maxit)
      fe <- full_eval(mech, bs$m); it_full <- 0L; m_fin <- bs$m
      if (abs(fe$auc_ratio / tg - 1) > tol) {                      # 200,000명으로 이분법 계속(선별 구간 안에서)
        bf <- bisect_log(function(m) full_eval(mech, m, cmax = FALSE)$auc_ratio, tg, br, rb, tol, maxit)
        m_fin <- bf$m; it_full <- bf$iter; fe <- full_eval(mech, m_fin)
      }
      res[[length(res) + 1]] <- data.table(model = model, mechanism = mech, direction = dirn, target = tg, reachable = TRUE, n_crossings = n_cross,
        multiplier = m_fin, auc_ratio = fe$auc_ratio, auc_se_log = fe$auc_se_log, cmax_ratio = fe$cmax_ratio, cmax_se_log = fe$cmax_se_log,
        iter_screen = bs$iter, iter_full = it_full, end_multiplier = end_m, end_auc_ratio = NA_real_, end_cmax_ratio = NA_real_, n_subjects = N)
      cat(sprintf("[%s]   target %.2f: m = %.5f, AUC %.5f (screen iters %d, full %d), Cmax %.4f\n", format(Sys.time(), "%H:%M:%S"), tg, m_fin, fe$auc_ratio, bs$iter, it_full, fe$cmax_ratio))
    }
  }
}
res <- rbindlist(res); scans <- rbindlist(scans)
res[, within_tol := reachable & abs(auc_ratio / target - 1) <= tol]
res[, overall_true_equiv := reachable & auc_ratio >= 0.80 & auc_ratio <= 1.25 & cmax_ratio >= 0.80 & cmax_ratio <= 1.25]
tag <- paste(what, collapse = "_")
fwrite(res, file.path(out_dir, sprintf("inversion_%s_%s.csv", model, tag))); fwrite(scans, file.path(out_dir, sprintf("inversion_scan_%s_%s.csv", model, tag)))
if (any(res$reachable & !res$within_tol)) warning("허용 오차 밖 역산 결과가 있습니다")
print(res[, .(mechanism, direction, target, reachable, multiplier, auc_ratio, cmax_ratio, end_auc_ratio)], digits = 4)
append_run_log(logfile, "done")
