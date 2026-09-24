#!/usr/bin/env Rscript
# §3-3 무작위 다변량 제품 공간 (config/oc_design.yaml random_space, D-040).
# 제품 20,000개: 라틴 하이퍼큐브, 기전별 로그 균등(F, ka, ke, Vmax, Km, V2). 두 모델에 같은 제품 집합.
# stage truth : 제품별 참 AUC0-inf 비·참 Cmax 비 — 모든 제품에 공통인 1,000명(공통 난수), 짝지은 추정과 MC 표준오차
# stage trials: 제품 i = 시험 i(대조군도 시험마다 새로 추출), B0, arm당 117명, pooled t, 평가변수 6종
# 사용법: Rscript scripts/32_oc_random_space.R <k2016|k2020> <truth|trials> [cores]
source("R/00_setup.R"); source_project()
args <- commandArgs(trailingOnly = TRUE); model <- args[1]; stage <- args[2]; cores <- if (length(args) >= 3) as.integer(args[3]) else 2L
oc <- read_cfg("oc_design.yaml"); design <- read_cfg("trial_design.yaml"); rs <- oc$random_space
out_dir <- proj_path("results", "oc"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
variant <- if (model == "k2016") "base" else "struct2020"
rv <- resolve_variant(variant, design); p <- rv$p; dose <- as.numeric(oc$estimand$dose_mg)
n_prod <- as.integer(rs$n_products); dims <- names(rs$ranges)
if (length(args) >= 4) { n_prod <- as.integer(args[4]); out_dir <- proj_path("results", "oc_test"); dir.create(out_dir, showWarnings = FALSE) }   # 시험 실행용

# 제품 집합(모델과 무관, 결정적)
prod_f <- file.path(out_dir, "random_products.csv")
lhs <- with_seed(as.integer(rs$lhs_seed), {
  m <- sapply(dims, function(d) { r <- as.numeric(unlist(rs$ranges[[d]])); u <- (sample.int(n_prod) - runif(n_prod)) / n_prod; exp(log(r[1]) + u * (log(r[2]) - log(r[1]))) })
  data.table(product = seq_len(n_prod), m)
})
if (!file.exists(prod_f)) fwrite(lhs, prod_f) else stopifnot(isTRUE(all.equal(fread(prod_f), lhs, tolerance = 1e-12, check.attributes = FALSE)))
mults <- function(i) as.list(unlist(lhs[i, ..dims]))
logfile <- start_run_log(sprintf("oc_random_%s_%s", model, stage), master_seed = as.integer(rs$trial_master_seed), run_mode = "final",
                         extra = list(model = model, stage = stage, n_products = n_prod, truth_subjects = rs$truth_subjects, cores = cores))

if (stage == "truth") {
  f <- file.path(out_dir, sprintf("random_truth_%s.csv", model))
  NT <- as.integer(rs$truth_subjects)
  subj <- truth_subjects(NT, p, design, derive_seed(as.integer(rs$truth_seed), "rs_truth", model))
  ref <- truth_metrics(individual_params(p, subj), dose, p$model_id)
  done <- if (file.exists(f)) fread(f, select = "product")$product else integer(0)
  block <- 20L
  one_block <- function(ids) {
    ip <- rbindlist(lapply(seq_along(ids), function(k) individual_params(apply_multipliers(p, mults(ids[k])), subj)[, id := id + (k - 1L) * NT]))
    tm <- truth_metrics(ip, dose, p$model_id, chunk = 25000L)
    rbindlist(lapply(seq_along(ids), function(k) {
      te <- tm[id > (k - 1L) * NT & id <= k * NT][, id := id - (k - 1L) * NT]
      r <- truth_ratio(ref, te)
      data.table(product = ids[k], true_auc_ratio = r$auc_ratio, auc_se_log = r$auc_se_log, true_cmax_ratio = r$cmax_ratio, cmax_se_log = r$cmax_se_log, n_subjects = r$n)
    }))
  }
  for (b0 in seq(1L, n_prod, by = 1000L)) {
    ids <- setdiff(b0:min(b0 + 999L, n_prod), done); if (!length(ids)) next
    t0 <- Sys.time()
    blocks <- split(ids, ceiling(seq_along(ids) / block))
    res <- parallel::mclapply(blocks, one_block, mc.cores = cores, mc.preschedule = FALSE)
    bad <- vapply(res, function(x) inherits(x, "try-error"), logical(1)); if (any(bad)) stop("참값 계산 실패: ", as.character(res[bad][[1]]))
    fwrite(rbindlist(res), f, append = file.exists(f))
    msg <- sprintf("truth products %d–%d done in %s", min(ids), max(ids), format(Sys.time() - t0)); cat(msg, "\n"); append_run_log(logfile, msg)
  }
} else if (stage == "trials") {
  f <- file.path(out_dir, sprintf("random_trials_be_%s.csv.gz", model))
  done <- if (file.exists(f)) unique(fread(f, select = "trial")$trial) else integer(0)
  MS <- as.integer(rs$trial_master_seed)
  invisible(get_model(p$model_id))
  one <- function(j) { sc <- list(P = list(code = "P", T_multipliers = mults(j)))
    run_trial_oc(j, p, design, sc, MS, rv$wt_spec, model_id = p$model_id)$be }
  for (b0 in seq(1L, n_prod, by = 1000L)) {
    ids <- setdiff(b0:min(b0 + 999L, n_prod), done); if (!length(ids)) next
    t0 <- Sys.time()
    res <- parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE)
    bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1)); if (any(bad)) stop("시험 실패: ", paste(ids[bad], collapse = ","))
    be <- rbindlist(res)[, `:=`(GMR = signif(GMR, 7), CI_lower = signif(CI_lower, 7), CI_upper = signif(CI_upper, 7))]
    fwrite(be, f, append = file.exists(f))
    msg <- sprintf("trials products %d–%d done in %s", min(ids), max(ids), format(Sys.time() - t0)); cat(msg, "\n"); append_run_log(logfile, msg)
  }
}
append_run_log(logfile, "done")
