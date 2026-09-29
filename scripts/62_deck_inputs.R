#!/usr/bin/env Rscript
# 62_deck_inputs.R: 결과보고 슬라이드의 분포 그림에 쓸 작은 요약 파일을 만든다(results/deck_inputs/). 새 모의 없음.
# 입력은 로컬에 저장된 대상자 수준 결과(용량 때문에 git에 넣지 않은 .rds)와, 대표 대상자 곡선은 보고서 그림 2-4(scripts/34b)와 같은
# 대상자·같은 모델 풀이다. 각 요약은 커밋된 요약 결과와 대조해 같을 때만 쓴다(중앙값·백분위·비율, 허용 오차 1e-9).
# 산출(results/deck_inputs/):
#   rep_profiles.csv        2016 모델, 60–90 kg(base) LLOQ 도달일 25·50·75백분위 대표 대상자 3명의 참 농도 곡선(0.05일 간격, 0.01 mg/L 이상)
#   rep_subjects.csv        대표 대상자: 백분위, LLOQ 도달(투여 후 일), 절벽 시작(1일 정의), 도달 연구일
#   cliff_lloq_day_hist.csv 모델별(base) 참 농도가 LLOQ에 닿는 연구일의 1일 간격 분포(대상자 비율 %)
#   coverage_hist.csv       모델별(B0, 60–90 kg) 채혈 구간 커버리지와 관측 대 참 비(AUClast; NCA AUCinf 규칙 B; 규칙 A 세트 (i))의 분포(구간 비율 %)
#   strata_armdiff_hist.csv 동일 제품 시험(두 모델, 2,000회)의 arm 간 무거운 층 비율 차이(시험군 - 대조군, %p) 1%p 간격 분포: AUClast 분석군, 규칙 A 세트 (i)·(iii)
#   provenance.csv          입력 파일 SHA-256, 행 수, 대조 결과
source("R/00_setup.R"); source_project()
out_dir <- proj_path("results", "deck_inputs"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff; design <- read_cfg("trial_design.yaml"); N <- as.integer(cf$n_subjects); SEED <- as.integer(cf$seed)
sha <- function(f) digest::digest(file = f, algo = "sha256")
PROV <- list(); chk <- function(what, a, b, tol = 1e-9) { ok <- isTRUE(all(abs(a - b) <= tol)); PROV[[length(PROV) + 1L]] <<- data.table(check = what, deck_value = paste(signif(a, 12), collapse = ";"), committed_value = paste(signif(b, 12), collapse = ";"), ok = ok)
  if (!ok) stop("mismatch with committed summary: ", what) }

# ---- 절벽: LLOQ 도달일 분포와 대표 대상자 ----------------------------------------------------------------------------------------------------
cf_rds <- proj_path("results", "cliff", "cliff_subjects.rds"); cs_all <- readRDS(cf_rds); cs_sum <- fread(proj_path("results", "cliff", "cliff_summary.csv"))
H <- rbindlist(lapply(c("k2016", "k2020"), function(m) {
  x <- cs_all[model == m & weight == "base"]; s <- cs_sum[model == m & weight == "base"]
  chk(sprintf("%s base: n", m), nrow(x), s$n); d <- x[!is.na(t_lloq), t_lloq + 1]
  chk(sprintf("%s base: median study day at LLOQ", m), median(d), s$lloq_studyday_median)
  chk(sprintf("%s base: 5th/95th percentile study day at LLOQ", m), quantile(d, c(0.05, 0.95), names = FALSE), c(s$lloq_studyday_p05, s$lloq_studyday_p95))
  chk(sprintf("%s base: median cliff length (1-day definition)", m), median(x$len1, na.rm = TRUE), s$len1_median)
  b <- floor(d); data.table(pk_model = m, day_bin = b)[, .N, by = .(pk_model, day_bin)][, pct := 100 * N / length(d)][order(day_bin)]
}))
fwrite(H, file.path(out_dir, "cliff_lloq_day_hist.csv"))
base_cs <- cs_all[model == "k2016" & weight == "base" & !is.na(t_lloq)]
qs_ <- quantile(base_cs$t_lloq, c(0.25, 0.5, 0.75))
rep_ids <- vapply(qs_, function(q) base_cs$id[which.min(abs(base_cs$t_lloq - q))], numeric(1))   # scripts/34b와 같은 선택
p16 <- load_params("k2016")
subj16 <- with_seed(derive_seed(SEED, "k2016", "base", "subj"), make_subjects(N, p16, cliff_weight_spec(design, cf$weights$base), 0.5, 0))
prof <- solve_model(individual_params(p16, subj16[id %in% rep_ids]), CJ(id = rep_ids, time = seq(0.05, 70, by = 0.05)), design$dose_mg, p16$model_id)
prof <- as.data.table(prof)[, .(id, time, C)]
# 풀이 곡선의 LLOQ 교차 시각이 저장된 t_lloq와 같은지(0.05일 격자 해상도 안) 확인
LLOQ <- study_lloq()
for (i in rep_ids) { pr <- prof[id == i]; tc <- pr[C >= LLOQ, max(time)]; chk(sprintf("representative subject %d: LLOQ crossing within the 0.05-day grid", i), abs(tc - base_cs[id == i, t_lloq]) <= 0.05, TRUE, 0) }
fwrite(prof[C > 0.01], file.path(out_dir, "rep_profiles.csv"))
fwrite(data.table(id = rep_ids, percentile = c(25, 50, 75), t_lloq = base_cs[match(rep_ids, id), t_lloq], start1 = base_cs[match(rep_ids, id), start1],
                  study_day_lloq = base_cs[match(rep_ids, id), t_lloq] + 1, len1 = base_cs[match(rep_ids, id), len1]), file.path(out_dir, "rep_subjects.csv"))

# ---- 커버리지와 관측 대 참 비(시험 모집단 B0) --------------------------------------------------------------------------------------------------
ind_f <- c(k2016 = proj_path("results", "individual", "nca_base_20000.rds"), k2020 = proj_path("results", "individual", "nca_struct2020_20000.rds"))
tcv <- fread(proj_path("results", "trialpop", "tp_coverage_individual.csv"))
MET <- c(window = "window coverage (true AUC0-tlast / true AUC0-inf)", auclast = "observed-to-true, AUClast (all subjects)",
         aucinf_B = "observed-to-true, AUCinf rule B (lambda-z estimable)", aucinf_A = "observed-to-true, AUCinf rule A (meeting set (i))")
C <- rbindlist(lapply(names(ind_f), function(m) {
  x <- readRDS(ind_f[[m]])[schedule == "B0"]; lz <- x$lambda_ok %in% TRUE; ai <- crit_ok(x, "i")
  v <- list(window = x$coverage_true, auclast = x$AUClast / x$AUCinf_true, aucinf_B = x$AUCinf[lz] / x$AUCinf_true[lz], aucinf_A = x$AUCinf[ai] / x$AUCinf_true[ai])
  rbindlist(lapply(names(v), function(k) { y <- v[[k]][is.finite(v[[k]])]; ref <- tcv[pk_model == m & metric == MET[[k]]]
    chk(sprintf("%s %s: n, median, 5th, 95th, min, max", m, k), c(length(y), median(y), quantile(y, c(0.05, 0.95), names = FALSE), min(y), max(y)), c(ref$n, ref$median, ref$p05, ref$p95, ref$min, ref$max))
    br <- c(seq(0.5, 1.5, by = 0.02), Inf); lo <- pmax(y, 0.5)   # 0.5 미만은 첫 구간, 1.5 초과는 마지막 구간(개수로 보고)
    b <- cut(lo, breaks = c(0.5 - 1e-12, br[-1]), right = FALSE, labels = FALSE)
    data.table(pk_model = m, metric = k, bin_lo = br[b], n_bin = 1L)[, .(n = .N), by = .(pk_model, metric, bin_lo)][, pct := 100 * n / length(y)][order(bin_lo)]
  }))
}))
fwrite(C, file.path(out_dir, "coverage_hist.csv"))

# ---- 층 균형: 동일 제품 시험의 arm 간 무거운 층 비율 차이(시험군 - 대조군, %p) 분포(scripts/55와 같은 식) -------------------------------------
cnt_f <- c(k2016 = proj_path("results", "trialpop", "trialpop_counts_k2016.csv.gz"), k2020 = proj_path("results", "trialpop", "trialpop_counts_k2020.csv.gz"))
CNT <- rbindlist(lapply(names(cnt_f), function(m) fread(cnt_f[[m]])[, pk_model := m]))
tsc <- fread(proj_path("results", "trialpop", "tp_strata_composition.csv"))
comp <- function(d, col) { w <- dcast(d, pk_model + trial + scenario ~ stratum, value.var = col); setnames(w, c("1", "2"), c("a1", "a2")); w[, h := 100 * a2 / (a1 + a2)][, .(pk_model, trial, scenario, h)] }
REFC <- CNT[scenario == "REF"]; TC <- CNT[scenario == "S00"]
SX <- c(auclast = "n_auclast", i = "n_i", iii = "n_iii")
AD <- rbindlist(lapply(names(SX), function(k) {
  z <- comp(TC, SX[[k]])[comp(REFC, SX[[k]])[, .(pk_model, trial, h_R = h)], on = c("pk_model", "trial")][, armdiff := h - h_R]
  for (m in names(cnt_f)) { v <- z[pk_model == m, armdiff]; ref <- tsc[pk_model == m & scenario == "S00" & analysis_set == k]
    chk(sprintf("%s S00 %s: n, median, 5th, 95th, share above 5 points", m, k), c(length(v), median(v), quantile(v, c(0.05, 0.95), names = FALSE), 100 * mean(abs(v) > 5)),
        c(ref$n_trials, ref$armdiff_median, ref$armdiff_p05, ref$armdiff_p95, ref$armdiff_abs_gt5_pct)) }
  z[, .(pk_model, analysis_set = k, trial, armdiff)]
}))
fwrite(AD[, .(n = .N), by = .(pk_model, analysis_set, bin = floor(armdiff))][, pct := 100 * n / sum(n), by = .(pk_model, analysis_set)][order(pk_model, analysis_set, bin)],
       file.path(out_dir, "strata_armdiff_hist.csv"))

fwrite(rbind(data.table(input = c(cf_rds, ind_f), sha256 = vapply(c(cf_rds, ind_f), sha, ""), note = "local subject-level result (git-ignored; produced by scripts/34 and scripts/10)"),
             data.table(input = proj_path("results", "cliff", "cliff_summary.csv"), sha256 = sha(proj_path("results", "cliff", "cliff_summary.csv")), note = "committed summary used for the checks"),
             data.table(input = proj_path("results", "trialpop", "tp_coverage_individual.csv"), sha256 = sha(proj_path("results", "trialpop", "tp_coverage_individual.csv")), note = "committed summary used for the checks"),
             data.table(input = c(cnt_f, proj_path("results", "trialpop", "tp_strata_composition.csv")), sha256 = vapply(c(cnt_f, proj_path("results", "trialpop", "tp_strata_composition.csv")), sha, ""), note = "committed per-trial counts and summary")),
       file.path(out_dir, "provenance.csv"))
fwrite(rbindlist(PROV), file.path(out_dir, "checks.csv"))
cat(sprintf("deck inputs written to %s: %d checks against committed summaries, all equal\n", out_dir, length(PROV)))
