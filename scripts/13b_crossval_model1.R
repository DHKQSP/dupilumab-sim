#!/usr/bin/env Rscript
# Model 1(자체 IIV·잔차) 교차검증 (검토 의견 통합본 §0-2). 검토자 독립 구현은 transit을 1차 흡수로 근사 → ±3% 이내면 일치.
# 20,000명, Syneos 일정(B0), 채혈 편차 없음. 비교용 2016 주 모델 60–90 kg 포함.
source("R/00_setup.R"); source_project()
design <- read_cfg("trial_design.yaml"); MASTER_SEED <- 20260923L
out_dir <- proj_path("results", "crossval"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- start_run_log("crossval_model1", master_seed = MASTER_SEED, run_mode = "final")
w6090 <- weight_spec_from_design(design, "base"); w78 <- list(mean = 78, sd = 11, trunc = c(55, 110))
m1 <- load_params("k2020"); m16 <- load_params("k2016")
run <- function(p, wt, tag) { x <- run_individual_population(20000, p, design, "B0", MASTER_SEED, wt, jitter = FALSE, model_id = p$model_id, tag = tag)$nca
  data.table(AUClast_mean = mean(x$AUClast, na.rm = TRUE), AUClast_geo = geo_mean(x$AUClast), AUClast_logcv = log_cv_pct(x$AUClast),
             Cmax_mean = mean(x$Cmax, na.rm = TRUE), Cmax_logcv = log_cv_pct(x$Cmax), tlast_median = median(x$tlast, na.rm = TRUE), tlast_p95 = q95(x$tlast),
             extrap_true_median = median(x$pct_extrap_true, na.rm = TRUE), extrap_true_p95 = q95(x$pct_extrap_true), extrap_nca_median = median(x$pct_extrap, na.rm = TRUE),
             extrap_gt20_pct = 100 * mean(x$pct_extrap > 20, na.rm = TRUE),
             reliable_pct = 100 * mean(x$lambda_ok & x$adj_r2 >= 0.80 & x$pct_extrap <= 20, na.rm = TRUE),   # 검토자 기준값과 같은 정의(adj R²·외삽만)
             reliable_new_def_pct = 100 * mean(x$reliable)) }   # 새 엔진 정의(span ratio 포함, D-039) — 기준값 없음
a <- run(m1, w6090, "xv_m1_6090"); b <- run(m1, w78, "xv_m1_78"); c16 <- run(m16, w6090, "xv_k16_6090")
long <- function(x, cond) melt(x[, condition := cond], id.vars = "condition", variable.name = "metric", value.name = "sim", variable.factor = FALSE)
sim <- rbind(long(a, "Model 1, 60–90 kg"), long(b, "Model 1, 약 78 kg (N(78,11) 55–110)"), long(c16, "2016 주 모델, 60–90 kg"))
ref <- rbind(data.table(condition = "Model 1, 60–90 kg", metric = c("AUClast_mean", "AUClast_geo", "AUClast_logcv", "Cmax_mean", "Cmax_logcv", "tlast_median", "tlast_p95", "extrap_true_median", "extrap_true_p95", "extrap_nca_median", "extrap_gt20_pct", "reliable_pct"),
                        ref = c(637, 587, 43.2, 34.9, 34.1, 35, 56, 0.61, 3.29, 2.79, 0.62, 91.4)),
             data.table(condition = "Model 1, 약 78 kg (N(78,11) 55–110)", metric = c("AUClast_mean", "Cmax_mean"), ref = c(610, 33.7)),
             data.table(condition = "2016 주 모델, 60–90 kg", metric = c("AUClast_mean", "AUClast_geo", "AUClast_logcv", "Cmax_mean", "reliable_pct"), ref = c(587, 546, 39.9, 38.2, 86.2)))
cmp <- merge(ref, sim, by = c("condition", "metric"))
cmp[, rel_diff_pct := 100 * (sim / ref - 1)]
cmp[, agree_3pct := abs(rel_diff_pct) <= 3]
cmp[metric %in% c("extrap_true_median", "extrap_true_p95", "extrap_nca_median", "extrap_gt20_pct"), note := "작은 백분율 지표: 상대 차이는 참고"]
cmp[metric %in% c("tlast_median", "tlast_p95"), note := "경과일(연구일 − 1)"]
cmp <- rbind(cmp, data.table(condition = "Model 1, 약 78 kg (N(78,11) 55–110)", metric = c("AUClast_ratio_vs_obs544", "Cmax_ratio_vs_obs32"), ref = c(610 / 544, 33.7 / 32), sim = c(b$AUClast_mean / 544, b$Cmax_mean / 32)), fill = TRUE)
fwrite(cmp, file.path(out_dir, "crossval_model1.csv")); print(cmp, digits = 4)
append_run_log(logfile, "done")
