#!/usr/bin/env Rscript
# §5-3 판정 불일치 분류·소비자 위험, §5-4 fallback 비용, §5-5 AUCinf 처리 규칙 비교, §5-6 표본 수 교차점검.
# 입력: results/trials5000/products5000_be_raw_base.csv.gz(주 모델 5,000회 이상), results/individual/nca_*_20000.rds, results/trials/products_per_endpoint_struct2020.csv 등.
source("R/00_setup.R"); source_project()
design <- read_cfg("trial_design.yaml"); lim <- as.numeric(unlist(design$be$limits))
out_dir <- proj_path("results", "fallback"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
be <- read_raw(proj_path("results", "trials5000", "products5000_be_raw_base.csv"))
st <- summarize_trials_ci(be); w <- st$wide
reps <- be[, .(n_trials = uniqueN(trial)), by = scenario]
truth <- st$per_endpoint[endpoint == "AUCinf_true", .(scenario, true_ratio = GMR_mean)]
truth[, true_inside := true_ratio >= lim[1] & true_ratio <= lim[2]]

# 5-3 불일치 분류
dis <- merge(st$concordance[, .(scenario, n_trials, last_pass_inf_fail, last_pass_inf_fail_lo, last_pass_inf_fail_hi)], truth, by = "scenario")
dis[, classification := fifelse(true_inside, "AUCinf 위음성(참값 비가 80–125% 안)", "참값 비 범위 밖: 소비자 위험 표 참조")]
setorder(dis, -last_pass_inf_fail)
fwrite(dis, file.path(out_dir, "discordance_classification.csv"))
cr <- rbindlist(lapply(truth[true_inside == FALSE, scenario], function(s_) { x <- w[scenario == s_]
  j2 <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax), nrow(x)); j3 <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_reliable), nrow(x))
  d <- paired_prop_diff_ci(x$pass_AUClast & x$pass_Cmax, x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_reliable)
  data.table(scenario = s_, true_ratio = truth[scenario == s_, true_ratio], n_trials = nrow(x), joint_last_cmax = j2$est, joint_last_cmax_lo = j2$lo, joint_last_cmax_hi = j2$hi,
             joint_3 = j3$est, joint_3_lo = j3$lo, joint_3_hi = j3$hi, only_inf_fails = d$est, only_inf_fails_lo = d$lo, only_inf_fails_hi = d$hi) }))
fwrite(cr, file.path(out_dir, "consumer_risk.csv"))

# 5-4 fallback 비용: 동시 통과율 (i) AUClast+Cmax, (ii) +AUCinf 신뢰군, (iii) +AUCinf 산출 가능 전체
fb <- rbindlist(lapply(c("S00", "KE110", "F097"), function(s_) { x <- w[scenario == s_]
  i <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax), nrow(x)); ii <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_reliable), nrow(x))
  iii <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_all), nrow(x))
  d2 <- paired_prop_diff_ci(x$pass_AUClast & x$pass_Cmax, x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_reliable)
  d3 <- paired_prop_diff_ci(x$pass_AUClast & x$pass_Cmax, x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_all)
  data.table(scenario = s_, true_ratio = truth[scenario == s_, true_ratio], n_trials = nrow(x), i_last_cmax = i$est, i_lo = i$lo, i_hi = i$hi,
             ii_plus_inf_rel = ii$est, ii_lo = ii$lo, ii_hi = ii$hi, iii_plus_inf_all = iii$est, iii_lo = iii$lo, iii_hi = iii$hi,
             cost_ii_pp = d2$est, cost_ii_lo = d2$lo, cost_ii_hi = d2$hi, cost_iii_pp = d3$est, cost_iii_lo = d3$lo, cost_iii_hi = d3$hi) }))
fwrite(fb, file.path(out_dir, "fallback_cost.csv"))

# 5-5 AUCinf 처리 규칙: (A) 신뢰군, (B) 산출 가능 전체, (C) 미충족은 AUClast 대입
rules <- c(A = "AUCinf_reliable", B = "AUCinf_all", C = "AUCinf_subC", `Reference: AUClast` = "AUClast")
rr <- rbindlist(lapply(names(rules), function(r) { ep <- rules[[r]]
  n <- be[endpoint == ep & scenario == "S00", .(n_R = mean(n_R), n_T = mean(n_T))]
  s00 <- be[endpoint == ep & scenario == "S00"]; pw <- wilson_ci(sum(s00$pass), nrow(s00))
  bias <- rbindlist(lapply(c("VM125", "F090"), function(s_) { g <- mc_mean_ci(be[endpoint == ep & scenario == s_, GMR]); tr <- truth[scenario == s_, true_ratio]
    data.table(scenario = s_, GMR = g$est, GMR_lo = g$lo, GMR_hi = g$hi, true_ratio = tr, bias_pct = 100 * (g$est / tr - 1)) }))
  agree <- rbindlist(lapply(c("S00", "VM125", "F090"), function(s_) { x <- w[scenario == s_]; v <- x[[paste0("pass_", ep)]] == x$pass_AUClast; a <- wilson_ci(sum(v, na.rm = TRUE), sum(!is.na(v)))
    data.table(scenario = s_, agree_with_AUClast = a$est, agree_lo = a$lo, agree_hi = a$hi) }))
  data.table(rule = r, endpoint = ep, n_R_mean = n$n_R, n_T_mean = n$n_T, pass_S00 = pw$est, pass_S00_lo = pw$lo, pass_S00_hi = pw$hi,
             bias_VM125_pct = bias[scenario == "VM125", bias_pct], GMR_VM125 = bias[scenario == "VM125", GMR], bias_F090_pct = bias[scenario == "F090", bias_pct], GMR_F090 = bias[scenario == "F090", GMR],
             agree_S00 = agree[scenario == "S00", agree_with_AUClast], agree_VM125 = agree[scenario == "VM125", agree_with_AUClast], agree_VM125_lo = agree[scenario == "VM125", agree_lo], agree_VM125_hi = agree[scenario == "VM125", agree_hi],
             agree_F090 = agree[scenario == "F090", agree_with_AUClast]) }))
fwrite(rr, file.path(out_dir, "aucinf_rules.csv"))

# 5-6 표본 수 교차점검: log-scale SD·CV (20,000명, B0), 경험적 검정력(arm당 117명)
sd_tab <- rbindlist(lapply(c(k2016 = "base", k2020 = "struct2020"), function(v) { x <- readRDS(proj_path("results", "individual", sprintf("nca_%s_20000.rds", v)))[schedule == "B0"]
  data.table(variant = v, AUClast_logsd = sd(log(x$AUClast), na.rm = TRUE), AUClast_logcv = log_cv_pct(x$AUClast), Cmax_logsd = sd(log(x$Cmax), na.rm = TRUE), Cmax_logcv = log_cv_pct(x$Cmax)) }))
sd_tab <- rbind(sd_tab, data.table(variant = "Cohen 2022 역산(AUClast 90% CI 0.96–1.28, n 62/63)", AUClast_logsd = 0.49, AUClast_logcv = 100 * sqrt(exp(0.49^2) - 1)), fill = TRUE)
sd_tab <- rbind(sd_tab, data.table(variant = "Syneos 제안", AUClast_logcv = 43), fill = TRUE)
fwrite(sd_tab, file.path(out_dir, "sample_size_logsd.csv"))
pw_row <- function(wide, s_, model, nrep) { x <- wide[scenario == s_]; if (!nrow(x)) return(NULL); a <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax), nrow(x)); b <- wilson_ci(sum(x$pass_AUClast & x$pass_Cmax & x$pass_AUCinf_reliable), nrow(x))
  data.table(model = model, scenario = s_, n_trials = nrow(x), power_last_cmax = a$est, lo = a$lo, hi = a$hi, power_3 = b$est, lo3 = b$lo, hi3 = b$hi) }
pw <- rbind(pw_row(w, "S00", "k2016", 0), pw_row(w, "F097", "k2016", 0))
b20 <- proj_path("results", "trials", "products_be_raw_struct2020.csv")
if (file.exists(b20) || file.exists(paste0(b20, ".gz"))) { w20 <- summarize_trials_ci(read_raw(b20))$wide; pw <- rbind(pw, pw_row(w20, "S00", "k2020", 0), pw_row(w20, "F097", "k2020", 0)) }
fwrite(pw, file.path(out_dir, "empirical_power.csv"))
print(dis); print(cr); print(fb); print(rr); print(sd_tab); print(pw)
