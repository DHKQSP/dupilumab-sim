#!/usr/bin/env Rscript
# results/key_numbers_en.md: key numbers in English with their source files (review 2026-09-24 section 6).
# Rules as summary_en.md: English only, no em-dash; proportions with 95% intervals where they come from trials.
source("R/00_setup.R"); source_project()
R <- function(...) { f <- proj_path("results", ...); if (file.exists(f) || file.exists(paste0(f, ".gz"))) read_raw(f) else NULL }
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2); f3 <- function(x) formatC(x, format = "f", digits = 3)
ci <- function(e, lo, hi, d = 1) sprintf("%s%% (95%% CI %s to %s)", formatC(e, format = "f", digits = d), formatC(lo, format = "f", digits = d), formatC(hi, format = "f", digits = d))
out <- c("# Key numbers (dupilumab biosimilar Phase 1 PK simulation)", "", sprintf("Generated %s. Each line names its source file under results/.", format(Sys.Date())), "")
add <- function(section, ...) out <<- c(out, if (!is.null(section)) c(paste0("## ", section), ""), ...)
ML <- c(k2016 = "2016 model", k2020 = "Model 1", base = "2016 model", struct2020 = "Model 1")

ev <- R("nca_engine", "engine_validation_summary.csv"); dr <- R("nca_engine", "dropout_reasons_B0.csv"); ed <- R("nca_engine", "engine_difference_individual_B0.csv")
if (!is.null(ev)) add("NCA engine (Phoenix WinNonlin-compatible)",
  sprintf("- Own engine vs NonCompart 0.8.4 vs PKNCA 0.12.1: lambda-z points identical in %s of %s profile comparisons; largest relative parameter difference %s (nca_engine/engine_validation_summary.csv).",
          format(sum(ev$lz_points_identical), big.mark = ","), format(sum(ev$n_profiles), big.mark = ","), formatC(max(ev$max_rel_diff), format = "e", digits = 1)),
  if (!is.null(dr)) sprintf("- AUCinf reliability failure at B0 (any flag or lambda-z not estimable): %s%% (2016 model) and %s%% (Model 1); span ratio below 2 flags %s%% and %s%% (nca_engine/dropout_reasons_B0.csv).",
          f1(dr[model == "base", any_flag_or_fail_pct]), f1(dr[model == "struct2020", any_flag_or_fail_pct]), f1(dr[model == "base", flag_span_pct]), f1(dr[model == "struct2020", flag_span_pct])),
  if (!is.null(ed)) sprintf("- Reliability rate, previous vs new engine: %s%% to %s%% (2016 model), %s%% to %s%% (Model 1) (nca_engine/engine_difference_individual_B0.csv).",
          f1(ed[model == "base"][1, reliable_pct]), f1(ed[model == "base"][2, reliable_pct]), f1(ed[model == "struct2020"][1, reliable_pct]), f1(ed[model == "struct2020"][2, reliable_pct])), "")

cs <- R("cliff", "cliff_summary.csv"); cp <- R("cliff", "cliff_points.csv")
if (!is.null(cs)) { b <- cs[model == "k2016" & weight == "base"]; n1 <- cp[model == "k2016" & weight == "base" & timing == "nominal"]
  add("Sampling cliff (2016 model, 60 to 90 kg, 20,000 subjects)",
      sprintf("- True LLOQ reached at study Day %s (5th to 95th percentile %s to %s); after Day 58 in %s%% (cliff/cliff_summary.csv).", f1(b$lloq_studyday_median), f1(b$lloq_studyday_p05), f1(b$lloq_studyday_p95), f1(b$lloq_after_day58_pct)),
      sprintf("- Cliff length: %s days (1-day definition), %s days (2-day definition); cliff starts at %s mg/L (median).", f2(b$len1_median), f2(b$len2_median), f2(b$c_start1_median)),
      sprintf("- Two or more samples in the cliff at nominal days, fixed schedules, 1-day definition: at most %s%%; three or more with daily Day 29 to 57 sampling: %s%% (cliff/cliff_points.csv).",
              f1(max(n1[definition_day == 1 & schedule != "daily_29_57", pct_ge2])), f1(n1[definition_day == 1 & schedule == "daily_29_57", pct_ge3])), "") }

bt <- R("oc", "boundary_type1.csv"); pw <- R("oc", "power.csv"); rk <- R("oc", "random_space_risks.csv"); pr <- R("oc", "prereg.csv")
if (!is.null(bt)) {
  x <- bt[config == "P2"][which.max(pass_pct)]; g <- bt[config == "G2"][which.max(pass_pct)]
  add("Operating characteristics (pre-specified design)",
      if (!is.null(pr)) sprintf("- Design pre-registered in commit %s (oc/prereg.csv).", substr(pr$prereg_commit, 1, 7)),
      sprintf("- Largest boundary type I error of P2 (AUClast + Cmax): %s, %s, %s, true AUC0-inf ratio %s (oc/boundary_type1.csv, 10,000 trials).", ci(x$pass_pct, x$lo, x$hi, 2), ML[[x$model]], x$mechanism, f2(x$target)),
      sprintf("- Largest boundary type I error of G2 (AUCinf + Cmax): %s, %s, %s, true ratio %s.", ci(g$pass_pct, g$lo, g$hi, 2), ML[[g$model]], g$mechanism, f2(g$target)),
      sprintf("- Boundary scenarios with P2 above 5%%: %d of %d.", nrow(bt[config == "P2" & pass_pct > 5]), nrow(bt[config == "P2"])))
  if (!is.null(pw)) { s0 <- pw[scenario == "S00" & config %in% c("P2", "F3A", "G2")]
    out <<- c(out, sprintf("- Power for identical products (10,000 trials): %s (oc/power.csv).", paste(sprintf("%s %s %s", ML[s0$model], c(P2 = "P2", F3A = "F3-A", G2 = "G2")[s0$config], ci(s0$pass_pct, s0$lo, s0$hi)), collapse = "; "))) }
  if (!is.null(rk)) { r <- rk[truth == "AUC0-inf" & scope == "전체" & config %in% c("P2", "G2")]
    out <<- c(out, sprintf("- Random product space, consumer risk (truth outside, 20,000 products): %s (oc/random_space_risks.csv).",
                          paste(sprintf("%s %s %s", ML[r[startsWith(metric, "소비자")]$model], r[startsWith(metric, "소비자")]$config, ci(r[startsWith(metric, "소비자")]$pct, r[startsWith(metric, "소비자")]$lo, r[startsWith(metric, "소비자")]$hi, 2)), collapse = "; "))) }
  out <- c(out, "")
}

dec <- R("trials", "schedule_decision_base.csv"); k <- R("individual200k", "criterion_d_200k_base.csv")
if (!is.null(dec)) { d3 <- dec[schedule == "D3"]
  add("Sampling density decision (2016 model)",
      sprintf("- Final schedule B0. Candidates meeting any pre-specified criterion: %s (trials/schedule_decision_base.csv).", if (any(dec$recommend)) paste(dec[recommend == TRUE, schedule], collapse = ", ") else "none"),
      sprintf("- D3 versus B0: AUClast CI width change %s%% (positive = narrower), reliability gain %s percentage points, NCA extrapolation above 20%% ratio %s%s.", f2(100 * d3$a_mean_width_rel_decrease), f2(d3$c_reliable_gain_pp), f3(d3$d_extrap20_ratio),
              if (!is.null(k)) sprintf("; 200,000 subjects %s (95%% CI %s to %s)", f3(k[schedule == "D3", extrap_gt20_ratio]), f3(k[schedule == "D3", d_ratio_boot_lo]), f3(k[schedule == "D3", d_ratio_boot_hi])) else ""), "") }

p1 <- R("rationale", "pillar1_coverage_B0.csv")
if (!is.null(p1)) { x <- p1[group == "전체"]
  add("Coverage of total exposure by AUClast (B0, 60 to 90 kg, 20,000 subjects per model)",
      sprintf("- True extrapolated share median %s%% (2016) and %s%% (Model 1); 95th percentile %s%% and %s%%; true coverage below 80%%: %s and %s (rationale/pillar1_coverage_B0.csv).",
              f2(x[model == "k2016", extrap_true_median]), f2(x[model == "k2020", extrap_true_median]), f2(x[model == "k2016", extrap_true_p95]), f2(x[model == "k2020", extrap_true_p95]),
              x[model == "k2016", coverage_lt80_pct_ci], x[model == "k2020", coverage_lt80_pct_ci]), "") }

txt <- paste(out, collapse = "\n")
if (grepl("—", txt)) stop("key_numbers_en.md contains an em-dash")
if (grepl("[가-힣]", txt)) { bad <- regmatches(txt, gregexpr("[^\n]*[가-힣][^\n]*", txt))[[1]]; stop("key_numbers_en.md contains Korean text: ", paste(head(bad, 3), collapse = " || ")) }
writeLines(txt, proj_path("results", "key_numbers_en.md"), useBytes = TRUE)
cat(txt, "\n")
