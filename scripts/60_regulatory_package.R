#!/usr/bin/env Rscript
# FDA 검토·실사 대비 패키지(regulatory/) 생성. 새 모의 없음: 커밋된 결과·config에서 문서를 만든다.
#  1) 부록 표(regulatory/tables/*.csv) — regulatory/src/build_tables.R
#  2) M&S 보고서(ICH M15 구조)와 FDA 예상 질의응답을 HTML·DOCX로 렌더링 — 모든 수치는 결과 파일에서 읽고 출처를 기록
#  3) 추적표 regulatory/traceability.csv(수치 → 파일·행 조건·열·SHA-256), 영문 규칙 검사(한글·en/em dash·U+2212 금지)
#  4) SHA-256 목록 regulatory/manifest_sha256.csv(프로그램·config·결과·문서), 검토자 안내 regulatory/README.md
# 문서에는 생성 시점의 커밋 해시와 작업 트리 청결 여부를 적는다. 제출용 확정본은 청결한 트리에서 생성한다.
# 사용법: Rscript scripts/60_regulatory_package.R
source("R/00_setup.R"); source(proj_path("regulatory", "src", "reg_helpers.R"))
reg_dir <- proj_path("regulatory"); src_dir <- file.path(reg_dir, "src"); tab_dir <- file.path(reg_dir, "tables")
git <- function(...) suppressWarnings(system2("git", c("-C", PROJ_ROOT, ...), stdout = TRUE, stderr = FALSE))
commit <- git("rev-parse", "HEAD")[1]
dirty <- git("status", "--porcelain", "--untracked-files=no")
dirty <- dirty[!grepl("^.. regulatory/", dirty)]          # regulatory/ 산출물 자체의 변경은 제외(생성 대상)
clean <- length(dirty) == 0
cat(sprintf("source commit %s, tracked tree outside regulatory/ clean: %s\n", commit, clean))

# 1) tables
sys.source(file.path(src_dir, "build_tables.R"), envir = new.env())

# 2) documents
docs <- c("MS_report", "FDA_questions", "SAP_text_proposals")
out_name <- c(MS_report = "MS_report", FDA_questions = "FDA_questions", SAP_text_proposals = "sap_text_proposals_en")   # 출력 파일 이름
doc_formats <- list(MS_report = c("html", "docx"), FDA_questions = c("html", "docx"), SAP_text_proposals = c("html", "md"))
VERSION <- "1.0.1 (draft for sponsor review)"
pars <- list(source_commit = commit, tree_clean = clean, version = VERSION)
for (d in docs) {
  for (fmt in doc_formats[[d]]) {
    of <- switch(fmt, html = rmarkdown::html_document(toc = TRUE, toc_depth = 3, number_sections = FALSE, self_contained = TRUE),
                 docx = rmarkdown::word_document(toc = TRUE, toc_depth = 3), md = rmarkdown::md_document(variant = "gfm", toc = TRUE, toc_depth = 2))
    rmarkdown::render(file.path(src_dir, paste0(d, ".Rmd")), output_format = of, output_file = paste0(out_name[[d]], ".", fmt), output_dir = reg_dir,
                      intermediates_dir = src_dir, knit_root_dir = src_dir, envir = list2env(list(REG_PARAMS = pars)), quiet = TRUE)
  }
  if ("md" %in% doc_formats[[d]]) check_english(file.path(reg_dir, paste0(out_name[[d]], ".md")))
  # English check on the rendered text (HTML without tags and embedded images)
  h <- paste(readLines(file.path(reg_dir, paste0(out_name[[d]], ".html")), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  h <- gsub("data:[a-z/+]+;base64,[A-Za-z0-9+/=]+", "", h); h <- gsub("<script[^>]*>.*?</script>", "", h, perl = TRUE); h <- gsub("<style[^>]*>.*?</style>", "", h, perl = TRUE)
  tf <- tempfile(fileext = ".txt"); writeLines(gsub("<[^>]+>", " ", h), tf, useBytes = TRUE); check_english(tf)
  cat("rendered", d, "\n")
}

# 3) traceability
tr <- rbindlist(lapply(docs, function(d) fread(file.path(tab_dir, sprintf("trace_%s.csv", d)))), fill = TRUE)
tr[, source_commit := commit]
fwrite(tr, file.path(reg_dir, "traceability.csv"))
unlink(file.path(tab_dir, sprintf("trace_%s.csv", docs)))
tf <- tempfile(fileext = ".csv"); fwrite(tr[, .(document, section, item, source_file, locator, value_printed)], tf); check_english(tf)

# 4) manifest: programs, configuration, dependency lock, tests, every result file cited, figures and the package itself
tracked <- git("ls-files")
keep <- tracked[grepl("^(R|scripts|config|tests|regulatory)/", tracked) | tracked %in% c("renv.lock", "renv/settings.json", "renv/activate.R", ".Rprofile", "SPEC.md", "DECISIONS.md", "README.md")]
cited <- unique(tr$source_file[file.exists(proj_path(tr$source_file))])
rendered <- unlist(lapply(docs, function(d) paste0(out_name[[d]], ".", doc_formats[[d]])))
files <- sort(unique(c(keep, cited, file.path("regulatory", c(rendered, "traceability.csv")),
                       file.path("regulatory", "tables", list.files(tab_dir)), file.path("regulatory", "figures", list.files(file.path(reg_dir, "figures"))))))
files <- files[file.exists(proj_path(files)) & !grepl("^regulatory/(manifest_sha256\\.csv|README\\.md|release_notes\\.md)$", files)]
mf <- data.table(path = files, bytes = file.size(proj_path(files)), sha256 = vapply(files, function(f) digest::digest(file = proj_path(f), algo = "sha256"), ""),
                 git_tracked = files %in% tracked, role = fifelse(grepl("^regulatory/", files), "regulatory document", fifelse(grepl("^results/", files), "result cited in the documents",
                                                             fifelse(grepl("^(R|scripts)/", files), "program", fifelse(grepl("^config/", files), "configuration", fifelse(grepl("^tests/", files), "test", "project record"))))))
fwrite(mf, file.path(reg_dir, "manifest_sha256.csv"))

# 5) reviewer guide
n_vals <- nrow(tr[value_printed != "(table)" & value_printed != "(figure)"])
readme <- c(
  "# Regulatory package: reviewer and inspector guide",
  "",
  "Modeling and simulation supporting AUC0-last and Cmax as co-primary endpoints, and the sampling schedule, of a single-dose PK similarity study of a proposed dupilumab biosimilar.",
  "",
  sprintf("Version %s. Generated %s UTC from commit `%s` (tracked files outside regulatory/ unchanged: %s). Status: draft for sponsor review; not approved.", sub(" .*$", "", VERSION), format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC"), commit, if (clean) "yes" else "no"),
  "",
  "## Read in this order",
  "",
  "| Document | Purpose |",
  "|---|---|",
  "| `MS_report.docx` / `MS_report.html` | Modeling and Simulation Report (ICH M15 structure): question of interest, context of use, model risk, methods, credibility evidence, results, discussion, appendices A to H. |",
  "| `FDA_questions.docx` / `FDA_questions.html` | Anticipated FDA questions with evidence-based responses and pointers to the report. |",
  "| `sap_text_proposals_en.md` | Proposed statistical analysis plan text for the PK analyses (primary endpoints and analysis model options, AUC0-inf as secondary endpoint, fallback, ADA, BLQ and AUC rules, sample size), with the supporting numbers. |",
  "| `release_notes.md` | Changes in this version. |",
  "| `tables/assumptions_register.csv` | Assumptions, their impact and the actions required before submission. |",
  "| `tables/prespecification_register.csv` | What was fixed before results, what was added afterwards (dates and commits). |",
  "| `tables/verification_qc.csv` | Verification and QC activities, criteria, results and evidence. |",
  "| `tables/parameter_provenance.csv` | Every model parameter with its published source. |",
  "| `tables/software_environment.csv` | Software versions (renv.lock) and roles. |",
  sprintf("| `traceability.csv` | %s printed values (plus tables and figures), each with its source file, row filter, column, raw value and SHA-256. |", format(n_vals, big.mark = ",")),
  "| `manifest_sha256.csv` | SHA-256 of every program, configuration, cited result and document in this package. |",
  "| `tables/label_translation.csv` | Mapping of Korean category labels found in some analysis outputs to the English codes used in the documents. |",
  "",
  "## Integrity and audit trail",
  "",
  "- **Version control.** Every change to code, configuration and results is recorded in the git history of this repository. The operating-characteristic design was committed before any result (see Appendix D of the report).",
  "- **Seeds and logs.** Every simulation run writes a log in `logs/` with the derived seeds, replicate counts, the SHA-256 of its configuration files and the session information.",
  "- **Integrity of this package.** It is identified by the source commit above and by `manifest_sha256.csv`. For a submission, archive this commit (for example as a signed tag) and the manifest in the sponsor's validated document management system.",
  "- **Attribution and human accountability.** The git history attributes changes to the AI coding assistant that prepared them under the sponsor's direction (commit author `Claude`). Human accountability is established by the signature page of the report after independent QC. Git by itself is not an electronic-signature system under 21 CFR Part 11.",
  "",
  "## Reproduce",
  "",
  "1. Restore the locked environment: R 4.3.3, then `renv::restore()`.",
  "2. Run `bash scripts/run_all.sh` to regenerate all results. The simulations are long; `scripts/run_round6_assembly.sh` regenerates the summaries from the stored results without new simulation.",
  "3. Run `Rscript scripts/38_repro_check.R` to check the committed results against pre-specified tolerances (`config/repro_check.yaml`).",
  "4. Run `Rscript scripts/60_regulatory_package.R` to regenerate this package. Every number is read from the result files, and generation stops if a stated premise is contradicted.",
  "",
  "## Before submission (sponsor actions)",
  "",
  "- Resolve the items of `tables/assumptions_register.csv` marked for action, in particular the LLOQ of the sponsor's assay and the protocol sampling times and windows.",
  "- Carry out independent QC of the report against `traceability.csv`, and sign the signature page.",
  "- Re-check the wording of the FDA and ICH references against their current versions.",
  "- Regenerate the package from a clean tree at the final commit, and archive it with its manifest."
)
writeLines(readme, file.path(reg_dir, "README.md"))
check_english(file.path(reg_dir, "README.md"))
cat(sprintf("regulatory package: %d traced values, %d files in the manifest\n", n_vals, nrow(mf)))

# 6) release notes (regulatory/release_notes.md; the release workflow requires the first line "# v<version>")
tag <- paste0("v", sub(" .*$", "", VERSION))
t1 <- fread(proj_path("results", "oc_models", "type1_models.csv")); de <- fread(proj_path("results", "oc_models", "p2_decomposition_models.csv"))
rv <- fread(proj_path("results", "oc_models", "m0_reverification.csv")); ex <- fread(proj_path("results", "oc_models", "expectations_check.csv"))
nn <- fread(proj_path("results", "sample_size", "ss_table_n_needed.csv")); tp <- fread(proj_path("results", "sample_size", "ss_table_power.csv"))
li <- fread(proj_path("results", "lloq", "lloq_individual_table.csv")); lt <- fread(proj_path("results", "lloq", "lloq_trial_type1.csv"))
f1 <- function(x) formatC(round(x + sign(x) * 1e-9, 1) + 0, format = "f", digits = 1); f2 <- function(x) formatC(round(x + sign(x) * 1e-9, 2) + 0, format = "f", digits = 2)
cls <- function(a) { x <- t1[analysis_model == a & config == "P2"]; sprintf("%d conservative, %d nominal, %d exceeding (%s%% to %s%%)", sum(x$class == "conservative"), sum(x$class == "nominal"), sum(x$class == "exceeding"), f2(min(x$pass_pct)), f2(max(x$pass_pct))) }
ov <- de[p2_pct > 5]
ov_l <- if (nrow(ov)) paste(sprintf("%s %s under %s: %s%% [%s, %s], %s, %s trials", ifelse(ov$pk_model == "k2016", "primary model", "2020 model"), ov$scenario, ov$analysis_model, f2(ov$p2_pct), f2(ov$lo), f2(ov$hi), ov$class, format(ov$n_trials, big.mark = ",")), collapse = "; ") else "none"
p117 <- function(cv_, a) f1(tp[input_model == "k2016" & cv == cv_ & gmr == 0.95 & n == 117 & analysis_model == a, analytic_pct])
nn_ <- function(cv_, a, tg) nn[cv == cv_ & gmr == 0.95 & analysis_model == a & target_pct == tg, n_evaluable_per_arm]
lr <- function(a) { x <- lt[model == a & resid == "fixed" & config == "P2", pass_pct]; sprintf("%s%% to %s%%", f2(min(x)), f2(max(x))) }
lk <- li[model == "k2016" & resid == "fixed"]
cri <- fread(proj_path("results", "criteria", "criteria_individual.csv")); cg <- fread(proj_path("results", "criteria", "criteria_g2_type1.csv")); rc <- fread(proj_path("results", "criteria", "restart_identity_check.csv"))
vc <- fread(proj_path("results", "criteria", "criteria_check_vs_v10.csv")); acb <- fread(proj_path("results", "atopic", "atopic_criteria_by_band.csv")); acf <- fread(proj_path("results", "atopic", "atopic_coverage_failing.csv"))
frr <- function(dist, s_) { x <- cri[distribution == dist & set == s_, fail_pct]; sprintf("%s%% to %s%%", f1(min(x)), f1(max(x))) }
g2n <- function(cf) sum(cg[analysis_model == "M0" & config == cf, pass_pct] > 5)
atX <- f1(acb[variant == "base" & distribution == "primary" & band == "all" & set == "i", fail_pct]); atY <- f1(acb[variant == "base" & distribution == "primary" & band == "all" & set == "iii", fail_pct])
atZ <- floor(100 * min(acf[distribution == "primary" & set == "iii", window_fail_min]))
tpf <- fread(proj_path("results", "trialpop", "tp_failure_by_set.csv")); tpa <- fread(proj_path("results", "trialpop", "tp_retained_per_arm.csv")); tad <- fread(proj_path("results", "trialpop", "tp_arm_difference.csv"))
tsi <- fread(proj_path("results", "trialpop", "tp_strata_individual.csv")); tid <- fread(proj_path("results", "trialpop", "tp_identity_check.csv")); tga <- fread(proj_path("results", "trialpop", "tp_gmr_agreement.csv"))
r2 <- function(x, f = f1, u = "%") if (identical(f(min(x)), f(max(x)))) sprintf("%s%s", f(min(x)), u) else sprintf("%s%s to %s%s", f(min(x)), u, f(max(x)), u)
rn <- c(
  paste("#", tag), "",
  sprintf("Regulatory package version %s, generated %s UTC from commit `%s`. Status: draft for sponsor review; not approved. Previous package: document version 0.9 (commit ed9e8ba, not tagged).",
          sub(" .*$", "", VERSION), format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC"), commit), "",
  "## Pre-registration", "",
  "- Analyses of this version were registered in `config/prereg_20260926.yaml` and committed before their results: analysis model (commit 521645a), numeric criteria for the expectation checks (c5b304b, before any result was read), LLOQ sensitivity and sample size (211e8ed), lambda-z criteria convention and adult atopic dermatitis body weight (325e63b), trial-population analyses of the correction directive (bae3574).", "",
  "## Results added", "",
  sprintf("- Analysis model (report Section 5.9). Boundary type I error of AUC0-last + Cmax in 16 cells: pooled t-test (M0) %s; ANOVA with the randomization weight stratum (M1) %s. Cells above 5%% (point estimate): %s.", cls("M0"), cls("M1"), ov_l),
  sprintf("- Expectations recorded before the results: %s.", paste(sprintf("%s (%s)", ex$expectation, ifelse(ex$consistent, "consistent", "not consistent")), collapse = "; ")),
  sprintf("- LLOQ (report Section 5.10). The study LLOQ is set in `config/assay.yaml` (single source). Between 0.02 and 0.5 mg/L: window coverage of AUC0-last below 80%% in %s%% to %s%% of subjects; AUC0-inf reliability, criteria (i), %s%% (0.02 mg/L) to %s%% (0.5 mg/L); boundary type I error of AUC0-last + Cmax %s (M0) and %s (M1) in three boundary scenarios.",
          f2(min(lk$coverage_lt80_pct)), f2(max(lk$coverage_lt80_pct)), f1(lk[abs(lloq - 0.02) < 1e-9, reliable_no_span_pct]), f1(lk[abs(lloq - 0.5) < 1e-9, reliable_no_span_pct]), lr("M0"), lr("M1")),
  sprintf("- Sample size (report Section 5.11). At a true GMR of 0.95, 117 evaluable subjects per arm give P2 power %s%% (M0) and %s%% (M1) at CV 43%%, %s%% and %s%% at CV 50%%. Evaluable per arm for 90%% power: %d (M0) and %d (M1) at CV 43%%, %d and %d at CV 50%%.",
          p117(43, "M0"), p117(43, "M1"), p117(50, "M0"), p117(50, "M1"), nn_(43, "M0", 90), nn_(43, "M1", 90), nn_(50, "M0", 90), nn_(50, "M1", 90)),
  sprintf("- Lambda-z criteria convention (report Sections 3.4, 5.3, 5.4). Phoenix WinNonlin applies no reliability criteria unless the user enters them; adjusted R-squared 0.80 is the lower end of the conventional range. Study population without a reliable AUC0-inf: %s (0.80), %s (0.90), %s (0.90 and span 3). AUC0-inf + Cmax with rule A exceeds 5%% in %d, %d, %d and %d of 16 boundary scenarios under criteria sets (i) to (iv) (M0).",
          frr("study_60_90", "i"), frr("study_60_90", "iii"), frr("study_60_90", "iv"), g2n("G2_A_i"), g2n("G2_A_ii"), g2n("G2_A_iii"), g2n("G2_A_iv")),
  sprintf("- Trial population (report Sections 5.2 and 5.3; healthy adults, 60 to 90 kg, weight-stratified randomization, two models). Without a reliable AUC0-inf: %s (criteria set (i)) and %s (set (iii)); retained per arm of 117 under rule A: median %s (set (i)) and %s (set (iii)). The heavier stratum fails more often (set (i) %s percentage points), and the share failing differs between arms when the products differ (Vmax x1.25: test minus reference %s percentage points, set (i)). The geometric mean of trial AUC0-last GMRs is within %s of the true AUC0-inf ratio in every scenario examined.",
          r2(tpf[set == "i", fail_pct]), r2(tpf[set == "iii", fail_pct]), paste(unique(tpa[set == "i", retained_median]), collapse = " to "), paste(unique(tpa[set == "iii", retained_median]), collapse = " to "),
          r2(tsi[set == "i", diff_pp], f2, ""), r2(tad[scenario == "VM125" & set == "i", diff_mean], f2, ""), formatC(max(tga$abs_diff), format = "f", digits = 3)),
  "- Robustness across body weight (report Appendix I): uniform weight bands, a population with a mean of 100 kg and an adult atopic dermatitis body-weight distribution are robustness checks, not evidence for the proposal (sponsor principle, correction directive).",
  "- Proposed statistical analysis plan text for the PK analyses: `sap_text_proposals_en.md` (criteria convention added to Section 2).",
  "- Anticipated FDA questions: Q1, Q7, Q11, Q12 and Q16 updated; Q17 (analysis model and stratification), Q18 (assay LLOQ), Q19 (is 0.80 a Phoenix WinNonlin criterion) and Q20 (patients outside the study weight range: validity judged in the trial population, wider weight range a robustness check) added.",
  "- Terms: window coverage (true AUC0-tlast / true AUC0-inf) and the observed-to-true ratio (observed AUC0-last or NCA AUC0-inf / true AUC0-inf) are defined and used consistently.",
  "- Report tables and figures of Section 5 renumbered in order (a duplicated table number in version 1.0 corrected).", "",
  "## Decisions", "",
  "- Confirmed: AUC0-inf as a secondary endpoint (Best Fit lambda-z; adjusted R-squared at least 0.80 and extrapolation at most 20%, no span criterion; two analysis sets with excluded subjects and their body weight; substitution rule as sensitivity); disclosure of the boundary cell above 5%; fallback if AUC0-inf is required as co-primary (all subjects with an estimable lambda-z, flagged subjects listed; exclusion and substitution as sensitivity analyses); CV 43% as the base assumption and 50% as sensitivity.",
  "- Proposed (sponsor to confirm): adjusted R-squared at least 0.80 as the reliability definition for the secondary AUC0-inf, with the count at 0.90 reported alongside.",
  "- Pending (sponsor): primary analysis model (M0 or M1); LLOQ of the validated assay; target power and sample size.", "",
  "## Verification", "",
  sprintf("- The regenerated pooled t-test results equal the stored results in %s rows (results/oc_models/m0_reverification.csv). LLOQ re-censoring at 0.078 mg/L reproduces the stored individual, cliff and trial results. Automated tests pass.", format(sum(rv$n_rows), big.mark = ",")),
  sprintf("- Regeneration restarted to add the criteria-set endpoints: %s rows of the trials completed before the restart are identical after it (results/criteria/restart_identity_check.csv). AUC0-inf + Cmax under rules A, B and C from the new files equals the version 1.0 values in all %d cells with equal trial counts (results/criteria/criteria_check_vs_v10.csv).",
          format(sum(rc$rows_matched), big.mark = ","), sum(vc$same_n)),
  sprintf("- Trial-population regeneration: per-arm counts equal the section1 n_R and n_T in all %s compared rows (results/trialpop/tp_identity_check.csv).", format(sum(tid$rows_compared), big.mark = ",")), "",
  "## Release assets", "",
  sprintf("- `dupilumab-sim-%s-source.tar.gz`: repository at the tag (git archive).", tag),
  sprintf("- `dupilumab-sim-%s-regulatory.zip`: the regulatory/ folder.", tag),
  sprintf("- `dupilumab-sim-%s-manifest_sha256.csv`: SHA-256 of every program, configuration, cited result and document.", tag),
  "- `SHA256SUMS.txt`: SHA-256 of the release assets.")
writeLines(rn, file.path(reg_dir, "release_notes.md")); check_english(file.path(reg_dir, "release_notes.md"))
cat("release notes written for", tag, "\n")
