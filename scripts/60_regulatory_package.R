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
docs <- c("MS_report", "FDA_questions")
pars <- list(source_commit = commit, tree_clean = clean, version = "0.9 (draft for sponsor review)")
for (d in docs) {
  for (fmt in c("html", "docx")) {
    of <- switch(fmt, html = rmarkdown::html_document(toc = TRUE, toc_depth = 3, number_sections = FALSE, self_contained = TRUE),
                 docx = rmarkdown::word_document(toc = TRUE, toc_depth = 3))
    rmarkdown::render(file.path(src_dir, paste0(d, ".Rmd")), output_format = of, output_file = paste0(d, ".", fmt), output_dir = reg_dir,
                      intermediates_dir = src_dir, knit_root_dir = src_dir, envir = list2env(list(REG_PARAMS = pars)), quiet = TRUE)
  }
  # English check on the rendered text (HTML without tags and embedded images)
  h <- paste(readLines(file.path(reg_dir, paste0(d, ".html")), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
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
files <- sort(unique(c(keep, cited, file.path("regulatory", c(paste0(docs, ".html"), paste0(docs, ".docx"), "traceability.csv")),
                       file.path("regulatory", "tables", list.files(tab_dir)), file.path("regulatory", "figures", list.files(file.path(reg_dir, "figures"))))))
files <- files[file.exists(proj_path(files)) & !grepl("^regulatory/(manifest_sha256\\.csv|README\\.md)$", files)]
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
  sprintf("Generated %s UTC from commit `%s` (tracked files outside regulatory/ unchanged: %s). Status: draft for sponsor review; not approved.", format(Sys.time(), "%Y-%m-%d %H:%M", tz = "UTC"), commit, if (clean) "yes" else "no"),
  "",
  "## Read in this order",
  "",
  "| Document | Purpose |",
  "|---|---|",
  "| `MS_report.docx` / `MS_report.html` | Modeling and Simulation Report (ICH M15 structure): question of interest, context of use, model risk, methods, credibility evidence, results, discussion, appendices A to H. |",
  "| `FDA_questions.docx` / `FDA_questions.html` | Anticipated FDA questions with evidence-based responses and pointers to the report. |",
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
