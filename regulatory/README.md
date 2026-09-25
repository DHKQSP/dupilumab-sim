# Regulatory package: reviewer and inspector guide

Modeling and simulation supporting AUC0-last and Cmax as co-primary endpoints, and the sampling schedule, of a single-dose PK similarity study of a proposed dupilumab biosimilar.

Generated 2026-09-25 08:19 UTC from commit `ed9e8ba4f2961b345d77f0dd32788d740cbd2763` (tracked files outside regulatory/ unchanged: yes). Status: draft for sponsor review; not approved.

## Read in this order

| Document | Purpose |
|---|---|
| `MS_report.docx` / `MS_report.html` | Modeling and Simulation Report (ICH M15 structure): question of interest, context of use, model risk, methods, credibility evidence, results, discussion, appendices A to H. |
| `FDA_questions.docx` / `FDA_questions.html` | Anticipated FDA questions with evidence-based responses and pointers to the report. |
| `tables/assumptions_register.csv` | Assumptions, their impact and the actions required before submission. |
| `tables/prespecification_register.csv` | What was fixed before results, what was added afterwards (dates and commits). |
| `tables/verification_qc.csv` | Verification and QC activities, criteria, results and evidence. |
| `tables/parameter_provenance.csv` | Every model parameter with its published source. |
| `tables/software_environment.csv` | Software versions (renv.lock) and roles. |
| `traceability.csv` | 223 printed values (plus tables and figures), each with its source file, row filter, column, raw value and SHA-256. |
| `manifest_sha256.csv` | SHA-256 of every program, configuration, cited result and document in this package. |
| `tables/label_translation.csv` | Mapping of Korean category labels found in some analysis outputs to the English codes used in the documents. |

## Integrity and audit trail

- **Version control.** Every change to code, configuration and results is recorded in the git history of this repository. The operating-characteristic design was committed before any result (see Appendix D of the report).
- **Seeds and logs.** Every simulation run writes a log in `logs/` with the derived seeds, replicate counts, the SHA-256 of its configuration files and the session information.
- **Integrity of this package.** It is identified by the source commit above and by `manifest_sha256.csv`. For a submission, archive this commit (for example as a signed tag) and the manifest in the sponsor's validated document management system.
- **Attribution and human accountability.** The git history attributes changes to the AI coding assistant that prepared them under the sponsor's direction (commit author `Claude`). Human accountability is established by the signature page of the report after independent QC. Git by itself is not an electronic-signature system under 21 CFR Part 11.

## Reproduce

1. Restore the locked environment: R 4.3.3, then `renv::restore()`.
2. Run `bash scripts/run_all.sh` to regenerate all results. The simulations are long; `scripts/run_round6_assembly.sh` regenerates the summaries from the stored results without new simulation.
3. Run `Rscript scripts/38_repro_check.R` to check the committed results against pre-specified tolerances (`config/repro_check.yaml`).
4. Run `Rscript scripts/60_regulatory_package.R` to regenerate this package. Every number is read from the result files, and generation stops if a stated premise is contradicted.

## Before submission (sponsor actions)

- Resolve the items of `tables/assumptions_register.csv` marked for action, in particular the LLOQ of the sponsor's assay and the protocol sampling times and windows.
- Carry out independent QC of the report against `traceability.csv`, and sign the signature page.
- Re-check the wording of the FDA and ICH references against their current versions.
- Regenerate the package from a clean tree at the final commit, and archive it with its manifest.
