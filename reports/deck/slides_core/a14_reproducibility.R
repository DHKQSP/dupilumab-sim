# A14 별첨 · 재현성. [모의]
# 카드 ① 핵심 덱 입력(scripts/63, 사전 등록 config/prereg_20260929.yaml section7, D-063): results/core_deck/provenance.csv의 대조 행(equal 열)과 입력 SHA-256 행.
#   곡선 모양 변형과 대표 대상자 모집단을 원래 시드 규칙으로 다시 만들어 커밋된 요약·저장된 NCA와 대조했다(허용 오차는 사전 등록).
#   사전 등록이 계산 전에 커밋되었는지는 scripts/63 실행 로그(logs/core_deck_inputs_*.log; git_sha, git_code_clean, 사전 등록 파일 SHA-256)로 전제 검사한다.
#   로그는 d* 출처로 쓸 수 없는 위치라(results/config/regulatory 밖) 수치를 읽지 않고 전제로만 쓴다(커밋 해시는 식별자).
# 카드 ② v1.0 결과 재현(results/repro/*; 결과보고 덱 A6과 같은 논리): 새 GitHub 러너, renv.lock 복원, 로컬과 유효숫자 8자리 일치. 항목은 v1.0 결과뿐이다.
# 카드 ③ 자동 시험(results/ci/*; A6과 같은 전제): 실행 72~101의 성공·취소·실패와 실패 원인, 고친 실행.
# 본문: 사전 등록 먼저, 시드 규칙(R/seeds.R derive_seed), renv.lock(regulatory/tables/software_environment.csv).

a14_field <- function(k, item = paste("GitHub reproducibility run:", k)) {
  x <- as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value)
  dderived(item, "repro/repro_github_run.csv", sprintf("field=='%s'", k), x, x)
}
a14_rx <- function(k, rx, item) {
  x <- as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value); m <- regmatches(x, regexec(rx, x))[[1]]
  premise(length(m) >= 2, sprintf("repro_github_run.csv field %s matches %s", k, rx))
  dderived(item, "repro/repro_github_run.csv", sprintf("field=='%s' :: regex %s", k, rx), m[2], m[2])
}
a14_ver <- function(sw, item) { SE <- "regulatory/tables/software_environment.csv"; x <- as.character(row1(SE, sprintf("Software=='%s'", sw))$Version)
  dderived(item, SE, sprintf("Software=='%s' :: Version", sw), x, x) }

slide_A14 <- function() {
  PV <- "core_deck/provenance.csv"; P29 <- "config/prereg_20260929.yaml"
  CB <- "ci/ci_failures_before.csv"; CS <- "ci/ci_failure_summary.csv"; CA <- "ci/ci_runs_after_fix.csv"; CV <- "ci/ci_runs_v101.csv"; CD <- "ci/ci_runs_deck.csv"
  RI <- "repro/repro_github_items.csv"; RC <- "repro/repro_check.csv"
  deck_slide("A14", tag = "sim")

  # ---- ① 핵심 덱 입력: 대조 행과 입력 행 ----
  pv <- rows(PV); inp <- startsWith(pv$check, "input ")
  premise(all(is.na(pv$equal[inp])) && all(pv$equal[!inp] %in% TRUE), "every identity check of scripts/63 is equal; input rows carry SHA-256 only")
  premise(all(grepl("^[0-9a-f]{64}$", pv$value[inp])), "input rows hold SHA-256 digests")
  RX <- list(cs = "grepl('^curve-shape .* regenerated vs curve_shape_B0', check)", nca = "grepl('regenerated B0 NCA equals stored', check)",
             sub = "grepl('true AUC at observed tlast|re-solved true concentration', check)",
             sum = "grepl('vs (tp_coverage_individual|lloq_individual_table|individual_resid12|tp_failure_by_set)', check)")
  cnt <- vapply(RX, function(w) nrow(rows(PV, w)), 1L)
  premise(sum(cnt) == sum(!inp) && all(cnt > 0), "the identity checks split into curve-shape, representative-population NCA, representative subjects and stored summaries")
  pr <- .read(P29)$section7
  seed_of <- function(s_) { m <- regmatches(s_, regexec("master seed ([0-9]+)", s_))[[1]]; premise(length(m) == 2, "master seed in the prereg text"); as.numeric(m[2]) }
  sd_rep <- seed_of(pr$representative_subjects$regeneration); sd_cs <- seed_of(pr$case_coverage$cases[[3]]$source)
  premise(sd_rep == sd_cs && sd_rep == .read("config/oc_design.yaml")$trials$master_seed, "same master seed for the curve-shape and representative-subject regeneration (the original seed rule)")
  # 사전 등록이 계산 전에 커밋됐다: 마지막 실행 로그(이름의 UTC 시각 순)가 읽은 사전 등록 SHA-256 = 현재 파일, 그 파일은 커밋된 상태(커밋 전 변경 목록에 없음)
  lg <- sort(list.files(proj_path("logs"), pattern = "^core_deck_inputs_[0-9]{8}T[0-9]{6}\\.log$", full.names = TRUE)); premise(length(lg) >= 1, "a run log of scripts/63")
  ll <- readLines(lg[length(lg)], warn = FALSE, encoding = "UTF-8")
  premise(any(ll == sprintf("master_seed: %s", format(sd_rep, scientific = FALSE))) && any(grepl("^\\[[0-9:]+\\] done$", ll)), "the latest scripts/63 run used the registered master seed and finished")
  sh <- sub("^  prereg_20260929.yaml: ", "", grep("^  prereg_20260929.yaml: ", ll, value = TRUE))
  premise(length(sh) == 1 && identical(sh, digest::digest(file = proj_path(P29), algo = "sha256")), "the registration read by the latest scripts/63 run is the current config/prereg_20260929.yaml")
  unc <- sub("^  uncommitted: [A-Z?]+ ", "", grep("^  uncommitted: ", ll, value = TRUE)); clean63 <- any(ll == "git_code_clean: TRUE")
  premise(clean63 == (length(unc) == 0) && !any(grepl("prereg_20260929", unc, fixed = TRUE)), "the registration was committed when scripts/63 ran; the clean flag matches the list of uncommitted files")
  cm63 <- substr(sub("^git_sha: ", "", grep("^git_sha: [0-9a-f]{40}$", ll, value = TRUE)), 1, 7); premise(length(cm63) == 1 && nchar(cm63) == 7, "commit of the scripts/63 run")
  d63 <- readLines(proj_path("DECISIONS.md"), warn = FALSE, encoding = "UTF-8")
  premise(any(grepl("앞선 세 번의 실행은 결과 전에 멈추거나 대체됐다", d63, fixed = TRUE) & grepl("사전 등록에 없는 조건이라 뺌", d63, fixed = TRUE) & grepl("세 로그는 지웠고", d63, fixed = TRUE)),
          "DECISIONS D-063: three earlier runs stopped or were superseded (one on a premise not in the registration, dropped); their logs were deleted")
  seed <- dderived("master seed of the regeneration (original seed rule of scripts/10 and scripts/21)", P29,
                   "section7.representative_subjects.regeneration :: regex master seed ([0-9]+)", sd_rep, format(sd_rep, scientific = FALSE))
  f1 <- list(neq = dcount(PV, "equal == TRUE", "scripts/63 identity checks equal to the committed results"),
             nchk = dcount(PV, "!startsWith(check, 'input ')", "scripts/63 identity checks"),
             nin = dcount(PV, "startsWith(check, 'input ')", "scripts/63 inputs recorded by SHA-256"),
             ncs = dcount(PV, RX$cs, "curve-shape draws regenerated and compared (base and four variants)"),
             nind = dcfg("trial_design.yaml", c("mc", "n_individual"), "subjects per model in the individual populations", function(x) fnum(as.numeric(x), 0, big = TRUE)),
             seed = seed, regd = dcfg("prereg_20260929.yaml", "registered_on", "registration date of the core-deck pre-registration", function(x) as.character(x)))

  # ---- ② v1.0 결과 재현(새 GitHub 러너) ----
  ri <- rows(RI); rc <- rows(RC)
  rgf <- function(k) as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value)
  premise(all(ri$pass_github) && all(ri$identical_to_local_8sig) && all(rc$pass) && nrow(ri) == nrow(rc), "every reproducibility item passes locally and on the clean runner, identical to 8 significant digits")
  premise(as.numeric(rgf("n_items")) == nrow(ri) && rgf("n_pass_github") == rgf("n_items") && rgf("n_identical_to_local_8sig") == rgf("n_items") && rgf("conclusion") == "success",
          "clean-runner run: all items pass and are identical to the local values")
  premise(grepl("GitHub-hosted", rgf("runner")) && rgf("run_number") == "1" && rgf("run_attempt") == "1", "first run of the repro workflow on a GitHub-hosted runner")
  premise(grepl("^hit", rgf("renv_cache")), "renv restored the locked library on the runner")
  premise(all(grepl("^(a[0-9]_(oc_trials|P2_pass_rate)|b[0-9]_(true_AUCinf_ratio_Vmax|screening_consistency)|c1_extrap_true_median)", ri$item)), "reproducibility items: v1.0 operating-characteristic trials, boundary multiplier search, extrapolation median")
  PR <- "regulatory/tables/prespecification_register.csv"
  t_run <- as.POSIXct(rgf("job_started_utc"), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"); preg <- rows(PR)[grepl("^pre-registered", Status)]
  t_reg <- min(as.POSIXct(sub(" UTC.*", "", preg$`Date (evidence)`), tz = "UTC", format = "%Y-%m-%d %H:%M"))
  premise(!is.na(t_run) && !is.na(t_reg) && t_run < t_reg, "the reproducibility run precedes the v1.0.1 pre-registrations (caption: v1.0 results only)")
  f2 <- list(pg = sprintf("%s/%s", a14_field("n_pass_github"), a14_field("n_items")),
             ni = a14_field("n_items", "GitHub reproducibility run: items (all pass)"),
             when = a14_rx("job_started_utc", "^([0-9]{4}-[0-9]{2}-[0-9]{2})T", "GitHub reproducibility run: job start date (UTC)"),
             cm = a14_rx("commit", "^([0-9a-f]{7})", "GitHub reproducibility run: commit (short hash)"),
             sig = dderived("significant digits in column name identical_to_local_8sig", RI, "column name identical_to_local_8sig", 8, "8"))

  # ---- ③ 자동 시험(GitHub Actions tests 워크플로) ----
  cb <- rows(CB); ca <- rows(CA, "workflow=='tests'"); cv <- rows(CV); cd <- rows(CD, "workflow=='tests'")
  runs <- rbind(cb[, .(run = run_number, conc = "failure")], ca[, .(run = run_number, conc = conclusion)], cv[, .(run = run_number, conc = conclusion)], cd[, .(run = run_number, conc = conclusion)])
  premise(!anyDuplicated(runs$run) && identical(sort(runs$run), seq_len(max(runs$run))) && all(runs$conc %in% c("success", "failure", "cancelled")), "tests runs 1 to the latest recorded run are each recorded once")
  premise(max(cb$run_number) + 1 == min(ca$run_number) && max(ca$run_number) + 1 == min(cv$run_number) && max(cv$run_number) + 1 == min(cd$run_number), "the four CI files are contiguous")
  premise(all(cb$category == "package install/build"), "runs before the fix failed in environment setup (notes)")
  fl <- cv[conclusion == "failure"]; premise(nrow(fl) > 0 && all(diff(fl$run_number) == 1) && all(grepl("^strict-skip rule", fl$note)), "v1.0.1 failures are consecutive and due to the strict-skip rule")
  fxr <- cv[grepl("^fix of runs", note)]; premise(nrow(fxr) == 1 && fxr$conclusion == "success" && fxr$run_number == max(fl$run_number) + 1, "the run after the strict-skip failures fixed them")
  dfl <- cd[conclusion == "failure"]; premise(nrow(dfl) == 1 && grepl("^environment setup", dfl$note), "one environment-setup failure in the deck record")
  dfx <- cd[grepl("^fix of run", note)]; premise(nrow(dfx) == 1 && dfx$conclusion == "success" && dfx$run_number == dfl$run_number + 1, "the run after it fixed it")
  premise(runs[run == max(run), conc] == "success", "the latest recorded run succeeded (card)")
  cnt3 <- function(k) { a_ <- nrow(ca[conclusion == k]); v_ <- nrow(cv[conclusion == k]); d_ <- nrow(cd[conclusion == k])
    dderived(sprintf("tests runs after the fix with conclusion %s (after_fix + v101 + deck)", k), CV,
             sprintf("count(conclusion=='%s') in ci_runs_after_fix.csv (workflow tests) + ci_runs_v101.csv + ci_runs_deck.csv (workflow tests)", k), c(a_, v_, d_), as.character(a_ + v_ + d_)) }
  dsrc("tests runs after the fix (with ci_runs_v101.csv)", CA, "(table)"); dsrc("tests runs after the fix (with ci_runs_v101.csv)", CD, "(table)")
  f3 <- list(fx1 = dint(CA, sprintf("workflow=='tests' & run_number==%d", min(ca$run_number)), "run_number", "first tests run after the environment fix"),
             lastd = dext(CD, "workflow=='tests'", "run_number", max, 0, "", "latest tests run recorded"),
             ns = cnt3("success"), nc = cnt3("cancelled"), nf = cnt3("failure"),
             rf = drange(CV, "conclusion=='failure'", "run_number", 0, "", "failed runs (strict-skip rule), run numbers"),
             fx2 = dint(CV, "grepl('^fix of runs', note)", "run_number", "run that fixed the strict-skip failures"),
             rd = dint(CD, "workflow=='tests' & conclusion=='failure'", "run_number", "failed run of the deck record (environment setup)"),
             fx3 = dint(CD, "workflow=='tests' & grepl('^fix of run', note)", "run_number", "run that fixed the environment-setup failure"))
  premise(as.numeric(f3$nf) == nrow(fl) + nrow(dfl), "failures after the fix = strict-skip failures + the environment-setup failure")

  y0 <- core_title(tx("A14.title", list(neq = f1$neq, ni = f2$ni)), tx("A14.kicker"))

  # ---- 카드 3개 ----
  C <- DK$txt$A14$cards; gap <- 0.22; cw <- (GEO$CW - 2 * gap) / 3; ch <- 2.95; cy <- y0 + 0.08
  v2 <- function(a, b) list(list(a, PAL$ink, 40), list(b, PAL$ink2, 20))
  core_card(C$c1$head, v2(sprintf("%s/%s", f1$neq, f1$nchk), C$c1$unit), fill(C$c1$label, f1), c(GEO$ML, cy, cw, ch), bg = PAL$tint_grey)
  core_card(C$c2$head, v2(f2$pg, C$c2$unit), fill(C$c2$label, f2), c(GEO$ML + cw + gap, cy, cw, ch), bg = PAL$tint_grey)
  core_card(C$c3$head, list(list(f3$ns, PAL$ink, 40), list(C$c3$ok, PAL$ink2, 20), list(f3$nf, PAL$ink, 40), list(C$c3$fail, PAL$ink2, 20)),
            fill(C$c3$label, f3), c(GEO$ML + 2 * (cw + gap), cy, cw, ch), bg = PAL$tint_grey)
  deck_visual(c(GEO$ML, cy, GEO$CW, ch))

  # ---- 본문(사전 등록, 시드, renv.lock)과 캡션 ----
  fb <- list(regd = f1$regd, seed = f1$seed, rv = a14_ver("R", "R version (renv.lock)"), rx = a14_ver("rxode2", "rxode2 version (renv.lock)"))
  capy <- core_caption(tx("A14.caption", list()), GEO$BODY_BOTTOM, size = 14)
  body <- tx("A14.body", fb); bh <- core_body_h(body, gap_pt = 6)
  premise(cy + ch + 0.2 <= capy - 0.08 - bh, "body fits between the cards and the caption")
  core_body(body, capy - 0.08, gap_pt = 6)

  run63 <- if (clean63) tx("A14.run_clean", list(cm63 = cm63)) else tx("A14.run_dirty", list(cm63 = cm63, unc = paste(unc, collapse = ", ")))
  deck_notes(tx("A14.notes", c(f1, f2, f3, fb, list(
    cm63 = cm63, run63 = run63, ncsv = dcount(PV, "grepl('^curve-shape (vmax|km)', check)", "curve-shape variants regenerated (without the base draw)"),
    nnca = dcount(PV, RX$nca, "regenerated representative-population NCA column checks (two models)"),
    ncol = { x <- nrow(rows(PV, RX$nca)) / 2; dderived("NCA columns compared per model", PV, sprintf("count of rows [%s] / 2 models", RX$nca), x, fnum(x, 0)) },
    nsub = dcount(PV, RX$sub, "representative-subject checks (true AUC at tlast, re-solved concentration)"),
    nsum = dcount(PV, RX$sum, "identity checks against stored summaries"),
    r1 = drange(CB, "TRUE", "run_number", 0, "", "runs before the environment fix, run numbers"),
    a = dint(CS, "cause=='no renv.lock in commit'", "n_runs", "runs failed: no renv.lock in commit"),
    b = dint(CS, "cause=='xml2 missing for testthat::JunitReporter'", "n_runs", "runs failed: xml2 missing for JunitReporter"),
    nl = dcount(RC, "pass==TRUE", "reproducibility items passing locally"),
    os = a14_rx("runner", "^(ubuntu-[0-9.]+) ", "GitHub runner image"),
    renv = a14_ver("renv", "renv version (renv.lock)")))))
  deck_end()
}
