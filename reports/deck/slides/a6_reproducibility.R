# A6 부록: 재현성. 자동 시험 기록(GitHub Actions tests 워크플로), 깨끗한 환경 재현 워크플로, 추적표 규모, 사전 명시·등록 이력. [모의]
# 자동 시험: results/ci/ci_failures_before.csv(실행 1~71), ci_runs_after_fix.csv(72~77, repro 1), ci_runs_v101.csv(78~98; GitHub Actions API, 2026-09-29 읽음).
# 재현: results/repro/repro_github_run.csv(field/value 문자열), repro_github_items.csv, repro_check.csv(로컬).
# 추적 규모: regulatory/tables/trace_summary.csv(scripts/60이 덱 행을 합치기 전에 쓰는 문서별 수; 덱은 regulatory/traceability.csv와
#   manifest_sha256.csv를 인용하지 않는다: 그 둘은 덱을 만든 뒤 다시 쓰이므로 인용하면 덱이 낡은 것으로 판정된다).
# 사전 명시·등록 이력: regulatory/tables/prespecification_register.csv(영문). 커밋 해시는 표 칸의 추적 문자열로만 쓴다.
# 이 슬라이드는 아토피 관련 말·파일을 쓰지 않는다(등록부의 해당 행은 표에서 뺐다; 부록 A5).

a6_field <- function(k, item = paste("GitHub reproducibility run:", k)) {
  x <- as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value)
  dderived(item, "repro/repro_github_run.csv", sprintf("field=='%s'", k), x, x)
}
a6_rx <- function(k, rx, item, nobreak_hyphen = FALSE) {
  x <- as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value); m <- regmatches(x, regexec(rx, x))[[1]]
  premise(length(m) >= 2, sprintf("repro_github_run.csv field %s matches %s", k, rx))
  dderived(item, "repro/repro_github_run.csv", sprintf("field=='%s' :: regex %s", k, rx), m[2], if (nobreak_hyphen) gsub("-", "\u2011", m[2]) else m[2])
}
# 등록부 한 행: 등록 시각·첫 결과를 "09-25 10:00 (521645a)" 꼴로(추적 문자열; 연도는 노트에 전체 날짜로)
a6_reg <- function(pat, col, item) {
  PR <- "regulatory/tables/prespecification_register.csv"; w <- sprintf("grepl('%s', Item)", pat)
  x <- as.character(row1(PR, w)[[col]]); m <- regmatches(x, regexec("^([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}) UTC \\(commit ([0-9a-f]{7})\\)$", x))[[1]]
  premise(length(m) == 3, sprintf("register date format [%s] %s", w, col))
  list(p = dderived(item, PR, sprintf("%s :: %s, month-day time (commit)", w, col), x, sprintf("%s (%s)", substr(m[2], 6, 16), m[3])), t = as.POSIXct(m[2], tz = "UTC", format = "%Y-%m-%d %H:%M"))
}

slide_A6 <- function() {
  CB <- "ci/ci_failures_before.csv"; CS <- "ci/ci_failure_summary.csv"; CA <- "ci/ci_runs_after_fix.csv"; CV <- "ci/ci_runs_v101.csv"
  RI <- "repro/repro_github_items.csv"; RC <- "repro/repro_check.csv"; TS <- "regulatory/tables/trace_summary.csv"; PR <- "regulatory/tables/prespecification_register.csv"
  deck_slide("A6", tag = "sim")
  L <- DK$txt$A6

  # ---- 자동 시험 기록: 전제 ----
  cb <- rows(CB); ca <- rows(CA, "workflow=='tests'"); cv <- rows(CV)
  runs <- rbind(cb[, .(run = run_number, conc = "failure")], ca[, .(run = run_number, conc = conclusion)], cv[, .(run = run_number, conc = conclusion)])
  premise(!anyDuplicated(runs$run) && identical(sort(runs$run), seq_len(max(runs$run))) && all(runs$conc %in% c("success", "failure", "cancelled")),
          "tests workflow runs 1 to the latest are each recorded once (success, failure or cancelled)")
  premise(all(nzchar(cb$failed_step)) && all(cb$category == "package install/build"), "runs before the fix failed in environment setup (package install/build)")
  fl <- cv[conclusion == "failure"]; premise(nrow(fl) > 0 && all(diff(fl$run_number) == 1) && all(grepl("^strict-skip rule", fl$note)), "v1.0.1 failures are consecutive and all due to the strict-skip rule")
  fxr <- cv[grepl("^fix of runs", note)]; premise(nrow(fxr) == 1 && fxr$conclusion == "success" && fxr$run_number == max(fl$run_number) + 1, "the run after the failures is the fix and succeeded")
  lastr <- cv[run_number == max(run_number)]; premise(lastr$conclusion == "success" && grepl("v1.0.1 package commit", lastr$note), "latest run succeeded on the v1.0.1 package commit")
  f1 <- ca[run_number == min(run_number)]; premise(f1$conclusion == "success" && grepl("first run after the fix", f1$note), "first run after the fix succeeded")
  premise(max(cb$run_number) + 1 == f1$run_number && max(ca$run_number) + 1 == min(cv$run_number), "the three CI files are contiguous")

  f <- list(
    r1 = drange(CB, "TRUE", "run_number", 0, "", "runs before the fix, run numbers"),
    a = dint(CS, "cause=='no renv.lock in commit'", "n_runs", "runs failed: no renv.lock in commit"),
    b = dint(CS, "cause=='xml2 missing for testthat::JunitReporter'", "n_runs", "runs failed: xml2 missing for JunitReporter"),
    fx1 = dint(CA, sprintf("workflow=='tests' & run_number==%d", f1$run_number), "run_number", "first run after the fix"),
    last = dint(CV, sprintf("run_number==%d", lastr$run_number), "run_number", "latest tests run (v1.0.1 package commit)"),
    rf = drange(CV, "conclusion=='failure'", "run_number", 0, "", "failed runs after the fix (strict-skip rule), run numbers"),
    fx2 = dint(CV, "grepl('^fix of runs', note)", "run_number", "run that fixed the strict-skip failures"))
  cnt <- function(k) { a_ <- nrow(ca[conclusion == k]); v_ <- nrow(cv[conclusion == k])
    dderived(sprintf("tests runs after the fix with conclusion %s (after_fix + v101)", k), CV, sprintf("count(conclusion=='%s') in ci_runs_after_fix.csv (workflow tests) + ci_runs_v101.csv", k), c(a_, v_), as.character(a_ + v_)) }
  dsrc("tests runs after the fix (with ci_runs_v101.csv)", CA, "(table)")
  f$ns <- cnt("success"); f$nc <- cnt("cancelled"); f$nf <- cnt("failure")

  # ---- 재현 워크플로 ----
  ri <- rows(RI); rc <- rows(RC)
  premise(all(ri$pass_github) && all(ri$identical_to_local_8sig) && all(rc$pass) && nrow(ri) == nrow(rc), "every reproducibility item passes locally and on the clean runner, identical to 8 significant digits")
  rgf <- function(k) as.character(row1("repro/repro_github_run.csv", sprintf("field=='%s'", k))$value)
  premise(as.numeric(rgf("n_items")) == nrow(ri) && rgf("n_pass_github") == rgf("n_items") && rgf("n_identical_to_local_8sig") == rgf("n_items") && rgf("conclusion") == "success",
          "clean-runner run: all items pass and are identical to the local values")
  g <- list(
    pg = sprintf("%s/%s", a6_field("n_pass_github"), a6_field("n_items")),
    ni = a6_field("n_items", "GitHub reproducibility run: items (all pass, title)"),
    nid = a6_field("n_identical_to_local_8sig"),
    sig = dderived("significant digits in column name identical_to_local_8sig", RI, "column name identical_to_local_8sig", 8, "8"),
    nl = dcount(RC, "pass==TRUE", "reproducibility items passing locally"),
    run = a6_field("run_number"),
    os = a6_rx("runner", "^(ubuntu-[0-9.]+) ", "GitHub runner image", nobreak_hyphen = TRUE),
    rv = a6_rx("r_version", "^R version ([0-9.]+) ", "R version on the clean runner"),
    nt = a6_rx("tests", "^([0-9]+) tests", "automated tests in the clean-environment run"),
    nft = a6_rx("tests", "^[0-9]+ tests, ([0-9]+) failed", "failed tests in the clean-environment run"),
    nsk = a6_rx("tests", " ([0-9]+) skipped \\(", "skipped tests in the clean-environment run"))

  # ---- 추적 규모: trace_summary.csv가 없으면(패키지 생성 스크립트를 다시 돌리기 전) 미산출로 표시한다(수치를 만들지 않는다) ----
  has_ts <- file.exists(proj_path(TS))
  if (has_ts) {
    ts <- rows(TS); premise(setequal(ts$document, c("MS_report", "FDA_questions", "SAP_text_proposals")), "trace summary covers the three regulatory documents")
    h <- list(tot = dderived("printed values in the three regulatory documents", TS, "sum(printed_values)", sum(ts$printed_values), fint(sum(ts$printed_values))),
              ms = dint(TS, "document=='MS_report'", "printed_values", "printed values, M&S report"),
              fq = dint(TS, "document=='FDA_questions'", "printed_values", "printed values, FDA questions"),
              sap = dint(TS, "document=='SAP_text_proposals'", "printed_values", "printed values, SAP text proposals"))
    hn <- list(sf_ms = dint(TS, "document=='MS_report'", "source_files", "source files, M&S report"),
               tb_ms = dint(TS, "document=='MS_report'", "tables", "table rows, M&S report"),
               fg_ms = dint(TS, "document=='MS_report'", "figures", "figure rows, M&S report"))
    tr_val <- h$tot; tr_lab <- tx("A6.card_trace", h); tr_note <- tx("A6.notes_trace", c(h, hn))
  } else {
    tr_val <- L$na; tr_lab <- tx("A6.card_trace_na"); tr_note <- tx("A6.notes_trace_na")
  }
  # 제목 전제: 엄격 건너뜀 수정 뒤 실행은 모두 성공
  premise(all(runs[run >= fxr$run_number, conc] == "success"), "every tests run from the strict-skip fix to the latest succeeded")

  deck_kicker(tx("A6.kicker")); deck_title(tx("A6.title", list(ni = g$ni, fx2 = f$fx2)))

  # ---- 왼쪽 위: 자동 시험 실행 띠 그림 ----
  XL <- GEO$ML; WL <- 7.15
  deck_text(tx("A6.ci_label"), c(XL, GEO$BODY_TOP, WL, 0.4), size = 16, bold = TRUE, color = PAL$ink2, label = "label_ci", gap_pt = 0)
  F <- L$fig
  CL <- c(success = F$success, failure = F$failure, cancelled = F$cancelled)
  runs[, lab := factor(CL[conc], levels = CL)]
  n1 <- max(cb$run_number); fr <- range(fl$run_number)
  ann <- data.table(x = c((1 + n1) / 2, mean(fr)), y = 1.62, lab = c(fill(F$before, list(a = nrow(cb[cause_key == "no renv.lock in commit"]), b = nrow(cb[cause_key == "xml2 missing for testthat::JunitReporter"]))), F$strict))
  premise(nrow(cb[cause_key == "no renv.lock in commit"]) == as.numeric(f$a) && nrow(cb[cause_key == "xml2 missing for testthat::JunitReporter"]) == as.numeric(f$b), "figure counts equal the failure summary")
  brk <- c(1, f1$run_number, fr[1], fxr$run_number, lastr$run_number)
  p <- ggplot(runs, aes(x = run, y = 1, fill = lab)) +
    geom_tile(width = 0.82, height = 0.7, colour = NA) +
    annotate("segment", x = 0.6, xend = n1 + 0.4, y = 1.43, yend = 1.43, colour = PAL$ink2, linewidth = 0.5) +
    annotate("segment", x = fr[1] - 0.4, xend = fr[2] + 0.4, y = 1.43, yend = 1.43, colour = PAL$ink2, linewidth = 0.5) +
    geom_text(data = ann, aes(x = x, y = y, label = lab), inherit.aes = FALSE, family = FONT, size = 4.1, colour = PAL$ink, vjust = 0) +
    scale_fill_manual(values = setNames(c(PAL$blue, PAL$orange, PAL$muted), CL), drop = FALSE) +
    scale_x_continuous(breaks = brk, labels = brk, expand = expansion(add = 0.8)) +
    scale_y_continuous(limits = c(0.6, 2.0), expand = expansion(0)) +
    labs(x = F$xlab, y = NULL) + theme_deck(13) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(), legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"), legend.key.size = grid::unit(0.9, "lines"))
  FY <- GEO$BODY_TOP + 0.42; FH <- 1.5
  deck_figure(p, "a6_ci_runs", c(XL, FY, WL, FH), src = c(CB, CA, CV))
  BY <- FY + FH + 0.04; BH <- 0.78
  deck_bullets(tx("A6.bullets", f), box = c(XL, BY, WL, BH), size = 16, gap_pt = 4)

  # ---- 왼쪽 아래: 사전 명시·등록 이력 표 ----
  REGROWS <- list(oc = "^Operating-characteristic design", s1 = "^Analysis-model re-judgement", s2 = "^Single study-LLOQ source", s3 = "^Sample-size table",
                  s4 = "^Criteria sets \\\\(iii\\\\) and \\\\(iv\\\\)", s6 = "^Trial-population analyses")
  rr <- lapply(names(REGROWS), function(k) { a_ <- a6_reg(REGROWS[[k]], "Date (evidence)", sprintf("register %s: registered (UTC, commit)", k))
    b_ <- a6_reg(REGROWS[[k]], "First results", sprintf("register %s: first results (UTC, commit)", k))
    premise(a_$t < b_$t, sprintf("register %s: registered before the first results", k)); c(a_$p, b_$p) })
  premise(identical(rows("oc/prereg.csv")$changed_since, FALSE), "oc_design.yaml unchanged since its pre-specification commit (results/oc/prereg.csv)")
  dsrc("oc_design.yaml unchanged since pre-specification", "oc/prereg.csv", "(table)")
  df <- data.frame(a = unlist(L$reg_rows[names(REGROWS)]), b = vapply(rr, `[`, "", 1), c = vapply(rr, `[`, "", 2), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("A6.table.head")
  TY <- BY + BH + 0.06
  deck_table(df, box = c(XL, TY, WL, GEO$BODY_BOTTOM - TY), widths = c(2.75, 2.2, 2.2), size = 12, label = "table_prereg")

  # ---- 오른쪽: 재현 워크플로 카드, 추적 규모 카드 ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  gap <- 0.14; CH1 <- 2.45
  deck_stat(g$pg, tx("A6.card_repro", g), c(XR, GEO$BODY_TOP, WR, CH1), color = PAL$blue, bg = PAL$tint_blue)
  deck_stat(tr_val, tr_lab, c(XR, GEO$BODY_TOP + CH1 + gap, WR, GEO$BODY_BOTTOM - GEO$BODY_TOP - CH1 - gap), color = PAL$ink2, bg = PAL$tint_grey)

  # ---- 노트 ----
  st <- rows(PR)
  cntp <- function(rx, it) dderived(it, PR, sprintf("count(grepl('%s', Status))", rx), sum(grepl(rx, st$Status)), as.character(sum(grepl(rx, st$Status))))
  n_ps <- sum(grepl("^pre-specified", st$Status)); n_pr <- sum(grepl("^pre-registered", st$Status)); n_ph <- sum(grepl("post hoc", st$Status))
  premise(sum(grepl("^pre-specified", st$Status) | grepl("^pre-registered", st$Status) | grepl("post hoc", st$Status)) == n_ps + n_pr + n_ph, "status groups do not overlap")
  deck_notes(tx("A6.notes", c(f, g, list(tr = tr_note,
    nreg = dcount(PR, "TRUE", "rows of the pre-specification register"),
    nps = cntp("^pre-specified", "register rows pre-specified"), npr = cntp("^pre-registered", "register rows pre-registered"), nph = cntp("post hoc", "register rows post hoc"),
    noth = dderived("register rows with other statuses (changes, amendments, decisions)", PR, "count of rows not matching ^pre-specified, ^pre-registered or post hoc", nrow(st) - n_ps - n_pr - n_ph, as.character(nrow(st) - n_ps - n_pr - n_ph)),
    js = a6_field("job_duration_s", "GitHub reproducibility run: job duration (s)"),
    rx = a6_rx("rxode2_version", "^([0-9.]+)$", "rxode2 version on the clean runner")))))
  deck_end()
}
