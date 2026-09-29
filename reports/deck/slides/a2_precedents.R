# A2 부록: 선례. 요청된 세 건(PKM14161, Cohen 2022, 토실리주맙 MSB11456 정맥 시험) 중 프로젝트에 기록된 내용만 적는다. [문헌+모의]
# PKM14161: 300 mg 검증 arm(Li 2020 Table 3)으로만 기록(평가변수·AUC0-inf 처리 기록 없음, 부분). Cohen 2022: 문헌(AUC0-inf 미산출).
# 토실리주맙 MSB11456 정맥 시험: 저장소 어디에도 기록이 없다(보고서, config, 문헌 발췌표, SPEC, DECISIONS 검색 결과 없음) → "프로젝트 기록 없음(확인 필요)".

a2_cfg <- function(section, match, key, item, fmt = function(x) format(x)) {
  y <- .read("config/design_clot2021.yaml"); lst <- if (section == "arms") y$arm_checks_300mg$arms else y$datasets
  hit <- Filter(function(z) all(vapply(names(match), function(k) identical(z[[k]], match[[k]]), TRUE)), lst)
  cond <- paste(sprintf("%s=='%s'", names(match), unlist(match)), collapse = " & ")
  premise(length(hit) == 1 && !is.null(hit[[1]][[key]]), sprintf("design_clot2021 %s [%s] %s", section, cond, key))
  x <- unlist(hit[[1]][[key]])
  dderived(item, "config/design_clot2021.yaml", sprintf("%s[%s].%s", if (section == "arms") "arm_checks_300mg.arms" else "datasets", cond, key), x, fmt(x))
}
a2_parse <- function(rel, where, col, rx, item) {
  r <- row1(rel, where); m <- regmatches(r[[col]], regexec(rx, r[[col]]))[[1]]; premise(length(m) == 2, sprintf("regex %s in %s [%s] %s", rx, rel, where, col))
  dderived(item, rel, sprintf("%s :: %s, regex %s", where, col, rx), as.numeric(m[2]), m[2])
}

slide_A2 <- function() {
  F16 <- "step1/step1f_300mg_arm_check.csv"; F20 <- "results/step1_k2020/step1f_300mg_arm_check.csv"; Q16 <- "step1/step1b_quant_gate.csv"; LSD <- "step1/step1d_cohen2022_logsd.csv"
  deck_slide("A2", tag = "litsim")
  deck_kicker(tx("A2.kicker")); deck_title(tx("A2.title"))
  # 전제: PKM14161은 2020 모델 개발 자료이고 2016 모델에는 외부 자료이며, 2016 모델 검증에서 두 arm 모두 통과
  y <- .read("config/design_clot2021.yaml")
  premise("PKM14161" %in% unlist(y$model_development_data$k2020) && !("PKM14161" %in% unlist(y$model_development_data$k2016)), "PKM14161 is in the 2020 model development data only")
  premise(all(rows(F16, "study=='PKM14161'")$pass_mean) && all(rows(F16, "study=='PKM14161'")$dev == "external"), "PKM14161 arms pass the 2016 model check as external data")
  # 토실리주맙 MSB11456: 근거 문서(보고서, SPEC, config, 문헌 결과)에 없고, 결정 기록에는 '기록이 없는 항목'으로만 나온다
  rx <- "MSB11456|tocilizumab|\ud1a0\uc2e4\ub9ac\uc8fc\ub9d9"; rd <- function(f) readLines(f, warn = FALSE, encoding = "UTF-8")
  ev_files <- c(proj_path("regulatory", "src", "MS_report.Rmd"), proj_path("regulatory", "src", "FDA_questions.Rmd"), proj_path("SPEC.md"),
                list.files(proj_path("config"), pattern = "\\.ya?ml$", full.names = TRUE), list.files(proj_path("results", "literature"), full.names = TRUE))
  premise(!any(grepl(rx, unlist(lapply(ev_files, rd)), ignore.case = TRUE)), "no record of the tocilizumab MSB11456 IV study in the report, SPEC, config or literature results")
  dl <- grep(rx, rd(proj_path("DECISIONS.md")), value = TRUE, ignore.case = TRUE)
  premise(all(grepl("\uae30\ub85d\uc774 \uc5c6", dl)), "DECISIONS mentions the tocilizumab study only as an item without a project record")

  lsd_lab <- rows("fallback/sample_size_logsd.csv", "startsWith(variant, 'Cohen')")$variant
  premise(length(lsd_lab) == 1 && grepl("90% CI", lsd_lab, fixed = TRUE), "the Cohen 2022 GMR interval is recorded as a 90% CI")
  dsrc("confidence level of the Cohen 2022 GMR interval (label)", "fallback/sample_size_logsd.csv", "(table)")
  f <- list(
    at = a2_cfg("arms", list(study = "PKM14161", arm = "test"), "auclast_mean", "PKM14161 test arm, published mean AUC0-last (Li 2020 Table 3)"),
    ar = a2_cfg("arms", list(study = "PKM14161", arm = "reference"), "auclast_mean", "PKM14161 reference arm, published mean AUC0-last (Li 2020 Table 3)"),
    tm = a2_cfg("arms", list(study = "PKM14161", arm = "test"), "tmax_median", "PKM14161 test arm, published median tmax (day)", function(x) fnum(x, 1)),
    d300 = f_dose(),
    q16 = sprintf("%s / %s", dv(F16, "study=='PKM14161' & arm=='test'", "ratio", 2, "", "PKM14161 test, sim/obs AUC0-last, 2016"), dv(F16, "study=='PKM14161' & arm=='reference'", "ratio", 2, "", "PKM14161 reference, sim/obs AUC0-last, 2016")),
    q20 = sprintf("%s / %s", dv(F20, "study=='PKM14161' & arm=='test'", "ratio", 2, "", "PKM14161 test, sim/obs AUC0-last, 2020"), dv(F20, "study=='PKM14161' & arm=='reference'", "ratio", 2, "", "PKM14161 reference, sim/obs AUC0-last, 2020")),
    tol = dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "AUC0-last mean tolerance of the validation (%)", num_fmt(0)),
    cd = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "dose_mg", "Cohen 2022 dose (mg)"),
    cc = a2_parse(Q16, "id=='Cohen2022_200mg_AI'", "presentation", "of ([0-9.]+) mg/mL", "Cohen 2022 concentration (mg/mL)"),
    n1 = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "n_subj", "Cohen 2022 autoinjector n"),
    n2 = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "n_subj", "Cohen 2022 prefilled syringe n"),
    gmr = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "gmr_auclast", "Cohen 2022 AUC0-last GMR", function(x) fnum(x, 2)),
    ci = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "gmr_ci", "Cohen 2022 AUC0-last GMR 90% CI", function(x) rng_fmt(x[1], x[2], 2)),
    lsd = dv(LSD, "gate_role=='external'", "obs_log_sd", 2, "", "Cohen 2022 implied log SD of AUC0-last"))

  # ---- 카드 세 개: 상태 표지 + 기록 내용 ---------------------------------------------------------------------------------------------
  gap <- 0.25; cw <- (GEO$CW - 2 * gap) / 3; yc <- GEO$BODY_TOP + 0.05; chip_h <- 0.42; ch <- 3.6
  cards <- list(pkm = list(bg = PAL$tint_grey, chip = PAL$ink2), cohen = list(bg = PAL$tint_blue, chip = PAL$blue), toci = list(bg = PAL$tint_orange, chip = PAL$orange))
  for (k in seq_along(cards)) {
    nm <- names(cards)[k]; x <- GEO$ML + (k - 1) * (cw + gap)
    deck_text(tx(sprintf("A2.%s.chip", nm)), c(x, yc, cw, chip_h), size = 16, bold = TRUE, color = "#ffffff", bg = cards[[nm]]$chip, geom = "roundRect", label = sprintf("chip_%s", nm), gap_pt = 0)
    deck_text(tx(sprintf("A2.%s.body", nm), f), c(x, yc + chip_h + 0.08, cw, ch), size = 16, bg = cards[[nm]]$bg, geom = "roundRect", label = sprintf("card_%s", nm), gap_pt = 6)
  }
  yb <- yc + chip_h + 0.08 + ch + 0.15
  deck_text(tx("A2.bottom"), c(GEO$ML, yb, GEO$CW, GEO$BODY_BOTTOM - yb), size = 16, bold = TRUE, label = "takeaway")
  deck_notes(tx("A2.notes", c(f, list(
    wc = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "weight_mean", "Cohen 2022 mean body weight (kg)", function(x) fnum(x, 1)),
    tmc = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "tmax_median", "Cohen 2022 median tmax (day)", function(x) fnum(x, 2)),
    ca16 = drange(Q16, "gate_role=='external' & grepl('Cohen', id)", "AUClast_ratio", 2, "", "Cohen 2022 sim/obs AUC0-last, 2016")))))
  deck_end()
}
