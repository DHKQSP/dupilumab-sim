# A2 부록: 선례. 요청된 세 건(PKM14161, Cohen 2022, 토실리주맙 MSB11456 정맥 시험) 중 프로젝트에 기록된 내용만 적는다. [문헌+모의]
# 제목은 MSB11456 기록의 상태(검색 발췌, 원문 미대조)를 앞에 밝힌다(premise로 상태 검사).
# PKM14161: 300 mg 검증 arm(Li 2020 Table 3)으로만 기록(평가변수·AUC0-inf 처리 기록 없음, 부분). Cohen 2022: 문헌(AUC0-inf 미산출).
# 토실리주맙 MSB11456 정맥 시험: config/literature_precedents.yaml(검색 발췌 두 건으로 확인, 원문은 분석 환경에서 접근 불가; D-062).

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
  premise(length(unique(vapply(Filter(function(z) identical(z$study, "PKM14161"), y$arm_checks_300mg$arms), function(z) z$tmax_median, 1))) == 1, "PKM14161: same published median tmax in both arms (card says 'both arms')")
  # 토실리주맙 MSB11456: 문헌 기록(config)의 상태가 '검색 발췌'인 동안은 카드와 노트에 원문 대조 필요를 적는다
  LP <- "literature_precedents.yaml"; mp <- .read(file.path("config", LP))$msb11456_iv
  premise(identical(mp$status, "search excerpt"), "MSB11456 record status is 'search excerpt' (card text says the full text was not compared)")
  premise(identical(mp$primary_endpoint, "AUC0-last") && grepl("AUC0-inf", mp$secondary_endpoints), "MSB11456: AUC0-last primary, AUC0-inf secondary (text)")
  premise(mp$n_test + mp$n_reference == mp$n_randomized, "MSB11456 arm sizes add up to the randomized number")
  rp <- function(k, lab) dcfg(LP, c("msb11456_iv", "ratio_pct", k), sprintf("MSB11456 IV study, geometric LS mean ratio (%%) and 90%% CI, %s", lab),
                              function(x) sprintf("%s%% (%s)", fnum(x[1], 2), rng_fmt(x[2], x[3], 2)))

  lsd_lab <- rows("fallback/sample_size_logsd.csv", "startsWith(variant, 'Cohen')")$variant
  premise(length(lsd_lab) == 1 && grepl("90% CI", lsd_lab, fixed = TRUE), "the Cohen 2022 GMR interval is recorded as a 90% CI")
  dsrc("confidence level of the Cohen 2022 GMR interval (label)", "fallback/sample_size_logsd.csv", "(table)")
  f <- list(
    at = a2_cfg("arms", list(study = "PKM14161", arm = "test"), "auclast_mean", "PKM14161 test arm, published mean AUC0-last (Li 2020 Table 3)"),
    ar = a2_cfg("arms", list(study = "PKM14161", arm = "reference"), "auclast_mean", "PKM14161 reference arm, published mean AUC0-last (Li 2020 Table 3)"),
    tm = a2_cfg("arms", list(study = "PKM14161", arm = "test"), "tmax_median", "PKM14161 published median tmax (day), same in both arms", function(x) fnum(x, 1)),
    d300 = f_dose(),
    q16t = dv(F16, "study=='PKM14161' & arm=='test'", "ratio", 2, "", "PKM14161 test, sim/obs AUC0-last, 2016"), q16r = dv(F16, "study=='PKM14161' & arm=='reference'", "ratio", 2, "", "PKM14161 reference, sim/obs AUC0-last, 2016"),
    q20t = dv(F20, "study=='PKM14161' & arm=='test'", "ratio", 2, "", "PKM14161 test, sim/obs AUC0-last, 2020"), q20r = dv(F20, "study=='PKM14161' & arm=='reference'", "ratio", 2, "", "PKM14161 reference, sim/obs AUC0-last, 2020"),
    tol = dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "AUC0-last mean tolerance of the validation (%)", num_fmt(0)),
    cd = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "dose_mg", "Cohen 2022 dose (mg)"),
    cc = a2_parse(Q16, "id=='Cohen2022_200mg_AI'", "presentation", "of ([0-9.]+) mg/mL", "Cohen 2022 concentration (mg/mL)"),
    n1 = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "n_subj", "Cohen 2022 autoinjector n"),
    n2 = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "n_subj", "Cohen 2022 prefilled syringe n"),
    ca1 = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "auclast_mean", "Cohen 2022 autoinjector, published mean AUC0-last (Table 2)"),
    ca2 = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "auclast_mean", "Cohen 2022 prefilled syringe, published mean AUC0-last (Table 2)"),
    gmr = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "gmr_auclast", "Cohen 2022 AUC0-last GMR", function(x) fnum(x, 2)),
    ci = a2_cfg("datasets", list(id = "Cohen2022_200mg_PFSS"), "gmr_ci", "Cohen 2022 AUC0-last GMR 90% CI", function(x) rng_fmt(x[1], x[2], 2)),
    lsd = dv(LSD, "gate_role=='external'", "obs_log_sd", 2, "", "Cohen 2022 implied log SD of AUC0-last"),
    tn = dcfg(LP, c("msb11456_iv", "n_randomized"), "MSB11456 IV study, randomized healthy adults"),
    tt = dcfg(LP, c("msb11456_iv", "n_test"), "MSB11456 IV study, test arm n"),
    trf = dcfg(LP, c("msb11456_iv", "n_reference"), "MSB11456 IV study, reference arm n"),
    tdose = dcfg(LP, c("msb11456_iv", "dose_mg_per_kg"), "MSB11456 IV study, dose (mg/kg)"),
    tday = dcfg(LP, c("msb11456_iv", "sampling_last_day"), "MSB11456 IV study, last sampling day"),
    tinf = dcfg(LP, c("msb11456_iv", "infusion_h"), "MSB11456 IV study, infusion duration (h)"),
    tal = rp("auc0_last", "AUC0-last"), tai = rp("auc0_inf", "AUC0-inf"), tcm = rp("cmax", "Cmax"))
  f$q16 <- sprintf("%s / %s", f$q16t, f$q16r); f$q20 <- sprintf("%s / %s", f$q20t, f$q20r)   # 노트: 시험 / 대조

  # ---- 카드 세 개(제목 순서: MSB11456, Cohen 2022, PKM14161): 상태 표지 + 기록 내용 ------------------------------------------------------
  # 폭은 글 양에 맞춰 나누고(세 카드의 추정 높이가 비슷하게), 카드 높이는 가장 긴 카드의 추정 높이에 맞춘다(빈 카드 공간을 줄임). 결론 문장은 카드 바로 아래에 둔다.
  cards <- list(toci = list(bg = PAL$tint_orange, chip = PAL$orange, w = 1.22), cohen = list(bg = PAL$tint_blue, chip = PAL$blue, w = 0.79), pkm = list(bg = PAL$tint_grey, chip = PAL$ink2, w = 0.99))
  gap <- 0.25; ws <- vapply(cards, function(z) z$w, 1); ws <- ws / sum(ws) * (GEO$CW - 2 * gap); yc <- GEO$BODY_TOP + 0.02; chip_h <- 0.42
  body <- lapply(names(cards), function(nm) tx(sprintf("A2.%s.body", nm), f)); names(body) <- names(cards)
  ch <- max(mapply(function(b, w) est_height(b, w, 16, gap_pt = 7, card = TRUE), body, ws)) + 0.08
  for (k in seq_along(cards)) {
    nm <- names(cards)[k]; x <- GEO$ML + sum(ws[seq_len(k - 1)]) + (k - 1) * gap
    deck_text(tx(sprintf("A2.%s.chip", nm)), c(x, yc, ws[k], chip_h), size = 16, bold = TRUE, color = "#ffffff", bg = cards[[nm]]$chip, geom = "roundRect", label = sprintf("chip_%s", nm), gap_pt = 0)
    deck_text(body[[nm]], c(x, yc + chip_h + 0.08, ws[k], ch), size = 16, bg = cards[[nm]]$bg, geom = "roundRect", label = sprintf("card_%s", nm), gap_pt = 7)
  }
  yb <- yc + chip_h + 0.08 + ch + 0.14
  deck_text(tx("A2.bottom"), c(GEO$ML, yb, GEO$CW, 0.45), size = 16, bold = TRUE, label = "takeaway", gap_pt = 0)
  deck_notes(tx("A2.notes", c(f, list(
    wc = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "weight_mean", "Cohen 2022 mean body weight (kg)", function(x) fnum(x, 1)),
    tmc = a2_cfg("datasets", list(id = "Cohen2022_200mg_AI"), "tmax_median", "Cohen 2022 median tmax (day)", function(x) fnum(x, 2)),
    ca16 = drange(Q16, "gate_role=='external' & grepl('Cohen', id)", "AUClast_ratio", 2, "", "Cohen 2022 sim/obs AUC0-last, 2016"),
    doi = dcfg(LP, c("msb11456_iv", "doi"), "MSB11456 IV study, DOI"), yr = dcfg(LP, c("msb11456_iv", "year"), "MSB11456 IV study, publication year"),
    pon = dcfg(LP, c("msb11456_iv", "published_online"), "MSB11456 IV study, published online"), pmid = dcfg(LP, c("msb11456_iv", "pmid"), "MSB11456 IV study, PubMed ID")))))
  deck_end()
}
