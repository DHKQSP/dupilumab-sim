# A12 별첨: 한계와 검증 상태(결과보고 덱 S25 한계 + S26 남은 결정을 핵심 덱 용어로 한 장에). [문헌+모의]
# 표 7행(첫 열 = 상태, 굵게): 모델 가정 2개(정량한계 아래 곡선 모양·Km 고정 → 민감도 분석으로만 확인(S8, 별첨 A8); 200 mg 제형 빠른 흡수 → 재현 안 됨(별첨 A1)),
# 검증 4개(실제 Phoenix WinNonlin 대조 미실시; 담당자 독립 QC 대기; 문헌·지침 인용은 검색 발췌로 원문 미대조(config/literature_core_deck.yaml,
# config/literature_precedents.yaml의 status); 핵심 덱 새 요약(scripts/63)은 원래 시드 규칙으로 다시 만들어 커밋된 결과와 일치(results/core_deck/provenance.csv)),
# 의뢰자 결정 대기 1개(LLOQ, 잠정 설계값, 주분석 모형, 목표 검정력; 결과보고 덱 S26). 상태는 config·결과 파일의 status 값과 전제(premise)로 확인한다.
# 캡션: 별첨 A13(규제 업무 포지션 페이퍼 7절 대응표)은 싣지 않음. 포지션 페이퍼와 직전 지시의 결과가 저장소에 없다(DECISIONS D-063).

# 문자열 열에서 정규식으로 숫자 하나를 읽는다(행들이 같은 값이어야 한다)
a12_parse <- function(rel, where, col, rx, item) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where))
  x <- unique(as.numeric(vapply(r[[col]], function(s) { m <- regmatches(s, regexec(rx, s))[[1]]; if (length(m) < 2) NA_character_ else m[2] }, "")))
  premise(length(x) == 1 && is.finite(x), sprintf("one number matching %s in %s [%s] %s", rx, rel, where, col))
  dderived(item, rel, sprintf("%s :: %s, regex %s", where, col, rx), x, format(x))
}
# 곡선 모양 민감도 변형의 배율(config/scenarios.yaml): 두 값이면 "a·×b"
a12_mult <- function(variants, par, item) {
  y <- .read("config/scenarios.yaml")$sensitivity_variants
  x <- sort(vapply(variants, function(v) as.numeric(y[[v]]$theta_multipliers[[par]]), 0)); premise(length(x) == 2 && all(is.finite(x)), paste("two multipliers of", par))
  dderived(item, "config/scenarios.yaml", sprintf("sensitivity_variants [%s] :: theta_multipliers.%s, both values", paste(variants, collapse = ", "), par), x,
           sprintf("%s·×%s", format(x[1]), format(x[2])))
}
# config 문자열 값(상태 등)을 그대로 적는다(추적 행에 남김)
a12_txt <- function(file, path, item) dcfg(file, path, item, function(x) as.character(x))

slide_A12 <- function() {
  CV <- "core_deck/coverage_by_case.csv"; PV <- "core_deck/provenance.csv"; INV <- "oc/inversion_all.csv"
  Q16 <- "step1/step1b_quant_gate.csv"; Q20 <- "results/step1_k2020/step1b_quant_gate.csv"
  EV <- "nca_engine/engine_validation_summary.csv"; QC <- "regulatory/tables/verification_qc.csv"
  LC <- "config/literature_core_deck.yaml"; LP <- "config/literature_precedents.yaml"
  LIf <- "lloq/lloq_individual_table.csv"; TPf <- "sample_size/ss_table_power.csv"; NNf <- "sample_size/ss_table_n_needed.csv"
  deck_slide("A12", tag = "litsim")
  L <- DK$txt$A12

  # ---- 전제: 모델 가정 --------------------------------------------------------------------------------------------------------------------
  P16 <- .read("config/params_typical.yaml"); P20 <- .read("config/params_k2020_model1.yaml"); PVv <- .read("config/params_variability.yaml")
  premise(isTRUE(P16$theta$Km$fixed) && isTRUE(P20$theta$Km$fixed) && P16$theta$Km$value == P20$theta$Km$value, "Km fixed at the same value in both models")
  premise(PVv$iiv_omega2$Km$omega2 == 0 && grepl("Km:\\s*\\{sd: 0, omega2: 0\\}", paste(readLines(proj_path("config", "params_k2020_model1.yaml"), encoding = "UTF-8"), collapse = "\n")),
          "no between-subject variability of Km in either model")
  CUR <- c(vmax080 = "vmax080_both", vmax125 = "vmax125_both", km05 = "km05_both", km10 = "km10_both")
  cur <- rows(CV, "case %in% c('vmax080','vmax125','km05','km10')")
  premise(nrow(cur) == 4 && all(mapply(function(cs, src) grepl(paste0("variant ", CUR[[cs]], " "), src, fixed = TRUE), cur$case, cur$source)), "curve-shape cases are the scripts/21 variants (both arms)")
  lim <- as.numeric(.read("config/trial_design.yaml")$be$limits); kinv <- rows(INV, "mechanism=='Km'")
  kx <- c(kinv[is.finite(end_auc_ratio), end_auc_ratio], kinv[reachable == TRUE, auc_ratio])
  premise(length(kx) > 0 && min(kx) > lim[1] && max(kx) < lim[2], "Km differences in the test arm never reach the equivalence limits (notes)")
  g16 <- rows(Q16, "gate_role=='gate'"); g20 <- rows(Q20, "gate_role=='gate'"); e16 <- rows(Q16, "gate_role=='external'")
  dose <- as.numeric(.read("config/trial_design.yaml")$dose_mg); tol <- as.numeric(.read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct)
  premise(all(abs(g16[dose_mg == dose, AUClast_ratio] - 1) <= tol / 100) && all(abs(g20[dose_mg == dose, AUClast_ratio] - 1) <= tol / 100), "study-dose AUClast within the exposure gate in both models (row 2: exposure reproduced)")
  premise(all(e16$dose_mg == 200) && all(!is.na(e16$tmax_obs_median)) && all(e16$tmax_sim_median > e16$tmax_obs_median), "200 mg data sets: simulated tmax later than observed (row 2: fast absorption not reproduced)")
  e20 <- rows(Q20, "gate_role=='external'")
  premise(all(e20$tmax_sim_median > e20$tmax_obs_median) && identical(sort(unique(e20$tmax_sim_median)), sort(unique(e16$tmax_sim_median))),
          "same in the 2020 model, with the same simulated tmax (row 2 gives one simulated value for both models)")

  # ---- 전제: 검증 상태 ---------------------------------------------------------------------------------------------------------------------
  ev <- rows(EV); premise(all(ev$pass) && setequal(unique(ev$comparison), c("this engine vs NonCompart", "PKNCA vs NonCompart", "this engine vs PKNCA")), "the NCA engine was compared with NonCompart and PKNCA only")
  premise(!any(grepl("phoenix|winnonlin", list.files(proj_path("results", "nca_engine"), recursive = TRUE), ignore.case = TRUE)), "no Phoenix output in results/nca_engine")
  qc <- row1(QC, "Activity=='Independent human QC of this report'"); premise(startsWith(qc$Result, "PENDING"), "independent human QC pending")
  lc <- .read(LC); lp <- .read(LP)
  st <- c(vapply(lc, function(e) e$status, ""), vapply(lp, function(e) e$status, ""))
  premise(length(lc) == 2 && length(lp) == 2 && all(st == "search excerpt"), "all four cited sources are search excerpts (status 'search excerpt')")
  premise(isFALSE(lc$fda_bla761055_clinpharm$page_verified), "FDA review page not verified")
  pv <- rows(PV)
  premise(nrow(rows(PV, "equal == FALSE")) == 0 && all(startsWith(pv[is.na(equal), check], "input ")) && nrow(rows(PV, "equal == TRUE")) > 0,
          "every identity check in provenance.csv is equal; the other rows are input hashes")
  # 직전 지시(포지션 페이퍼 연계)의 결과는 저장소에 없다(D-063): 결과·문서 폴더에 포지션 페이퍼 파일이 없다
  pp <- list.files(proj_path("results"), pattern = "position[_ -]?paper|포지션", recursive = TRUE, ignore.case = TRUE)
  pp2 <- list.files(proj_path("regulatory"), pattern = "position[_ -]?paper|포지션", recursive = TRUE, ignore.case = TRUE)
  premise(length(pp) == 0 && length(pp2) == 0 && any(grepl("A13", readLines(proj_path("DECISIONS.md"), warn = FALSE, encoding = "UTF-8"), fixed = TRUE)),
          "no position paper and no results of the earlier directive in the repository (appendix A13 omitted, DECISIONS D-063)")

  # ---- 전제: 의뢰자 결정 대기 ----------------------------------------------------------------------------------------------------------------
  TD <- .read("config/trial_design.yaml"); AS <- .read("config/assay.yaml")
  premise(AS$lloq_mg_L$status == "assumption", "assay LLOQ is a placeholder (status assumption)")
  premise(TD$day1_postdose_time$status == "assumption" && TD$sampling_window_days$status == "assumption" && TD$stratification$split_kg$status == "assumption" &&
          TD$weight$base$status == "assumption" && TD$weight$sex_ratio_male$status == "assumption", "design values are placeholders (status assumption)")
  cv_v <- as.numeric(.read("config/prereg_20260926.yaml")$section3$base_cv_pct); n_v <- as.numeric(TD$n_per_arm)
  wp <- function(am) sprintf("input_model=='k2016' & cv==%s & gmr==0.95 & n==%s & analysis_model=='%s'", cv_v, n_v, am)
  wn <- function(am, tg) sprintf("cv==%s & gmr==0.95 & analysis_model=='%s' & target_pct==%s", cv_v, am, tg)

  # ---- 제목 ------------------------------------------------------------------------------------------------------------------------------
  km <- dcfg("params_typical.yaml", c("theta", "Km", "value"), "Km (mg/L), fixed in both models", function(x) format(x))
  y0 <- core_title(tx("A12.title", list(km = km)), tx("A12.kicker"))

  # ---- 표 --------------------------------------------------------------------------------------------------------------------------------
  p05c <- 100 * min(cur$p05)
  lz <- { a <- sum(ev$lz_points_identical); b <- sum(ev$lz_points_identical + ev$lz_points_mismatch)
    dderived("lambda-z windows identical over all comparisons", EV, "sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", c(a, b), sprintf("%s / %s", fint(a), fint(b))) }
  f <- list(
    km = km, v = a12_mult(c("vmax080_both", "vmax125_both"), "Vmax", "Vmax multipliers, both arms, curve-shape cases"),
    k = a12_mult(c("km05_both", "km10_both"), "Km", "Km multipliers, both arms, curve-shape cases"),
    p05 = dderived("smallest 5th percentile of AUClast/AUCinf over the curve-shape cases, rounded down", CV, "case %in% c('vmax080','vmax125','km05','km10') :: floor(min(p05) x 1000) / 10", p05c, paste0(fnum(fl(p05c, 1), 1), "%")),
    d200 = drange(Q16, "gate_role=='external'", "dose_mg", 0, "", "dose of the 175 mg/mL external data sets (mg)"),
    c175 = a12_parse(Q16, "gate_role=='external'", "presentation", "([0-9.]+) mg/mL", "concentration of the 200 mg presentation (mg/mL)"),
    to = drange(Q16, "gate_role=='external' & !is.na(tmax_obs_median)", "tmax_obs_median", 1, "", "observed median tmax, 200 mg (175 mg/mL) data sets (day)"),
    ts = drange(Q16, "gate_role=='external'", "tmax_sim_median", 0, "", "simulated median tmax, 200 mg data sets, 2016 (day)"),
    d300 = f_dose(), lz = lz,
    qc = local({ .record("independent human QC status", QC, "Activity=='Independent human QC of this report' :: Result", qc$Result, qc$Result); qc$Result }),
    page = dcfg("literature_core_deck.yaml", c("fda_bla761055_clinpharm", "page"), "FDA BLA 761055 review page given in the directive (not verified)", num_fmt(0)),
    e80 = dcfg("literature_core_deck.yaml", c("ema_be_guideline_80pct", "coverage_min_pct"), "EMA BE guideline: minimum coverage of AUC0-t (%)", num_fmt(0)),
    nlit = dderived("cited sources with status 'search excerpt' (core deck and results deck literature files)", "config/literature_core_deck.yaml", "count of entries with status == 'search excerpt' in literature_core_deck.yaml and literature_precedents.yaml", length(st), as.character(length(st))),
    n34 = dcount(PV, "equal == TRUE", "identity checks of the regenerated core-deck summaries (all equal)"),
    lloq = f_lloq(), n_arm = f_n_arm(),
    tg90 = dint(NNf, wn("M1", 90), "target_pct", "target power (%)"), n90 = dint(NNf, wn("M1", 90), "n_evaluable_per_arm", "n per arm for the higher target, M1"),
    tg85 = dint(NNf, wn("M1", 85), "target_pct", "target power (%)"), n85 = dint(NNf, wn("M1", 85), "n_evaluable_per_arm", "n per arm for the lower target, M1"))
  dsrc("literature status (results deck citations)", LP, "(table)")
  H <- L$table
  df <- data.frame(a = tx("A12.table.status", f), b = tx("A12.table.item", f), c = tx("A12.table.detail", f), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- unlist(H$head)
  cap <- tx("A12.caption")
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  deck_table(df, box = c(GEO$ML, y0, GEO$CW, capy - 0.16 - y0), widths = c(1.95, 3.1, 7.18), size = 14, align_num = FALSE, label = "table_limits",
             highlight = 3:5, highlight_fill = PAL$tint_grey)   # 아직 하지 않은 검증(제목)
  deck_src_first(c(PV, QC, LC))

  # ---- 노트 ------------------------------------------------------------------------------------------------------------------------------
  deck_notes(tx("A12.notes", c(f, list(
    kmed = drange("curve_shape/curve_shape_B0.csv", "variant %in% c('km05_both','km10_both')", "extrap_true_median", 2, "%", "Km x0.5 and x10 (both arms), median true extrapolation"),
    im = local({ r <- rows(INV, "mechanism=='Km'"); x <- range(r$end_multiplier)
      dderived("Km multipliers searched in the test arm (range ends)", INV, "mechanism=='Km' :: range(end_multiplier)", x, sprintf("%s~×%s", format(x[1]), format(x[2]))) }),
    kr = dderived("true AUCinf ratio range over the Km search (range ends and reached rows)", INV, "mechanism=='Km' :: end_auc_ratio, auc_ratio", range(kx), sprintf("%s~%s", fnum(min(kx), 3), fnum(max(kx), 3))),
    lims = dderived("equivalence limits (%), printed without decimals", "config/trial_design.yaml", "be.limits :: x 100, 0 decimals", 100 * lim, rng_fmt(100 * lim[1], 100 * lim[2], 0)),
    ce = dcount(CV, "case %in% c('vmax080','vmax125','km05','km10') & pct_lt80 == 0", "curve-shape cases without any subject below 80%"),
    a16 = drange(Q16, "gate_role=='external'", "AUClast_ratio", 2, "", "AUClast sim/obs, 200 mg (175 mg/mL) data sets, 2016"),
    a20 = drange(Q20, "gate_role=='external'", "AUClast_ratio", 2, "", "AUClast sim/obs, 200 mg (175 mg/mL) data sets, 2020"),
    c150 = a12_parse(Q16, "gate_role=='gate' & dose_mg==300", "presentation", "([0-9.]+) mg/mL", "concentration of the study presentation (mg/mL)"),
    g300a = drange(Q16, sprintf("gate_role=='gate' & dose_mg==%s", dose), "AUClast_ratio", 2, "", "AUClast sim/obs, study-dose data sets, 2016"),
    g300b = drange(Q20, sprintf("gate_role=='gate' & dose_mg==%s", dose), "AUClast_ratio", 2, "", "AUClast sim/obs, study-dose data sets, 2020"),
    gate = dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "exposure gate: AUClast mean within +/- tolerance (%)", num_fmt(0)),
    gc16 = drange(Q16, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, study presentation, 2016"),
    gc20 = drange(Q20, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, study presentation, 2020"),
    mx = local({ x <- max(ev$max_rel_diff); e <- floor(log10(x)); m <- signif(x, 2) / 10^e
      sup <- chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(e))
      dderived("largest relative parameter difference", EV, "max(max_rel_diff), printed as mantissa x 10^exponent", x, sprintf("%s×10%s", format(m), sup)) }),
    s_fda = a12_txt("literature_core_deck.yaml", c("fda_bla761055_clinpharm", "status"), "status of the FDA BLA 761055 review citation"),
    s_ema = a12_txt("literature_core_deck.yaml", c("ema_be_guideline_80pct", "status"), "status of the EMA bioequivalence guideline citation"),
    s_msb = a12_txt("literature_precedents.yaml", c("msb11456_iv", "status"), "status of the MSB11456 study citation"),
    s_mab = a12_txt("literature_precedents.yaml", c("ema_2012_mab", "status"), "status of the EMA 2012 monoclonal antibody guideline citation"),
    rdate = a12_txt("literature_core_deck.yaml", c("fda_bla761055_clinpharm", "review_date"), "FDA BLA 761055 review date"),
    nin = dcount(PV, "startsWith(check, 'input ')", "input files hashed in provenance.csv"),
    nrep = dcount(PV, "equal == TRUE & grepl('subject', check)", "representative-subject checks in provenance.csv"),
    grid = f_lloq_grid(),
    cov = drange(LIf, "model=='k2016' & resid=='fixed'", "coverage_lt80_pct", 2, "%", "LLOQ grid, share below 80% AUClast/AUCinf, 2016 model"),
    thr = dderived("AUClast/AUCinf threshold in column name coverage_lt80_pct (percent)", LIf, "column name coverage_lt80_pct", 80, "80"),
    t1 = dcfg("trial_design.yaml", c("day1_postdose_time", "value"), "Day 1 post-dose sampling time (day after dose)", function(x) format(x)),
    wm = dcfg("trial_design.yaml", c("weight", "base", "mean"), "body weight mean (kg)", num_fmt(0)), wsd = dcfg("trial_design.yaml", c("weight", "base", "sd"), "body weight SD (kg)", num_fmt(0)),
    male = dcfg("trial_design.yaml", c("weight", "sex_ratio_male", "value"), "share male (%)", function(x) fnum(100 * as.numeric(x), 0)),
    split = f_split(), n_rand = f_n_rand(),
    cv = dcfg("prereg_20260926.yaml", c("section3", "base_cv_pct"), "protocol CV (%)", num_fmt(0)),
    g = dv(TPf, wp("M1"), "gmr", 2, "", "true GMR of the power row"),
    pw1 = dv(TPf, wp("M1"), "analytic_pct", 1, "%", "AUClast + Cmax power at the protocol n, M1"), pw0 = dv(TPf, wp("M0"), "analytic_pct", 1, "%", "AUClast + Cmax power at the protocol n, M0"),
    n90_0 = dint(NNf, wn("M0", 90), "n_evaluable_per_arm", "n per arm for the higher target, M0"), n85_0 = dint(NNf, wn("M0", 85), "n_evaluable_per_arm", "n per arm for the lower target, M0")))))
  deck_end()
}
