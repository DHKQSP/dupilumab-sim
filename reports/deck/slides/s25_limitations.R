# S25 한계: 모델 한계 4개(표: 정량한계 아래 곡선 모양·Km 고정, Cmax 과대예측, 200 mg 제형, ADA)와
# 검증 절차 한계 2개(카드: NCA 엔진의 실제 Phoenix 대조 미실시, 사람의 독립 QC 대기). [문헌+모의]
# Cmax 과대예측 "최대 27%"는 2016 모델의 600 mg(300 mg x 2) 자료에서 나온 값이다. 300 mg 자료만의 최대와 2020 모델 최대를 함께 적는다.
# Phoenix 대조 미실시는 보고서 한계 목록(6.5절, 부록 B)에 없는 새 항목이다(results/nca_engine에 Phoenix 출력 없음).

# 결과 열 이름에 담긴 기준값(예: coverage_lt80_pct의 80)을 읽는다
s25_col_thr <- function(rel, col, item) {
  premise(col %in% names(.read(rel)), sprintf("%s has column %s", rel, col))
  x <- as.numeric(sub("^[a-z_]*lt([0-9]+)_pct$", "\\1", col)); premise(is.finite(x), sprintf("threshold in column name %s", col))
  dderived(item, rel, sprintf("column name %s :: threshold", col), x, fnum(x, 0))
}
# 문자열 열에서 정규식으로 숫자 하나를 읽는다(행들이 같은 값이어야 한다). 원문에 한글이 있어도 추적 행에는 숫자와 영문 조건만 남긴다
s25_parse <- function(rel, where, col, rx, item) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where))
  x <- unique(as.numeric(vapply(r[[col]], function(s) { m <- regmatches(s, regexec(rx, s))[[1]]; if (length(m) < 2) NA_character_ else m[2] }, "")))
  premise(length(x) == 1 && is.finite(x), sprintf("one number matching %s in %s [%s] %s", rx, rel, where, col))
  dderived(item, rel, sprintf("%s :: %s, regex %s", where, col, rx), x, format(x))
}
# 곡선 모양 민감도 변형의 배율(config/scenarios.yaml). 슬라이드 21과 같은 표기: 값이 둘이면 "a·×b"(두 점만 모의), 셋 이상이면 "a~×b"(범위 끝)
s25_mult <- function(variants, par, item) {
  y <- .read("config/scenarios.yaml")$sensitivity_variants
  x <- sort(vapply(variants, function(v) as.numeric(y[[v]]$theta_multipliers[[par]]), 0))
  premise(all(is.finite(x)), paste("multipliers of", paste(variants, collapse = ", ")))
  p <- if (length(x) == 1L) format(x) else if (length(x) == 2L) sprintf("%s·×%s", format(x[1]), format(x[2])) else sprintf("%s~×%s", format(min(x)), format(max(x)))
  dderived(item, "config/scenarios.yaml", sprintf("sensitivity_variants [%s] :: theta_multipliers.%s, %s", paste(variants, collapse = ", "), par, if (length(x) == 2L) "both values" else "min and max"), x, p)
}

slide_S25 <- function() {
  CSH <- "curve_shape/curve_shape_B0.csv"; PC <- "rationale/pillar1_coverage_B0.csv"; INV <- "oc/inversion_all.csv"
  Q16 <- "step1/step1b_quant_gate.csv"; Q20 <- "results/step1_k2020/step1b_quant_gate.csv"; ADA <- "individual/individual_ada10_by_subgroup.csv"
  EV <- "nca_engine/engine_validation_summary.csv"; QC <- "regulatory/tables/verification_qc.csv"
  deck_slide("S25", tag = "litsim")
  deck_kicker(tx("S25.kicker")); deck_title(tx("S25.title"))

  # ---- 전제 검사 --------------------------------------------------------------------------------------------------------------
  P16 <- .read("config/params_typical.yaml"); P20 <- .read("config/params_k2020_model1.yaml"); PV <- .read("config/params_variability.yaml")
  premise(isTRUE(P16$theta$Km$fixed) && isTRUE(P20$theta$Km$fixed) && P16$theta$Km$value == P20$theta$Km$value, "Km fixed at the same value in both models")
  premise(PV$iiv_omega2$Km$omega2 == 0 && grepl("Km:\\s*\\{sd: 0, omega2: 0\\}", paste(readLines(proj_path("config", "params_k2020_model1.yaml"), encoding = "UTF-8"), collapse = "\n")), "no between-subject variability of Km in either model")
  cs <- rows(CSH); ok <- cs[stress_test == FALSE]
  premise(length(unique(round(ok$coverage_lt80_pct, 2))) == 1, "share below 80% window coverage is the same in all non-stress variants (table: all variants)")
  # 창 포착률 5백분위수: 노출 확인 범위(±15%) 밖인 비스트레스 변형(Vmax x0.8)을 빼도 같은 값(노트)
  gok <- ok[abs(AUClast_ratio_vs_obs544 - 1) <= 0.15]
  premise(identical(sort(setdiff(ok$variant, gok$variant)), "vmax080_both") &&
          max(c(ok$extrap_true_p95, rows(PC, "group=='all'")$extrap_true_p95)) == max(c(gok$extrap_true_p95, rows(PC, "group=='all'")$extrap_true_p95)),
          "smallest 5th-percentile coverage unchanged without the non-stress variant outside the exposure check (Vmax x0.8)")
  km_v <- ok[grepl("^km", variant), variant]; vm_v <- ok[grepl("^vmax", variant), variant]
  premise(length(km_v) >= 2 && length(vm_v) >= 2 && all(cs[stress_test == TRUE, grepl("^vmax", variant)]), "curve-shape variants: Km and Vmax (both arms), the stress test is a Vmax variant")
  g16 <- rows(Q16, "gate_role=='gate'"); g20 <- rows(Q20, "gate_role=='gate'")
  premise(g16[which.max(Cmax_ratio), dose_mg] == 600 && g20[which.max(Cmax_ratio), dose_mg] == 600, "the largest Cmax ratio is in a 600 mg data set in both models")
  premise(all(g16[dose_mg == 600]$Cmax_ratio > 1) && all(g20[dose_mg == 600]$Cmax_ratio > 1), "600 mg Cmax over-predicted in both models (text: over-prediction)")
  premise(all(g16[dose_mg == 600, presentation] == "2 x 300 mg"), "600 mg given as two 300 mg injections")
  dose <- as.numeric(.read("config/trial_design.yaml")$dose_mg); premise(all(g16[dose_mg != 600, dose_mg] == dose), "the other study-presentation data sets are at the study dose")
  # 시험 제형(300 mg) 재현: tmax(관측이 있는 자료) 일치, AUC0-last 모의/관측이 두 모델 모두 노출 gate(±tol%) 안
  tol <- as.numeric(.read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct); premise(is.finite(tol), "exposure gate tolerance in config")
  premise(all(abs(g16[dose_mg == dose, AUClast_ratio] - 1) <= tol / 100) && all(abs(g20[dose_mg == dose, AUClast_ratio] - 1) <= tol / 100), "study-dose AUC0-last within the exposure gate in both models (row 3: exposure reproduced)")
  t3 <- g16[dose_mg == dose & !is.na(tmax_obs_median)]; t3b <- g20[dose_mg == dose & !is.na(tmax_obs_median)]
  premise(nrow(t3) >= 1 && all(t3$tmax_obs_median == t3$tmax_sim_median) && all(t3b$tmax_obs_median == t3b$tmax_sim_median), "study-dose tmax reproduced in both models (row 3)")
  # Km 역산: 시험군 Km을 바꿔도 참 AUC0-inf 비가 동등 한계에 닿지 않는다(행 1: AUC로 검출되지 않음)
  lim <- as.numeric(.read("config/trial_design.yaml")$be$limits); kinv <- rows(INV, "mechanism=='Km'")
  kx <- c(kinv[is.finite(end_auc_ratio), end_auc_ratio], kinv[reachable == TRUE, auc_ratio])
  premise(length(kx) > 0 && min(kx) > lim[1] && max(kx) < lim[2] && !any(kinv[reachable == TRUE, target] <= lim[1] | kinv[reachable == TRUE, target] >= lim[2]), "Km differences never reach the equivalence limits (row 1: not detected by AUC)")
  premise(isTRUE(cs[variant == "vmax050_both", stress_test]) && abs(cs[variant == "vmax050_both", AUClast_ratio_vs_obs544] - 1) > tol / 100, "the Vmax x0.5 stress test is outside the exposure gate (row 1: inconsistent with observed exposure)")
  qc <- row1(QC, "Activity=='Independent human QC of this report'"); premise(startsWith(qc$Result, "PENDING"), "independent human QC pending")
  ev <- rows(EV); premise(all(ev$pass) && setequal(unique(ev$comparison), c("this engine vs NonCompart", "PKNCA vs NonCompart", "this engine vs PKNCA")), "the NCA engine was compared with NonCompart and PKNCA only")
  premise(!any(grepl("phoenix|winnonlin", list.files(proj_path("results", "nca_engine")), ignore.case = TRUE)), "no Phoenix output in results/nca_engine")

  # ---- 수치 -------------------------------------------------------------------------------------------------------------------
  thr <- s25_col_thr(CSH, "coverage_lt80_pct", "window coverage threshold (%) of the column coverage_lt80_pct")
  cov_min95 <- 100 - max(c(ok$extrap_true_p95, rows(PC, "group=='all'")$extrap_true_p95))
  dsrc("both-model coverage used in the smallest 5th-percentile window coverage", PC, "(table)")
  f <- list(
    km = dcfg("params_typical.yaml", c("theta", "Km", "value"), "Km (mg/L), fixed in both models", function(x) format(x)),
    k = s25_mult(km_v, "Km", "Km multipliers, both arms, curve-shape variants"),
    v = s25_mult(vm_v, "Vmax", "Vmax multipliers, both arms, non-stress curve-shape variants"),
    thr = thr,
    lt80 = dext(CSH, "stress_test==FALSE", "coverage_lt80_pct", max, 2, "%", "share below 80% window coverage, non-stress curve-shape variants (all equal)"),
    cov95 = dderived("smallest 95th-percentile coverage (100 - max 95th percentile of true extrapolation) over both models and non-stress curve-shape variants", CSH,
                     "stress_test==FALSE :: 100 - max(extrap_true_p95), with rationale/pillar1_coverage_B0.csv group=='all'", cov_min95, fnum(floor(cov_min95 * 10) / 10, 1)),
    im = local({ r <- rows(INV, "mechanism=='Km'"); x <- range(r$end_multiplier)
      dderived("Km multipliers searched in the test arm (range ends)", INV, "mechanism=='Km' :: range(end_multiplier)", x, sprintf("%s~×%s", format(x[1]), format(x[2]))) }),
    sv = s25_mult(cs[stress_test == TRUE, variant], "Vmax", "Vmax multiplier of the stress test"),
    gate = dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "exposure gate: AUClast mean within +/- tolerance (%)", num_fmt(0)),
    kr = local({ r <- rows(INV, "mechanism=='Km' & is.finite(end_auc_ratio)"); r2 <- rows(INV, "mechanism=='Km' & reachable==TRUE"); x <- range(c(r$end_auc_ratio, r2$auc_ratio))
      dderived("true AUC0-inf ratio range over Km x0.01 to x100 (range ends and reached rows)", INV, "mechanism=='Km' :: end_auc_ratio, auc_ratio", x, sprintf("%s~%s", fnum(x[1], 3), fnum(x[2], 3))) }),
    gc16 = drange(Q16, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2016"),
    gc20 = drange(Q20, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2020"),
    c27 = dderived("largest Cmax over-prediction (%), 2016 model, study presentation", Q16, "gate_role=='gate' :: (max(Cmax_ratio) - 1) x 100", 100 * (max(g16$Cmax_ratio) - 1), fnum(100 * (max(g16$Cmax_ratio) - 1), 0)),
    c300 = local({ r <- rows(Q16, "gate_role=='gate' & dose_mg==300"); x <- 100 * (max(r$Cmax_ratio) - 1)
      dderived("largest Cmax over-prediction (%), 2016 model, 300 mg data sets", Q16, "gate_role=='gate' & dose_mg==300 :: (max(Cmax_ratio) - 1) x 100", x, fnum(x, 0)) }),
    c20 = dderived("largest Cmax over-prediction (%), 2020 model, study presentation", Q20, "gate_role=='gate' :: (max(Cmax_ratio) - 1) x 100", 100 * (max(g20$Cmax_ratio) - 1), fnum(100 * (max(g20$Cmax_ratio) - 1), 0)),
    c20b = local({ r <- rows(Q20, "gate_role=='gate' & dose_mg==300"); x <- 100 * (max(r$Cmax_ratio) - 1)
      dderived("largest Cmax over-prediction (%), 2020 model, 300 mg data sets", Q20, "gate_role=='gate' & dose_mg==300 :: (max(Cmax_ratio) - 1) x 100", x, fnum(x, 0)) }),
    d600 = dint(Q16, sprintf("id=='%s'", g16[which.max(Cmax_ratio), id]), "dose_mg", "dose of the data set with the largest Cmax ratio (mg)"),
    d300 = f_dose(),
    d200 = drange(Q16, "gate_role=='external'", "dose_mg", 0, "", "dose of the 175 mg/mL external data sets (mg)"),
    c175 = s25_parse(Q16, "gate_role=='external'", "presentation", "([0-9.]+) mg/mL", "concentration of the 200 mg presentation (mg/mL)"),
    c150 = s25_parse(Q16, "gate_role=='gate' & dose_mg==300", "presentation", "([0-9.]+) mg/mL", "concentration of the study presentation (mg/mL)"),
    v2 = s25_parse(Q16, "gate_role=='gate' & dose_mg==300", "presentation", "^([0-9.]+) mL of", "volume of the study presentation (mL)"),
    to = drange(Q16, "gate_role=='external' & !is.na(tmax_obs_median)", "tmax_obs_median", 1, "", "observed median tmax, 200 mg (175 mg/mL) data sets (day)"),
    ts = drange(Q16, "gate_role=='external'", "tmax_sim_median", 0, "", "simulated median tmax, 200 mg data sets, 2016 (day)"),
    t300 = drange(Q16, "gate_role=='gate' & !is.na(tmax_obs_median)", "tmax_obs_median", 1, "", "observed median tmax, study presentation (day)"),
    a16 = drange(Q16, "gate_role=='external'", "AUClast_ratio", 2, "", "AUC0-last sim/obs, 200 mg (175 mg/mL) data sets, 2016"),
    a20 = drange(Q20, "gate_role=='external'", "AUClast_ratio", 2, "", "AUC0-last sim/obs, 200 mg (175 mg/mL) data sets, 2020"),
    frac = dcfg("trial_design.yaml", c("ada_sensitivity", "fraction"), "ADA-like subgroup share (%)", function(x) fnum(100 * x, 0)),
    onset = dcfg("trial_design.yaml", c("ada_sensitivity", "onset_day"), "ADA-like onset (day after dose)", num_fmt(0)),
    kem = dcfg("trial_design.yaml", c("ada_sensitivity", "ke_multiplier"), "ADA-like ke multiplier", function(x) format(x)),
    e1 = dv(ADA, "ada==1 & schedule=='B0'", "extrap_true_median", 2, "", "ADA-like subgroup, median true extrapolation (%)"),
    e0 = dv(ADA, "ada==0 & schedule=='B0'", "extrap_true_median", 2, "", "other subjects, median true extrapolation (%)"),
    alt = dv(ADA, "ada==1 & schedule=='B0'", "coverage_lt80_pct", 2, "", "ADA-like subgroup, window coverage below 80% (%)"),
    lz = local({ a <- sum(ev$lz_points_identical); b <- sum(ev$lz_points_identical + ev$lz_points_mismatch)
      dderived("lambda-z windows identical over all comparisons", EV, "sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", c(a, b), sprintf("%s / %s", fint(a), fint(b))) }),
    mx = local({ x <- max(ev$max_rel_diff); e <- floor(log10(x)); m <- signif(x, 2) / 10^e   # 표시: 4.9×10⁻¹³ (지수는 위첨자)
      sup <- chartr("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079", as.character(e))
      dderived("largest relative parameter difference", EV, "max(max_rel_diff), printed as mantissa x 10^exponent", x, sprintf("%s\u00d710%s", format(m), sup)) }),
    qc = local({ .record("independent human QC status", QC, "Activity=='Independent human QC of this report' :: Result", qc$Result, qc$Result); qc$Result }))

  # ---- 왼쪽: 모델 한계 표 -------------------------------------------------------------------------------------------------------
  H <- DK$txt$S25$table
  df <- data.frame(a = tx("S25.table.col1", f), b = tx("S25.table.col2", f), c = tx("S25.table.col3", f), check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- tx("S25.table.head")
  if (nzchar(Sys.getenv("S25_DEBUG"))) { data.table::fwrite(df, Sys.getenv("S25_DEBUG")); writeLines(c(tx("S25.card_nca", f), "", tx("S25.card_qc", f)), paste0(Sys.getenv("S25_DEBUG"), ".cards")) }
  tw <- 8.4
  deck_text(tx("S25.left_label"), c(GEO$ML, GEO$BODY_TOP - 0.04, tw, 0.45), size = 16, bold = TRUE, color = PAL$ink2, label = "label_model")
  tsz <- 13.5; tws <- c(1.95, 3.35, 3.1); ty <- GEO$BODY_TOP + 0.43
  deck_table(df, box = c(GEO$ML, ty, tw, GEO$BODY_BOTTOM - ty), widths = tws, size = tsz, align_num = FALSE)
  # 표의 렌더링 높이 추정(deck_table과 같은 줄 수 x 표 줄 높이 1.19 x 글자 크기 + 행마다 셀 여백 8 pt): 오른쪽 카드의 아래 끝을 표의 아래 선에 맞춘다
  twi <- tws / sum(tws) * tw
  nl <- function(v, w, b) vapply(as.character(v), function(s) est_lines(s, w - 0.14, tsz, b), 1L)
  tl <- max(mapply(function(v, w) nl(v, w, TRUE), names(df), twi)) + sum(apply(sapply(seq_along(twi), function(j) nl(df[[j]], twi[j], j == 1)), 1, max))
  tab_bot <- min(GEO$BODY_BOTTOM, ty + tl * tsz * 1.19 / 72 + (nrow(df) + 1) * 8 / 72)

  # ---- 오른쪽: 검증 절차 한계 카드 두 개 ------------------------------------------------------------------------------------------
  xr <- GEO$ML + tw + 0.25; wr <- GEO$W - GEO$MR - xr
  deck_text(tx("S25.right_label"), c(xr, GEO$BODY_TOP - 0.04, wr, 0.45), size = 16, bold = TRUE, color = PAL$ink2, label = "label_verif")
  # 두 카드 높이는 추정 글 높이에 비례해 나눈다(위 끝 = 표 위 끝, 아래 끝 = 표 아래 선)
  y0 <- ty; cg <- 0.18; cn <- tx("S25.card_nca", f); cq <- tx("S25.card_qc", f)
  he <- c(est_height(cn, wr, 16, 5, card = TRUE), est_height(cq, wr, 16, 5, card = TRUE)); ch <- he / sum(he) * (tab_bot - y0 - cg)
  deck_text(cn, c(xr, y0, wr, ch[1]), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "card_nca", gap_pt = 5)
  deck_text(cq, c(xr, y0 + ch[1] + cg, wr, ch[2]), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "card_qc", gap_pt = 5)

  deck_notes(tx("S25.notes", c(f, list(
    pl = local({ col <- "extrap_true_p95"; premise(col %in% names(.read(CSH)), "percentile column"); x <- as.numeric(sub("^.*_p([0-9]+)$", "\\1", col))
      dderived("percentile level of the column extrap_true_p95", CSH, sprintf("column name %s :: percentile", col), x, fnum(x, 0)) }),
    kmed = drange(CSH, "stress_test==FALSE & grepl('^km', variant)", "extrap_true_median", 2, "%", "Km variants, median true extrapolation"),
    kp95 = drange(CSH, "stress_test==FALSE & grepl('^km', variant)", "extrap_true_p95", 2, "%", "Km variants, 95th percentile of true extrapolation"),
    vp95 = drange(CSH, "stress_test==FALSE & grepl('^vmax', variant)", "extrap_true_p95", 2, "%", "Vmax variants, 95th percentile of true extrapolation"),
    sp95 = dv(CSH, "variant=='vmax050_both'", "extrap_true_p95", 2, "%", "stress test, 95th percentile of true extrapolation"),
    slt = dv(CSH, "variant=='vmax050_both'", "coverage_lt80_pct", 3, "%", "stress test, window coverage below 80%"),
    tr = local({ r <- rows(INV, "mechanism=='Km' & reachable==TRUE"); premise(nrow(r) == 2 && all(abs(r$target - 1.05) < 1e-9) && setequal(r$model, c("k2016", "k2020")), "the only reachable Km target is the same in both models")
      dv(INV, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "target", 2, "", "only reachable Km target ratio") }),
    r16 = dv(INV, "mechanism=='Km' & reachable==TRUE & model=='k2016'", "multiplier", 0, "", "Km multiplier reaching the target, 2016"),
    r20 = dv(INV, "mechanism=='Km' & reachable==TRUE & model=='k2020'", "multiplier", 0, "", "Km multiplier reaching the target, 2020"),
    ce16 = drange(Q16, "gate_role=='external'", "Cmax_ratio", 2, "", "Cmax sim/obs, 200 mg data sets, 2016"),
    ce20 = drange(Q20, "gate_role=='external'", "Cmax_ratio", 2, "", "Cmax sim/obs, 200 mg data sets, 2020"),
    tl1 = dv(ADA, "ada==1 & schedule=='B0'", "tlast_median", 1, "", "ADA-like subgroup, median tlast (day)"),
    tl0 = dv(ADA, "ada==0 & schedule=='B0'", "tlast_median", 1, "", "other subjects, median tlast (day)"),
    ng = dcount(Q16, "gate_role=='gate'", "study-presentation data sets"),
    g300a = drange(Q16, sprintf("gate_role=='gate' & dose_mg==%s", dose), "AUClast_ratio", 2, "", "AUC0-last sim/obs, study-dose data sets, 2016"),
    g300b = drange(Q20, sprintf("gate_role=='gate' & dose_mg==%s", dose), "AUClast_ratio", 2, "", "AUC0-last sim/obs, study-dose data sets, 2020"),
    v80 = s25_mult("vmax080_both", "Vmax", "Vmax multiplier of the non-stress variant outside the exposure check"),
    n_arm = f_n_arm()))))
  deck_end()
}
