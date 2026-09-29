# S23 통계분석계획(SAP) 제안 요지: 1차(AUC0-last + Cmax, 주분석 M1(의뢰자 결정 전), M0 민감도), 이차(AUC0-inf, 기준 세트 (i), 두 분석군, adjusted R² 0.90 인원 병기),
# 민감도(규칙 C 세트 (i), 항약물항체 상태별), FDA 요구 시 대응안(규칙 B 공동 1차). 카드 네 개(도식). [문헌+모의]
# 출처: regulatory/sap_text_proposals_en.md 1~7절의 문안(수치는 같은 결과 파일에서 SAP 문서와 같은 행 조건으로 읽는다).
# 주의: SAP 문안은 주분석 모형을 "[M0 or M1, to be selected by the sponsor]"로 두고 7절 상태가 "pending (sponsor)"다. 덱은 M1을 기준으로 하되 확인 대기로 적는다.
# 규칙 C 수치는 세트 (i)(criteria_g2_type1.csv G2_C_i = type1_models.csv G2_Ci). 보고서 요약·6.3절의 세트 (ii) 값(G2_Cii)과 섞지 않는다.
S23_T1 <- "oc_models/type1_models.csv"; S23_CG <- "criteria/criteria_g2_type1.csv"
s23_w <- function(am, cf, extra = "") sprintf("analysis_model=='%s' & config=='%s'%s", am, cf, extra)
# 한 구성의 최댓값 칸(행 조건에 모델·시나리오를 넣어 SAP 문서의 max_t1과 같은 locator)
s23_max <- function(rel, am, cf, ci = FALSE, item = sprintf("%s %s maximum", am, cf)) {
  r <- rows(rel, s23_w(am, cf)); premise(nrow(r) == 16, sprintf("%s %s: 16 boundary cells", am, cf)); rr <- r[which.max(pass_pct)]
  w <- s23_w(am, cf, sprintf(" & pk_model=='%s' & scenario=='%s'", rr$pk_model, rr$scenario))
  if (ci) dci(rel, w, "pass_pct", "lo", "hi", 2, "%", item) else dv(rel, w, "pass_pct", 2, "%", item)
}

slide_S23 <- function() {
  T1 <- S23_T1; CG <- S23_CG; TPF <- "trialpop/tp_failure_by_set.csv"
  deck_slide("S23", tag = "litsim")

  # ---- 전제 ----
  for (am in c("M0", "M1")) {
    a <- rows(CG, s23_w(am, "G2_C_i")); b <- rows(T1, s23_w(am, "G2_Ci"))
    premise(nrow(a) == 16 && nrow(b) == 16 && isTRUE(all.equal(max(a$pass_pct), max(b$pass_pct))), sprintf("rule C set (i) maximum identical in the criteria file and type1_models (%s)", am))
  }
  mx1 <- rows(CG, s23_w("M1", "G2_C_i")); mx1 <- mx1[which.max(pass_pct)]; p2 <- rows(T1, s23_w("M1", "P2")); p2 <- p2[which.max(pass_pct)]
  premise(mx1$pk_model == p2$pk_model && mx1$scenario == p2$scenario, "rule C (i) and P2 are highest in the same boundary cell under M1 (notes)")
  premise(nrow(rows(T1, s23_w("M1", "F3B", " & class=='conservative'"))) == 16, "fallback F3B conservative in all 16 cells under M1 (text)")
  # 대응안의 규칙 A·C 민감도(F3A, F3C)는 세트 (ii) 탈락 정의(AUCinf_A, AUCinf_C)로 계산되었고 세트 (i) 세 평가변수 구성은 없다(노트)
  OC <- .read("config/oc_design.yaml")$configurations; CB <- "criteria/criteria_bias.csv"
  premise(identical(unlist(OC$F3A$endpoints), c("AUClast", "Cmax", "AUCinf_A")) && identical(unlist(OC$F3C$endpoints), c("AUClast", "Cmax", "AUCinf_C")) &&
          all(rows(CB, "endpoint=='AUCinf_A'")$config == "A_ii") && all(rows(CB, "endpoint=='AUCinf_C'")$config == "C_ii") &&
          !any(grepl("^F3.*i$", unique(rows(T1)$config))), "fallback rule A and C sensitivity use set (ii); no set (i) three-endpoint configuration (notes)")
  premise(nrow(rows(CG, s23_w("M1", "G2_C_i", " & class=='nominal'"))) >= 1 && rows(CG, s23_w("M1", "G2_C_i"))[which.max(pass_pct)]$class == "nominal", "rule C set (i) M1 largest cell is nominal (text)")

  f <- list(r1 = drange(T1, s23_w("M1", "P2"), "pass_pct", 2, "%", "M1 P2 boundary type I error range"),
            r0 = drange(T1, s23_w("M0", "P2"), "pass_pct", 2, "%", "M0 P2 boundary type I error range"),
            c1 = dcount(T1, s23_w("M1", "P2", " & class=='conservative'"), "M1 P2 cells classified conservative"),
            e1 = dcount(T1, s23_w("M1", "P2", " & class=='exceeding'"), "M1 P2 cells classified exceeding"),
            r2i = f_set("i", "r2"), exi = f_set("i", "extrap"), r2iii = f_set("iii", "r2"), exiii = f_set("iii", "extrap"),
            fi = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
            fiii = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models"),
            gc0 = s23_max(CG, "M0", "G2_C_i", item = "M0 rule C set (i) largest boundary value"),
            gc1 = s23_max(CG, "M1", "G2_C_i", item = "M1 rule C set (i) largest boundary value"),
            fb1 = drange(T1, s23_w("M1", "F3B"), "pass_pct", 2, "%", "M1 F3B boundary type I error range"),
            gbn = dcount(T1, s23_w("M1", "G2_B", " & pass_pct > 5"), "M1 G2_B cells with point estimate above 5%"),
            gbl = dcount(T1, s23_w("M1", "G2_B", " & lo > 5"), "M1 G2_B cells with Wilson lower bound above 5%"),
            ncell = dcount(T1, s23_w("M1", "G2_B"), "boundary cells per configuration (M1)"),
            nom = f_nominal())
  deck_kicker(tx("S23.kicker")); deck_title(tx("S23.title", f))

  # ---- 카드 네 개(2 x 2): 1차, 이차 / 민감도, 대응안 ----
  gap <- 0.22; cw <- (GEO$CW - gap) / 2; y1 <- GEO$BODY_TOP; hh <- 0.46
  chs <- c(2.6, GEO$BODY_BOTTOM - y1 - gap - 2.6)   # 위 줄(1차, 이차)이 문단이 많아 더 높다
  cards <- list(list(k = "primary", fill = PAL$tint_orange, col = PAL$orange), list(k = "secondary", fill = PAL$tint_blue, col = PAL$blue),
                list(k = "sens", fill = PAL$tint_grey, col = PAL$ink), list(k = "fallback", fill = PAL$tint_grey, col = PAL$ink2))
  for (i in seq_along(cards)) {
    cd <- cards[[i]]; rw <- (i - 1) %/% 2 + 1; x <- GEO$ML + ((i - 1) %% 2) * (cw + gap); y <- y1 + (rw - 1) * (chs[1] + gap); ch <- chs[rw]
    deck_text(" ", c(x, y, cw, ch), size = 16, bg = cd$fill, geom = "roundRect", label = sprintf("card_%s", cd$k))
    deck_text(tx(sprintf("S23.cards.%s.head", cd$k), f), c(x + 0.12, y + 0.08, cw - 0.24, hh), size = 18, bold = TRUE, color = cd$col, label = sprintf("head_%s", cd$k))
    deck_text(tx(sprintf("S23.cards.%s.body", cd$k), f), c(x + 0.12, y + 0.08 + hh, cw - 0.24, ch - hh - 0.14), size = 16, label = sprintf("body_%s", cd$k), gap_pt = 5)
  }

  # ---- 노트 ----
  SS <- "oc_models/sd_se_models.csv"; PW <- "oc_models/power_models.csv"; RS <- "reliability/reliability_two_flag_sets_summary.csv"; EX <- "sap/sap_exclusion_support.csv"
  AD <- "individual/individual_ada10_by_subgroup.csv"
  ex_ <- function(m, s_, col, d, it) dv(EX, sprintf("model=='%s' & set=='%s'", m, s_), col, d, "", it)
  ad_ <- function(g, col, d, it) dv(AD, sprintf("ada==%d & schedule=='B0'", g), col, d, "", it)
  pw_ <- function(am, cf = "P2") dci(PW, sprintf("pk_model=='k2016' & scenario=='S00' & analysis_model=='%s' & config=='%s'", am, cf), "pass_pct", "lo", "hi", 1, "%", sprintf("power k2016 S00 %s %s", am, cf))
  deck_notes(tx("S23.notes", c(f, list(
    split = f_split(), wt = f_wt_range(),
    sd0 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M0' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE M0 AUCinf_true boundary cells"),
    sd1 = drange(SS, "endpoint=='AUCinf_true' & analysis_model=='M1' & !scenario %in% c('S00','F097')", "sd_se_ratio", 3, "", "SD/SE M1 AUCinf_true boundary cells"),
    pw1 = pw_("M1"), pw0 = pw_("M0"), pwb1 = pw_("M1", "F3B"),
    p2m1 = s23_max(T1, "M1", "P2", TRUE, "M1 P2 maximum"), p2m0 = s23_max(T1, "M0", "P2", TRUE, "M0 P2 maximum"),
    lz = drange(RS, "variant %in% c('base','struct2020')", "lambda_ok_pct", 1, "%", "lambda-z estimable, two models"),
    ri = drange(RS, "variant %in% c('base','struct2020')", "reliable_i_pct", 1, "%", "reliability criteria (i), two models"),
    wd16 = ex_("k2016", "criteria_i", "wt_diff", 1, "criteria (i) exclusion k2016 wt_diff"), wd20 = ex_("k2020", "criteria_i", "wt_diff", 1, "criteria (i) exclusion k2020 wt_diff"),
    gr16 = ex_("k2016", "criteria_i", "aucinf_true_gm_ratio_excluded_to_retained", 2, "criteria (i) exclusion k2016 true AUC0-inf GM ratio"),
    gr20 = ex_("k2020", "criteria_i", "aucinf_true_gm_ratio_excluded_to_retained", 2, "criteria (i) exclusion k2020 true AUC0-inf GM ratio"),
    exb16 = paste0(ex_("k2016", "rule_B", "excluded_pct", 1, "rule B exclusion k2016"), "%"), exb20 = paste0(ex_("k2020", "rule_B", "excluded_pct", 1, "rule B exclusion k2020"), "%"),
    gc0ci = s23_max(CG, "M0", "G2_C_i", TRUE, "M0 rule C set (i) largest with CI"), gc1ci = s23_max(CG, "M1", "G2_C_i", TRUE, "M1 rule C set (i) largest with CI"),
    gco1 = dcount(CG, s23_w("M1", "G2_C_i", " & pass_pct > 5"), "M1 G2_C_i cells above 5% (point)"),
    gcl1 = dcount(CG, s23_w("M1", "G2_C_i", " & lo > 5"), "M1 G2_C_i cells with Wilson lower bound above 5%"),
    gcii0 = s23_max(CG, "M0", "G2_C_ii", item = "M0 rule C set (ii) largest (not used)"), gcii1 = s23_max(CG, "M1", "G2_C_ii", item = "M1 rule C set (ii) largest (not used)"),
    ada = dcfg("trial_design.yaml", c("other_schedules_out_of_scope", "ADA_days"), "ADA sampling days after the pre-dose sample (non-zero entries)", function(x) paste(fnum(x[x > 0], 0), collapse = ", ")),
    adaf = dcfg("trial_design.yaml", c("ada_sensitivity", "fraction"), "ADA-like subgroup share (%)", function(x) paste0(fnum(100 * x, 0), "%")),
    onset = dcfg("trial_design.yaml", c("ada_sensitivity", "onset_day"), "ADA-like subgroup onset (day)", num_fmt(0)),
    kem = dcfg("trial_design.yaml", c("ada_sensitivity", "ke_multiplier"), "ADA-like subgroup linear elimination multiplier", num_fmt(0)),
    tl1 = ad_(1L, "tlast_median", 1, "ADA-like subgroup median tlast (days)"), tl0 = ad_(0L, "tlast_median", 1, "other subjects median tlast (days)"),
    ex1 = paste0(ad_(1L, "extrap_true_median", 2, "ADA-like subgroup median true extrapolation (%)"), "%"),
    ex0 = paste0(ad_(0L, "extrap_true_median", 2, "other subjects median true extrapolation (%)"), "%"),
    fb0 = drange(T1, s23_w("M0", "F3B"), "pass_pct", 2, "%", "M0 F3B boundary type I error range"),
    fbm1 = s23_max(T1, "M1", "F3B", TRUE, "M1 F3B maximum"),
    gbm1 = s23_max(T1, "M1", "G2_B", TRUE, "M1 G2_B maximum"),
    gbn0 = dcount(T1, s23_w("M0", "G2_B", " & pass_pct > 5"), "M0 G2_B cells with point estimate above 5%"),
    gbl0 = dcount(T1, s23_w("M0", "G2_B", " & lo > 5"), "M0 G2_B cells with Wilson lower bound above 5%"),
    fa1 = drange(T1, s23_w("M1", "F3A"), "pass_pct", 2, "%", "M1 F3A (set ii) boundary type I error range"),
    fc1 = drange(T1, s23_w("M1", "F3C"), "pass_pct", 2, "%", "M1 F3C (set ii) boundary type I error range")))))
  deck_end()
}
