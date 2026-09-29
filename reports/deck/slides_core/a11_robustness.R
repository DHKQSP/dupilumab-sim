# A11 별첨: 견고성 점검(근거 아님). 체중 범위·아토피 피부염 성인 체중 분포(결과보고 덱 A5)와 가정 변경 요약(결과보고 덱 S21).
# 핵심 덱에서 results/atopic/, results/weight_generalization/ 파일과 "아토피", "환자 모집단"이라는 말을 쓰는 유일한 슬라이드다
# (check_deck 4: text/ko_core/00_common.yaml atopic_allowed_slides A11). 근거는 시험 모집단(건강인, 체중 층화, B0)뿐이다.
# 왼쪽 그림(두 패널 나란히, 체중 축 공유, 회색 띠 = 연구 범위, 점선 = 개발 자료 범위 밖일 수 있는 구간; 파랑 = AUClast/AUCinf 5번째 백분위수(참값), 주황 = AUCinf 신뢰 기준 미달; 모델은 표식·선 모양):
#   (a) 체중 구간 40~150 kg(구간 안 균일 분포): weight_bands_B0_abcd.csv. 이 파일의 신뢰 비율은 세트 (i)뿐(reliable_rsq_extrap_pct; scripts/24).
#   (b) 아토피 피부염 성인 가정 체중 분포(주 분포): 세트 (iii)(= 신뢰할 수 있는 AUCinf) 미확보 비율(atopic_criteria_by_band.csv),
#       AUClast/AUCinf 5번째 백분위수(atopic_individual_pillar1.csv). 점의 x = 구간 체중 중앙값(a) 또는 평균(b).
# 오른쪽 표: 결과보고 덱 S21의 가정 다섯 개를 핵심 덱 용어로 요약(1종 오류는 분석 모형 M1, 신뢰 기준은 있는 곳에서 세트 (iii);
#   정량한계 시험은 세트 (ii)만 계산됨). 위 띠: 근거 아님. 아래 캡션: 세트 정의, 분석 모형, 제목 값의 뜻.
# 제목 수치: 세트 (iii) 미달 환자의 AUClast/AUCinf 최솟값(세 모델 변형, 내림; atopic_coverage_failing.csv window_fail_min).
# weight_bands_B0_abcd.csv의 dev_range_note 열은 한글이다. 행 조건에는 영문 코드(model, band)만 쓴다.

# 곡선 모양 민감도 변형의 배율(config/scenarios.yaml): 두 값이면 "a·×b"
a11_mult <- function(variants, par, item) {
  y <- .read("config/scenarios.yaml")$sensitivity_variants
  x <- sort(vapply(variants, function(v) as.numeric(y[[v]]$theta_multipliers[[par]]), 0)); premise(length(x) == 2 && all(is.finite(x)), paste("two multipliers of", par))
  dderived(item, "config/scenarios.yaml", sprintf("sensitivity_variants [%s] :: theta_multipliers.%s, both values", paste(variants, collapse = ", "), par), x,
           sprintf("%s·×%s", format(x[1]), format(x[2])))
}
a11_lo <- function(b) as.numeric(sub("-.*", "", b))

a11_fig <- function(WB, ABN, AIP, L, ML, k130) {
  wb <- rows(WB, "model %in% c('a','b')")
  a <- wb[, .(panel = "wt", pk = c(a = "k2016", b = "k2020")[model], x = WT_median, cov = 100 - extrap_true_p95, fail = 100 - reliable_rsq_extrap_pct)]
  V2 <- c(base = "k2016", struct2020 = "k2020")
  fb <- rows(ABN, "distribution=='primary' & set=='iii' & variant %in% c('base','struct2020') & band %in% c('below 60','60-90','above 90-100','above 100')")
  cb <- rows(AIP, "distribution=='primary' & variant %in% c('base','struct2020') & band %in% c('below 60','60-90','above 90-100','above 100')")
  premise(nrow(fb) == 8 && nrow(cb) == 8, "atopic: four weight bands x two models in both files")
  b <- merge(fb[, .(variant, band, fail = fail_pct)], cb[, .(variant, band, x = wt_mean, cov = 100 - extrap_true_p95)], by = c("variant", "band"))
  b <- b[, .(panel = "at", pk = V2[variant], x, cov, fail)]
  d <- melt(rbind(a, b), id.vars = c("panel", "pk", "x"), variable.name = "metric", value.name = "y")
  PN <- c(wt = L$panel_wt, at = L$panel_at)
  d[, panel := factor(PN[panel], levels = PN)]
  d[, grp := paste(pk, metric)]
  # 패널마다 x 범위: (a)는 150 kg 구간까지, (b)는 자료(가장 무거운 구간 평균)까지만(끝값 글자 자리 포함)
  xr_ <- d[, .(xmax = max(x)), by = panel]; premise(xr_[panel == PN[["at"]], xmax] < 130 && xr_[panel == PN[["wt"]], xmax] > 130, "panel (b) data stop below 130 kg, panel (a) data reach the heaviest band")
  blank <- data.table(panel = factor(PN[c("wt", "wt", "at", "at")], levels = PN), x = c(40, 168, 40, xr_[panel == PN[["at"]], xmax] + 34), y = 50)   # 끝값 글자가 패널 끝에 닿지 않게
  # 끝값 표시: 가장 무거운 구간의 미달 비율(두 모델)
  endv <- d[metric == "fail", .SD[which.max(x)], by = .(panel, pk)]
  # 지표 이름(직접 표시, 두 줄): 선과 겹치지 않는 높이. 파랑 이름은 파랑 선 아래; 주황 이름은 (a) 주황 선 위, (b) 주황 선 아래
  lab <- data.table(panel = factor(PN, levels = PN), cov = L$cov, fail = c(L$fail_i, L$fail_iii), ycov = c(87, 87), yfail = c(58, 22))
  premise(max(d[metric == "fail" & panel == PN[["wt"]], y]) < 35 && min(d[metric == "fail" & panel == PN[["at"]], y]) > 25 && min(d[metric == "cov", y]) > 90,
          "direct labels sit clear of the lines")
  shade <- data.table(panel = factor(PN, levels = PN))
  p <- ggplot(d, aes(x, y)) +
    geom_blank(data = blank, aes(x, y), inherit.aes = FALSE) +
    geom_rect(data = shade, aes(xmin = 60, xmax = 90, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = PAL$tint_grey) +
    geom_vline(data = shade[1], aes(xintercept = k130), linetype = "dotted", colour = PAL$ink2, linewidth = 0.7) +
    geom_line(aes(group = grp, colour = metric, linetype = pk), linewidth = 0.9) +
    geom_point(aes(colour = metric, shape = pk), size = 2.8) +
    geom_text(data = endv, aes(x, y, label = paste0(fnum(y, 1), "%"), colour = metric), hjust = -0.2, size = PT(14), family = FONT, show.legend = FALSE) +
    geom_text(data = lab, aes(x = 42, y = ycov, label = cov), inherit.aes = FALSE, hjust = 0, vjust = 1, size = PT(14), family = FONT, colour = PAL$blue, fontface = "bold", lineheight = 0.95) +
    geom_text(data = lab, aes(x = 42, y = yfail, label = fail), inherit.aes = FALSE, hjust = 0, vjust = 1, size = PT(14), family = FONT, colour = PAL$orange, fontface = "bold", lineheight = 0.95) +
    facet_wrap(~ panel, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = c(cov = PAL$blue, fail = PAL$orange), guide = "none") +
    scale_shape_manual(values = CORE_MODEL_SHAPE, breaks = names(CORE_MODEL_SHAPE), labels = unlist(ML[names(CORE_MODEL_SHAPE)]), name = NULL) +
    scale_linetype_manual(values = CORE_MODEL_LT, breaks = names(CORE_MODEL_LT), labels = unlist(ML[names(CORE_MODEL_LT)]), name = NULL) +   # breaks: 범례 순서 주 모델(2020) 먼저, 이름과 짝
    scale_x_continuous(breaks = c(60, 90, 120, 150), expand = expansion(add = c(1, 0))) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%"), expand = expansion(mult = c(0, 0.03))) +
    coord_cartesian(clip = "off") + labs(x = L$xlab, y = NULL) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.spacing.x = grid::unit(14, "pt"),
          strip.text = element_text(size = 16, face = "bold", colour = PAL$ink, hjust = 0, lineheight = 0.95, margin = margin(2, 0, 4, 0)), plot.margin = margin(4, 8, 4, 2))
  p
}

slide_A11 <- function() {
  WB <- "weight_generalization/weight_bands_B0_abcd.csv"; ABN <- "atopic/atopic_criteria_by_band.csv"; ACF <- "atopic/atopic_coverage_failing.csv"
  AIP <- "atopic/atopic_individual_pillar1.csv"; AW <- "atopic/atopic_weight_table.csv"; AB <- "atopic/atopic_trials_bias.csv"; ACC <- "atopic/atopic_trials_concordance.csv"
  AAF <- "atopic/atopic_trials_arm_failure.csv"; TPF <- "trialpop/tp_failure_by_set.csv"
  T1 <- "oc_models/type1_models.csv"; CG <- "criteria/criteria_g2_type1.csv"; CV <- "core_deck/coverage_by_case.csv"; LT <- "lloq/lloq_trial_type1.csv"
  TRS <- "trialpop/tp_residual_sensitivity.csv"
  deck_slide("A11", tag = "sim")
  L <- DK$txt$A11; ML <- DK$txt$common$models_short
  nom <- f_nominal(); nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  premise(abs(nomv - 5) < 1e-9, "nominal level 5% (the 'pass_pct > 5' filters)")
  V2C <- "pk_model=='k2020' & scenario=='V2_up_080'"

  # ---- 전제: 체중 구간·아토피 ----------------------------------------------------------------------------------------------------------
  wb <- rows(WB, "model %in% c('a','b')")
  premise(all(wb[model == "a", model_variant] == "base") && all(wb[model == "b", model_variant] == "struct2020"), "weight bands: model a = 2016 model (base), model b = 2020 model (Model 1)")
  for (m in c("a", "b")) { x <- wb[model == m][order(a11_lo(band))]
    premise(nrow(x) == 6 && all(diff(x$reliable_rsq_extrap_pct) < 0) && all(diff(x$extrap_true_p95) > 0), sprintf("heavier band: lower set (i) reliability and lower 5th-percentile coverage (model %s)", m)) }
  premise(all(wb$coverage_lt80_pct <= 0.01), "coverage below 80% at most 0.01% in every weight band")
  s24 <- paste(readLines(proj_path("scripts", "24_weight_generalization.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  cs1 <- .read("config/prereg_20260926.yaml")$section4$criteria_sets$i
  premise(grepl("reliable_rsq_extrap_pct = 100 * mean(x$lambda_ok & x$adj_r2 >= 0.80 & x$pct_extrap <= 20", s24, fixed = TRUE) && cs1$adj_r2_min == 0.8 && cs1$extrap_max_pct == 20 && is.null(cs1$span_ratio_min),
          "weight-band reliability column is criteria set (i) (adjusted R-squared 0.80, extrapolation 20%, no span condition)")
  bands <- c("below 60", "60-90", "above 90-100", "above 100")
  cb <- rows(ABN, "distribution=='primary' & set=='iii'")
  premise(setequal(unique(cb$variant), c("base", "struct2020", "k2016_bmi_vc0817")), "three model variants")
  for (vr in unique(cb$variant)) { y <- cb[variant == vr][match(bands, band), fail_pct]
    premise(length(y) == 4 && !anyNA(y) && all(diff(y) > 0), sprintf("atopic: set (iii) failing share rises across weight bands (%s)", vr)) }
  ap <- rows(AIP, "distribution=='primary'")
  premise(all(ap$coverage_lt80_pct < 0.05), "atopic: coverage below 80% under 0.05% in every band and variant")
  premise(nrow(unique(ap[band %in% bands, .(band, round(wt_mean, 9))])) == 4, "atopic: band mean weight identical over variants (same population)")
  acf <- rows(ACF, "distribution=='primary' & set=='iii'"); premise(nrow(acf) == 3, "three model variants in the coverage file")
  # 무거운 구간 표시: 모델 개발 자료 범위 밖일 수 있는 구간(dev_range_note가 있는 구간)
  fl_ <- wb[model == "a" & dev_range_note != ""]; premise(nrow(fl_) == 1 && fl_$band == wb[model == "a"][which.max(a11_lo(band)), band], "only the heaviest band is flagged as possibly outside the model development data")
  premise(all(wb[model == "b" & dev_range_note != "", band] == fl_$band), "the same band is flagged for the 2020 model")
  k130v <- a11_lo(fl_$band)
  k130 <- dderived("lower bound of the weight band flagged as possibly outside the model development data (kg)", WB, "model=='a' & dev_range_note != '' :: lower bound of band", k130v, fnum(k130v, 0))
  pa <- .read("config/population_atopic.yaml")
  premise(identical(pa$distributions$primary$dist, "lognormal") && identical(pa$status, "assumption"), "atopic weight distribution is a lognormal placeholder (status assumption)")

  # ---- 전제: 가정 변경 요약(결과보고 덱 S21) ---------------------------------------------------------------------------------------------------
  premise(nrow(rows(T1, sprintf("config=='P2' & pass_pct > 5 & !(%s)", V2C))) == 0, "AUClast + Cmax above 5% nowhere but the 2020 V2 cell, all analysis models")
  for (am in c("M0", "M1", "M2")) { r <- rows(T1, sprintf("analysis_model=='%s' & config=='P2' & class!='conservative'", am))
    premise(nrow(r) == 1 && r$pk_model == "k2020" && r$scenario == "V2_up_080", paste("the only non-conservative AUClast + Cmax cell is the 2020 V2 cell,", am)) }
  cls <- vapply(c("M0", "M1", "M2"), function(am) row1(T1, sprintf("%s & analysis_model=='%s' & config=='P2'", V2C, am))$class, "")
  premise(identical(unname(cls), c("nominal", "exceeding", "exceeding")), "V2 cell class: nominal under M0, exceeding under M1 and M2")
  premise(row1(T1, sprintf("analysis_model=='M1' & config=='P2' & %s", V2C))$pass_pct == max(rows(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020'")$pass_pct), "the 2020 maximum is the V2 cell")
  CUR <- c(vmax080 = "vmax080_both", vmax125 = "vmax125_both", km05 = "km05_both", km10 = "km10_both")
  cur <- rows(CV, "case %in% c('vmax080','vmax125','km05','km10')")
  premise(nrow(cur) == 4 && all(mapply(function(cs, src) grepl(paste0("variant ", CUR[[cs]], " "), src, fixed = TRUE), cur$case, cur$source)), "curve-shape cases are the scripts/21 variants")
  premise(all(cur$pct_lt80 == 0), "no subject below 80% in the curve-shape cases")
  lt <- rows(LT)
  premise(nrow(rows(LT, "config=='P2' & pass_pct > 5")) == 0, "AUClast + Cmax at or below 5% in every LLOQ cell")
  premise(length(unique(lt$scenario)) == 3 && all(lt$n_trials == 5000) && !"pk_model" %in% names(lt), "LLOQ trials: 3 boundary cells, 5,000 trials, 2016 model only (no PK model column)")
  premise(!any(grepl("iii", unique(lt$config))) && "G2_Aii" %in% lt$config, "LLOQ trials: AUCinf rule A only with set (ii) (no set (iii))")
  premise(nrow(rows(LT, "config=='G2_Aii' & resid=='fixed' & model=='M1' & pass_pct <= 5")) == 0, "AUCinf (set ii) + Cmax above 5% in every LLOQ cell (M1)")
  rf <- list.files(proj_path("results"), pattern = "resid12", recursive = TRUE)
  premise(!any(grepl("^(oc|oc_models|criteria|lloq)/|type1", rf)), "no boundary type I error was computed with the 12% residual")
  premise(row1(TRS, "variant=='resid12' & set=='iii'")$fail_pct < row1(TRS, "variant=='k2016' & set=='iii'")$fail_pct && row1(TRS, "variant=='resid12' & set=='iii'")$fail_pct > 20,
          "smaller residual: fewer subjects without a reliable AUCinf, but still more than 20% (notes: failures remain)")
  premise(setequal(unique(lt$scenario), c("F_down_080", "Vmax_down_125", "Vmax_up_080")), "LLOQ boundary cells: absorbed amount down, target-mediated elimination down and up (notes)")
  for (m in c("k2016", "k2020")) premise(nrow(rows(CG, sprintf("analysis_model=='M1' & config=='G2_A_iii' & pk_model=='%s' & pass_pct > 5", m))) >= 1, paste("AUCinf (set iii) + Cmax above 5% in both models (notes: exceedance remains),", m))

  # ---- 제목 ----------------------------------------------------------------------------------------------------------------------------
  # 제목: 그린 자료(두 모델)의 신뢰 기준 미달 최대와 AUClast/AUCinf 80% 미만 비율 최대(체중 구간과 아토피 구간)
  fpl <- rows(ABN, "distribution=='primary' & set=='iii' & variant %in% c('base','struct2020') & band %in% c('below 60','60-90','above 90-100','above 100')")
  wmax_a <- max(100 - wb$reliable_rsq_extrap_pct); hi_v <- max(fpl$fail_pct); premise(hi_v > wmax_a, "the largest plotted failing share is in panel (b) (title)")
  cpl <- rows(AIP, "distribution=='primary' & variant %in% c('base','struct2020') & band %in% c('below 60','60-90','above 90-100','above 100')")
  lt_v <- max(c(wb$coverage_lt80_pct, cpl$coverage_lt80_pct))
  hi <- dderived("largest plotted share without a reliable AUCinf (set iii, atopic weight bands, two models)", ABN,
                 "distribution=='primary' & set=='iii' & variant in (base, struct2020) & band in 4 bands :: max(fail_pct)", hi_v, paste0(fnum(hi_v, 1), "%"))
  dsrc("largest share below 80% AUClast/AUCinf, plotted weight and atopic bands", AIP, "(table)")
  lt80 <- dderived("largest share below 80% AUClast/AUCinf over the plotted weight bands (models a, b) and atopic bands (two models), rounded up", WB,
                 "max(coverage_lt80_pct) over weight_bands_B0_abcd.csv [model in a, b] and atopic_individual_pillar1.csv [primary, base and struct2020, 4 bands]", lt_v, paste0(fnum(cl(lt_v, 2), 2), "%"))
  th80 <- dderived("AUClast/AUCinf threshold in column name coverage_lt80_pct (percent)", WB, "column name coverage_lt80_pct", 80, "80")
  zmin <- 100 * min(acf$window_fail_min)
  z <- dderived("smallest AUClast/AUCinf (true) among patients without a reliable AUCinf (set iii), three model variants, rounded down to 0.1 (%)", ACF,
                "distribution=='primary' & set=='iii' :: floor(1000 * min(window_fail_min)) / 10", zmin, paste0(fnum(fl(zmin, 1), 1), "%"))
  hib <- fpl[which.max(fail_pct)]; premise(hib$band == "above 100" && hib$variant == "base", "the largest plotted set (iii) failing share is the above-100 kg band of the assumed atopic distribution, 2016 model (title)")
  bk <- unlist(pa$bands_kg); premise(length(bk) == 3 && hib$band == sprintf("above %s", format(bk[3])), "the heaviest atopic band starts at the last configured band limit (title)")
  premise(any(grepl("접근할 수 없었다", readLines(proj_path("config", "population_atopic.yaml"), warn = FALSE, encoding = "UTF-8"), fixed = TRUE)) && grepl("^placeholder", pa$distributions$primary$source),
          "atopic weight distribution: placeholder from literature summaries; phase 3 raw data not accessible (caption)")
  k100t <- dcfg("population_atopic.yaml", "bands_kg", "lower weight bound of the heaviest atopic band (kg)", function(x) fnum(as.numeric(x[[3]]), 0))
  y0 <- core_title(tx("A11.title", list(hi = hi, lt = lt80, th = th80, k100 = k100t)), tx("A11.kicker"))

  # ---- 근거 아님 띠 -----------------------------------------------------------------------------------------------------------------------
  wt <- f_wt_range()
  FW <- 6.1; XT <- GEO$ML + FW + 0.25; TW <- GEO$W - GEO$MR - XT
  ne <- tx("A11.not_evidence", list(wt = wt))
  nh <- est_height(ne, FW - 0.2, SZ$body, 0) + 0.1
  deck_text(ne, c(GEO$ML, y0, FW, nh), size = SZ$body, bg = PAL$tint_grey, geom = "roundRect", label = "text_not_evidence", gap_pt = 0)
  yc <- y0 + nh + 0.1

  # ---- 캡션(아래) -------------------------------------------------------------------------------------------------------------------------
  cap <- tx("A11.caption", list(k = k130, nom = nom, r2i = f_set("i", "r2"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), lt = lt80, th = th80,
                                m = dcfg("population_atopic.yaml", c("distributions", "primary", "mean"), "atopic placeholder distribution, mean (kg)", num_fmt(0)),
                                s = dcfg("population_atopic.yaml", c("distributions", "primary", "sd"), "atopic placeholder distribution, SD (kg)", num_fmt(0))))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)

  # ---- 오른쪽 표: 가정 변경 요약 ---------------------------------------------------------------------------------------------------------------
  pmax_ <- function(m) dext(T1, sprintf("analysis_model=='M1' & config=='P2' & pk_model=='%s'", m), "pass_pct", max, 2, "%", sprintf("largest boundary type I error, AUClast + Cmax, M1, %s", m))
  p05c <- 100 * min(cur$p05)
  R <- list(
    list(p16 = pmax_("k2016"), p20 = pmax_("k2020"), n = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii'", "boundary cells, AUCinf (set iii) + Cmax, M1"),
         k = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii' & pass_pct > 5", "AUCinf (set iii) + Cmax cells above 5%, M1"), nom = nom),
    list(v = a11_mult(c("vmax080_both", "vmax125_both"), "Vmax", "Vmax multipliers, both arms, curve-shape cases"),
         km = a11_mult(c("km05_both", "km10_both"), "Km", "Km multipliers, both arms, curve-shape cases"),
         th = dderived("AUClast/AUCinf threshold in column name pct_lt80 (percent)", CV, "column name pct_lt80", 80, "80"),
         lt = drange(CV, "case %in% c('vmax080','vmax125','km05','km10')", "pct_lt80", 0, "%", "share below 80% AUClast/AUCinf, curve-shape cases"),
         p05 = dderived("smallest 5th percentile of AUClast/AUCinf over the curve-shape cases, rounded down", CV, "case %in% c('vmax080','vmax125','km05','km10') :: floor(min(p05) x 1000) / 10", p05c, paste0(fnum(fl(p05c, 1), 1), "%"))),
    list(lg = f_lloq_grid(), nc = dderived("boundary cells in the LLOQ trials", LT, "length(unique(scenario))", length(unique(lt$scenario)), as.character(length(unique(lt$scenario)))),
         p1 = drange(LT, "model=='M1' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid, AUClast + Cmax, M1"),
         nt = dint(LT, "model=='M1' & resid=='fixed' & config=='P2' & scenario=='F_down_080' & abs(lloq - 0.078) < 1e-9", "n_trials", "LLOQ trials per cell"),
         g1 = drange(LT, "model=='M1' & resid=='fixed' & config=='G2_Aii'", "pass_pct", 1, "%", "LLOQ grid, AUCinf (set ii, rule A) + Cmax, M1")),
    list(s24 = dv(TRS, "variant=='k2016' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model"), s12 = dv(TRS, "variant=='resid12' & set=='iii'", "sigma_prop_pct", 0, "%", "proportional residual, variant"),
         f24 = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "without a reliable AUCinf (set iii), residual 24.2%"), f12 = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "without a reliable AUCinf (set iii), residual 12%")),
    list(nom = nom, c = dcount(T1, "analysis_model=='M1' & config=='P2' & class=='conservative'", "AUClast + Cmax cells conservative, M1")))
  premise(all(vapply(c("M0", "M2"), function(am) nrow(rows(T1, sprintf("analysis_model=='%s' & config=='P2' & class=='conservative'", am))), 1L) == as.integer(R[[5]]$c)), "same number of conservative cells under M0, M1 and M2")
  H <- L$table
  df <- data.frame(a = vapply(1:5, function(i) fill(H$cond[[i]], R[[i]]), ""), b = vapply(1:5, function(i) fill(H$result[[i]], R[[i]]), ""), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- unlist(H$head)
  deck_table(df, box = c(XT, y0, TW, capy - 0.1 - y0), widths = c(1.7, TW - 1.7), size = 14, align_num = FALSE, label = "table_robust")

  # ---- 왼쪽 그림 -------------------------------------------------------------------------------------------------------------------------
  nb <- dint(WB, "model=='a' & band=='40-60'", "n", "subjects per weight band"); premise(all(wb$n == wb$n[1]), "same number of subjects in every weight band")
  F <- L$fig
  p <- a11_fig(WB, ABN, AIP, F, ML, k130v)
  deck_figure(p, "a11_weight_atopic", c(GEO$ML, yc, FW, capy - 0.1 - yc), src = c(WB, ABN, AIP))

  # ---- 노트 ------------------------------------------------------------------------------------------------------------------------------
  wv <- function(m, b, col, it) dv(WB, sprintf("model=='%s' & band=='%s'", m, b), col, 1, "%", it)
  fi <- function(m, b) { r <- row1(WB, sprintf("model=='%s' & band=='%s'", m, b)); x <- 100 - r$reliable_rsq_extrap_pct
    dderived(sprintf("weight band %s, model %s: failing set (i) = 100 - reliable_rsq_extrap_pct", b, m), WB, sprintf("model=='%s' & band=='%s' :: 100 - reliable_rsq_extrap_pct", m, b), x, paste0(fnum(x, 1), "%")) }
  c5 <- function(m, b) { r <- row1(WB, sprintf("model=='%s' & band=='%s'", m, b)); x <- 100 - r$extrap_true_p95
    dderived(sprintf("weight band %s, model %s: 5th percentile AUClast/AUCinf = 100 - extrap_true_p95", b, m), WB, sprintf("model=='%s' & band=='%s' :: 100 - extrap_true_p95", m, b), x, paste0(fnum(x, 1), "%")) }
  at <- function(vr, b, it) dv(ABN, sprintf("variant=='%s' & distribution=='primary' & band=='%s' & set=='iii'", vr, b), "fail_pct", 1, "%", it)
  AWS <- "distribution=='primary' & source=='simulated (20,000 subjects)'"
  aw <- function(col, it) dv(AW, AWS, col, 1, "%", sprintf("primary weight distribution, simulated, %s", it))
  k60 <- dderived("weight threshold in column name below_60 (kg)", AW, "column name below_60", 60, "60"); k100 <- dderived("weight threshold in column name above_100 (kg)", AW, "column name above_100", 100, "100")
  wtr <- as.numeric(unlist(.read("config/trial_design.yaml")$weight$inclusion_kg))
  premise(all(c("below_60", "above_90", "above_100", "outside_60_90") %in% names(rows(AW))) && identical(wtr, c(60, 90)), "column thresholds below_60 / above_90 are the study weight range")
  od <- rbindlist(lapply(c("a", "b", "d"), function(m) rows(sprintf("weight_generalization/obese_trials_dropout_%s.csv", m))))
  dsrc("heavier dropouts in trials of the heavier population", sprintf("weight_generalization/obese_trials_dropout_%s.csv", c("b", "d")), "(table)")
  dwr <- range(od$wt_dropout_mean - od$wt_retained_mean); premise(all(dwr > 0), "failing subjects heavier than retained in every scenario and arm")
  nr_ <- .read("config/nca_rules.yaml")$standard$reliability; cs2 <- .read("config/prereg_20260926.yaml")$section4$criteria_sets$ii
  premise(nr_$adj_r2_min == cs2$adj_r2_min && nr_$extrap_max_pct == cs2$extrap_max_pct && nr_$span_ratio_min == cs2$span_ratio_min, "engine reliability flag used for the dropout comparison equals criteria set (ii)")
  bias <- function(ep) drange(AB, sprintf("analysis_model=='M0' & endpoint=='%s'", ep), "bias_pct", 2, "%", sprintf("atopic trials, bias of %s GMR, M0", ep))
  pf <- rows(TPF, "set=='iii'")
  premise(abs(row1(ABN, "variant=='base' & distribution=='primary' & band=='all' & set=='iii'")$fail_pct - pf[pk_model == "k2016", fail_pct]) < 1 &&
          abs(row1(ABN, "variant=='struct2020' & distribution=='primary' & band=='all' & set=='iii'")$fail_pct - pf[pk_model == "k2020", fail_pct]) < 1,
          "atopic failing share (set iii) within 1 point of the study population in both models (notes: similar)")
  deck_notes(tx("A11.notes", list(
    wt = wt, nb = nb, k130 = k130, nom = nom, z = z, hi = hi, lt = lt80, th = th80, k100t = k100t, p1n = R[[3]]$p1,
    a1 = fi("a", "40-60"), a6 = fi("a", "130-150"), b1 = fi("b", "40-60"), b6 = fi("b", "130-150"),
    ca1 = c5("a", "40-60"), ca6 = c5("a", "130-150"), cb1 = c5("b", "40-60"), cb6 = c5("b", "130-150"),
    wlt = dext(WB, "model %in% c('a','b')", "coverage_lt80_pct", max, 2, "%", "largest share below 80% AUClast/AUCinf over weight bands, two models"),
    dw = dderived("dropouts heavier than retained subjects (kg)", "weight_generalization/obese_trials_dropout_a.csv", "range of wt_dropout_mean - wt_retained_mean over models a, b, d", dwr, rng_fmt(dwr[1], dwr[2], 1)),
    m = dcfg("population_atopic.yaml", c("distributions", "primary", "mean"), "atopic placeholder distribution, mean (kg)", num_fmt(0)),
    s = dcfg("population_atopic.yaml", c("distributions", "primary", "sd"), "atopic placeholder distribution, SD (kg)", num_fmt(0)),
    tr = dcfg("population_atopic.yaml", c("distributions", "primary", "trunc"), "atopic placeholder distribution, truncation (kg)", function(x) rng_fmt(x[1], x[2], 0)),
    out = aw("outside_60_90", "share outside 60 to 90 kg"), a100 = aw("above_100", "share above 100 kg"), k60 = k60, k100 = k100,
    f16 = at("base", "all", "atopic, 2016 model, set iii failing, all"), f20 = at("struct2020", "all", "atopic, 2020 model, set iii failing, all"),
    f3 = drange(ABN, "distribution=='primary' & band=='all' & set=='iii'", "fail_pct", 1, "%", "atopic, three variants, set iii failing, all"),
    t16 = dv(TPF, "pk_model=='k2016' & set=='iii'", "fail_pct", 1, "%", "study population, set iii failing, 2016"), t20 = dv(TPF, "pk_model=='k2020' & set=='iii'", "fail_pct", 1, "%", "study population, set iii failing, 2020"),
    lo16 = at("base", "below 60", "atopic, 2016 model, set iii failing, below 60 kg"), hi16 = at("base", "above 100", "atopic, 2016 model, set iii failing, above 100 kg"),
    lo20 = at("struct2020", "below 60", "atopic, 2020 model, set iii failing, below 60 kg"), hi20 = at("struct2020", "above 100", "atopic, 2020 model, set iii failing, above 100 kg"),
    sens = drange(ABN, "distribution!='primary' & band=='all' & set=='iii'", "fail_pct", 1, "%", "sensitivity weight distributions, set iii failing, three variants"),
    wmed = drange(ACF, "distribution=='primary' & set=='iii'", "window_fail_median", 1, "%", "set iii failing patients: median AUClast/AUCinf, three variants", scale = 100),
    w05 = drange(ACF, "distribution=='primary' & set=='iii'", "window_fail_p05", 1, "%", "set iii failing patients: 5th percentile AUClast/AUCinf, three variants", scale = 100),
    wm16 = dv(ACF, "variant=='base' & distribution=='primary' & set=='iii'", "window_fail_min", 1, "%", "set iii failing patients: minimum AUClast/AUCinf, 2016 model", scale = 100),
    wm20 = dv(ACF, "variant=='struct2020' & distribution=='primary' & set=='iii'", "window_fail_min", 1, "%", "set iii failing patients: minimum AUClast/AUCinf, 2020 model", scale = 100),
    alt = dext(AIP, "distribution=='primary' & band=='all'", "coverage_lt80_pct", max, 3, "%", "atopic, share below 80% AUClast/AUCinf, three variants"),
    bl = bias("AUClast"), bai = bias("AUCinf_Aiii"),
    ag = drange(ACC, "analysis_model=='M0' & config=='A_iii'", "agree_pct", 1, "%", "atopic trials, AUClast + Cmax vs AUCinf (set iii, rule A) + Cmax decision agreement, M0"),
    ad = drange(AAF, "scenario=='VM125' & set=='iii'", "diff_mean", 2, "", "atopic trials, test minus reference failing set iii, VM125 (points)"),
    ntr = dcfg("prereg_20260926.yaml", c("section5", "trial", "trials"), "atopic trials per scenario (pre-registered)", function(x) fint(as.numeric(x))),
    m0 = dv(T1, sprintf("%s & analysis_model=='M0' & config=='P2'", V2C), "pass_pct", 2, "%", "AUClast + Cmax, 2020 V2 cell, M0"),
    m1 = dci(T1, sprintf("%s & analysis_model=='M1' & config=='P2'", V2C), "pass_pct", "lo", "hi", 2, "%", "AUClast + Cmax, 2020 V2 cell, M1"),
    m2 = dv(T1, sprintf("%s & analysis_model=='M2' & config=='P2'", V2C), "pass_pct", 2, "%", "AUClast + Cmax, 2020 V2 cell, M2"),
    n20 = f_reps("ext"), n10 = f_reps("boundary"), npm = dcount(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2016'", "boundary cells per model"),
    imax = dext(CG, "analysis_model=='M1' & config=='G2_A_iii'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii) + Cmax, M1"),
    p0 = drange(LT, "model=='M0' & resid=='fixed' & config=='P2'", "pass_pct", 2, "%", "LLOQ grid, AUClast + Cmax, M0"),
    nt = dint(LT, "model=='M1' & resid=='fixed' & config=='P2' & scenario=='F_down_080' & abs(lloq - 0.078) < 1e-9", "n_trials", "LLOQ trials per cell"),
    i24 = dv(TRS, "variant=='k2016' & set=='i'", "fail_pct", 1, "%", "set (i) failing, residual 24.2%"), i12 = dv(TRS, "variant=='resid12' & set=='i'", "fail_pct", 1, "%", "set (i) failing, residual 12%"),
    r2i = f_set("i", "r2"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"))))
  deck_end()
}
