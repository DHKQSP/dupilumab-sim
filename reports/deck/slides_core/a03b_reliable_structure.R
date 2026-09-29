# 별첨 A3b ① 세부: 신뢰할 수 있는 AUCinf 미달은 구조적이고 무작위가 아니다(지시 2026-09-29). 결과 파일에 있는 곳은 모두 세트 (iii)(이 덱의 정의)로 쓰고,
# 세트 (iii)이 없는 정량한계 격자(세트 (i)·(ii)만 계산)만 세트 (i)로 적고 그렇게 밝힌다. 시험 모집단(건강인, 체중 층화 배정, B0)만.
# 결과보고 덱 S11(구조적 원인), S12(비무작위 탈락), S13(처리군 의존 탈락)의 자료 논리를 따른다.
#  그림 왼쪽: 비례 잔차를 절반으로(2016 모델 구조, 추정값 -> 12%) 줄였을 때 세트 (iii)·(iv) 미달(Wilson 95% 구간; trialpop/tp_residual_sensitivity.csv).
#  그림 오른쪽: 세트 (iii) 미달 비율의 arm 간 차이(시험군 - 대조군, 시험 2,000회 평균; 평균의 95% 구간은 표식보다 좁다): 동일 제품과
#              참 AUCinf 비 0.80 경계 칸 다섯 개(trialpop/tp_arm_difference.csv). 모델은 표식 모양으로 구분한다.
#  요점: 구조적(잔차; 정량한계 lloq/lloq_individual_table.csv, 세트 (i)), 비무작위(tp_characteristics.csv, tp_strata_individual.csv), arm마다 다름(그림).
#  주의: arm 간 차이와 층 구성은 인원 수 기반이라 분석 모형(M0/M1)과 무관하다.
A3B_SC <- c("S00", "ka_down_080", "Vmax_up_080", "F_down_080", "ke_up_080", "V2_up_080")
a3b_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?%?$", p)) p else paste0("+", p)
a3b_in <- function(x) paste0("c(", paste0("'", x, "'", collapse = ","), ")")

slide_A3b <- function() {
  TRS <- "trialpop/tp_residual_sensitivity.csv"; TAD <- "trialpop/tp_arm_difference.csv"; TAF <- "trialpop/tp_arm_difference_fit.csv"
  TCH <- "trialpop/tp_characteristics.csv"; TSI <- "trialpop/tp_strata_individual.csv"; TSC <- "trialpop/tp_strata_composition.csv"
  LI <- "lloq/lloq_individual_table.csv"; LW <- "reliability/reliability_lz_window_by_schedule.csv"
  deck_slide("A3b", tag = "sim")
  L <- DK$txt$A3b; ML <- DK$txt$common$models_short; MOD <- c("k2016", "k2020")
  lq0 <- .read("config/assay.yaml")$lloq_mg_L$value; lqmin <- min(unlist(.read("config/assay.yaml")$lloq_sensitivity_mg_L))
  wl <- function(m, rv, lq) sprintf("model=='%s' & resid=='%s' & abs(lloq - %s) < 1e-9", m, rv, format(lq))

  # ---- 전제: 문장이 기대는 사실 ----
  tr <- rows(TRS); g <- function(vr, s_, col = "fail_pct") tr[variant == vr & set == s_][[col]]
  premise(nrow(tr[variant %in% c("k2016", "resid12") & set %in% c("iii", "iv")]) == 4, "residual sensitivity: sets (iii) and (iv) at two residual levels (2016 structure)")
  ratio <- g("resid12", "iii", "sigma_prop_pct") / g("k2016", "iii", "sigma_prop_pct"); premise(ratio > 0.45 && ratio < 0.55, "the proportional residual is roughly halved")
  premise(abs(g("resid12", "iv") - g("k2016", "iv")) < 1 && g("resid12", "iv", "fail_lo") < g("k2016", "iv", "fail_hi") && g("k2016", "iv", "fail_lo") < g("resid12", "iv", "fail_hi"),
          "set (iv) failure unchanged: difference below 1 point and overlapping Wilson intervals")
  premise(g("resid12", "iii", "fail_hi") < g("k2016", "iii", "fail_lo") && g("resid12", "iii", "fail_lo") > 20, "set (iii) falls with the smaller residual but more than a fifth still fails (text: remains)")
  premise(g("resid12", "iii", "adj_r2_median") > g("k2016", "iii", "adj_r2_median"), "median adjusted R-squared rises with the smaller residual (notes)")
  premise(g("resid12", "iv") - g("resid12", "iii") > g("k2016", "iv") - g("k2016", "iii"), "the span-only share rises with the smaller residual (notes)")
  lzb <- rows(LW, "variant %in% c('base','struct2020') & schedule=='B0'"); csb <- rows("cliff/cliff_summary.csv", "weight=='base'")
  premise(all(lzb[variant == "base", upper_median] + 1 < csb[model == "k2016", lloq_studyday_median]) && all(lzb[variant == "struct2020", upper_median] + 1 < csb[model == "k2020", lloq_studyday_median]),
          "median end of the lambda-z window (study day) lies before the median LLOQ day in both models (notes: window before the cliff)")
  premise(all(abs(lzb$span_median - as.numeric(.read("config/prereg_20260926.yaml")$section4$criteria_sets$iv$span_ratio_min)) < 0.5), "median span ratio close to the set (iv) threshold (notes)")
  LK <- c(k2016 = "lloq/lloq_individual_k2016.csv", k2020 = "lloq/lloq_individual_k2020.csv")
  mk <- function(m, lq, col) row1(LK[[m]], sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lq)))[[col]]
  for (m in names(LK)) premise(mk(m, lqmin, "tlast_median") > mk(m, lq0, "tlast_median") && mk(m, lqmin, "flag_rsq_pct") > mk(m, lq0, "flag_rsq_pct"),
                               sprintf("%s: at the lowest LLOQ the last quantifiable time is later and adjusted R-squared fails more often (notes)", m))
  lil <- rows(LI, sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lqmin)))
  premise(nrow(lil) == 2 && all(lil$d_rel_i_pp < 0) && all(lil$d_rel_i_hi < 0), "set (i) reliability lower at the lowest LLOQ than at the study LLOQ, both models (paired CI below 0)")
  premise(all(rows(LI, sprintf("resid=='scaled' & abs(lloq - %s) < 1e-9", format(lqmin)))$d_rel_i_hi < 0), "also with the additive residual scaled to the LLOQ (notes)")
  ch <- rows(TCH, "set=='iii'"); premise(nrow(ch) == 2 && all(ch$true_aucinf_gmr_hi < 1) && all(ch$auclast_gmr_hi < 1) && all(ch$wt_diff_lo > 0),
                                         "set (iii) failing subjects: lower true AUCinf and AUClast (upper 95% limit below 1), heavier (lower limit above 0), both models")
  si <- rows(TSI, "set=='iii'"); premise(nrow(si) == 2 && all(si$diff_lo > 0), "set (iii): heavier stratum fails more often (Newcombe lower limit above 0), both models")
  sc <- rows(TSC, "scenario=='S00' & analysis_set %in% c('auclast','iii')")
  gtc <- grep("^armdiff_abs_gt[0-9]+_pct$", names(sc), value = TRUE); premise(length(gtc) == 1, "one threshold share column in tp_strata_composition.csv")
  thr_v <- as.numeric(sub("^armdiff_abs_gt([0-9]+)_pct$", "\\1", gtc))
  premise(all(sc[analysis_set == "auclast"][[gtc]] == 0) && all(sc[analysis_set == "iii"][[gtc]] > 0) && all(abs(sc$rand_armdiff_abs_max - 100 / as.numeric(.read("config/trial_design.yaml")$n_per_arm)) < 1e-9),
          "AUClast analysis set never exceeds the threshold; set (iii) analysis set does; randomization bound = one subject")
  a <- copy(rows(TAD, sprintf("set=='iii' & scenario %%in%% %s", a3b_in(A3B_SC))))
  premise(nrow(a) == 2 * length(A3B_SC), "identical products and five boundary cells at 0.80, two models")
  premise(all(abs(a[scenario == "S00", diff_mean]) < 0.5) && all(a[scenario == "S00", diff_lo < 0 & diff_hi > 0]), "identical products: arm difference near zero, 95% interval includes zero (set iii, two models)")
  premise(all(abs(a[scenario != "S00", auc_ratio] - 0.8) < 0.005), "boundary cells: true AUCinf ratio 0.80")
  premise(all(a[scenario == "ka_down_080", diff_lo] > 0) && all(a[scenario == "V2_up_080", diff_hi] < 0), "ka decreased: test arm fails more; V2 increased: test arm fails less (set iii, two models)")
  cw <- max(a$diff_hi - a$diff_lo); premise(cw < 1, "every 95% interval of the mean arm difference is narrower than 1 point (hidden by the markers)")
  aa <- rows(TAD, "set=='iii' & scenario!='S00'"); mxv <- max(abs(aa$diff_mean))
  premise(abs(mxv - max(abs(a[scenario != "S00", diff_mean]))) < 1e-12, "the largest absolute set (iii) arm difference over all test products is in the plotted cells (title)")

  # ---- 제목 ----
  f <- list(r3b = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model with residual 12%"),
            mx = dderived("largest absolute test minus reference difference in the share failing set (iii) over the test products of both models (points)", TAD,
                          "set=='iii' & scenario!='S00' :: max(abs(diff_mean))", mxv, fnum(mxv, 1)))
  y0 <- core_title(tx("A3b.title", f), tx("A3b.kicker"))

  # ---- 요점(아래) ----
  s24 <- dv(TRS, "variant=='k2016' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model")
  s12 <- dv(TRS, "variant=='resid12' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, variant")
  s0 <- max(abs(a[scenario == "S00", diff_mean])); s0p <- fnum(s0, 1); premise(as.numeric(s0p) >= s0, "printed bound for identical products is not below the largest absolute mean difference")
  bl <- tx("A3b.bullets", list(
    s24 = s24, s12 = s12, r3a = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model (residual as estimated)"), r3b = f$r3b,
    dlq = drange(LI, sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lqmin)), "d_rel_i_pp", 1, "", "set (i) reliability decrease, lowest LLOQ vs study LLOQ, two models (points)", scale = -1),
    gmr = drange(TCH, "set=='iii'", "true_aucinf_gmr", 2, "", "true AUCinf ratio, failing to retained subjects, set (iii), two models"),
    split = f_split(), sd = drange(TSI, "set=='iii'", "diff_pp", 1, "", "stratum difference, set iii, diff_pp"),
    s0 = dderived("largest absolute mean arm difference, identical products, set (iii), two models (points)", TAD, "scenario=='S00' & set=='iii' :: max(abs(diff_mean))", s0, s0p),
    rng = { x <- range(a[scenario != "S00", diff_mean])
      dderived("mean arm difference in the share failing set (iii), boundary cells at 0.80, range over mechanisms and models (points)", TAD,
               sprintf("set=='iii' & scenario %%in%% %s & scenario!='S00' :: range(diff_mean)", a3b_in(A3B_SC)), x, sprintf("%s~%s", fnum(x[1], 1), a3b_signed(fnum(x[2], 1)))) }))
  BH <- est_height(bl, GEO$CW, 18, 6, indent = 0.3) + 0.04
  BY <- GEO$BODY_BOTTOM - BH
  deck_bullets(bl, box = c(GEO$ML, BY, GEO$CW, BH), size = 18, gap_pt = 6)

  # ---- 그림 왼쪽: 잔차 절반 ----
  FL <- L$fig
  r <- copy(tr[variant %in% c("k2016", "resid12") & set %in% c("iii", "iv")])
  rl <- c(k2016 = fill(FL$res$k2016, list(s = fnum(g("k2016", "iii", "sigma_prop_pct"), 1))), resid12 = fill(FL$res$resid12, list(s = fnum(g("resid12", "iii", "sigma_prop_pct"), 1))))
  sl <- c(iii = FL$sets$iii, iv = fill(FL$sets$iv, list(v = fnum(as.numeric(.read("config/prereg_20260926.yaml")$section4$criteria_sets$iv$span_ratio_min), 0))))
  r[, res := factor(rl[variant], levels = rl)][, sx := factor(sl[set], levels = sl)]
  pd <- position_dodge(width = 0.8)
  p1 <- ggplot(r, aes(res, fail_pct, fill = res)) +
    geom_col(width = 0.72) +
    geom_errorbar(aes(ymin = fail_lo, ymax = fail_hi), width = 0.2, colour = PAL$ink2, linewidth = 0.5) +
    geom_text(aes(y = fail_hi, label = fnum(fail_pct, 1)), vjust = -0.45, size = PT(15), family = FONT, colour = PAL$ink) +
    facet_wrap(~sx, nrow = 1) +
    scale_fill_manual(values = c(PAL$orange, ORANGE_LIGHT), guide = "none") +
    scale_y_continuous(limits = c(0, 85), breaks = seq(0, 80, 20), labels = function(v) paste0(v, "%"), expand = expansion(mult = 0)) +
    labs(x = FL$xres, y = NULL, subtitle = FL$sub1) + theme_core(16) +
    theme(panel.grid.major.x = element_blank(), strip.text = element_text(hjust = 0.5, size = 15, lineheight = 0.95), axis.text.x = element_text(size = 14, colour = PAL$ink, lineheight = 0.95),
          panel.spacing.x = grid::unit(10, "pt"))

  # ---- 그림 오른쪽: arm 간 미달 차이(세트 (iii)) ----
  a[, row := length(A3B_SC) + 1 - match(scenario, A3B_SC)][, y := row + fifelse(pk_model == "k2016", 0.15, -0.15)]
  a[, mod := factor(unlist(ML[pk_model]), levels = unlist(ML))]
  lab_y <- unique(a[, .(row, scenario, key = paste(mechanism, direction, sep = "_"))])[order(row)]
  lab_y[, lab := vapply(seq_len(.N), function(j) if (scenario[j] == "S00") FL$s00 else FL$mech[[key[j]]], "")]
  xr <- range(a$diff_mean); xl <- c(floor(xr[1] / 2) * 2 - 1, ceiling(xr[2] / 2) * 2)
  xv <- xl[2] + 1.2                                                           # 값 열(두 모델)
  vt <- a[, .(lab = paste(vapply(fnum(diff_mean[match(MOD, pk_model)], 1), a3b_signed, ""), collapse = " / ")), by = row]
  nt <- unique(a$n_trials); premise(length(nt) == 1, "same number of trials in every plotted cell")
  r80 <- mean(a[scenario != "S00", auc_ratio])
  p2 <- ggplot(a, aes(diff_mean, y)) +
    geom_hline(yintercept = length(A3B_SC) - 0.5, colour = PAL$grid, linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = diff_mean, yend = y, colour = mod), linewidth = 0.8, alpha = 0.45, show.legend = FALSE) +
    geom_point(aes(shape = mod, colour = mod), size = 3.4) +
    geom_text(data = vt, aes(x = xv, y = row, label = lab), inherit.aes = FALSE, hjust = 0, size = PT(15), family = FONT, colour = PAL$ink) +
    annotate("text", x = xv, y = length(A3B_SC) + 0.85, label = FL$vhead, hjust = 0, size = PT(14), family = FONT, colour = PAL$ink2) +
    scale_shape_manual(values = unname(CORE_MODEL_SHAPE), name = NULL) + scale_colour_manual(values = unname(CORE_MODEL_COL), name = NULL) +
    scale_y_continuous(breaks = lab_y$row, labels = lab_y$lab, limits = c(0.5, length(A3B_SC) + 1.2), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = c(xl[1], xv + 6.2), breaks = seq(-10, 15, 5), labels = function(v) ifelse(v > 0, sprintf("+%g", v), sprintf("%g", v)), expand = expansion(add = 0)) +
    coord_cartesian(clip = "off") +
    labs(x = fill(FL$xlab, list(n = fint(nt))), y = NULL, subtitle = fill(FL$sub2, list(r = fnum(r80, 2)))) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 15, colour = PAL$ink))
  FH <- BY - 0.1 - y0; WLF <- 3.75                                          # 두 그림(패널 높이를 서로 맞추지 않는다)
  deck_figure(p1, "a3b_residual", c(GEO$ML, y0, WLF, FH), src = TRS)
  deck_figure(p2, "a3b_arm_difference", c(GEO$ML + WLF + 0.2, y0, GEO$CW - WLF - 0.2, FH), src = TAD)

  # ---- 노트 ----
  lci <- function(m, rv) dci(LI, wl(m, rv, lqmin), "d_rel_i_pp", "d_rel_i_lo", "d_rel_i_hi", 2, "%p", sprintf("set (i) reliability change with paired 95%% CI, lowest LLOQ vs study LLOQ, %s, residual %s", m, rv))
  ad <- function(sc_) drange(TAD, sprintf("scenario=='%s' & set=='iii'", sc_), "diff_mean", 1, "", sprintf("test minus reference failing, %s, set iii, diff_mean", sc_))
  wci <- function(m) dci(TCH, sprintf("pk_model=='%s' & set=='iii'", m), "wt_diff_kg", "wt_diff_lo", "wt_diff_hi", 2, " kg", sprintf("failing minus retained body weight with Welch 95%% CI, set iii, %s", m))
  gci <- function(m) dci(TCH, sprintf("pk_model=='%s' & set=='iii'", m), "true_aucinf_gmr", "true_aucinf_gmr_lo", "true_aucinf_gmr_hi", 3, "", sprintf("true AUCinf ratio failing to retained with 95%% CI, set iii, %s", m))
  sci <- function(m) dci(TSI, sprintf("pk_model=='%s' & set=='iii'", m), "diff_pp", "diff_lo", "diff_hi", 1, "%p", sprintf("stratum difference with Newcombe 95%% CI, set iii, %s", m))
  q <- function(m, col, it) dv(TSI, sprintf("pk_model=='%s' & set=='iii'", m), col, 1, "%", it)
  mech <- function(m, lq_, it) dv(LK[[m]], sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lq_)), "flag_rsq_pct", 1, "%", it)
  spd <- function(vr) { x <- g(vr, "iv") - g(vr, "iii")
    dderived(sprintf("meeting set (iii) but failing span 3 = set (iv) minus set (iii) failing, %s (points)", vr), TRS, sprintf("variant=='%s' :: fail_pct[set=='iv'] - fail_pct[set=='iii']", vr), x, fnum(x, 1)) }
  deck_notes(tx("A3b.notes", list(
    wt = f_wt_range(), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), sp4 = f_set("iv", "span"), r2i = f_set("i", "r2"),
    s24 = s24, s12 = s12, r3a = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model (residual as estimated)"), r3b = f$r3b,
    r4a = dv(TRS, "variant=='k2016' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, 2016 model (residual as estimated)"),
    r4b = dv(TRS, "variant=='resid12' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, 2016 model with proportional residual 12%"),
    ci4a = dspan(TRS, "variant=='k2016' & set=='iv'", "fail_lo", "fail_hi", 1, "", "set (iv) failing, Wilson 95% interval, 2016 model"),
    ci4b = dspan(TRS, "variant=='resid12' & set=='iv'", "fail_lo", "fail_hi", 1, "", "set (iv) failing, Wilson 95% interval, residual 12%"),
    ar2a = dv(TRS, "variant=='k2016' & set=='i'", "adj_r2_median", 3, "", "median adjusted R-squared, 2016 model"),
    ar2b = dv(TRS, "variant=='resid12' & set=='i'", "adj_r2_median", 3, "", "median adjusted R-squared, residual 12%"),
    m20 = dv(TRS, "variant=='k2020' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2020 model"),
    spm = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "span_median", 2, "", "median span ratio, B0, two models"),
    lq = f_lloq(), lqm = dv(LI, wl("k2016", "fixed", lqmin), "lloq", 2, "", "lowest LLOQ in the sensitivity grid (mg/L)"),
    c16 = lci("k2016", "fixed"), c20 = lci("k2020", "fixed"), cs16 = lci("k2016", "scaled"), cs20 = lci("k2020", "scaled"),
    g16 = gci("k2016"), g20 = gci("k2020"), al = drange(TCH, "set=='iii'", "auclast_gmr", 2, "", "AUClast ratio failing to retained, set (iii), two models"),
    w16 = wci("k2016"), w20 = wci("k2020"),
    l16 = q("k2016", "fail_light_pct", "failing, lighter stratum, set iii, k2016"), h16 = q("k2016", "fail_heavy_pct", "failing, heavier stratum, set iii, k2016"),
    l20 = q("k2020", "fail_light_pct", "failing, lighter stratum, set iii, k2020"), h20 = q("k2020", "fail_heavy_pct", "failing, heavier stratum, set iii, k2020"),
    d16 = sci("k2016"), d20 = sci("k2020"), split = f_split(),
    thr = dderived(sprintf("threshold in column name %s (points)", gtc), TSC, sprintf("column name %s", gtc), thr_v, fnum(thr_v, 0)),
    gt = drange(TSC, "scenario=='S00' & analysis_set=='iii'", gtc, 1, "%", "share of identical-product trials with a between-arm heavier-stratum difference above 5 points, set (iii), two models"),
    pr = paste0(dspan(TSC, "scenario=='S00' & analysis_set=='iii'", "armdiff_p05", "armdiff_p95", 1, "", "armdiff, S00, iii: 5th to 95th percentile, two models"), "%p"),
    rb = dext(TSC, "analysis_set=='auclast'", "rand_armdiff_abs_max", max, 1, "", "largest between-arm difference in the heavier-stratum share at randomization (points)"),
    ntr = dint(TSC, "pk_model=='k2016' & scenario=='S00' & analysis_set=='auclast'", "n_trials", "regenerated trials per scenario"),
    s00 = ad("S00"), ka = ad("ka_down_080"), vm = ad("Vmax_up_080"), fd = ad("F_down_080"), ke = ad("ke_up_080"), v2 = ad("V2_up_080"),
    vmd = ad("Vmax_down_125"), fu = ad("F_up_125"), ked = ad("ke_down_125"),
    r80 = drange(TAD, sprintf("set=='iii' & scenario %%in%% %s & scenario!='S00'", a3b_in(A3B_SC)), "auc_ratio", 2, "", "true AUCinf ratio, plotted boundary cells, two models"),
    r125 = drange(TAD, "grepl('_125$', scenario) & set=='iii'", "auc_ratio", 2, "", "true AUCinf ratio, boundary cells at 1.25, two models"),
    w = dderived("widest 95% CI of the mean arm difference, set (iii), plotted cells, two models (points)", TAD, sprintf("set=='iii' & scenario %%in%% %s :: max(diff_hi - diff_lo)", a3b_in(A3B_SC)), cw, fnum(cw, 2)),
    fr = drange(TAD, "scenario=='S00' & set=='iii'", "fail_R_mean", 1, "%", "reference-arm mean share failing, set iii, two models"),
    r2f = drange(TAF, "set=='iii'", "r_squared", 2, "", "arm difference per 10% lower true ratio, set iii, r_squared"),
    spa = spd("k2016"), spb = spd("resid12"),
    f0 = mech("k2016", lq0, "adjusted R-squared below 0.80, 2016 model, study LLOQ (%)"), f1 = mech("k2016", lqmin, "adjusted R-squared below 0.80, 2016 model, lowest LLOQ (%)"),
    g0 = mech("k2020", lq0, "adjusted R-squared below 0.80, 2020 model, study LLOQ (%)"), g1 = mech("k2020", lqmin, "adjusted R-squared below 0.80, 2020 model, lowest LLOQ (%)"),
    z = drange(TSC, "scenario=='S00' & analysis_set=='auclast'", gtc, 1, "%", "share of identical-product trials with a between-arm heavier-stratum difference above 5 points, AUClast analysis set, two models"))))
  deck_end()
}
