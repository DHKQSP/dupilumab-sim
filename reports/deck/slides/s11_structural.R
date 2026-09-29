# S11 논거 ①-3 구조적 원인. 시험 모집단(건강인 60~90 kg, B0)만.
#  왼쪽: 비례 잔차를 절반으로 줄여도(2016 모델 24.2% → 12%) 세트 (iv)(span ≥ 3) 탈락은 그대로다. 세트 (ii)도 span 조건(≥ 2)이 있으나
#        잔차를 줄이면 탈락이 준다(그림에 표시). 세트 (iii) 충족·span 미달 몫은 오히려 는다(요점). results/trialpop/tp_residual_sensitivity.csv
#        (12% 변형은 2016 구조에만 있다; 대응 CI 없음, Wilson 구간만). span 미달의 원인: λz 창(B0) 대 반감기(reliability_lz_window_by_schedule.csv),
#        tp_failure_reasons.csv.
#  오른쪽: LLOQ가 낮을수록 세트 (i) AUC0-inf 신뢰 충족이 떨어진다(두 모델, 잔차 두 변형 모두). results/lloq/lloq_individual_table.csv.
#        세트 (ii)는 2016 모델·잔차 추정값에서 성립하지 않는다(표시). 세트 (iii)·(iv)는 LLOQ 격자에서 계산되지 않았다(노트).
#        기전 열(마지막 정량 시점, adjusted R² 미달): lloq_individual_k2016.csv, lloq_individual_k2020.csv.
s11_signed <- function(p) if (grepl("^-", p) || grepl("^0(\\.0+)?($|%)", p)) p else paste0("+", p)          # 양수에 "+"(인쇄값은 d* 함수 그대로)
s11_nb <- function(p) gsub("(^|~|\\s)-(?=[0-9])", "\\1\u2011", p, perl = TRUE)     # 표시층: 음수 부호에서 줄바꿈 방지(U+2011)

slide_S11 <- function() {
  TRS <- "trialpop/tp_residual_sensitivity.csv"; TPR <- "trialpop/tp_failure_reasons.csv"; LW <- "reliability/reliability_lz_window_by_schedule.csv"
  LI <- "lloq/lloq_individual_table.csv"; LK <- c(k2016 = "lloq/lloq_individual_k2016.csv", k2020 = "lloq/lloq_individual_k2020.csv")
  deck_slide("S11", tag = "sim")
  L <- DK$txt$S11
  sets <- c("i", "ii", "iii", "iv")
  lq0 <- .read("config/assay.yaml")$lloq_mg_L$value; lqmin <- min(unlist(.read("config/assay.yaml")$lloq_sensitivity_mg_L))
  wl <- function(m, rv, lq) sprintf("model=='%s' & resid=='%s' & abs(lloq - %s) < 1e-9", m, rv, format(lq))

  # ---- 전제: 문장이 기대는 사실 ----
  tr <- rows(TRS); g <- function(vr, s_, col = "fail_pct") tr[variant == vr & set == s_][[col]]
  premise(nrow(tr[variant %in% c("k2016", "resid12")]) == 8, "residual sensitivity: four sets x two residual levels (2016 structure)")
  ratio <- g("resid12", "i", "sigma_prop_pct") / g("k2016", "i", "sigma_prop_pct"); premise(ratio > 0.45 && ratio < 0.55, "the proportional residual is roughly halved (title)")
  premise(abs(g("resid12", "iv") - g("k2016", "iv")) < 1 && g("resid12", "iv", "fail_lo") < g("k2016", "iv", "fail_hi") && g("k2016", "iv", "fail_lo") < g("resid12", "iv", "fail_hi"),
          "set (iv) failure unchanged: difference below 1 point and overlapping Wilson intervals (title)")
  premise(all(sapply(c("i", "ii", "iii"), function(s_) g("resid12", s_, "fail_hi") < g("k2016", s_, "fail_lo"))), "sets (i) to (iii) fall with the smaller residual (figure)")
  premise(g("resid12", "i", "adj_r2_median") > g("k2016", "i", "adj_r2_median"), "median adjusted R-squared rises with the smaller residual (notes)")
  li <- rows(LI); premise(nrow(li) == 24, "LLOQ table: two models x two residual variants x six LLOQs")
  mono <- li[order(lloq), .(up = all(diff(reliable_no_span_pct) > 0)), by = .(model, resid)]
  premise(nrow(mono) == 4 && all(mono$up), "set (i) reliability rises with the LLOQ in both models under both residual variants (lower LLOQ, lower reliability)")
  premise(all(rows(LI, sprintf("abs(lloq - %s) < 1e-9", format(lqmin)))$d_rel_i_hi < 0), "set (i) reliability at the lowest LLOQ below the study LLOQ in all four cases (paired CI)")
  r2 <- row1(LI, wl("k2016", "fixed", lqmin)); premise(r2$d_rel_ii_lo < 0 && r2$d_rel_ii_hi > 0 && !li[model == "k2016" & resid == "fixed"][order(lloq), all(diff(reliable_pct) > 0)],
                                                       "set (ii), 2016 model, residual as estimated: no decrease at the lowest LLOQ (CI includes 0) and not monotone")
  mk <- function(m, lq, col) row1(LK[[m]], sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lq)))[[col]]
  for (m in names(LK)) premise(mk(m, lqmin, "tlast_median") > mk(m, lq0, "tlast_median") && mk(m, lqmin, "flag_rsq_pct") > mk(m, lq0, "flag_rsq_pct"),
                               sprintf("%s: at the lowest LLOQ the last quantifiable time is later and adjusted R-squared fails more often (mechanism)", m))
  lzb <- rows(LW, "variant %in% c('base','struct2020') & schedule=='B0'"); csb <- rows("cliff/cliff_summary.csv", "weight=='base'")
  premise(all(lzb[variant == "base", upper_median] + 1 < csb[model == "k2016", lloq_studyday_median]) && all(lzb[variant == "struct2020", upper_median] + 1 < csb[model == "k2020", lloq_studyday_median]),
          "median end of the lambda-z window (study day) lies before the median LLOQ day in both models (text: window before the cliff)")
  sa <- .read("config/params_variability.yaml")$residual$sigma_add$value
  premise(sa == .read("config/params_k2020_model1.yaml")$residual$sigma_add$value, "same additive residual SD in both models")

  # ---- 제목 ----
  f <- list(a = dv(TRS, "variant=='k2016' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, 2016 model (residual as estimated)"),
            b = dv(TRS, "variant=='resid12' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, 2016 model with proportional residual 12%"))
  deck_kicker(tx("S11.kicker")); deck_title(tx("S11.title", f))

  # ---- 배치 ----
  WP <- (GEO$CW - 0.35) / 2; XL <- GEO$ML; XR <- XL + WP + 0.35
  HH <- 0.42; FY <- GEO$BODY_TOP + HH + 0.02; FH <- 2.38; BY <- FY + FH + 0.08; BHt <- GEO$BODY_BOTTOM - BY
  s24 <- dv(TRS, "variant=='k2016' & set=='i'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model")
  s12 <- dv(TRS, "variant=='resid12' & set=='i'", "sigma_prop_pct", 1, "%", "proportional residual, sensitivity variant")

  # ---- 왼쪽: 잔차 절반 ----
  deck_text(tx("S11.head_left"), c(XL, GEO$BODY_TOP, WP, HH), size = 17, bold = TRUE, color = PAL$blue, label = "text_head_left")
  d <- copy(tr[variant %in% c("k2016", "resid12")])
  rl <- c(k2016 = fill(L$fig$res$k2016, list(s = fnum(g("k2016", "i", "sigma_prop_pct"), 1))), resid12 = fill(L$fig$res$resid12, list(s = fnum(g("resid12", "i", "sigma_prop_pct"), 1))))
  d[, res := factor(rl[variant], levels = rl)]
  d[, sx := factor(unlist(DK$txt$common$sets[set]), levels = unlist(DK$txt$common$sets[sets]))]
  pd <- position_dodge(width = 0.8)
  spv <- function(s_) .read("config/prereg_20260926.yaml")$section4$criteria_sets[[s_]]$span_ratio_min
  premise(is.null(spv("i")) && is.null(spv("iii")) && !is.null(spv("ii")) && !is.null(spv("iv")), "sets (ii) and (iv) carry a span condition, (i) and (iii) do not (figure labels)")
  premise(g("resid12", "ii", "fail_hi") < g("k2016", "ii", "fail_lo"), "set (ii), also with a span condition, falls with the smaller residual (figure)")
  # 두 잔차 수준은 같은 2016 모델이므로 모델 색(파랑/주황)을 쓰지 않고 파랑 계열 명암으로 구분한다. 음영은 잔차와 무관하게 그대로인 세트 (iv)
  p1 <- ggplot(d, aes(x = sx, y = fail_pct, fill = res)) + geom_blank() +
    annotate("rect", xmin = 3.52, xmax = 4.48, ymin = 0, ymax = 100, fill = PAL$tint_grey) +
    annotate("text", x = 4, y = 97, label = fill(L$fig$iv_note, list(v = fnum(spv("iv"), 0))), vjust = 1, size = 3.9, family = FONT, colour = PAL$ink2, fontface = "bold") +
    annotate("text", x = 2, y = g("k2016", "ii", "fail_hi") + 16, label = fill(L$fig$iv_note, list(v = fnum(spv("ii"), 0))), vjust = 0, size = 3.7, family = FONT, colour = PAL$ink2) +
    geom_col(position = pd, width = 0.74, colour = "white", linewidth = 0.5) +
    geom_errorbar(aes(ymin = fail_lo, ymax = fail_hi), position = pd, width = 0.2, linewidth = 0.45, colour = PAL$ink2) +
    geom_text(aes(y = fail_hi, label = fnum(fail_pct, 1)), position = pd, vjust = -0.5, size = 3.8, family = FONT, colour = PAL$ink) +
    scale_fill_manual(values = c(PAL$blue, "#9cc3ef")) + scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL, subtitle = L$fig$ylab1) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          panel.grid.major.x = element_blank(), plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))
  deck_figure(p1, "s11_residual_halved", c(XL, FY, WP, FH), src = TRS)

  spd <- function(vr) { x <- g(vr, "iv") - g(vr, "iii")
    dderived(sprintf("meeting set (iii) but failing span 3 = set (iv) minus set (iii) failing, %s (points)", vr), TRS,
             sprintf("variant=='%s' :: fail_pct[set=='iv'] - fail_pct[set=='iii']", vr), x, fnum(x, 1)) }
  premise(abs(g("k2016", "iv") - g("k2016", "iii") - row1(TPR, "pk_model=='k2016' & set=='iv' & reason=='span ratio below threshold'")$hierarchical_pct) < 1e-6,
          "set (iv) minus set (iii) equals the hierarchical span share (derivation)")
  bl <- list(s24 = s24, s12 = s12,
             r3a = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model (residual as estimated)"),
             r3b = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model with residual 12%"),
             r4a = f$a, r4b = f$b, spa = spd("k2016"), spb = spd("resid12"), sp4 = f_set("iv", "span"),
             win = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "window_median", 1, "", "median lambda-z window, B0, two models (days)"),
             hl = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "HL_median", 1, "", "median half-life from lambda-z, B0, two models (days)"),
             any = drange(TPR, "set=='iv' & reason=='span ratio below threshold'", "any_pct", 1, "%", "span ratio below 3 (any flag), two models"),
             spm = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "span_median", 2, "", "median span ratio, B0, two models"))
  premise(all(abs(rows(LW, "variant %in% c('base','struct2020') & schedule=='B0'")$span_median - as.numeric(spv("iv"))) < 0.5), "median span ratio lies just above the set (iv) threshold (bullet: close to 3)")
  premise(g("resid12", "iv") - g("resid12", "iii") > g("k2016", "iv") - g("k2016", "iii"), "the span-only share rises with the smaller residual (bullet)")
  deck_bullets(tx("S11.bullets_left", bl), box = c(XL, BY, WP, BHt), size = 16, gap_pt = 7, label = "body_left")

  # ---- 오른쪽: LLOQ 역설 ----
  deck_text(tx("S11.head_right"), c(XR, GEO$BODY_TOP, WP, HH), size = 17, bold = TRUE, color = PAL$blue, label = "text_head_right")
  q <- copy(li); q[, model_ := factor(model_lab()[model], levels = model_lab())]; q[, rv := factor(unlist(L$fig$resid[resid]), levels = unlist(L$fig$resid[c("fixed", "scaled")]))]
  q[, lx := log10(lloq)]
  grid_ <- sort(unique(q$lloq)); yl <- range(q$reliable_no_span_pct)
  # 끝 숫자: 가장 낮은 LLOQ에서 연구 LLOQ 대비 변화(%p, 자료 그대로). 겹치지 않게 세로로 벌리고 연결선(선 모양 = 잔차 변형)으로 잇는다
  yr <- c(floor(yl[1]) - 2.5, ceiling(yl[2]) + 2.5)
  el_ <- q[abs(lloq - lqmin) < 1e-9, .(model_, rv, y = reliable_no_span_pct, lab = paste0(vapply(fnum(d_rel_i_pp, 2), s11_signed, ""), "%p"))][order(y)]
  sep <- diff(yr) * 0.15; el_[, yl_ := y]
  for (k in seq_len(nrow(el_))[-1]) if (el_$yl_[k] - el_$yl_[k - 1] < sep) el_$yl_[k] <- el_$yl_[k - 1] + sep
  el_[, yl_ := yl_ - (mean(yl_) - mean(y))]; if (min(el_$yl_) < yr[1] + 0.5) el_[, yl_ := yl_ + (yr[1] + 0.5 - min(yl_))]
  lx0 <- log10(lqmin)
  p2 <- ggplot(q, aes(x = lx, y = reliable_no_span_pct, colour = model_, shape = model_, linetype = rv, group = interaction(model_, rv))) +
    geom_vline(xintercept = log10(lq0), colour = PAL$ink2, linewidth = 0.5, linetype = "22") +
    annotate("text", x = log10(lq0) + 0.035, y = yl[2] + 1.3, label = fill(L$fig$study, list(v = fnum(lq0, 3))), hjust = 1, vjust = 0, size = 3.8, family = FONT, colour = PAL$ink2) +
    geom_line(data = q[resid == "fixed"], linewidth = 1.05) + geom_point(data = q[resid == "fixed"], size = 2.6) +
    geom_line(data = q[resid == "scaled"], linewidth = 0.6) + geom_point(data = q[resid == "scaled"], size = 1.7) +          # 비례 변형은 가늘게 위에
    geom_segment(data = el_, aes(x = lx0 - 0.03, xend = lx0 - 0.11, y = y, yend = yl_, colour = model_, linetype = rv), inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = el_, aes(x = lx0 - 0.13, y = yl_, label = lab, colour = model_), inherit.aes = FALSE, hjust = 0, size = 3.6, family = FONT, show.legend = FALSE) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = unname(MODEL_SHAPE)) + scale_linetype_manual(values = c("solid", "42")) +
    scale_x_reverse(breaks = log10(grid_[abs(grid_ - lq0) > 1e-9]), labels = vapply(grid_[abs(grid_ - lq0) > 1e-9], function(x) format(x), ""),     # 연구 LLOQ는 점선 글자로
                    expand = expansion(add = c(0.06, 0.42))) +
    scale_y_continuous(limits = yr) +
    labs(x = L$fig$xlab2, y = NULL, subtitle = L$fig$ylab2) + theme_deck(12) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 1), linetype = guide_legend(order = 2, override.aes = list(colour = PAL$ink2, linewidth = 0.8))) +
    theme(legend.position = "top", legend.justification = "left", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = grid::unit(0, "pt"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"), legend.key.width = grid::unit(2.2, "lines"),
          plot.subtitle = element_text(colour = PAL$ink2, size = 12, margin = margin(0, 0, 2, 0)))
  deck_figure(p2, "s11_lloq_reliability", c(XR, FY, WP, FH), src = LI)

  lq <- f_lloq(); lqm <- dv(LI, wl("k2016", "fixed", lqmin), "lloq", 2, "", "lowest LLOQ in the sensitivity grid (mg/L)")
  dd <- function(m, rv) s11_nb(dv(LI, wl(m, rv, lqmin), "d_rel_i_pp", 2, "", sprintf("set (i) reliability change, lowest LLOQ vs study LLOQ, %s, residual %s (points)", m, rv)))
  mech <- function(m, lq_, col, d_, it, unit = "%") dv(LK[[m]], sprintf("resid=='fixed' & abs(lloq - %s) < 1e-9", format(lq_)), col, d_, unit, it)
  br <- list(lq = lq, lqm = lqm, d16 = dd("k2016", "fixed"), d20 = dd("k2020", "fixed"), s16 = dd("k2016", "scaled"), s20 = dd("k2020", "scaled"),
             sa = dcfg("params_variability.yaml", c("residual", "sigma_add", "value"), "additive residual SD (mg/L), both models", num_fmt(2)),
             f0 = mech("k2016", lq0, "flag_rsq_pct", 1, "adjusted R-squared below 0.80, 2016 model, study LLOQ (%)"),
             f1 = mech("k2016", lqmin, "flag_rsq_pct", 1, "adjusted R-squared below 0.80, 2016 model, lowest LLOQ (%)"),
             ii16 = s11_nb(s11_signed(dci(LI, wl("k2016", "fixed", lqmin), "d_rel_ii_pp", "d_rel_ii_lo", "d_rel_ii_hi", 2, "%p", "set (ii) reliability change with paired 95% CI, lowest LLOQ vs study LLOQ, 2016 model, residual as estimated"))))
  deck_bullets(tx("S11.bullets_right", br), box = c(XR, BY, WP, BHt), size = 16, gap_pt = 7, label = "body_right")

  # ---- 노트 ----
  lci <- function(m, rv) dci(LI, wl(m, rv, lqmin), "d_rel_i_pp", "d_rel_i_lo", "d_rel_i_hi", 2, "%p", sprintf("set (i) reliability change with paired 95%% CI, lowest LLOQ vs study LLOQ, %s, residual %s", m, rv))
  rel <- function(m, rv, lq_, col = "reliable_no_span_pct") dv(LI, wl(m, rv, lq_), col, 1, "%", sprintf("%s, %s, residual %s, LLOQ %s", col, m, rv, format(lq_)))
  lqmax <- max(unlist(.read("config/assay.yaml")$lloq_sensitivity_mg_L))
  deck_notes(tx("S11.notes", list(
    wt = f_wt_range(), n = dint("trialpop/tp_failure_by_set.csv", "pk_model=='k2016' & set=='i'", "n", "subjects per model, trial population"), r2i = f_set("i", "r2"),
    s24 = s24, s12 = s12, s20m = dv(TRS, "variant=='k2020' & set=='i'", "sigma_prop_pct", 1, "%", "proportional residual, 2020 model"),
    r1a = dv(TRS, "variant=='k2016' & set=='i'", "fail_pct", 1, "%", "set (i) failing, 2016 model"), r1b = dv(TRS, "variant=='resid12' & set=='i'", "fail_pct", 1, "%", "set (i) failing, residual 12%"),
    r2a = dv(TRS, "variant=='k2016' & set=='ii'", "fail_pct", 1, "%", "set (ii) failing, 2016 model"), r2b = dv(TRS, "variant=='resid12' & set=='ii'", "fail_pct", 1, "%", "set (ii) failing, residual 12%"),
    r3a = bl$r3a, r3b = bl$r3b, r4a = f$a, r4b = f$b,
    ci4a = dspan(TRS, "variant=='k2016' & set=='iv'", "fail_lo", "fail_hi", 1, "", "set (iv) failing, Wilson 95% interval, 2016 model"),
    ci4b = dspan(TRS, "variant=='resid12' & set=='iv'", "fail_lo", "fail_hi", 1, "", "set (iv) failing, Wilson 95% interval, residual 12%"),
    ar2a = dv(TRS, "variant=='k2016' & set=='i'", "adj_r2_median", 3, "", "median adjusted R-squared, 2016 model"),
    ar2b = dv(TRS, "variant=='resid12' & set=='i'", "adj_r2_median", 3, "", "median adjusted R-squared, residual 12%"),
    m20iv = dv(TRS, "variant=='k2020' & set=='iv'", "fail_pct", 1, "%", "set (iv) failing, 2020 model"), m20iii = dv(TRS, "variant=='k2020' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2020 model"),
    spa = bl$spa, spb = bl$spb, sp4 = bl$sp4, sp2 = f_set("ii", "span"),
    h4 = drange(TPR, "set=='iv' & reason=='span ratio below threshold'", "hierarchical_pct", 1, "%", "span below 3 (hierarchical), two models"),
    any = bl$any, win = bl$win, hl = bl$hl,
    w3 = drange(LW, "variant %in% c('base','struct2020') & schedule=='B0'", "min_3pt_nominal_window_days", 0, "", "shortest three-point nominal window after Day 22, B0 (days)"),
    spm = bl$spm,
    lq = lq, lqm = lqm, lqx = dv(LI, wl("k2016", "fixed", lqmax), "lloq", 1, "", "highest LLOQ in the sensitivity grid (mg/L)"),
    c16 = lci("k2016", "fixed"), c20 = lci("k2020", "fixed"), cs16 = lci("k2016", "scaled"), cs20 = lci("k2020", "scaled"),
    i16a = rel("k2016", "fixed", lqmin), i16b = rel("k2016", "fixed", lq0), i16c = rel("k2016", "fixed", lqmax),
    i20a = rel("k2020", "fixed", lqmin), i20b = rel("k2020", "fixed", lq0), i20c = rel("k2020", "fixed", lqmax),
    ii16 = br$ii16,
    ii20 = dv(LI, wl("k2020", "fixed", lqmin), "d_rel_ii_pp", 2, "%p", "set (ii) reliability change, lowest LLOQ, 2020 model, residual as estimated"),
    iis16 = dv(LI, wl("k2016", "scaled", lqmin), "d_rel_ii_pp", 2, "%p", "set (ii) reliability change, lowest LLOQ, 2016 model, residual scaled"),
    iis20 = dv(LI, wl("k2020", "scaled", lqmin), "d_rel_ii_pp", 2, "%p", "set (ii) reliability change, lowest LLOQ, 2020 model, residual scaled"),
    sa = br$sa, f0 = br$f0, f1 = br$f1,
    t0 = mech("k2016", lq0, "tlast_median", 1, "median last quantifiable time, 2016 model, study LLOQ (days after dose)", ""),
    t1 = mech("k2016", lqmin, "tlast_median", 1, "median last quantifiable time, 2016 model, lowest LLOQ (days after dose)", ""),
    q0 = mech("k2016", lq0, "quant_at_last_pct", 1, "quantifiable at the last planned sample, 2016 model, study LLOQ (%)"),
    q1 = mech("k2016", lqmin, "quant_at_last_pct", 1, "quantifiable at the last planned sample, 2016 model, lowest LLOQ (%)"),
    fs0 = mech("k2016", lq0, "flag_span_pct", 1, "span ratio below 2, 2016 model, study LLOQ (%)"),
    fs1 = mech("k2016", lqmin, "flag_span_pct", 1, "span ratio below 2, 2016 model, lowest LLOQ (%)"),
    g0 = mech("k2020", lq0, "flag_rsq_pct", 1, "adjusted R-squared below 0.80, 2020 model, study LLOQ (%)"),
    g1 = mech("k2020", lqmin, "flag_rsq_pct", 1, "adjusted R-squared below 0.80, 2020 model, lowest LLOQ (%)"),
    nli = dint(LI, wl("k2016", "fixed", lq0), "n", "subjects per model and LLOQ"))))
  deck_end()
}
