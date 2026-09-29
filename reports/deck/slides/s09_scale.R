# S09 논거 ①-1 규모: 기준 세트별 미달(λz 산출 불가 / 산출되나 미달), arm당 유지자 수, 공개 SAP 근거. 시험 모집단만.
slide_S09 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TPR <- "trialpop/tp_failure_reasons.csv"; TPA <- "trialpop/tp_retained_per_arm.csv"; TRS <- "trialpop/tp_residual_sensitivity.csv"
  deck_slide("S09", tag = "litsim")
  f <- list(iii = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models"),
            i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
            n_arm = f_n_arm(), wt = f_wt_range(), r2iii = f_set("iii", "r2"))
  deck_kicker(tx("S09.kicker")); deck_title(tx("S09.title", f))
  # 표: 세트별 정의(외삽 기준은 네 세트 공통이라 머리글에), 미달(두 모델 범위), 그중 λz 산출되나 미달, arm당 유지자 중앙값(5~95백분위)
  sets <- c("i", "ii", "iii", "iv")
  exs <- vapply(sets, function(s_) .read("config/prereg_20260926.yaml")$section4$criteria_sets[[s_]]$extrap_max_pct, 0)
  premise(length(unique(exs)) == 1, "all four criteria sets share the extrapolation limit (table header)")
  setdef <- function(s_) {
    r2 <- f_set(s_, "r2"); sp <- if (s_ %in% c("ii", "iv")) f_set(s_, "span") else NA
    fill(if (is.na(sp)) tx_raw("S09.table.def") else tx_raw("S09.table.def_span"), list(r2 = r2, sp = sp))
  }
  ret <- function(s_) sprintf("%s (%s)", drange(TPA, sprintf("set=='%s'", s_), "retained_median", 0, "", sprintf("retained per arm, median, set (%s)", s_)),
                               dspan(TPA, sprintf("set=='%s'", s_), "retained_p05", "retained_p95", 0, "", sprintf("retained per arm, 5th to 95th percentile, set (%s)", s_)))
  df <- data.frame(a = vapply(sets, function(s_) sprintf("%s  %s", DK$txt$common$sets[[s_]], setdef(s_)), ""),
                   b = vapply(sets, function(s_) drange(TPF, sprintf("set=='%s'", s_), "fail_pct", 1, "%", sprintf("failing, set (%s)", s_)), ""),
                   c = vapply(sets, function(s_) drange(TPF, sprintf("set=='%s'", s_), "est_fail_pct", 1, "%", sprintf("estimable but failing, set (%s)", s_)), ""),
                   d = vapply(sets, ret, ""), check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- tx("S09.table.head", list(n = f$n_arm, ex = f_set("i", "extrap")))
  TH <- 1.9
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, 7.75, TH), widths = c(2.9, 1.3, 1.35, 2.2), highlight = 3)
  # 표 아래 설명: 값의 범위(두 모델), 열별 수준(대상자 대 시험), span 정의
  tpf <- rows(TPF); tpa <- rows(TPA)
  premise(length(unique(tpf$n)) == 1 && length(unique(tpa$n_trials)) == 1, "same number of subjects and of identical-product trials in both models and all sets (caption)")
  n_ind <- dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"); n_tr <- dint(TPA, "pk_model=='k2016' & set=='i'", "n_trials", "identical-product trials per model")
  cap <- tx("S09.table.caption", list(n_ind = n_ind, n_tr = n_tr)); CH <- est_height(nobreak(cap), 7.75, 16)
  deck_text(cap, c(GEO$ML, GEO$BODY_TOP + TH + 0.06, 7.75, CH), size = 16, color = PAL$ink2, label = "caption_table")

  # 요점: λz 산출 불가 / 산출되나 미달의 사유, 세트 (i)이 하한인 이유와 공개 SAP, 잔차 크기
  tr_ <- rows(TPR)
  hs <- function(m, s_, r) tr_[pk_model == m & set == s_ & reason == r, hierarchical_pct]
  RR <- "adjusted R-squared below threshold"; SP <- "span ratio below threshold"
  fb <- rows(TPF)
  for (m in c("k2016", "k2020")) {
    for (s_ in c("i", "iii")) premise(hs(m, s_, RR) >= 0.9 * fb[pk_model == m & set == s_, est_fail_pct], sprintf("%s set (%s): adjusted R-squared is most of the estimable failures (bullet)", m, s_))
    for (s_ in c("ii", "iv")) premise(hs(m, s_, SP) > 0, sprintf("%s set (%s): span failures add to the other reasons (bullet)", m, s_))
    premise(fb[pk_model == m][order(fail_pct)]$set[1] == "i", sprintf("%s: set (i) fails least of the four sets (bullet: lower bound)", m))
  }
  r2s <- vapply(sets, function(s_) as.numeric(.read("config/prereg_20260926.yaml")$section4$criteria_sets[[s_]]$adj_r2_min), 0)
  premise(r2s[["i"]] == min(r2s) && is.null(.read("config/prereg_20260926.yaml")$section4$criteria_sets$i$span_ratio_min), "set (i) has the lowest adjusted R-squared and no span condition (most lenient)")
  premise(hs("k2020", "iv", SP) > hs("k2020", "iv", RR), "2020 model set (iv): span failures exceed adjusted R-squared failures (notes)")
  sap <- .read("config/prereg_20260926.yaml")$section6$s2_1_failure_by_set$public_saps
  ncts <- regmatches(sap, gregexpr("NCT[0-9]{8}", sap))[[1]]
  premise(length(ncts) == 3 && all(grepl("0\\.90", strsplit(sap, "NCT[0-9]{8}")[[1]][-1])), "three public SAPs, each with adjusted R-squared 0.90 (prereg section6)")
  premise(grepl("web-search excerpt", sub(paste0(".*", ncts[3]), "", sap)), "the last SAP was verified from a web-search excerpt (notes)")
  premise(abs(as.numeric(.read("config/prereg_20260926.yaml")$section4$criteria_sets$iii$adj_r2_min) - 0.90) < 1e-9, "set (iii) threshold equals the SAP threshold")
  nsap <- dderived("public SAPs with adjusted R-squared 0.90 (count of NCT identifiers)", "config/prereg_20260926.yaml", "section6.s2_1_failure_by_set.public_saps :: count of NCT identifiers", length(ncts), as.character(length(ncts)))
  sp4h <- drange(TPR, sprintf("set=='iv' & reason=='%s'", SP), "hierarchical_pct", 1, "%", "set (iv) span below 3 (hierarchical), two models")
  lz <- drange(TPF, "set=='i'", "lz_pct", 2, "%", "lambda-z not estimable (all sets), two models")
  BY <- GEO$BODY_TOP + TH + 0.06 + CH + 0.14
  deck_bullets(tx("S09.bullets", list(lz = lz, r2iii = f$r2iii, sp4h = sp4h, nsap = nsap, ncts = paste(ncts, collapse = ", "))),
               box = c(GEO$ML, BY, 7.75, GEO$BODY_BOTTOM - BY), size = 16, gap_pt = 8)
  # 그림: 세트별 미달(두 모델 나란히), Wilson 95% 구간. λz 산출 불가/산출되나 미달의 분해는 표에 둔다
  d <- rows(TPF)[, .(pk_model, set, fail_pct, fail_lo, fail_hi)]
  premise(max(d$fail_hi) < 75, "every Wilson upper bound below 75% (y axis ends at 80 with room for the value labels)")
  L <- DK$txt$S09$fig
  d[, set := factor(set, levels = sets, labels = unlist(DK$txt$common$sets[sets]))]; d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  p <- ggplot(d, aes(x = set, y = fail_pct, fill = model)) + geom_col(position = position_dodge(width = 0.78), width = 0.72, colour = "white", linewidth = 0.6) +
    geom_errorbar(aes(ymin = fail_lo, ymax = fail_hi), position = position_dodge(width = 0.78), width = 0.2, linewidth = 0.5, colour = PAL$ink2) +
    # 막대 값: 흰 바탕(테두리 없음)으로 가로 격자선이 글자 뒤에서 끊기게 한다
    geom_label(aes(y = fail_hi, label = fnum(fail_pct, 1)), position = position_dodge(width = 0.78), vjust = -0.3, size = 4.3, family = FONT, colour = PAL$ink,
               fill = "white", label.size = 0, label.padding = unit(0.08, "lines"), label.r = unit(0, "lines")) +
    scale_fill_manual(values = unname(MODEL_COL)) + scale_y_continuous(limits = c(0, 80), breaks = seq(0, 75, by = 25), expand = expansion(mult = c(0, 0.01))) +
    labs(x = NULL, y = NULL, subtitle = L$ylab) + theme_deck(14) +
    theme(panel.grid.major.x = element_blank(), legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(colour = PAL$ink2, size = 13, lineheight = 1.1, margin = margin(0, 0, 4, 0)))
  FX <- GEO$ML + 7.75 + 0.3
  deck_figure(p, "s09_failure_by_set", c(FX, GEO$BODY_TOP, GEO$W - GEO$MR - FX, 5.0), src = TPF)
  hr <- function(s_, r, it) drange(TPR, sprintf("set=='%s' & reason=='%s'", s_, r), "hierarchical_pct", 1, "%", it)
  deck_notes(tx("S09.notes", list(i = f$i, iii = f$iii, iv = drange(TPF, "set=='iv'", "fail_pct", 1, "%", "set (iv) failing"),
                                  r4_24 = dv(TRS, "variant=='k2016' & set=='iv'", "fail_pct", 1, "%", "set (iv), residual 24.2%"),
                                  r4_12 = dv(TRS, "variant=='resid12' & set=='iv'", "fail_pct", 1, "%", "set (iv), residual 12%"),
                                  n_ind = n_ind, n_tr = n_tr,
                                  r24 = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model (residual 24.2%)"),
                                  r12 = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model with residual 12%"),
                                  s24 = dv(TRS, "variant=='k2016' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model"),
                                  s12 = dv(TRS, "variant=='resid12' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, variant"),
                                  n = f$n_arm, nsap = nsap, sp4 = f_set("iv", "span"), sp4h = sp4h,
                                  r2h_i = hr("i", RR, "set (i) adjusted R-squared below 0.80 (hierarchical), two models"),
                                  r2h_iii = hr("iii", RR, "set (iii) adjusted R-squared below 0.90 (hierarchical), two models"),
                                  sp2h = hr("ii", SP, "set (ii) span below 2 (hierarchical), two models"),
                                  sp4h20 = dv(TPR, sprintf("pk_model=='k2020' & set=='iv' & reason=='%s'", SP), "hierarchical_pct", 1, "%", "set (iv) span below 3 (hierarchical), 2020 model"),
                                  r2h4_20 = dv(TPR, sprintf("pk_model=='k2020' & set=='iv' & reason=='%s'", RR), "hierarchical_pct", 1, "%", "set (iv) adjusted R-squared below 0.90 (hierarchical), 2020 model"))))
  deck_end()
}
