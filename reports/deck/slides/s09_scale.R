# S09 논거 ①-1 규모: 기준 세트별 미달(λz 산출 불가 / 산출되나 미달), arm당 잔류 인원, 공개 SAP 근거. 시험 모집단만.
slide_S09 <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TPA <- "trialpop/tp_retained_per_arm.csv"; TRS <- "trialpop/tp_residual_sensitivity.csv"
  deck_slide("S09", tag = "litsim")
  f <- list(iii = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models"),
            i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
            n_arm = f_n_arm(), wt = f_wt_range(), r2iii = f_set("iii", "r2"))
  deck_kicker(tx("S09.kicker")); deck_title(tx("S09.title", f))
  # 표: 세트별 정의, 미달(두 모델), λz 산출 불가, 산출되나 미달, arm당 잔류 중앙값(5~95백분위)
  setdef <- function(s_) {
    r2 <- f_set(s_, "r2"); ex <- f_set(s_, "extrap"); sp <- if (s_ %in% c("ii", "iv")) f_set(s_, "span") else NA
    fill(if (is.na(sp)) tx_raw("S09.table.def") else tx_raw("S09.table.def_span"), list(r2 = r2, ex = ex, sp = sp))
  }
  ret <- function(s_) sprintf("%s (%s)", drange(TPA, sprintf("set=='%s'", s_), "retained_median", 0, "", sprintf("retained per arm, median, set (%s)", s_)),
                               dspan(TPA, sprintf("set=='%s'", s_), "retained_p05", "retained_p95", 0, "", sprintf("retained per arm, 5th to 95th percentile, set (%s)", s_)))
  sets <- c("i", "ii", "iii", "iv")
  df <- data.frame(a = vapply(sets, function(s_) sprintf("%s  %s", DK$txt$common$sets[[s_]], setdef(s_)), ""),
                   b = vapply(sets, function(s_) drange(TPF, sprintf("set=='%s'", s_), "fail_pct", 1, "%", sprintf("failing, set (%s)", s_)), ""),
                   c = vapply(sets, function(s_) drange(TPF, sprintf("set=='%s'", s_), "est_fail_pct", 1, "%", sprintf("estimable but failing, set (%s)", s_)), ""),
                   d = vapply(sets, ret, ""), check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- tx("S09.table.head", list(n = f$n_arm))
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, 7.55, 2.75), widths = c(3.1, 1.45, 1.45, 1.55), highlight = 3)
  lz <- drange(TPF, "set=='i'", "lz_pct", 2, "%", "lambda-z not estimable (all sets), two models")
  deck_bullets(tx("S09.bullets", list(lz = lz, r2iii = f$r2iii, r24 = dv(TRS, "variant=='k2016' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model (residual 24.2%)"),
                                      r12 = dv(TRS, "variant=='resid12' & set=='iii'", "fail_pct", 1, "%", "set (iii) failing, 2016 model with residual 12%"),
                                      s24 = dv(TRS, "variant=='k2016' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, 2016 model"),
                                      s12 = dv(TRS, "variant=='resid12' & set=='iii'", "sigma_prop_pct", 1, "%", "proportional residual, variant"))),
               box = c(GEO$ML, 4.55, 7.55, 2.2), size = 16)
  # 그림: 세트별 미달(두 모델 나란히), Wilson 95% 구간. λz 산출 불가/산출되나 미달의 분해는 표에 둔다
  d <- rows(TPF)[, .(pk_model, set, fail_pct, fail_lo, fail_hi)]
  L <- DK$txt$S09$fig
  d[, set := factor(set, levels = sets, labels = unlist(DK$txt$common$sets[sets]))]; d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  p <- ggplot(d, aes(x = set, y = fail_pct, fill = model)) + geom_col(position = position_dodge(width = 0.78), width = 0.72, colour = "white", linewidth = 0.6) +
    geom_errorbar(aes(ymin = fail_lo, ymax = fail_hi), position = position_dodge(width = 0.78), width = 0.2, linewidth = 0.5, colour = PAL$ink2) +
    geom_text(aes(y = fail_hi, label = fnum(fail_pct, 1)), position = position_dodge(width = 0.78), vjust = -0.55, size = 4.3, family = FONT, colour = PAL$ink) +
    scale_fill_manual(values = unname(MODEL_COL)) + scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = L$ylab) + theme_deck(14) + theme(panel.grid.major.x = element_blank())
  deck_figure(p, "s09_failure_by_set", c(8.3, GEO$BODY_TOP, GEO$W - GEO$MR - 8.3, 5.0), src = TPF)
  deck_notes(tx("S09.notes", list(i = f$i, iii = f$iii, iv = drange(TPF, "set=='iv'", "fail_pct", 1, "%", "set (iv) failing"),
                                  r4_24 = dv(TRS, "variant=='k2016' & set=='iv'", "fail_pct", 1, "%", "set (iv), residual 24.2%"),
                                  r4_12 = dv(TRS, "variant=='resid12' & set=='iv'", "fail_pct", 1, "%", "set (iv), residual 12%"),
                                  n_ind = dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"),
                                  n_tr = dint(TPA, "pk_model=='k2016' & set=='i'", "n_trials", "identical-product trials per model"),
                                  sp4 = f_set("iv", "span"))))
  deck_end()
}
