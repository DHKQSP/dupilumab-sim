# S9 ③ 판정 안정성(지시 §2, 그림 4-3): 경계 1종 오류 점도표, 두 구성만(M1: 체중 층 포함).
#   AUClast + Cmax = P2(results/oc_models/type1_models.csv), AUCinf + Cmax = 규칙 A 세트 (iii)(신뢰할 수 있는 AUCinf 미달 제외; criteria_g2_type1.csv G2_A_iii).
# 사전 등록 7e: G2_A_iii(M1)가 있으므로 세트 (i) 대체는 쓰지 않는다(전제로 확인). 기전 이름은 짧은 한국어(문구 파일), 목표 비는 자료에서.
slide_S9 <- function() {
  T1 <- "oc_models/type1_models.csv"; CG <- "criteria/criteria_g2_type1.csv"
  deck_slide("S9", tag = "sim")
  a <- copy(rows(T1, "analysis_model=='M1' & config=='P2'")); b <- copy(rows(CG, "analysis_model=='M1' & config=='G2_A_iii'"))
  premise(nrow(a) == 16 && nrow(b) == 16 && setequal(paste(a$pk_model, a$scenario), paste(b$pk_model, b$scenario)), "16 boundary cells in both configurations (no fallback to set (i))")
  nom <- f_nominal(); nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  f <- list(max = dext(CG, "analysis_model=='M1' & config=='G2_A_iii'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii, rule A) + Cmax, M1"),
            n = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells"),
            k = dcount(CG, "analysis_model=='M1' & config=='G2_A_iii' & pass_pct > 5", "AUCinf (set iii) + Cmax cells above 5%, M1"),
            z = dcount(T1, "analysis_model=='M1' & config=='P2' & pass_pct <= 5", "AUClast + Cmax cells at or below 5%, M1"),
            p2max = dext(T1, "analysis_model=='M1' & config=='P2'", "pass_pct", max, 1, "%", "largest boundary type I error, AUClast + Cmax, M1"), nom = nom)
  premise(abs(nomv - 5) < 1e-9, "nominal level 5% (the 'pass_pct > 5' filters)")
  y0 <- core_title(tx("S9.title", f), tx("S9.kicker"))
  L <- DK$txt$S9$fig; ML <- DK$txt$common$models_short
  SC <- c("F_down_080", "F_up_125", "ke_up_080", "ke_down_125", "Vmax_up_080", "Vmax_down_125", "V2_up_080", "ka_down_080")
  d <- rbind(a[, .(pk_model, scenario, mechanism, direction, target, pass_pct, cfg = "last")], b[, .(pk_model, scenario, mechanism, direction, target, pass_pct, cfg = "inf")])
  premise(setequal(unique(d$scenario), SC), "eight boundary scenarios per model")
  d[, row := length(SC) + 1 - match(scenario, SC)][, y := row + c(last = 0.17, inf = -0.17)[cfg]]
  d[, mod := factor(unlist(ML[pk_model]), levels = unlist(ML))][, cfgf := factor(cfg, levels = c("last", "inf"), labels = c(L$cfg$last, L$cfg$inf))]
  lab_y <- unique(d[, .(row, mechanism, direction, target)])[order(row)]
  premise(nrow(lab_y) == length(SC), "one label per row (same in both models)")
  lab_y[, lab := vapply(seq_len(.N), function(i) fill(L$cell, list(mech = L$mech[[paste(mechanism[i], direction[i], sep = "_")]], tgt = fnum(target[i], 2))), "")]
  mx <- rbind(d[cfg == "inf", .SD[which.max(pass_pct)], by = .(pk_model)], d[cfg == "last" & pass_pct > nomv])   # 표시: AUCinf 구성의 모델별 최댓값, AUClast 구성의 명목 초과 칸
  xmax <- ceiling(max(d$pass_pct) / 5) * 5 + 2
  ka <- unique(d[mechanism == "ka", .(mod, row)]); premise(all(d[mechanism == "ka", pass_pct] == 0), "ka-down cells: 0% in both configurations (label)")
  p <- ggplot(d, aes(pass_pct, y, colour = cfgf, shape = cfgf)) +
    geom_hline(yintercept = seq(1.5, length(SC) - 0.5, 1), colour = PAL$grid, linewidth = 0.35) +
    geom_vline(xintercept = nomv, linetype = "22", colour = PAL$ink2, linewidth = 0.7) +
    geom_segment(aes(x = 0, xend = pass_pct, yend = y), linewidth = 0.9, alpha = 0.5) +
    geom_point(size = 3.6) +
    geom_text(data = mx, aes(label = paste0(fnum(pass_pct, 1), "%")), hjust = -0.35, size = PT(15), family = FONT, fontface = "bold", show.legend = FALSE) +
    geom_label(data = ka, aes(x = nomv + 0.4, y = row + 0.02, label = L$ka), inherit.aes = FALSE, hjust = 0, size = PT(14), family = FONT, colour = PAL$ink2,
               fill = "white", label.size = 0, label.padding = grid::unit(0.05, "lines")) +
    facet_wrap(~ mod, nrow = 1) +
    scale_colour_manual(values = setNames(c(PAL$blue, PAL$orange), c(L$cfg$last, L$cfg$inf)), name = NULL) +
    scale_shape_manual(values = setNames(c(15, 18), c(L$cfg$last, L$cfg$inf)), name = NULL) +   # 네모·마름모: 모델 표식(●/▲)과 겹치지 않게
    scale_y_continuous(breaks = lab_y$row, labels = lab_y$lab, limits = c(0.5, length(SC) + 0.5), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, 5), expand = expansion(add = c(0.3, 0))) +
    labs(x = fill(L$xlab, list(v = nom)), y = NULL) + theme_core(16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), axis.text.y = element_text(size = 14, colour = PAL$ink), panel.spacing.x = grid::unit(24, "pt"))
  f$kb <- dcount(CG, "analysis_model=='M1' & config=='G2_B' & pass_pct > 5", "AUCinf (all lambda-z estimable) + Cmax cells above 5%, M1")
  f$maxb <- dext(CG, "analysis_model=='M1' & config=='G2_B'", "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (all lambda-z estimable) + Cmax, M1")
  premise(all(rows("criteria/criteria_bias.csv", "analysis_model=='M1' & scenario=='Vmax_up_080' & endpoint %in% c('AUCinf_Aiii','AUCinf_B')")$bias_dir == "toward_1"), "Vmax up: NCA AUCinf GMR biased toward 1 (notes)")
  body <- tx("S9.body", f)
  bt <- unlist(.read("config/oc_design.yaml")$boundary_targets)
  cap <- tx("S9.caption", list(tgt = dderived("boundary true AUCinf ratios (both)", "config/oc_design.yaml", "boundary_targets :: both values", bt, paste(fnum(bt, 2), collapse = ", ")),
                               reps = f_reps("boundary"), ext = dint(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials in the extended cell")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "s9_type1", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(T1, CG))
  mxr <- function(rel, cf, pk) dv(rel, sprintf("analysis_model=='M1' & config=='%s' & pk_model=='%s' & pass_pct==max(pass_pct[analysis_model=='M1' & config=='%s' & pk_model=='%s'])", cf, pk, cf, pk), "pass_pct", 2, "%", sprintf("largest cell %s %s M1", cf, pk))
  deck_notes(tx("S9.notes", c(f, list(
    i16 = mxr(CG, "G2_A_iii", "k2016"), i20 = mxr(CG, "G2_A_iii", "k2020"),
    ci = dci(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "pass_pct", "lo", "hi", 2, "%", "AUClast + Cmax, 2020 model V2 up cell, M1"),
    ref = drange(T1, "analysis_model=='M1' & config=='AUCinf_true_only'", "pass_pct", 2, "%", "reference judged with the true AUCinf, M1"),
    i0 = dcount(CG, "analysis_model=='M0' & config=='G2_A_iii' & pass_pct > 5", "AUCinf (set iii) + Cmax cells above 5%, M0"),
    p0 = dcount(T1, "analysis_model=='M0' & config=='P2' & pass_pct > 5", "AUClast + Cmax cells above 5%, M0"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"),
    bi16 = dv("criteria/criteria_bias.csv", "analysis_model=='M1' & pk_model=='k2016' & scenario=='Vmax_up_080' & endpoint=='AUCinf_Aiii'", "bias_pct", 1, "%", "AUCinf GMR bias toward 1, set (iii) rule A, 2016 Vmax up"),
    bi20 = dv("criteria/criteria_bias.csv", "analysis_model=='M1' & pk_model=='k2020' & scenario=='Vmax_up_080' & endpoint=='AUCinf_Aiii'", "bias_pct", 1, "%", "AUCinf GMR bias toward 1, set (iii) rule A, 2020 Vmax up"),
    bb16 = dv("criteria/criteria_bias.csv", "analysis_model=='M1' & pk_model=='k2016' & scenario=='Vmax_up_080' & endpoint=='AUCinf_B'", "bias_pct", 1, "%", "AUCinf GMR bias toward 1, all estimable, 2016 Vmax up"),
    bb20 = dv("criteria/criteria_bias.csv", "analysis_model=='M1' & pk_model=='k2020' & scenario=='Vmax_up_080' & endpoint=='AUCinf_B'", "bias_pct", 1, "%", "AUCinf GMR bias toward 1, all estimable, 2020 Vmax up"),
    ke = dci(CG, "analysis_model=='M1' & config=='G2_A_iii' & pk_model=='k2016' & scenario=='ke_up_080'", "pass_pct", "lo", "hi", 2, "%", "AUCinf (set iii) + Cmax, 2016 ke up cell, M1"),
    ext = dint(T1, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials in the extended cell")))))
  deck_end()
}
