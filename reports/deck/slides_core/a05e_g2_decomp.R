# A5e 가이드라인 기본 구성의 1 쪽 치우침 분해(v1.2, 지시 2026-09-29 "S9 재설계" §5 "기전 패널(외삽 대 탈락 편향)", D-065):
# 경계 16칸의 1 쪽 치우침(%): AUClast, AUCinf(B = λz 산출 전원: 외삽), AUCinf(A = 세트 (iii) 충족자만: 외삽 + 기준 미달자 제외), B → A 화살표 = 제외(탈락) 몫.
#   results/oc_curves/oc_bias_relation.csv(bias_pct = 로그 기하평균비 평균의 참 AUCinf 비 대비 1 쪽 치우침, %). 2종 오류 분해는 별첨 A5f.
slide_A5e <- function() {
  BR <- "oc_curves/oc_bias_relation.csv"
  deck_slide("A5e", tag = "sim")
  L <- DK$txt$A5e; ML <- DK$txt$common$models_short; MS <- names(ML)
  br <- rows(BR, "endpoint %in% c('AUClast','AUCinf_B','AUCinf_Aiii')"); premise(nrow(br) == 16 * 3, "16 boundary cells x three endpoints")
  for (m in MS) { v <- br[pk_model == m & scenario == "Vmax_up_080"]
    premise(v[endpoint == "AUCinf_B", bias_pct] > 0 && v[endpoint == "AUCinf_Aiii", bias_pct] > v[endpoint == "AUCinf_B", bias_pct], sprintf("%s Vmax up: rule B biased toward 1 and rule A more (title)", m)) }
  vb <- function(m, ep) dv(BR, sprintf("pk_model=='%s' & scenario=='Vmax_up_080' & endpoint=='%s'", m, ep), "bias_pct", 1, "%", sprintf("bias toward 1, %s, Vmax up, %s (%%)", ep, m))
  f <- list(b20 = vb("k2020", "AUCinf_B"), b16 = vb("k2016", "AUCinf_B"), a20 = vb("k2020", "AUCinf_Aiii"), a16 = vb("k2016", "AUCinf_Aiii"))
  y0 <- core_title(tx("A5e.title", f), tx("A5e.kicker"))

  # ---- 왼쪽: 치우침 분해 ----
  SC <- c("F_down_080", "F_up_125", "ke_up_080", "ke_down_125", "Vmax_up_080", "Vmax_down_125", "V2_up_080", "ka_down_080")
  premise(setequal(unique(br$scenario), SC), "eight boundary scenarios per model")
  LF <- DK$txt$A5d$fig                                                  # 칸 이름은 별첨 A5d 점도표와 같다
  br[, row := length(SC) + 1 - match(scenario, SC)][, y := row + c(AUClast = 0.2, AUCinf_B = 0, AUCinf_Aiii = -0.2)[endpoint]]
  br[, mod := factor(unlist(ML[pk_model]), levels = unlist(ML))]
  EP <- c("AUClast", "AUCinf_B", "AUCinf_Aiii"); le <- unlist(L$fig$ep[EP]); br[, epf := factor(endpoint, levels = EP, labels = le)]
  lab_y <- unique(br[, .(row, mechanism, direction, target)])[order(row)]
  lab_y[, lab := vapply(seq_len(.N), function(i) fill(LF$cell, list(mech = LF$mech[[paste(mechanism[i], direction[i], sep = "_")]], tgt = fnum(target[i], 2))), "")]
  ar <- dcast(br[endpoint %in% c("AUCinf_B", "AUCinf_Aiii")], mod + row ~ endpoint, value.var = "bias_pct")
  p1 <- ggplot(br, aes(bias_pct, y)) +
    geom_hline(yintercept = seq(1.5, length(SC) - 0.5, 1), colour = PAL$grid, linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.6) +
    geom_segment(data = ar, aes(x = AUCinf_B, xend = AUCinf_Aiii, y = row, yend = row - 0.2), inherit.aes = FALSE, colour = PAL$orange, linewidth = 0.7,
                 arrow = grid::arrow(length = grid::unit(5, "pt"), type = "closed")) +
    geom_point(aes(colour = epf, shape = epf, fill = epf), size = 3.0, stroke = 1.0) +
    facet_wrap(~ mod, nrow = 1) +
    scale_colour_manual(values = setNames(c(PAL$blue, PAL$orange, PAL$orange), le), name = NULL) +
    scale_fill_manual(values = setNames(c(PAL$blue, "white", PAL$orange), le), name = NULL) +
    scale_shape_manual(values = setNames(c(22, 23, 23), le), name = NULL) +
    scale_y_continuous(breaks = lab_y$row, labels = lab_y$lab, limits = c(0.5, length(SC) + 0.5), expand = expansion(mult = 0)) +
    scale_x_continuous(labels = function(v) paste0(v, "%")) +
    labs(x = L$fig$xlab, y = NULL) + theme_core(16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), axis.text.y = element_text(size = 14, colour = PAL$ink), panel.spacing.x = grid::unit(24, "pt"))
  dsrc("bias decomposition by boundary cell", c(BR), "(figure)")

  # ---- 본문 ----
  premise(all(br[endpoint == "AUCinf_Aiii" & mechanism %in% c("F", "Vmax"), bias_pct] > 0), "AUCinf rule A biased toward 1 at every F and Vmax boundary cell (body line 2)")
  premise(all(rows("oc_curves/oc_curves_pass.csv", "analysis_model=='M1' & kind=='boundary' & config=='G2A_iii' & mechanism %in% c('F','Vmax')")$pass_pct > 5),
          "AUCinf (A) + Cmax type I above 5% at every F and Vmax boundary cell (body line 2)")
  body <- c(tx("A5e.body1", f), tx("A5e.body2", list(nom = f_nominal())))
  cap <- tx("A5e.caption", list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), reps = f_reps("boundary")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  deck_figure(p1, "a5e_bias", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(BR))

  # ---- 노트 ----
  deck_notes(tx("A5e.notes", c(f, list(
    v2a = dv(BR, "pk_model=='k2020' & scenario=='V2_up_080' & endpoint=='AUCinf_Aiii'", "bias_pct", 1, "%", "bias toward 1, AUCinf rule A, 2020 V2 up (%)"),
    v2b = dv(BR, "pk_model=='k2020' & scenario=='V2_up_080' & endpoint=='AUCinf_B'", "bias_pct", 1, "%", "bias toward 1, AUCinf rule B, 2020 V2 up (%)"),
    ka20 = dv(BR, "pk_model=='k2020' & scenario=='ka_down_080' & endpoint=='AUCinf_Aiii'", "bias_pct", 1, "%", "bias toward 1, AUCinf rule A, 2020 ka down (%)"),
    ka16 = dv(BR, "pk_model=='k2016' & scenario=='ka_down_080' & endpoint=='AUCinf_Aiii'", "bias_pct", 1, "%", "bias toward 1, AUCinf rule A, 2016 ka down (%)"),
    lv20 = dv(BR, "pk_model=='k2020' & scenario=='Vmax_up_080' & endpoint=='AUClast'", "bias_pct", 1, "%", "bias toward 1, AUClast, 2020 Vmax up (%)"),
    lv16 = dv(BR, "pk_model=='k2016' & scenario=='Vmax_up_080' & endpoint=='AUClast'", "bias_pct", 1, "%", "bias toward 1, AUClast, 2016 Vmax up (%)"),
    reps = f_reps("boundary"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap")))))
  deck_end()
}
