# A5d 가이드라인 기본 구성(AUCinf + Cmax)의 1종 오류 원인(v1.2, 지시 2026-09-29 "S9 재설계" §5; v1.1 본문 S9 점도표를 옮기고 관계식 그림을 더함, D-065):
# 왼쪽 = 경계 16칸 점도표(두 모델 패널, 칸마다 AUClast + Cmax, AUCinf(규칙 A) + Cmax, AUCinf(규칙 B) + Cmax; results/oc_curves/oc_curves_pass.csv, M1).
# 오른쪽 = 편향과 1종 오류의 관계: 경계 칸 x 단일 지표(AUClast, AUCinf 규칙 A·B, 참 AUCinf) 64점의 모의 통과율 대 b/SD(1 쪽 치우침 / 시험 간 표준편차)와
#   정규 근사 곡선 100 x Φ(b/SD - t)(results/oc_curves/oc_bias_relation.csv; scripts/65, 사전 등록 뒤 추가한 기술 요약).
slide_A5d <- function() {
  PF <- "oc_curves/oc_curves_pass.csv"; BR <- "oc_curves/oc_bias_relation.csv"; PD <- "oc_models/p2_decomposition_models.csv"; T1 <- "oc_models/type1_models.csv"
  deck_slide("A5d", tag = "sim")
  L <- DK$txt$A5d$fig; ML <- DK$txt$common$models_short; MS <- names(ML)
  CF <- c("P2", "G2A_iii")                                             # 점도표는 v1.1 S9와 같은 두 구성(규칙 B는 본문 1줄과 별첨 A5c)
  d <- rows(PF, "analysis_model=='M1' & kind=='boundary' & config %in% c('P2','G2A_iii')")
  premise(nrow(d) == 16 * 2, "16 boundary cells x two configurations")
  nom <- f_nominal(); nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  premise(abs(nomv - 5) < 1e-9, "nominal level 5% (the 'pass_pct > 5' filters)")
  W <- function(cf) sprintf("analysis_model=='M1' & kind=='boundary' & config=='%s'", cf)
  f <- list(max = dext(PF, W("G2A_iii"), "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (set iii, rule A) + Cmax, M1"),
            n = dcount(PF, W("P2"), "boundary cells"),
            k = dcount(PF, paste(W("G2A_iii"), "& pass_pct > 5"), "AUCinf (set iii) + Cmax cells above 5%, M1"),
            kb = dcount(PF, paste(W("G2B"), "& pass_pct > 5"), "AUCinf (all lambda-z estimable) + Cmax cells above 5%, M1"),
            maxb = dext(PF, W("G2B"), "pass_pct", max, 1, "%", "largest boundary type I error, AUCinf (all lambda-z estimable) + Cmax, M1"),
            z = dcount(PF, paste(W("P2"), "& pass_pct <= 5"), "AUClast + Cmax cells at or below 5%, M1"),
            p2max = dext(PF, W("P2"), "pass_pct", max, 1, "%", "largest boundary type I error, AUClast + Cmax, M1"), nom = nom)
  for (m in MS) premise(d[config == "G2A_iii" & pk_model == m][which.max(pass_pct), mechanism] == "Vmax", sprintf("%s: largest AUCinf + Cmax boundary cell is target-mediated elimination (title)", m))

  # ---- 관계식 자료(전제: 예측과 모의가 가깝다, 참 AUCinf는 치우침이 없다) ----
  br <- rows(BR, "TRUE"); premise(nrow(br) == 16 * 4, "64 boundary cell x single-endpoint rows")
  dmax <- max(abs(br$pred_pass_pct - br$sim_pass_pct)); premise(dmax < 1.5, "normal approximation within 1.5 percentage points of the simulated pass in every row (body line 2)")
  dcur <- max(abs(100 * pnorm(br$b_over_sd - median(br$t_q)) - br$sim_pass_pct))   # 그린 곡선(모든 점의 t 중앙값, SE = SD)과의 최대 차이(검토 M4)
  premise(all(abs(br[endpoint == "AUCinf_true", b_over_sd]) < 0.05), "true AUCinf: no bias at any boundary cell (figure)")
  premise(all(br[endpoint == "AUCinf_Aiii" & scenario == "Vmax_up_080", b_over_sd] > 0.5), "AUCinf rule A at the Vmax up cell: bias above half a standard deviation in both models (title)")
  y0 <- core_title(tx("A5d.title", f), tx("A5d.kicker"))

  # ---- 왼쪽: 점도표 ----
  SC <- c("F_down_080", "F_up_125", "ke_up_080", "ke_down_125", "Vmax_up_080", "Vmax_down_125", "V2_up_080", "ka_down_080")
  premise(setequal(unique(d$code), SC), "eight boundary scenarios per model")
  d[, row := length(SC) + 1 - match(code, SC)][, y := row + c(P2 = 0.18, G2A_iii = -0.18)[config]]
  d[, mod := factor(unlist(ML[pk_model]), levels = unlist(ML))][, cfgf := factor(config, levels = CF, labels = unlist(L$cfg[CF]))]
  lab_y <- unique(d[, .(row, mechanism, direction, target)])[order(row)]
  premise(nrow(lab_y) == length(SC), "one label per row (same in both models)")
  lab_y[, lab := vapply(seq_len(.N), function(i) fill(L$cell, list(mech = L$mech[[paste(mechanism[i], direction[i], sep = "_")]], tgt = fnum(target[i], 2))), "")]
  mx <- rbind(d[config == "G2A_iii", .SD[which.max(pass_pct)], by = .(pk_model)], d[config == "P2" & pass_pct > nomv])
  xmax <- ceiling(max(d$pass_pct) / 5) * 5 + 6
  ka <- unique(d[mechanism == "ka", .(mod, row)]); premise(all(d[mechanism == "ka", pass_pct] == 0), "ka-down cells: 0% in all three configurations (label)")
  labs3 <- unlist(L$cfg[CF])
  p1 <- ggplot(d, aes(pass_pct, y, colour = cfgf, shape = cfgf, fill = cfgf)) +
    geom_hline(yintercept = seq(1.5, length(SC) - 0.5, 1), colour = PAL$grid, linewidth = 0.35) +
    geom_vline(xintercept = nomv, linetype = "22", colour = PAL$ink2, linewidth = 0.7) +
    geom_segment(aes(x = 0, xend = pass_pct, yend = y), linewidth = 0.8, alpha = 0.45) +
    geom_point(size = 3.0, stroke = 1.0) +
    geom_text(data = mx, aes(label = paste0(fnum(pass_pct, 1), "%")), hjust = -0.3, size = PT(14), family = FONT, fontface = "bold", show.legend = FALSE) +
    geom_label(data = ka, aes(x = nomv + 0.4, y = row + 0.02, label = L$ka), inherit.aes = FALSE, hjust = 0, size = PT(14), family = FONT, colour = PAL$ink2,
               fill = "white", label.size = 0, label.padding = grid::unit(0.05, "lines")) +
    facet_wrap(~ mod, nrow = 1) +
    scale_colour_manual(values = setNames(c(PAL$blue, PAL$orange), labs3), name = NULL) +
    scale_fill_manual(values = setNames(c(PAL$blue, PAL$orange), labs3), name = NULL) +
    scale_shape_manual(values = setNames(c(22, 23), labs3), name = NULL) +
    scale_y_continuous(breaks = lab_y$row, labels = lab_y$lab, limits = c(0.5, length(SC) + 0.5), expand = expansion(mult = 0)) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, 10), expand = expansion(add = c(0.3, 0))) +
    labs(x = fill(L$xlab, list(v = nom)), y = NULL) + theme_core(16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.justification = "left",
          legend.margin = margin(0, 0, 0, 0), axis.text.y = element_text(size = 14, colour = PAL$ink), panel.spacing.x = grid::unit(16, "pt"))

  # ---- 오른쪽: 편향과 1종 오류의 관계 ----
  EP <- c("AUClast", "AUCinf_Aiii", "AUCinf_B", "AUCinf_true"); le <- unlist(L$rel$ep[EP])
  br[, epf := factor(endpoint, levels = EP, labels = le)]
  tbar <- median(br$t_q); xs <- seq(min(br$b_over_sd) - 0.1, max(br$b_over_sd) + 0.1, length.out = 200)
  cur <- data.table(x = xs, y = 100 * pnorm(xs - tbar))
  hl <- br[endpoint == "AUCinf_Aiii" & scenario == "Vmax_up_080"][, lab := unlist(ML[pk_model])]
  hl[, lx := b_over_sd + fifelse(pk_model == "k2020", 0.30, -0.30)][, ly := sim_pass_pct + fifelse(pk_model == "k2020", -8, 9)][, hj := fifelse(pk_model == "k2020", 0, 1)]
  p2 <- ggplot(br, aes(b_over_sd, sim_pass_pct)) +
    geom_hline(yintercept = nomv, linetype = "22", colour = PAL$ink2, linewidth = 0.6) + geom_vline(xintercept = 0, colour = PAL$grid, linewidth = 0.6) +
    geom_line(data = cur, aes(x, y), inherit.aes = FALSE, colour = PAL$ink2, linewidth = 0.8) +
    geom_point(aes(colour = epf, shape = epf, fill = epf), size = 2.8, stroke = 1.0) +
    geom_segment(data = hl, aes(x = lx, y = ly, xend = b_over_sd, yend = sim_pass_pct), inherit.aes = FALSE, colour = PAL$ink2, linewidth = 0.4) +
    geom_label(data = hl, aes(x = lx, y = ly, label = lab, hjust = hj), inherit.aes = FALSE, size = PT(14), family = FONT, colour = PAL$ink, fill = "white", label.size = 0,
               label.padding = grid::unit(0.1, "lines")) +
    annotate("text", x = min(xs) + 0.02, y = max(br$sim_pass_pct) * 0.98, label = L$rel$formula, hjust = 0, vjust = 1, size = PT(16), family = FONT, colour = PAL$ink) +
    scale_colour_manual(values = setNames(c(PAL$blue, PAL$orange, PAL$orange, PAL$muted), le), name = NULL) +
    scale_fill_manual(values = setNames(c(PAL$blue, PAL$orange, "white", PAL$muted), le), name = NULL) +
    scale_shape_manual(values = setNames(c(22, 23, 23, 21), le), name = NULL) +
    scale_y_continuous(labels = function(v) paste0(v, "%"), expand = expansion(mult = c(0.02, 0.04))) +
    labs(x = L$rel$xlab, y = NULL, subtitle = L$rel$sub) + theme_core(16) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank())
  dsrc("boundary dot plot and bias relation", c(PF, BR), "(figure)")

  # ---- 본문·캡션 ----
  body <- c(tx("A5d.body1", f),
            tx("A5d.body2", list(nr = dcount(BR, "TRUE", "boundary cell x single-endpoint rows"),
                                 dmax = dderived("largest |predicted - simulated| single-endpoint boundary pass, rounded up (percentage points)", BR, "all rows :: max(abs(pred_pass_pct - sim_pass_pct)), rounded up", dmax, paste0(fnum(cl(dmax, 1), 1), "%p")),
                                 dcur = dderived("largest |drawn curve - simulated| single-endpoint boundary pass, rounded up (percentage points)", BR, "all rows :: max(abs(100 pnorm(b_over_sd - median(t_q)) - sim_pass_pct)), rounded up", dcur, paste0(fnum(cl(dcur, 1), 1), "%p")))))
  bt <- unlist(.read("config/oc_design.yaml")$boundary_targets)
  cap <- tx("A5d.caption", list(tgt = dderived("boundary true AUCinf ratios (both)", "config/oc_design.yaml", "boundary_targets :: both values", bt, paste(fnum(bt, 2), collapse = ", ")),
                               reps = f_reps("boundary"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), ci = f_ci_level()))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  fw <- GEO$CW * 0.56; gap <- 0.25
  deck_figure(p1, "a5d_type1", c(GEO$ML, y0, fw, by - 0.08 - y0), src = c(PF))
  deck_figure(p2, "a5d_relation", c(GEO$ML + fw + gap, y0, GEO$CW - fw - gap, by - 0.08 - y0), src = c(BR))
  deck_visual(c(GEO$ML, y0, GEO$CW, by - 0.08 - y0))

  # ---- 노트 ----
  rel <- function(m, sc, ep, col, d_, item, unit = "%") dv(BR, sprintf("pk_model=='%s' & scenario=='%s' & endpoint=='%s'", m, sc, ep), col, d_, unit, item)
  mxr <- function(cf, pk) dv(PF, sprintf("%s & pk_model=='%s' & pass_pct==max(pass_pct[analysis_model=='M1' & kind=='boundary' & config=='%s' & pk_model=='%s'])", W(cf), pk, cf, pk), "pass_pct", 2, "%", sprintf("largest cell %s %s M1", cf, pk))
  deck_notes(tx("A5d.notes", c(f, list(
    i20 = mxr("G2A_iii", "k2020"), i16 = mxr("G2A_iii", "k2016"),
    ci = dci(PF, paste(W("P2"), "& pk_model=='k2020' & code=='V2_up_080'"), "pass_pct", "lo", "hi", 2, "%", "AUClast + Cmax, 2020 model V2 up cell, M1"),
    ke = dci(PF, paste(W("G2A_iii"), "& pk_model=='k2016' & code=='ke_up_080'"), "pass_pct", "lo", "hi", 2, "%", "AUCinf (set iii) + Cmax, 2016 ke up cell, M1"),
    ext = dint(PF, paste(W("P2"), "& pk_model=='k2020' & code=='V2_up_080'"), "n_trials", "trials in the extended cell"),
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"),
    b20 = rel("k2020", "Vmax_up_080", "AUCinf_Aiii", "bias_pct", 1, "AUCinf (rule A) bias toward 1, 2020 Vmax up (%)"),
    b16 = rel("k2016", "Vmax_up_080", "AUCinf_Aiii", "bias_pct", 1, "AUCinf (rule A) bias toward 1, 2016 Vmax up (%)"),
    s20 = rel("k2020", "Vmax_up_080", "AUCinf_Aiii", "b_over_sd", 2, "b/SD, AUCinf rule A, 2020 Vmax up", ""),
    s16 = rel("k2016", "Vmax_up_080", "AUCinf_Aiii", "b_over_sd", 2, "b/SD, AUCinf rule A, 2016 Vmax up", ""),
    q20 = rel("k2020", "Vmax_up_080", "AUCinf_Aiii", "pred_pass_pct", 1, "predicted single-endpoint pass, AUCinf rule A, 2020 Vmax up"),
    m20 = rel("k2020", "Vmax_up_080", "AUCinf_Aiii", "sim_pass_pct", 1, "simulated single-endpoint pass, AUCinf rule A, 2020 Vmax up"),
    q16 = rel("k2016", "Vmax_up_080", "AUCinf_Aiii", "pred_pass_pct", 1, "predicted single-endpoint pass, AUCinf rule A, 2016 Vmax up"),
    m16 = rel("k2016", "Vmax_up_080", "AUCinf_Aiii", "sim_pass_pct", 1, "simulated single-endpoint pass, AUCinf rule A, 2016 Vmax up"),
    bb20 = rel("k2020", "Vmax_up_080", "AUCinf_B", "bias_pct", 1, "AUCinf (rule B) bias toward 1, 2020 Vmax up (%)"),
    bb16 = rel("k2016", "Vmax_up_080", "AUCinf_B", "bias_pct", 1, "AUCinf (rule B) bias toward 1, 2016 Vmax up (%)"),
    bl = rel("k2020", "V2_up_080", "AUClast", "bias_pct", 1, "AUClast bias toward 1 relative to the true AUCinf ratio, 2020 V2 up (%)"),
    tq = drange(BR, "TRUE", "t_q", 3, "", "t quantile (0.95, n - 3), all rows"),
    tr = drange(BR, "endpoint=='AUCinf_true'", "sim_pass_pct", 1, "%", "simulated boundary pass with the true AUCinf, all 16 cells"),
    naw = dcount(PD, "analysis_model=='M1' & auclast_bias_dir=='away_from_1'", "boundary cells where the NCA AUClast GMR lies beyond the true AUCinf ratio, M1"),
    i0 = dcount("criteria/criteria_g2_type1.csv", "analysis_model=='M0' & config=='G2_A_iii' & pass_pct > 5", "AUCinf (set iii) + Cmax cells above 5%, M0"),
    p0 = dcount(T1, "analysis_model=='M0' & config=='P2' & pass_pct > 5", "AUClast + Cmax cells above 5%, M0")))))
  deck_end()
}
