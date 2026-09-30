# S7 ③ 판정 성능(v1.2, 지시 2026-09-29 "S9 재설계: 1종·2종 오류", config/prereg_20260929_oc.yaml section8, D-064):
# 왼쪽 65% = 주 모델(2020)·표적 매개 소실(Vmax) 기전의 운용특성 곡선(참 AUCinf 비 0.80~1.25, 로그 눈금; P2, F3A_iii, G2A_iii, 이상적 계단),
# 오른쪽 35% = 요약표(구성 3 × 1종 오류 최대(경계 16칸), 2종 오류(동일 제품), 2종 오류(참 비 0.95~1.05); 각 칸은 두 모델 범위).
# 제목·본문 문구는 사전 등록 규칙(oc_wording_rules.csv: T, 4a, 4b, protect, G2_type1)으로 고른다. 분석 모형 M1.
slide_S7 <- function() {
  PF <- "oc_curves/oc_curves_pass.csv"; T1 <- "oc_curves/oc_type1_summary.csv"; T2 <- "oc_curves/oc_type2_summary.csv"; WR <- "oc_curves/oc_wording_rules.csv"; PD <- "oc_curves/oc_paired.csv"
  deck_slide("S7", tag = "sim")
  W <- rows(WR, "TRUE"); premise(nrow(W) == 2, "wording rules for both models")
  allr <- function(k) all(W[[k]])
  L <- DK$txt$S7
  # ---- 제목(규칙) ----
  g2c <- if (allr("rule_4a_keep_type2") && allr("rule_G2_type1")) "g2_both" else if (allr("rule_G2_type1")) "g2_type1" else "g2_none"
  f3c <- if (allr("rule_4b_small_cost")) L$title$f3_small else if (g2c %in% c("g2_both", "g2_type1")) L$title$f3_cost_ellipsis else L$title$f3_cost   # "A하면 X가, B하면 Y가 는다"(2줄 제한)
  tt <- c(if (allr("rule_T")) L$title$prefix, paste0(f3c, ", ", L$title[[g2c]]))
  y0 <- core_title(paste(tt, collapse = "\n"), tx("S7.kicker"))

  # ---- 왼쪽: 운용특성 곡선(2020 모델, Vmax) ----
  M <- "k2020"; CF <- c("P2", "F3A_iii", "G2A_iii")
  d <- rows(PF, sprintf("pk_model=='%s' & analysis_model=='M1' & (mechanism=='Vmax' | code=='S00') & config %%in%% c('P2','F3A_iii','G2A_iii')", M))
  premise(nrow(d) == 9 * 3, "nine targets (0.80-1.25) x three configurations for the 2020 model, Vmax")
  d[, cfg := factor(config, levels = CF, labels = unlist(L$fig$cfg[CF]))]
  lims <- c(0.80, 1.25); xl <- c(0.70, 1.43)
  shade <- data.table(xmin = c(xl[1], lims[2]), xmax = c(lims[1], xl[2]))
  ideal <- data.table(x = c(xl[1], lims[1], lims[1], lims[2], lims[2], xl[2]), y = c(0, 0, 100, 100, 0, 0))
  cols <- setNames(c(PAL$blue, PAL$ink, PAL$orange), unlist(L$fig$cfg[CF])); shp <- setNames(c(16, 15, 18), unlist(L$fig$cfg[CF])); lty <- setNames(c("solid", "22", "solid"), unlist(L$fig$cfg[CF]))
  p <- ggplot(d, aes(target, pass_pct, colour = cfg, shape = cfg, linetype = cfg)) +
    annotate("rect", xmin = shade$xmin, xmax = shade$xmax, ymin = -Inf, ymax = Inf, fill = PAL$tint_orange) +
    annotate("rect", xmin = lims[1], xmax = lims[2], ymin = -Inf, ymax = Inf, fill = PAL$tint_grey) +
    geom_path(data = ideal, aes(x, y), inherit.aes = FALSE, linetype = "12", colour = PAL$ink2, linewidth = 0.6) +
    geom_line(linewidth = 1.0) + geom_point(size = 3.4) +
    annotate("label", x = sqrt(xl[1] * lims[1]), y = 55, label = L$fig$out, size = PT(14), family = FONT, colour = PAL$orange, fill = PAL$tint_orange, label.size = 0, lineheight = 0.95) +
    annotate("label", x = sqrt(lims[2] * xl[2]), y = 55, label = L$fig$out, size = PT(14), family = FONT, colour = PAL$orange, fill = PAL$tint_orange, label.size = 0, lineheight = 0.95) +
    annotate("label", x = 1, y = 40, label = L$fig$inn, size = PT(14), family = FONT, colour = PAL$ink2, fill = PAL$tint_grey, label.size = 0, lineheight = 0.95) +
    scale_colour_manual(values = cols, name = NULL) + scale_shape_manual(values = shp, name = NULL) + scale_linetype_manual(values = lty, name = NULL) +
    scale_x_log10(limits = xl, breaks = c(0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.11, 1.18, 1.25), labels = function(v) formatC(v, format = "f", digits = 2), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = L$fig$xlab, y = NULL, subtitle = fill(L$fig$sub, list(model = DK$txt$common$models[[M]]))) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  dsrc("operating characteristic curve, 2020 model, Vmax", c(PF), "(figure)")

  # ---- 오른쪽: 요약표(각 칸 = 두 모델 범위) ----
  t1v <- function(cf) dv(T1, sprintf("config=='%s' & scope=='both'", cf), "max_pct", 1, "%", sprintf("largest boundary type I error over 16 cells, %s, M1", cf))
  t2id <- function(cf) { r <- rows(PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s'", cf)); x <- range(100 - r$pass_pct)
    dderived(sprintf("type II error, identical product, %s, M1, range over two models", cf), PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s' :: range(100 - pass_pct)", cf), x, rng_fmt(x[1], x[2], 1, "%")) }
  t2nr <- function(cf) { r <- rows(T2, sprintf("config=='%s' & t2group=='near_095_105' & scope=='both'", cf)); premise(nrow(r) == 1, "one near-one row")
    dderived(sprintf("type II error, true ratio 0.95 and 1.05, %s, M1, range over mechanisms and models", cf), T2, sprintf("config=='%s' & t2group=='near_095_105' & scope=='both' :: min_pct, max_pct", cf),
             c(r$min_pct, r$max_pct), rng_fmt(r$min_pct, r$max_pct, 0, "%")) }
  nt <- as.numeric(unlist(.read("config/prereg_20260929_oc.yaml")$section8$reps$near_one$targets))
  ab <- list(a = dderived("near-one targets, lower (true AUCinf ratio)", "config/prereg_20260929_oc.yaml", "section8.reps.near_one.targets :: first", nt[1], fnum(nt[1], 2)),
             b = dderived("near-one targets, upper (true AUCinf ratio)", "config/prereg_20260929_oc.yaml", "section8.reps.near_one.targets :: second", nt[2], fnum(nt[2], 2)))
  n16 <- dcount(PF, "analysis_model=='M1' & config=='P2' & kind=='boundary'", "boundary cells, both models")
  mxc <- rows(T2, "t2group=='near_095_105' & scope=='both' & config %in% c('P2','F3A_iii','G2A_iii')")$max_cell
  premise(length(unique(mxc)) == 1 && grepl(" ka_up_", mxc[1]), "near-one type II maximum is the same absorption-rate cell for all three configurations (table note)")
  kac <- rows(PF, sprintf("config=='P2' & analysis_model=='M1' & code=='%s'", sub("^\\S+ ", "", mxc[1])))
  premise(all(kac$cmax_in_limits), "that cell's true Cmax ratio is inside the limits (it stays in the type II cells by the pre-registered rule)")
  tnote <- fill(L$tab$note, list(cr = drange(PF, sprintf("config=='P2' & analysis_model=='M1' & code=='%s'", kac$code[1]), "cmax_ratio", 2, "", "true Cmax ratio of the near-one maximum cell, two models"),
                                 t = dderived("target of the near-one maximum cell", PF, sprintf("config=='P2' & analysis_model=='M1' & code=='%s' :: target", kac$code[1]), kac$target[1], fnum(kac$target[1], 2))))
  tab <- data.table(a = unlist(L$tab$rows[CF]), b = vapply(CF, t1v, ""), c = vapply(CF, t2id, ""), d = paste0(vapply(CF, t2nr, ""), L$tab$dagger))
  setnames(tab, vapply(unlist(L$tab$head), function(h) fill(h, c(ab, list(n = n16))), ""))
  # ---- 본문(규칙) ----
  red <- { r <- rows(PD, "comparison=='F3A_iii - P2' & kind=='boundary'"); x <- max(-r$diff_pass_pp)
    dderived("largest type I error reduction by adding AUCinf (F3A_iii vs P2), 16 boundary cells, M1 (percentage points)", PD, "comparison=='F3A_iii - P2' & kind=='boundary' :: max(-diff_pass_pp)", x, paste0(fnum(x, 1), "%p")) }
  inc_id <- { r <- rows(PD, "comparison=='F3A_iii - P2' & scenario=='S00'"); x <- range(-r$diff_pass_pp)
    dderived("type II error increase by adding AUCinf (F3A_iii vs P2), identical product, two models (percentage points)", PD, "comparison=='F3A_iii - P2' & scenario=='S00' :: range(-diff_pass_pp)", x, rng_fmt(x[1], x[2], 1, "%p")) }
  inc_max <- { x <- max(W$F3_t2_inc_max); dderived("largest type II error increase by adding AUCinf over the type II cells (F3A_iii vs P2), both models (percentage points)", WR, "TRUE :: max(F3_t2_inc_max)", x, paste0(fnum(x, 0), "%p")) }
  # 1종 오류가 줄어드는 칸(사용자 지시 2026-09-29: 말초 분포 극단 칸 외에도 여러 칸이면 본문에 밝힌다): 감소 >= 보호 규칙 기준(1.0%p) 칸 수와 그중 P2가 이미 명목 이하인 칸 수
  thp_v <- 1.0; thp <- dderived("wording rule protect threshold (percentage points)", "config/prereg_20260929_oc.yaml", "section8.wording.rule_protect :: threshold in the text", thp_v, "1.0%p")
  rb1 <- merge(rows(PD, "comparison=='F3A_iii - P2' & kind=='boundary'")[, .(pk_model, scenario, red = -diff_pass_pp)],
               rows(PF, "analysis_model=='M1' & kind=='boundary' & config=='P2'")[, .(pk_model, scenario = code, p2 = pass_pct)], by = c("pk_model", "scenario"))
  premise(nrow(rb1) == 16, "16 boundary cells with paired reductions")
  k1 <- dderived("boundary cells where adding AUCinf (F3A_iii) lowers type I error by at least 1.0 point, both models, M1", PD, "comparison=='F3A_iii - P2' & kind=='boundary' :: count(-diff_pass_pp >= 1.0)", sum(rb1$red >= thp_v), fnum(sum(rb1$red >= thp_v), 0))
  k2 <- dderived("of those, cells where AUClast + Cmax is already at or below 5%, M1", PF, "P2 boundary cells with F3A_iii reduction >= 1.0 point :: count(pass_pct <= 5)", sum(rb1$red >= thp_v & rb1$p2 <= 5), fnum(sum(rb1$red >= thp_v & rb1$p2 <= 5), 0))
  b1 <- if (allr("rule_4b_small_cost")) fill(L$body$f3_small, list(x = red, y = inc_id)) else if (allr("rule_protect_almost_none")) fill(L$body$f3_none, list(x = red, y = inc_id, z = inc_max)) else fill(L$body$f3_some, list(x = red, y = inc_id, z = inc_max, k1 = k1, k2 = k2, thp = thp, nom = f_nominal()))
  b2 <- if (allr("rule_4a_keep_type2") && allr("rule_G2B_type1")) L$body$g2_both else if (allr("rule_G2B_type1")) L$body$g2_type1 else L$body$g2_none
  body <- c(b1, b2)
  rn <- dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "near_one", "trials"), "trials per cell, true ratio 0.95 and 1.05", function(x) fnum(as.numeric(x), 0, TRUE))
  ro <- dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "other", "trials"), "trials per cell, true ratio 0.85, 0.90, 1.11 and 1.18", function(x) fnum(as.numeric(x), 0, TRUE))
  cap <- tx("S7.caption", list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), ro = ro))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  fw <- GEO$CW * 0.58; tw <- GEO$CW - fw - 0.25
  deck_figure(p, "s7_oc_k2020_vmax", c(GEO$ML, y0, fw, by - 0.08 - y0), src = c(PF))
  th <- deck_table_h(tab, tw, c(1.72, 0.92, 1.05, 1.20), 14) + 0.04
  ty <- y0 + 0.42
  deck_text(fill(L$tab$title, list(n = n16)), c(GEO$ML + fw + 0.25, y0, tw, 0.38), size = 14, bold = TRUE, color = PAL$ink2, label = "caption_tabtitle")
  deck_table(tab, c(GEO$ML + fw + 0.25, ty, tw, th), widths = c(1.72, 0.92, 1.05, 1.20), size = 14, highlight = 1, highlight_fill = PAL$tint_blue, label = "table_oc")
  nh <- est_height(tnote, tw, 14, 0) + 0.04
  deck_text(tnote, c(GEO$ML + fw + 0.25, ty + th + 0.08, tw, nh), size = 14, color = PAL$ink2, label = "caption_tabnote")
  deck_visual(c(GEO$ML, y0, GEO$CW, by - 0.08 - y0))

  # ---- 노트 ----
  c4 <- rows("oc_curves/oc_rule4c_cells.csv", "TRUE")
  premise(rb1[which.max(red), pk_model] == "k2020" && rb1[which.max(red), scenario] == "V2_up_080" && rb1[p2 > 5, .N] == 1 && rb1[p2 > 5, scenario] == "V2_up_080", "largest reduction is the 2020 V2 cell, the only P2 cell above 5% (notes)")
  deck_notes(tx("S7.notes", list(x = red, z = inc_max, nom = f_nominal(),
    ruleT = L$yn[[as.character(allr("rule_T"))]], r4a = L$yn[[as.character(allr("rule_4a_keep_type2"))]], r4b = L$yn[[as.character(allr("rule_4b_small_cost"))]],
    rpr = L$yn[[as.character(allr("rule_protect_almost_none"))]],
    p2max = t1v("P2"), g2bmax = t1v("G2B"), f3bmax = t1v("F3B"),
    n4c = dderived("cells at true ratio 0.90-0.95 where G2 passes more often than P2 (paired interval above 0)", "oc_curves/oc_rule4c_cells.csv", "count of rows", nrow(c4), fnum(nrow(c4), 0)),
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), rb = f_reps("boundary"), rn = rn, ro = ro, a = ab$a, b = ab$b,
    rx = dint("oc_models/type1_models.csv", "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials in the extended cell"),
    c = dderived("type II group (c), lower target", "config/prereg_20260929_oc.yaml", "section8.statistics.type2 :: (c) lower", 0.90, "0.90"),
    d = dderived("type II group (c), upper target", "config/prereg_20260929_oc.yaml", "section8.statistics.type2 :: (c) upper", 1.11, "1.11"),
    thT = dderived("wording rule T threshold: P2 type I maximum (%)", "config/prereg_20260929_oc.yaml", "section8.wording.rule_T :: threshold in the text", 6.0, "6.0%"),
    th4 = dderived("wording rules 4a and 4b threshold (percentage points)", "config/prereg_20260929_oc.yaml", "section8.wording.rule_4a and rule_4b :: threshold in the text", 2.0, "2.0%p"),
    thp = thp,
    c4a = dderived("rule 4c range, lower", "config/prereg_20260929_oc.yaml", "section8.wording.rule_4c :: lower end", 0.90, "0.90"),
    c4b = dderived("rule 4c range, upper", "config/prereg_20260929_oc.yaml", "section8.wording.rule_4c :: upper end", 0.95, "0.95"))))
  deck_end()
}
