# A5c 운용특성 곡선 전체(v1.2, 지시 2026-09-29 "S9 재설계" §5 "새 별첨: 모든 OC 곡선"; config/prereg_20260929_oc.yaml section8, D-064):
# 행 = PK 모델(위 2020 주 모델), 열 = 바꾼 기전 5개(F, ke, Vmax, V2, ka), 구성 5개(P2, F3A_iii, F3B, G2A_iii, G2B). 분석 모형 M1.
# 표식: 구성 계열은 색(AUClast + Cmax 파랑, 세 지표 검정, AUCinf + Cmax 주황), 규칙은 선과 채움(규칙 A 실선·채운 점, 규칙 B 파선·빈 점).
# 회색 띠 = 참 Cmax 비가 동등 한계 밖인 칸(모든 구성이 Cmax에서 실패). 동일 제품 칸(S00)은 모든 기전의 참 비 1.00 점이다.
slide_A5c <- function() {
  PF <- "oc_curves/oc_curves_pass.csv"; T1 <- "oc_curves/oc_type1_summary.csv"; T2 <- "oc_curves/oc_type2_summary.csv"; UR <- "oc_curves/oc_unreachable_cells.csv"
  deck_slide("A5c", tag = "sim")
  L <- DK$txt$A5c; CF <- c("P2", "F3A_iii", "F3B", "G2A_iii", "G2B"); MK <- c("F", "ke", "Vmax", "V2", "ka"); MS <- names(DK$txt$common$models_short)
  d0 <- rows(PF, "analysis_model=='M1' & config %in% c('P2','F3A_iii','F3B','G2A_iii','G2B')")
  s0 <- d0[code == "S00"]; premise(nrow(s0) == 2 * length(CF), "identity cell for both models and five configurations")
  d <- rbind(d0[code != "S00"], rbindlist(lapply(MK, function(k) copy(s0)[, mechanism := k])))
  premise(setequal(unique(d$mechanism), MK), "five mechanisms")
  ur <- rows(UR, "TRUE"); nt <- length(unlist(.read("config/prereg_20260929_oc.yaml")$section8$cells$targets))
  cnt <- d[config == "P2", .N, by = .(pk_model, mechanism)][ur[, .(nu = .N), by = .(pk_model, mechanism)], on = c("pk_model", "mechanism"), nu := i.nu]
  premise(all(cnt$N + fcoalesce(cnt$nu, 0L) == nt), "every mechanism x model has all nine targets, reachable or listed as unreachable")

  # ---- 제목 전제 ----
  w <- dcast(d0, pk_model + code ~ config, value.var = "pass_pct")
  premise(all(w$F3A_iii <= pmin(w$P2, w$F3B, w$G2A_iii, w$G2B) + 1e-9), "three endpoints (rule A) pass least often in every cell (title)")
  premise(all(w$F3A_iii <= w$P2 + 1e-9 & w$F3B <= w$P2 + 1e-9), "adding AUCinf never passes more often than AUClast + Cmax (body line 1)")
  for (cf in c("G2A_iii", "G2B")) for (m in MS) {
    b <- d0[config == cf & pk_model == m & kind == "boundary"]; premise(b[which.max(pass_pct), mechanism] == "Vmax", sprintf("%s %s: largest boundary pass is a target-mediated elimination cell (title)", cf, m)) }
  y0 <- core_title(tx("A5c.title"), tx("A5c.kicker"))

  # ---- 그림 ----
  lab <- unlist(L$fig$cfg[CF])
  d[, cfg := factor(config, levels = CF, labels = lab)][, mech := factor(mechanism, levels = MK, labels = unlist(L$fig$mech[MK]))]
  d[, mod := factor(unlist(DK$txt$common$models_short[pk_model]), levels = unlist(DK$txt$common$models_short))]
  band <- unique(d[cmax_in_limits == FALSE & code != "S00", .(mod, mech, target)])
  lims <- c(0.80, 1.25); xl <- c(0.775, 1.29)
  cols <- setNames(c(PAL$blue, PAL$ink, PAL$ink, PAL$orange, PAL$orange), lab)
  fills <- setNames(c(PAL$blue, PAL$ink, "white", PAL$orange, "white"), lab)
  shp <- setNames(c(21, 22, 22, 23, 23), lab); lty <- setNames(c("solid", "solid", "22", "solid", "22"), lab)
  p <- ggplot(d, aes(target, pass_pct, colour = cfg, fill = cfg, shape = cfg, linetype = cfg, group = cfg)) +
    geom_rect(data = band, aes(xmin = target / 1.03, xmax = target * 1.03, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = PAL$grid, colour = NA) +
    geom_vline(xintercept = lims, linetype = "12", colour = PAL$ink2, linewidth = 0.5) +
    geom_line(linewidth = 0.75) + geom_point(size = 2.4, stroke = 0.9) +
    facet_grid(mod ~ mech, switch = "y") +
    scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = fills, name = NULL) + scale_shape_manual(values = shp, name = NULL) +
    scale_linetype_manual(values = lty, name = NULL) +
    scale_x_log10(limits = xl, breaks = c(0.80, 1.00, 1.25), labels = function(v) formatC(v, format = "f", digits = 2), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 50), labels = function(v) paste0(v, "%"), expand = expansion(mult = c(0.02, 0.04))) +
    labs(x = L$fig$xlab, y = NULL) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.key.width = grid::unit(30, "pt"),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.spacing.x = grid::unit(14, "pt"), panel.spacing.y = grid::unit(20, "pt"),
          strip.text = element_text(size = 15, colour = PAL$ink, face = "bold", hjust = 0), strip.placement = "outside",
          strip.text.y.left = element_text(size = 15, colour = PAL$ink, face = "bold", angle = 90, hjust = 0.5))
  dsrc("operating characteristic curves, all mechanisms, both models", c(PF, UR), "(figure)")

  # ---- 본문 ----
  t1 <- function(cf, sc = "both") dv(T1, sprintf("config=='%s' & scope=='%s'", cf, sc), "max_pct", 1, "%", sprintf("largest boundary type I error, %s, %s, M1", cf, sc))
  t2id <- function(cf) { r <- rows(PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s'", cf)); x <- range(100 - r$pass_pct)
    dderived(sprintf("type II error, identical product, %s, M1, range over two models", cf), PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s' :: range(100 - pass_pct)", cf), x, rng_fmt(x[1], x[2], 1, "%")) }
  g2b <- rows(PF, "code=='S00' & analysis_model=='M1' & config=='G2B'"); p2 <- rows(PF, "code=='S00' & analysis_model=='M1' & config=='P2'")
  premise(all(100 - g2b[order(pk_model)]$pass_pct <= 100 - p2[order(pk_model)]$pass_pct), "rule B AUCinf + Cmax: identical-product type II not above AUClast + Cmax in either model (body line 2)")
  for (m in MS) premise(row1(T1, sprintf("config=='G2B' & scope=='%s'", m))$max_pct > row1(T1, sprintf("config=='P2' & scope=='%s'", m))$max_pct,
                        sprintf("rule B AUCinf + Cmax: larger boundary type I maximum than AUClast + Cmax, %s (body line 2)", m))
  ka <- d0[mechanism == "ka" & config == "P2" & code != "S00"]; premise(any(!ka$cmax_in_limits) & all(d0[mechanism != "ka" & code != "S00" & !cmax_in_limits, mechanism] %in% "V2"),
                                                                    "Cmax outside the limits: ka cells and the V2 cells only (body line 3)")
  body <- tx("A5c.body")
  rn <- dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "near_one", "trials"), "trials per cell, true ratio 0.95 and 1.05", function(x) fnum(as.numeric(x), 0, TRUE))
  ro <- dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "other", "trials"), "trials per cell, true ratio 0.85, 0.90, 1.11 and 1.18", function(x) fnum(as.numeric(x), 0, TRUE))
  nr <- as.numeric(unlist(.read("config/prereg_20260929_oc.yaml")$section8$reps$near_one$targets))
  cap <- tx("A5c.caption", list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), rx = dint(PF, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & code=='V2_up_080'", "n_trials", "trials in the extended cell"), ro = ro))
  ab <- list(a = dderived("near-one targets, lower (true AUCinf ratio)", "config/prereg_20260929_oc.yaml", "section8.reps.near_one.targets :: first", nr[1], fnum(nr[1], 2)),
             b = dderived("near-one targets, upper (true AUCinf ratio)", "config/prereg_20260929_oc.yaml", "section8.reps.near_one.targets :: second", nr[2], fnum(nr[2], 2)))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "a5c_oc_all", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(PF, UR))

  # ---- 노트 ----
  urs <- ur[order(match(pk_model, MS), mechanism, target)]; premise(setequal(unique(urs$mechanism), c("V2", "ka")), "unreachable cells are V2 and ka cells only (notes)")
  deck_notes(tx("A5c.notes", list(
    g2 = t1("G2A_iii"), g2_20 = t1("G2A_iii", "k2020"), g2_16 = t1("G2A_iii", "k2016"), g2b_20 = t1("G2B", "k2020"), g2b_16 = t1("G2B", "k2016"),
    f3b = t1("F3B"), g2id = t2id("G2A_iii"), f3bid = t2id("F3B"), f3 = t1("F3A_iii"), p2 = t1("P2"), f3id = t2id("F3A_iii"), p2id = t2id("P2"), g2bid = t2id("G2B"), g2b = t1("G2B"),
    nur = dcount(UR, "TRUE", "unreachable mechanism x target cells, both models"),
    v2u = { x <- sort(urs[mechanism == "V2" & pk_model == "k2020", target]); premise(identical(x, sort(urs[mechanism == "V2" & pk_model == "k2016", target])), "same unreachable V2 targets in both models (notes)")
      dderived("unreachable V2 targets, both models", UR, "mechanism=='V2' :: target", x, paste(fnum(x, 2), collapse = ", ")) },
    ka20 = { x <- sort(urs[mechanism == "ka" & pk_model == "k2020", target]); dderived("unreachable ka targets, 2020 model", UR, "mechanism=='ka' & pk_model=='k2020' :: target", x, paste(fnum(x, 2), collapse = ", ")) },
    ka16 = { x <- sort(urs[mechanism == "ka" & pk_model == "k2016", target]); dderived("unreachable ka targets, 2016 model", UR, "mechanism=='ka' & pk_model=='k2016' :: target", x, paste(fnum(x, 2), collapse = ", ")) },
    kac = drange(PF, "analysis_model=='M1' & config=='P2' & mechanism=='ka'", "cmax_ratio", 2, "", "true Cmax ratio, ka cells, both models"),
    v2c = dv(PF, "analysis_model=='M1' & config=='P2' & pk_model=='k2020' & code=='V2_down_118'", "cmax_ratio", 2, "", "true Cmax ratio, 2020 model V2 cell at 1.18"),
    c4a = dderived("rule 4c range, lower", "config/prereg_20260929_oc.yaml", "section8.wording.rule_4c :: lower end", 0.90, "0.90"),
    c4b = dderived("rule 4c range, upper", "config/prereg_20260929_oc.yaml", "section8.wording.rule_4c :: upper end", 0.95, "0.95"),
    n4c = dcount("oc_curves/oc_rule4c_cells.csv", "TRUE", "rule 4c cells (G2 passes more often than P2 at true ratio 0.90-0.95, paired interval above 0)"),
    c4det = core_c4det(L$fig$c4cell, L$fig$mech),
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), ci = f_ci_level(), lim = { lim_v <- as.numeric(unlist(.read("config/trial_design.yaml")$be$limits))
      dderived("equivalence limits in percent (definition text)", "config/trial_design.yaml", "be.limits :: x 100, printed as a range", 100 * lim_v, rng_fmt(100 * lim_v[1], 100 * lim_v[2], 2)) }, rb = f_reps("boundary"), rn = rn, ro = ro, a = ab$a, b = ab$b)))
  deck_end()
}
