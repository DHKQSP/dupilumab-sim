# S17 논거 ③-2 시험 수준: 시나리오별 시험 GMR의 기하평균(exp(mean_log_gmr))을 참 AUC0-inf 비와 비교. 평가변수: 참 AUC0-inf(모델 적분, 잔차 없음),
# AUC0-last, NCA AUC0-inf 규칙 A (i)·A (ii)·B·C (i). 분석 모형 M1(이 덱의 기준), M0 병기(노트). 경계 16칸 + 제품 시나리오(S00 동일 제품, F097 = F ×0.97).
# 자료: results/oc_models/sd_se_models.csv(시험 모집단 재생성 시험, 칸당 10,000회; 2020 모델 V2 증가는 20,000회),
#       results/criteria/criteria_bias.csv(같은 시험의 기준 세트 (iii)·(iv) 편향), results/rationale/pillar2_products_B0.csv와
#       results/trials/products_concordance*.csv(v1.0 제품 시나리오 시험의 AUC0-last 대 NCA 판정 일치율, 로그 GMR 상관).
# 끝점 코드: AUCinf_A = 규칙 A 세트 (ii), AUCinf_Ai = 규칙 A 세트 (i), AUCinf_Ci = 규칙 C 세트 (i)(R/oc.R, R/oc_models.R).
# 미산출: AUC0-last 판정과 이상적 AUC0-inf(모델 적분값) 판정의 시험별 일치율은 프로젝트에서 계산하지 않았다.
s17_maxabs <- function(rel, where, col, d, unit, item) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- max(abs(r[[col]]))
  dderived(item, rel, sprintf("max(abs(%s)) over %d rows [%s]", col, nrow(r), where), x, paste0(fnum(x, d), unit))
}
# 두 파일(두 모델)에 걸친 범위: 최솟값과 최댓값을 각각 그 값이 있는 파일로 기록한다
s17_range2 <- function(rels, where, col, d, unit, item) {
  v <- lapply(rels, function(r) { x <- rows(r, where)[[col]]; premise(length(x) >= 1, sprintf("%s [%s] matched no rows", r, where)); x })
  mn <- vapply(v, min, 1); mx <- vapply(v, max, 1); i <- which.min(mn); j <- which.max(mx)
  a <- dderived(paste(item, "(minimum over both models)"), rels[i], sprintf("min over rows [%s] :: %s", where, col), mn[i], paste0(fnum(mn[i], d), unit))
  b <- dderived(paste(item, "(maximum over both models)"), rels[j], sprintf("max over rows [%s] :: %s", where, col), mx[j], paste0(fnum(mx[j], d), unit))
  rng_fmt(mn[i], mx[j], d, unit)
}
s17_signed <- function(p) if (!startsWith(p, "-")) paste0("+", p) else p
slide_S17 <- function() {
  SSf <- "oc_models/sd_se_models.csv"; CB <- "criteria/criteria_bias.csv"
  EP <- c("AUCinf_true", "AUClast", "AUCinf_Ai", "AUCinf_A", "AUCinf_B", "AUCinf_Ci")
  BND <- "!scenario %in% c('S00','F097')"
  wb <- function(ep, am = "M1") sprintf("analysis_model=='%s' & endpoint=='%s' & %s", am, ep, BND)
  deck_slide("S17", tag = "sim")
  L <- DK$txt$S17
  premise(nrow(rows(SSf, wb("AUClast"))) == 16 && all(vapply(EP, function(e) nrow(rows(SSf, wb(e))) == 16, TRUE)), "16 boundary cells per endpoint under M1")
  premise(nrow(rows(CB, "analysis_model=='M1' & config %in% c('A_iii','A_iv')")) == 32, "criteria_bias.csv holds the 16 boundary cells for sets (iii) and (iv), M1")
  lim <- unlist(.read("config/trial_design.yaml")$be$limits)
  lo <- dcfg("trial_design.yaml", c("be", "limits"), "lower equivalence limit", function(x) fnum(x[1], 2))
  hi <- dcfg("trial_design.yaml", c("be", "limits"), "upper equivalence limit", function(x) fnum(x[2], 2))

  # ---- 제목: 편향 절댓값 최대(AUC0-last; 규칙 A (i)·(ii)·B). 규칙 A (iii)·(iv)의 최대는 그림·표 밖이라 요점(하위 항목)에 둔다 ----
  nb <- dcount(SSf, wb("AUClast"), "boundary cells (both models), M1")
  f <- list(a = s17_maxabs(SSf, wb("AUClast"), "bias_pct", 2, "%", "AUC0-last: largest absolute bias of the geometric mean trial GMR vs the true AUC0-inf ratio, 16 boundary cells, M1"),
            b = s17_maxabs(SSf, sprintf("analysis_model=='M1' & endpoint %%in%% c('AUCinf_Ai','AUCinf_A','AUCinf_B') & %s", BND), "bias_pct", 2, "%",
                           "NCA AUC0-inf rules A (i), A (ii), B: largest absolute bias, 16 boundary cells, M1"),
            c = s17_maxabs(CB, "analysis_model=='M1' & config %in% c('A_iii','A_iv')", "bias_pct", 2, "%",
                           "NCA AUC0-inf rule A criteria sets (iii) and (iv): largest absolute bias, 16 boundary cells, M1"))
  deck_kicker(tx("S17.kicker")); deck_title(tx("S17.title", f))

  # ---- 왼쪽: 경계 16칸의 편향(동등 쪽 = 동등 판정이 쉬워지는 방향을 양수로), 평가변수별, 두 모델 ----
  d <- copy(rows(SSf, sprintf("analysis_model=='M1' & endpoint %%in%% c(%s) & %s", paste(sprintf("'%s'", EP), collapse = ","), BND)))
  premise(all(d$true_ratio != 1), "boundary cells have a true ratio different from 1")
  premise(all(abs(d[true_ratio > 1]$true_ratio - lim[2]) < 0.01) && all(abs(d[true_ratio < 1]$true_ratio - lim[1]) < 0.01), "boundary cells sit at the lower or the upper equivalence limit")
  d[, tw := ifelse(true_ratio < 1, bias_pct, -bias_pct)]
  EL <- unlist(L$ep[EP]); d[, ep := factor(EL[endpoint], levels = rev(EL))]
  # 겹침 방지: 같은 평가변수·모델의 8칸을 값 순서로 세 줄에 번갈아 놓는다(세로 위치만 어긋남, 가로 위치가 값)
  setorder(d, endpoint, pk_model, tw); d[, k := (seq_len(.N) - 1L) %% 3L, by = .(endpoint, pk_model)]
  d[, y := as.numeric(ep) + ifelse(pk_model == "k2016", 0.2, -0.2) + (k - 1) * 0.085]
  d[, model := factor(model_lab()[pk_model], levels = model_lab())]
  xr <- range(d$tw); xl <- c(floor(xr[1]) - 0.5, ceiling(xr[2]) + 0.5)
  ex_ <- rbind(d[endpoint == "AUClast"][which.max(abs(bias_pct))], d[endpoint %in% c("AUCinf_Ai", "AUCinf_A", "AUCinf_B")][which.max(abs(bias_pct))])
  premise(nrow(ex_) == 2 && abs(abs(ex_$bias_pct[1]) - max(abs(d[endpoint == "AUClast"]$bias_pct))) < 1e-12, "figure labels: the two title extremes")
  ex_[, lab := paste0(ifelse(tw > 0, "+", ""), fnum(tw, 2), "%")]
  ex_[, `:=`(hj = ifelse(tw < 0, 1, 0.5), xt = ifelse(tw < 0, tw - 0.25, tw), yt = ifelse(tw < 0, y, y + 0.32))]
  FW <- 6.05; FH <- 2.72
  p <- ggplot(d, aes(x = tw, y = y, colour = model, shape = model)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = PAL$tint_orange, alpha = 0.7) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_point(size = 2.1, alpha = 0.75, stroke = 0.6) +
    # 제목의 두 극단(AUC0-last 절댓값 최대, NCA 규칙 A (i)·(ii)·B 절댓값 최대)을 그림에 수치로 표시(자료에서 바로)
    geom_label(data = ex_, aes(x = xt, y = yt, label = lab), inherit.aes = FALSE, hjust = ex_$hj, size = 3.7, family = FONT, colour = PAL$ink,
               fill = "white", label.size = 0, label.padding = grid::unit(1, "pt")) +
    annotate("text", x = xl[2], y = length(EL) + 0.66, label = L$fig$right, hjust = 1, vjust = 1, size = 3.9, family = FONT, colour = PAL$orange, fontface = "bold") +
    annotate("text", x = xl[1], y = length(EL) + 0.66, label = L$fig$left, hjust = 0, vjust = 1, size = 3.9, family = FONT, colour = PAL$ink2) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = c(16, 17)) +
    scale_y_continuous(breaks = seq_along(levels(d$ep)), labels = levels(d$ep), limits = c(0.55, length(EL) + 0.7), expand = expansion(0)) +
    scale_x_continuous(limits = xl, breaks = seq(ceiling(xl[1] / 2) * 2, floor(xl[2] / 2) * 2, by = 2)) +
    labs(x = fill(L$fig$xlab, list(hi = hi)), y = NULL) + theme_deck(12) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          axis.title.x = element_text(hjust = 1, margin = margin(4, 0, 0, 0)),
          panel.grid.major.y = element_line(colour = PAL$grid, linewidth = 0.35), axis.text.y = element_text(colour = PAL$ink, size = 12))
  deck_figure(p, "s17_bias_strip", c(GEO$ML, GEO$BODY_TOP, FW, FH), src = SSf)

  # ---- 왼쪽 아래: 요점 ----
  v2 <- dv(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2020' & scenario=='V2_up_080'", "bias_pct", 2, "%", "AUC0-last bias, 2020 model V2 up (true ratio 0.80), M1")
  aw <- dcount(SSf, sprintf("%s & bias_dir=='away_from_1'", wb("AUClast")), "AUC0-last cells biased away from 1, M1")
  premise(as.integer(aw) == as.integer(nb) - 1L && nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1' & pk_model=='k2020' & scenario=='V2_up_080'", wb("AUClast")))) == 1,
          "AUC0-last is biased away from 1 in all boundary cells but the 2020 model V2 cell (M1)")
  premise(all(vapply(c("AUCinf_Ai", "AUCinf_A", "AUCinf_B"), function(e) nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1'", wb(e)))) > as.integer(nb) / 2, TRUE)),
          "NCA rules A (i), A (ii) and B are biased toward 1 in most boundary cells (M1)")
  tw_n <- function(ep, am = "M1") dcount(SSf, sprintf("%s & bias_dir=='toward_1'", wb(ep, am)), sprintf("cells biased toward 1 (Monte Carlo interval excludes 0), %s, %s", ep, am))
  # 기준 세트 (iii)·(iv)의 규칙 A: 같은 시험의 동등 쪽 편향 칸 수(M1, criteria_bias.csv)
  tw_cb <- function(cf) dcount(CB, sprintf("analysis_model=='M1' & config=='%s' & bias_dir=='toward_1'", cf), sprintf("cells biased toward 1, rule A criteria set config %s, M1", cf))
  premise(all(vapply(c("A_i", "A_ii", "B"), function(cf) nrow(rows(CB, sprintf("analysis_model=='M1' & config=='%s' & bias_dir=='toward_1'", cf))) ==
                       nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1'", wb(c(A_i = "AUCinf_Ai", A_ii = "AUCinf_A", B = "AUCinf_B")[[cf]])))), TRUE)),
          "criteria_bias.csv and sd_se_models.csv agree on the toward-1 counts of rules A (i), A (ii), B (same trials, M1)")
  kci <- tw_n("AUCinf_Ci")
  premise(as.integer(kci) == 1L && nrow(rows(SSf, sprintf("%s & bias_dir=='toward_1'", wb("AUClast")))) == 1L, "rule C (i) and AUC0-last each have one cell biased toward 1 (M1)")
  BY <- GEO$BODY_TOP + FH + 0.08
  deck_bullets(tx("S17.bullets", list(aw = aw, n = nb, v2 = s17_signed(v2), k_ai = tw_n("AUCinf_Ai"), k_aii = tw_n("AUCinf_A"), k_b = tw_n("AUCinf_B"), ci = kci,
                                       k_aiii = tw_cb("A_iii"), k_aiv = tw_cb("A_iv"), c = f$c)), box = c(GEO$ML, BY, FW, GEO$BODY_BOTTOM - BY), size = 16)

  # ---- 오른쪽: 표 제목, 시나리오별 기하평균비 표(M1) ----
  XR <- GEO$ML + FW + 0.3; WR <- GEO$W - GEO$MR - XR
  gm <- function(ep, m, sc) {
    r <- row1(SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & pk_model=='%s' & scenario=='%s'", ep, m, sc)); x <- exp(r$mean_log_gmr)
    dderived(sprintf("geometric mean of trial GMRs, %s, %s, %s, M1", ep, m, sc), SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & pk_model=='%s' & scenario=='%s' :: exp(mean_log_gmr)", ep, m, sc), x, fnum(x, 4))
  }
  tr <- function(m, sc) dv(SSf, sprintf("analysis_model=='M1' & endpoint=='AUClast' & pk_model=='%s' & scenario=='%s'", m, sc), "true_ratio", 4, "", sprintf("true AUC0-inf ratio (common virtual subjects), %s, %s", m, sc))
  SC <- list(c("k2016", "Vmax_up_080"), c("k2020", "ka_down_080"), c("k2016", "F097"))
  premise(all(vapply(SC, function(s) all(rows(SSf, sprintf("analysis_model=='M1' & pk_model=='%s' & scenario=='%s'", s[1], s[2]))$n_trials == 10000), TRUE)), "the three table scenarios have 10,000 trials each (caption)")
  ntr <- dint(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2016' & scenario=='F_down_080'", "n_trials", "trials per boundary cell")
  npop <- dcfg("oc_design.yaml", c("estimand", "population", "n_subjects"), "common virtual subjects for the true ratio", function(x) fnum(as.numeric(x), 0, big = TRUE))
  f97 <- dcfg("prereg_20260926.yaml", c("section1", "scenarios", "products", "F097", "F"), "F multiplier of the F097 product scenario, test arm", function(x) fnum(as.numeric(x), 2))
  CH <- 0.69
  deck_text(tx("S17.caption", list(ntr = ntr, npop = npop)), c(XR, GEO$BODY_TOP, WR, CH), size = 16, label = "text_caption", color = PAL$ink2, gap_pt = 0)
  df <- data.frame(a = unlist(L$tab_ep[EP]), b = vapply(EP, tw_n, ""),
                   c = vapply(EP, gm, "", m = SC[[1]][1], sc = SC[[1]][2]), d = vapply(EP, gm, "", m = SC[[2]][1], sc = SC[[2]][2]), e = vapply(EP, gm, "", m = SC[[3]][1], sc = SC[[3]][2]),
                   stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("S17.table.head", list(n = nb, t1 = tr(SC[[1]][1], SC[[1]][2]), t2 = tr(SC[[2]][1], SC[[2]][2]), t3 = tr(SC[[3]][1], SC[[3]][2]), f97 = f97))
  TY <- GEO$BODY_TOP + CH + 0.03; TH <- 2.59
  deck_table(df, box = c(XR, TY, WR, TH), widths = c(1.5, 0.93, 1.15, 1.15, 1.15), size = 12, highlight = 2, label = "table_gm")

  # ---- 오른쪽 아래: 미산출(판정 일치율) ----
  CONC <- c("trials/products_concordance.csv", "trials/products_concordance_struct2020.csv")
  cor <- s17_range2(CONC, "schedule=='B0'", "cor_logGMR_last_true", 2, "", "correlation of trial log GMR, AUC0-last vs true AUC0-inf (trial subjects), v1.0 product scenarios")
  agr_b <- s17_range2(CONC, "schedule=='B0'", "agree_last_infall", 1, "%", "decision agreement AUC0-last vs NCA AUC0-inf rule B (lambda-z estimable), v1.0 product scenarios")
  n_conc <- local({ x <- vapply(CONC, function(r) { v <- unique(rows(r, "schedule=='B0'")$n_trials); premise(length(v) == 1, paste("one trial count per product scenario in", r)); v }, 1)
    premise(length(unique(x)) == 1, "same trial count in both concordance files")
    for (r in CONC) dderived("trials per product scenario, concordance file", r, "schedule=='B0' :: unique(n_trials)", x[[r]], fint(x[[r]])); fint(x[[1]]) })
  agr <- drange("rationale/pillar2_products_B0.csv", "TRUE", "agree", 1, "%", "decision agreement AUC0-last vs NCA AUC0-inf rule A (ii), product scenarios, two models")
  na_ <- local({ x <- range(rows("rationale/pillar2_products_B0.csv", "TRUE")$n_trials)
    dderived("trials per product scenario, agreement with rule A (ii), range over scenarios and models", "rationale/pillar2_products_B0.csv", "TRUE :: range(n_trials)", x,
             if (x[1] == x[2]) fint(x[1]) else sprintf("%s~%s", fint(x[1]), fint(x[2]))) })
  NY <- TY + TH + 0.07
  premise(local({ r <- range(rows("rationale/pillar2_products_B0.csv", "TRUE")$n_trials); x <- as.numeric(gsub(",", "", n_conc)); x >= r[1] && x <= r[2] }),
          "the concordance files' trial count lies inside the pillar2 range, so one range covers all three quantities (text)")
  deck_text(tx("S17.na", list(cor = cor, agr = agr, agr_b = agr_b, na = na_)), c(XR, NY, WR, GEO$BODY_BOTTOM - NY), size = 16, label = "text_na", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 4)

  # ---- 노트 ----
  rng <- function(ep, am = "M1") drange(SSf, wb(ep, am), "bias_pct", 2, "%", sprintf("bias range, 16 boundary cells, %s, %s", ep, am))
  rcb <- function(cf) drange(CB, sprintf("analysis_model=='M1' & config=='%s'", cf), "bias_pct", 2, "%", sprintf("bias range, 16 boundary cells, criteria set config %s, M1", cf))
  gs00 <- function(ep) { r <- rows(SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & scenario=='S00'", ep)); x <- range(exp(r$mean_log_gmr))
    dderived(sprintf("identical products: geometric mean of trial GMRs, two models, %s, M1", ep), SSf, sprintf("analysis_model=='M1' & endpoint=='%s' & scenario=='S00' :: range(exp(mean_log_gmr))", ep), x, rng_fmt(x[1], x[2], 4)) }
  deck_notes(tx("S17.notes", list(
    ntr = ntr, npop = npop, f97 = f97, lo = lo, hi = hi,
    ntr_v2 = dint(SSf, "analysis_model=='M1' & endpoint=='AUClast' & pk_model=='k2020' & scenario=='V2_up_080'", "n_trials", "trials, 2020 model V2 cell"),
    n_arm = f_n_arm(), wt = f_wt_range(),
    r_true = rng("AUCinf_true"), r_last = rng("AUClast"), r_last0 = rng("AUClast", "M0"), r_ai = rng("AUCinf_Ai"), r_aii = rng("AUCinf_A"), r_b = rng("AUCinf_B"), r_ci = rng("AUCinf_Ci"),
    r_aiii = rcb("A_iii"), r_aiv = rcb("A_iv"), r_cii = rcb("C_ii"), r_ciii = rcb("C_iii"), r_civ = rcb("C_iv"),
    k_true = tw_n("AUCinf_true"), k_true0 = tw_n("AUCinf_true", "M0"), k_last = tw_n("AUClast"), k_last0 = tw_n("AUClast", "M0"),
    k_ai = tw_n("AUCinf_Ai"), k_ai0 = tw_n("AUCinf_Ai", "M0"), k_aii = tw_n("AUCinf_A"), k_aii0 = tw_n("AUCinf_A", "M0"),
    k_b = tw_n("AUCinf_B"), k_b0 = tw_n("AUCinf_B", "M0"), k_ci = kci, k_ci0 = tw_n("AUCinf_Ci", "M0"),
    s_true = gs00("AUCinf_true"), s_last = gs00("AUClast"), s_ai = gs00("AUCinf_Ai"), s_b = gs00("AUCinf_B"),
    v2 = s17_signed(v2), a = f$a, b = f$b, c = f$c,
    g_last = gm("AUClast", "k2016", "Vmax_up_080"), g_ai = gm("AUCinf_Ai", "k2016", "Vmax_up_080"), g_b = gm("AUCinf_B", "k2016", "Vmax_up_080"),
    g_true = gm("AUCinf_true", "k2016", "Vmax_up_080"), t_v = tr("k2016", "Vmax_up_080"), n = nb,
    agr = agr, agr_b = agr_b, n_conc = n_conc, cor = cor, na = na_)))
  deck_end()
}
