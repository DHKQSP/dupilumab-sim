# S07 방법 · 모델 적격성: 연구 제형(300 mg, 600 mg = 2 x 300 mg) 검증 gate, 200 mg 1.14 mL 제형 외부 점검과 흡수 진단(채택 안 함),
# 300 mg 단일 arm(Li 2020)의 모의/관측 산술평균 비(두 모델), 완전 외부 재판정, PKM12350과 가정 체중.
# 그림: 자료별 AUC0-last 모의/관측(두 모델)과 허용 범위. 단일 arm 비는 산술평균 비다(관측 기하평균 미보고). 연구 간 이질성 통계량은 계산된 적이 없다.
# 문자열 열의 k번째 수(예: 제형 "1.14 mL of 175 mg/mL")를 출처와 함께 읽는다
s07_num <- function(rel, where, col, k, n_expected, item) {
  r <- row1(rel, where); rx <- "[0-9]+(\\.[0-9]+)?"; m <- regmatches(r[[col]], gregexpr(rx, r[[col]]))[[1]]
  premise(length(m) == n_expected, sprintf("%s [%s] %s holds %d numbers", rel, where, col, n_expected))
  dderived(item, rel, sprintf("%s :: %s, number %d of %d (regex %s)", where, col, k, n_expected, rx), as.numeric(m[k]), m[k])
}
slide_S07 <- function() {
  Q16 <- "step1/step1b_quant_gate.csv"; Q20 <- "results/step1_k2020/step1b_quant_gate.csv"
  A16 <- "step1/step1f_300mg_arm_check.csv"; A20 <- "step1_k2020/step1f_300mg_arm_check.csv"
  H16 <- "step1/step1h_external_only_gate.csv"; H20 <- "results/step1_k2020/step1h_external_only_gate.csv"
  WS <- "step1/step1i_arm_weight_sensitivity.csv"; AB <- "step1/appendix_absorption_diagnostic.csv"
  DSM <- "step1/appendix_dose_shape_model.csv"; DSO <- "step1/appendix_dose_shape_observed.csv"; GS <- "step1/step1_status.csv"
  deck_slide("S07", tag = "litsim")
  q16 <- rows(Q16); q20 <- rows(Q20); a16 <- rows(A16); a20 <- rows(A20)
  premise(all(q16[gate_role == "gate", pass_mean & pass_cv]) && all(q20[gate_role == "gate", pass_mean & pass_cv]), "every study-presentation data set passes mean and CV criteria in both models")
  premise(all(rows(GS)[role == "gate", status] == "PASS") && all(rows("results/step1_k2020/step1_status.csv")[role == "gate", status] == "PASS"), "gate overall PASS in both models")
  premise(all(q16[gate_role == "external", tmax_sim_median > tmax_obs_median]), "200 mg: simulated tmax later than observed")
  tol <- dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "validation criterion: simulated AUC0-last mean within +-x% of observed (%)", num_fmt(0))
  d200 <- drange(Q16, "gate_role=='external'", "dose_mg", 0, "", "dose of the other presentation (mg)")
  deck_kicker(tx("S07.kicker")); deck_title(tx("S07.title", list(tol = tol, d200 = d200)))
  to <- drange(Q16, "gate_role=='external'", "tmax_obs_median", 1, "", "observed median tmax, other presentation (days)")
  ts <- drange(Q16, "gate_role=='external'", "tmax_sim_median", 0, "", "simulated median tmax, other presentation (days)")

  # ---- 그림: 자료별 AUC0-last 모의/관측, 두 모델, 허용 범위 ----
  L <- DK$txt$S07$fig
  tolv <- .read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct / 100
  arm_lab <- function(st, ar) sprintf("%s %s", st, unlist(L$arm)[ar])
  g1 <- merge(q16[gate_role == "gate", .(id, dose_mg, k2016 = AUClast_ratio)], q20[, .(id, k2020 = AUClast_ratio)], by = "id", sort = FALSE)
  g1[, lab := sprintf("%s mg %s", fnum(dose_mg, 0), vapply(id, function(z) L$ds[[z]], ""))]
  g2 <- merge(a16[, .(study, arm, k2016 = ratio)], a20[, .(study, arm, k2020 = ratio)], by = c("study", "arm"), sort = FALSE)
  g2[, lab := arm_lab(study, arm)]
  g3 <- merge(q16[gate_role == "external", .(id, dose_mg, k2016 = AUClast_ratio)], q20[, .(id, k2020 = AUClast_ratio)], by = "id", sort = FALSE)
  g3[, lab := vapply(id, function(z) L$ds[[z]], "")]
  premise(nrow(g1) == 5 && nrow(g2) == 6 && nrow(g3) == 4 && !anyNA(c(g1$k2020, g2$k2020, g3$k2020)), "5 study-presentation sets, 6 single arms, 4 external 200 mg arms, both models")
  premise(all(g2$k2020 > g2$k2016), "2020 model ratio above 2016 model ratio in every single arm")
  hdr <- c(fill(L$g1, list()), fill(L$g2, list(dose = fnum(.read("config/design_clot2021.yaml")$variability_checks$li2020$dose_mg, 0))), fill(L$g3, list(dose = fnum(unique(g3$dose_mg), 0))))
  G <- list(g1, g2, g3); y <- 0; items <- list(); heads <- list()
  for (k in seq_along(G)) { y <- y - 1; heads[[k]] <- data.table(y = y, lab = hdr[k]); g <- G[[k]]
    for (i in seq_len(nrow(g))) { y <- y - 1; items[[length(items) + 1L]] <- data.table(y = y, lab = g$lab[i], k2016 = g$k2016[i], k2020 = g$k2020[i]) }
    y <- y - 0.35 }
  I <- rbindlist(items); H <- rbindlist(heads)
  P <- melt(I, id.vars = c("y", "lab"), measure.vars = c("k2016", "k2020"), variable.name = "model", value.name = "ratio")
  P[, model := factor(model_lab()[as.character(model)], levels = model_lab())]
  xl <- c(0.7, 1.35); FW <- 5.5
  # 허용 범위 음영과 기준선 1은 자료 행이 있는 구간에만 그린다(묶음 제목 줄은 비워 제목이 음영을 가리지 않게)
  I[, grp := rep(seq_along(G), vapply(G, nrow, 1L))]
  GB <- I[, .(y0 = min(y) - 0.5, y1 = max(y) + 0.5), by = grp]
  premise(max(g3$k2016, g3$k2020) < 1 + tolv, "200 mg points lie inside the band, left of the tmax note (right of the band)")
  p <- ggplot() +
    geom_rect(data = GB, aes(xmin = 1 - tolv, xmax = 1 + tolv, ymin = y0, ymax = y1), fill = PAL$tint_grey) +
    geom_segment(data = GB, aes(x = 1, xend = 1, y = y0, yend = y1), colour = PAL$muted, linewidth = 0.5) +
    geom_segment(data = I, aes(x = k2016, xend = k2020, y = y, yend = y), colour = PAL$muted, linewidth = 0.6) +
    geom_point(data = P, aes(ratio, y, colour = model, shape = model), size = 3) +
    geom_text(data = H, aes(x = xl[1], y = y, label = lab), hjust = 0, vjust = 0.5, size = 4.3, fontface = "bold", family = FONT, colour = PAL$ink) +
    annotate("text", x = 1, y = max(H$y) + 1.25, label = fill(L$band, list(tol = tol)), hjust = 0.5, vjust = 0.5, size = 3.9, family = FONT, colour = PAL$ink2) +
    # tmax 설명은 허용 범위 음영 오른쪽 흰 바탕에 둔다(음영 경계와 세로 격자선을 가리지 않게 흰 상자)
    annotate("label", x = 1 + tolv + 0.012, y = mean(I[grp == length(G), y]), label = fill(L$tmax, list(to = to, ts = ts)), hjust = 0, vjust = 0.5, size = 3.9, family = FONT,
             colour = PAL$ink2, lineheight = 0.95, fill = "white", label.size = 0, label.padding = unit(0.12, "lines")) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = unname(MODEL_SHAPE)) +
    scale_y_continuous(breaks = I$y, labels = I$lab, limits = c(min(I$y) - 0.5, max(H$y) + 1.75), expand = expansion(0)) +
    scale_x_continuous(limits = xl, breaks = seq(0.7, 1.3, by = 0.1), expand = expansion(0)) +
    labs(x = L$xlab, y = NULL) + theme_deck(13) +
    # 범례는 패널 오른쪽 끝에 맞춘다: 두 모델 이름이 패널 폭보다 길어 가운데 정렬이면 그림 오른쪽 끝에서 잘린다(왼쪽 y축 글자 위 빈 곳으로 넘친다)
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 12, colour = PAL$ink), legend.margin = margin(0, 0, 0, 0),
          legend.justification = "right")
  deck_figure(p, "s07_validation_ratios", c(GEO$ML, GEO$BODY_TOP, FW, GEO$BODY_BOTTOM - GEO$BODY_TOP), src = c(Q16, Q20, A16, A20))

  # ---- 요점 ----
  fct <- list(
    tol = tol, d200 = d200, dose = f_dose(),
    gr16 = drange(Q16, "gate_role=='gate'", "AUClast_ratio", 2, "", "AUC0-last sim/obs range, 2016, study presentation"),
    gr20 = drange(Q20, "gate_role=='gate'", "AUClast_ratio", 2, "", "AUC0-last sim/obs range, 2020, study presentation"),
    gc16 = drange(Q16, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2016"),
    gc20 = drange(Q20, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2020"),
    vol = s07_num(Q16, "id=='PKM14271_200mg_test'", "presentation", 1, 2, "other presentation: volume (mL)"),
    conc = s07_num(Q16, "id=='PKM14271_200mg_test'", "presentation", 2, 2, "other presentation: concentration (mg/mL)"),
    to = to, ts = ts,
    r16 = drange(Q16, "gate_role=='external'", "AUClast_ratio", 2, "", "AUC0-last sim/obs, other presentation, 2016"),
    r20 = drange(Q20, "gate_role=='external'", "AUClast_ratio", 2, "", "AUC0-last sim/obs, other presentation, 2020"),
    km = drange(DSM, "ka_multiplier %in% c(1.5, 2)", "ka_multiplier", 1, "", "absorption diagnostic: ka multipliers"),
    tk = drange(AB, "condition %in% c('cohen_like','pkm14271_like') & ka_multiplier %in% c(1.5, 2)", "tmax_sim_median", 0, "", "absorption diagnostic: simulated tmax at ka x1.5 to x2 (days)"),
    s1 = dv(DSM, "ka_multiplier==1", "ratio_per_mg", 2, "", "dose-normalized 200:300 AUC0-last ratio, simulated, ka x1"),
    s2 = drange(DSM, "ka_multiplier %in% c(1.5, 2)", "ratio_per_mg", 2, "", "dose-normalized 200:300 AUC0-last ratio, simulated, ka x1.5 to x2"),
    o = dv(DSO, "!grepl('PKM14271', basis)", "ratio_per_mg", 2, "", "dose-normalized 200:300 AUC0-last ratio, observed, four arms n-weighted"),
    a16 = drange(A16, "TRUE", "ratio", 2, "", "300 mg single arms: mean AUC0-last sim/obs, 2016"),
    a20 = drange(A20, "TRUE", "ratio", 2, "", "300 mg single arms: mean AUC0-last sim/obs, 2020"),
    n16 = dcount(H16, "TRUE", "fully external items, 2016"), f16 = dcount(H16, "pass_mean==FALSE", "fully external items failing, 2016"),
    n20 = dcount(H20, "TRUE", "fully external items, 2020"), f20 = dcount(H20, "pass_mean==FALSE", "fully external items failing, 2020"),
    x16 = dv(H16, "item=='PKM12350 test'", "ratio", 2, "", "PKM12350 test arm sim/obs, 2016 (fully external judgement)"),
    w0 = dcfg("design_clot2021.yaml", c("arm_checks_300mg", "weight", "mean"), "assumed mean weight of the single arms (kg)", num_fmt(0)),
    p1 = dv(A16, "study=='PKM12350' & arm=='test'", "auclast_obs", 1, "", "PKM12350 test arm observed mean AUC0-last"),
    p2 = dv(A16, "study=='PKM12350' & arm=='reference'", "auclast_obs", 1, "", "PKM12350 reference arm observed mean AUC0-last"),
    pmx = dext(A16, "TRUE", "auclast_obs", max, 1, "", "largest observed 300 mg single-arm mean AUC0-last"))
  # Cmax 비: 최대 코호트(이름은 그림 문구에서, 용량은 자료에서), 과대 예측이 아닌 값은 많아야 하나
  qg <- q16[gate_role == "gate"]; cid <- qg[which.max(Cmax_ratio), id]
  fct$cmx <- sprintf("%s mg %s", dv(Q16, sprintf("id=='%s'", cid), "dose_mg", 0, "", "dose of the cohort with the largest Cmax sim/obs, 2016"), sub(" \\(.*\\)$", "", L$ds[[cid]]))
  premise(sum(c(qg$Cmax_ratio, q20[gate_role == "gate", Cmax_ratio]) <= 1) <= 1, "Cmax over-predicted in all but at most one study-presentation cohort over both models (text: mostly over-predicted)")
  premise(nrow(a16) == 6, "six 300 mg single arms (text: lowest of six)")
  premise(length(list.files(proj_path("results"), pattern = "heterogen", recursive = TRUE, ignore.case = TRUE)) == 0, "no between-study heterogeneity statistic in the result files (text: not computed)")
  eo <- .read("config/design_clot2021.yaml")$external_only_gate
  premise(!setequal(unlist(eo$k2016), unlist(eo$k2020)), "fully external item lists differ by model (text: different numbers of items)")
  fl <- rbind(rows(H16, "pass_mean==FALSE"), rows(H20, "pass_mean==FALSE"))
  premise(all(fl$type == "Li 2020 arm"), "every fully external failure is a Li 2020 single arm (no reported weight)")
  premise(identical(a16[order(auclast_obs)][1:2, study], c("PKM12350", "PKM12350")), "PKM12350 arms have the two lowest observed means of the six 300 mg arms")
  wsr <- rows(WS, "model=='k2016' & study=='PKM12350' & arm=='test'"); premise(nrow(wsr) >= 2 && all(diff(wsr[order(weight_mean), ratio]) < 0), "PKM12350 test ratio falls as the assumed weight rises (2016 model)")
  wlo <- min(wsr$weight_mean); whi <- max(wsr$weight_mean)
  premise(abs(wsr[weight_mean == wlo, ratio] - 1) > tolv && abs(wsr[weight_mean == whi, ratio] - 1) <= tolv, "PKM12350 test fails the gate mean criterion at the lowest assumed weight and passes at the highest (text: pass or fail depends on the assumed weight)")
  fct$wl <- dv(WS, sprintf("model=='k2016' & study=='PKM12350' & arm=='test' & weight_mean==%s", wlo), "weight_mean", 0, "", "assumed weight, lowest (kg)")
  fct$wh <- dv(WS, sprintf("model=='k2016' & study=='PKM12350' & arm=='test' & weight_mean==%s", whi), "weight_mean", 0, "", "assumed weight, highest (kg)")
  fct$rl <- dv(WS, sprintf("model=='k2016' & study=='PKM12350' & arm=='test' & weight_mean==%s", wlo), "ratio", 2, "", "PKM12350 test sim/obs at the lowest assumed weight, 2016")
  fct$rh <- dv(WS, sprintf("model=='k2016' & study=='PKM12350' & arm=='test' & weight_mean==%s", whi), "ratio", 2, "", "PKM12350 test sim/obs at the highest assumed weight, 2016")
  deck_bullets(tx("S07.bullets", fct), box = c(GEO$ML + FW + 0.25, GEO$BODY_TOP, GEO$W - GEO$MR - (GEO$ML + FW + 0.25), GEO$BODY_BOTTOM - GEO$BODY_TOP), size = 16, gap_pt = 4)

  # ---- 노트 ----
  deck_notes(tx("S07.notes", list(
    tol = tol, cvr = dcfg("design_clot2021.yaml", c("gate", "log_cv_range_pct"), "validation criterion: simulated log-scale CV range (%)", function(x) rng_fmt(x[1], x[2], 0)),
    cv16 = drange(Q16, "gate_role=='gate'", "AUClast_sim_logcv", 0, "", "simulated log-scale CV range, 2016"),
    cv20 = drange(Q20, "gate_role=='gate'", "AUClast_sim_logcv", 0, "", "simulated log-scale CV range, 2020"),
    z = drange(Q16, "gate_role=='external'", "AUClast_z", 1, "", "z of AUC0-last, other presentation, 2016"),
    m16 = dv("step1/step1f_300mg_arm_check_summary.csv", "TRUE", "model_mean_pooled", 0, "", "pooled simulated mean AUC0-last, 300 mg, 2016"),
    m20 = dv("step1_k2020/step1f_300mg_arm_check_summary.csv", "TRUE", "model_mean_pooled", 0, "", "pooled simulated mean AUC0-last, 300 mg, 2020"),
    rng = dspan("step1/step1f_300mg_arm_check_summary.csv", "TRUE", "arm_min", "arm_max", 1, "", "observed 300 mg single-arm mean range"),
    x16 = fct$x16, w0 = fct$w0,
    x78 = dv(WS, sprintf("model=='k2016' & study=='PKM12350' & arm=='test' & weight_mean==%s", .read("config/design_clot2021.yaml")$arm_checks_300mg$weight$mean), "ratio", 2, "", "PKM12350 test sim/obs at the assumed weight, weight-sensitivity run, 2016"),
    f20 = fct$f20, n20 = fct$n20, n16 = fct$n16, f16 = fct$f16, cmx = fct$cmx,
    cmin20 = dext(Q20, "gate_role=='gate'", "Cmax_ratio", min, 3, "", "smallest Cmax sim/obs, 2020, study presentation"))))
  deck_end()
}
