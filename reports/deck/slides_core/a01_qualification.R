# A1 별첨 · 모델 적격성(S3 본문 "모델은 연구 제형 300 mg의 관측 AUClast를 ±15% 안에서 재현한다"의 근거). [문헌+모의]
# 그림: 자료별 AUClast 모의/관측(평균의 비), 두 모델(표식 모양), gate 허용 범위 ±tol 음영(config/design_clot2021.yaml gate.auclast_mean_tol_pct).
#   위 = 연구 제형(300 mg 2 mL of 150 mg/mL, 600 mg = 2 x 300 mg; gate 판정 5개), 아래 = 200 mg 1.14 mL 제형(외부 점검 4개, gate 제외, 회색).
#   속이 빈 표식 = 그 모델의 개발 자료(config/design_clot2021.yaml datasets dev_k2016/dev_k2020가 internal 또는 partial).
# 오른쪽 본문: 연구 제형(300 mg 부분집합은 S3와 같은 범위), 200 mg 제형, 완전 외부 재판정(step1h, step1_status: 두 모델 FAIL, D-030).
# 자료 논리는 결과보고 덱 S07(step1b_quant_gate.csv, step1f, step1h)과 A1(literature_numeric.csv, step1g)을 그대로 쓴다. 세부(Cmax, 200 mg 흡수,
# 300 mg 단일 arm, 가정 체중 민감도(step1i), 문헌 평균비)는 노트에 둔다.

# 문자열 열의 k번째 수(예: 제형 "2 mL of 150 mg/mL")를 출처와 함께 읽는다
a1_num <- function(rel, where, col, k, n_expected, item) {
  r <- row1(rel, where); rx <- "[0-9]+(\\.[0-9]+)?"; m <- regmatches(r[[col]], gregexpr(rx, r[[col]]))[[1]]
  premise(length(m) == n_expected, sprintf("%s [%s] %s holds %d numbers", rel, where, col, n_expected))
  dderived(item, rel, sprintf("%s :: %s, number %d of %d (regex %s)", where, col, k, n_expected, rx), as.numeric(m[k]), m[k])
}
# 사전 명시 등록부 한 행의 커밋(Date (evidence) 열)
a1_commit <- function(item_rx, item) {
  PR <- "regulatory/tables/prespecification_register.csv"; w <- sprintf("grepl('%s', Item)", item_rx)
  r <- row1(PR, w); m <- regmatches(r[["Date (evidence)"]], regexec("commit ([0-9a-f]{7})", r[["Date (evidence)"]]))[[1]]
  premise(length(m) == 2, sprintf("commit hash in register row [%s]", w))
  dderived(item, PR, sprintf("%s :: Date (evidence), regex commit ([0-9a-f]{7})", w), m[2], m[2])
}

slide_A1 <- function() {
  Q16 <- "step1/step1b_quant_gate.csv"; Q20 <- "step1_k2020/step1b_quant_gate.csv"
  A16 <- "step1/step1f_300mg_arm_check.csv"; A20 <- "step1_k2020/step1f_300mg_arm_check.csv"
  H16 <- "step1/step1h_external_only_gate.csv"; H20 <- "step1_k2020/step1h_external_only_gate.csv"
  LN <- "literature/literature_numeric.csv"; PR <- "regulatory/tables/prespecification_register.csv"
  deck_slide("A1", tag = "litsim")
  L <- DK$txt$A1$fig; ML <- DK$txt$common$models_short
  q16 <- rows(Q16); q20 <- rows(Q20)
  tolv <- .read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct

  # ---- 전제: 두 모델 같은 자료, 연구 제형 전부 합격(평균·CV), gate 종합 PASS, 200 mg 제형은 외부 점검 ----
  premise(identical(q16$id, q20$id) && identical(q16$gate_role, q20$gate_role), "same data sets and roles in both models")
  for (q in list(q16, q20)) {
    premise(all(q[gate_role == "gate", pass_mean & pass_cv]), "every study-presentation data set passes the mean and CV criteria")
    premise(all(abs(q[gate_role == "gate", AUClast_ratio] - 1) * 100 <= tolv), "every study-presentation AUClast ratio lies within the tolerance (title)")
    premise(all(q[gate_role == "external", tmax_sim_median] > q[gate_role == "external", tmax_obs_median]), "200 mg presentation: simulated tmax later than observed (text: fast absorption not reproduced)")
  }
  for (s_ in c("step1/step1_status.csv", "step1_k2020/step1_status.csv")) premise(all(rows(s_)[role == "gate", status] == "PASS"), sprintf("%s: every gate criterion and the overall gate PASS", s_))
  premise(all(q16[gate_role == "gate" & dose_mg == 300, presentation] == "2 mL of 150 mg/mL") && all(q16[gate_role == "gate" & dose_mg == 600, presentation] == "2 x 300 mg") &&
          setequal(unique(q16[gate_role == "gate", dose_mg]), c(300, 600)), "study presentation: 300 mg = 2 mL of 150 mg/mL, 600 mg = two 300 mg injections")
  premise(all(startsWith(q16[gate_role == "external", presentation], "1.14 mL of 175 mg/mL")) && all(q16[gate_role == "external", dose_mg] == 200), "external rows: 200 mg, 1.14 mL of 175 mg/mL")
  premise(as.numeric(.read("config/trial_design.yaml")$dose_mg) == 300, "study dose is the 300 mg of the study presentation")
  reg <- rows(PR)
  premise(nrow(reg[grepl("^Scope of model validation restricted to the study presentation", Item) & grepl("^post hoc", Status)]) == 1, "register: restricting the gate to the study presentation was post hoc (text: after the first comparison)")
  premise(nrow(reg[grepl("^Model validation criteria", Item) & Status == "pre-specified"]) == 1 && nrow(reg[grepl("^Models and parameter values", Item) & Status == "pre-specified"]) == 1,
          "register: validation criteria and parameter values pre-specified (caption: criteria fixed before the comparison, parameters not adjusted)")

  tol <- dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "gate tolerance for the simulated AUClast mean (%)", num_fmt(0))
  g_all <- c(q16[gate_role == "gate", AUClast_ratio], q20[gate_role == "gate", AUClast_ratio])
  rng <- dderived("simulated / observed AUClast mean, study presentation, range over both models", Q16,
                  sprintf("gate_role=='gate' :: range(AUClast_ratio) over %s and %s", Q16, Q20), range(g_all), rng_fmt(min(g_all), max(g_all), 2))
  y0 <- core_title(tx("A1.title", list(tol = tol)), tx("A1.kicker"))

  # ---- 그림: 자료별 모의/관측, 두 모델(표식), 허용 범위 음영 ----
  premise(all(q16$id %in% names(L$ds)), "a label for every data set")
  d <- rbind(q16[, .(id, gate_role, dose_mg, pk = "k2016", ratio = AUClast_ratio)], q20[, .(id, gate_role, dose_mg, pk = "k2020", ratio = AUClast_ratio)])
  d[, row := match(id, q16$id)][, mod := factor(unlist(ML[pk]), levels = unlist(ML))]
  # 모델 개발 자료 표시(config/design_clot2021.yaml datasets의 dev_k2016/dev_k2020: internal 또는 partial이면 속이 빈 표식)
  ds_cfg <- .read("config/design_clot2021.yaml")$datasets; dv_st <- rbindlist(lapply(ds_cfg, function(z) data.table(id = z$id, k2016 = z$dev_k2016, k2020 = z$dev_k2020)))
  premise(all(d$id %in% dv_st$id) && all(unlist(dv_st[, .(k2016, k2020)]) %in% c("external", "internal", "partial")), "development-data status recorded for every data set")
  d[, dev := mapply(function(i, m) dv_st[id == i][[m]] != "external", id, pk)]
  premise(!any(d[gate_role == "external", dev]), "no 200 mg row is development data (only gate rows are marked)")
  premise(setequal(d[dev == TRUE & pk == "k2016", id], "clot300_nonasian_pooled") && setequal(d[dev == TRUE & pk == "k2020", id], c("clot300_nonasian_pooled", "clot300_japanese_TDU12265", "clot600_japanese")) &&
          all(unlist(dv_st[id == "clot300_nonasian_pooled", .(k2016, k2020)]) == "partial"), "development data: pooled non-Asian partly for both models, Japanese rows for the 2020 model (caption)")
  dsrc("development-data status of the validation data sets (hollow markers)", "config/design_clot2021.yaml")
  d[, lab := vapply(seq_len(.N), function(i) if (gate_role[i] == "gate") fill(L$gate_lab, list(dose = fnum(dose_mg[i], 0), ds = L$ds[[id[i]]])) else L$ds[[id[i]]], "")]
  w <- dcast(d, id + gate_role + row + lab ~ pk, value.var = "ratio")[, vtxt := sprintf("%s / %s", fnum(k2016, 2), fnum(k2020, 2))]
  XL <- c(0.78, 1.20); XV <- 1.225
  premise(min(d$ratio) > XL[1] && max(d$ratio) < XL[2], "every ratio inside the plotted range")
  d200 <- unique(q16[gate_role == "external", dose_mg])
  mk <- function(role, sub, xlab) {
    dd <- d[gate_role == role]; ww <- w[gate_role == role]; lv <- rev(ww[order(row), lab])
    dd[, labf := factor(lab, levels = lv)]; ww[, labf := factor(lab, levels = lv)]
    col <- if (role == "gate") PAL$blue else PAL$ink2
    p <- ggplot() +
      annotate("rect", xmin = 1 - tolv / 100, xmax = 1 + tolv / 100, ymin = -Inf, ymax = Inf, fill = PAL$tint_blue) +
      geom_vline(xintercept = 1, colour = PAL$muted, linewidth = 0.5) +
      geom_segment(data = ww, aes(x = k2016, xend = k2020, y = labf, yend = labf), colour = PAL$muted, linewidth = 0.7) +
      geom_point(data = dd, aes(ratio, labf, shape = mod, fill = dev), colour = col, size = 3.8, stroke = 1.1) +
      geom_text(data = ww, aes(x = XV, y = labf, label = vtxt), hjust = 0, size = PT(15), family = FONT, colour = PAL$ink) +
      annotate("text", x = XV, y = Inf, label = L$vals, hjust = 0, vjust = -0.5, size = PT(14), family = FONT, colour = PAL$ink2) +
      scale_shape_manual(values = setNames(c(21, 24), unlist(ML)), name = NULL, drop = FALSE) +                       # 21/24 = CORE_MODEL_SHAPE의 원·삼각형(속 채움 조절용)
      scale_fill_manual(values = c(`FALSE` = col, `TRUE` = "white"), breaks = "TRUE", labels = L$dev, name = NULL, drop = FALSE) +
      scale_x_continuous(breaks = c(1 - tolv / 100, 1, 1 + tolv / 100), labels = function(v) fnum(v, 2)) +
      scale_y_discrete(expand = expansion(add = 0.6)) +
      coord_cartesian(xlim = XL, clip = "off") +
      labs(x = xlab, y = NULL, subtitle = sub) + theme_core(16) +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 15, colour = PAL$ink),
            plot.subtitle = element_text(size = 16, face = "bold", colour = if (role == "gate") PAL$blue else PAL$ink2, margin = margin(0, 0, 4, 0)),
            plot.title.position = "plot", plot.margin = margin(2, 96, 2, 4))
    p
  }
  # 범례(모델 = 표식 모양)는 위 그림에만(제목 줄 아래 왼쪽; 오른쪽 값 열 머리글과 떨어지게)
  # 범례(모델 = 표식 모양, 속이 빈 표식 = 개발 자료)는 그림 아래(오른쪽 값 열 머리글과 겹치지 않게)
  p1 <- mk("gate", L$g_gate, NULL) + guides(shape = guide_legend(order = 1, override.aes = list(colour = PAL$ink, fill = PAL$ink)),
                                             fill = guide_legend(order = 2, override.aes = list(shape = 21, colour = PAL$ink, size = 3.8))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- mk("external", fill(L$g_ext, list(dose = fnum(d200, 0))), L$xlab) + guides(shape = "none", fill = "none")
  p <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(nrow(w[gate_role == "gate"]), nrow(w[gate_role == "external"]))) +
    patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
                                                       legend.box.spacing = grid::unit(2, "pt"), legend.spacing.x = grid::unit(2, "pt"), legend.text = element_text(size = 14, colour = PAL$ink))

  # ---- 본문 두 줄, 캡션 ----
  to <- drange(Q16, "gate_role=='external'", "tmax_obs_median", 1, "", "observed median tmax, 200 mg presentation (days)")
  ts <- drange(Q16, "gate_role=='external'", "tmax_sim_median", 0, "", "simulated median tmax, 200 mg presentation, 2016 model (days)")
  premise(identical(sort(q16[gate_role == "external", tmax_sim_median]), sort(q20[gate_role == "external", tmax_sim_median])), "same simulated tmax in both models (text states one value)")
  fb <- list(conc = a1_num(Q16, "id=='clot300_chinese'", "presentation", 2, 2, "study presentation: concentration (mg/mL)"),
             vol = a1_num(Q16, "id=='clot300_chinese'", "presentation", 1, 2, "study presentation: volume (mL)"),
             d600 = drange(Q16, "gate_role=='gate' & dose_mg > 300", "dose_mg", 0, "", "higher study-presentation dose (mg)"),
             n = dcount(Q16, "gate_role=='gate'", "study-presentation data sets"),
             r16 = drange(Q16, "gate_role=='gate'", "AUClast_ratio", 2, "", "AUClast sim/obs range, 2016 model, study presentation"),
             r20 = drange(Q20, "gate_role=='gate'", "AUClast_ratio", 2, "", "AUClast sim/obs range, 2020 model, study presentation"),
             d200 = drange(Q16, "gate_role=='external'", "dose_mg", 0, "", "dose of the other presentation (mg)"),
             conc2 = a1_num(Q16, "id=='PKM14271_200mg_test'", "presentation", 2, 2, "other presentation: concentration (mg/mL)"),
             vol2 = a1_num(Q16, "id=='PKM14271_200mg_test'", "presentation", 1, 2, "other presentation: volume (mL)"),
             to = to, ts = ts,
             d300 = drange(Q16, "gate_role=='gate' & dose_mg==300", "dose_mg", 0, "", "study-presentation single-injection dose (mg)"),
             n3 = dcount(Q16, "gate_role=='gate' & dose_mg==300", "study-presentation 300 mg data sets"),
             r16_3 = drange(Q16, "gate_role=='gate' & dose_mg==300", "AUClast_ratio", 2, "", "AUClast sim/obs range, 2016 model, 300 mg study presentation (as S3)"),
             r20_3 = drange(Q20, "gate_role=='gate' & dose_mg==300", "AUClast_ratio", 2, "", "AUClast sim/obs range, 2020 model, 300 mg study presentation (as S3)"))
  # 배치: 왼쪽 그림(본문 높이 전체), 오른쪽 요점 두 문단, 아래 캡션
  h16 <- rows(H16); h20 <- rows(H20); fl_ <- rbind(h16[pass_mean == FALSE], h20[pass_mean == FALSE])
  premise(nrow(fl_) > 0 && all(fl_$type == "Li 2020 arm"), "every fully external failure is a Li 2020 single arm (no reported weight; caption)")
  for (s_ in c("step1/step1_status.csv", "step1_k2020/step1_status.csv")) premise(startsWith(rows(s_)[role == "gate(외부만)", status], "FAIL"), sprintf("%s: fully external re-judgement FAIL (caption)", s_))
  premise(isTRUE(.read("config/design_clot2021.yaml")$arm_checks_300mg$weight$status == "assumption"), "Li 2020 arm weight is an assumption (caption)")
  w0 <- dcfg("design_clot2021.yaml", c("arm_checks_300mg", "weight", "mean"), "assumed mean weight of the single arms (kg)", num_fmt(0))
  xo <- list(n16 = dcount(H16, "TRUE", "fully external items, 2016 model"), f16 = dcount(H16, "pass_mean==FALSE", "fully external items failing, 2016 model"),
             n20 = dcount(H20, "TRUE", "fully external items, 2020 model"), f20 = dcount(H20, "pass_mean==FALSE", "fully external items failing, 2020 model"))
  # 배치: 왼쪽 = 그림 + 그 아래 캡션(그림 폭), 오른쪽 = 본문 세 문단(연구 제형, 200 mg 제형, 완전 외부 재판정)
  FW <- 6.75; XR <- GEO$ML + FW + 0.3; WR <- GEO$W - GEO$MR - XR
  cap <- tx("A1.caption", list(tol = tol, d200 = fb$d200))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, x = GEO$ML, width = FW, size = 14)
  body <- tx("A1.body", c(fb, xo, list(w0 = w0))); bh <- core_body_h(body, WR, gap_pt = 12)
  premise(y0 + 0.1 + bh <= GEO$BODY_BOTTOM, "right-column body fits above the footer")
  deck_text(body, c(XR, y0 + 0.1, WR, bh), size = SZ$body, label = "body", gap_pt = 12)
  deck_figure(p, "a1_qualification", c(GEO$ML, y0, FW, capy - 0.12 - y0), src = c(Q16, Q20))

  # ---- 노트: 기준, 자료별 값, 200 mg, Cmax, 단일 arm, 완전 외부 재판정, 문헌 평균비 ----
  pr2 <- function(id_, item) sprintf("%s / %s", dv(Q16, sprintf("id=='%s'", id_), "AUClast_ratio", 2, "", paste(item, "2016 model")), dv(Q20, sprintf("id=='%s'", id_), "AUClast_ratio", 2, "", paste(item, "2020 model")))
  qg <- q16[gate_role == "gate"]; cid <- qg[which.max(Cmax_ratio), id]
  premise(sum(c(qg$Cmax_ratio, q20[gate_role == "gate", Cmax_ratio]) <= 1) <= 1, "Cmax over-predicted in all but at most one study-presentation cohort over both models (notes: mostly over-predicted)")
  premise(sum(q20[gate_role == "gate", Cmax_ratio] < 1) == 1 && all(q16[gate_role == "gate", Cmax_ratio] > 1), "2020 model: exactly one study-presentation Cmax ratio below 1, none in the 2016 model (notes)")
  eo <- .read("config/design_clot2021.yaml")$external_only_gate
  premise(!setequal(unlist(eo$k2016), unlist(eo$k2020)), "fully external item lists differ by model (notes: different numbers of items)")
  premise(nrow(rows(A16)) == 6 && nrow(rows(A20)) == 6, "six 300 mg single arms")
  # 문헌 평균비: 일치 판정 허용 폭은 결과 파일 판정 열의 '±n%p'에서 읽고, 판정을 수치로 다시 계산해 결과 파일과 같은지 본다
  W <- list(c300 = 'source=="Clot 2021 Table 3" & dose_mg==300', c600 = 'source=="Clot 2021 Table 3" & dose_mg==600', pkm = "grepl('PKM12350', source)")
  ltol <- local({ s <- row1(LN, W$c300)$judgment; m <- regmatches(s, regexec("±([0-9]+)%p", s))[[1]]; premise(length(m) == 2, "tolerance in the judgment label")
    x <- as.numeric(m[2]); list(x = x, p = dderived("agreement tolerance of the NCA mean ratio (points), from the judgment label", LN, sprintf("%s :: judgment, regex +/-([0-9]+)%%p", W$c300), x, fnum(x, 0))) })
  agree <- function(k) { r <- row1(LN, W[[k]]); max(abs(c(r$nca_mean_ratio_k2016, r$nca_mean_ratio_k2020) - r$lit_mean_ratio)) * 100 <= ltol$x + 1e-12 }
  premise(agree("c300") && agree("pkm") && !agree("c600"), "both 300 mg literature rows agree within the tolerance, the 600 mg row does not")
  premise(startsWith(row1(LN, W$c300)$judgment, "일치") && startsWith(row1(LN, W$pkm)$judgment, "일치") && !startsWith(row1(LN, W$c600)$judgment, "일치"), "recomputed agreement equals the stored judgment")
  r600 <- row1(LN, W$c600)
  premise(all(c(r600$true_mean_ratio_k2016, r600$true_mean_ratio_k2020, r600$nca_mean_ratio_k2016, r600$nca_mean_ratio_k2020) < r600$lit_mean_ratio),
          "600 mg: simulated true and NCA mean ratios below the published NCA ratio in both models (notes: conservative direction)")
  lv <- function(k, col, item) dv(LN, W[[k]], col, 1, "%", item, scale = 100)
  two <- function(k, stem, item) sprintf("%s / %s", lv(k, paste0(stem, "_k2016"), paste(item, "2016 model")), lv(k, paste0(stem, "_k2020"), paste(item, "2020 model")))
  wmin <- function(rel, m) { r <- rows(rel, "study=='PKM12350'"); g <- r[, .(pass = all(within_15)), by = weight_mean][order(weight_mean)]
    x <- min(g$weight_mean[g$pass]); premise(any(g$pass) && all(g$pass[g$weight_mean >= x]), sprintf("%s: PKM12350 arms within 15%% from one assumed weight upward", m))
    dderived(sprintf("lowest assumed mean weight with both PKM12350 arms within 15%%, %s (kg)", m), rel, "study=='PKM12350' :: min(weight_mean) with all(within_15)", x, fnum(x, 0)) }
  ws <- list(ws16 = wmin("step1/step1i_arm_weight_sensitivity.csv", "k2016"), ws20 = wmin("step1_k2020/step1i_arm_weight_sensitivity.csv", "k2020"))
  deck_notes(tx("A1.notes", c(ws, list(
    tol = tol, cvr = dcfg("design_clot2021.yaml", c("gate", "log_cv_range_pct"), "validation criterion: simulated log-scale CV range (%)", function(x) rng_fmt(x[1], x[2], 0)),
    cv16 = drange(Q16, "gate_role=='gate'", "AUClast_sim_logcv", 0, "", "simulated log-scale CV range, 2016 model"),
    cv20 = drange(Q20, "gate_role=='gate'", "AUClast_sim_logcv", 0, "", "simulated log-scale CV range, 2020 model"),
    x1 = pr2("clot300_nonasian_pooled", "AUClast sim/obs, non-Asian pooled 300 mg,"), x2 = pr2("clot300_japanese_TDU12265", "AUClast sim/obs, Japanese 300 mg,"),
    x3 = pr2("clot300_chinese", "AUClast sim/obs, Chinese 300 mg,"), x4 = pr2("clot600_chinese", "AUClast sim/obs, Chinese 600 mg,"), x5 = pr2("clot600_japanese", "AUClast sim/obs, Japanese 600 mg,"),
    d200 = fb$d200, d600 = fb$d600, to = to, ts = ts, rng = rng,
    d300 = fb$d300,
    e16 = drange(Q16, "gate_role=='external'", "AUClast_ratio", 2, "", "AUClast sim/obs, 200 mg presentation, 2016 model"),
    e20 = drange(Q20, "gate_role=='external'", "AUClast_ratio", 2, "", "AUClast sim/obs, 200 mg presentation, 2020 model"),
    z = drange(Q16, "gate_role=='external'", "AUClast_z", 1, "", "z of AUClast, 200 mg presentation, 2016 model"),
    cmt = a1_commit("^Scope of model validation restricted", "commit of the post hoc decision restricting the gate to the study presentation"),
    c16 = drange(Q16, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2016 model, study presentation"),
    c20 = drange(Q20, "gate_role=='gate'", "Cmax_ratio", 2, "", "Cmax sim/obs range, 2020 model, study presentation"),
    cmx = sprintf("%s mg %s", dv(Q16, sprintf("id=='%s'", cid), "dose_mg", 0, "", "dose of the cohort with the largest Cmax sim/obs, 2016 model"), L$ds[[cid]]),
    cmin20 = dext(Q20, "gate_role=='gate'", "Cmax_ratio", min, 3, "", "smallest Cmax sim/obs, 2020 model, study presentation"),
    w0 = w0,
    a16 = drange(A16, "TRUE", "ratio", 2, "", "300 mg single arms: mean AUClast sim/obs, 2016 model"),
    a20 = drange(A20, "TRUE", "ratio", 2, "", "300 mg single arms: mean AUClast sim/obs, 2020 model"),
    n16 = xo$n16, f16 = xo$f16, n20 = xo$n20, f20 = xo$f20,
    l300 = lv("c300", "lit_mean_ratio", "Clot 2021 300 mg published mean AUClast/AUCinf"), n300 = two("c300", "nca_mean_ratio", "Clot 2021 300 mg simulated NCA mean ratio,"),
    t300 = two("c300", "true_mean_ratio", "Clot 2021 300 mg simulated true mean ratio,"),
    lpk = lv("pkm", "lit_mean_ratio", "PKM12350 control arm published mean AUClast/AUCinf"), npk = two("pkm", "nca_mean_ratio", "PKM12350 simulated NCA mean ratio,"),
    l600 = lv("c600", "lit_mean_ratio", "Clot 2021 600 mg published mean AUClast/AUCinf"), n600 = two("c600", "nca_mean_ratio", "Clot 2021 600 mg simulated NCA mean ratio,"),
    t600 = two("c600", "true_mean_ratio", "Clot 2021 600 mg simulated true mean ratio,"), lt = ltol$p))))
  deck_end()
}
