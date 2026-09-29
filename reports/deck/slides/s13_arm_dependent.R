# S13 논거 ① 처리군 의존 탈락: 시험군 - 대조군 기준 미달 비율 차이(시험 2,000회 평균, 95% 구간)를 시나리오별로(동일 제품 S00, 제품 시나리오,
# 경계 칸) 세트 (i)·(iii), 두 모델로 그린다. 동일 제품에서는 약 0, 제품이 다르면 기전에 따라 방향·크기가 다르게 벌어진다.
# 수치: results/trialpop/tp_arm_difference.csv(scripts/55), tp_arm_difference_fit.csv(서술용 선형 적합), tp_identity_check.csv(재생성 일치).
# 주의: 인원 수 기반 결과라 분석 모형(M0/M1)과 무관하다(M1로 표시하지 않는다). 세트 (iii)에서는 ke 시나리오도 CI가 0을 포함하지 않는다.
S13_PROD <- c("F097", "KE110", "KE120", "VM125")
s13_mult <- function(sc) {
  y <- .read("config/scenarios.yaml")$scenarios[[sc]]$T_multipliers; premise(length(y) == 1, sprintf("scenario %s has one multiplier", sc))
  dcfg("scenarios.yaml", c("scenarios", sc, "T_multipliers", names(y)), sprintf("scenario %s multiplier on %s", sc, names(y)), num_fmt(2))
}
slide_S13 <- function() {
  TAD <- "trialpop/tp_arm_difference.csv"; TAF <- "trialpop/tp_arm_difference_fit.csv"; TID <- "trialpop/tp_identity_check.csv"
  deck_slide("S13", tag = "sim")
  L <- DK$txt$S13
  a <- rows(TAD, "set %in% c('i','iii')")
  BND <- setdiff(unique(a$scenario), c("S00", S13_PROD))
  g <- function(sc, s_) a[scenario == sc & set == s_]

  # ---- 전제(부호·크기·같은 참 비) ----
  premise(all(abs(a[scenario == "S00", diff_mean]) < 0.5) && all(a[scenario == "S00", diff_lo < 0 & diff_hi > 0]), "identical products: arm difference near zero and 95% interval includes zero (sets i and iii, two models)")
  premise(all(g("VM125", "i")$diff_lo > 0) && all(g("ka_down_080", "i")$diff_lo > 0) && all(g("Vmax_up_080", "i")$diff_lo > 0), "Vmax x1.25, ka decreased and Vmax increased: test arm fails more (set i)")
  premise(all(a[scenario == "V2_up_080", diff_hi] < 0), "V2 increased: test arm fails less (sets i and iii)")
  premise(length(BND) == 8 && all(abs(a[grepl("_080$", scenario), auc_ratio] - 0.8) < 0.005) && all(abs(a[grepl("_125$", scenario), auc_ratio] - 1.25) < 0.005) &&
            all(grepl("_(080|125)$", BND)), "boundary cells: true AUC0-inf ratio 0.80 or 1.25")
  premise(all(g("KE110", "iii")$diff_hi < 0) && all(g("KE120", "iii")$diff_hi < 0), "ke scenarios under set iii: 95% interval excludes zero (negative)")
  ri <- a[set == "i" & scenario != "S00"]; k <- which.max(abs(ri$diff_mean))
  premise(ri$scenario[k] == "ka_down_080", "largest set i difference is the ka decreased boundary cell")
  f_ <- rows(TAF, "set %in% c('i','iii')")
  premise(f_[pk_model == "k2020" & set == "i", slope_lo < 0 & slope_hi > 0] && all(f_[set == "iii", slope_lo < 0 & slope_hi > 0]) && f_[pk_model == "k2016" & set == "i", slope_hi < 0],
          "fit slope CI includes zero for k2020 set i and both models set iii only")

  # ---- 제목 ----
  zero <- dderived("reference value: no between-arm difference in the share failing (points)", TAD, "reference line diff_mean = 0", 0, "0")
  mx <- dderived("largest absolute test minus reference difference in the share failing set i over the test scenarios of both models (points)", TAD,
                 "set=='i' & scenario!='S00' :: max(abs(diff_mean))", abs(ri$diff_mean[k]), fnum(abs(ri$diff_mean[k]), 1))
  deck_kicker(tx("S13.kicker")); deck_title(tx("S13.title", list(zero = zero, mx = mx)))

  # ---- 왼쪽: 시나리오별 arm 간 차이(세트 (i)·(iii), 두 모델) ----
  ad <- function(sc, s_, d = 2) drange(TAD, sprintf("scenario=='%s' & set=='%s'", sc, s_), "diff_mean", d, "", sprintf("test minus reference failing, %s, set %s, diff_mean", sc, s_))
  # 음의 차이는 크기로("시험군이 덜 탈락"): 전제에서 부호 확인
  adn <- function(sc, s_, d = 2) { premise(all(g(sc, s_)$diff_mean < 0), sprintf("%s set %s: test arm fails less", sc, s_))
    drange(TAD, sprintf("scenario=='%s' & set=='%s'", sc, s_), "diff_mean", d, "", sprintf("reference minus test failing (magnitude), %s, set %s, -diff_mean", sc, s_), scale = -1) }
  r80 <- drange(TAD, "grepl('_080$', scenario) & set=='i'", "auc_ratio", 2, "", "true AUC0-inf ratio, boundary cells at 0.80, two models")
  r125 <- drange(TAD, "grepl('_125$', scenario) & set=='i'", "auc_ratio", 2, "", "true AUC0-inf ratio, boundary cells at 1.25, two models")
  mult <- vapply(S13_PROD, s13_mult, "")
  G <- L$fig$groups
  grp_lab <- c(s00 = G$s00, prod = G$prod, b080 = fill(G$b080, list(r = r80)), b125 = fill(G$b125, list(r = r125)))
  dd <- copy(a)
  dd[, grp := fifelse(scenario == "S00", "s00", fifelse(scenario %in% S13_PROD, "prod", fifelse(grepl("_080$", scenario), "b080", "b125")))]
  # 행 순서: 무리마다 세트 (i) 두 모델 평균 차이가 큰 것부터(제품 시나리오는 참 비 큰 것부터)
  ordk <- dd[set == "i", .(m = mean(diff_mean), r = mean(auc_ratio)), by = .(scenario, grp)]
  ordk[, key := fifelse(grp == "prod", -r, -m)]
  ordk <- ordk[order(match(grp, names(grp_lab)), key)]; ordk[, row := rev(seq_len(.N))]
  MECH <- c(F = "F", ke = "ke", Vmax = "Vmax", ka = "ka", V2 = "V2")
  ordk[, lab := vapply(seq_len(.N), function(j) { sc <- scenario[j]; r_ <- a[scenario == sc & set == "i"]
    if (sc == "S00") L$fig$s00
    else if (sc %in% S13_PROD) fill(L$fig$prod, list(mech = MECH[[r_$mechanism[1]]], mult = mult[[sc]], r = fnum(mean(r_$auc_ratio), 2)))
    else fill(L$fig$bnd, list(mech = MECH[[r_$mechanism[1]]], dir = L$fig$dir[[r_$direction[1]]])) }, "")]
  dd <- ordk[, .(scenario, row, lab)][dd, on = "scenario"]
  dd[, y := row + fifelse(pk_model == "k2016", 0.17, -0.17)]
  dd[, grp := factor(grp_lab[grp], levels = grp_lab)]; dd[, model := factor(model_lab()[pk_model], levels = model_lab())]
  dd[, setl := factor(unlist(DK$txt$common$sets)[set], levels = unlist(DK$txt$common$sets[c("i", "iii")]))]
  labmap <- setNames(ordk$lab, ordk$row)
  ntr_n <- unique(a$n_trials); premise(length(ntr_n) == 1, "same number of trials in every scenario")
  p <- ggplot(dd, aes(x = diff_mean, y = y, colour = model, shape = model)) +
    geom_vline(xintercept = 0, colour = PAL$ink2, linewidth = 0.5) +
    geom_errorbarh(aes(xmin = diff_lo, xmax = diff_hi), height = 0, linewidth = 0.7) +
    geom_point(size = 2.6) +
    facet_grid(grp ~ setl, scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = unname(MODEL_SHAPE)) +
    scale_y_continuous(breaks = ordk$row, labels = function(b) unname(labmap[as.character(b)]), expand = expansion(add = 0.5)) +
    scale_x_continuous(breaks = seq(-10, 15, by = 5)) +
    labs(x = fill(L$fig$xlab, list(n = fint(ntr_n))), y = NULL) + theme_deck(13) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(2, "pt"),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12, face = "bold", colour = PAL$ink2, lineheight = 0.95),
          strip.text.x = element_text(size = 13, face = "bold"), panel.grid.major.y = element_blank(), panel.spacing.y = grid::unit(6, "pt"),
          panel.spacing.x = grid::unit(14, "pt"), panel.background = element_rect(fill = "#fafaf8", colour = NA), axis.text.y = element_text(size = 12, colour = PAL$ink))
  XL <- GEO$ML; WF <- 7.1
  deck_figure(p, "s13_arm_difference", c(XL, GEO$BODY_TOP, WF, GEO$BODY_BOTTOM - GEO$BODY_TOP), src = TAD)

  # ---- 오른쪽: 요점 ----
  XR <- XL + WF + 0.25; WR <- GEO$W - GEO$MR - XR
  deck_bullets(tx("S13.bullets", list(
    zero = zero, s00i = ad("S00", "i"), s00iii = ad("S00", "iii"), vmm = mult[["VM125"]], vm = ad("VM125", "i"), ka = ad("ka_down_080", "i"),
    r80 = r80, v2 = adn("V2_up_080", "i"), v2iii = adn("V2_up_080", "iii"),
    r2i = drange(TAF, "set=='i'", "r_squared", 2, "", "arm difference per 10% lower true ratio, set i, r_squared"),
    r2iii = drange(TAF, "set=='iii'", "r_squared", 2, "", "arm difference per 10% lower true ratio, set iii, r_squared"))),
    box = c(XR, GEO$BODY_TOP, WR, GEO$BODY_BOTTOM - GEO$BODY_TOP), size = 16, gap_pt = 9)

  # ---- 노트 ----
  idr <- rows(TID); premise(nrow(idr) == 2 && all(idr$ok), "trial-population regeneration equals section1 in both models")
  nid <- dderived("regenerated per-arm counts identical to section1 (rows compared, both models)", TID, "sum(rows_compared), all ok", sum(idr$rows_compared), fint(sum(idr$rows_compared)))
  deck_notes(tx("S13.notes", list(
    wt = f_wt_range(), ntr = dint(TAD, "pk_model=='k2016' & scenario=='S00' & set=='i'", "n_trials", "regenerated trials per scenario"), nid = nid,
    fr_i = drange(TAD, "scenario=='S00' & set=='i'", "fail_R_mean", 2, "%", "reference-arm mean share failing, set i, two models"),
    fr_iii = drange(TAD, "scenario=='S00' & set=='iii'", "fail_R_mean", 2, "%", "reference-arm mean share failing, set iii, two models"),
    f97 = mult[["F097"]], k110 = mult[["KE110"]], k120 = mult[["KE120"]], vm125 = mult[["VM125"]],
    pr = drange(TAD, "scenario %in% c('F097','KE110','KE120','VM125') & set=='i'", "auc_ratio", 2, "", "true AUC0-inf ratio, product scenarios, two models"),
    r80 = r80, r125 = r125, nb = dcount(TAD, "pk_model=='k2016' & set=='i' & !scenario %in% c('S00','F097','KE110','KE120','VM125')", "boundary cells per model"),
    p_s00 = dspan(TAD, "scenario=='S00' & set=='i'", "diff_p05", "diff_p95", 2, "", "per-trial arm difference, identical products, set i, 5th to 95th percentile, two models"),
    zero = zero, ke120 = ad("KE120", "i"),
    ke_iii = drange(TAD, "scenario %in% c('KE110','KE120') & set=='iii'", "diff_mean", 2, "", "reference minus test failing (magnitude), ke scenarios, set iii, two models", scale = -1),
    nsc = dint(TAF, "pk_model=='k2016' & set=='i'", "n_scenarios", "scenario means per linear fit"))))
  deck_end()
}
