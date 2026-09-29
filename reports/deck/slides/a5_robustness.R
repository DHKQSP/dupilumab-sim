# A5 부록: 견고성: 체중 범위와 아토피 모집단(보고서 부록 I 요약). 평가변수 선택의 근거가 아니다(근거는 시험 모집단뿐). [모의]
# 덱에서 results/atopic/, results/weight_generalization/ 파일과 "아토피"라는 말을 쓰는 유일한 슬라이드다(check_deck 4: atopic_allowed_slides A5).
# 수치는 보고서 부록 I의 사실 청크(MS_report.Rmd 189~211행, 1360~1370행, 734~738행)와 같은 파일·행 조건·열로 읽는다(검사 8 대조).
# 그림: 체중 구간(40~150 kg, 구간 안 균일 분포, 구간당 20,000명; 구간 폭은 같지 않다)의 세트 (i) 신뢰 비율(reliable_rsq_extrap_pct), 2016 모델(model a)과 2020 모델(model b).
# weight_bands_B0_abcd.csv의 dev_range_note 열은 한글이다. 행 조건에는 영문 코드(model, band)만 쓴다.

slide_A5 <- function() {
  WB <- "weight_generalization/weight_bands_B0_abcd.csv"; ABN <- "atopic/atopic_criteria_by_band.csv"; ACF <- "atopic/atopic_coverage_failing.csv"
  AW <- "atopic/atopic_weight_table.csv"; AB <- "atopic/atopic_trials_bias.csv"; ACC <- "atopic/atopic_trials_concordance.csv"; AAF <- "atopic/atopic_trials_arm_failure.csv"
  deck_slide("A5", tag = "sim")
  L <- DK$txt$A5

  # ---- 전제 ----
  wb <- rows(WB, "model %in% c('a','b')")
  premise(all(wb[model == "a", model_variant] == "base") && all(wb[model == "b", model_variant] == "struct2020"), "weight bands: model a = 2016 model (base), model b = 2020 Model 1")
  lo_of <- function(b) as.numeric(sub("-.*", "", b))
  for (m in c("a", "b")) { x <- wb[model == m][order(lo_of(band))]
    premise(nrow(x) == 6 && all(diff(x$reliable_rsq_extrap_pct) < 0), sprintf("reliability under set (i) falls with every heavier band (model %s)", m)) }
  bands <- c("below 60", "60-90", "above 90-100", "above 100")
  cb <- rows(ABN, "distribution=='primary' & set %in% c('i','iii')")
  for (vr in unique(cb$variant)) for (s_ in c("i", "iii")) { y <- cb[variant == vr & set == s_][match(bands, band), fail_pct]
    premise(length(y) == 4 && !anyNA(y) && all(diff(y) > 0), sprintf("atopic: failing share rises across weight bands (%s, set %s)", vr, s_)) }
  premise(length(unique(cb$variant)) == 3, "three model variants")

  # ---- 제목·kicker ----
  acf <- rows(ACF, "distribution=='primary' & set=='iii'"); premise(nrow(acf) == 3, "three model variants in the coverage file")
  z <- paste0(dderived("smallest coverage of the true AUC0-inf by the sampling window in patients failing criteria set (iii), three model variants, rounded down to 0.1 (%)", ACF,
                       "distribution=='primary' & set=='iii' :: floor(1000 * min(window_fail_min)) / 10", 100 * min(acf$window_fail_min), fnum(floor(1000 * min(acf$window_fail_min)) / 10, 1)), "%")
  deck_kicker(tx("A5.kicker")); deck_title(tx("A5.title", list(z = z)), box = c(GEO$ML, GEO$TITLE_TOP, 8.6, GEO$TITLE_H))   # 두 절이 한 줄씩(쉼표 뒤에서 줄바꿈), 오른쪽 위 태그와 떨어지게

  # ---- 근거 아님 표시 ----
  wt <- f_wt_range()
  CH <- 0.84   # 두 줄 문구: 아래 여백이 위 여백과 비슷하게(검사 9 tight_bottom)
  deck_text(tx("A5.not_evidence", list(wt = wt)), c(GEO$ML, GEO$BODY_TOP, GEO$CW, CH), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "text_not_evidence", gap_pt = 0)

  # ---- 왼쪽: 체중 구간 그림 + 요점 ----
  Y0 <- GEO$BODY_TOP + CH + 0.12; XL <- GEO$ML; WL <- 6.9
  d <- copy(wb)[, .(model, band, y = reliable_rsq_extrap_pct)]
  d[, pk := c(a = "k2016", b = "k2020")[model]]; d[, mlab := factor(model_lab()[pk], levels = model_lab())]
  ord <- unique(d[order(lo_of(band)), band]); d[, x := match(band, ord)]
  xl <- gsub("-", "~", ord)
  st <- match(c("60-75", "75-90"), ord); premise(!anyNA(st) && diff(st) == 1, "study-range bands 60-75 and 75-90 are adjacent")
  lab <- d[x %in% c(1, length(ord))]
  F <- L$fig
  # 모델 개발 자료 범위 밖일 수 있는 구간(dev_range_note가 있는 구간): 점선과 표시 문구
  fl <- wb[model == "a" & dev_range_note != ""]; premise(nrow(fl) == 1 && lo_of(fl$band) == 130 && fl$band == ord[length(ord)], "only the heaviest band (from 130 kg) is flagged as possibly outside the model development data")
  premise(all(wb[model == "b" & dev_range_note != "", band] == fl$band), "the same band is flagged for the 2020 model")
  k130 <- dderived("lower bound of the weight band flagged as possibly outside the model development data (kg)", WB, "model=='a' & dev_range_note != '' :: lower bound of band", lo_of(fl$band), fnum(lo_of(fl$band), 0))
  xf <- match(fl$band, ord)
  p <- ggplot(d, aes(x = x, y = y, colour = mlab, shape = mlab, linetype = mlab)) +
    annotate("rect", xmin = min(st) - 0.5, xmax = max(st) + 0.5, ymin = -Inf, ymax = Inf, fill = PAL$tint_grey) +
    annotate("text", x = mean(st), y = 64, label = fill(F$study, list(wt = wt)), family = FONT, size = 12 / .pt, colour = PAL$ink2, vjust = 0) +
    annotate("segment", x = xf - 0.5, xend = xf - 0.5, y = 60, yend = 102, colour = PAL$muted, linewidth = 0.5, linetype = "dotted") +
    annotate("text", x = xf - 0.42, y = 101, label = fill(F$extrap, list(k = k130)), family = FONT, size = 12 / .pt, colour = PAL$ink2, hjust = 0, vjust = 1, lineheight = 0.95) +   # 눈금 글자(12 pt)와 같은 크기
    geom_line(linewidth = 0.9) + geom_point(size = 3) +
    geom_text(data = lab, aes(label = fnum(y, 1)), vjust = ifelse(lab$pk == "k2016", 1.9, -1.0), size = 4.3, family = FONT, show.legend = FALSE) +
    scale_colour_manual(values = unname(MODEL_COL)) + scale_shape_manual(values = unname(MODEL_SHAPE)) + scale_linetype_manual(values = unname(MODEL_LT)) +
    scale_x_continuous(breaks = seq_along(ord), labels = xl, expand = expansion(add = c(0.35, 0.75))) +   # 오른쪽: 외삽 표시 문구 자리
    scale_y_continuous(limits = c(60, 102), breaks = seq(60, 90, 10), expand = expansion(mult = c(0, 0))) +   # 100 눈금선 없음: 첫 구간 값 표시(2020 모델)가 눈금선에 걸리지 않게
    labs(x = F$xlab, y = NULL, subtitle = F$ylab) + theme_deck(13) +
    theme(legend.position = "top", legend.justification = "left", legend.key.width = grid::unit(2.2, "lines"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"), panel.grid.major.x = element_blank(), plot.subtitle = element_text(colour = PAL$ink2, size = 13, margin = margin(0, 0, 2, 0)))
  FH <- 2.49
  deck_figure(p, "a5_reliability_by_weight_band", c(XL, Y0, WL, FH), src = WB)
  # 그림 값의 추적(끝 구간 값)
  hi <- dv(WB, sprintf("model=='a' & band=='%s'", ord[1]), "reliable_rsq_extrap_pct", 1, "%", "weight bands, 2016 model, reliable under set (i), lightest band")
  lo <- dv(WB, sprintf("model=='a' & band=='%s'", ord[length(ord)]), "reliable_rsq_extrap_pct", 1, "%", "weight bands, 2016 model, reliable under set (i), heaviest band")
  hi20 <- dv(WB, sprintf("model=='b' & band=='%s'", ord[1]), "reliable_rsq_extrap_pct", 1, "%", "weight bands, 2020 Model 1, reliable under set (i), lightest band")
  lo20 <- dv(WB, sprintf("model=='b' & band=='%s'", ord[length(ord)]), "reliable_rsq_extrap_pct", 1, "%", "weight bands, 2020 Model 1, reliable under set (i), heaviest band")
  wa <- rows(WB, "model=='a'")
  wmed <- dderived("true extrapolation median range across weight bands", WB, "model=='a' :: range(extrap_true_median)", range(wa$extrap_true_median),
                   rng_fmt(min(wa$extrap_true_median), max(wa$extrap_true_median), 2, "%"))
  wcov <- dderived("largest coverage below 80% across weight bands (%)", WB, "model=='a' :: max(coverage_lt80_pct)", max(wa$coverage_lt80_pct), paste0(fnum(max(wa$coverage_lt80_pct), 2), "%"))
  thr <- dderived("window coverage threshold in column name coverage_lt80_pct (percent)", WB, "column name coverage_lt80_pct", 80, "80")
  brng <- dderived("uniform weight bands, lowest and highest bound (kg)", WB, "model=='a' :: min(lower bound of band), max(upper bound of band)", c(min(lo_of(wa$band)), max(as.numeric(sub(".*-", "", wa$band)))),
                   rng_fmt(min(lo_of(wa$band)), max(as.numeric(sub(".*-", "", wa$band))), 0))
  nb <- dint(WB, sprintf("model=='a' & band=='%s'", ord[1]), "n", "subjects per weight band")
  premise(all(wa$n == wa$n[1]), "same number of subjects in every weight band")
  od <- rbindlist(lapply(c("a", "b", "d"), function(m) rows(sprintf("weight_generalization/obese_trials_dropout_%s.csv", m))))
  dsrc("heavier dropouts in trials of the heavier population", sprintf("weight_generalization/obese_trials_dropout_%s.csv", c("b", "d")), "(table)")
  dwr <- range(od$wt_dropout_mean - od$wt_retained_mean); premise(all(dwr > 0), "failing subjects heavier than retained in every scenario and arm")
  dw <- dderived("dropouts heavier than retained subjects (kg)", "weight_generalization/obese_trials_dropout_a.csv", "range of wt_dropout_mean - wt_retained_mean over models a, b, d", dwr,
                 rng_fmt(dwr[1], dwr[2], 1))
  YB <- Y0 + FH + 0.08
  # 탈락자 = 엔진의 기본 신뢰 표시(reliable)를 못 넘은 대상자(R/summarize.R summarize_dropout). 엔진 표시는 세트 (ii)와 같다(config/nca_rules.yaml)
  nr_ <- .read("config/nca_rules.yaml")$standard$reliability; cs_ <- .read("config/prereg_20260926.yaml")$section4$criteria_sets$ii
  premise(nr_$adj_r2_min == cs_$adj_r2_min && nr_$extrap_max_pct == cs_$extrap_max_pct && nr_$span_ratio_min == cs_$span_ratio_min, "engine reliability flag used for the dropout comparison equals criteria set (ii)")
  premise(setequal(unique(od$model), c("a", "b", "d")), "dropout comparison covers the three model variants (a, b, d)")
  premise(length(unique(od$scenario)) > 1 && setequal(unique(od$arm), c("R", "T")), "dropout weight range spans several scenarios and both arms (bullet: model variants, scenarios, arms)")
  deck_bullets(tx("A5.bullets", list(brng = brng, hi = hi, lo = lo, wmed = wmed, thr = thr, wcov = wcov, dw = dw)),
               box = c(XL, YB, WL, GEO$BODY_BOTTOM - YB), size = 16, gap_pt = 4)

  # ---- 오른쪽: 아토피 피부염 성인 체중 분포(자리표시자) 카드 세 개 ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  deck_text(tx("A5.atopic_label"), c(XR, Y0, WR, 0.4), size = 16, bold = TRUE, color = PAL$ink2, label = "label_atopic", gap_pt = 0)
  at <- function(s_) dv(ABN, sprintf("variant=='base' & distribution=='primary' & band=='all' & set=='%s'", s_), "fail_pct", 1, "%", sprintf("atopic population, 2016 model, share failing criteria set (%s)", s_))
  atr <- function(s_) drange(ABN, sprintf("distribution=='primary' & band=='all' & set=='%s'", s_), "fail_pct", 1, "%", sprintf("atopic, three variants, set %s, band all, fail_pct", s_))
  AWS <- "distribution=='primary' & source=='simulated (20,000 subjects)'"
  aw <- function(col, it) dv(AW, AWS, col, 1, "%", sprintf("primary weight distribution, simulated, %s", it))
  wtr <- as.numeric(unlist(.read("config/trial_design.yaml")$weight$inclusion_kg))
  k60 <- dderived("weight threshold in column name below_60 (kg)", AW, "column name below_60", 60, "60"); k90 <- dderived("weight threshold in column name above_90 (kg)", AW, "column name above_90", 90, "90")
  premise(all(c("below_60", "above_90", "outside_60_90") %in% names(rows(AW))) && identical(wtr, c(60, 90)), "column thresholds below_60 / above_90 are the study weight range")
  pa <- .read("config/population_atopic.yaml")$distributions$primary
  premise(identical(pa$dist, "lognormal") && identical(.read("config/population_atopic.yaml")$status, "assumption"), "atopic weight distribution is a lognormal placeholder (status assumption)")
  cards <- list(
    list(v = sprintf("%s / %s", at("i"), at("iii")), l = tx("A5.card_fail", list(ri = atr("i"), riii = atr("iii"))), col = PAL$blue, bg = PAL$tint_blue),   # 2016 모델 값: 그림의 2016 모델 색
    list(v = paste0("≥ ", z), l = tx("A5.card_cov", list(med = drange(ACF, "distribution=='primary' & set=='iii'", "window_fail_median", 1, "%", "patients failing set iii: window_fail_median x 100, three variants", scale = 100))),
         col = PAL$ink, bg = PAL$tint_grey),   # 세 모델 변형의 최솟값: 모델 색 없음
    list(v = aw("outside_60_90", "share outside 60 to 90 kg"),
         l = tx("A5.card_wt", list(wt = wt, k60 = k60, k90 = k90, b60 = aw("below_60", "share below 60 kg"), a90 = aw("above_90", "share above 90 kg"),
                                   m = dcfg("population_atopic.yaml", c("distributions", "primary", "mean"), "atopic placeholder distribution, mean (kg)", num_fmt(0)),
                                   s = dcfg("population_atopic.yaml", c("distributions", "primary", "sd"), "atopic placeholder distribution, SD (kg)", num_fmt(0)))),
         col = PAL$ink2, bg = PAL$tint_grey))
  cy <- Y0 + 0.4; gap <- 0.07; chh <- (GEO$BODY_BOTTOM - cy - 2 * gap) / 3
  for (k in seq_along(cards)) deck_stat(cards[[k]]$v, cards[[k]]$l, c(XR, cy + (k - 1) * (chh + gap), WR, chh), color = cards[[k]]$col, bg = cards[[k]]$bg, value_size = 26)

  # ---- 노트 ----
  p130 <- local({ x <- pa; sl <- sqrt(log(1 + (x$sd / x$mean)^2)); ml <- log(x$mean) - sl^2 / 2; tr <- as.numeric(unlist(x$trunc))
    kb <- lo_of(fl$band); pp <- 100 * (plnorm(tr[2], ml, sl) - plnorm(kb, ml, sl)) / (plnorm(tr[2], ml, sl) - plnorm(tr[1], ml, sl))
    dderived(sprintf("share of the primary distribution above %s kg (analytic, truncated)", kb), "config/population_atopic.yaml", sprintf("distributions.primary :: lognormal share above %s kg within the truncation", kb), pp, fnum(pp, 1)) })
  bias <- function(ep) drange(AB, sprintf("analysis_model=='M0' & endpoint=='%s'", ep), "bias_pct", 2, "%", sprintf("atopic trials, bias of %s GMR, M0", ep))
  deck_notes(tx("A5.notes", list(
    wt = wt, hi20 = hi20, lo20 = lo20, p130 = p130, k130 = k130, nb = nb, k60 = k60, thr = thr, wcov = wcov, dw = dw,
    tr = dcfg("population_atopic.yaml", c("distributions", "primary", "trunc"), "atopic placeholder distribution, truncation (kg)", function(x) rng_fmt(x[1], x[2], 0)),
    a100 = aw("above_100", "share above 100 kg"),
    k100 = dderived("weight threshold in column name above_100 (kg)", AW, "column name above_100", 100, "100"),
    sens_i = drange(ABN, "distribution!='primary' & band=='all' & set=='i'", "fail_pct", 1, "%", "sensitivity weight distributions, set i"),
    sens_iii = drange(ABN, "distribution!='primary' & band=='all' & set=='iii'", "fail_pct", 1, "%", "sensitivity weight distributions, set iii"),
    b60i = dv(ABN, "variant=='base' & distribution=='primary' & band=='below 60' & set=='i'", "fail_pct", 1, "%", "atopic, base, set i, band below 60"),
    b100i = dv(ABN, "variant=='base' & distribution=='primary' & band=='above 100' & set=='i'", "fail_pct", 1, "%", "atopic, base, set i, band above 100"),
    bl = bias("AUClast"), bai = bias("AUCinf_Ai"),
    ag = drange(ACC, "analysis_model=='M0' & config=='A_i'", "agree_pct", 1, "%", "atopic trials, P2 and G2 A_i decision agreement, M0"),
    ad = drange(AAF, "scenario=='VM125' & set=='i'", "diff_mean", 2, "", "atopic trials, test minus reference failing set i, VM125 (points)"),
    ntr = dcfg("prereg_20260926.yaml", c("section5", "trial", "trials"), "atopic trials per scenario (pre-registered)", function(x) fint(as.numeric(x))))))
  deck_end()
}
