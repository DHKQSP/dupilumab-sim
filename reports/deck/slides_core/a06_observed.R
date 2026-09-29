# 별첨 A6 실측 자료(300 mg): 저장소에 있는 공개 실측값만 쓴다.
#  왼쪽 그림: 공개 비구획 평균 AUClast/AUCinf(평균 AUClast ÷ 평균 AUCinf) 대 같은 조건 모의의 비구획 평균비와 참 평균비(두 모델).
#    PKM12350 대조군·시험군(FDA BLA 761055 임상약리 리뷰 Table 4.2.c), Clot 2021 Table 3 300 mg 코호트. 자료: results/literature/literature_numeric.csv
#    (시험군 값은 같은 파일 note 열 '시험군 488 대 554(88.1%)'에서 정규식으로 읽는다; AUCinf 산출 대상자 수가 다를 수 있어 일치 판정 밖).
#    일치 기준(±n%p)은 결과 파일 판정 문구에서 읽고, 판정을 수치로 다시 계산해 결과 파일 판정과 같은지 검사한다(결과보고 덱 A1과 같은 방식).
#  오른쪽 그림: 300 mg arm별 관측 평균 AUClast(PKM12350·PKM14161 시험약/대조약) 대 모의 평균(두 모델)과 검증 허용 범위(관측 ±15%).
#    자료: results/step1/step1f_300mg_arm_check.csv, results/step1_k2020/step1f_300mg_arm_check.csv(출처 표기는 config/design_clot2021.yaml arm_checks_300mg.source).
#    PKM14161은 2020 모델 개발 자료(config dev_k2020: internal)라 †로 표시한다.
#  Table 4.2.c의 GMR 행·표준편차는 SPEC.md 10절(서지) 문장에만 있고 결과 파일에 없어 싣지 않는다(캡션). FDA 리뷰 Table 4.2.b, TYENNE BMER(BLA 761275),
#  MSB11456 피하 초록은 저장소에 없어 뺐다(캡션, 별첨 A12).
slide_A6 <- function() {
  LN <- "literature/literature_numeric.csv"; A16 <- "step1/step1f_300mg_arm_check.csv"; A20 <- "step1_k2020/step1f_300mg_arm_check.csv"; DC <- "config/design_clot2021.yaml"
  deck_slide("A6", tag = "litsim")
  L <- DK$txt$A6
  W <- list(pkm = "grepl('PKM12350', source)", clot = "source=='Clot 2021 Table 3' & dose_mg==300")
  premise(nrow(rows(LN, W$pkm)) == 1 && nrow(rows(LN, W$clot)) == 1, "one PKM12350 control-arm row and one Clot 2021 300 mg row")
  premise(all(vapply(W, function(w) row1(LN, w)$dose_mg == 300, TRUE)) && as.numeric(.read("config/trial_design.yaml")$dose_mg) == 300, "both literature rows are at the study dose")
  # 일치 기준(±n%p): 판정 문구에서 읽는다
  tol <- local({ s <- row1(LN, W$clot)$judgment; m <- regmatches(s, regexec("±([0-9]+)%p", s))[[1]]; premise(length(m) == 2, "tolerance in the judgment label")
    x <- as.numeric(m[2]); list(x = x, p = dderived("agreement tolerance of the NCA mean ratio (points), from the judgment label", LN, sprintf("%s :: judgment, regex +/-([0-9]+)%%p", W$clot), x, fnum(x, 0))) })
  for (k in names(W)) { r <- row1(LN, W[[k]]); d <- 100 * (c(r$nca_mean_ratio_k2016, r$nca_mean_ratio_k2020) - r$lit_mean_ratio)
    premise(all(abs(d) <= tol$x + 1e-12) && startsWith(r$judgment, "일치"), sprintf("%s: simulated NCA mean ratio within the tolerance in both models, as in the stored judgment", k))
    premise(all(c(r$true_mean_ratio_k2016, r$true_mean_ratio_k2020) > pmax(r$lit_mean_ratio, r$nca_mean_ratio_k2016, r$nca_mean_ratio_k2020)),
            sprintf("%s: simulated true mean ratio above the published and simulated NCA ratios (caption: NCA ratio usually lower)", k)) }
  tr <- { x <- 100 * unlist(lapply(W, function(w) { r <- row1(LN, w); c(r$true_mean_ratio_k2016, r$true_mean_ratio_k2020) }))
    dderived("simulated true mean AUClast/AUCinf ratio, both 300 mg literature conditions, two models, rounded", LN, "PKM12350 and Clot 300 mg rows :: range(true_mean_ratio_k2016, true_mean_ratio_k2020) x 100",
             range(x), rng_fmt(min(x), max(x), 0, "%")) }
  dose <- dint(LN, W$clot, "dose_mg", "dose of the literature conditions (mg)")
  y0 <- core_title(tx("A6.title", list(tol = tol$p, tr = tr)), tx("A6.kicker", list(dose = dose)))

  # ---- PKM12350 시험군의 공개 평균비(결과 파일 note 열: '시험군 488 대 554(88.1%)') ----
  pnote <- row1(LN, W$pkm)$note; mt <- regmatches(pnote, regexec("([0-9]+) 대 ([0-9]+)\\(([0-9.]+)%\\)", pnote))[[1]]
  premise(length(mt) == 4 && abs(100 * as.numeric(mt[2]) / as.numeric(mt[3]) - as.numeric(mt[4])) < 0.05, "PKM12350 test arm: published AUClast and AUCinf means and their ratio in the note agree")
  premise(grepl("대상자 수가 다를 가능성", pnote) && grepl("확인 불가", pnote), "note: the test-arm AUCinf mean may use a different number of subjects, not verifiable (figure label, body)")
  tv_ <- as.numeric(mt[4])
  pt <- dderived("PKM12350 test arm published NCA mean ratio AUClast/AUCinf (from the note)", LN, "grepl('PKM12350', source) :: note, pattern 'AUClast-mean (Korean vs) AUCinf-mean(ratio%)', third number", tv_, paste0(fnum(tv_, 1), "%"))
  r_pk <- row1(LN, W$pkm); premise(abs(tv_ - 100 * r_pk$lit_mean_ratio) > tol$x && all(abs(tv_ - 100 * c(r_pk$nca_mean_ratio_k2016, r_pk$nca_mean_ratio_k2020)) > tol$x), "test arm ratio outside the agreement tolerance (title restricted to the control arm and Clot)")

  # ---- 캡션·본문(아래) ----
  src <- dcfg("design_clot2021.yaml", c("arm_checks_300mg", "source"), "source of the 300 mg single-arm AUClast means", function(x) as.character(x))
  capy <- core_caption(tx("A6.caption"), GEO$BODY_BOTTOM, size = 14)
  gt <- dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "validation criterion: simulated mean AUClast within +-x% of observed (%)", num_fmt(0))
  tolv <- .read(DC)$gate$auclast_mean_tol_pct / 100
  a16 <- rows(A16, "study %in% c('PKM12350','PKM14161')"); a20 <- rows(A20, "study %in% c('PKM12350','PKM14161')")
  premise(nrow(a16) == 4 && nrow(a20) == 4, "four arms (PKM12350, PKM14161; test and reference) in both models")
  premise(all(abs(c(a16[study == "PKM14161", ratio], a20[study == "PKM14161", ratio]) - 1) <= tolv) &&
          all(c(a16[study == "PKM12350", ratio], a20[study == "PKM12350", ratio]) > 1), "PKM14161 arms within the tolerance; PKM12350 arms simulated above observed (body)")
  arm_cfg <- .read(DC)$arm_checks_300mg
  premise(all(vapply(arm_cfg$arms, function(z) if (z$study == "PKM14161") z$dev_k2020 == "internal" && z$dev_k2016 == "external" else TRUE, TRUE)) &&
          all(vapply(arm_cfg$arms, function(z) if (z$study == "PKM12350") z$dev_k2020 == "external" && z$dev_k2016 == "external" else TRUE, TRUE)),
          "PKM14161 is 2020-model development data, PKM12350 external to both models (dagger, body)")
  premise(identical(arm_cfg$weight$status, "assumption"), "Li 2020 arm weights not reported: assumed for both studies (body)")
  dsrc("development-data status of the 300 mg arms", DC)
  rr <- function(st) { w <- sprintf("study=='%s'", st)
    drange(A16, w, "ratio", 2, "", sprintf("simulated / observed mean AUClast, %s arms, 2016", st)); drange(A20, w, "ratio", 2, "", sprintf("simulated / observed mean AUClast, %s arms, 2020", st))
    x <- range(c(rows(A16, w)$ratio, rows(A20, w)$ratio))
    dderived(sprintf("simulated / observed mean AUClast, %s arms, two models", st), A16, sprintf("%s :: range(ratio) over the step1 and step1_k2020 files", w), x, rng_fmt(x[1], x[2], 2)) }
  body <- tx("A6.body", list(tol = gt, r1 = rr("PKM14161"), r2 = rr("PKM12350"),
                             wt = dcfg("design_clot2021.yaml", c("arm_checks_300mg", "weight", "mean"), "assumed mean weight of the single arms (kg)", num_fmt(0))))
  by <- core_body(body, capy - 0.06)

  # ---- 왼쪽 그림: 평균 AUClast/AUCinf ----
  FL <- L$fig; ML <- DK$txt$common$models_short
  d <- rbindlist(lapply(names(W), function(k) { r <- row1(LN, W[[k]])
    data.table(k = k, type = c("lit", "nca", "nca", "tv", "tv"), model = c(NA, "k2016", "k2020", "k2016", "k2020"),
               v = 100 * c(r$lit_mean_ratio, r$nca_mean_ratio_k2016, r$nca_mean_ratio_k2020, r$true_mean_ratio_k2016, r$true_mean_ratio_k2020)) }))
  d <- rbind(d, data.table(k = "pkm", type = "litt", model = NA, v = tv_))              # PKM12350 시험군 공개 값(일치 판정 밖)
  dsrc("figure a6: published and simulated mean ratios", LN)
  # 세로 위치: 자료마다 머리글 한 줄 + 문헌(PKM12350은 대조군·시험군)·모의 비구획·모의 참값(두 모델은 위아래로 어긋나게)
  TY <- list(pkm = c(lit = 1, litt = 2, nca = 3, tv = 4), clot = c(lit = 1, nca = 2, tv = 3)); blk <- c(pkm = 0, clot = 5.3); OFF <- 0.2
  d[, y := -(blk[k] + mapply(function(k_, t_) TY[[k_]][[t_]], k, type)) + ifelse(is.na(model), 0, ifelse(model == "k2016", OFF, -OFF))]
  d[, shp := ifelse(is.na(model), "lit", model)][, ctype := ifelse(type == "litt", "lit", type)]
  lab <- d[, .(v = if (type[1] == "tv") min(v) else max(v), y = -(blk[k[1]] + TY[[k[1]]][[type[1]]]), lab = paste0(paste(fnum(v, 1), collapse = " / "), "%")), by = .(k, type)]
  lab[type == "litt", lab := fill(FL$litt_lab, list(v = lab))]
  lab[, hj := ifelse(type == "tv", 1.12, -0.12)]                       # 참값은 점 왼쪽(오른쪽 끝 100% 근처라 잘리지 않게)
  hd <- data.table(k = names(W), y = -blk[names(W)], lab = vapply(names(W), function(k) fill(FL$coh[[k]], list(dose = dose)), ""))
  band <- d[type == "lit", .(k, lo = v - tol$x, hi = v + tol$x, y0 = -(blk[k] + vapply(k, function(k_) max(TY[[k_]]), 1) + 0.45), y1 = -(blk[k] + 0.55))]
  yt <- rbindlist(lapply(names(W), function(k) data.table(y = -(blk[[k]] + TY[[k]]), lab = unlist(FL$type[[k]][names(TY[[k]])]))))
  XL <- c(floor(min(c(band$lo, d$v))) - 0.6, 100.4)
  premise(max(lab[type %in% c("lit", "nca"), v]) < 99, "value labels right of the literature and NCA points stay inside the panel")
  p1 <- ggplot(d) +
    geom_rect(data = band, aes(xmin = lo, xmax = pmin(hi, 100), ymin = y0, ymax = y1), fill = PAL$tint_grey, inherit.aes = FALSE) +   # AUClast/AUCinf는 100%를 넘지 않으므로 100%에서 자른다
    geom_vline(xintercept = 100, colour = PAL$muted, linewidth = 0.4) +
    geom_point(aes(v, y, shape = shp, colour = ctype), size = 3.4) +
    geom_text(data = lab, aes(v, y, label = lab, hjust = hj), size = PT(14), family = FONT, colour = PAL$ink) +
    geom_text(data = hd, aes(XL[1], y, label = lab), hjust = 0, size = PT(15), family = FONT, fontface = "bold", colour = PAL$ink) +
    scale_shape_manual(values = c(lit = 15, k2016 = 16, k2020 = 17), labels = c(lit = FL$lit_shape, k2016 = ML$k2016, k2020 = ML$k2020), name = NULL,
                       breaks = c("lit", "k2016", "k2020")) +
    scale_colour_manual(values = c(lit = PAL$ink, nca = PAL$orange, tv = PAL$blue), labels = c(nca = FL$c_nca, tv = FL$c_tv), name = NULL, breaks = c("nca", "tv")) +
    guides(shape = guide_legend(order = 1, override.aes = list(colour = PAL$ink)), colour = guide_legend(order = 2, override.aes = list(shape = 15, size = 4.2))) +
    scale_x_continuous(breaks = seq(ceiling(XL[1] / 2) * 2, 100, 2), labels = function(v) paste0(v, "%")) +
    scale_y_continuous(breaks = yt$y, labels = yt$lab) +
    coord_cartesian(xlim = XL, ylim = c(min(yt$y) - 0.6, 0.45), expand = FALSE, clip = "off") +
    labs(x = FL$xlab, y = NULL, subtitle = fill(FL$band, list(tol = tol$p))) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box = "horizontal", legend.spacing.x = grid::unit(2, "pt"),
          legend.box.spacing = grid::unit(2, "pt"), legend.text = element_text(size = 14, colour = PAL$ink), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 14, colour = PAL$ink),
          plot.subtitle = element_text(size = 14, colour = PAL$ink2, margin = margin(0, 0, 2, 0)), plot.margin = margin(2, 6, 2, 2))
  FW <- 6.35
  deck_figure(p1, "a6_ratio_literature", c(GEO$ML, y0, FW, by - 0.08 - y0), src = LN)

  # ---- 오른쪽 그림: arm별 관측 평균 AUClast 대 모의 평균 ----
  ar <- rbind(a16[, .(study, arm, obs = auclast_obs, sim = auclast_sim, model = "k2016")], a20[, .(study, arm, obs = auclast_obs, sim = auclast_sim, model = "k2020")])
  premise(all(a16$auclast_obs == a20[match(paste(a16$study, a16$arm), paste(study, arm)), auclast_obs]), "same observed means in both model files")
  ar[, lab := sprintf("%s %s%s", study, FL$arm[[arm]], if (study == "PKM14161") FL$dev_mark else ""), by = .(study, arm)]
  ord <- rev(unique(ar[order(study, -xtfrm(arm))]$lab))
  ar[, lab := factor(lab, levels = ord)]
  ob <- unique(ar[, .(lab, obs)])
  dsrc("figure a6: 300 mg arm means", c(A16, A20))
  XR <- GEO$ML + FW + 0.3; WR2 <- GEO$W - GEO$MR - XR
  p2 <- ggplot() +
    geom_segment(data = ob, aes(x = obs * (1 - tolv), xend = obs * (1 + tolv), y = lab, yend = lab), colour = PAL$tint_grey, linewidth = 9) +
    geom_point(data = ob, aes(obs, lab), shape = 15, size = 3.6, colour = PAL$ink) +
    geom_text(data = ob, aes(obs, lab, label = fnum(obs, 0)), vjust = -1.25, size = PT(14), family = FONT, colour = PAL$ink) +
    geom_point(data = ar, aes(sim, lab, shape = model), colour = PAL$blue, size = 3.4) +
    scale_shape_manual(values = c(k2016 = 16, k2020 = 17), labels = c(k2016 = paste(ML$k2016, FL$sim), k2020 = paste(ML$k2020, FL$sim)), name = NULL) +
    scale_x_continuous(limits = c(400, 720), breaks = seq(400, 700, 100), expand = expansion(0)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.8))) +
    labs(x = FL$xlab2, y = NULL, subtitle = fill(FL$band2, list(tol = gt, src = src))) + theme_core(16) +
    theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 14, colour = PAL$ink), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 14, colour = PAL$ink),
          plot.subtitle = element_text(size = 14, colour = PAL$ink2, margin = margin(0, 0, 2, 0), lineheight = 1.0))
  premise(all(ar$obs * (1 - tolv) > 400 & ar$sim < 720 & ar$obs * (1 + tolv) < 720), "arm values inside the x range")
  deck_figure(p2, "a6_arm_auclast", c(XR, y0, WR2, by - 0.08 - y0), src = c(A16, A20))

  # ---- 노트 ----
  v <- function(k, col, it) dv(LN, W[[k]], col, 1, "%", it, scale = 100)
  dd <- function(k) { r <- row1(LN, W[[k]]); x <- 100 * (c(r$nca_mean_ratio_k2016, r$nca_mean_ratio_k2020) - r$lit_mean_ratio)
    dderived(sprintf("%s: simulated NCA minus published mean ratio (points), 2016 and 2020", k), LN, sprintf("%s :: (nca_mean_ratio_k2016, nca_mean_ratio_k2020 - lit_mean_ratio) x 100", W[[k]]), x,
             sprintf("%s / %s", fnum(x[1], 1), fnum(x[2], 1))) }
  W600 <- "source=='Clot 2021 Table 3' & dose_mg==600"; W200 <- "source=='Clot 2021 Table 3' & dose_mg==200"
  r600 <- row1(LN, W600); premise(all(c(r600$nca_mean_ratio_k2016, r600$nca_mean_ratio_k2020, r600$true_mean_ratio_k2016, r600$true_mean_ratio_k2020) < r600$lit_mean_ratio),
                                  "600 mg: simulated NCA and true mean ratios below the published ratio (notes)")
  am <- function(f_, st, ar_, it) dv(f_, sprintf("study=='%s' & arm=='%s'", st, ar_), "ratio", 2, "", it)
  deck_notes(tx("A6.notes", list(
    dose = dose, tol = tol$p,
    pl = v("pkm", "lit_mean_ratio", "PKM12350 control arm, published NCA mean ratio"), pn16 = v("pkm", "nca_mean_ratio_k2016", "PKM12350, simulated NCA mean ratio, 2016"),
    pn20 = v("pkm", "nca_mean_ratio_k2020", "PKM12350, simulated NCA mean ratio, 2020"), pt16 = v("pkm", "true_mean_ratio_k2016", "PKM12350, simulated true mean ratio, 2016"),
    pt20 = v("pkm", "true_mean_ratio_k2020", "PKM12350, simulated true mean ratio, 2020"), pd = dd("pkm"),
    cl = v("clot", "lit_mean_ratio", "Clot 300 mg, published NCA mean ratio"), cn16 = v("clot", "nca_mean_ratio_k2016", "Clot 300 mg, simulated NCA mean ratio, 2016"),
    cn20 = v("clot", "nca_mean_ratio_k2020", "Clot 300 mg, simulated NCA mean ratio, 2020"), ct16 = v("clot", "true_mean_ratio_k2016", "Clot 300 mg, simulated true mean ratio, 2016"),
    ct20 = v("clot", "true_mean_ratio_k2020", "Clot 300 mg, simulated true mean ratio, 2020"), cd = dd("clot"),
    pw = dv("step1/step1g_literature_coverage.csv", "grepl('PKM12350', source)", "weight_mean", 0, "", "PKM12350 simulated cohort mean weight (assumed, kg)"),
    cw = dv("step1/step1g_literature_coverage.csv", "grepl('Clot', source) & dose_mg==300", "weight_mean", 1, "", "Clot 300 mg simulated cohort mean weight (kg)"),
    ptest = dderived("PKM12350 test arm published AUClast and AUCinf means (from the note)", LN, "grepl('PKM12350', source) :: note, regex groups 1 and 2", as.numeric(mt[2:3]), sprintf("%s / %s", mt[2], mt[3])), pt = pt,
    d600 = dint(LN, W600, "dose_mg", "dose of the 600 mg literature row (mg)"), l600 = dv(LN, W600, "lit_mean_ratio", 1, "%", "Clot 600 mg, published NCA mean ratio", scale = 100),
    t600 = { x <- 100 * c(r600$true_mean_ratio_k2016, r600$true_mean_ratio_k2020)
      dderived("Clot 600 mg, simulated true mean ratio, two models", LN, sprintf("%s :: range(true_mean_ratio_k2016, true_mean_ratio_k2020) x 100", W600), range(x), rng_fmt(min(x), max(x), 1, "%")) },
    d200 = dint(LN, W200, "dose_mg", "dose of the 200 mg literature row (mg)"),
    o1 = dv(A16, "study=='PKM12350' & arm=='test'", "auclast_obs", 1, "", "PKM12350 test arm observed mean AUClast"),
    o2 = dv(A16, "study=='PKM12350' & arm=='reference'", "auclast_obs", 1, "", "PKM12350 reference arm observed mean AUClast"),
    o3 = dv(A16, "study=='PKM14161' & arm=='test'", "auclast_obs", 0, "", "PKM14161 test arm observed mean AUClast"),
    o4 = dv(A16, "study=='PKM14161' & arm=='reference'", "auclast_obs", 0, "", "PKM14161 reference arm observed mean AUClast"),
    s16a = dv(A16, "study=='PKM12350' & arm=='test'", "auclast_sim", 0, "", "simulated mean AUClast, PKM12350 schedule, 2016"),
    s16b = dv(A16, "study=='PKM14161' & arm=='test'", "auclast_sim", 0, "", "simulated mean AUClast, PKM14161 schedule, 2016"),
    s20a = dv(A20, "study=='PKM12350' & arm=='test'", "auclast_sim", 0, "", "simulated mean AUClast, PKM12350 schedule, 2020"),
    s20b = dv(A20, "study=='PKM14161' & arm=='test'", "auclast_sim", 0, "", "simulated mean AUClast, PKM14161 schedule, 2020"),
    q1 = am(A16, "PKM12350", "test", "PKM12350 test sim/obs, 2016"), q2 = am(A16, "PKM12350", "reference", "PKM12350 reference sim/obs, 2016"),
    q3 = am(A20, "PKM12350", "test", "PKM12350 test sim/obs, 2020"), q4 = am(A20, "PKM12350", "reference", "PKM12350 reference sim/obs, 2020"),
    gt = gt)))
  deck_end()
}
