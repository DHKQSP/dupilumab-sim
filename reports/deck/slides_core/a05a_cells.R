# A5a ③ 세부(별첨): 경계 16칸 표. S9 점도표의 칸별 값. 두 표(2016 모델, 2020 모델; 표마다 8행, 14pt): 조건(S9 fig.mech 짧은 이름 → 목표 참 AUCinf 비),
# AUClast + Cmax(P2, results/oc_models/type1_models.csv)와 AUCinf + Cmax(규칙 A 세트 (iii), results/criteria/criteria_g2_type1.csv G2_A_iii)의 경계 1종 오류(M1),
# Wilson 95% 분류 표시(▲ 초과: 하한 > 5%, ○ 명목: 구간이 5% 포함, 표시 없음: 보수적 = 상한 < 5%; 결과 파일 class 열, 전제로 lo/hi와 대조).
# 시험 수: 칸마다 10,000회, 연장 칸(2020 모델 말초 분포 증가)은 20,000회(†). 자료 논리는 결과보고 덱 s15_mechanism.R, s14_rules.R.
A5A_SC <- c("F_down_080", "F_up_125", "ke_up_080", "ke_down_125", "Vmax_up_080", "Vmax_down_125", "V2_up_080", "ka_down_080")   # S9 그림과 같은 위에서 아래 순서

slide_A5a <- function() {
  T1 <- "oc_models/type1_models.csv"; CG <- "criteria/criteria_g2_type1.csv"
  WP <- "analysis_model=='M1' & config=='P2'"; WG <- "analysis_model=='M1' & config=='G2_A_iii'"
  deck_slide("A5a", tag = "sim")
  L <- DK$txt$A5a; LS <- DK$txt$S9$fig                                   # 기전 이름·칸 이름 틀은 S9와 같은 문구(text/ko_core/s09.yaml fig)
  ML <- DK$txt$common$models

  # ---- 전제 ----
  a <- rows(T1, WP); b <- rows(CG, WG)
  premise(nrow(a) == 16 && nrow(b) == 16 && setequal(paste(a$pk_model, a$scenario), paste(b$pk_model, b$scenario)), "16 boundary cells in both configurations")
  for (pk in c("k2016", "k2020")) premise(setequal(a[pk_model == pk, scenario], A5A_SC) && setequal(b[pk_model == pk, scenario], A5A_SC), sprintf("eight boundary scenarios, %s", pk))
  nomv <- 100 * (1 - .read("config/trial_design.yaml")$be$ci_level) / 2
  wcls <- function(d) ifelse(d$hi < nomv, "conservative", ifelse(d$lo > nomv, "exceeding", "nominal"))
  premise(all(wcls(a) == a$class) && all(wcls(b) == b$class), "class column equals the Wilson classification against the nominal level (marks)")
  lab_chk <- unique(rbind(a, b, fill = TRUE)[, .(scenario, mechanism, direction, target)]); premise(nrow(lab_chk) == length(A5A_SC), "one condition label per scenario (same in both models and files)")
  base_n <- unique(a[!(pk_model == "k2020" & scenario == "V2_up_080"), n_trials]); ext <- a[n_trials != base_n[1]]
  premise(length(base_n) == 1 && nrow(ext) == 1 && ext$pk_model == "k2020" && ext$scenario == "V2_up_080" &&
            b[pk_model == "k2020" & scenario == "V2_up_080", n_trials] == ext$n_trials && length(unique(b[!(pk_model == "k2020" & scenario == "V2_up_080"), n_trials])) == 1,
          "every cell has the pre-specified trial count except the extended 2020 V2 cell, the same in both files (dagger)")
  premise(all(a[scenario == "ka_down_080", pass_pct] == 0) && all(b[scenario == "ka_down_080", pass_pct] == 0) &&
            all(rows(T1, "grepl('^ka_', scenario)")$cmax_ratio < min(unlist(.read("config/trial_design.yaml")$be$limits))), "ka-down cells: 0% in both configurations, true Cmax ratio outside the limits (caption)")
  pe_ <- a[class == "exceeding"]; premise(nrow(pe_) == 1 && pe_$pk_model == "k2020" && pe_$scenario == "V2_up_080" && b[pk_model == "k2020" & scenario == "V2_up_080", class] == "conservative",
                                         "the one exceeding AUClast + Cmax cell is the 2020 V2 up cell, where the AUCinf configuration is conservative (body, notes)")
  top <- b[, .SD[which.max(pass_pct)], by = pk_model]
  premise(nrow(top) == 2 && length(unique(top$scenario)) == 1, "AUCinf configuration: the highest cell is the same condition in both models (body)")

  # ---- 제목 ----
  f <- list(n = dcount(T1, WP, "boundary cells, M1 P2"), nom = f_nominal(),
            e = dcount(CG, paste(WG, "& class=='exceeding'"), "AUCinf (set iii) + Cmax cells classified exceeding (Wilson lower bound above 5%), M1"),
            c = dcount(T1, paste(WP, "& class=='conservative'"), "AUClast + Cmax cells classified conservative (Wilson upper bound below 5%), M1"))
  y0 <- core_title(tx("A5a.title", f), tx("A5a.kicker", f))

  # ---- 표 두 개(모델별 8행) ----
  mk <- function(pk) {
    rel_lab <- function(sc) { r <- a[pk_model == pk & scenario == sc]
      s <- fill(LS$cell, list(mech = LS$mech[[paste(r$mechanism, r$direction, sep = "_")]], tgt = dv(T1, sprintf("%s & pk_model=='%s' & scenario=='%s'", WP, pk, sc), "target", 2, "", "target true AUCinf ratio of the boundary cell")))
      if (r$n_trials != base_n) paste0(s, L$table$ext) else s }
    val <- function(rel, w, sc, item) { ww <- sprintf("%s & pk_model=='%s' & scenario=='%s'", w, pk, sc)
      paste0(dv(rel, ww, "pass_pct", 2, "%", item), L$table$mark[[row1(rel, ww)$class]]) }
    df <- data.frame(a = vapply(A5A_SC, rel_lab, ""),
                     b = vapply(A5A_SC, function(sc) val(T1, WP, sc, sprintf("boundary type I error, AUClast + Cmax, %s, M1", pk)), ""),
                     c = vapply(A5A_SC, function(sc) val(CG, WG, sc, sprintf("boundary type I error, AUCinf (set iii) + Cmax, %s, M1", pk)), ""), stringsAsFactors = FALSE)
    names(df) <- tx("A5a.table.head"); df
  }
  gap <- 0.30; tw <- (GEO$CW - gap) / 2
  cap <- tx("A5a.caption", list(nom = f$nom, reps = f_reps("boundary"),
                                ext = dint(T1, sprintf("%s & pk_model=='k2020' & scenario=='V2_up_080'", WP), "n_trials", "trials in the extended cell (2020 model, V2 up)"),
                                zero = dext(T1, paste(WP, "& grepl('^ka_', scenario)"), "pass_pct", max, 0, "%", "AUClast + Cmax in the ka-down cells, M1")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  mech_top <- LS$mech[[paste(top$mechanism[1], top$direction[1], sep = "_")]]
  pn_ <- b[pass_pct > nomv & class != "exceeding"]
  premise(nrow(pn_) == 1 && pn_$pk_model == "k2016" && pn_$scenario == "ke_up_080" && pn_$class == "nominal", "AUCinf configuration: the only cell above 5% (point) that is not exceeding is the 2016 linear-elimination-up cell, nominal (body)")
  bf <- list(nom = f$nom, mke = LS$mech[["ke_up"]], k = dcount(CG, paste(WG, "& pass_pct > 5"), "AUCinf (set iii) + Cmax cells above 5% (point), M1 (S9 count)"),
             ke = dv(CG, sprintf("%s & pk_model=='k2016' & scenario=='ke_up_080'", WG), "pass_pct", 2, "%", "AUCinf (set iii) + Cmax, 2016 linear elimination up cell, M1"))
  by <- core_body(tx("A5a.body", bf), capy - 0.06)
  hh <- 0.39; th <- by - 0.10 - y0 - hh
  for (k in 1:2) {
    pk <- c("k2016", "k2020")[k]; x <- GEO$ML + (k - 1) * (tw + gap)
    deck_text(ML[[pk]], c(x, y0, tw, hh), size = 16, bold = TRUE, label = sprintf("label_%s", pk), gap_pt = 0)
    hl <- if (pk == pe_$pk_model) match(pe_$scenario, A5A_SC) else NULL             # AUClast 구성의 초과 칸(별첨 A5b에서 분해)
    deck_table(mk(pk), box = c(x, y0 + hh, tw, th), widths = c(2.3, 1.2, 2.45), size = 14, label = sprintf("table_%s", pk), highlight = hl, highlight_fill = PAL$tint_blue)
  }
  dsrc("boundary cell tables", c(T1, CG), "(table)")

  # ---- 노트 ----
  cnt <- function(rel, w, cls, item) dcount(rel, sprintf("%s & class=='%s'", w, cls), item)
  deck_notes(tx("A5a.notes", c(f, bf[c("mke", "ke")], list(mech = mech_top, 
    pn = cnt(T1, WP, "nominal", "AUClast + Cmax cells nominal, M1"), pe = cnt(T1, WP, "exceeding", "AUClast + Cmax cells exceeding, M1"),
    gc = cnt(CG, WG, "conservative", "AUCinf (set iii) + Cmax cells conservative, M1"), gn = cnt(CG, WG, "nominal", "AUCinf (set iii) + Cmax cells nominal, M1"),
    gk = dcount(CG, paste(WG, "& pass_pct > 5"), "AUCinf (set iii) + Cmax cells above 5% (point), M1"),
    pk_ = dcount(T1, paste(WP, "& pass_pct > 5"), "AUClast + Cmax cells above 5% (point), M1"),
    i16 = dci(CG, sprintf("%s & pk_model=='k2016' & scenario=='%s'", WG, top[pk_model == "k2016", scenario]), "pass_pct", "lo", "hi", 2, "%", "largest cell with interval, AUCinf (set iii) + Cmax, 2016, M1"),
    i20 = dci(CG, sprintf("%s & pk_model=='k2020' & scenario=='%s'", WG, top[pk_model == "k2020", scenario]), "pass_pct", "lo", "hi", 2, "%", "largest cell with interval, AUCinf (set iii) + Cmax, 2020, M1"),
    p2ci = dci(T1, sprintf("%s & pk_model=='k2020' & scenario=='V2_up_080'", WP), "pass_pct", "lo", "hi", 2, "%", "AUClast + Cmax, 2020 V2 up cell, M1"),
    g2v = dci(CG, sprintf("%s & pk_model=='k2020' & scenario=='V2_up_080'", WG), "pass_pct", "lo", "hi", 2, "%", "AUCinf (set iii) + Cmax, 2020 V2 up cell, M1"),
    cmr = drange(T1, paste(WP, "& grepl('^ka_', scenario)"), "cmax_ratio", 3, "", "true Cmax ratio, ka-down cells"),
    npm = dcount(T1, paste(WP, "& pk_model=='k2016'"), "boundary cells per PK model"), r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), reps = f_reps("boundary"),
    m0g = dcount(CG, "analysis_model=='M0' & config=='G2_A_iii' & class=='exceeding'", "AUCinf (set iii) + Cmax cells exceeding, M0"),
    m0p = dcount(T1, "analysis_model=='M0' & config=='P2' & class=='conservative'", "AUClast + Cmax cells conservative, M0")))))
  deck_end()
}
