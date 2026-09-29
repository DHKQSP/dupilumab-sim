# S8 ② 보강(지시 §2, 그림 4-2): 케이스별 창 포착률(AUClast/AUCinf, 참값) 분포. 행 = 기본 두 모델, 곡선 모양(Vmax, Km; 두 군 모두), 정량한계, 비례 잔차,
# 체중 층(results/core_deck/coverage_by_case.csv, scripts/63). 표시: 중앙값 점, 5~95백분위 막대, 최솟값 표식, EMA 80% 기준선, 기본 조건 최솟값 주석,
# 문헌 두 점(비구획 평균비: Clot 2021 300 mg, FDA 리뷰 Table 4.2.c PKM12350 대조군; literature_numeric.csv). 제목 = 사전 등록 7b 규칙(title_rule.csv).
# 케이스 이름의 수치(영문 label 열에서 정규식으로 읽고 추적 행을 남긴다)
s8_lab <- function(case, rx, item) { rel <- "core_deck/coverage_by_case.csv"; s <- row1(rel, sprintf("case=='%s'", case))$label; m <- regmatches(s, regexec(rx, s))[[1]]
  premise(length(m) == 2, sprintf("case label %s matches %s", case, rx)); dderived(item, rel, sprintf("case=='%s' :: label, regex '%s'", case, rx), as.numeric(m[2]), m[2]) }
slide_S8 <- function() {
  CV <- "core_deck/coverage_by_case.csv"; TR <- "core_deck/title_rule.csv"; LN <- "literature/literature_numeric.csv"; LC <- "config/literature_core_deck.yaml"
  deck_slide("S8", tag = "litsim")
  tr <- row1(TR, "TRUE"); premise(!isTRUE(tr$all_cases_min_ge_84), "title rule: not every case has a minimum of at least 84% (use the 'otherwise' wording)")
  rw <- rows(CV, "role=='row'"); premise(nrow(rw) == tr$n_cases, "eleven case rows")
  premise(abs(tr$max_pct_lt80 - max(rw$pct_lt80)) < 1e-12 && abs(tr$min_p05 - min(rw$p05)) < 1e-12, "title rule values equal the case table")
  f <- list(bmin = { x <- 100 * min(tr$base_min_k2016, tr$base_min_k2020)
              dderived("window coverage, base case minimum over both models, rounded down", TR, "TRUE :: floor(min(base_min_k2016, base_min_k2020) x 100)", x, paste0(fnum(fl(x, 0), 0), "%")) },
            lt80 = dderived("largest share below 80% over all cases, rounded up", TR, "TRUE :: ceiling(max_pct_lt80, 3 decimals)", tr$max_pct_lt80, paste0(fnum(cl(tr$max_pct_lt80, 3), 3), "%")),
            p05 = { x <- 100 * tr$min_p05; dderived("smallest 5th percentile over all cases, rounded down", TR, "TRUE :: floor(min_p05 x 100)", x, paste0(fnum(fl(x, 0), 0), "%")) },
            e80 = dcfg("literature_core_deck.yaml", c("ema_be_guideline_80pct", "coverage_min_pct"), "EMA BE guideline: minimum coverage of AUC0-t (%)", num_fmt(0)))
  y0 <- core_title(tx("S8.title", f), tx("S8.kicker"))
  L <- DK$txt$S8$fig

  # ---- 행 이름: 영문 label 열에서 수치를 정규식으로 읽어 한국어 틀에 넣는다(하드코딩 금지) ----
  num1 <- function(s, rx) { m <- regmatches(s, regexec(rx, s))[[1]]; premise(length(m) == 2, sprintf("label '%s' matches %s", s, rx)); m[2] }
  lab_of <- function(case, label) switch(case,
    k2016_base = L$rows$k2016_base, k2020_base = L$rows$k2020_base,
    vmax080 = , vmax125 = fill(L$rows$vmax, list(v = num1(label, "x([0-9.]+)"))), km05 = , km10 = fill(L$rows$km, list(v = num1(label, "x([0-9.]+)"))),
    lloq002 = , lloq05 = fill(L$rows$lloq, list(v = num1(label, "LLOQ ([0-9.]+)"))), resid12 = fill(L$rows$resid, list(v = num1(label, "residual ([0-9.]+)%"))),
    wt60_75 = fill(L$rows$wt_lo, list(v = num1(label, "60-([0-9.]+) kg"))), wt75_90 = fill(L$rows$wt_hi, list(v = num1(label, ">([0-9.]+)-90"))))
  grp_of <- c(k2016_base = "base", k2020_base = "base", vmax080 = "curve", vmax125 = "curve", km05 = "curve", km10 = "curve", lloq002 = "lloq", lloq05 = "lloq",
              resid12 = "resid", wt60_75 = "wt", wt75_90 = "wt")
  d <- copy(rw)[, .(case, label, median, p05, p95, min)][, lab := mapply(lab_of, case, label)][, grp := grp_of[case]]
  premise(!anyNA(d$grp), "every case has a group")
  lit <- rbind(data.table(case = "lit_clot", v = row1(LN, "grepl('^Clot 2021', source) & dose_mg==300")$lit_mean_ratio, lab = L$rows$lit_clot),
               data.table(case = "lit_pkm", v = row1(LN, "grepl('PKM12350', source)")$lit_mean_ratio, lab = L$rows$lit_pkm))[, grp := "lit"]
  G <- c("base", "curve", "lloq", "resid", "wt", "lit")
  d[, grp := factor(grp, levels = G, labels = unlist(L$groups[G]))]; lit[, grp := factor(grp, levels = G, labels = unlist(L$groups[G]))]
  ord <- c(rev(lit$lab), rev(d$lab)); d[, lab := factor(lab, levels = ord)]; lit[, lab := factor(lab, levels = ord)]
  e80 <- .read(LC)$ema_be_guideline_80pct$coverage_min_pct / 100
  xlo <- floor(min(c(d$min, e80)) * 100 - 5) / 100                     # 기준선 왼쪽에 'EMA 기준' 표지 자리
  bl <- d[case %in% c("k2016_base", "k2020_base")]
  p <- ggplot() +
    geom_vline(xintercept = e80, colour = PAL$ink2, linewidth = 0.8, linetype = "22") +
    geom_segment(data = d, aes(x = p05, xend = p95, y = lab, yend = lab), colour = PAL$blue, linewidth = 4.2, alpha = 0.45) +
    geom_point(data = d, aes(median, lab), shape = 21, fill = PAL$blue, colour = "white", size = 4.2, stroke = 0.7) +
    geom_point(data = d, aes(min, lab), shape = 124, colour = PAL$ink, size = 6) +
    geom_label(data = d, aes(min, lab, label = paste0(fnum(100 * min, 1), "%")), hjust = -0.18, size = PT(14), family = FONT, colour = PAL$ink2,
               fill = "white", label.size = 0, label.padding = grid::unit(0.06, "lines")) +
    geom_point(data = lit, aes(v, lab), shape = 22, fill = PAL$ink, colour = "white", size = 4.0) +
    geom_text(data = lit, aes(v, lab, label = paste0(fnum(100 * v, 1), "%")), hjust = 1.35, size = PT(14), family = FONT, colour = PAL$ink2) +
    geom_text(data = data.table(grp = factor(levels(d$grp)[1], levels = levels(d$grp)), x = e80), aes(x = x, y = Inf, label = fill(L$ema, list(v = fnum(100 * e80, 0)))),
              hjust = 1.05, vjust = 1.2, size = PT(14), family = FONT, colour = PAL$ink2) +
    facet_grid(rows = vars(grp), scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(limits = c(xlo, 1), breaks = seq(0.8, 1, 0.05), labels = function(v) paste0(round(100 * v), "%"), expand = expansion(add = c(0, 0.004))) +
    labs(x = L$xlab, y = NULL) + theme_core(16) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = 14, colour = PAL$ink2, face = "plain"), panel.spacing.y = grid::unit(3, "pt"),
          panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 14, colour = PAL$ink), panel.grid.minor = element_blank())
  body <- tx("S8.body", list(e80 = f$e80))
  cap <- tx("S8.caption", list())
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  deck_figure(p, "s8_case_ranges", c(GEO$ML, y0, GEO$CW, by - 0.08 - y0), src = c(CV, LN, LC))

  # 노트: 사전 등록 규칙과 판정, 케이스별 최솟값
  mn <- function(case, item) dv(CV, sprintf("case=='%s'", case), "min", 1, "%", item, scale = 100)
  deck_notes(tx("S8.notes", list(
    m16 = mn("k2016_base", "minimum, base 2016"), m20 = mn("k2020_base", "minimum, base 2020"), mlq = mn("lloq05", "minimum, LLOQ 0.5"), mv = mn("vmax125", "minimum, Vmax x1.25"),
    mk = mn("km05", "minimum, Km x0.5"), mr = mn("resid12", "minimum, proportional residual 12%"), mcb = mn("curve_base", "minimum, curve-shape base draw (note row)"),
    p05v = dv(CV, "case=='vmax125'", "p05", 1, "%", "5th percentile, Vmax x1.25", scale = 100), lt = dv(CV, "case=='lloq05'", "pct_lt80", 3, "%", "share below 80%, LLOQ 0.5"),
    nsub = dint(CV, "case=='lloq05'", "n", "subjects, LLOQ 0.5 case"), e80 = f$e80,
    nlt = { r <- row1(CV, "case=='lloq05'"); x <- r$pct_lt80 * r$n / 100; dderived("subjects below 80%, LLOQ 0.5 case", CV, "case=='lloq05' :: pct_lt80 x n / 100", x, fnum(x, 0)) },
    t84 = dderived("title rule threshold (%)", "config/prereg_20260929.yaml", "section7.title_rule.all_cases_min_ge_84 :: threshold", 84, "84"),
    lq = s8_lab("lloq05", "LLOQ ([0-9.]+)", "LLOQ of the LLOQ case (mg/L)"), vm = s8_lab("vmax125", "x([0-9.]+)", "Vmax multiplier"), km = s8_lab("km05", "x([0-9.]+)", "Km multiplier"),
    rs = s8_lab("resid12", "residual ([0-9.]+)%", "proportional residual (%)"),
    cd = dv(LN, "grepl('^Clot 2021', source) & dose_mg==300", "dose_mg", 0, "", "Clot 2021 dose (mg)"),
    e20 = dcfg("literature_core_deck.yaml", c("ema_be_guideline_80pct", "observations_share_pct"), "EMA BE guideline: share of observations below 80% that triggers discussion (%)", num_fmt(0)),
    clot = dv(LN, "grepl('^Clot 2021', source) & dose_mg==300", "lit_mean_ratio", 1, "%", "Clot 2021 300 mg published NCA mean ratio", scale = 100),
    clots = dv(LN, "grepl('^Clot 2021', source) & dose_mg==300", "true_mean_ratio_k2016", 1, "%", "simulated true mean ratio, Clot 300 mg cohort, 2016", scale = 100),
    pkm = dv(LN, "grepl('PKM12350', source)", "lit_mean_ratio", 1, "%", "PKM12350 control arm published NCA mean ratio", scale = 100),
    pkms = dv(LN, "grepl('PKM12350', source)", "true_mean_ratio_k2016", 1, "%", "simulated true mean ratio, PKM12350 cohort, 2016", scale = 100))))
  deck_end()
}
