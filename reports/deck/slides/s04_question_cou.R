# S04 질문, 사용 맥락, 모델 위험(ICH M15 요소)과 이 분석이 답하는 범위. [문헌] 보고서 2절 표를 한국어로 줄이고, 근거 범위(시험 모집단)와 답하지 않는 것을 옆에 둔다.
# 질문 2의 구간(Day a~b)은 config의 추가 채혈 일정(D1, D4)이 채혈을 더하는 B0 구간에서 만든다(문구에 숫자를 쓰지 않음).
s04_qwin <- function() {
  sch <- .read("config/trial_design.yaml")$schedules; b0 <- unlist(sch$B0$days)
  br <- function(s_) { add <- setdiff(unlist(sch[[s_]]$days), b0); c(max(b0[b0 < min(add)]), min(b0[b0 > max(add)])) }
  w <- br("D1"); premise(identical(w, br("D4")), "D1 and D4 add samples within the same B0 interval (question 2 window)")
  sd <- w + 1
  list(lo = dderived("question 2 window start (study day): B0 sample before the samples added in D1 and D4", "config/trial_design.yaml",
                     "schedules.D1.days not in B0 :: max(B0 days below the added days) + 1", sd[1], fnum(sd[1], 0)),
       hi = dderived("question 2 window end (study day): B0 sample after the samples added in D1 and D4", "config/trial_design.yaml",
                     "schedules.D1.days not in B0 :: min(B0 days above the added days) + 1", sd[2], fnum(sd[2], 0)))
}
# Km 0.01~100배에서 참 AUC0-inf 비 범위(보고서 5.6절 km_rng와 같은 계산: 범위 끝 행과 도달 행)
s04_km <- function() {
  INV <- "oc/inversion_all.csv"; r <- rows(INV, "mechanism=='Km' & is.finite(end_auc_ratio)"); r2 <- rows(INV, "mechanism=='Km' & reachable==TRUE")
  x <- range(c(r$end_auc_ratio, r2$auc_ratio)); lim <- unlist(.read("config/trial_design.yaml")$be$limits)
  premise(x[1] > lim[1] && x[2] < lim[2], "Km changes keep the true AUC0-inf ratio inside the equivalence limits (text: Km difference not visible in AUC)")
  dderived("true AUC0-inf ratio range over Km changes (range ends and reached rows)", INV, "mechanism=='Km' :: range of end_auc_ratio, auc_ratio", x, rng_fmt(x[1], x[2], 3))
}

# f_study_days()와 같은 계산. config의 채혈일이 정수·실수 혼합이라 yaml이 목록으로 읽으므로 unlist한다(공용 f_study_days는 이 경우 멈춤: 수정 요청)
s04_study_days <- function(schedule = "B0", which = c("all", "last", "n")) {
  which <- match.arg(which); d <- unlist(.read("config/trial_design.yaml")$schedules[[schedule]]$days); premise(length(d) > 0, paste("schedule", schedule))
  sd <- d + 1; p <- switch(which, all = paste(fnum(sd[sd == round(sd)], 0), collapse = ", "), last = fnum(max(sd), 0), n = as.character(length(d)))
  dderived(sprintf("schedule %s, %s (study day = days after dose + 1)", schedule, which), "config/trial_design.yaml", sprintf("schedules.%s.days :: %s", schedule, which), sd, p)
}

slide_S04 <- function() {
  deck_slide("S04", tag = "lit")
  deck_kicker(tx("S04.kicker")); deck_title(tx("S04.title"))
  q <- s04_qwin()
  f <- list(d_lo = q$lo, d_hi = q$hi, wt = f_wt_range(), dose = f_dose(), n = f_n_arm(), nr = f_n_rand(),
            ns = s04_study_days("B0", "n"), last = s04_study_days("B0", "last"), lloq = f_lloq())

  # 왼쪽: ICH M15 요소 표(보고서 2절)
  keys <- c("q1", "q2", "cou", "risk", "cred")
  df <- data.frame(a = vapply(keys, function(k) tx(sprintf("S04.table.rows.%s", k), f)[1], ""),
                   b = vapply(keys, function(k) tx(sprintf("S04.table.rows.%s", k), f)[2], ""), check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- tx("S04.table.head")
  # 표 너비: 오른쪽 카드의 가장 긴 줄(주분석 M1 ..., M0 함께 제시)이 한 줄에 들어가도록 오른쪽 열을 넓힌다
  tw <- 6.6
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, tw, GEO$BODY_BOTTOM - GEO$BODY_TOP), widths = c(1.55, 5.05), size = 14, align_num = FALSE)

  # 오른쪽: 근거 범위(시험 모집단) 카드와 답하지 않는 것(카드 높이는 마지막 줄 아래 여백이 위 여백과 같도록)
  xr <- GEO$ML + tw + 0.25; wr <- GEO$CW - tw - 0.25; hs <- 2.58
  deck_text(tx("S04.scope", f), c(xr, GEO$BODY_TOP, wr, hs), size = 16, bg = PAL$tint_blue, geom = "roundRect", label = "scope", gap_pt = 3)
  yn <- GEO$BODY_TOP + hs + 0.10
  deck_text(tx("S04.not.head"), c(xr, yn, wr, 0.44), size = 18, bold = TRUE, label = "not_head")
  deck_bullets(tx("S04.not.bullets"), c(xr, yn + 0.42, wr, GEO$BODY_BOTTOM - yn - 0.42), size = 16, gap_pt = 5)

  INV <- "oc/inversion_all.csv"
  ada <- function(k, item, fmt) dcfg("trial_design.yaml", c("ada_sensitivity", k), item, fmt)
  deck_notes(tx("S04.notes", c(f, list(ada_f = ada("fraction", "ADA-like subgroup: share of subjects", function(x) paste0(fnum(100 * x, 0), "%")),
                                        ada_d = ada("onset_day", "ADA-like subgroup: onset (days)", num_fmt(0)),
                                        ada_k = ada("ke_multiplier", "ADA-like subgroup: ke multiplier", num_fmt(0)),
                                        km = s04_km(), km_lo = dext(INV, "mechanism=='Km'", "end_multiplier", min, 2, "", "smallest Km multiplier examined"),
                                        km_hi = dext(INV, "mechanism=='Km'", "end_multiplier", max, 0, "", "largest Km multiplier examined")))))
  deck_end()
}
