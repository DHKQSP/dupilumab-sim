# S03 배경: 규제 기대(단회 투여 PK 시험의 1차 평가변수 AUC0-inf)와 비교 유효성 시험 면제 흐름, 두필루맙에서의 문제. [문헌]
# 기관·연도 카드 세 개(EMA 2012, FDA 2016, FDA 2025 초안) + 아래 두필루맙 문제 요점. 수치는 설계 상수(동등성 한계, CI 수준)뿐.

slide_S03 <- function() {
  deck_slide("S03", tag = "lit")
  deck_kicker(tx("S03.kicker")); deck_title(tx("S03.title"))
  f <- list(ci = f_ci_level(), lim = f_limits())

  # 위: 지침 카드 세 개(연도순)
  y0 <- GEO$BODY_TOP + 0.05; gw <- 0.25; cw <- (GEO$CW - 2 * gw) / 3; chh <- 2.75
  cards <- c("ema", "fda16", "fda25"); fills <- c(PAL$tint_blue, PAL$tint_blue, PAL$tint_orange)
  for (k in seq_along(cards)) {
    x <- GEO$ML + (k - 1) * (cw + gw)
    deck_text(tx(sprintf("S03.cards.%s", cards[k]), f), c(x, y0, cw, chh), size = 17, bg = fills[k], geom = "roundRect", label = sprintf("card_%s", cards[k]), gap_pt = 7)
  }

  deck_text(tx("S03.ema_caveat"), c(GEO$ML + 0.02, y0 + chh - 0.45, cw - 0.04, 0.42), size = 16, color = PAL$muted, label = "ema_caveat")

  # 아래: 두필루맙에서의 문제(왼쪽 표지, 오른쪽 요점)
  yb <- y0 + chh + 0.25; hb <- GEO$BODY_BOTTOM - yb; lw <- 2.3
  deck_text(tx("S03.problem.label"), c(GEO$ML, yb, lw, 0.95), size = 18, bold = TRUE, bg = PAL$tint_grey, geom = "roundRect", label = "problem_label")
  deck_bullets(tx("S03.problem.bullets"), c(GEO$ML + lw + 0.2, yb, GEO$CW - lw - 0.2, hb), size = 16, gap_pt = 6)

  deck_notes(tx("S03.notes", f))
  deck_end()
}
