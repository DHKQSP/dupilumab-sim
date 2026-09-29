# S03 배경: 규제 기대(단회 투여 PK 시험의 1차 평가변수 AUC0-inf)와 비교 유효성 시험 면제 흐름, 두필루맙에서의 문제. [문헌]
# 기관·연도 카드 세 개(EMA 2012, FDA 2016, FDA 2025 초안) + 아래 두필루맙 문제 요점. 수치는 설계 상수(동등성 한계, CI 수준)뿐.
# 원문과 직접 대조하지 못한 인용(EMA 2012, FDA 2016)은 카드 머리에 같은 방식으로 '원문 대조 전'을 적고, 카드 아래 한 줄에 사유와 계획, 자세한 사정은 노트에 둔다.

# 글자 없는 카드 배경(16pt 빈 글상자): 표지 글상자를 카드 세로 가운데에 따로 놓기 위해
s03_panel <- function(box, fill, label) deck_text(" ", box, size = 16, bg = fill, geom = "roundRect", label = label)

slide_S03 <- function() {
  deck_slide("S03", tag = "lit")
  deck_kicker(tx("S03.kicker")); deck_title(tx("S03.title"))
  f <- list(ci = f_ci_level(), lim = f_limits())

  # 위: 지침 카드 세 개(연도순)
  y0 <- GEO$BODY_TOP + 0.05; gw <- 0.25; cw <- (GEO$CW - 2 * gw) / 3; chh <- 2.62
  cards <- c("ema", "fda16", "fda25"); fills <- c(PAL$tint_blue, PAL$tint_blue, PAL$tint_orange)
  for (k in seq_along(cards)) {
    x <- GEO$ML + (k - 1) * (cw + gw)
    deck_text(tx(sprintf("S03.cards.%s", cards[k]), f), c(x, y0, cw, chh), size = 16, bg = fills[k], geom = "roundRect", label = sprintf("card_%s", cards[k]), gap_pt = 7)
  }

  # 카드 바로 아래: '원문 대조 전' 표시의 사유와 계획 한 줄(자세한 사정은 노트)
  yn <- y0 + chh + 0.04; hn <- 0.39
  deck_text(tx("S03.verify_note"), c(GEO$ML, yn, GEO$CW, hn), size = 16, color = PAL$ink2, label = "text_verify")

  # 아래: 두필루맙에서의 문제(왼쪽 표지 카드는 요점 묶음 높이에 맞추고 글은 세로 가운데, 오른쪽 요점)
  yb <- yn + hn + 0.12; hb <- GEO$BODY_BOTTOM - yb; lw <- 1.9
  bl <- tx("S03.problem.bullets"); bw <- GEO$CW - lw - 0.2
  hb_used <- min(hb, est_height(bl, bw, 16, 6, indent = 0.3) + 0.08)
  s03_panel(c(GEO$ML, yb, lw, hb_used), PAL$tint_grey, "problem_panel")
  hl <- 0.82
  deck_text(tx("S03.problem.label"), c(GEO$ML + 0.02, yb + (hb_used - hl) / 2, lw - 0.04, hl), size = 18, bold = TRUE, label = "problem_label")
  deck_bullets(bl, c(GEO$ML + lw + 0.2, yb, bw, hb), size = 16, gap_pt = 6)

  deck_notes(tx("S03.notes", f))
  deck_end()
}
