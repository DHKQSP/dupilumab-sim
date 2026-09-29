# S01 표지: 시험, 제안, 보고 성격, 버전·날짜·커밋. 어두운 배경.
slide_S01 <- function() {
  deck_slide("S01", tag = "none", dark = TRUE)
  deck_text(tx("S01.kicker"), c(GEO$ML, 1.35, 9, 0.5), size = 18, bold = TRUE, color = "#8fb8ea", label = "kicker")
  deck_title(tx("S01.title"), size = 36, box = c(GEO$ML, 1.9, 11.6, 1.45))
  deck_text(tx("S01.subtitle"), c(GEO$ML, 3.42, 11.6, 1.56), size = 22, color = "#d5dde6", label = "subtitle")
  deck_text(tx("S01.scope", list(wt = f_wt_range(), dose = f_dose())), c(GEO$ML, 5.1, 11.6, 0.6), size = 18, color = "#aab6c2", label = "scope")
  # 커밋 표기는 바닥글과 같은 문자열(작업 트리가 커밋과 다르면 '*')
  deck_text(tx("S01.meta", list(version = DK$version, date = DK$date, commit = paste0(DK$commit, if (isTRUE(DK$clean)) "" else "*"))), c(GEO$ML, 6.05, 8, 0.5), size = 16, color = "#aab6c2", label = "meta")
  deck_notes(tx("S01.notes"))
  deck_end()
}
