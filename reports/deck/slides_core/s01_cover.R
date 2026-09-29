# S1 표지(지시 2026-09-29 §2): 제목·부제는 지시 문구 그대로. 제목의 약어(PK, AUCinf)는 같은 슬라이드의 용어 줄에서 풀이한다
# (check.abbrev_anywhere_slides). 어두운 배경.
slide_S1 <- function() {
  deck_slide("S1", tag = "none", dark = TRUE)
  deck_text(tx("S1.kicker"), c(GEO$ML, 1.20, 11.6, 0.45), size = 18, bold = TRUE, color = "#8fb8ea", label = "kicker")
  deck_title(tx("S1.title"), size = 36, box = c(GEO$ML, 1.75, 12.0, 1.50))
  deck_text(tx("S1.subtitle", list(wt = f_wt_range(), dose = f_dose())), c(GEO$ML, 3.45, 12.0, 1.00), size = 22, color = "#d5dde6", label = "subtitle")
  deck_text(tx("S1.glossary"), c(GEO$ML, 4.75, 12.0, 1.05), size = 18, color = "#aab6c2", label = "glossary", gap_pt = 2)
  deck_text(tx("S1.meta", list(version = DK$version, date = DK$date)), c(GEO$ML, 6.15, 8, 0.45), size = 16, color = "#aab6c2", label = "caption_meta")
  deck_notes(tx("S1.notes"))
  deck_end()
}
