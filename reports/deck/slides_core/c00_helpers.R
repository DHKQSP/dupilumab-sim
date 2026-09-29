# c00_helpers.R: 핵심 덱(v1.1) 슬라이드 공용 도우미. 파일 이름 순서로 가장 먼저 읽힌다(build_deck.R --deck core).
# 배치 원칙(지시 2026-09-29 §1): 제목(결론 문장, 28pt, 2줄 이내) → 주 그림(본문 영역의 60% 이상) → 본문(18pt, 3줄 이내) → 각주(14pt 이상).

# 제목: 줄 수에 맞춘 높이로 두고 본문 시작 위치를 돌려준다(한 줄 제목이면 그림이 위로 올라간다)
core_title <- function(s, kicker = NULL) {
  if (!is.null(kicker)) deck_kicker(kicker)
  n <- est_lines(s, GEO$CW - 0.1 - 0.2, SZ$title, TRUE); h <- n * SZ$title * 1.2 * 1.05 / 72 + 0.14
  deck_title(s, box = c(GEO$ML, GEO$TITLE_TOP, GEO$CW - 0.1, h))
  GEO$TITLE_TOP + h + 0.16
}
# 본문(18pt): 문단 벡터. 아래 끝(bottom)에 맞춰 높이를 정하고 위쪽 y를 돌려준다
core_body_h <- function(paras, width = GEO$CW, size = SZ$body, gap_pt = 4) est_height(paras, width, size, gap_pt) + 0.02
core_body <- function(paras, bottom, x = GEO$ML, width = GEO$CW, size = SZ$body, gap_pt = 4, label = "body") {
  h <- core_body_h(paras, width, size, gap_pt); deck_text(paras, c(x, bottom - h, width, h), size = size, label = label, gap_pt = gap_pt); bottom - h
}
# 각주·캡션(14~16pt, 본문 줄 수에 세지 않음)
core_caption <- function(paras, bottom, x = GEO$ML, width = GEO$CW, size = 14, color = PAL$ink2, label = "caption", gap_pt = 2) {
  h <- est_height(paras, width, size, gap_pt) + 0.02; deck_text(paras, c(x, bottom - h, width, h), size = size, color = color, label = label, gap_pt = gap_pt); bottom - h
}
# 카드: 머리말(18pt 굵게) + 큰 수치 줄(한 줄 이상; 줄마다 조각 list(text, color[, size]) 목록) + 설명(18pt)
core_card <- function(head, value, label, box, bg, value_size = 48, head_color = PAL$ink) {
  if (!is.list(value[[1]][[1]])) value <- list(value)                  # 한 줄이면 줄 목록으로 감싼다
  vsz <- function(v) if (length(v) >= 3) v[[3]] else value_size
  hv <- sum(vapply(value, function(ln) max(vapply(ln, vsz, 1)), 1)) * 1.2 * 0.95 / 72 + 0.08 * length(value)
  fit_check("card", c(head, label), c(box[1], box[2], box[3], box[4] - hv), SZ$body, gap_pt = 6, card = TRUE)
  for (ln in value) { wv <- sum(vapply(ln, function(v) text_w(v[[1]], vsz(v), vsz(v) >= 30), 1))
    if (wv > box[3] - 2 * CARD_INS[["lr"]]) stop(sprintf("%s card value line '%s' %.2f in wider than %.2f in", REG$sec, paste(vapply(ln, `[[`, "", 1), collapse = ""), wv, box[3] - 2 * CARD_INS[["lr"]]), call. = FALSE) }
  p0 <- para(head, SZ$body, head_color, TRUE, "left", gap_pt = 6)
  pv <- lapply(value, function(ln) do.call(fpar, c(lapply(ln, function(v) ftext(nobreak(v[[1]]), ftp(vsz(v), v[[2]], vsz(v) >= 30))),
                                                   list(fp_p = fp_par(text.align = "left", padding.bottom = 6, line_spacing = 0.95)))))
  p2 <- para(label, SZ$body, PAL$ink, FALSE, "left", gap_pt = 0)
  DK$x <- ph_with(DK$x, do.call(block_list, c(list(p0), pv, list(p2))), location = loc(box, "stat", bg = bg, geom = "roundRect", ln = no_line())); invisible(NULL)
}
# 도식 상자(글자 있음): 머리말 + 내용 줄. 글자 없는 상자는 deck_box
core_diag_box <- function(head, lines, box, bg = PAL$tint_grey, head_color = PAL$ink, size = SZ$body, label = "diag") {
  fit_check(label, c(head, lines), box, size, gap_pt = 4, card = TRUE)
  ps <- c(list(para(head, size, head_color, TRUE, "left", gap_pt = 6)), lapply(lines, function(z) para(z, size, PAL$ink, FALSE, "left", gap_pt = 4)))
  DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(box, label, bg = bg, geom = "roundRect", ln = no_line())); invisible(NULL)
}
# 화살표(글자 없는 도형: 검사 6은 빈 글자 도형을 건너뛴다)
core_arrow <- function(box, fill = PAL$muted, geom = "rightArrow") deck_box(box, fill = fill, geom = geom, label = "arrow")

# 인쇄 값 보조: '이상' 하한은 내림, '이하' 상한은 올림(사전 등록 7b, D-063)
fl <- function(x, d) floor(x * 10^d + 1e-9) / 10^d
cl <- function(x, d) ceiling(x * 10^d - 1e-9) / 10^d
# 두 모델 창 포착률(기본 조건) 최솟값의 내림: '84% 이상'
core_cov_floor <- function(d = 0) {
  f <- "core_deck/coverage_by_case.csv"; w <- "case %in% c('k2016_base','k2020_base')"; r <- rows(f, w); premise(nrow(r) == 2, "two base cases")
  x <- 100 * min(r$min); dderived("window coverage, base case, smallest subject over both models, rounded down", f, sprintf("%s :: floor(min(min) x 100)", w), x, paste0(fnum(fl(x, d), d), "%"))
}
# 공개 SAP의 adjusted R² 기준: 사전 등록(2026-09-26) 6절 public_saps 문장에서 정규식으로 읽는다(결과보고 덱 A4와 같은 locator)
core_sap <- function(nct, rx, item) {
  P <- c("section6", "s2_1_failure_by_set", "public_saps")
  x <- .read("config/prereg_20260926.yaml")[[P[1]]][[P[2]]][[P[3]]]; m <- regmatches(x, regexec(sprintf("%s \\(adjusted R-squared %s ([0-9.]+)", nct, rx), x))[[1]]
  premise(length(m) == 2, sprintf("public SAP threshold for %s (%s)", nct, rx))
  v <- as.numeric(m[2]); dderived(item, "config/prereg_20260926.yaml", sprintf("%s :: regex '%s \\(adjusted R-squared %s ([0-9.]+)'", paste(P, collapse = "."), nct, rx), v, fnum(v, 2))
}
