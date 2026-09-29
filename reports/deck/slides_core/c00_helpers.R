# c00_helpers.R: 핵심 덱(v1.1) 슬라이드 공용 도우미. 파일 이름 순서로 가장 먼저 읽힌다(build_deck.R --deck core).
# 배치 원칙(지시 2026-09-29 §1): 제목(결론 문장, 28pt, 2줄 이내) → 주 그림(본문 영역의 60% 이상) → 본문(18pt, 3줄 이내) → 각주(14pt 이상).

# 제목: 줄 수에 맞춘 높이로 두고 본문 시작 위치를 돌려준다(한 줄 제목이면 그림이 위로 올라간다)
# 제목 안의 "\n"은 줄을 나누는 자리다(LibreOffice는 한글 옆의 줄바꿈 금지 문자를 무시해 '0.005% / 이하'처럼 갈라지므로 끊을 곳을 직접 정한다).
core_title <- function(s, kicker = NULL, size = SZ$title, box_x = GEO$ML, top = GEO$TITLE_TOP, width = GEO$CW - 0.1, color = NULL) {
  if (!is.null(kicker)) deck_kicker(kicker)
  ln <- strsplit(s, "\n", fixed = TRUE)[[1]]
  n <- sum(vapply(ln, est_lines, 1L, width = width - 0.2, size = size, bold = TRUE)); h <- n * size * 1.2 * 1.05 / 72 + 0.14
  if (length(ln) == 1L) deck_title(s, size = size, box = c(box_x, top, width, h)) else {
    DK$cur$title <- strip_markup(paste(ln, collapse = " ")); fit_check("title", ln, c(box_x, top, width, h), size, gap_pt = 0, bold = TRUE, line = 1.05)
    DK$cur$title_lines <- n
    if (!is.null(LIMITS$title_lines) && n > LIMITS$title_lines && isTRUE(DK$strict)) stop(sprintf("%s: title needs %d lines (limit %d)", REG$sec, n, LIMITS$title_lines), call. = FALSE)
    col <- color %||% (if (DK$cur$dark) PAL$dark_ink else PAL$ink)
    DK$x <- ph_with(DK$x, do.call(block_list, lapply(ln, function(z) para(z, size, col, bold = TRUE, gap_pt = 0, line = 1.05))), location = loc(c(box_x, top, width, h), "title"))
  }
  top + h + 0.16
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
  pv <- lapply(value, function(ln) do.call(fpar, c(lapply(ln, function(v) ftext(nobreak(v[[1]]), ftp(vsz(v), v[[2]], vsz(v) >= 30))),
                                                   list(fp_p = fp_par(text.align = "left", padding.bottom = 6, line_spacing = 0.95)))))
  labs_ <- strsplit(label, "\n", fixed = TRUE)[[1]]
  p2 <- lapply(seq_along(labs_), function(i) para(labs_[i], SZ$body, PAL$ink, FALSE, "left", gap_pt = if (i < length(labs_)) 4 else 0))
  # 두 도형: 카드 바탕 + 위에 머리말(위 정렬, 카드끼리 줄 맞춤), 머리말 아래 영역에 수치·설명(같은 바탕색, 위 정렬). 설명의 "\n"은 문단 나눔
  hh <- SZ$body * 1.2 * 1.1 / 72 + 2 * CARD_INS[["tb"]] + 0.04
  DK$x <- ph_with(DK$x, para(head, SZ$body, head_color, TRUE, "left", gap_pt = 0), location = loc(box, "card", bg = bg, geom = "roundRect", ln = no_line()))
  DK$x <- ph_with(DK$x, do.call(block_list, c(pv, p2)), location = loc(c(box[1], box[2] + hh, box[3], box[4] - hh - 0.06), "cardval", bg = bg, geom = "rect", ln = no_line()))   # 위 정렬: 카드끼리 큰 수치 줄 맞춤
  invisible(NULL)
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

# ---- 대표 대상자 음영 그림(S7: 2016 모델, 별첨 A4a: 2020 모델) ------------------------------------------------------------------
# 패널마다: 선형 눈금 참 농도 곡선(Day 1~70), Day 1~관측 tlast 아래 파랑 음영(AUClast 쪽), tlast 이후 주황 음영(외삽), 채혈점(정량 = 채운 점,
# 정량한계 미만 = 빈 점), 정량한계 선, tlast 세로선, 위에 큰 글자 'AUClast/AUCinf = xx.x%'. 오른쪽 위 삽입도(로그 눈금, tlast 7일 전 ~ Day 60):
# 참 곡선, 비구획 외삽선 Clast x exp(-λz (t - tlast)) 점선, 비구획 외삽 면적 빗금 대 참 외삽 면적 주황, 두 면적의 배수.
# 입력: results/core_deck/rep_subjects.csv, rep_profiles.csv, rep_obs.csv(scripts/63). 글자는 모두 14pt 이상.
core_hatch <- function(x0, x1, u_top, uf, xr, ur, aspect, n = 26) {
  s <- (ur / xr) * aspect                                            # 화면에서 45도가 되는 기울기(로그 단위 / 일)
  xs <- seq(x0, x1, length.out = 400); a_rng <- c(uf - s * x1, max(u_top(xs)) - s * x0)
  rbindlist(lapply(seq(a_rng[1], a_rng[2], length.out = n), function(a) {
    u <- a + s * xs; ok <- u >= uf & u <= u_top(xs); if (sum(ok) < 2) return(NULL)
    i <- range(which(ok)); data.table(x = xs[i[1]], xend = xs[i[2]], y = 10^u[i[1]], yend = 10^u[i[2]]) }))
}
core_shaded_panels <- function(model, L) {
  RSf <- "core_deck/rep_subjects.csv"; RPf <- "core_deck/rep_profiles.csv"; ROf <- "core_deck/rep_obs.csv"
  S <- rows(RSf, sprintf("model=='%s'", model)); premise(nrow(S) == 3 && setequal(S$role, c("median", "p05", "min")), sprintf("%s: three representative subjects", model))
  S <- S[match(c("median", "p05", "min"), role)]
  lloq_v <- .read("config/assay.yaml")$lloq_mg_L$value
  mk <- function(i) {
    s <- S[i]; pr <- rows(RPf, sprintf("model=='%s' & id==%d", model, s$id)); ob <- rows(ROf, sprintf("model=='%s' & id==%d", model, s$id))
    pr <- rbind(data.table(study_day = 1, conc = 0), pr[, .(study_day, conc)])
    tl <- s$study_day_tlast; ctl <- approx(pr$study_day, pr$conc, tl)$y
    pa <- rbind(pr[study_day < tl], data.table(study_day = tl, conc = ctl)); pb <- rbind(data.table(study_day = tl, conc = ctl), pr[study_day > tl])
    q <- ob[blq == FALSE]; b <- ob[blq == TRUE]
    ymax <- max(c(pr$conc, q$conc_obs)) * 1.15
    head <- fill(L$head, list(v = fnum(100 * s$coverage_true, 1)))
    p <- ggplot() +
      geom_ribbon(data = pa, aes(study_day, ymin = 0, ymax = conc), fill = PAL$blue, alpha = 0.30) +
      geom_ribbon(data = pb, aes(study_day, ymin = 0, ymax = conc), fill = PAL$orange, alpha = 0.85) +
      geom_line(data = pr, aes(study_day, conc), colour = PAL$ink, linewidth = 0.8) +
      geom_hline(yintercept = lloq_v, colour = PAL$ink2, linetype = "22", linewidth = 0.4) +
      geom_vline(xintercept = tl, colour = PAL$ink2, linetype = "42", linewidth = 0.5) +
      geom_point(data = q, aes(study_day, conc_obs), shape = 21, fill = PAL$blue, colour = "white", size = 2.3, stroke = 0.4) +
      geom_point(data = b, aes(study_day, 0), shape = 21, fill = "white", colour = PAL$blue, size = 2.1, stroke = 0.8) +
      scale_x_continuous(limits = c(1, 70), breaks = c(1, 15, 29, 43, 57, 70), expand = expansion(mult = 0)) +
      scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.02, 0))) +
      coord_cartesian(clip = "off") +
      labs(x = L$xlab, y = NULL, title = head, subtitle = fill(L$sub, list(role = L$role[[s$role]]))) + theme_core(16) +
      theme(plot.title = element_text(size = 20, face = "bold", colour = PAL$blue, margin = margin(0, 0, 2, 0)), plot.subtitle = element_text(size = 15, colour = PAL$ink2),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), plot.margin = margin(4, 10, 4, 4))
    # ---- 삽입도(로그 눈금) ----
    x0 <- max(1, tl - 7); x1 <- 60; fl_ <- 0.01
    tail <- pr[study_day >= x0 & study_day <= x1 & conc > 0]
    top <- max(tail$conc) * 2.5; ur <- log10(top) - log10(fl_); xr <- x1 - x0
    ti <- pb[study_day <= x1 & conc >= fl_]
    ip <- ggplot() + geom_ribbon(data = ti, aes(study_day, ymin = fl_, ymax = conc), fill = PAL$orange, alpha = 0.85)
    if (isTRUE(s$lambda_ok)) {
      u_top <- function(x) pmin(log10(s$Clast) - s$lambda_z * (x - tl) / log(10), log10(top))
      hz <- core_hatch(tl, x1, u_top, log10(fl_), xr, ur, aspect = 1.9 / 1.5)
      nl <- data.table(study_day = seq(tl, x1, by = 0.1))[, conc := s$Clast * exp(-s$lambda_z * (study_day - tl))]
      r <- s$extrap_ratio_nca_to_true; lab <- fill(L$ratio, list(r = fnum(r, if (r < 1) 2 else 1)))
      ip <- ip + geom_segment(data = hz, aes(x = x, xend = xend, y = y, yend = yend), colour = PAL$ink2, linewidth = 0.35) +
        geom_line(data = nl[conc >= fl_], aes(study_day, conc), colour = PAL$ink, linetype = "22", linewidth = 0.7)
    } else lab <- L$no_lz
    ip <- ip + geom_line(data = pr[study_day >= x0 & study_day <= x1 & conc >= fl_], aes(study_day, conc), colour = PAL$ink, linewidth = 0.8) +
      geom_hline(yintercept = lloq_v, colour = PAL$ink2, linetype = "22", linewidth = 0.4) +
      scale_y_log10(limits = c(fl_, top), breaks = c(0.01, 0.1, 1, 10), labels = function(x) formatC(x, format = "fg"), expand = expansion(mult = 0)) +
      scale_x_continuous(limits = c(x0, x1), breaks = pretty(c(x0, x1), 3), expand = expansion(mult = 0)) +
      labs(x = NULL, y = NULL, title = lab) + theme_core(16) +
      theme(plot.title = element_text(size = 14, face = "bold", colour = PAL$ink, margin = margin(0, 0, 2, 0), lineheight = 0.95), panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", colour = PAL$grid, linewidth = 0.6), plot.margin = margin(3, 6, 2, 3))
    p + patchwork::inset_element(ip, left = 0.45, bottom = 0.24, right = 1.0, top = 1.0, align_to = "panel")
  }
  patchwork::wrap_plots(lapply(1:3, mk), nrow = 1)
}
