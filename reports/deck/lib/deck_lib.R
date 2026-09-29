# deck_lib.R: 결과보고 슬라이드(v1.0.1) 공용 함수. officer로 16:9 슬라이드를 절대 위치로 그린다.
# 원칙
#  - 수치: 모두 결과·config 파일에서 읽고 regulatory/src/reg_helpers.R의 .record로 출처를 남긴다(아래 d* 함수). 문구 파일에는 숫자를 쓰지 않는다.
#  - 문구: text/<lang>/*.yaml. 자리표시자 {name}은 슬라이드 코드가 만든 수치 문자열로 바꾼다. 같은 코드로 영문판을 만들 수 있다.
#  - 형식: 본문 16pt 이상, 표 12pt 이상, 슬라이드당 요점 5개 이하, 표 본문 6행 이하. 어기면 빌드가 멈춘다.
#  - 넘침: 글상자마다 Pretendard 글꼴 폭으로 줄 수를 추정해 상자 높이를 넘으면 빌드가 멈춘다(렌더링 육안 확인과 별도).
suppressPackageStartupMessages({ library(officer); library(flextable); library(ggplot2); library(data.table) })
source(proj_path("regulatory", "src", "reg_helpers.R"))

# 추적 행은 영문 추적표(regulatory/traceability.csv)에 합쳐지므로 한글을 넣지 않는다(단위 "명", "일" 등은 문구 파일에)
.record_reg <- .record
.record <- function(item, rel, locator, raw, printed) {
  if (grepl("[\u3131-\u318e\uac00-\ud7a3]", paste(item, locator, printed), perl = TRUE))
    stop(sprintf("%s: trace item, locator and printed value must not contain Korean (put Korean units in the text file): %s | %s | %s", REG$sec, item, locator, printed), call. = FALSE)
  .record_reg(item, rel, locator, raw, printed)
}
DK <- new.env()
GEO <- list(W = 13.333, H = 7.5, ML = 0.55, MR = 0.55, KICK_TOP = 0.28, TITLE_TOP = 0.52, TITLE_H = 1.0, BODY_TOP = 1.65, BODY_BOTTOM = 6.80,
            FOOT_TOP = 7.03, FOOT_H = 0.30, TAG_W = 1.15, TAG_H = 0.34)
GEO$CW <- GEO$W - GEO$ML - GEO$MR                     # content width
PAL <- list(blue = "#2a78d6", orange = "#eb6834", green = "#1baf7a", ink = "#0b0b0b", ink2 = "#52514e", muted = "#8a8984", grid = "#e6e5e0",
            surface = "#ffffff", tint_blue = "#eaf2fc", tint_orange = "#fdf0ea", tint_grey = "#f3f3f1", dark = "#16202b", dark_ink = "#f4f6f8")
FONT <- "Pretendard"
SZ <- list(title = 26, kicker = 13, body = 18, body_min = 16, small = 16, table = 13, table_min = 12, foot = 9, tag = 11, stat = 40, stat_label = 16)
LIMITS <- list(bullets = 5L, table_rows = 6L)
# 덱 종류(D-063): results = 결과보고 덱 v1.0.1(32장, 기술 백업), core = 핵심 덱 v1.1(정량 근거 중심). 기본값은 results(기존 빌드와 같은 결과).
# core는 글자 하한(제목 28, 본문 18, 표·그림·캡션 14pt), 본문 슬라이드의 본문 줄 수(3줄)와 주 그림 면적(본문 영역의 60%)을 빌드에서 강제한다.
GEO0 <- GEO; SZ0 <- SZ; LIMITS0 <- LIMITS
DECK_PROFILES <- list(
  results = list(text_suffix = "", slides = "slides", fig_dir = "figures", base = "dupilumab_endpoint_results_%s", trace = "deck_traceability.csv", meta = "deck_meta.csv",
                 doc = "deck", footer = "full", geo = list(), sz = list(), limits = list()),
  core = list(text_suffix = "_core", slides = "slides_core", fig_dir = "figures_core", base = "dupilumab_AUCinf_core_deck_%s", trace = "core_deck_traceability.csv", meta = "core_deck_meta.csv",
              doc = "deck_core", footer = "version_source",
              geo = list(KICK_TOP = 0.30, TITLE_TOP = 0.66, TITLE_H = 1.10, BODY_TOP = 1.86, BODY_BOTTOM = 6.84, FOOT_TOP = 7.06, TAG_W = 1.30, TAG_H = 0.40),
              sz = list(title = 28, kicker = 16, body = 18, body_min = 18, small = 18, table = 14, table_min = 14, foot = 10, tag = 14, stat = 44, stat_label = 18, caption_min = 14, fig_min = 14),
              limits = list(bullets = 3L, table_rows = 8L, body_lines = 3L, fig_share = 0.60, title_lines = 2L))
)

# ---- 초기화 ------------------------------------------------------------------------------------------------------------------------
deck_init <- function(lang = "ko", strict = TRUE, fig_dir = NULL, deck = "results") {
  prof <- DECK_PROFILES[[deck]]; if (is.null(prof)) stop("unknown deck: ", deck)
  DK$deck <- deck; DK$prof <- prof
  GEO <<- modifyList(GEO0, prof$geo); GEO$CW <<- GEO$W - GEO$ML - GEO$MR; SZ <<- modifyList(SZ0, prof$sz); LIMITS <<- modifyList(LIMITS0, prof$limits)
  if (is.null(fig_dir)) fig_dir <- proj_path("reports", "deck", prof$fig_dir)
  DK$lang <- lang; DK$strict <- strict; DK$fig_dir <- fig_dir; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  DK$txt <- deck_text_load(paste0(lang, prof$text_suffix)); DK$meta <- list(); DK$fit <- list(); DK$n <- 0L; DK$total <- NA_integer_; DK$headline <- list()
  git <- function(...) suppressWarnings(system2("git", c("-C", PROJ_ROOT, ...), stdout = TRUE, stderr = FALSE))
  DK$commit <- git("rev-parse", "--short=7", "HEAD")[1]; DK$date <- git("log", "-1", "--format=%cs")[1]   # 날짜 = 소스 커밋 날짜(빌드마다 같음)
  dirty <- git("status", "--porcelain", "--untracked-files=no")
  dirty <- dirty[!grepl("^.. reports/deck/(dupilumab_endpoint_results|deck_traceability|deck_meta|figures/|dupilumab_AUCinf_core_deck|core_deck_|figures_core/)", dirty)]   # 덱 산출물 자체의 변경은 제외(생성 대상)
  DK$clean <- length(dirty) == 0
  DK$version <- DK$txt$common$version
  reg_init(prof$doc); DK$x <- read_pptx(proj_path("reports", "deck", "template_16x9.pptx"))
  s <- slide_size(DK$x); stopifnot(abs(s$width - GEO$W) < 0.01, abs(s$height - GEO$H) < 0.01)
  invisible(DK)
}
deck_text_load <- function(lang) {
  d <- proj_path("reports", "deck", "text", lang); fs <- sort(list.files(d, pattern = "\\.ya?ml$", full.names = TRUE))
  out <- list(); for (f in fs) { y <- yaml::read_yaml(f); dup <- intersect(names(y), names(out)); if (length(dup)) stop("duplicate text keys in ", f, ": ", paste(dup, collapse = ", ")); out <- c(out, y) }
  out
}

# ---- 문구 --------------------------------------------------------------------------------------------------------------------------
# tx("S09.title", facts) -> 문자열(또는 목록이면 문자 벡터). {name}을 facts[[name]]로 바꾸고, 남은 자리표시자나 없는 키는 오류.
tx_raw <- function(key) { p <- strsplit(key, ".", fixed = TRUE)[[1]]; y <- DK$txt; for (k in p) { if (is.null(y[[k]])) stop("missing text key: ", key, call. = FALSE); y <- y[[k]] }; y }
tx <- function(key, facts = list()) {
  y <- tx_raw(key); if (is.list(y)) y <- unlist(y)
  vapply(y, function(s) fill(s, facts, key), "", USE.NAMES = FALSE)
}
fill <- function(s, facts, key = "") {
  keys <- unique(regmatches(s, gregexpr("\\{[A-Za-z0-9_.]+\\}", s))[[1]])
  for (k in keys) { nm <- substr(k, 2, nchar(k) - 1); if (is.null(facts[[nm]])) stop(sprintf("text %s: no fact for %s", key, k), call. = FALSE)
    s <- gsub(k, as.character(facts[[nm]]), s, fixed = TRUE) }
  if (grepl("[{}]", s)) stop(sprintf("text %s: unresolved braces in '%s'", key, s), call. = FALSE)
  if (grepl("–|—|−", s)) stop(sprintf("text %s: en-dash, em-dash or U+2212 in '%s'", key, s), call. = FALSE)
  if (grepl("(^|[^A-Za-z])(NA|NaN|Inf)([^A-Za-z]|$)", s)) stop(sprintf("text %s: NA/NaN/Inf in '%s'", key, s), call. = FALSE)
  s
}

# ---- 수치(출처 기록) -----------------------------------------------------------------------------------------------------------------
# 모든 함수는 인쇄 문자열을 돌려주고 .record로 추적 행을 남긴다(document = "deck", section = 슬라이드 ID).
rng_fmt <- function(a, b, d, unit = "") {
  A <- fnum(a, d); B <- fnum(b, d)
  if (A == B) return(paste0(A, unit))
  if (DK$lang == "ko") sprintf("%s~%s%s", A, B, unit) else sprintf("%s%s to %s%s", A, unit, B, unit)
}
ci_fmt <- function(e, l, h, d, unit = "%") if (DK$lang == "ko") sprintf("%s%s (95%% CI %s~%s)", fnum(e, d), unit, fnum(l, d), fnum(h, d)) else sprintf("%s%s (95%% CI %s to %s)", fnum(e, d), unit, fnum(l, d), fnum(h, d))
dv <- function(rel, where, col, d = 1, unit = "", item = col, scale = 1, big = FALSE) {
  r <- row1(rel, where); x <- r[[col]] * scale; p <- paste0(fnum(x, d, big), unit)
  .record(item, rel, sprintf("%s :: %s%s", where, col, if (scale != 1) sprintf(" x %s", scale) else ""), x, p); p
}
dint <- function(rel, where, col, item = col, unit = "") { r <- row1(rel, where); x <- r[[col]]; p <- paste0(fint(x), unit); .record(item, rel, sprintf("%s :: %s", where, col), x, p); p }
drange <- function(rel, where, col, d = 1, unit = "", item = col, scale = 1) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- range(r[[col]] * scale)
  p <- rng_fmt(x[1], x[2], d, unit); .record(item, rel, sprintf("range over %d rows [%s] :: %s%s", nrow(r), where, col, if (scale != 1) sprintf(" x %s", scale) else ""), x, p); p
}
# 여러 행의 lo 열 최솟값 ~ hi 열 최댓값(예: 두 모델의 5~95백분위 범위)
dspan <- function(rel, where, lo_col, hi_col, d = 1, unit = "", item = paste(lo_col, hi_col)) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- c(min(r[[lo_col]]), max(r[[hi_col]]))
  p <- rng_fmt(x[1], x[2], d, unit); .record(item, rel, sprintf("over %d rows [%s] :: min(%s), max(%s)", nrow(r), where, lo_col, hi_col), x, p); p
}
dext <- function(rel, where, col, fun = max, d = 1, unit = "", item = col, scale = 1) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- fun(r[[col]] * scale)
  p <- paste0(fnum(x, d), unit); .record(item, rel, sprintf("%s over %d rows [%s] :: %s", deparse(substitute(fun)), nrow(r), where, col), x, p); p
}
dci <- function(rel, where, est = "pass_pct", lo = "lo", hi = "hi", d = 2, unit = "%", item = est, scale = 1) {
  r <- row1(rel, where); e <- r[[est]] * scale; l <- r[[lo]] * scale; h <- r[[hi]] * scale; p <- ci_fmt(e, l, h, d, unit)
  .record(item, rel, sprintf("%s :: %s [%s, %s]", where, est, lo, hi), c(e, l, h), p); p
}
dcount <- function(rel, where, item = "count") { r <- rows(rel, where); p <- as.character(nrow(r)); .record(item, rel, sprintf("count of rows [%s]", where), nrow(r), p); p }
dcfg <- function(file, path, item = paste(path, collapse = "."), fmt = function(x) paste(format(x), collapse = ", ")) {
  y <- .read(file.path("config", file)); for (k in path) y <- y[[k]]
  premise(!is.null(y), sprintf("config/%s: %s not found", file, paste(path, collapse = "."))); p <- fmt(unlist(y))
  .record(item, file.path("config", file), paste(path, collapse = "."), unlist(y), p); p
}
# 파생값(차이·비 등): 계산식을 locator에 적는다
dderived <- function(item, rel, locator, raw, printed) { .record(item, rel, locator, raw, printed); printed }
# 표·그림 전체 출처
dsrc <- function(item, rels, what = "(figure)") { for (r_ in rels) .record(item, r_, "whole object", NA, what); invisible(NULL) }
# 핵심 수치(표지 요약·논리 도식): M&S 보고서 추적표·key_numbers_en.md와의 일치를 check_deck.R가 반드시 확인한다
headline <- function(p) { DK$headline[[length(DK$headline) + 1L]] <- list(section = REG$sec, printed = p, idx = length(REG$trace)); p }

# ---- 글자 폭과 넘침 추정 ------------------------------------------------------------------------------------------------------------------
strip_markup <- function(s) gsub("\\*\\*|__", "", s)
text_w <- function(s, size, bold = FALSE) systemfonts::string_width(s, family = FONT, size = size, res = 72, bold = bold) / 72
est_lines <- function(s, width, size, bold = FALSE) {
  if (grepl("\n", s, fixed = TRUE)) return(sum(vapply(strsplit(s, "\n", fixed = TRUE)[[1]], est_lines, 1L, width = width, size = size, bold = bold)))   # 명시적 줄바꿈
  # 굵게 표시(**..**, __..__)한 부분은 굵은 글꼴 폭으로 잰다: 단어마다 (굵은 부분, 보통 부분) 폭을 더한다
  seg <- strsplit(s, "\\*\\*|__", perl = TRUE)[[1]]; if (!length(seg)) return(1L)
  pieces <- list(); cur_w <- character(0); cur_b <- logical(0)
  for (k in seq_along(seg)) { b_ <- bold || (k %% 2 == 0); tok <- strsplit(seg[k], " ", fixed = TRUE)[[1]]
    if (!length(tok)) next
    for (j in seq_along(tok)) { if (j > 1) { pieces[[length(pieces) + 1L]] <- list(w = cur_w, b = cur_b); cur_w <- character(0); cur_b <- logical(0) }
      if (nzchar(tok[j])) { cur_w <- c(cur_w, tok[j]); cur_b <- c(cur_b, b_) } }
    if (endsWith(seg[k], " ")) { pieces[[length(pieces) + 1L]] <- list(w = cur_w, b = cur_b); cur_w <- character(0); cur_b <- logical(0) } }
  pieces[[length(pieces) + 1L]] <- list(w = cur_w, b = cur_b)
  pieces <- Filter(function(p) length(p$w) > 0, pieces); if (!length(pieces)) return(1L)
  wlen <- vapply(pieces, function(p) sum(mapply(function(w, b) text_w(w, size, b), p$w, p$b)), 0)
  sp <- text_w(" ", size, bold); lines <- 1L; cur <- 0
  for (ww in wlen) {
    if (ww > width) { chunks <- ceiling(ww / width); lines <- lines + (cur > 0) + chunks - 1L; cur <- ww - (chunks - 1) * width; next }
    add <- if (cur == 0) ww else sp + ww
    if (cur + add > width) { lines <- lines + 1L; cur <- ww } else cur <- cur + add }
  as.integer(lines)
}
# 문단 목록의 높이(인치): 줄 높이 = 1.2 x 줄 간격 배수 x 글자 크기(Pretendard 렌더링 실측: 줄 간격 1.1에서 1.32배, 1.05에서 1.26배,
# 표 1.0에서 1.19배), 문단 간격 gap_pt, 상하 여백 0.10 in(글상자 기본 여백 0.05 in x 2). 3% 넘게 넘치면 실패(실제 넘침은 check_deck 9가 PDF로 확인)
# 채운 카드(배경색이 있고 높이 CARD_MIN_H 이상인 글상자)는 안쪽 여백이 좌우 0.12 in, 상하 0.08 in(pp_bullets가 bodyPr에 적는다)
CARD_MIN_H <- 0.6; CARD_INS <- c(lr = 0.12, tb = 0.08)
est_height <- function(paras, width, size, gap_pt = 6, indent = 0, bold = FALSE, line = 1.1, card = FALSE) {
  inner <- width - (if (card) 2 * CARD_INS[["lr"]] else 0.2) - indent
  n <- vapply(paras, function(p) est_lines(p, inner, size, bold), 1L)
  gaps <- if (length(gap_pt) > 1L) sum(head(gap_pt, -1L)) else max(0, length(paras) - 1) * gap_pt   # 문단별 간격(마지막 문단 뒤는 세지 않음)
  sum(n) * size * 1.2 * line / 72 + gaps / 72 + (if (card) 2 * CARD_INS[["tb"]] else 0.10)
}
fit_check <- function(label, paras, box, size, gap_pt = 6, indent = 0, bold = FALSE, line = 1.1, card = FALSE) {
  h <- est_height(paras, box[3], size, gap_pt, indent, bold, line, card)
  DK$fit[[length(DK$fit) + 1L]] <- data.table(slide = REG$sec, shape = label, est_h = h, box_h = box[4], ratio = h / box[4])
  if (h > box[4] * 1.03 && isTRUE(DK$strict)) stop(sprintf("%s %s: estimated text height %.2f in exceeds box %.2f in", REG$sec, label, h, box[4]), call. = FALSE)
  invisible(h)
}

# ---- 글자 조각 --------------------------------------------------------------------------------------------------------------------------
ftp <- function(size = SZ$body, color = PAL$ink, bold = FALSE, italic = FALSE)
  fp_text(font.size = size, color = color, bold = bold, italic = italic, font.family = FONT, eastasia.family = FONT, hansi.family = FONT, cs.family = FONT)
# "**굵게**", "__강조색__" 표기를 run으로 나눈다
runs <- function(s, size, color = PAL$ink, accent = PAL$blue, bold = FALSE) {
  parts <- regmatches(s, gregexpr("\\*\\*[^*]+\\*\\*|__[^_]+__|[^*_]+|[*_]", s))[[1]]
  lapply(parts, function(p) {
    if (grepl("^\\*\\*.*\\*\\*$", p)) ftext(substr(p, 3, nchar(p) - 2), ftp(size, color, TRUE))
    else if (grepl("^__.*__$", p)) ftext(substr(p, 3, nchar(p) - 2), ftp(size, accent, TRUE))
    else ftext(p, ftp(size, color, bold)) })
}
BUL <- function(level) sprintf("⁣B%d⁣", level)   # 글머리표 표지(pp_bullets가 내어쓰기 글머리표로 바꾼다)
# 표시층: AUC0-inf 같은 용어가 붙임표에서 줄바꿈되지 않도록 줄바꿈 없는 붙임표(U+2011)로 바꾼다(문구 파일은 그대로)
# 한 덩어리로 읽는 말이 줄 끝에서 갈라지지 않게 한다(LibreOffice는 숫자와 한글, 숫자와 단위 사이에서도 줄을 바꾼다):
#  숫자 뒤 공백 + 단위/모델, 'Day'·부등호 뒤 공백 → 줄바꿈 없는 공백(U+00A0); 숫자·%와 붙은 한글 사이 → 단어 결합자(U+2060)
# 음수 부호로 쓴 줄바꿈 없는 붙임표(U+2011; 일부 슬라이드 코드가 넣음)를 ASCII '-'로 바꾼다. 앞이 숫자인 붙임표(날짜 2026-09-25 등)는 그대로 둔다
neg_ascii <- function(s) gsub("(?<![0-9A-Za-z])\u2011(?=[0-9])", "-", s, perl = TRUE)
nobreak <- function(s) {
  s <- neg_ascii(s)                                                                    # 음수 부호는 ASCII '-'로 통일(아래에서 숫자와 붙인다)
  s <- gsub("AUC0-(inf|last|tlast)", "AUC0\u2011\\1", s, perl = TRUE)
  s <- gsub("([0-9%)]) (?=(kg|mg|mL|L/|mg/|mg·|day|h\\b|모델|명|칸|일|회|배|개|시간|점|mg/kg)(?![A-Za-z]))", "\\1\u00a0\u2060", s, perl = TRUE)   # LibreOffice는 NBSP 뒤 한글에서 줄을 바꾸므로 U+2060을 덧붙인다
  s <- gsub("(Day|≥|≤|>|<|×) (?=[0-9])", "\\1\u00a0\u2060", s, perl = TRUE)
  s <- gsub("(?<![A-Za-z])(SD|CI|GMR|CV|Wilson|span|R²) (?=[0-9(])", "\\1\u00a0", s, perl = TRUE)   # 통계량 이름과 값
  s <- gsub("([0-9]%) (?=CI)", "\\1\u00a0", s, perl = TRUE)                                       # '90% CI'
  s <- gsub("(?<=[0-9%A-Za-z)\u00b2])(?=[\uac00-\ud7a3])", "\u2060", s, perl = TRUE)   # 숫자·영문·닫는 괄호와 바로 붙은 한글(단위, 조사)
  s <- gsub("(?<=^|[\\s(~,;:/])-(?=[0-9])", "-\u2060", s, perl = TRUE)                   # 음수 부호와 숫자
  # 한글과 여는 괄호는 붙이지 않는다: LibreOffice는 괄호 뒤 결합자를 무시하고, 붙인 덩어리가 길면 한글 단어 가운데서 줄을 바꾼다(시험 렌더링 확인)
  s <- gsub("(세트|규칙|모델|기준|분석군) (?=[(A-D])", "\\1\u00a0\u2060", s, perl = TRUE)          # '세트 (i)', '규칙 A
  # 한 단위로 읽는 용어 쌍: 'Wilson 하한', '참 AUC0-inf', 'Km 0.01', '(iv) 65.9%'.
  # 한글 단어끼리의 쌍('창 포착률', '1종 오류')은 붙이지 않는다: LibreOffice는 한글 옆의 NBSP·U+2060·U+202F를 줄바꿈 금지로 보지 않고
  # 그 앞뒤에서 줄을 바꾸며, 줄 첫머리에 공백이 보인다(2026-09-29 렌더링 시험: 16칸 중 10~12칸 분리)
  s <- gsub("(Wilson) (?=하한)", "\\1\u00a0\u2060", s, perl = TRUE)
  s <- gsub("(참|모의|관측) (?=(AUC|Cmax))", "\\1\u00a0", s, perl = TRUE)
  s <- gsub("(?<![A-Za-z])(Km|Vmax|V2|ka|ke|Vc|CL|Q) (?=[×0-9])", "\\1\u00a0\u2060", s, perl = TRUE)
  s <- gsub("(\\((?:i|ii|iii|iv)\\)) (?=-?[0-9][0-9.,]*%)", "\\1\u00a0\u2060", s, perl = TRUE)   # 세트 표지와 뒤따르는 백분율만('세트 (i) 2016 모델'은 붙이지 않음)
  s
}
para <- function(s, size = SZ$body, color = PAL$ink, bold = FALSE, align = "left", bullet = NA, gap_pt = 6, line = 1.1, accent = PAL$blue) {
  r <- runs(nobreak(s), size, color, accent, bold)
  if (!is.na(bullet)) r <- c(list(ftext(BUL(bullet), ftp(size, color))), r)
  do.call(fpar, c(r, list(fp_p = fp_par(text.align = align, padding.bottom = gap_pt, line_spacing = line))))
}

# ---- 슬라이드 뼈대 -----------------------------------------------------------------------------------------------------------------------
loc <- function(box, label, bg = NULL, geom = NULL, ln = NULL) ph_location(left = box[1], top = box[2], width = box[3], height = box[4], newlabel = label, bg = bg, geom = geom, ln = ln)
no_line <- function() sp_line(color = "#ffffff", lwd = 0, lty = "solid")
deck_slide <- function(id, tag = c("sim", "lit", "litsim", "none"), dark = FALSE) {
  tag <- match.arg(tag); sec(id); DK$n <- DK$n + 1L
  DK$x <- add_slide(DK$x, layout = "Blank", master = "Office Theme")
  DK$cur <- list(id = id, n = DK$n, dark = dark, bullets = 0L, shapes = 0L, title = NA_character_, table_rows = 0L, body_lines = 0L, visual = 0, title_lines = NA_integer_)
  if (dark) DK$x <- ph_with(DK$x, fpar(ftext(" ", ftp(8))), location = loc(c(0, 0, GEO$W, GEO$H), "background", bg = PAL$dark, geom = "rect", ln = no_line()))
  if (tag != "none") {
    lab <- DK$txt$common$tags[[tag]]
    b <- c(GEO$W - GEO$MR - GEO$TAG_W, GEO$KICK_TOP - 0.08, GEO$TAG_W, GEO$TAG_H - 0.04)   # 제목 첫 줄과 겹치지 않게 위로
    DK$x <- ph_with(DK$x, fpar(ftext(lab, ftp(SZ$tag, if (dark) PAL$dark_ink else PAL$ink2, TRUE)), fp_p = fp_par(text.align = "center")),
                    location = loc(b, "tag", bg = if (dark) "#2a3644" else PAL$tint_grey, geom = "roundRect", ln = no_line()))
  }
  invisible(NULL)
}
deck_kicker <- function(s) {
  b <- c(GEO$ML, GEO$KICK_TOP, GEO$CW - GEO$TAG_W - 0.2, 0.34)
  DK$x <- ph_with(DK$x, para(s, SZ$kicker, if (DK$cur$dark) "#8fb8ea" else PAL$blue, bold = TRUE, gap_pt = 0), location = loc(b, "kicker")); invisible(NULL)
}
deck_title <- function(s, size = SZ$title, box = c(GEO$ML, GEO$TITLE_TOP, GEO$CW - 0.1, GEO$TITLE_H)) {
  DK$cur$title <- strip_markup(s); fit_check("title", s, box, size, bold = TRUE, line = 1.05)
  DK$cur$title_lines <- est_lines(s, box[3] - 0.2, size, TRUE)
  if (!is.null(LIMITS$title_lines) && DK$cur$title_lines > LIMITS$title_lines && isTRUE(DK$strict))
    stop(sprintf("%s: title needs %d lines (limit %d)", REG$sec, DK$cur$title_lines, LIMITS$title_lines), call. = FALSE)
  DK$x <- ph_with(DK$x, para(s, size, if (DK$cur$dark) PAL$dark_ink else PAL$ink, bold = TRUE, gap_pt = 0, line = 1.05), location = loc(box, "title")); invisible(NULL)
}
# items: 문자 벡터. 앞에 "- "가 붙은 항목은 2단계(요점 수에 세지 않음)
deck_bullets <- function(items, box, size = SZ$body, gap_pt = 8, label = "body", color = NULL) {
  stopifnot(size >= SZ$body_min); lvl <- ifelse(startsWith(items, "- "), 1L, 0L); items <- sub("^- ", "", items)
  n0 <- sum(lvl == 0L); DK$cur$bullets <- DK$cur$bullets + n0
  if (DK$cur$bullets > LIMITS$bullets) stop(sprintf("%s: %d bullet points (limit %d)", REG$sec, DK$cur$bullets, LIMITS$bullets), call. = FALSE)
  # 하위 글머리표는 윗 요점에 붙인다: 하위 항목 앞 문단과 하위 항목끼리는 간격 4 pt(마지막 하위 항목 뒤는 원래 간격)
  gp <- rep(gap_pt, length(items)); nx <- c(lvl[-1L], 0L); gp[nx == 1L] <- min(gap_pt, 4)
  fit_check(label, items, box, size, gp, indent = 0.3)
  if (startsWith(label, "body")) deck_count_lines(items, box[3] - 0.2 - 0.3, size)
  col <- color %||% (if (DK$cur$dark) PAL$dark_ink else PAL$ink)
  ps <- lapply(seq_along(items), function(i) para(items[i], if (lvl[i] == 1L) max(SZ$body_min, size - 2) else size, col, bullet = lvl[i], gap_pt = gp[i]))
  DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(box, label)); invisible(NULL)
}
deck_text <- function(s, box, size = SZ$body, bold = FALSE, color = NULL, align = "left", label = "text", bg = NULL, gap_pt = 6, geom = NULL) {
  stopifnot(size >= SZ$body_min || grepl("^(caption|label|axis|kicker|foot)", label))
  if (!is.null(SZ$caption_min) && size < SZ$caption_min) stop(sprintf("%s %s: %gpt below the %gpt minimum", REG$sec, label, size, SZ$caption_min), call. = FALSE)
  fit_check(label, s, box, size, gap_pt, bold = bold, card = !is.null(bg) && box[4] >= CARD_MIN_H)
  if (startsWith(label, "body")) deck_count_lines(s, box[3] - (if (!is.null(bg) && box[4] >= CARD_MIN_H) 2 * CARD_INS[["lr"]] else 0.2), size, bold)
  col <- color %||% (if (DK$cur$dark) PAL$dark_ink else PAL$ink)
  ps <- lapply(s, function(z) para(z, size, col, bold, align, gap_pt = gap_pt))
  DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(box, label, bg = bg, geom = geom, ln = if (!is.null(bg)) no_line() else NULL)); invisible(NULL)
}
# 큰 수치 한 개와 설명(수치는 반드시 d* 함수 결과)
deck_stat <- function(value, label, box, color = PAL$blue, bg = PAL$tint_blue, label_size = SZ$stat_label, value_size = SZ$stat) {
  fit_check("stat_label", label, c(box[1], box[2], box[3], box[4] - value_size * 1.25 / 72), label_size, gap_pt = 0, card = box[4] >= CARD_MIN_H)
  ps <- list(para(value, value_size, color, TRUE, "left", gap_pt = 2, line = 0.9), para(label, label_size, PAL$ink, FALSE, "left", gap_pt = 0))   # 큰 수치 줄은 줄 간격 0.9(글자 위 여백이 커 가운데 정렬해도 위가 넓어 보인다)
  DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(box, "stat", bg = bg, geom = "roundRect", ln = no_line())); invisible(NULL)
}
deck_box <- function(box, fill = PAL$tint_grey, geom = "roundRect", label = "shape", line_col = NULL) {
  DK$x <- ph_with(DK$x, fpar(ftext(" ", ftp(SZ$body_min))), location = loc(box, label, bg = fill, geom = geom, ln = if (is.null(line_col)) no_line() else sp_line(color = line_col, lwd = 1.5)))
  invisible(NULL)
}
# 표: df는 문자열 data.frame(수치는 d* 결과). 본문 6행 이하.
# 표 높이 추정(인치): 셀마다 줄 수(Pretendard 폭) x 줄 높이 + 행마다 위아래 여백(pad pt x 2). 셀 여백 0.14 in(실제 넘침은 검사 9가 PDF로 확인)
deck_table_h <- function(df, width, widths = NULL, size = SZ$table, bold_col1 = TRUE, pad = 4) {
  df <- as.data.frame(df); for (j in seq_along(df)) df[[j]] <- nobreak(as.character(df[[j]])); names(df) <- nobreak(names(df))
  if (is.null(widths)) widths <- rep(width / ncol(df), ncol(df)) else widths <- widths / sum(widths) * width
  nl <- function(v, w, b = FALSE) vapply(as.character(v), function(s) est_lines(s, w - 0.14, size, b), 1L)
  hdr <- max(mapply(function(v, w) max(nl(v, w, TRUE)), names(df), widths))
  bod <- apply(matrix(sapply(seq_along(widths), function(j) nl(df[[j]], widths[j], bold_col1 && j == 1)), nrow = nrow(df)), 1, max)
  (hdr + sum(bod)) * size * 1.2 / 72 + (nrow(df) + 1) * 2 * pad / 72
}
# pad: 셀 위아래 여백(pt). group_end: 주면 행 사이 가는 선 대신 그 행들 아래에만 굵은 선(묶음 구분)
deck_table <- function(df, box, widths = NULL, size = SZ$table, header_fill = PAL$tint_blue, bold_col1 = TRUE, align_num = TRUE, highlight = NULL, highlight_fill = PAL$tint_orange, label = "table", align_cols = NULL,
                       pad = 4, group_end = NULL) {
  stopifnot(size >= SZ$table_min)
  if (nrow(df) > LIMITS$table_rows) stop(sprintf("%s: table with %d body rows (limit %d)", REG$sec, nrow(df), LIMITS$table_rows), call. = FALSE)
  DK$cur$table_rows <- max(DK$cur$table_rows, nrow(df))
  df <- as.data.frame(df); for (j in seq_along(df)) df[[j]] <- nobreak(as.character(df[[j]])); names(df) <- nobreak(names(df))
  ft <- flextable(df)
  ft <- font(ft, fontname = FONT, part = "all", eastasia.family = FONT, hansi.family = FONT, cs.family = FONT)
  ft <- fontsize(ft, size = size, part = "all"); ft <- color(ft, color = PAL$ink, part = "all")
  ft <- bold(ft, part = "header"); ft <- bg(ft, bg = header_fill, part = "header")
  if (bold_col1) ft <- bold(ft, j = 1, part = "body")
  if (!is.null(highlight)) ft <- bg(ft, i = highlight, bg = highlight_fill, part = "body")
  ft <- border_remove(ft)
  ft <- if (is.null(group_end)) hline(ft, border = fp_border_default(color = PAL$grid, width = 0.75), part = "body") else hline(ft, i = group_end, border = fp_border_default(color = PAL$ink2, width = 1.5), part = "body")
  ft <- hline_bottom(ft, border = fp_border_default(color = PAL$ink2, width = 1), part = "header"); ft <- hline_top(ft, border = fp_border_default(color = PAL$ink2, width = 1), part = "header")
  ft <- padding(ft, padding.top = pad, padding.bottom = pad, padding.left = 5, padding.right = 5, part = "all")
  if (align_num && ncol(df) > 1) ft <- align(ft, j = 2:ncol(df), align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "all"); ft <- valign(ft, valign = "center", part = "all")
  if (!is.null(align_cols)) { stopifnot(length(align_cols) == ncol(df)); for (j in seq_along(align_cols)) ft <- align(ft, j = j, align = align_cols[j], part = "all") }   # 열별 정렬
  if (is.null(widths)) widths <- rep(box[3] / ncol(df), ncol(df)) else widths <- widths / sum(widths) * box[3]
  ft <- width(ft, width = widths)
  h_est <- deck_table_h(df, box[3], widths, size, bold_col1, pad)
  DK$fit[[length(DK$fit) + 1L]] <- data.table(slide = REG$sec, shape = label, est_h = h_est, box_h = box[4], ratio = h_est / box[4])
  if (h_est > box[4] * 1.03 && isTRUE(DK$strict)) stop(sprintf("%s %s: estimated table height %.2f in exceeds box %.2f in", REG$sec, label, h_est, box[4]), call. = FALSE)
  DK$x <- ph_with(DK$x, ft, location = loc(box, label)); invisible(NULL)
}
# 그림: 슬라이드 자리 크기 그대로 결과 파일에서 다시 그린다(늘리거나 줄이지 않음). src: 그림 자료 파일(추적 기록)
# 본문 줄 수(core: 본문 슬라이드 3줄 이내). 이름이 body로 시작하는 글상자만 센다(카드·캡션·그림 글자 제외)
deck_count_lines <- function(paras, width, size, bold = FALSE) {
  DK$cur$body_lines <- DK$cur$body_lines + sum(vapply(strip_markup(paras), est_lines, 1L, width = width, size = size, bold = bold)); invisible(NULL) }
# 주 시각 요소(그림, 도식 묶음, 카드 묶음)의 면적 등록(core: 본문 영역의 60% 이상)
deck_visual <- function(box) { DK$cur$visual <- max(DK$cur$visual, box[3] * box[4]); invisible(NULL) }
# 그림 안 글자 크기 하한(core: 14pt). 테마의 모든 글자 요소와 글자 레이어(geom_text/label, annotate)를 검사한다. patchwork는 하위 그림까지
gg_leaves <- function(p) {
  if (inherits(p, "patchwork")) { out <- list(); for (q in p$patches$plots) out <- c(out, gg_leaves(q))
    q0 <- p; q0$patches <- NULL; class(q0) <- setdiff(class(q0), "patchwork")
    if (length(q0$layers)) out <- c(out, list(q0)); out }                # wrap_plots의 바탕(빈 ggplot, 그리는 글자 없음)은 건너뛴다
  else if (inherits(p, "ggplot")) list(p) else list()
}
fig_text_check <- function(p, name, min_pt = SZ$fig_min) {
  if (is.null(min_pt)) return(invisible(NULL))
  EL <- c("axis.text.x.bottom", "axis.text.x.top", "axis.text.y.left", "axis.text.y.right", "axis.title.x.bottom", "axis.title.x.top", "axis.title.y.left", "axis.title.y.right",
          "legend.text", "legend.title", "strip.text.x.top", "strip.text.y.right", "plot.title", "plot.subtitle", "plot.caption", "plot.tag")
  bad <- character(0)
  for (q in gg_leaves(p)) {
    th <- ggplot2:::plot_theme(q)
    for (el in EL) { e <- ggplot2::calc_element(el, th); if (inherits(e, "element_text") && !is.null(e$size) && e$size < min_pt - 1e-6) bad <- c(bad, sprintf("%s %.1fpt", el, e$size)) }
    b <- ggplot2::ggplot_build(q)
    for (i in seq_along(q$layers)) if (inherits(q$layers[[i]]$geom, c("GeomText", "GeomLabel"))) {
      sz <- b$data[[i]]$size * ggplot2::.pt; if (length(sz) && any(sz < min_pt - 1e-6)) bad <- c(bad, sprintf("text layer %d %.1fpt", i, min(sz))) }
  }
  if (length(bad)) stop(sprintf("%s figure %s: text below %gpt: %s", REG$sec, name, min_pt, paste(unique(bad), collapse = ", ")), call. = FALSE)
  invisible(NULL)
}
deck_figure <- function(p, name, box, src, dpi = 220, visual = TRUE) {
  dsrc(sprintf("figure %s", name), src); fig_text_check(p, name); if (visual) deck_visual(box)
  f <- file.path(DK$fig_dir, sprintf("%s_%s.png", DK$lang, name))
  ragg::agg_png(f, width = box[3], height = box[4], units = "in", res = dpi, background = "white"); print(p); grDevices::dev.off()
  DK$x <- ph_with(DK$x, external_img(f, width = box[3], height = box[4]), location = loc(box, sprintf("figure_%s", name)), use_loc_size = TRUE); invisible(f)
}
deck_notes <- function(s) { s <- neg_ascii(s); DK$x <- set_notes(DK$x, value = paste(s, collapse = "\n"), location = notes_location_type("body")); DK$cur$notes <- paste(s, collapse = "\n"); invisible(NULL) }
# 슬라이드 끝: 바닥글(결과 파일 ID, 버전, 커밋, 쪽 번호)과 메타 기록
# 바닥글 출처 목록에서 맨 앞에 둘 파일(슬라이드의 핵심 출처가 설정 파일이라 '외 N개'에 가려질 때)
deck_src_first <- function(files) { DK$cur$src_first <- files; invisible(NULL) }
deck_end <- function() {
  tr <- rbindlist(REG$trace, fill = TRUE); src <- if (nrow(tr)) unique(tr[section == REG$sec, source_file]) else character(0)
  ids <- unique(c(basename(src[!grepl("^config/", src)]), basename(src[grepl("^config/", src)])))   # 결과 파일 먼저, 설정 파일은 뒤
  if (length(DK$cur$src_first)) ids <- unique(c(intersect(basename(DK$cur$src_first), ids), ids))   # 슬라이드가 지정한 핵심 출처를 맨 앞에(deck_src_first)
  lab <- DK$txt$common$footer
  # 출처 목록은 한 줄에 들어가는 만큼만 적고 나머지는 "외 N개"로 줄인다(전체 목록은 추적표에 있음)
  rw <- if (identical(DK$prof$footer, "version_source")) 1.2 else 3.0          # 오른쪽(버전) 칸 폭
  wmax <- GEO$CW - rw - 0.1 - 0.25
  mk <- function(k) paste0(lab$source, " ", if (k < length(ids)) sprintf("%s %s", paste(ids[seq_len(k)], collapse = ", "), fill(lab$more, list(n = length(ids) - k))) else paste(ids, collapse = ", "))
  k <- min(length(ids), 5L); while (k > 1L && text_w(mk(k), SZ$foot) > wmax) k <- k - 1L
  left <- if (length(ids)) mk(k) else ""
  right <- if (identical(DK$prof$footer, "version_source")) sprintf("%s%s", DK$version, if (DK$clean) "" else "*") else sprintf("%s  |  %s %s%s  |  %d", DK$version, lab$commit, DK$commit, if (DK$clean) "" else "*", DK$cur$n)
  col <- if (DK$cur$dark) "#9aa6b2" else PAL$muted
  DK$x <- ph_with(DK$x, fpar(ftext(left, ftp(SZ$foot, col))), location = loc(c(GEO$ML, GEO$FOOT_TOP, GEO$CW - rw - 0.1, GEO$FOOT_H), "footer_src"))
  DK$x <- ph_with(DK$x, fpar(ftext(right, ftp(SZ$foot, col)), fp_p = fp_par(text.align = "right")), location = loc(c(GEO$W - GEO$MR - rw, GEO$FOOT_TOP, rw, GEO$FOOT_H), "footer_ver"))
  m <- data.table(n = DK$cur$n, id = DK$cur$id, title = DK$cur$title, bullets = DK$cur$bullets, table_rows = DK$cur$table_rows,
                  sources = paste(src, collapse = "; "), notes = DK$cur$notes %||% "")
  if (!is.null(LIMITS$body_lines)) {                                   # core: 본문 슬라이드의 본문 줄 수, 주 시각 요소 면적 비, 제목 줄 수
    body_area <- GEO$CW * (GEO$BODY_BOTTOM - GEO$BODY_TOP); share <- DK$cur$visual / body_area
    main <- DK$cur$id %in% unlist(DK$txt$common$main_limits)
    m[, `:=`(body_lines = DK$cur$body_lines, fig_share = round(share, 3), title_lines = DK$cur$title_lines, main_limits = main)]
    if (main && isTRUE(DK$strict)) {
      if (DK$cur$body_lines > LIMITS$body_lines) stop(sprintf("%s: %d body lines (limit %d)", REG$sec, DK$cur$body_lines, LIMITS$body_lines), call. = FALSE)
      if (share < LIMITS$fig_share) stop(sprintf("%s: main visual covers %.0f%% of the body area (minimum %.0f%%)", REG$sec, 100 * share, 100 * LIMITS$fig_share), call. = FALSE)
    }
  }
  DK$meta[[length(DK$meta) + 1L]] <- m
  invisible(NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 그림 테마(보고서와 같은 검증 팔레트; 색 외에 선 모양·표식으로도 구분) -------------------------------------------------------------
theme_deck <- function(base = 14) {
  theme_minimal(base_size = base, base_family = FONT) + theme(
    plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = PAL$grid, linewidth = 0.35), axis.text = element_text(colour = PAL$ink2, size = base - 1),
    axis.title = element_text(colour = PAL$ink2, size = base), legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(colour = PAL$ink2, size = base - 1), strip.text = element_text(face = "bold", colour = PAL$ink, size = base),
    plot.margin = margin(6, 12, 6, 6))
}
# 핵심 덱(core) 그림 테마: 모든 글자 14pt 이상(fig_text_check가 확인). 색은 두 가지 뜻으로만 쓴다: AUClast 파랑, 외삽·AUCinf 주황.
# 모델 구분은 색이 아니라 표식 모양·선 모양·명도(검정 계열)로 한다(CORE_MODEL_*).
theme_core <- function(base = 16) {
  theme_deck(base) + theme(axis.text = element_text(colour = PAL$ink2, size = base - 2), axis.title = element_text(colour = PAL$ink2, size = base - 1),
                           legend.text = element_text(colour = PAL$ink2, size = base - 2), strip.text = element_text(face = "bold", colour = PAL$ink, size = base, hjust = 0),
                           plot.title = element_text(size = base, face = "bold", colour = PAL$ink), plot.subtitle = element_text(size = base - 1, colour = PAL$ink2),
                           plot.caption = element_text(size = base - 2, colour = PAL$ink2, hjust = 0), plot.tag = element_text(size = base - 2, colour = PAL$ink2))
}
PT <- function(pt) pt / ggplot2::.pt                                      # 그림 글자 크기(pt) → ggplot size(mm)
CORE_MODEL_COL <- c(k2016 = PAL$ink, k2020 = "#7a7974"); CORE_MODEL_SHAPE <- c(k2016 = 16, k2020 = 17); CORE_MODEL_LT <- c(k2016 = "solid", k2020 = "22")
ORANGE_LIGHT <- "#f5b99c"; BLUE_LIGHT <- "#9cc3ee"
MODEL_COL <- c(k2016 = PAL$blue, k2020 = PAL$orange); MODEL_SHAPE <- c(k2016 = 16, k2020 = 17); MODEL_LT <- c(k2016 = "solid", k2020 = "22")
model_lab <- function() { m <- DK$txt$common$models; c(k2016 = m$k2016, k2020 = m$k2020) }

# ---- 저장과 후처리 -------------------------------------------------------------------------------------------------------------------------
# 저장 후 처리: 글머리표 표지를 PowerPoint 글머리표(내어쓰기)로 바꾸고, 자리표시자 표지 제거, 한국어 어절 줄바꿈, 모서리 반지름 통일
pp_bullets <- function(path) {
  tmp <- tempfile("pp_"); dir.create(tmp); utils::unzip(path, exdir = tmp)
  sl <- list.files(file.path(tmp, "ppt", "slides"), pattern = "^slide[0-9]+\\.xml$", full.names = TRUE)
  ns <- c(a = "http://schemas.openxmlformats.org/drawingml/2006/main", p = "http://schemas.openxmlformats.org/presentationml/2006/main")
  nfix <- 0L
  for (f in sl) {
    doc <- xml2::read_xml(f); ts <- xml2::xml_find_all(doc, ".//a:t", ns)
    for (t in ts) {
      s <- xml2::xml_text(t); m <- regmatches(s, regexec("^⁣B([0-9])⁣", s))[[1]]
      if (!length(m)) next
      lvl <- as.integer(m[2]); xml2::xml_text(t) <- sub("^⁣B[0-9]⁣", "", s)
      r <- xml2::xml_parent(t); p <- xml2::xml_parent(r)
      if (!nzchar(xml2::xml_text(t))) xml2::xml_remove(r)
      ppr <- xml2::xml_find_first(p, "./a:pPr", ns)
      if (inherits(ppr, "xml_missing")) { xml2::xml_add_child(p, "a:pPr", .where = 0); ppr <- xml2::xml_find_first(p, "./a:pPr", ns) }
      marL <- if (lvl == 0) 285750L else 571500L; ind <- -228600L
      xml2::xml_set_attr(ppr, "marL", marL); xml2::xml_set_attr(ppr, "indent", ind)
      for (old in xml2::xml_find_all(ppr, "./a:buNone|./a:buChar|./a:buFont|./a:buClr", ns)) xml2::xml_remove(old)
      bc <- xml2::xml_add_child(ppr, "a:buClr"); xml2::xml_add_child(bc, "a:srgbClr", val = if (lvl == 0) "2A78D6" else "8A8984")
      xml2::xml_add_child(ppr, "a:buFont", typeface = "Arial"); xml2::xml_add_child(ppr, "a:buChar", char = if (lvl == 0) "●" else "◦")
      nfix <- nfix + 1L
    }
    # 자리표시자 표지(p:ph)를 지워 보통 도형으로 만든다: LibreOffice(PDF)는 자리표시자의 도형 모양(roundRect, rightArrow)을 무시한다.
    # 위치·채우기·글자 서식은 도형에 모두 명시되어 있어 모양만 바뀐다.
    for (ph in xml2::xml_find_all(doc, ".//p:nvPr/p:ph", ns)) xml2::xml_remove(ph)
    # 모양이 없는 도형(자리표시자에서 사각형을 물려받던 것)에는 사각형을 명시한다(없으면 채우기가 그려지지 않는다)
    for (sppr in xml2::xml_find_all(doc, ".//p:sp/p:spPr[not(a:prstGeom) and not(a:custGeom)]", ns)) {
      xf <- xml2::xml_find_first(sppr, "./a:xfrm", ns); if (inherits(xf, "xml_missing")) next
      g <- xml2::xml_add_sibling(xf, "a:prstGeom", prst = "rect", .where = "after"); xml2::xml_add_child(g, "a:avLst")
    }
    # 한국어는 어절 단위로 줄을 바꾼다(PowerPoint의 '한글 단어 잘림 허용' 끔)
    for (ppr in xml2::xml_find_all(doc, ".//a:pPr", ns)) xml2::xml_set_attr(ppr, "eaLnBrk", "0")
    # 채운 카드의 안쪽 여백(좌우 0.12 in, 상하 0.08 in): 글이 카드 가장자리에 붙지 않게(높이 CARD_MIN_H 이상, 글자가 있는 도형)
    for (sp in xml2::xml_find_all(doc, ".//p:sp[p:spPr/a:solidFill and p:txBody]", ns)) {
      nm <- xml2::xml_attr(xml2::xml_find_first(sp, ".//p:cNvPr", ns), "name"); if (nm %in% c("background", "tag")) next
      ext <- xml2::xml_find_first(sp, "./p:spPr/a:xfrm/a:ext", ns); if (inherits(ext, "xml_missing") || as.numeric(xml2::xml_attr(ext, "cy")) < CARD_MIN_H * 914400) next
      if (!nzchar(trimws(xml2::xml_text(xml2::xml_find_first(sp, "./p:txBody", ns))))) next
      bp <- xml2::xml_find_first(sp, "./p:txBody/a:bodyPr", ns)
      xml2::xml_set_attrs(bp, c(xml2::xml_attrs(bp), lIns = as.character(round(CARD_INS[["lr"]] * 914400)), rIns = as.character(round(CARD_INS[["lr"]] * 914400)),
                                tIns = as.character(round(CARD_INS[["tb"]] * 914400)), bIns = as.character(round(CARD_INS[["tb"]] * 914400))))
    }
    # 큰 수치 카드(stat)와 낮은 채운 도형(칩·띠, 높이 CARD_MIN_H 미만)은 글을 세로 가운데에 둔다(위아래 여백을 같게; 위 정렬이면 아래 가장자리에 붙는다)
    for (sp in xml2::xml_find_all(doc, ".//p:sp[p:spPr/a:solidFill and p:txBody]", ns)) {
      nm <- xml2::xml_attr(xml2::xml_find_first(sp, ".//p:cNvPr", ns), "name"); if (nm %in% c("background", "tag")) next
      ext <- xml2::xml_find_first(sp, "./p:spPr/a:xfrm/a:ext", ns); if (inherits(ext, "xml_missing")) next
      if (!(nm == "stat" || startsWith(nm, "cmid") || as.numeric(xml2::xml_attr(ext, "cy")) < CARD_MIN_H * 914400)) next   # cmid*: 핵심 덱의 세로 가운데 상자
      if (!nzchar(trimws(xml2::xml_text(xml2::xml_find_first(sp, "./p:txBody", ns))))) next
      xml2::xml_set_attr(xml2::xml_find_first(sp, "./p:txBody/a:bodyPr", ns), "anchor", "ctr")
    }
    # 둥근 사각형의 모서리 반지름을 0.1 in로 통일(기본값은 짧은 변의 16.7%라 큰 카드가 지나치게 둥글다)
    for (sp in xml2::xml_find_all(doc, ".//p:sp[p:spPr/a:prstGeom[@prst='roundRect']]", ns)) {
      ext <- xml2::xml_find_first(sp, "./p:spPr/a:xfrm/a:ext", ns); if (inherits(ext, "xml_missing")) next
      m <- min(as.numeric(xml2::xml_attr(ext, "cx")), as.numeric(xml2::xml_attr(ext, "cy"))); if (!is.finite(m) || m <= 0) next
      av <- xml2::xml_find_first(sp, "./p:spPr/a:prstGeom/a:avLst", ns)
      for (g in xml2::xml_children(av)) xml2::xml_remove(g)
      xml2::xml_add_child(av, "a:gd", name = "adj", fmla = sprintf("val %d", as.integer(round(min(50000, 91440 / m * 100000)))))
    }
    xml2::write_xml(doc, f)
  }
  file.remove(path); old <- setwd(tmp); on.exit(setwd(old)); zip::zip(path, files = list.files(".", recursive = TRUE, all.files = TRUE), mode = "mirror"); setwd(old)
  invisible(nfix)
}
deck_write <- function(path, trace_path = NULL, meta_path = NULL) {
  print(DK$x, target = path); n <- pp_bullets(path)
  if (!is.null(trace_path)) {
    tr <- rbindlist(REG$trace, fill = TRUE); files <- unique(tr$source_file)
    sh <- vapply(files, function(f) digest::digest(file = proj_path(f), algo = "sha256"), ""); tr[, source_sha256 := sh[source_file]]
    tr[, slide_n := match(section, vapply(DK$meta, function(m) m$id, ""))]; tr[, source_commit := DK$commit]
    fwrite(tr, trace_path)
  }
  if (!is.null(meta_path)) { m <- rbindlist(DK$meta); fwrite(m, meta_path); fwrite(rbindlist(DK$fit), sub("\\.csv$", "_fit.csv", meta_path))
    fwrite(rbindlist(lapply(DK$headline, as.data.table)), sub("\\.csv$", "_headline.csv", meta_path)) }
  invisible(n)
}
