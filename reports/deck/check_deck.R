#!/usr/bin/env Rscript
# check_deck.R: 결과보고 슬라이드 자동 검사. 하나라도 실패하면 0이 아닌 코드로 끝난다.
#  1 숫자 추적: 슬라이드(본문·표·발표자 노트)의 모든 숫자가 그 슬라이드의 추적 행(deck_traceability.csv)에 있다(허용 문맥 제외).
#  2 빈 값: 자리표시자, NA/NaN/Inf, TODO가 없다.        3 대시: en-dash, em-dash, U+2212가 없고 음수 부호는 ASCII '-'다(글자와 슬라이드·노트 XML 전체, 글머리표 문자 포함).
#  4 아토피: 부록 A5 밖(본문 슬라이드)에 아토피 모집단 말이나 그 결과 파일 출처가 없다.
#  5 약어: 첫 등장 슬라이드에 풀이가 있다.               6 형식: 본문 16pt, 표 12pt 이상; 요점 5개, 표 본문 6행 이하.
#  7 추적 완결: 추적 행의 출처 파일이 있고 SHA-256이 현재 파일과 같다.
#  8 대조: 같은 파일·같은 행 조건·같은 열의 M&S 보고서 추적 행과 원값이 같고, 소수 자릿수가 같으면 인쇄값도 같다.
#          key_numbers_en.md에서 같은 파일을 인용한 줄의 숫자와 비교한다. 핵심 수치(headline)는 보고서나 key numbers에서 반드시 확인되어야 한다.
#  9 렌더링 배치: render_deck.py로 만든 PDF의 글줄 위치를 pptx 도형 상자와 대조한다(check_render.py: 넘침, 겹침, 가장자리 여백).
# 사용법: Rscript reports/deck/check_deck.R [pptx] [--lang ko] [--pdf 경로(기본: pptx와 같은 이름)]
source("R/00_setup.R"); suppressPackageStartupMessages({ library(data.table); library(xml2) })
args <- commandArgs(trailingOnly = TRUE); lang <- if ("--lang" %in% args) args[match("--lang", args) + 1] else "ko"
pptx <- if (length(args) && !startsWith(args[1], "--")) args[1] else proj_path("reports", "deck", sprintf("dupilumab_endpoint_results_v1.0.1%s.pptx", if (lang == "ko") "" else paste0("_", lang)))
base <- sub("\\.pptx$", "", pptx)
tr_f <- if (file.exists(paste0(base, "_traceability.csv"))) paste0(base, "_traceability.csv") else proj_path("reports", "deck", if (lang == "ko") "deck_traceability.csv" else sprintf("deck_traceability_%s.csv", lang))
me_f <- if (file.exists(paste0(base, "_meta.csv"))) paste0(base, "_meta.csv") else proj_path("reports", "deck", if (lang == "ko") "deck_meta.csv" else sprintf("deck_meta_%s.csv", lang))
tr <- fread(tr_f, encoding = "UTF-8", colClasses = list(character = c("value_raw", "value_printed", "locator"))); meta <- fread(me_f, encoding = "UTF-8")
hl_f <- sub("\\.csv$", "_headline.csv", me_f); hl <- if (file.exists(hl_f) && file.size(hl_f) > 0) fread(hl_f, encoding = "UTF-8", colClasses = list(character = "printed")) else data.table()
txt <- list(); for (f in sort(list.files(proj_path("reports", "deck", "text", lang), pattern = "\\.ya?ml$", full.names = TRUE))) txt <- c(txt, yaml::read_yaml(f))
CK <- txt$check
res <- list(); add <- function(check, slide, status, detail = "") res[[length(res) + 1L]] <<- data.table(check = check, slide = slide, status = status, detail = detail)

# ---- pptx 글자 추출(슬라이드 순서, 도형 이름별) -------------------------------------------------------------------------------------------
tmp <- tempfile("ck_"); dir.create(tmp); utils::unzip(pptx, exdir = tmp)
ns <- c(a = "http://schemas.openxmlformats.org/drawingml/2006/main", p = "http://schemas.openxmlformats.org/presentationml/2006/main",
        r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
pres <- read_xml(file.path(tmp, "ppt", "presentation.xml")); prel <- read_xml(file.path(tmp, "ppt", "_rels", "presentation.xml.rels"))
rid <- xml_attr(xml_find_all(pres, ".//p:sldIdLst/p:sldId", ns), "r:id", ns)
relmap <- setNames(xml_attr(xml_children(prel), "Target"), xml_attr(xml_children(prel), "Id"))
slide_files <- file.path(tmp, "ppt", relmap[rid])
stopifnot(length(slide_files) == nrow(meta))
para_text <- function(node) gsub("\u00a0", " ", gsub("\u2060", "", paste(vapply(xml_find_all(node, ".//a:p", ns), function(p) paste(xml_text(xml_find_all(p, ".//a:t", ns)), collapse = ""), ""), collapse = "\n")))   # 표시용 결합 문자 제거
S <- rbindlist(lapply(seq_along(slide_files), function(i) {
  d <- read_xml(slide_files[i]); sh <- xml_find_all(d, "//p:spTree/p:sp | //p:spTree/p:graphicFrame", ns)
  rows_ <- lapply(sh, function(s) {
    nm <- xml_attr(xml_find_first(s, ".//p:cNvPr", ns), "name"); kind <- xml_name(s)
    szs <- as.integer(xml_attr(xml_find_all(s, ".//a:rPr[@sz] | .//a:endParaRPr[@sz]", ns), "sz"))
    off <- xml_find_first(s, "./p:spPr/a:xfrm/a:off | ./p:xfrm/a:off", ns)
    data.table(n = i, shape = nm, kind = kind, text = para_text(s), min_sz = if (length(szs)) min(szs) else NA_integer_,
               y = if (inherits(off, "xml_missing")) NA_real_ else as.numeric(xml_attr(off, "y")), x = if (inherits(off, "xml_missing")) NA_real_ else as.numeric(xml_attr(off, "x")))
  })
  out <- rbindlist(rows_)
  # 발표자 노트
  rf <- file.path(dirname(slide_files[i]), "_rels", paste0(basename(slide_files[i]), ".rels"))
  if (file.exists(rf)) { rr <- read_xml(rf); tg <- xml_attr(xml_children(rr), "Target"); nt <- tg[grepl("notesSlide", tg)]
    if (length(nt)) { nd <- read_xml(normalizePath(file.path(dirname(slide_files[i]), nt))); body <- xml_find_all(nd, "//p:sp[.//p:ph[@type='body']]", ns)
      out <- rbind(out, data.table(n = i, shape = "notes", kind = "notes", text = paste(vapply(body, para_text, ""), collapse = "\n"), min_sz = NA_integer_, y = NA_real_, x = NA_real_)) } }
  out
}))
meta <- meta[order(n)]; S[, id := meta$id[n]]   # n = 슬라이드의 pptx 안 순서(부분 빌드에서도 meta와 같은 순서)
if (!nrow(S) || anyNA(S$id)) stop("slide text extraction failed")
for (sid in meta$id) if (!nrow(S[id == sid & kind != "notes" & nzchar(text)])) stop("no text extracted for ", sid)
vis <- S[!grepl("^footer_|^tag$|^background$", shape) & kind != "notes"][order(n, round(y / 45720), x)]   # 읽는 순서(위에서 아래, 0.05 in 단위로 같은 줄이면 왼쪽부터)
num_tokens <- function(s) {
  s <- gsub("‑", "-", s); for (rx in unlist(CK$number_contexts)) s <- gsub(rx, " ", s, perl = TRUE)
  m <- regmatches(s, gregexpr("(?<![A-Za-z0-9_.‑-])-?[0-9][0-9,]*(\\.[0-9]+)?|(?<=[~\\s(])-[0-9][0-9,]*(\\.[0-9]+)?", s, perl = TRUE))[[1]]
  unique(gsub(",", "", m))
}
norm_tok <- function(x) { x <- sub("^\\+", "", x); x }

# ---- 1 숫자 추적 ---------------------------------------------------------------------------------------------------------------------
for (sid in meta$id) {
  allowed <- unique(unlist(lapply(tr[section == sid, value_printed], num_tokens)))
  allowed <- c(allowed, sub("^-", "", allowed))
  body <- S[id == sid & !grepl("^footer_|^tag$|^background$", shape)]
  miss <- character(0)
  for (k in seq_len(nrow(body))) { tk <- num_tokens(body$text[k]); bad <- setdiff(sub("^-", "", tk), sub("^-", "", allowed))
    if (length(bad)) miss <- c(miss, sprintf("%s: %s", body$shape[k], paste(bad, collapse = " "))) }
  if (length(miss)) add("1 numbers traced", sid, "FAIL", paste(miss, collapse = " | ")) else add("1 numbers traced", sid, "pass", sprintf("%d traced values", nrow(tr[section == sid & !value_printed %in% c("(figure)", "(table)")])))
}
# ---- 2, 3 빈 값·대시 ---------------------------------------------------------------------------------------------------------------------
for (sid in meta$id) {
  t_ <- paste(S[id == sid, text], collapse = "\n")
  if (grepl("[{}]|(^|[^A-Za-z])(NA|NaN|Inf|TODO|XXX)([^A-Za-z]|$)", t_)) add("2 empty values", sid, "FAIL", regmatches(t_, regexpr(".{0,40}([{}]|NA|NaN|Inf|TODO|XXX).{0,40}", t_))) else add("2 empty values", sid, "pass")
  # 글자뿐 아니라 슬라이드·노트 XML 전체(글머리표 문자 a:buChar, 대체 텍스트 등)를 본다
  i_ <- match(sid, meta$id); rf_ <- file.path(dirname(slide_files[i_]), "_rels", paste0(basename(slide_files[i_]), ".rels"))
  xml_ <- paste(readLines(slide_files[i_], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (file.exists(rf_)) { tg_ <- xml_attr(xml_children(read_xml(rf_)), "Target"); for (nt_ in tg_[grepl("notesSlide", tg_)])
    xml_ <- paste(xml_, paste(readLines(normalizePath(file.path(dirname(slide_files[i_]), nt_)), warn = FALSE, encoding = "UTF-8"), collapse = "\n")) }
  if (grepl("[–—−]", t_)) add("3 dashes", sid, "FAIL", "en-dash, em-dash or U+2212")
  else if (grepl("(?<![0-9A-Za-z])\u2011(?=[0-9])", t_, perl = TRUE)) add("3 dashes", sid, "FAIL", "non-breaking hyphen U+2011 used as a minus sign (use ASCII '-')")
  else if (grepl("[–—−]", xml_)) add("3 dashes", sid, "FAIL", "en-dash, em-dash or U+2212 in the slide XML (bullet character or attribute)") else add("3 dashes", sid, "pass")
}
# ---- 4 아토피 ---------------------------------------------------------------------------------------------------------------------------
ok_sl <- unlist(CK$atopic_allowed_slides)
for (sid in setdiff(meta$id, ok_sl)) {
  t_ <- tolower(paste(S[id == sid, text], collapse = "\n")); hit <- unlist(CK$atopic_terms)[vapply(unlist(CK$atopic_terms), function(w) grepl(tolower(w), t_, fixed = TRUE), TRUE)]
  src <- unique(tr[section == sid & grepl(paste(unlist(CK$atopic_sources), collapse = "|"), source_file), source_file])
  if (length(hit) || length(src)) add("4 atopic outside A5", sid, "FAIL", paste(c(hit, src), collapse = ", ")) else add("4 atopic outside A5", sid, "pass")
}
# ---- 5 약어 첫 등장 ---------------------------------------------------------------------------------------------------------------------
ab <- CK$abbreviations
for (a in names(ab)) {
  rx <- sprintf("(?<![A-Za-z])%s(?![A-Za-z])", gsub("([&.])", "\\\\\\1", a))   # PKM12350, PKNCA 같은 이름 안의 글자는 약어가 아니다
  first <- NA_character_; for (sid in meta$id) { t_ <- paste(vis[id == sid, text], collapse = "\n"); if (grepl(rx, gsub("‑", "-", t_), perl = TRUE)) { first <- sid; break } }
  if (is.na(first)) next
  t_ <- gsub("‑", "-", paste(vis[id == first, text], collapse = "\n"))
  # 풀이는 읽는 순서로 첫 등장 앞이나 바로 뒤(괄호 풀이)에 있어야 한다
  pa <- regexpr(rx, t_, perl = TRUE)[1]; pe <- regexpr(ab[[a]], t_, perl = TRUE)[1]
  if (pe > 0 && pe <= pa + nchar(a) + 3) add("5 abbreviation first use", first, "pass", a)
  else add("5 abbreviation first use", first, "FAIL", if (pe < 0) sprintf("%s first used without '%s'", a, ab[[a]]) else sprintf("%s first used before its expansion '%s' (reading order)", a, ab[[a]]))
}
# ---- 6 형식 ---------------------------------------------------------------------------------------------------------------------------
for (k in seq_len(nrow(S))) {
  r <- S[k]; if (is.na(r$min_sz) || grepl("^footer_|^tag$|^background$", r$shape) || r$kind == "notes" || !nzchar(trimws(r$text))) next   # 글자 없는 도형(배경·화살표) 제외
  lim <- if (r$kind == "graphicFrame") 1200L else if (r$shape %in% c("kicker")) 1300L else if (grepl("^figure_", r$shape)) 0L else 1600L
  if (r$min_sz < lim) add("6 font size", r$id, "FAIL", sprintf("%s %.1fpt < %.0fpt", r$shape, r$min_sz / 100, lim / 100))
}
for (k in seq_len(nrow(meta))) { m <- meta[k]
  if (m$bullets > 5 || m$table_rows > 6) add("6 limits", m$id, "FAIL", sprintf("bullets %d, table rows %d", m$bullets, m$table_rows)) else add("6 limits", m$id, "pass", sprintf("bullets %d, table rows %d", m$bullets, m$table_rows)) }
# ---- 7 추적 완결 ------------------------------------------------------------------------------------------------------------------------
u <- unique(tr[, .(source_file, source_sha256)])
for (k in seq_len(nrow(u))) { f <- proj_path(u$source_file[k])
  if (!file.exists(f)) add("7 trace sources", "all", "FAIL", paste("missing", u$source_file[k]))
  else if (digest::digest(file = f, algo = "sha256") != u$source_sha256[k]) add("7 trace sources", "all", "FAIL", paste("changed since build", u$source_file[k])) }
add("7 trace sources", "all", "pass", sprintf("%d source files, %d trace rows", nrow(u), nrow(tr)))
# ---- 8 대조: M&S 보고서 추적표, key_numbers_en.md -----------------------------------------------------------------------------------------------------
rt <- fread(proj_path("regulatory", "traceability.csv"), encoding = "UTF-8", colClasses = list(character = c("value_raw", "value_printed", "locator")))
rt <- rt[!value_printed %in% c("(figure)", "(table)")]; dt <- tr[!value_printed %in% c("(figure)", "(table)")]
norm_p <- function(p) { p <- gsub("~", " to ", p); p <- gsub("%", "", p); gsub("\\s+", " ", trimws(p)) }
dec <- function(p) { m <- regmatches(p, gregexpr("[0-9]+\\.([0-9]+)", p))[[1]]; if (length(m)) max(nchar(sub("^[0-9]+\\.", "", m))) else 0L }
mm <- merge(dt[, .(section, item, source_file, locator, raw_d = value_raw, p_d = value_printed)], rt[, .(source_file, locator, raw_r = value_raw, p_r = value_printed, doc_r = document)],
            by = c("source_file", "locator"), allow.cartesian = TRUE)
nr <- 0L; np <- 0L
if (nrow(mm)) for (k in seq_len(nrow(mm))) { r <- mm[k]
  if (!identical(r$raw_d, r$raw_r)) { add("8 report cross-check", r$section, "FAIL", sprintf("%s [%s]: raw %s vs report %s", r$source_file, r$locator, r$raw_d, r$raw_r)); next }
  nr <- nr + 1L
  if (dec(r$p_d) == dec(r$p_r) && !identical(num_tokens(r$p_d), num_tokens(r$p_r))) add("8 report cross-check", r$section, "FAIL", sprintf("%s: printed '%s' vs report '%s'", r$item, r$p_d, r$p_r)) else np <- np + 1L }
add("8 report cross-check", "all", "pass", sprintf("%d deck values share file, filter and column with a report value; %d identical raw values; %d printed values consistent", nrow(mm), nr, np))
kn <- readLines(proj_path("results", "key_numbers_en.md"), encoding = "UTF-8", warn = FALSE)
kn_num <- list(); for (l in kn) { cit <- regmatches(l, gregexpr("[A-Za-z0-9_]+/[A-Za-z0-9_]+\\.csv", l))[[1]]; for (c_ in basename(cit)) kn_num[[c_]] <- unique(c(kn_num[[c_]], num_tokens(sub("\\([^()]*\\.csv[^()]*\\)\\.?$", "", l)))) }
kn_state <- function(p, f) { b <- basename(f); if (is.null(kn_num[[b]])) return("file not cited"); tk <- sub("^-", "", num_tokens(p)); kk <- sub("^-", "", kn_num[[b]])
  if (all(tk %in% kk)) return("matched")
  cl <- vapply(tk, function(t) any(abs(as.numeric(kk) - as.numeric(t)) <= 0.5 * 10^-min(dec(t), 3) + 1e-9, na.rm = TRUE), TRUE); if (all(cl)) "consistent after rounding" else "not in key numbers" }
dt[, kn := mapply(kn_state, value_printed, source_file)]
add("8 key numbers cross-check", "all", "pass", paste(sprintf("%s: %d", names(table(dt$kn)), as.integer(table(dt$kn))), collapse = "; "))
if (nrow(hl)) for (k in seq_len(nrow(hl))) { h <- hl[k]; r <- tr[h$idx]
  in_rep <- nrow(mm[section == h$section & p_d == h$printed]) > 0 || any(norm_p(rt$value_printed) == norm_p(h$printed))
  st <- kn_state(h$printed, r$source_file)
  if (in_rep || st %in% c("matched", "consistent after rounding")) add("8 headline", h$section, "pass", sprintf("%s (report: %s; key numbers: %s)", h$printed, in_rep, st))
  else add("8 headline", h$section, "FAIL", sprintf("%s not found in the report trace or key numbers (%s)", h$printed, r$source_file)) }

# ---- 9 렌더링 배치: PDF(LibreOffice)의 글줄이 자기 도형 밖으로 나가거나 다른 도형·그림과 겹치지 않는다 ----------------------------------------
pdf <- if ("--pdf" %in% args) args[match("--pdf", args) + 1] else paste0(base, ".pdf")
if (!("--pdf" %in% args) && !file.exists(pdf)) {   # render_deck.py의 기본 출력(<PNG 폴더>/<이름>.pdf)도 찾는다
  cand <- list.files(dirname(pptx), pattern = paste0("^", basename(base), "\\.pdf$"), recursive = TRUE, full.names = TRUE)
  if (length(cand)) pdf <- cand[which.max(file.mtime(cand))] }
if (!file.exists(pdf) || file.mtime(pdf) < file.mtime(pptx)) add("9 rendered layout", "all", "FAIL", sprintf("PDF missing or older than the pptx (%s); run render_deck.py", basename(pdf))) else {
  rc <- paste0(base, "_render.csv")
  system2("python3", c(proj_path("reports", "deck", "check_render.py"), shQuote(pptx), shQuote(pdf), shQuote(rc)), stdout = FALSE)
  rr <- fread(rc, encoding = "UTF-8", colClasses = list(character = c("text", "detail", "shape")))
  for (i in seq_len(nrow(meta))) { x <- rr[slide == i & kind != "unmatched"]
    if (nrow(x)) add("9 rendered layout", meta$id[i], "FAIL", paste(sprintf("%s [%s] %s: %s", x$kind, x$shape, x$detail, x$text), collapse = " | "))
    else add("9 rendered layout", meta$id[i], "pass", sprintf("%d unmatched lines", nrow(rr[slide == i & kind == "unmatched"]))) }
}

R <- rbindlist(res); out <- sub("\\.pptx$", "_checks.csv", pptx); if (!grepl("_checks\\.csv$", out)) out <- paste0(out, "_checks.csv")
fwrite(R, out); fwrite(dt[, .(section, item, source_file, value_printed, kn)], sub("_checks\\.csv$", "_keynumbers_xref.csv", out))
cat(sprintf("checks: %d pass, %d FAIL (%s)\n", sum(R$status == "pass"), sum(R$status == "FAIL"), out))
if (any(R$status == "FAIL")) { print(R[status == "FAIL"], nrows = 200); quit(status = 1) }
