#!/usr/bin/env Rscript
# make_template.R: 결과보고 슬라이드의 16:9 템플릿(reports/deck/template_16x9.pptx)을 officer 기본 템플릿에서 만든다.
# 바꾸는 것은 두 가지뿐이다: 슬라이드 크기(13.333 x 7.5 in = 12192000 x 6858000 EMU)와 테마 글꼴(라틴·동아시아 모두 Pretendard).
# 슬라이드는 모두 Blank 레이아웃에 절대 위치로 그리므로 다른 레이아웃의 자리표시자 위치는 쓰지 않는다.
# 사용법: Rscript reports/deck/make_template.R
source("R/00_setup.R")
src <- file.path(find.package("officer"), "template", "template.pptx")
out <- proj_path("reports", "deck", "template_16x9.pptx")
tmp <- tempfile("tpl_"); dir.create(tmp)
utils::unzip(src, exdir = tmp)
pres <- file.path(tmp, "ppt", "presentation.xml")
x <- readLines(pres, warn = FALSE, encoding = "UTF-8")
x <- sub('<p:sldSz cx="9144000" cy="6858000" type="screen4x3"/>', '<p:sldSz cx="12192000" cy="6858000"/>', x, fixed = TRUE)
stopifnot(any(grepl('<p:sldSz cx="12192000" cy="6858000"/>', x, fixed = TRUE)))
writeLines(x, pres, useBytes = TRUE)
th <- file.path(tmp, "ppt", "theme", "theme1.xml")
y <- readLines(th, warn = FALSE, encoding = "UTF-8")
y <- gsub('<a:latin typeface="Calibri"/>', '<a:latin typeface="Pretendard"/>', y, fixed = TRUE)
y <- gsub('<a:ea typeface=""/>', '<a:ea typeface="Pretendard"/>', y, fixed = TRUE)
stopifnot(sum(lengths(regmatches(y, gregexpr('typeface="Pretendard"', y)))) >= 4)
writeLines(y, th, useBytes = TRUE)
if (file.exists(out)) file.remove(out)
old <- setwd(tmp); on.exit(setwd(old))
zip::zip(out, files = list.files(".", recursive = TRUE, all.files = TRUE), mode = "mirror")
setwd(old)
tpl <- officer::read_pptx(out); s <- officer::slide_size(tpl)
stopifnot(abs(s$width - 13.333) < 0.01, abs(s$height - 7.5) < 0.01)
cat("template written:", out, sprintf("(%.3f x %.3f in)\n", s$width, s$height))
