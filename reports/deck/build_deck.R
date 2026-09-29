#!/usr/bin/env Rscript
# build_deck.R: 결과보고 슬라이드(v1.0.1)를 만든다. 새 모의 없음: 커밋된 결과·config 파일만 읽는다.
# 사용법: Rscript reports/deck/build_deck.R [--lang ko] [--only S05,S09] [--out 경로] [--no-strict]
#   전체 빌드(--only 없음)는 reports/deck/dupilumab_endpoint_results_v1.0.1.pptx와 추적표(deck_traceability.csv),
#   슬라이드 목록(deck_meta.csv), 넘침 추정(deck_meta_fit.csv), 핵심 수치 목록(deck_meta_headline.csv)을 쓴다.
#   --only는 개발용: 지정한 슬라이드만 --out(기본 임시 경로)에 쓰고 추적표는 같은 폴더에 둔다.
source("R/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, default = NULL) { i <- match(k, args); if (is.na(i)) default else args[i + 1] }
lang <- opt("--lang", "ko"); only <- opt("--only"); strict <- !("--no-strict" %in% args)
source(proj_path("reports", "deck", "lib", "deck_lib.R"))
source(proj_path("reports", "deck", "lib", "deck_facts_common.R"))
for (f in sort(list.files(proj_path("reports", "deck", "slides"), pattern = "\\.R$", full.names = TRUE)))   # 부분 빌드에서는 다른 슬라이드 파일의 오류로 멈추지 않는다
  tryCatch(source(f), error = function(e) if (is.null(only)) stop(e) else message("skipped (does not source): ", basename(f), ": ", conditionMessage(e)))
deck_init(lang, strict = strict)
order <- unlist(DK$txt$common$order)
sel <- if (is.null(only)) order else intersect(order, strsplit(only, ",")[[1]])
if (!length(sel)) stop("no slide selected")
DK$total <- length(order)
for (id in sel) {
  fn <- get0(paste0("slide_", id), mode = "function")
  if (is.null(fn)) stop("slide function missing: slide_", id)
  if (!is.null(only)) DK$n <- match(id, order) - 1L          # 부분 빌드에서도 쪽 번호는 전체 순서
  fn()
  if (is.null(DK$meta[[length(DK$meta)]]) || DK$meta[[length(DK$meta)]]$id != id) stop("slide ", id, " did not call deck_end()")
}
base <- sprintf("dupilumab_endpoint_results_%s", DK$version)
if (is.null(only)) {
  out <- opt("--out", proj_path("reports", "deck", paste0(base, if (lang == "ko") "" else paste0("_", lang), ".pptx")))
  tr <- proj_path("reports", "deck", if (lang == "ko") "deck_traceability.csv" else sprintf("deck_traceability_%s.csv", lang))
  me <- proj_path("reports", "deck", if (lang == "ko") "deck_meta.csv" else sprintf("deck_meta_%s.csv", lang))
} else {
  out <- opt("--out", file.path(tempdir(), paste0(base, "_partial.pptx"))); dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE); tr <- sub("\\.pptx$", "_traceability.csv", out); me <- sub("\\.pptx$", "_meta.csv", out)
}
n <- deck_write(out, trace_path = tr, meta_path = me)
cat(sprintf("deck written: %s (%d slides, %d bullet paragraphs, %d traced values, commit %s%s)\n", out, length(sel), n,
            sum(!rbindlist(REG$trace)$value_printed %in% c("(figure)", "(table)")), DK$commit, if (DK$clean) "" else ", tree not clean"))
