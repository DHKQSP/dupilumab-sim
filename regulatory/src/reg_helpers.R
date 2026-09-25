# Helpers for the regulatory (FDA-facing) documents in regulatory/.
# Every number printed in a regulatory document is read from a committed result or config file through the functions below,
# which record the source file, the row filter, the column, the raw value and the printed value (traceability, Appendix G).
# English only in the outputs. Plain ASCII hyphen and the word "to" are used for ranges (no en-dash, em-dash or U+2212).
suppressPackageStartupMessages({ library(data.table) })

REG <- new.env()
reg_init <- function(doc) { REG$doc <- doc; REG$sec <- ""; REG$trace <- list(); REG$cache <- list(); invisible(NULL) }
sec <- function(s) { REG$sec <- s; invisible("") }

.rel_path <- function(rel) if (grepl("^(config|results|regulatory)/", rel)) rel else file.path("results", rel)
# Some analysis outputs carry Korean category labels (the analysis team's working language). They are mapped to the English
# codes below when read, so that filters, tables and the traceability file are in English. The mapping is written to
# regulatory/tables/label_translation.csv.
LABEL_EN <- c("생산자 위험(범위 안 불통과)" = "producer", "소비자 위험(범위 밖 통과)" = "consumer", "경계 근처" = "near_boundary", "전체" = "all",
              "AUC0-inf·Cmax 모두" = "AUC0-inf and Cmax", "동일" = "identical", "경계" = "boundary", "범위 밖" = "outside", "범위 안" = "inside",
              "체중 60-75" = "weight 60-75", "체중 >75-90" = "weight >75-90",
              "Theoph (12명)" = "Theoph (12 profiles)", "Indometh (6명, 혈관외 규칙)" = "Indometh (6 profiles, extravascular rules)",
              "두필루맙 모의 1,000개 (2016 500 + Model 1 500)" = "simulated dupilumab (1,000 profiles: 500 per model)",
              "자체 vs NonCompart" = "this engine vs NonCompart", "자체 vs PKNCA" = "this engine vs PKNCA")
.translate <- function(d) {
  for (cn in names(d)) if (is.character(d[[cn]])) { x <- d[[cn]]; hit <- x %in% names(LABEL_EN); if (any(hit)) { x[hit] <- unname(LABEL_EN[x[hit]]); set(d, j = cn, value = x) } }
  d
}
.read <- function(rel) {
  p <- .rel_path(rel)
  if (is.null(REG$cache[[p]])) {
    fp <- proj_path(p)
    if (!file.exists(fp)) stop("regulatory document: missing source file ", p, call. = FALSE)
    REG$cache[[p]] <- if (grepl("\\.ya?ml$", p)) yaml::read_yaml(fp) else .translate(fread(fp, encoding = "UTF-8"))
  }
  REG$cache[[p]]
}
.record <- function(item, rel, locator, raw, printed) {
  REG$trace[[length(REG$trace) + 1L]] <- data.table(document = REG$doc, section = REG$sec, item = item, source_file = .rel_path(rel),
                                                    locator = locator, value_raw = paste(format(raw, digits = 15), collapse = "; "), value_printed = printed)
  invisible(printed)
}
# lookups used inline return their text visibly (knitr prints nothing for an invisible inline value)
.vis <- function(x) { force(x); x }   # force first: evaluating the promise would otherwise leave R_Visible off
premise <- function(ok, msg) if (!isTRUE(all(ok))) stop("regulatory document premise failed: ", msg, call. = FALSE)

# rounding: half away from zero (as in the analysis outputs, R/summarize.R round_half_away)
fnum <- function(x, d = 1, big = FALSE) { x <- round(x + sign(x) * 1e-9, d) + 0; formatC(x, format = "f", digits = d, big.mark = if (big) "," else "") }
fint <- function(x) formatC(round(x), format = "d", big.mark = ",")

# rows of a result table matching a filter given as text (recorded verbatim)
rows <- function(rel, where = "TRUE") { d <- .read(rel); d[eval(parse(text = where), envir = d)] }
row1 <- function(rel, where) { r <- rows(rel, where); premise(nrow(r) == 1, sprintf("%s [%s] must match exactly one row (matched %d)", rel, where, nrow(r))); r }

# single value
v <- function(rel, where, col, d = 1, item = col, big = FALSE, scale = 1) {
  r <- row1(rel, where); x <- r[[col]] * scale
  .vis(.record(item, rel, sprintf("%s :: %s%s", where, col, if (scale != 1) sprintf(" x %s", scale) else ""), x, fnum(x, d, big)))
}
vi <- function(rel, where, col, item = col) { r <- row1(rel, where); x <- r[[col]]; .vis(.record(item, rel, sprintf("%s :: %s", where, col), x, fint(x))) }
vtxt <- function(rel, where, col, item = col) { r <- row1(rel, where); x <- as.character(r[[col]]); .vis(.record(item, rel, sprintf("%s :: %s", where, col), x, x)) }
# estimate with 95% interval: "5.18% (95% CI 4.88 to 5.50)"
vci <- function(rel, where, est = "pass_pct", lo = "lo", hi = "hi", d = 2, item = est, unit = "%", scale = 1, style = c("paren", "bracket")) {
  style <- match.arg(style); r <- row1(rel, where); e <- r[[est]] * scale; l <- r[[lo]] * scale; h <- r[[hi]] * scale
  txt <- if (style == "paren") sprintf("%s%s (95%% CI %s to %s)", fnum(e, d), unit, fnum(l, d), fnum(h, d)) else sprintf("%s%s [%s, %s]", fnum(e, d), unit, fnum(l, d), fnum(h, d))
  .vis(.record(item, rel, sprintf("%s :: %s [%s, %s]", where, est, lo, hi), c(e, l, h), txt))
}
# range over matching rows: "a to b"
vrange <- function(rel, where, col, d = 1, item = col, unit = "", scale = 1) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); x <- range(r[[col]] * scale)
  .vis(.record(item, rel, sprintf("range over %d rows [%s] :: %s", nrow(r), where, col), x, sprintf("%s%s to %s%s", fnum(x[1], d), unit, fnum(x[2], d), unit)))
}
# the row with the largest (or smallest) value of col among rows matching where; returns list(text, row); the recorded
# filter identifies that row by model and scenario
vci_ext <- function(rel, where, col = "pass_pct", lo = "lo", hi = "hi", d = 2, item = col, fun = which.max) {
  r <- rows(rel, where); premise(nrow(r) >= 1, sprintf("%s [%s] matched no rows", rel, where)); k <- fun(r[[col]]); rr <- r[k]
  w2 <- sprintf("%s & model=='%s' & scenario=='%s'", where, rr$model, rr$scenario)
  list(text = vci(rel, w2, col, lo, hi, d, item), row = rr)
}
vmax <- function(rel, where, col, d = 1, item = col) { r <- rows(rel, where); x <- max(r[[col]]); .vis(.record(item, rel, sprintf("max over %d rows [%s] :: %s", nrow(r), where, col), x, fnum(x, d))) }
vcount <- function(rel, where, item = "count") { r <- rows(rel, where); .vis(.record(item, rel, sprintf("count of rows [%s]", where), nrow(r), as.character(nrow(r)))) }
# config value by path (character vector of keys)
cfgv <- function(file, path, item = paste(path, collapse = "."), fmt = function(x) paste(format(x), collapse = ", ")) {
  y <- .read(file.path("config", file)); for (k in path) y <- y[[k]]
  premise(!is.null(y), sprintf("config/%s: %s not found", file, paste(path, collapse = ".")))
  .vis(.record(item, file.path("config", file), paste(path, collapse = "."), unlist(y), fmt(unlist(y))))
}
# table-level provenance (every value in a rendered table comes from these files)
tab_src <- function(item, rels, locator = "whole table") { for (r_ in rels) .record(item, r_, locator, NA, "(table)"); invisible(NULL) }
fig <- function(rel, caption, item = caption) {
  p <- .rel_path(rel); premise(file.exists(proj_path(p)), paste("figure missing:", p)); .record(item, rel, "figure", NA, "(figure)")
  knitr::include_graphics(proj_path(p))
}

# labels used in the regulatory documents (internal codes are not shown)
MODEL_LAB <- c(k2016 = "Kovalenko 2016 (primary model)", k2020 = "Kovalenko 2020 Model 1 (structural sensitivity model)")
MODEL_SHORT <- c(k2016 = "2016 model", k2020 = "2020 Model 1")
MECH_LAB <- c(F = "bioavailability (F)", ka = "absorption rate constant (ka)", ke = "linear elimination rate constant (ke)",
              Vmax = "maximum target-mediated elimination rate (Vmax)", Km = "Michaelis-Menten constant (Km)", V2 = "peripheral volume (V2)")
CONFIG_LAB <- c(P2 = "AUC0-last + Cmax (proposed)", AUClast_only = "AUC0-last alone (reference)",
                G2_Aii = "AUC0-inf + Cmax, rule A, criteria (ii) (guidance default as pre-specified)", G2_Ai = "AUC0-inf + Cmax, rule A, criteria (i)",
                G2_B = "AUC0-inf + Cmax, rule B", G2_Cii = "AUC0-inf + Cmax, rule C, criteria (ii)", G2_Ci = "AUC0-inf + Cmax, rule C, criteria (i)",
                F3A = "AUC0-last + AUC0-inf + Cmax, rule A, criteria (ii)", F3B = "AUC0-last + AUC0-inf + Cmax, rule B", F3C = "AUC0-last + AUC0-inf + Cmax, rule C, criteria (ii)",
                AUCinf_Aii = "AUC0-inf alone, rule A, criteria (ii)", AUCinf_Ai = "AUC0-inf alone, rule A, criteria (i)", AUCinf_B = "AUC0-inf alone, rule B",
                AUCinf_Cii = "AUC0-inf alone, rule C, criteria (ii)", AUCinf_Ci = "AUC0-inf alone, rule C, criteria (i)")
scen_lab <- function(mechanism, direction, target, multiplier = NULL) {
  s <- sprintf("%s %s, true AUC0-inf ratio %s", mechanism, ifelse(direction == "up", "increased", "decreased"), formatC(as.numeric(target), format = "f", digits = 2))
  if (!is.null(multiplier)) s <- sprintf("%s (x%s)", s, formatC(as.numeric(multiplier), format = "f", digits = 2))
  s
}

# write the traceability rows of this document
reg_trace_write <- function(path) {
  tr <- rbindlist(REG$trace, fill = TRUE)
  files <- unique(tr$source_file)
  sh <- vapply(files, function(f) digest::digest(file = proj_path(f), algo = "sha256"), "")
  tr[, source_sha256 := sh[source_file]]
  fwrite(tr, path); invisible(tr)
}

# English-output check: no Korean, no en-dash, em-dash or U+2212 in a generated text file
check_english <- function(path) {
  x <- readLines(path, encoding = "UTF-8", warn = FALSE)
  bad <- grep("[–—−]|[가-힣]", x, value = TRUE)
  if (length(bad)) stop("non-English or dash characters in ", path, ": ", paste(substr(head(bad, 3), 1, 160), collapse = " || "), call. = FALSE)
  invisible(TRUE)
}
