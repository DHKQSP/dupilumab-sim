# 별첨 A3a ① 세부: 기준 세트별 결과(지시 2026-09-29). 결과보고 덱 S09(slides/s09_scale.R)의 자료 논리를 따른다. 시험 모집단(건강인, B0)만.
#  표(세트 (i)~(iv) x 두 모델, 8행): 미달(대상자 수준, Wilson 95% 구간), λz 산출 불가, λz가 산출되나 미달과 그 사유(계층 분류: adjusted R² 미달,
#  외삽 초과, span 미달; λz 산출 불가와 세 사유의 합 = 미달), arm당 충족 인원(동일 제품 시험, 중앙값과 5~95백분위). 세트 (iii) 행 강조(이 덱의 정의).
#  자료: results/trialpop/tp_failure_by_set.csv, tp_failure_reasons.csv, tp_retained_per_arm.csv. 세트 정의: config/prereg_20260926.yaml section4.
a3a_SETS <- c("i", "ii", "iii", "iv")
a3a_R <- c(lz = "lambda-z not estimable", r2 = "adjusted R-squared below threshold", ex = "extrapolation above 20%", sp = "span ratio below threshold")
a3a_w <- function(m, s_) sprintf("pk_model=='%s' & set=='%s'", m, s_)
# deck_table과 같은 식의 표 높이 추정(상자 높이를 표에 맞춘다)
a3a_table_h <- function(df, width, widths, size) {
  widths <- widths / sum(widths) * width
  nl <- function(v, w, b = FALSE) vapply(nobreak(as.character(v)), function(s) est_lines(s, w - 0.14, size, b), 1L)
  hdr <- max(mapply(function(v, w) max(nl(v, w, TRUE)), names(df), widths))
  bod <- apply(matrix(sapply(seq_along(widths), function(j) nl(df[[j]], widths[j], j == 1)), nrow = nrow(df)), 1, max)
  (hdr + sum(bod)) * size * 1.2 / 72 + (nrow(df) + 1) * 8 / 72
}

# deck_table과 같은 표(글꼴·크기·높이 검사·행 수 등록)에 세트 사이의 굵은 선을 더한다: 같은 세트의 두 모델 행 사이에는 선을 긋지 않고,
# 세트가 바뀌는 곳(group_end 행 아래)에 진한 선을 그어 모델 쌍을 묶어 보이게 한다(시각 검토 V, A3a).
a3a_table <- function(df, box, widths, size, highlight, group_end, label = "table") {
  stopifnot(size >= SZ$table_min)
  if (nrow(df) > LIMITS$table_rows) stop(sprintf("%s: table with %d body rows (limit %d)", REG$sec, nrow(df), LIMITS$table_rows), call. = FALSE)
  DK$cur$table_rows <- max(DK$cur$table_rows, nrow(df))
  df <- as.data.frame(df); for (j in seq_along(df)) df[[j]] <- nobreak(as.character(df[[j]])); names(df) <- nobreak(names(df))
  ft <- flextable(df)
  ft <- font(ft, fontname = FONT, part = "all", eastasia.family = FONT, hansi.family = FONT, cs.family = FONT)
  ft <- fontsize(ft, size = size, part = "all"); ft <- color(ft, color = PAL$ink, part = "all")
  ft <- bold(ft, part = "header"); ft <- bg(ft, bg = PAL$tint_blue, part = "header"); ft <- bold(ft, j = 1, part = "body")
  ft <- bg(ft, i = highlight, bg = PAL$tint_orange, part = "body")
  ft <- border_remove(ft)
  ft <- hline(ft, i = group_end, border = fp_border_default(color = PAL$ink2, width = 1.5), part = "body")
  ft <- hline_bottom(ft, border = fp_border_default(color = PAL$ink2, width = 1), part = "header"); ft <- hline_top(ft, border = fp_border_default(color = PAL$ink2, width = 1), part = "header")
  ft <- padding(ft, padding.top = 4, padding.bottom = 4, padding.left = 5, padding.right = 5, part = "all")
  ft <- align(ft, j = 2:ncol(df), align = "center", part = "all"); ft <- align(ft, j = 1, align = "left", part = "all"); ft <- valign(ft, valign = "center", part = "all")
  w <- widths / sum(widths) * box[3]; ft <- width(ft, width = w)
  h_est <- a3a_table_h(df, box[3], widths, size)
  DK$fit[[length(DK$fit) + 1L]] <- data.table(slide = REG$sec, shape = label, est_h = h_est, box_h = box[4], ratio = h_est / box[4])
  if (h_est > box[4] * 1.03 && isTRUE(DK$strict)) stop(sprintf("%s %s: estimated table height %.2f in exceeds box %.2f in", REG$sec, label, h_est, box[4]), call. = FALSE)
  DK$x <- ph_with(DK$x, ft, location = loc(box, label)); invisible(NULL)
}

slide_A3a <- function() {
  TPF <- "trialpop/tp_failure_by_set.csv"; TPR <- "trialpop/tp_failure_reasons.csv"; TRA <- "trialpop/tp_retained_per_arm.csv"; PR <- "config/prereg_20260926.yaml"
  deck_slide("A3a", tag = "litsim")
  L <- DK$txt$A3a; ML <- DK$txt$common$models_short; MOD <- c("k2016", "k2020")
  cs <- .read(PR)$section4$criteria_sets

  # ---- 전제: 표와 문장이 기대는 사실 ----
  exs <- vapply(a3a_SETS, function(s_) as.numeric(cs[[s_]]$extrap_max_pct), 0)
  premise(length(unique(exs)) == 1, "all four criteria sets share the extrapolation limit (table header)")
  fb <- rows(TPF); fr <- rows(TPR); ta <- rows(TRA)
  for (m in MOD) for (s_ in a3a_SETS) {
    x <- fb[pk_model == m & set == s_]; h <- fr[pk_model == m & set == s_]
    premise(nrow(x) == 1 && nrow(h) == 4 && setequal(h$reason, a3a_R), sprintf("%s set (%s): one failure row and four reasons", m, s_))
    premise(abs(sum(h$hierarchical_pct) - x$fail_pct) < 1e-6, sprintf("%s set (%s): hierarchical reasons add up to the failing share (table)", m, s_))
    premise(abs(h[reason == a3a_R[["lz"]], hierarchical_pct] - x$lz_pct) < 1e-9 && abs(x$fail_pct - x$lz_pct - x$est_fail_pct) < 1e-6,
            sprintf("%s set (%s): not estimable + estimable but failing = failing", m, s_))
    premise(is.null(cs[[s_]]$span_ratio_min) == (h[reason == a3a_R[["sp"]], hierarchical_pct] == 0), sprintf("%s set (%s): span failures only where the set has a span condition", m, s_))
    premise(h[reason == a3a_R[["r2"]], hierarchical_pct] > 0.9 * x$est_fail_pct || !is.null(cs[[s_]]$span_ratio_min),
            sprintf("%s set (%s): without a span condition, adjusted R-squared is most of the estimable failures (bullet)", m, s_))
    premise(nrow(ta[pk_model == m & set == s_]) == 1 && ta[pk_model == m & set == s_, retained_p95] <= as.numeric(.read("config/trial_design.yaml")$n_per_arm), "retained per arm at most the evaluable count")
  }
  for (m in MOD) premise(fb[pk_model == m][order(fail_pct)]$set[1] == "i" && fb[pk_model == m][order(-fail_pct)]$set[1] == "iv", sprintf("%s: set (i) fails least and set (iv) most (text)", m))
  r2s <- vapply(a3a_SETS, function(s_) as.numeric(cs[[s_]]$adj_r2_min), 0)
  premise(r2s[["i"]] == min(r2s) && is.null(cs$i$span_ratio_min), "set (i) has the lowest adjusted R-squared and no span condition (most lenient)")
  premise(is.null(cs$iii$span_ratio_min) && r2s[["iii"]] == max(r2s), "set (iii) = the deck definition: highest adjusted R-squared, no span condition")
  sap <- .read(PR)$section6$s2_1_failure_by_set$public_saps; ncts <- regmatches(sap, gregexpr("NCT[0-9]{8}", sap))[[1]]
  SRX <- c(NCT04117607 = "at least", NCT04441905 = "at least", NCT04700163 = "above")
  premise(setequal(ncts, names(SRX)) && all(vapply(names(SRX), function(z) core_sap(z, SRX[[z]], sprintf("public SAP %s: adjusted R-squared %s", z, SRX[[z]])) == f_set("iii", "r2"), TRUE)),
          "three public SAPs, each with the set (iii) adjusted R-squared threshold")

  # ---- 제목 ----
  f <- list(iii = drange(TPF, "set=='iii'", "fail_pct", 0, "%", "share without a reliable AUCinf, set (iii), two models"),
            i = drange(TPF, "set=='i'", "fail_pct", 0, "%", "share failing set (i), two models"))
  y0 <- core_title(tx("A3a.title", f), tx("A3a.kicker"))

  # ---- 표: 세트 x 모델 8행 ----
  mk_row <- function(s_, m, k) {
    w <- a3a_w(m, s_); x <- row1(TPF, w); r <- row1(TRA, w)
    setc <- if (k == 1) fill(L$table$set, list(s = L$sets[[s_]], r2 = f_set(s_, "r2")))
            else if (!is.null(cs[[s_]]$span_ratio_min)) fill(L$table$span, list(v = f_set(s_, "span"))) else if (s_ == "iii") L$table$deck else ""
    fail <- dderived(sprintf("share failing with Wilson 95%% interval, set (%s), %s", s_, m), TPF, sprintf("%s :: fail_pct (fail_lo to fail_hi)", w),
                     c(x$fail_pct, x$fail_lo, x$fail_hi), sprintf("%s%% (%s~%s)", fnum(x$fail_pct, 1), fnum(x$fail_lo, 1), fnum(x$fail_hi, 1)))
    rs <- function(k_) { if (k_ == "sp" && is.null(cs[[s_]]$span_ratio_min)) return(L$table$no_span)
      dv(TPR, sprintf("%s & reason=='%s'", w, a3a_R[[k_]]), "hierarchical_pct", 1, "%", sprintf("failing reason %s (hierarchical), set (%s), %s", k_, s_, m)) }
    ret <- dderived(sprintf("evaluable subjects per arm meeting set (%s), median (5th to 95th percentile), %s", s_, m), TRA, sprintf("%s :: retained_median (retained_p05 to retained_p95)", w),
                    c(r$retained_median, r$retained_p05, r$retained_p95), sprintf("%s (%s~%s)", fnum(r$retained_median, 0), fnum(r$retained_p05, 0), fnum(r$retained_p95, 0)))
    c(setc, ML[[m]], fail, dv(TPF, w, "lz_pct", 1, "%", sprintf("lambda-z not estimable, %s", m)),
      dv(TPF, w, "est_fail_pct", 1, "%", sprintf("estimable lambda-z but failing, set (%s), %s", s_, m)), rs("r2"), rs("ex"), rs("sp"), ret)
  }
  M <- do.call(rbind, lapply(a3a_SETS, function(s_) rbind(mk_row(s_, "k2016", 1), mk_row(s_, "k2020", 2))))
  df <- as.data.frame(M, stringsAsFactors = FALSE); names(df) <- tx("A3a.table.head", list(ex = f_set("i", "extrap"), n = f_n_arm()))
  WD <- c(2.1, 1.05, 1.95, 0.95, 1.05, 1.0, 0.95, 1.05, 1.75)
  TH <- a3a_table_h(df, GEO$CW, WD, 14) + 0.04
  a3a_table(df, box = c(GEO$ML, y0, GEO$CW, TH), widths = WD, size = 14, highlight = 5:6, group_end = c(2, 4, 6, 8), label = "table_sets")

  # ---- 요점과 캡션 ----
  RR <- a3a_R[["r2"]]; SP <- a3a_R[["sp"]]
  bl <- tx("A3a.bullets", list(
    lz = drange(TPF, "set=='iii'", "lz_pct", 1, "%", "lambda-z not estimable, two models"),
    r2h = drange(TPR, sprintf("set=='iii' & reason=='%s'", RR), "hierarchical_pct", 1, "%", "set (iii) adjusted R-squared below 0.90 (hierarchical), two models"),
    nsap = dderived("public SAPs with the set (iii) adjusted R-squared threshold (count of NCT identifiers)", PR, "section6.s2_1_failure_by_set.public_saps :: count of NCT identifiers", length(ncts), as.character(length(ncts))),
    r2 = f_set("iii", "r2")))
  cap <- tx("A3a.caption", list(n = dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"), na = f_n_arm(),
                                ntr = dint(TRA, "pk_model=='k2016' & set=='i'", "n_trials", "identical-product trials per model")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  BY <- y0 + TH + 0.08
  deck_bullets(bl, box = c(GEO$ML, BY, GEO$CW, capy - 0.04 - BY), size = 18, gap_pt = 6)

  # ---- 노트 ----
  deck_notes(tx("A3a.notes", list(
    wt = f_wt_range(), n = dint(TPF, "pk_model=='k2016' & set=='i'", "n", "subjects per model"), ntr = dint(TRA, "pk_model=='k2016' & set=='i'", "n_trials", "identical-product trials per model"),
    na = f_n_arm(), ex = f_set("i", "extrap"), r2i = f_set("i", "r2"), r2 = f_set("iii", "r2"), sp2 = f_set("ii", "span"), sp4 = f_set("iv", "span"),
    i = drange(TPF, "set=='i'", "fail_pct", 1, "%", "trial population, set (i) failing, two models"),
    ii = drange(TPF, "set=='ii'", "fail_pct", 1, "%", "set (ii) failing, two models"),
    iii = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "trial population, set (iii) failing, two models"),
    iv = drange(TPF, "set=='iv'", "fail_pct", 1, "%", "set (iv) failing, two models"),
    sp4a = drange(TPR, sprintf("set=='iv' & reason=='%s'", SP), "any_pct", 1, "%", "span ratio below 3 (any flag), two models"),
    sp4h = drange(TPR, sprintf("set=='iv' & reason=='%s'", SP), "hierarchical_pct", 1, "%", "set (iv) span below 3 (hierarchical), two models"),
    exa = drange(TPR, "set=='iii' & reason=='extrapolation above 20%'", "any_pct", 1, "%", "extrapolation above 20% (any flag), two models"),
    lam = drange(TRA, "set=='lambda'", "retained_median", 0, "", "median subjects per arm with an estimable lambda-z, two models"),
    rmin = dext(TRA, "set=='iii'", "retained_min", min, 0, "", "smallest subjects per arm meeting set (iii) over all trials, two models"),
    ri = drange(TRA, "set=='i'", "retained_median", 0, "", "median evaluable subjects per arm meeting set (i), two models"),
    ri3 = drange(TRA, "set=='iii'", "retained_median", 0, "", "median evaluable subjects per arm with a reliable AUCinf, set (iii), two models"),
    ncts = paste(ncts, collapse = ", "),
    s3 = core_sap("NCT04700163", "above", "public SAP NCT04700163: adjusted R-squared above"))))
  deck_end()
}
