# A5f 2종 오류 분해(v1.2, 지시 2026-09-29 "S9 재설계" §2 "2종 오류 분해(분석 세트 축소 대 AUCinf 잡음·편향; 규칙 A(iii) 대 규칙 B 비교)"와 §5 "2종 오류의 분석 세트 기여", D-065):
# 표 = AUCinf(A) + Cmax - AUClast + Cmax의 2종 오류 증가 = [AUCinf(B) + Cmax - AUClast + Cmax](AUCinf 추정) + [AUCinf(A) + Cmax - AUCinf(B) + Cmax](분석 세트 축소),
#   모델 x 칸 묶음(동일 제품, 참 비 0.95·1.05, 0.90·1.11; 참 Cmax 비가 한계 안인 칸)의 값 또는 범위, 군당 분석 대상자(AUCinf(A)) 중앙값.
#   results/oc_curves/oc_decomposition.csv(사전 등록 section8 statistics.decomposition; 같은 시험의 쌍대 차이, 두 몫의 합 = 합계).
slide_A5f <- function() {
  DC <- "oc_curves/oc_decomposition.csv"; PR <- "config/prereg_20260929_oc.yaml"
  deck_slide("A5f", tag = "sim")
  L <- DK$txt$A5f; ML <- DK$txt$common$models_short; MS <- names(ML)
  dc <- rows(DC, "config=='G2'"); premise(nrow(dc) > 0 && all(dc$reduction_lo > 0), "analysis-set reduction part significantly above 0 in every type II cell, AUCinf + Cmax (title)")
  premise(all(rows(DC, "config=='F3'")$reduction_lo > 0), "analysis-set reduction part significantly above 0 in every type II cell, three endpoints (notes)")
  premise(all(dc$reduction_pp > dc$estimation_pp), "AUCinf + Cmax: analysis-set part larger than the estimation part in every type II cell (title)")
  premise(all(abs(dc$total_pp - dc$estimation_pp - dc$reduction_pp) < 1e-6), "the two parts add up to the total (table)")
  premise(all(dc[total_pp > 0, reduction_pp >= 0.5 * total_pp]), "analysis-set part is at least half of every positive total (body line 1)")
  premise(dc[, mean(estimation_pp < 0)] > 0.5, "estimation part negative in most cells (body line 2)")
  nA <- rows(DC, "config=='G2' & scenario=='S00'")$n_arm_Aiii; nL <- rows(DC, "config=='G2' & scenario=='S00'")$n_arm_AUClast; premise(length(unique(nL)) == 1, "same AUClast analysis count in both models")
  f <- list(nl = dint(DC, "config=='G2' & scenario=='S00' & pk_model=='k2020'", "n_arm_AUClast", "median AUClast analysis count per arm, identical product"),
            na = dderived("median AUCinf (set iii) analysis count per arm, identical product, range over two models", DC, "config=='G2' & scenario=='S00' :: range(n_arm_Aiii)", range(nA), rng_fmt(min(nA), max(nA), 0)))
  y0 <- core_title(tx("A5f.title", f), tx("A5f.kicker"))

  # ---- 표 ----
  nt <- as.numeric(unlist(.read(PR)$section8$reps$near_one$targets)); premise(isTRUE(all.equal(nt, c(0.95, 1.05))), "near-one targets 0.95 and 1.05")
  mt <- c(0.90, 1.11)
  tg <- function(v, what) dderived(sprintf("type II group target (%s)", what), PR, sprintf("section8.statistics.type2 :: %s", what), v, fnum(v, 2))
  tt <- list(a = tg(nt[1], "(b) lower"), b = tg(nt[2], "(b) upper"), c = tg(mt[1], "(c) lower"), d = tg(mt[2], "(c) upper"))
  Wg <- function(m, g_) sprintf("config=='G2' & pk_model=='%s' & %s", m, switch(g_, identity = "scenario=='S00'", near = "(abs(target - 0.95) < 1e-9 | abs(target - 1.05) < 1e-9)",
                                                                           mid = "(abs(target - 0.90) < 1e-9 | abs(target - 1.11) < 1e-9)"))
  premise(all(vapply(MS, function(m) nrow(rows(DC, Wg(m, "identity"))) + nrow(rows(DC, Wg(m, "near"))) + nrow(rows(DC, Wg(m, "mid"))) == nrow(dc[pk_model == m]), TRUE)),
          "every AUCinf + Cmax type II cell is in one of the three groups")
  cell <- function(m, g_, col, d_ = 1) { w <- Wg(m, g_); r <- rows(DC, w); x <- range(r[[col]])
    if (g_ == "identity") dv(DC, w, col, d_, "", sprintf("type II decomposition %s, %s, identical product (percentage points)", col, m))
    else dderived(sprintf("type II decomposition %s, %s, %s cells, range (percentage points)", col, m, g_), DC, sprintf("%s :: range(%s)", w, col), x, rng_fmt(x[1], x[2], d_)) }
  nn <- function(m, g_) { w <- Wg(m, g_); r <- rows(DC, w); x <- range(r$n_arm_Aiii)
    dderived(sprintf("median AUCinf (set iii) analysis count per arm, %s, %s cells, range", m, g_), DC, sprintf("%s :: range(n_arm_Aiii)", w), x, if (x[1] == x[2]) fnum(x[1], 1) else rng_fmt(x[1], x[2], 1)) }
  ncell <- function(m, g_) dcount(DC, Wg(m, g_), sprintf("type II cells, %s, %s", m, g_))
  G <- c("identity", "near", "mid")
  tab <- rbindlist(lapply(MS, function(m) rbindlist(lapply(G, function(g_) data.table(
    a = fill(L$tab$rows[[g_]], c(list(m = ML[[m]]), tt)), b = cell(m, g_, "total_pp"), c = cell(m, g_, "reduction_pp"), d = cell(m, g_, "estimation_pp"), e = nn(m, g_))))))
  setnames(tab, unlist(L$tab$head))

  # ---- 본문·캡션 ----
  body <- c(tx("A5f.body1"), tx("A5f.body2"))
  cap <- tx("A5f.caption", list(r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)
  by <- core_body(body, capy - 0.06)
  wd <- c(3.8, 1.9, 2.1, 2.1, 2.3); wd <- wd * GEO$CW / sum(wd)
  th <- deck_table_h(tab, GEO$CW, wd, 16) + 0.04; ht <- est_height(L$tab$title, GEO$CW, 16, 0) + 0.06
  if (y0 + ht + 0.04 + th >= by - 0.1) message(sprintf("A5f table: y0 %.2f ht %.2f th %.2f by %.2f", y0, ht, th, by)); premise(y0 + ht + 0.04 + th < by - 0.1, "table fits above the body")
  deck_text(L$tab$title, c(GEO$ML, y0, GEO$CW, ht), size = 16, bold = TRUE, color = PAL$ink2, label = "caption_tabtitle")
  deck_table(tab, c(GEO$ML, y0 + ht + 0.04, GEO$CW, th), widths = wd, size = 16, group_end = 3, label = "table_decomp")
  deck_visual(c(GEO$ML, y0 + ht + 0.04, GEO$CW, th))

  # ---- 노트 ----
  f3 <- function(col) { r <- rows(DC, "config=='F3'"); x <- range(r[[col]])
    dderived(sprintf("three endpoints (F3) decomposition %s, all type II cells, range (percentage points)", col), DC, sprintf("config=='F3' :: range(%s)", col), x, rng_fmt(x[1], x[2], 1)) }
  deck_notes(tx("A5f.notes", c(f, list(
    nb = { r <- rows(DC, "config=='G2' & scenario=='S00'"); x <- range(r$n_arm_B); dderived("median AUCinf (rule B) analysis count per arm, identical product, range over two models", DC, "config=='G2' & scenario=='S00' :: range(n_arm_B)", x, rng_fmt(x[1], x[2], 1)) },
    f3t = f3("total_pp"), f3r = f3("reduction_pp"), f3e = f3("estimation_pp"),
    emin = dext(DC, "config=='G2'", "estimation_pp", min, 1, "", "most negative estimation part, AUCinf + Cmax (percentage points)"),
    emax = dext(DC, "config=='G2'", "estimation_pp", max, 1, "", "largest estimation part, AUCinf + Cmax (percentage points)"),
    rlo = dext(DC, "config=='G2'", "reduction_lo", min, 2, "", "smallest lower 95% limit of the analysis-set part, AUCinf + Cmax (percentage points)"),
    rn = dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "near_one", "trials"), "trials per cell, true ratio 0.95 and 1.05", function(x) fnum(as.numeric(x), 0, TRUE)),
    ro = dcfg("prereg_20260929_oc.yaml", c("section8", "reps", "other", "trials"), "trials per cell, true ratio 0.85, 0.90, 1.11 and 1.18", function(x) fnum(as.numeric(x), 0, TRUE)),
    rid = dint(DC, "config=='G2' & scenario=='S00' & pk_model=='k2020'", "n_trials", "trials, identical product"),
    n20n = ncell("k2020", "near"), n20m = ncell("k2020", "mid"), n16n = ncell("k2016", "near"), n16m = ncell("k2016", "mid")), tt)))
  deck_end()
}
