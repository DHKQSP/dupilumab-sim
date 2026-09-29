# A4 부록: NCA 엔진 규칙과 신뢰 기준 세트 정의. [문헌+모의]
# 규칙: config/nca_rules.yaml(standard). 기준 세트: config/prereg_20260926.yaml section4$criteria_sets(f_set). 세트의 지위(사전 명시·사후·사전 등록)는
# regulatory/tables/prespecification_register.csv(영문)에서 premise로 확인한다. Phoenix Lambda Z Acceptance Criteria가 사용자 선택 항목이라는 사실은
# 보고서 3.4절과 prereg section4 facts_for_documents(Certara Phoenix 8.2 NCA 문서)의 서술이다(문헌).
# 엔진 검증: results/nca_engine/engine_validation_summary.csv. 이 파일의 dataset·comparison 열은 한글이며 .read()가 LABEL_EN으로 영문 코드로 바꾼다
# (예: "simulated dupilumab (1,000 profiles: 500 per model)", "this engine vs NonCompart"). 행 조건은 영문 코드로만 쓴다(추적표에 한글 금지).
# 참조 구현 버전은 regulatory/tables/software_environment.csv(renv.lock에서 만든 영문 표)에서 읽는다.

# 공개 SAP의 adjusted R² 기준: 사전 등록 6절(s2_1_failure_by_set$public_saps, 확인 방법 포함)의 문장에서 정규식으로 읽는다
A4_SAP <- c("section6", "s2_1_failure_by_set", "public_saps")
a4_sap <- function(nct, rx, item) {
  x <- .read("config/prereg_20260926.yaml")[[A4_SAP[1]]][[A4_SAP[2]]][[A4_SAP[3]]]; m <- regmatches(x, regexec(sprintf("%s \\(adjusted R-squared %s ([0-9.]+)", nct, rx), x))[[1]]
  premise(length(m) == 2, sprintf("public SAP threshold for %s (%s) in prereg %s", nct, rx, paste(A4_SAP, collapse = ".")))
  v <- as.numeric(m[2]); dderived(item, "config/prereg_20260926.yaml", sprintf("%s :: regex '%s \\(adjusted R-squared %s ([0-9.]+)'", paste(A4_SAP, collapse = "."), nct, rx), v, fnum(v, 2))
}
# 상대 차이를 가수x10의 위첨자 지수로(코드 표기 e-NN 꼴 대신)
a4_sci <- function(x) { e <- floor(log10(abs(x))); m <- signif(x / 10^e, 2); if (m >= 10) { m <- m / 10; e <- e + 1 }
  sup <- c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079", "-" = "\u207b")
  paste0(format(m), "\u00d710", paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }

a4_txt <- function(rel, where, col, item) {
  x <- as.character(row1(rel, where)[[col]]); premise(!grepl("[ㄱ-ㆎ가-힣]", x), sprintf("ASCII text in %s [%s] %s", rel, where, col))
  dderived(item, rel, sprintf("%s :: %s", where, col), x, x)
}

slide_A4 <- function() {
  NR <- "nca_rules.yaml"; EV <- "nca_engine/engine_validation_summary.csv"; SE <- "regulatory/tables/software_environment.csv"; PR <- "regulatory/tables/prespecification_register.csv"
  deck_slide("A4", tag = "litsim")
  L <- DK$txt$A4
  nr <- .read(file.path("config", NR))$standard; pr <- .read("config/prereg_20260926.yaml")$section4

  # ---- 전제: 규칙 문구와 config가 같다 ----
  premise(identical(nr$auc_method, "linear_up_log_down") && identical(nr$blq$pre_first_quant, "zero") && identical(nr$blq$intermediate, "missing") &&
            identical(nr$blq$post_tlast, "exclude") && identical(nr$blq$two_consecutive_blq, "truncate"), "BLQ and AUC rules in nca_rules.yaml match the slide text")
  premise(identical(nr$lambda_z$selection, "max_adj_r2") && identical(nr$lambda_z$tie_rule, "more_points") && isTRUE(nr$lambda_z$after_tmax_only) && isTRUE(nr$lambda_z$exclude_cmax) &&
            identical(nr$lambda_z$positive_slope_windows, "excluded_before_selection") && identical(nr$aucinf, "AUClast + Clast_obs / lambda_z"), "lambda-z Best Fit and AUCINF_obs rules match the slide text")
  premise(grepl("제외", nr$aucinf_rules$A) && grepl("포함", nr$aucinf_rules$B) && grepl("AUClast", nr$aucinf_rules$C), "rules A (exclude), B (include), C (substitute AUClast) in nca_rules.yaml")
  cs <- pr$criteria_sets
  premise(nr$reliability$adj_r2_min == cs$ii$adj_r2_min && nr$reliability$extrap_max_pct == cs$ii$extrap_max_pct && nr$reliability$span_ratio_min == cs$ii$span_ratio_min,
          "engine reliability flag in nca_rules.yaml equals criteria set (ii)")
  premise(is.null(cs$i$span_ratio_min) && is.null(cs$iii$span_ratio_min) && !is.null(cs$ii$span_ratio_min) && !is.null(cs$iv$span_ratio_min), "span only in sets (ii) and (iv)")
  premise(any(grepl("optional user entries", unlist(pr$facts_for_documents))) && any(grepl("flagged .*not excluded", unlist(pr$facts_for_documents))) &&
            any(grepl("Certara Phoenix 8.2", unlist(pr$facts_for_documents))), "Phoenix Lambda Z Acceptance Criteria: optional, flag only (prereg section4 facts)")
  # 세트의 지위: 사전 명시 등록부(영문)
  reg <- rows(PR)
  premise(nrow(reg[grepl("^Reliability criteria set \\(i\\)", Item) & Status == "post hoc" & grepl("set \\(ii\\) is the pre-specified definition", `How reported`)]) == 1,
          "register: set (i) post hoc, set (ii) the pre-specified definition")
  premise(nrow(reg[grepl("^Criteria sets \\(iii\\) and \\(iv\\)", Item) & grepl("^pre-registered \\(section4\\)", Status)]) == 1, "register: sets (iii) and (iv) pre-registered in section 4")
  dsrc("status of the criteria sets (pre-specified, post hoc, pre-registered)", PR, "(table)")

  f <- list(r2i = f_set("i", "r2"), r2iii = f_set("iii", "r2"))
  # 공개 SAP 세 건의 기준(카드): 사전 등록 6절 문장. NCT04117607·NCT04441905는 "at least"(≥), NCT04700163은 "above"(>)이고 셋 모두 세트 (iii)의 값과 같다
  sap <- list(s1 = a4_sap("NCT04117607", "at least", "public SAP NCT04117607: adjusted R-squared at least"),
              s2 = a4_sap("NCT04441905", "at least", "public SAP NCT04441905: adjusted R-squared at least"),
              s3 = a4_sap("NCT04700163", "above", "public SAP NCT04700163: adjusted R-squared above"))
  premise(all(as.numeric(unlist(sap)) == cs$iii$adj_r2_min), "the three public SAP thresholds equal the adjusted R-squared of criteria set (iii)")
  ps <- .read("config/prereg_20260926.yaml")$section6$s2_1_failure_by_set$public_saps
  premise(grepl("NCT04117607 \\([^)]*span", ps) && !grepl("NCT04441905 \\([^)]*span", ps) && !grepl("NCT04700163 \\([^)]*span", ps), "only NCT04117607 has a span requirement")
  premise(grepl("NCT04117607 \\([^)]*; verified\\)", ps) && grepl("NCT04441905 \\([^)]*; verified\\)", ps) && grepl("NCT04700163 \\(.*verified from a web-search excerpt", ps),
          "verification status of the three public SAPs (two verified, NCT04700163 from a web-search excerpt)")
  premise(grepl("lower end of the conventional range", pr$reporting) && grepl("set \\(i\\)", pr$reporting) && any(grepl("^0\\.80 is a statistical-analysis-plan convention", unlist(pr$facts_for_documents))),
          "0.80 (set i) reported as the lower end of the conventional range (prereg section4)")
  premise(cs$i$adj_r2_min == 0.8, "set (i) adjusted R-squared is the 0.80 named in prereg section4")
  deck_kicker(tx("A4.kicker")); deck_title(tx("A4.title"), box = c(GEO$ML, GEO$TITLE_TOP, GEO$CW - GEO$TAG_W - 0.2, GEO$TITLE_H))   # 제목이 오른쪽 위 태그 아래로 들어가지 않게

  # ---- 왼쪽: NCA 엔진 규칙(요점 5개: BLQ, AUC0-last, λz, AUC0-inf, 미달자 처리 규칙 A·B·C) + 엔진 검증 카드 ----
  XL <- GEO$ML; WL <- 6.05
  deck_text(tx("A4.rules_label"), c(XL, GEO$BODY_TOP, WL, 0.4), size = 16, bold = TRUE, color = PAL$ink2, label = "label_rules", gap_pt = 0)
  mp <- dcfg(NR, c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0))
  tol <- dcfg(NR, c("standard", "lambda_z", "tie_tolerance"), "lambda-z adjusted R-squared tie tolerance", num_fmt(4))
  BH <- 3.05
  z0 <- dderived("BLQ before the first quantifiable value is set to zero", file.path("config", NR), "standard.blq.pre_first_quant == 'zero'", nr$blq$pre_first_quant, "0")
  deck_bullets(tx("A4.rules", list(mp = mp, tol = tol, z = z0)), box = c(XL, GEO$BODY_TOP + 0.42, WL, BH), size = 16, gap_pt = 5)

  ev <- rows(EV); premise(nrow(ev) == 9 && all(ev$pass) && all(ev$lz_points_mismatch == 0) && all(ev$na_mismatch == 0), "every NCA engine comparison passed (9 rows)")
  lz_i <- sum(ev$lz_points_identical); lz_t <- sum(ev$lz_points_identical + ev$lz_points_mismatch)
  lz <- dderived("lambda-z windows identical over all comparisons", EV, "sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", c(lz_i, lz_t), sprintf("%s/%s", fint(lz_i), fint(lz_t)))
  mx <- max(ev$max_rel_diff); mxp <- dderived("largest relative parameter difference", EV, "max(max_rel_diff)", mx, a4_sci(mx))
  WSIM <- "comparison=='this engine vs NonCompart' & dataset=='simulated dupilumab (1,000 profiles: 500 per model)'"
  WTH <- "comparison=='this engine vs NonCompart' & dataset=='Theoph (12 profiles)'"
  WIN <- "comparison=='this engine vs NonCompart' & dataset=='Indometh (6 profiles, extravascular rules)'"
  g <- list(lz = lz, mx = mxp,
            nc = a4_txt(SE, "Software=='NonCompart'", "Version", "NonCompart version (renv.lock)"),
            pk = a4_txt(SE, "Software=='PKNCA'", "Version", "PKNCA version (renv.lock)"),
            nth = dint(EV, WTH, "n_profiles", "Theoph profiles"), nin = dint(EV, WIN, "n_profiles", "Indometh profiles"),
            nsim = dint(EV, WSIM, "n_profiles", "simulated profiles"))
  premise(nrow(ev) / length(unique(ev$dataset)) == 3 && length(unique(ev$comparison)) == 3, "three pairwise comparisons per data set (three pairs on the slide)")
  np_ <- ev[, .(n = n_profiles[1]), by = dataset]; premise(all(ev[, uniqueN(n_profiles), by = dataset]$V1 == 1), "same profile count in every comparison of a data set")
  g$ntot <- dderived("profiles per engine comparison (Theoph + Indometh + simulated)", EV, "sum over datasets of n_profiles (one comparison each)", sum(np_$n), fint(sum(np_$n)))
  g$npair <- dderived("engine pairs compared per data set", EV, "uniqueN(comparison)", length(unique(ev$comparison)), as.character(length(unique(ev$comparison))))
  premise(lz_t == sum(np_$n) * length(unique(ev$comparison)), "compared windows = profiles x engine pairs")
  cy <- GEO$BODY_TOP + 0.42 + BH + 0.1
  deck_stat(g$lz, tx("A4.engine", g), c(XL, cy, WL, GEO$BODY_BOTTOM - cy), color = PAL$blue, bg = PAL$tint_blue, value_size = 32)

  # ---- 오른쪽 위: 기준 세트 표 ----
  XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR
  deck_text(tx("A4.sets_label"), c(XR, GEO$BODY_TOP, WR, 0.4), size = 16, bold = TRUE, color = PAL$ink2, label = "label_sets", gap_pt = 0)
  sets <- c("i", "ii", "iii", "iv")
  sp <- function(s_) if (s_ %in% c("ii", "iv")) f_set(s_, "span") else L$none
  df <- data.frame(a = unlist(DK$txt$common$sets[sets]), b = vapply(sets, f_set, "", what = "r2"), c = vapply(sets, function(s_) paste0(f_set(s_, "extrap"), "%"), ""),
                   d = vapply(sets, sp, ""), e = unlist(L$status[sets]), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("A4.table.head")
  TY <- GEO$BODY_TOP + 0.42; TH <- 1.9
  deck_table(df, box = c(XR, TY, WR, TH), widths = c(0.8, 1.42, 0.96, 1.12, 1.58), size = 13, label = "table_sets")   # 강조 행 없음(모든 세트를 같은 지위로 보인다)

  # ---- 오른쪽 아래: Phoenix 사실 ----
  PY <- TY + TH + 0.15
  deck_text(tx("A4.phoenix", c(f, sap)), c(XR, PY, WR, GEO$BODY_BOTTOM - PY), size = 16, bg = PAL$tint_orange, geom = "roundRect", label = "text_phoenix", gap_pt = 6)

  deck_notes(tx("A4.notes", c(f, g, sap, list(mp = mp, tol = tol, ex = f_set("i", "extrap"), sp2 = f_set("ii", "span"), sp4 = f_set("iv", "span")))))
  deck_end()
}
