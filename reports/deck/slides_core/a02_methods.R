# A2 별첨 · 방법 세부: 비구획 분석(NCA) 규칙(config/nca_rules.yaml standard), λz 신뢰 기준 세트 (i)~(iv)(config/prereg_20260926.yaml section4 criteria_sets, f_set),
# Phoenix Lambda Z Acceptance Criteria가 사용자 선택 항목이라는 사실(사전 등록 4절 facts_for_documents, Certara Phoenix 8.2 NCA 문서)과 공개 SAP 기준(core_sap),
# 엔진 검증(results/nca_engine/engine_validation_summary.csv; 행 조건은 LABEL_EN 영문 코드). 자료 논리는 결과보고 덱 A4·S06과 같다. [문헌+모의]
# 배치: 왼쪽 = 규칙 요점 3개, 오른쪽 = 기준 세트 표(14pt, 4행, 등록 상태 포함), 아래 = Phoenix·공개 SAP 세 건 상자(전체 폭, 회색), 캡션(엔진 검증, 시험 모의 방법).

# 상대 차이를 가수x10의 위첨자 지수로(코드 표기 e-NN 꼴 대신; 결과보고 덱 A4와 같은 표기)
a2_sci <- function(x) { e <- floor(log10(abs(x))); m <- signif(x / 10^e, 2); if (m >= 10) { m <- m / 10; e <- e + 1 }
  sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹", "-" = "⁻")
  paste0(format(m), "×10", paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }
# 공개 SAP 문장(사전 등록 6절)에서 NCT04117607의 span 요건(반감기 배수)을 읽는다
a2_sap_span <- function() {
  P <- c("section6", "s2_1_failure_by_set", "public_saps"); x <- .read("config/prereg_20260926.yaml")[[P[1]]][[P[2]]][[P[3]]]
  m <- regmatches(x, regexec("NCT04117607 \\(adjusted R-squared at least ([0-9.]+), span at least ([0-9.]+) half-lives", x))[[1]]
  premise(length(m) == 3, "NCT04117607: adjusted R-squared and span in prereg section6 public_saps")
  list(r2 = as.numeric(m[2]), span = as.numeric(m[3]))
}

slide_A2 <- function() {
  NR <- "nca_rules.yaml"; EV <- "nca_engine/engine_validation_summary.csv"; SE <- "regulatory/tables/software_environment.csv"; PR <- "regulatory/tables/prespecification_register.csv"
  deck_slide("A2", tag = "litsim")
  L <- DK$txt$A2
  nr <- .read(file.path("config", NR))$standard; pr <- .read("config/prereg_20260926.yaml")$section4; cs <- pr$criteria_sets

  # ---- 전제: 슬라이드 문구와 config가 같다(결과보고 덱 A4와 같은 검사) ----
  premise(identical(nr$auc_method, "linear_up_log_down") && identical(nr$blq$pre_first_quant, "zero") && identical(nr$blq$intermediate, "missing") &&
            identical(nr$blq$post_tlast, "exclude") && identical(nr$blq$two_consecutive_blq, "truncate"), "BLQ and AUC rules in nca_rules.yaml match the slide text")
  premise(identical(nr$lambda_z$selection, "max_adj_r2") && identical(nr$lambda_z$tie_rule, "more_points") && isTRUE(nr$lambda_z$after_tmax_only) && isTRUE(nr$lambda_z$exclude_cmax) &&
            identical(nr$lambda_z$positive_slope_windows, "excluded_before_selection") && identical(nr$aucinf, "AUClast + Clast_obs / lambda_z"), "lambda-z Best Fit and AUCinf rules match the slide text")
  premise(nr$reliability$adj_r2_min == cs$ii$adj_r2_min && nr$reliability$extrap_max_pct == cs$ii$extrap_max_pct && nr$reliability$span_ratio_min == cs$ii$span_ratio_min,
          "engine reliability flag in nca_rules.yaml equals criteria set (ii) (table: engine default flag)")
  premise(is.null(cs$i$span_ratio_min) && is.null(cs$iii$span_ratio_min) && !is.null(cs$ii$span_ratio_min) && !is.null(cs$iv$span_ratio_min), "span only in sets (ii) and (iv)")
  premise(length(unique(vapply(cs, function(z) z$extrap_max_pct, 1))) == 1, "same extrapolation limit in every set")
  fd <- unlist(pr$facts_for_documents)
  premise(any(grepl("optional user entries", fd)) && any(grepl("flagged \\(Accepted / Not_Accepted\\), not excluded", fd)) && any(grepl("no acceptance values are applied unless entered", fd)) &&
            any(grepl("Certara Phoenix 8.2", fd)), "Phoenix Lambda Z Acceptance Criteria: optional user entries; not applied unless entered; failing profiles flagged, not excluded (prereg section4 facts)")
  premise(grepl("lower end of the conventional range", pr$reporting) && grepl("set \\(i\\)", pr$reporting), "set (i) reported as the lower end of the conventional range (table)")
  sp <- a2_sap_span(); premise(sp$r2 == cs$iv$adj_r2_min && sp$span == cs$iv$span_ratio_min, "set (iv) equals the NCT04117607 criteria (table)")
  reg <- rows(PR)
  premise(nrow(reg[grepl("^Reliability criteria set \\(i\\)", Item) & Status == "post hoc" & grepl("set \\(ii\\) is the pre-specified definition", `How reported`)]) == 1,
          "register: set (i) post hoc, set (ii) the pre-specified definition (notes)")
  premise(nrow(reg[grepl("^Criteria sets \\(iii\\) and \\(iv\\)", Item) & grepl("^pre-registered \\(section4\\)", Status)]) == 1, "register: sets (iii) and (iv) pre-registered in section 4 (notes)")
  dsrc("status of the criteria sets (pre-specified, post hoc, pre-registered)", PR, "(table)")

  f <- list(r2i = f_set("i", "r2"), r2iii = f_set("iii", "r2"))
  y0 <- core_title(tx("A2.title", f), tx("A2.kicker"))

  # ---- 엔진 검증과 시험 모의 방법(캡션 두 줄) ----
  ev <- rows(EV); premise(nrow(ev) == 9 && all(ev$pass) && all(ev$lz_points_mismatch == 0) && all(ev$na_mismatch == 0), "every NCA engine comparison passed (9 rows)")
  premise(setequal(unique(unlist(strsplit(unique(ev$comparison), " vs ", fixed = TRUE))), c("this engine", "NonCompart", "PKNCA")), "engines: this engine, NonCompart, PKNCA")
  premise(!any(grepl("phoenix", list.files(proj_path("results/nca_engine"), recursive = TRUE), ignore.case = TRUE)), "no Phoenix output in results/nca_engine (caption: not compared with Phoenix itself)")
  np_ <- ev[, .(n = n_profiles[1]), by = dataset]; premise(all(ev[, uniqueN(n_profiles), by = dataset]$V1 == 1), "same profile count in every comparison of a data set")
  lz_i <- sum(ev$lz_points_identical); lz_t <- sum(ev$lz_points_identical + ev$lz_points_mismatch)
  premise(lz_t == sum(np_$n) * length(unique(ev$comparison)), "compared windows = profiles x engine pairs")
  g <- list(lz = dderived("lambda-z windows identical over all comparisons", EV, "sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", c(lz_i, lz_t), sprintf("%s/%s", fint(lz_i), fint(lz_t))),
            ntot = dderived("profiles per engine comparison (Theoph + Indometh + simulated)", EV, "sum over datasets of n_profiles (one comparison each)", sum(np_$n), fint(sum(np_$n))),
            mx = local({ x <- max(ev$max_rel_diff); dderived("largest relative parameter difference (scientific notation)", EV, "max(max_rel_diff) :: printed as mantissa x 10^exponent", x, a2_sci(x)) }))
  oc <- .read("config/oc_design.yaml")
  premise(grepl("기하평균비", oc$estimand$quantity) && grepl("적분", oc$estimand$auc) && grepl("공통 난수", oc$estimand$population$paired_crn), "truth: geometric mean ratio of model-integral values on common-random-number subjects (caption)")
  premise(grepl("이분법", oc$inversion$method) && identical(as.numeric(unlist(oc$boundary_targets)), as.numeric(unlist(.read("config/trial_design.yaml")$be$limits))), "boundary cells: multiplier found by bisection so that the true AUCinf ratio equals an equivalence limit (caption)")
  mth <- list(nt = dcfg("oc_design.yaml", c("estimand", "population", "n_subjects"), "common-random-number subjects of the true ratio", function(x) fnum(as.numeric(x), 0, big = TRUE)), reps = f_reps("boundary"))
  capy <- core_caption(c(tx("A2.caption", g), tx("A2.method", mth)), GEO$BODY_BOTTOM, size = 14)

  # ---- 왼쪽: NCA 규칙 요점(3개 + AUCinf 하위 요점) ----
  XL <- GEO$ML; WL <- 6.1; XR <- XL + WL + 0.3; WR <- GEO$W - GEO$MR - XR; LH <- 0.40
  mp <- dcfg(NR, c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0))
  tol <- dcfg(NR, c("standard", "lambda_z", "tie_tolerance"), "lambda-z adjusted R-squared tie tolerance", num_fmt(4))
  z0 <- dderived("BLQ before the first quantifiable value is set to zero", file.path("config", NR), "standard.blq.pre_first_quant == 'zero'", nr$blq$pre_first_quant, "0")
  deck_text(tx("A2.rules_label"), c(XL, y0, WL, LH), size = 16, bold = TRUE, color = PAL$ink2, label = "label_rules", gap_pt = 0)
  rl <- tx("A2.rules", list(mp = mp, tol = tol, z = z0))
  rh <- est_height(vapply(rl, function(s) sub("^- ", "", s), ""), WL, SZ$body, gap_pt = 5, indent = 0.3) + 0.06
  deck_bullets(rl, box = c(XL, y0 + LH + 0.02, WL, rh), size = SZ$body, gap_pt = 5, label = "body_rules")

  # ---- 오른쪽 위: 기준 세트 표(등록 상태 포함; 외삽 한도는 네 세트 공통이라 머리말에) ----
  deck_text(tx("A2.sets_label", list(ex = f_set("i", "extrap"))), c(XR, y0, WR, LH), size = 16, bold = TRUE, color = PAL$ink2, label = "label_sets", gap_pt = 0)
  sets <- c("i", "ii", "iii", "iv")
  spn <- function(s_) if (s_ %in% c("ii", "iv")) fill(L$span, list(x = f_set(s_, "span"))) else L$none   # span 요건은 세트 (ii), (iv)만(전제)
  use <- function(s_) if (s_ == "iv") fill(L$use$iv, list(ex = f_set("iv", "extrap"))) else L$use[[s_]]
  premise(!is.null(cs$i), "set (i) is part of the registered criteria sets (table: defined post hoc, then registered)")
  df <- data.frame(a = unlist(L$set_names[sets]), b = vapply(sets, function(s_) fill(L$ge, list(x = f_set(s_, "r2"))), ""),
                   d = vapply(sets, spn, ""), s = unlist(L$status[sets]), e = vapply(sets, use, ""), stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tx("A2.table.head")
  TW_ <- c(0.5, 1.03, 1.0, 1.39, 1.91); TY <- y0 + LH + 0.02
  TH <- deck_table_h(df, WR, TW_, size = 14) + 0.18                  # 렌더링 표 높이(추정보다 조금 크다)
  deck_table(df, box = c(XR, TY, WR, TH), widths = TW_, size = 14, label = "table_sets", highlight = 3, group_end = 4)   # 행 사이 선 없이 아래 굵은 선(표 끝 표시)

  # ---- 오른쪽 아래: 공개 SAP 세 건(사전 등록 6절 문장에서 읽는다) ----
  sap <- core_sap("NCT04117607", "at least", "public SAP NCT04117607: adjusted R-squared at least")
  s2 <- core_sap("NCT04441905", "at least", "public SAP NCT04441905: adjusted R-squared at least")
  s3 <- core_sap("NCT04700163", "above", "public SAP NCT04700163: adjusted R-squared above")
  premise(s2 == sap && s3 == sap && as.numeric(sap) == cs$iii$adj_r2_min, "three public SAPs use the set (iii) threshold value")
  psx <- .read("config/prereg_20260926.yaml")$section6$s2_1_failure_by_set$public_saps
  premise(grepl("NCT04117607 \\([^)]*; verified\\)", psx) && grepl("NCT04441905 \\([^)]*; verified\\)", psx) && grepl("NCT04700163 \\(.*web-search excerpt", psx) &&
          grepl("otherwise AUC\\(INF\\) and other terminal parameters excluded", psx), "SAP verification status and the NCT04700163 exclusion rule (table)")
  sp4 <- dderived("public SAP NCT04117607: span at least (half-lives)", "config/prereg_20260926.yaml", "section6.s2_1_failure_by_set.public_saps :: regex 'NCT04117607 \\(adjusted R-squared at least ([0-9.]+), span at least ([0-9.]+) half-lives' group 2",
                  sp$span, fnum(sp$span, 0))
  # ---- 아래 전체 폭: Phoenix 사실과 공개 SAP 세 건(회색 상자) ----
  S <- L$sap
  ph <- tx("A2.phoenix", list(s1 = fill(S$ge, list(x = sap)), sp = sp4, s2 = fill(S$ge, list(x = s2)), s3 = fill(S$gt, list(x = s3))))
  phh <- est_height(ph, GEO$CW, SZ$body, gap_pt = 4, card = TRUE) + 0.10
  PY <- max(y0 + LH + 0.02 + rh, TY + TH) + 0.10
  premise(PY + phh <= capy - 0.08, "Phoenix and SAP box fits between the columns and the caption")
  deck_text(ph, c(GEO$ML, PY, GEO$CW, phh), size = SZ$body, label = "body_phoenix", bg = PAL$tint_grey, geom = "roundRect", gap_pt = 4)

  deck_notes(tx("A2.notes", c(f, g, mth, list(mp = mp, tol = tol, z = z0, ex = f_set("i", "extrap"), sp2 = f_set("ii", "span"), sp4 = f_set("iv", "span"),
    s1 = sap, s2 = s2, s3 = s3,
    nc = dderived("NonCompart version (renv.lock)", SE, "Software=='NonCompart' :: Version", row1(SE, "Software=='NonCompart'")$Version, row1(SE, "Software=='NonCompart'")$Version),
    pk = dderived("PKNCA version (renv.lock)", SE, "Software=='PKNCA' :: Version", row1(SE, "Software=='PKNCA'")$Version, row1(SE, "Software=='PKNCA'")$Version),
    nth = dint(EV, "comparison=='this engine vs NonCompart' & dataset=='Theoph (12 profiles)'", "n_profiles", "Theoph profiles"),
    nin = dint(EV, "comparison=='this engine vs NonCompart' & dataset=='Indometh (6 profiles, extravascular rules)'", "n_profiles", "Indometh profiles"),
    nsim = dint(EV, "comparison=='this engine vs NonCompart' & dataset=='simulated dupilumab (1,000 profiles: 500 per model)'", "n_profiles", "simulated profiles")))))
  deck_end()
}
