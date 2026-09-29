# A1 부록: 문헌 정합성. 공개된 AUC0-last/AUC0-inf 평균비(문헌 NCA)와 같은 조건 모의(NCA, 참값; 두 모델)의 비교 표 4행,
# 아래에 원개발사 문헌의 정성 서술(프로젝트 문헌 발췌표). [문헌+모의]
# 결과 파일 literature_numeric.csv의 출처·판정 열에 한글이 있어 행은 영문 조건(Clot 행: source=="Clot 2021 Table 3", PKM12350 행: grepl)으로 고르고,
# 표에 보이는 자료 이름은 step1g_literature_coverage.csv의 영문 출처 열에서 읽는다. 일치 판정은 수치로 다시 계산하고 결과 파일의 판정과 같은지 검사한다.

a1_label <- function(rel, where, item) {
  r <- row1(rel, where); x <- as.character(r$source)
  premise(!grepl("[ㄱ-ㆎ가-힣]", x), sprintf("ASCII source label in %s [%s]", rel, where))
  .record(item, rel, sprintf("%s :: source", where), x, x); x
}

slide_A1 <- function() {
  LN <- "literature/literature_numeric.csv"; G16 <- "step1/step1g_literature_coverage.csv"; G20 <- "results/step1_k2020/step1g_literature_coverage.csv"
  deck_slide("A1", tag = "litsim")
  ln <- rows(LN); premise(nrow(ln) == 4 && sum(ln$source == "Clot 2021 Table 3") == 3 && sum(grepl("PKM12350", ln$source)) == 1, "literature table: three Clot 2021 doses and the PKM12350 control arm")
  # 판정 허용 폭(%p): 결과 파일 판정 문구의 '±n%p'에서 읽는다
  wj <- 'source=="Clot 2021 Table 3" & dose_mg==300'
  tol <- local({ s <- row1(LN, wj)$judgment; m <- regmatches(s, regexec("±([0-9]+)%p", s))[[1]]; premise(length(m) == 2, "tolerance in the judgment label")
    x <- as.numeric(m[2]); list(x = x, p = dderived("agreement tolerance of the NCA mean ratio (points), from the judgment label", LN, sprintf("%s :: judgment, regex +/-([0-9]+)%%p", wj), x, fnum(x, 0))) })
  W <- list(c200 = 'source=="Clot 2021 Table 3" & dose_mg==200', c300 = wj, c600 = 'source=="Clot 2021 Table 3" & dose_mg==600', pkm = "grepl('PKM12350', source)")
  WG <- list(c200 = "grepl('Clot', source) & dose_mg==200", c300 = "grepl('Clot', source) & dose_mg==300", c600 = "grepl('Clot', source) & dose_mg==600", pkm = "grepl('PKM12350', source)")
  agree_lab <- "일치"   # 결과 파일 판정 열의 '일치' 접두어(검사용)
  row_of <- function(k) {
    w <- W[[k]]; r <- row1(LN, w); v <- function(col, it) dv(LN, w, col, 1, "%", sprintf("%s, %s", k, it), scale = 100)
    d16 <- 100 * (r$nca_mean_ratio_k2016 - r$lit_mean_ratio); d20 <- 100 * (r$nca_mean_ratio_k2020 - r$lit_mean_ratio)
    agree <- abs(d16) <= tol$x + 1e-12 && abs(d20) <= tol$x + 1e-12
    premise(agree == startsWith(r$judgment, agree_lab), sprintf("recomputed agreement equals the stored judgment (%s)", k))
    dd <- dderived(sprintf("%s, simulated NCA minus published mean ratio (points), 2016 and 2020", k), LN, sprintf("%s :: (nca_mean_ratio_k2016, nca_mean_ratio_k2020 - lit_mean_ratio) x 100", w),
                   c(d16, d20), sprintf("%s / %s", fnum(d16, 1), fnum(d20, 1)))
    wt <- dv(G16, WG[[k]], "weight_mean", 1, "", sprintf("%s, mean body weight of the simulated cohort (kg)", k))
    premise(row1(G16, WG[[k]])$weight_mean == row1(G20, WG[[k]])$weight_mean, sprintf("same weight in both models (%s)", k))
    lab <- a1_label(G16, WG[[k]], sprintf("%s, source label", k))
    if (k == "pkm") {   # 공개 값은 PKM12350 대조군 arm뿐(결과 파일 출처 열에 '대조군'이 있다): 표의 자료 이름에 밝힌다
      ctrl <- tx("A1.table.ctrl")
      premise(grepl(sprintf("(PKM12350 %s)", ctrl), r$source, fixed = TRUE) && grepl("(PKM12350)", lab, fixed = TRUE), "PKM12350 literature row is the control arm")
      lab <- sub("(PKM12350)", sprintf("(PKM12350 %s)", ctrl), lab, fixed = TRUE)
    }
    interp <- if (agree) tx("A1.interp.agree", list(d = dd)) else {
      premise(r$nca_mean_ratio_k2016 < r$lit_mean_ratio && r$nca_mean_ratio_k2020 < r$lit_mean_ratio && r$true_mean_ratio_k2016 < r$lit_mean_ratio && r$true_mean_ratio_k2020 < r$lit_mean_ratio,
              sprintf("simulated NCA and true mean ratios below the published NCA ratio in both models (%s)", k))
      tx("A1.interp.lower", list(d = dd)) }
    if (k == "c200") interp <- paste0(interp, "\n", tx("A1.interp.ext"))   # 판정·차이 한 줄, 덧붙임 한 줄(가운데 정렬 두 줄이 고르게)
    if (k == "pkm") interp <- paste0(interp, "\n", tx("A1.interp.wt"))
    c(tx("A1.table.src", list(lab = lab, dose = dint(LN, w, "dose_mg", sprintf("%s, dose (mg)", k)), wt = wt)),
      v("lit_mean_ratio", "published mean AUC0-last/AUC0-inf ratio"),
      sprintf("%s / %s", v("nca_mean_ratio_k2016", "simulated NCA mean ratio, 2016"), v("nca_mean_ratio_k2020", "simulated NCA mean ratio, 2020")),
      sprintf("%s / %s", v("true_mean_ratio_k2016", "simulated true mean ratio, 2016"), v("true_mean_ratio_k2020", "simulated true mean ratio, 2020")),
      interp)
  }
  keys <- c("c300", "pkm", "c600", "c200")
  m <- t(vapply(keys, row_of, character(5))); df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- tx("A1.table.head", list(tol = tol$p))
  agree_of <- function(k) { r <- row1(LN, W[[k]]); abs(r$nca_mean_ratio_k2016 - r$lit_mean_ratio) * 100 <= tol$x + 1e-12 && abs(r$nca_mean_ratio_k2020 - r$lit_mean_ratio) * 100 <= tol$x + 1e-12 }
  premise(agree_of("c300") && agree_of("pkm") && !agree_of("c600") && row1(LN, W$pkm)$dose_mg == 300, "title: both 300 mg rows agree, the 600 mg row does not")
  r600 <- row1(LN, W$c600)
  premise(r600$true_mean_ratio_k2016 < r600$lit_mean_ratio && r600$true_mean_ratio_k2020 < r600$lit_mean_ratio &&
          r600$nca_mean_ratio_k2016 < r600$lit_mean_ratio && r600$nca_mean_ratio_k2020 < r600$lit_mean_ratio,
          "title: at 600 mg the simulated true mean ratio (and the simulated NCA ratio) is below the published NCA ratio in both models")
  deck_kicker(tx("A1.kicker")); premise(as.numeric(.read("config/trial_design.yaml")$dose_mg) == 300, "study dose equals the 300 mg literature rows")
  deck_title(tx("A1.title", list(tol = tol$p, dose = f_dose(), d600 = dint(LN, W$c600, "dose_mg", "dose of the lower-coverage literature row (mg)"))))
  # 해석 열: 판정과 차이 한 줄 + 덧붙임 한 줄(명시적 줄바꿈). 일치 기준과 반올림 전 계산은 머리글에, 하한의 뜻은 표 아래 설명에. 강조 행 없음(색만으로 뜻을 나타내지 않는다)
  th <- 2.8
  deck_table(df, box = c(GEO$ML, GEO$BODY_TOP, GEO$CW, th), widths = c(3.45, 1.05, 2.2, 2.2, 3.33), size = 13, align_num = TRUE, align_cols = c("left", "center", "center", "center", "left"))
  yc <- GEO$BODY_TOP + th + 0.08; hc <- 0.42
  deck_text(tx("A1.caption"), c(GEO$ML, yc, GEO$CW, hc), size = 16, color = PAL$ink2, label = "caption_lowerbound", gap_pt = 0)
  # 정성 서술(프로젝트 문헌 발췌표 literature_qualitative.csv; 보고서 1.1절)
  dsrc("qualitative statements of the originator literature", "literature/literature_qualitative.csv", "(table)")
  yq <- yc + hc + 0.1; hq <- GEO$BODY_BOTTOM - yq; wq <- (GEO$CW - 0.3) / 2
  deck_text(tx("A1.qual_label"), c(GEO$ML, yq, GEO$CW, 0.4), size = 16, bold = TRUE, color = PAL$ink2, label = "label_qual", gap_pt = 0)
  q <- tx("A1.qual")   # 왼쪽 두 개(각 두 줄), 오른쪽 세 개(두 줄 + 한 줄 + 한 줄)로 두 단의 높이를 맞춘다
  deck_bullets(q[1:2], c(GEO$ML, yq + 0.4, wq, hq - 0.4), size = 16, gap_pt = 4, label = "qual_left")
  deck_bullets(q[3:5], c(GEO$ML + wq + 0.3, yq + 0.4, wq, hq - 0.4), size = 16, gap_pt = 4, label = "qual_right")
  deck_notes(tx("A1.notes", list(
    tol = tol$p,
    n16 = dv(G16, "grepl('Clot', source) & dose_mg==300", "mean_ratio_nca_reliable", 1, "%", "Clot 300 mg, simulated NCA mean ratio in subjects meeting the criteria, 2016", scale = 100),
    t16 = dv(G16, "grepl('Clot', source) & dose_mg==300", "tlast_cohort_median", 1, "", "Clot 300 mg, median of simulated cohort median tlast, 2016 (day)"),
    t20 = dv(G20, "grepl('Clot', source) & dose_mg==300", "tlast_cohort_median", 1, "", "Clot 300 mg, median of simulated cohort median tlast, 2020 (day)"),
    lt = local({ s <- row1(LN, wj)$literature; m <- regmatches(s, regexec("tlast ([0-9]+) ", s))[[1]]; premise(length(m) == 2, "published tlast in the literature label")
      dderived("Clot 300 mg, published median tlast (day), from the literature label", LN, sprintf("%s :: literature, regex tlast ([0-9]+)", wj), as.numeric(m[2]), m[2]) }),
    pt = local({ s <- row1(LN, W$pkm)$note; m <- regmatches(s, gregexpr("[0-9]+", s))[[1]]; premise(length(m) >= 2, "test-arm numbers in the PKM12350 note")
      dderived("PKM12350 test arm published AUC0-t and AUC0-inf (from the note)", LN, "grepl('PKM12350', source) :: note, first two integers", as.numeric(m[1:2]), sprintf("%s / %s", m[1], m[2])) }),
    ns = dint(G16, "grepl('Clot', source) & dose_mg==300", "n", "simulated subjects per literature condition, 2016"))))
  deck_end()
}
