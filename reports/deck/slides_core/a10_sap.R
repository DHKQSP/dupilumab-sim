# A10 SAP·대응안(별첨): 카드 세 개 = (1) 제안 SAP(1차 AUClast + Cmax, 2차 AUCinf 두 분석군), (2) AUCinf 공동 1차 요구 시 대응안(AUClast 유지 + AUCinf 규칙 B 추가 =
# 세 지표 공동 1차, config F3B; 비용 = 동시 통과율 감소, results/fallback/fallback_cost.csv cost_iii_pp, 기준 검정력 empirical_power.csv), (3) AUCinf + Cmax만은 권하지 않음.
# 표 = 구성별 경계 1종 오류(M1: 점추정 > 5% 칸 수, 최대; type1_models.csv P2·F3B·G2_B·G2_Ai, criteria_g2_type1.csv G2_A_iii),
# AUCinf 분석에서 arm당 빠지는 인원(평가 가능 인원 - trialpop/tp_retained_per_arm.csv 중앙값), 빠지는 대상자의 참 AUCinf(남는 대상자 대비, tp_characteristics.csv true_aucinf_gmr).
# 규칙: A = 기준 미달자 제외(세트 (i) adjusted R² ≥ 0.80, (iii) 신뢰할 수 있는 AUCinf), B = λz 산출 전원, C = 미달자에 AUClast 대입. 자료 논리: 결과보고 덱 s23_sap.R, s24_questions.R.
slide_A10 <- function() {
  T1 <- "oc_models/type1_models.csv"; CG <- "criteria/criteria_g2_type1.csv"; FC <- "fallback/fallback_cost.csv"; EP <- "fallback/empirical_power.csv"
  TC <- "trialpop/tp_characteristics.csv"; TR <- "trialpop/tp_retained_per_arm.csv"; TPF <- "trialpop/tp_failure_by_set.csv"; PW <- "oc_models/power_models.csv"
  deck_slide("A10", tag = "sim")
  L <- DK$txt$A10
  w1 <- function(cf, extra = "") sprintf("analysis_model=='M1' & config=='%s'%s", cf, extra)
  n_arm <- as.numeric(.read("config/trial_design.yaml")$n_per_arm)

  # ---- 전제 ----
  for (cf in c("P2", "F3B", "G2_B", "G2_Ai")) premise(nrow(rows(T1, w1(cf))) == 16, sprintf("16 boundary cells, %s, M1", cf))
  premise(nrow(rows(CG, w1("G2_A_iii"))) == 16, "16 boundary cells, G2_A_iii, M1")
  OC <- .read("config/oc_design.yaml")$configurations
  premise(identical(unlist(OC$F3B$endpoints), c("AUClast", "Cmax", "AUCinf_B")) && grepl("λz", .read("config/oc_design.yaml")$aucinf_rules$B), "F3B = AUClast + Cmax + AUCinf rule B (all subjects with an estimable lambda-z)")
  premise(all(rows(T1, w1("F3B"))$class == "conservative"), "fallback F3B conservative in all 16 cells under M1 (card)")
  premise(all(rows(TC, "set %in% c('i','iii')")$true_aucinf_gmr_hi < 1), "failing subjects have a lower true AUCinf than retained subjects, sets (i) and (iii), both models (card: low-exposure subjects drop out)")
  premise(all(rows(TR, "set %in% c('i','iii','lambda')")$retained_median <= n_arm) && all(rows(TR, "set %in% c('i','iii','lambda')")$n_trials == rows(TR, "set=='iii'")$n_trials[1]), "retained per arm at most the evaluable count; one trial count")
  premise(all(rows(TC, "set=='iii'")$wt_diff_lo > 0), "set (iii) failing subjects are heavier than retained subjects, both models (notes)")
  fc <- rows(FC, "TRUE"); premise(all(fc$cost_iii_pp >= 0) && setequal(fc$scenario, c("S00", "KE110", "F097")) && all(fc$true_ratio > 0.94), "fallback cost: three near-equivalent products, non-negative loss (card: 'at most')")
  ep <- rows(EP, "model=='k2016'"); premise(all(abs(ep$power_last_cmax - fc[match(ep$scenario, fc$scenario), i_last_cmax]) < 1e-9), "empirical power and fallback cost come from the same trials (AUClast + Cmax column identical)")
  for (cf in c("G2_A_iii", "G2_B", "G2_Ai")) { r <- rows(if (cf == "G2_A_iii") CG else T1, w1(cf)); premise(sum(r$pass_pct > 5) >= 8, sprintf("%s: many cells above 5%% (card)", cf)) }

  # ---- 제목 ----
  f <- list(fb = dext(T1, w1("F3B"), "pass_pct", max, 2, "%", "largest boundary type I error, fallback AUClast + AUCinf (rule B) + Cmax, M1"), nom = f_nominal(),
            n = dcount(T1, w1("P2"), "boundary cells"))
  y0 <- core_title(tx("A10.title", f), tx("A10.kicker"))

  # ---- 캡션 ----
  cap <- tx("A10.caption", list(nom = f$nom, n_arm = f_n_arm(), ntr = dint(TR, "pk_model=='k2016' & set=='iii'", "n_trials", "simulated trials (retained per arm)"),
                               r2i = f_set("i", "r2"), exi = f_set("i", "extrap")))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)

  # ---- 표 ----
  ret <- function(s_, item) { r <- rows(TR, sprintf("set=='%s'", s_)); x <- n_arm - r$retained_median
    dderived(item, TR, sprintf("set=='%s' :: n_per_arm (%s) - retained_median, range over models", s_, fnum(n_arm, 0)), x, rng_fmt(min(x), max(x), 0)) }
  gm <- function(s_, item) drange(TC, sprintf("set=='%s'", s_), "true_aucinf_gmr", 2, "", item)
  cnt <- function(rel, cf) dcount(rel, w1(cf, " & pass_pct > 5"), sprintf("%s cells above 5%% (point), M1", cf))
  mx <- function(rel, cf) dext(rel, w1(cf), "pass_pct", max, 2, "%", sprintf("%s largest boundary type I error, M1", cf))
  U <- L$table$units
  RW <- list(P2 = list(T1, "P2", U$na, U$na), F3B = list(T1, "F3B", fill(U$lost_inf, list(k = ret("lambda", "subjects per arm without an estimable lambda-z (AUCinf only)"))), U$nc),
             Aiii = list(CG, "G2_A_iii", fill(U$lost, list(k = ret("iii", "subjects per arm without a reliable AUCinf, set (iii)"))), fill(U$gmr, list(g = gm("iii", "true AUCinf ratio failing to retained, set (iii)")))),
             B = list(T1, "G2_B", fill(U$lost, list(k = ret("lambda", "subjects per arm without an estimable lambda-z"))), U$nc),
             Ai = list(T1, "G2_Ai", fill(U$lost, list(k = ret("i", "subjects per arm failing set (i)"))), fill(U$gmr, list(g = gm("i", "true AUCinf ratio failing to retained, set (i)")))))
  df <- data.frame(a = unlist(L$table$rows[names(RW)]), b = vapply(RW, function(z) cnt(z[[1]], z[[2]]), ""), c = vapply(RW, function(z) mx(z[[1]], z[[2]]), ""),
                   d = vapply(RW, `[[`, "", 3), e = vapply(RW, `[[`, "", 4), stringsAsFactors = FALSE)
  names(df) <- tx("A10.table.head", f)
  # 높이: 표 추정 높이를 먼저 구해 카드에 나머지를 준다
  TWD <- c(4.55, 1.35, 1.2, 2.2, 2.93)
  th <- local({ w <- TWD / sum(TWD) * GEO$CW; nl <- function(v, w_, b = FALSE) vapply(as.character(v), function(s) est_lines(nobreak(s), w_ - 0.14, 14, b), 1L)
    hdr <- max(mapply(function(v, w_) max(nl(v, w_, TRUE)), names(df), w)); bod <- apply(matrix(sapply(seq_along(w), function(j) nl(df[[j]], w[j], j == 1)), nrow = nrow(df)), 1, max)
    (hdr + sum(bod)) * 14 * 1.2 / 72 + (nrow(df) + 1) * 8 / 72 + 0.04 })
  ty <- capy - 0.10 - th
  deck_table(df, box = c(GEO$ML, ty, GEO$CW, th), widths = TWD, size = 14, highlight = 3:5, highlight_fill = PAL$tint_orange, label = "table_configs")
  dsrc("configuration table", c(T1, CG, TR, TC), "(table)")

  # ---- 카드 세 개 ----
  gap <- 0.2; cw <- (GEO$CW - 2 * gap) / 3; ch <- ty - 0.18 - y0
  cost <- dext(FC, "TRUE", "cost_iii_pp", max, 2, "", "largest loss in joint pass rate when AUCinf (rule B) is added, near-equivalent products, 2016 model")
  C <- L$cards; fills <- c(PAL$tint_blue, PAL$tint_grey, PAL$tint_orange); hc <- c(PAL$blue, PAL$ink, PAL$orange)
  cf_k <- local({ k <- c(nrow(rows(CG, w1("G2_A_iii", " & pass_pct > 5"))), nrow(rows(T1, w1("G2_B", " & pass_pct > 5"))), nrow(rows(T1, w1("G2_Ai", " & pass_pct > 5"))))
    dderived("AUCinf + Cmax cells above 5% (point), M1, rule A set (iii), rule B, rule A set (i): range", CG, "analysis_model=='M1' & config in (G2_A_iii [criteria file], G2_B, G2_Ai [type1_models]) & pass_pct > 5 :: range of row counts", k, rng_fmt(min(k), max(k), 0)) })
  vals <- list(list(), list(cost = cost, fb = f$fb), list(k = cf_k, n = f$n, nom = f$nom))
  for (i in 1:3) {
    x <- GEO$ML + (i - 1) * (cw + gap); cd <- C[[i]]
    ps <- c(list(para(cd$head, SZ$body, hc[i], TRUE, "left", gap_pt = 6)), lapply(cd$lines, function(z) para(fill(z, vals[[i]]), SZ$body, PAL$ink, FALSE, "left", gap_pt = 4)))
    fit_check(sprintf("card_%d", i), c(cd$head, vapply(cd$lines, function(z) fill(z, vals[[i]]), "")), c(x, y0, cw, ch), SZ$body, gap_pt = 5, card = TRUE)
    DK$x <- ph_with(DK$x, do.call(block_list, ps), location = loc(c(x, y0, cw, ch), sprintf("card_%d", i), bg = fills[i], geom = "roundRect", ln = no_line()))
  }

  # ---- 노트 ----
  fcv <- function(s_, col, d, item) dv(FC, sprintf("scenario=='%s'", s_), col, d, "", item)
  pdiff <- local({ r <- rows(PW, "scenario %in% c('S00','F097') & analysis_model=='M1' & config %in% c('P2','F3B')")
    w <- dcast(r, pk_model + scenario ~ config, value.var = "pass_pct"); x <- w$P2 - w$F3B; premise(nrow(w) == 4 && all(x >= 0), "M1 power: F3B at most P2 in the four cells")
    dderived("P2 minus F3B power, S00 and F097, both models, M1 (points)", PW, "scenario in (S00, F097) & analysis_model=='M1' & config in (P2, F3B) :: range of P2 - F3B", x, rng_fmt(min(x), max(x), 2)) })
  deck_notes(tx("A10.notes", c(f, list(
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), r2i = f_set("i", "r2"),
    fb0 = dext(T1, "analysis_model=='M0' & config=='F3B'", "pass_pct", max, 2, "%", "largest boundary type I error, F3B, M0"),
    fbci = local({ r <- rows(T1, w1("F3B"))[which.max(pass_pct)]; dci(T1, w1("F3B", sprintf(" & pk_model=='%s' & scenario=='%s'", r$pk_model, r$scenario)), "pass_pct", "lo", "hi", 2, "%", "F3B largest cell with interval, M1") }),
    fa = dext(T1, w1("F3A"), "pass_pct", max, 2, "%", "largest boundary type I error, F3A (rule A set ii), M1"),
    fcx = dext(T1, w1("F3C"), "pass_pct", max, 2, "%", "largest boundary type I error, F3C (rule C set ii), M1"),
    fcn = dcount(T1, w1("F3C", " & pass_pct > 5"), "F3C cells above 5% (point), M1"),
    p0 = fcv("S00", "i_last_cmax", 2, "joint pass AUClast + Cmax, identical products (%)"), p3 = fcv("S00", "iii_plus_inf_all", 2, "joint pass with AUCinf rule B, identical products (%)"),
    c0 = dci(FC, "scenario=='S00'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 3, "", "loss adding AUCinf rule B, identical products (points)"),
    ck = dci(FC, "scenario=='KE110'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 2, "", "loss adding AUCinf rule B, KE110 (points)"),
    cf9 = dci(FC, "scenario=='F097'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 2, "", "loss adding AUCinf rule B, F097 (points)"),
    rk = fcv("KE110", "true_ratio", 3, "true AUCinf ratio KE110"), rf = fcv("F097", "true_ratio", 3, "true AUCinf ratio F097"),
    ca = drange(FC, "scenario %in% c('KE110','F097')", "cost_ii_pp", 2, "", "loss adding AUCinf rule A (set ii), KE110 and F097 (points)"),
    ntrf = drange(FC, "TRUE", "n_trials", 0, "", "trials per product in the fallback cost set"),
    ep = dv(EP, "model=='k2016' & scenario=='F097'", "power_last_cmax", 2, "%", "empirical power AUClast + Cmax, F097, 2016 model"),
    pd = pdiff,
    gb = dext(T1, w1("G2_B"), "pass_pct", max, 2, "%", "largest boundary type I error, G2_B, M1"), kb = dcount(T1, w1("G2_B", " & pass_pct > 5"), "G2_B cells above 5% (point), M1"), gbl = dcount(T1, w1("G2_B", " & lo > 5"), "G2_B cells with Wilson lower bound above 5%, M1"),
    gi = dext(CG, w1("G2_A_iii"), "pass_pct", max, 2, "%", "largest boundary type I error, G2_A_iii, M1"), gil = dcount(CG, w1("G2_A_iii", " & lo > 5"), "G2_A_iii cells with Wilson lower bound above 5%, M1"),
    fail3 = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), two models"),
    fail1 = drange(TPF, "set=='i'", "fail_pct", 1, "%", "share failing set (i), two models"),
    lz = drange(TPF, "set=='iii'", "lz_pct", 2, "%", "share without an estimable lambda-z, two models"),
    wd = drange(TC, "set=='iii'", "wt_diff_kg", 1, "", "body weight difference failing minus retained, set (iii), two models (kg)"),
    ret3 = drange(TR, "set=='iii'", "retained_median", 0, "", "median subjects per arm with a reliable AUCinf, two models"),
    n_arm = f_n_arm()))))
  deck_end()
}
