
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
  fc <- rows(FC, "TRUE"); premise(all(fc$cost_iii_pp >= 0) && setequal(fc$scenario, c("S00", "KE110", "F097")) && all(fc$true_ratio > 0.94), "fallback cost: three near-equivalent products, non-negative loss (card: 'at most')")
  ep <- rows(EP, "model=='k2016'"); premise(all(abs(ep$power_last_cmax - fc[match(ep$scenario, fc$scenario), i_last_cmax]) < 1e-9), "empirical power and fallback cost come from the same trials (AUClast + Cmax column identical)")
  for (cf in c("G2_A_iii", "G2_B", "G2_Ai")) { r <- rows(if (cf == "G2_A_iii") CG else T1, w1(cf)); premise(sum(r$pass_pct > 5) >= 8, sprintf("%s: many cells above 5%% (card)", cf)) }
  # 대응안이 규칙 A가 아니라 규칙 B인 이유(카드 2, 표): 규칙 A(세트 (ii) 미달자 제외)는 1종 오류는 낮지만 통과율 손실이 크고 노출이 낮은 대상자를 뺀다
  premise(identical(unlist(OC$F3A$endpoints), c("AUClast", "Cmax", "AUCinf_A")) && grepl("span ratio < 2", .read("config/oc_design.yaml")$aucinf_rules$A, fixed = TRUE), "F3A = AUClast + Cmax + AUCinf rule A with the set (ii) flags")
  premise(nrow(rows(T1, w1("F3A", " & pass_pct > 5"))) == 0, "F3A has no cell above 5% under M1")
  premise(all(fc$cost_ii_pp > fc$cost_iii_pp) && all(rows(TC, "set=='ii'")$true_aucinf_gmr_hi < 1), "rule A loses more joint passes than rule B and drops lower-exposure subjects (set ii, both models)")
  premise(nrow(rows("oc_curves/oc_type1_summary.csv", "config=='F3A_iii' & scope=='both'")) == 1, "three-endpoint configuration with rule A set (iii) computed in the v1.2 operating characteristics (notes)")
  premise(.read("config/trial_design.yaml")$be$method == "pooled_t", "fallback cost trials analysed with the pooled t-test (M0)")
  LPr <- .read("config/literature_precedents.yaml")
  premise(LPr$ema_2012_mab$primary_endpoint_single_dose == "AUC0-inf" && LPr$ema_2012_mab$subcutaneous_co_primary == "Cmax" && LPr$ema_2012_mab$status == "search excerpt",
          "EMA 2012 mAb guideline (search excerpt): AUC0-inf primary in single-dose studies, Cmax co-primary for subcutaneous use")
  premise(LPr$msb11456_iv$primary_endpoint == "AUC0-last" && grepl("intravenous", LPr$msb11456_iv$route) && LPr$msb11456_iv$status == "search excerpt", "MSB11456 (search excerpt): single intravenous dose, AUC0-last primary")

  # 대응안을 기본으로 두지 않는 이유(설명 상자·캡션·노트, 검토 2차): (1) 세 지표 판정은 AUClast + Cmax 판정의 부분집합이라 칸마다 통과율이 같거나 낮다(교집합);
  # AUCinf(규칙 B) + Cmax만으로는 여러 칸이 5%를 넘으므로 낮은 1종 오류는 AUCinf 정보가 아니라 지표 하나를 더 통과해야 하는 데서 온다. (2) 규칙 B는 λz가 산출되면
  # 세트 (iii) 미달이어도 AUCinf를 쓴다. (3) AUClast + Cmax의 5% 초과 한 칸은 AUClast 기하평균비가 참 비보다 1 쪽으로 치우친 2020 모델 V2 칸(별첨 A5b)
  pf3 <- merge(rows(T1, w1("P2"))[, .(pk_model, scenario, p2 = pass_pct)], rows(T1, w1("F3B"))[, .(pk_model, scenario, f3 = pass_pct)], by = c("pk_model", "scenario"))
  premise(nrow(pf3) == 16 && all(pf3$f3 <= pf3$p2 + 1e-9), "three-endpoint fallback passes no more often than AUClast + Cmax in every cell (intersection rule)")
  premise(nrow(rows(T1, w1("G2_B", " & pass_pct > 5"))) >= 8, "AUCinf (rule B) + Cmax alone exceeds 5% in many cells: the fallback's low error is not AUCinf information")
  p2x <- rows(T1, w1("P2", " & pass_pct > 5")); DEf <- "oc_models/p2_decomposition_models.csv"
  premise(nrow(p2x) == 1 && p2x$pk_model == "k2020" && p2x$scenario == "V2_up_080", "the one AUClast + Cmax cell above 5% is the 2020 V2 cell (M1)")
  de1 <- row1(DEf, "pk_model=='k2020' & mechanism=='V2' & analysis_model=='M1'")
  premise(de1$auclast_bias_dir == "toward_1" && de1$auclast_bias_lo_pct > 0 && abs(de1$p2_pct - p2x$pass_pct) < 1e-9, "2020 V2 cell: AUClast GMR biased toward 1 (interval above 0), same cell and value as type1_models")
  premise(!grepl("Rsq|span|Extrap", .read("config/oc_design.yaml")$aucinf_rules$B) && grepl("span ratio", .read("config/oc_design.yaml")$aucinf_rules$A, fixed = TRUE),
          "rule B applies no reliability flag: every subject with an estimable lambda-z is kept, including those failing set (iii)")

  # ---- 제목 ----
  f <- list(fb = dext(T1, w1("F3B"), "pass_pct", max, 2, "%", "largest boundary type I error, fallback AUClast + AUCinf (rule B) + Cmax, M1"), nom = f_nominal(),
            n = dcount(T1, w1("P2"), "boundary cells"))
  y0 <- core_title(tx("A10.title", f), tx("A10.kicker"))

  # ---- 캡션 ----
  fcr <- drange(FC, "TRUE", "true_ratio", 2, "", "true AUCinf ratio of the three fallback-cost products")
  pdiff <- local({ r <- rows(PW, "scenario %in% c('S00','F097') & analysis_model=='M1' & config %in% c('P2','F3B')")
    w <- dcast(r, pk_model + scenario ~ config, value.var = "pass_pct"); x <- w$P2 - w$F3B; premise(nrow(w) == 4 && all(x >= 0), "M1 power: F3B at most P2 in the four cells")
    dderived("P2 minus F3B power, S00 and F097, both models, M1 (points)", PW, "scenario in (S00, F097) & analysis_model=='M1' & config in (P2, F3B) :: range of P2 - F3B", x, rng_fmt(min(x), max(x), 2)) })
  wy <- list(p2m = dv(T1, w1("P2", " & pk_model=='k2020' & scenario=='V2_up_080'"), "pass_pct", 2, "%", "AUClast + Cmax, 2020 V2 cell, M1"),
             bias = dv(DEf, "pk_model=='k2020' & mechanism=='V2' & analysis_model=='M1'", "auclast_bias_pct", 2, "%", "AUClast GMR bias toward 1 against the true AUCinf ratio, 2020 V2 cell, M1"),
             efl = drange(TPF, "set=='iii'", "est_fail_pct", 0, "%", "subjects with an estimable lambda-z but failing set (iii), two models"))
  cap <- tx("A10.caption", c(wy, list(nom = f$nom, n = f$n, n_arm = f_n_arm(), r2i = f_set("i", "r2"), exi = f_set("i", "extrap"), sp2 = f_set("ii", "span"), fcr = fcr,
                                ey = dcfg("literature_precedents.yaml", c("ema_2012_mab", "year"), "EMA mAb biosimilar guideline year", num_fmt(0)))))
  capy <- core_caption(cap, GEO$BODY_BOTTOM, size = 14)

  # ---- 표 ----
  ret <- function(s_, item) { r <- rows(TR, sprintf("set=='%s'", s_)); x <- n_arm - r$retained_median
    dderived(item, TR, sprintf("set=='%s' :: n_per_arm (%s) - retained_median, range over models", s_, fnum(n_arm, 0)), x, rng_fmt(min(x), max(x), 0)) }
  gm <- function(s_, item) drange(TC, sprintf("set=='%s'", s_), "true_aucinf_gmr", 2, "", item)
  cnt <- function(rel, cf) dcount(rel, w1(cf, " & pass_pct > 5"), sprintf("%s cells above 5%% (point), M1", cf))
  mx <- function(rel, cf) dext(rel, w1(cf), "pass_pct", max, 2, "%", sprintf("%s largest boundary type I error, M1", cf))
  U <- L$table$units
  cst <- function(col, item) drange(FC, "TRUE", col, 2, "%p", item)
  RW <- list(P2 = list(T1, "P2", U$na, U$na, U$ref), F3B = list(T1, "F3B", fill(U$lost_inf, list(k = ret("lambda", "subjects per arm without an estimable lambda-z (AUCinf only)"))), U$nc,
                                                        cst("cost_iii_pp", "loss in joint pass rate adding AUCinf rule B, three near-equivalent products, 2016 model, M0")),
             F3A = list(T1, "F3A", fill(U$lost_inf, list(k = ret("ii", "subjects per arm failing set (ii) (AUCinf only)"))), fill(U$gmr, list(g = gm("ii", "true AUCinf ratio failing to retained, set (ii)"))),
                        cst("cost_ii_pp", "loss in joint pass rate adding AUCinf rule A (set ii), three near-equivalent products, 2016 model, M0")),
             Aiii = list(CG, "G2_A_iii", fill(U$lost, list(k = ret("iii", "subjects per arm without a reliable AUCinf, set (iii)"))), fill(U$gmr, list(g = gm("iii", "true AUCinf ratio failing to retained, set (iii)"))), U$nc),
             B = list(T1, "G2_B", fill(U$lost, list(k = ret("lambda", "subjects per arm without an estimable lambda-z"))), U$nc, U$nc),
             Ai = list(T1, "G2_Ai", fill(U$lost, list(k = ret("i", "subjects per arm failing set (i)"))), fill(U$gmr, list(g = gm("i", "true AUCinf ratio failing to retained, set (i)"))), U$nc))
  df <- data.frame(a = unlist(L$table$rows[names(RW)]), b = vapply(RW, function(z) cnt(z[[1]], z[[2]]), ""), c = vapply(RW, function(z) mx(z[[1]], z[[2]]), ""),
                   f = vapply(RW, `[[`, "", 5), d = vapply(RW, `[[`, "", 3), e = vapply(RW, `[[`, "", 4), stringsAsFactors = FALSE)
  names(df) <- tx("A10.table.head", f)
  # 높이: 표 추정 높이를 먼저 구해 카드에 나머지를 준다
  TWD <- c(3.9, 1.35, 0.95, 1.45, 2.2, 2.38)
  th <- deck_table_h(df, GEO$CW, TWD, 14, pad = 2) + 0.04
  ty <- capy - 0.10 - th
  deck_table(df, box = c(GEO$ML, ty, GEO$CW, th), widths = TWD, size = 14, highlight = 4:6, label = "table_configs", pad = 2)   # 셀 위아래 여백 2pt: 카드 세 개 + 표 + 캡션을 한 장에
  dsrc("configuration table", c(T1, CG, TR, TC, FC), "(table)")

  # ---- 설명 상자: 대응안을 기본으로 두지 않는 이유, AUCinf + Cmax만은 권하지 않음(표의 주황 행) ----
  cf_k <- local({ k <- c(nrow(rows(CG, w1("G2_A_iii", " & pass_pct > 5"))), nrow(rows(T1, w1("G2_B", " & pass_pct > 5"))), nrow(rows(T1, w1("G2_Ai", " & pass_pct > 5"))))
    dderived("AUCinf + Cmax cells above 5% (point), M1, rule A set (iii), rule B, rule A set (i): range", CG, "analysis_model=='M1' & config in (G2_A_iii [criteria file], G2_B, G2_Ai [type1_models]) & pass_pct > 5 :: range of row counts", k, rng_fmt(min(k), max(k), 0)) })
  # 대응안의 2종 오류 비용(v1.2, 지시 2026-09-29 "S9 재설계" §5: A10은 F3-B 대응안 유지, 2종 오류 비용 추가): S7과 같은 시험(results/oc_curves, M1)
  PF <- "oc_curves/oc_curves_pass.csv"; PDc <- "oc_curves/oc_paired.csv"
  t2i <- function(cf) { r <- rows(PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s'", cf)); x <- range(100 - r$pass_pct)
    dderived(sprintf("type II error, identical product, %s, M1, range over two models", cf), PF, sprintf("code=='S00' & analysis_model=='M1' & config=='%s' :: range(100 - pass_pct)", cf), x, rng_fmt(x[1], x[2], 1, "%")) }
  W2 <- "comparison=='F3B - P2' & (kind=='identity' | (cmax_in_limits & (abs(target - 0.95) < 1e-9 | abs(target - 1.05) < 1e-9 | abs(target - 0.90) < 1e-9 | abs(target - 1.11) < 1e-9)))"
  i2 <- rows(PDc, W2); premise(nrow(i2) > 0 && max(-i2$diff_pass_pp) < max(-rows(PDc, sub("F3B", "F3A_iii", W2))$diff_pass_pp), "fallback (rule B) type II increase smaller than with set (iii) AUCinf (why box)")
  f3binc <- dderived("largest type II error increase, F3B vs P2, type II cells, both models, M1 (percentage points)", PDc, paste(W2, ":: max(-diff_pass_pp)"), max(-i2$diff_pass_pp), paste0(fnum(max(-i2$diff_pass_pp), 1), "%p"))
  t2c <- list(f3bid = t2i("F3B"), p2id = t2i("P2"), f3aid = t2i("F3A_iii"), f3binc = f3binc)
  why <- tx("A10.why", c(wy, t2c, list(k = cf_k, n = f$n, nom = f$nom))); wyy <- y0 + 0.02; wh <- ty - 0.14 - wyy
  deck_text(why, c(GEO$ML, wyy, GEO$CW, wh), size = 16, bg = PAL$tint_grey, geom = "roundRect", label = "caption_why", gap_pt = 4)
  deck_visual(c(GEO$ML, wyy, GEO$CW, wh))

  # ---- 노트 ----
  fcv <- function(s_, col, d, item) dv(FC, sprintf("scenario=='%s'", s_), col, d, "", item)
  deck_notes(tx("A10.notes", c(f, list(
    r2 = f_set("iii", "r2"), ex = f_set("iii", "extrap"), r2i = f_set("i", "r2"),
    fb0 = dext(T1, "analysis_model=='M0' & config=='F3B'", "pass_pct", max, 2, "%", "largest boundary type I error, F3B, M0"),
    fbci = local({ r <- rows(T1, w1("F3B"))[which.max(pass_pct)]; dci(T1, w1("F3B", sprintf(" & pk_model=='%s' & scenario=='%s'", r$pk_model, r$scenario)), "pass_pct", "lo", "hi", 2, "%", "F3B largest cell with interval, M1") }),
    fa = dext(T1, w1("F3A"), "pass_pct", max, 2, "%", "largest boundary type I error, F3A (rule A set ii), M1"),
    fcx = dext(T1, w1("F3C"), "pass_pct", max, 2, "%", "largest boundary type I error, F3C (rule C set ii), M1"),
    fcn = dcount(T1, w1("F3C", " & pass_pct > 5"), "F3C cells above 5% (point), M1"),
    p0 = paste0(fcv("S00", "i_last_cmax", 2, "joint pass AUClast + Cmax, identical products (%)"), "%"), p3 = paste0(fcv("S00", "iii_plus_inf_all", 2, "joint pass with AUCinf rule B, identical products (%)"), "%"),
    c0 = dci(FC, "scenario=='S00'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 3, "%p", "loss adding AUCinf rule B, identical products (points)"),
    ck = dci(FC, "scenario=='KE110'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 2, "%p", "loss adding AUCinf rule B, KE110 (points)"),
    cf9 = dci(FC, "scenario=='F097'", "cost_iii_pp", "cost_iii_lo", "cost_iii_hi", 2, "%p", "loss adding AUCinf rule B, F097 (points)"),
    rk = fcv("KE110", "true_ratio", 3, "true AUCinf ratio KE110"), rf = fcv("F097", "true_ratio", 3, "true AUCinf ratio F097"),
    ca = drange(FC, "scenario %in% c('KE110','F097')", "cost_ii_pp", 2, "", "loss adding AUCinf rule A (set ii), KE110 and F097 (points)"),
    ntrf = local({ x <- range(rows(FC, "TRUE")$n_trials); dderived("trials per product in the fallback cost set", FC, "TRUE :: range(n_trials)", x, sprintf("%s~%s", fint(x[1]), fint(x[2]))) }),
    ep = dv(EP, "model=='k2016' & scenario=='F097'", "power_last_cmax", 2, "%", "empirical power AUClast + Cmax, F097, 2016 model"),
    pd = pdiff,
    gb = dext(T1, w1("G2_B"), "pass_pct", max, 2, "%", "largest boundary type I error, G2_B, M1"), kb = dcount(T1, w1("G2_B", " & pass_pct > 5"), "G2_B cells above 5% (point), M1"), gbl = dcount(T1, w1("G2_B", " & lo > 5"), "G2_B cells with Wilson lower bound above 5%, M1"),
    gi = dext(CG, w1("G2_A_iii"), "pass_pct", max, 2, "%", "largest boundary type I error, G2_A_iii, M1"), gil = dcount(CG, w1("G2_A_iii", " & lo > 5"), "G2_A_iii cells with Wilson lower bound above 5%, M1"),
    fail3 = drange(TPF, "set=='iii'", "fail_pct", 1, "%", "share without a reliable AUCinf, set (iii), two models"),
    fail1 = drange(TPF, "set=='i'", "fail_pct", 1, "%", "share failing set (i), two models"),
    lz = drange(TPF, "set=='iii'", "lz_pct", 2, "%", "share without an estimable lambda-z, two models"),
    ret3 = drange(TR, "set=='iii'", "retained_median", 0, "", "median subjects per arm with a reliable AUCinf, two models"),
    ntr = dint(TR, "pk_model=='k2016' & set=='iii'", "n_trials", "simulated trials (retained per arm)"), n_arm = f_n_arm(),
    p2m = wy$p2m, bias = wy$bias, efl = wy$efl, nom = f$nom, fb = f$fb, sp2 = f_set("ii", "span"),
    bci = dci(DEf, "pk_model=='k2020' & mechanism=='V2' & analysis_model=='M1'", "auclast_bias_pct", "auclast_bias_lo_pct", "auclast_bias_hi_pct", 2, "%", "AUClast GMR bias with 95% interval, 2020 V2 cell, M1"),
    p2x = dcount(T1, w1("P2", " & pass_pct > 5"), "AUClast + Cmax cells above 5% (point), M1"),
    f3a3 = dv("oc_curves/oc_type1_summary.csv", "config=='F3A_iii' & scope=='both'", "max_pct", 1, "%", "largest boundary type I error over 16 cells, F3A_iii, M1")), t2c)))
  deck_end()
}
