# A3 부록: 모델 파라미터와 출처. 두 PK 모델(2016 모델 주분석, 2020 Model 1 민감도)의 구조 파라미터·개체간 변동·잔차를 config/params_*.yaml에서 dcfg로 읽는다. [문헌]
# regulatory/tables/parameter_provenance.csv는 Value 열이 문자·숫자 혼합이라 수치 출처로 쓰지 않는다(보고서 부록 A의 표시용 표).
# 모델은 속도 상수(ke, k12, k21)로 모수화되어 있어 CL, Q, 말초 용적은 추정값이 아니다. 75 kg 환산값은 dderived로 계산식을 남긴다.
# 출처 열의 논문 표기(Table 2, Table 1, 보충 Table S2)는 문구 파일에 두고, config의 source 필드가 같은 표를 가리키는지 premise로 검사한다.

a3_v <- function(file, path, item, d) dcfg(file, path, item, num_fmt(d))

slide_A3 <- function() {
  PT <- "params_typical.yaml"; PV <- "params_variability.yaml"; P20 <- "params_k2020_model1.yaml"
  deck_slide("A3", tag = "lit")
  L <- DK$txt$A3
  pt <- .read(file.path("config", PT)); pv <- .read(file.path("config", PV)); p20 <- .read(file.path("config", P20))
  core <- c("Vc", "ke", "k12", "k21", "ka", "Vmax", "Km", "F")

  # ---- 전제: 모든 값은 논문 표에서 확정(재추정·보정 없음), 출처 표 번호가 문구와 같다 ----
  premise(all(vapply(core, function(k) identical(pt$theta[[k]]$status, "confirmed") && identical(p20$theta[[k]]$status, "confirmed"), TRUE)) &&
            identical(pv$status, "confirmed") && identical(pt$covariates$WT_on_Vc$theta_WT$status, "confirmed") && identical(p20$covariates$WT_on_Vc$theta_WT$status, "confirmed"),
          "every structural parameter of both models has status 'confirmed' (taken from the publications, not re-estimated)")
  premise(all(vapply(core, function(k) grepl("Kovalenko 2016 Table 2", pt$theta[[k]]$source, fixed = TRUE), TRUE)) && grepl("Kovalenko 2016 Table 2", pv$source, fixed = TRUE),
          "2016 model values and IIV from Kovalenko 2016 Table 2")
  premise(all(vapply(core, function(k) grepl("Kovalenko 2020 Table 1 Model 1", p20$theta[[k]]$source, fixed = TRUE), TRUE)) &&
            grepl("Supplementary Table 2", p20$iiv$source, fixed = TRUE) && grepl("Supplementary Table 2", p20$residual$sigma_prop$source, fixed = TRUE),
          "2020 Model 1 values from Table 1, IIV and residual error from Supplementary Table 2")
  premise(grepl("(V2", pt$theta$Vc$source, fixed = TRUE) && grepl("(k23)", pt$theta$k12$source, fixed = TRUE) && grepl("(k32)", pt$theta$k21$source, fixed = TRUE) &&
            grepl("(kcp)", p20$theta$k12$source, fixed = TRUE) && grepl("kpc", p20$theta$k21$source, fixed = TRUE) && grepl("(Vm)", p20$theta$Vmax$source, fixed = TRUE),
          "paper notation in the source column (2016: V2, k23, k32; 2020: kcp, kpc, Vm)")
  premise(isTRUE(pt$theta$Km$fixed) && isTRUE(p20$theta$Km$fixed), "Km fixed in both models")
  zero16 <- c("k12", "k21", "F", "Km"); zero20 <- c("k12", "k21", "F", "Km")
  premise(all(vapply(zero16, function(k) pv$iiv_omega2[[k]]$omega2 == 0, TRUE)) && all(vapply(zero20, function(k) p20$iiv$sd[[k]]$omega2 == 0, TRUE)),
          "no IIV on k12, k21, F and Km in either model ('none' in the table)")
  premise(all(vapply(c("Vc", "ke", "ka", "Vmax", "MTT"), function(k) abs(p20$iiv$sd[[k]]$sd^2 - p20$iiv$sd[[k]]$omega2) < 5e-4, TRUE)),
          "2020 IIV variances equal the squared SD of Supplementary Table 2")
  premise(pt$covariates$WT_on_Vc$WT_ref$value == p20$covariates$WT_on_Vc$WT_ref$value, "same reference weight in both models")

  # 제목의 '변환은 SD 제곱과 kpc 유도뿐': 두 모델의 파라미터 출처 문구 중 '유도'는 2020 k21(kpc)뿐이고, kpc = kcp / Mpc가 수치로 맞다
  srcs16 <- c(vapply(core, function(k) pt$theta[[k]]$source, ""), pt$covariates$WT_on_Vc$theta_WT$source, vapply(names(pv$iiv_omega2), function(k) pv$iiv_omega2[[k]]$source, ""))
  srcs20 <- setNames(c(vapply(core, function(k) p20$theta[[k]]$source, ""), p20$covariates$WT_on_Vc$theta_WT$source), c(core, "theta_WT"))
  premise(!any(grepl("유도", srcs16)) && identical(names(srcs20)[grepl("유도", srcs20)], "k21"), "the only derived structural value is the 2020 kpc (k21); the 2016 values are taken as published")
  mpc <- as.numeric(regmatches(p20$theta$k21$source, regexec("Mpc ([0-9.]+)", p20$theta$k21$source))[[1]][2])
  premise(is.finite(mpc) && abs(p20$theta$k12$value / mpc - p20$theta$k21$value) < 5e-3, "2020 kpc equals kcp / Mpc (notes)")
  premise(!any(grepl("WT", c(names(pv$iiv_omega2), names(p20$iiv$sd)))), "no IIV on the weight exponent in either model ('not applicable' in the table)")

  wr <- a3_v(PT, c("covariates", "WT_on_Vc", "WT_ref", "value"), "reference body weight for Vc (kg), both models", 0)
  deck_kicker(tx("A3.kicker")); deck_title(tx("A3.title"))

  # ---- 표 하나(두 모델을 나란히; 같은 행이 같은 파라미터라 행 높이가 저절로 맞는다): 파라미터 | 2016 값, 분산, 출처 | 2020 값, 분산, 출처 ----
  # 모델별 열 묶음 위에 모델 이름과 출처 표를 적는다. 머리글은 두 묶음이 같으므로 2020 쪽 이름 끝에 공백을 붙여 열 이름만 구분한다(보이는 글자는 같다).
  none <- L$none
  cols <- function(m) {
    if (m == "k2016") { f <- PT; iv <- function(k, d) a3_v(PV, c("iiv_omega2", k, "omega2"), sprintf("2016 model, IIV variance %s", k), d); nm <- "2016 model" }
    else { f <- P20; iv <- function(k, d) a3_v(P20, c("iiv", "sd", k, "omega2"), sprintf("2020 Model 1, IIV variance %s", k), d); nm <- "2020 Model 1" }
    wexp <- c("covariates", "WT_on_Vc", "theta_WT", "value")
    th <- function(k, d, it) a3_v(f, c("theta", k, "value"), sprintf("%s, %s", nm, it), d)
    dk12 <- if (m == "k2016") 4 else 3
    dvm <- if (m == "k2016") 3 else 2
    val <- c(th("Vc", 2, "Vc (L at 75 kg)"),
             a3_v(f, wexp, sprintf("%s, weight exponent on Vc", nm), 3),
             th("ke", 4, "ke (1/day)"),
             sprintf("%s / %s", th("k12", dk12, "k12 (1/day)"), th("k21", 3, "k21 (1/day)")),
             sprintf("%s / %s", th("ka", 3, "ka (1/day)"), th("F", 3, "F")),
             sprintf("%s / %s", th("Vmax", dvm, "Vmax (mg/L/day)"), th("Km", 2, "Km (mg/L), fixed")))
    # ω²는 유효 숫자 세 자리로 적는다(2020 ka, MTT는 SD 제곱이라 네 자리로 저장됨)
    om <- c(iv("Vc", 4), L$na, iv("ke", if (m == "k2016") 3 else 4), none, sprintf("%s / %s", iv("ka", 3), none), sprintf("%s / %s", iv("Vmax", 4), none))
    df <- data.frame(b = val, c = om, d = unlist(L$src[[m]]), stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- paste0(tx("A3.table.head"), if (m == "k2020") " " else ""); df
  }
  df <- cbind(data.frame(a = tx("A3.rows", list(wr = wr)), stringsAsFactors = FALSE), cols("k2016"), cols("k2020"))
  names(df)[1] <- tx("A3.table.param")
  wd <- c(2.05, 1.55, 1.45, 2.05, 1.55, 1.45, 2.05); wd <- wd / sum(wd) * GEO$CW
  ly <- GEO$BODY_TOP; lh <- 0.68; ty <- ly + lh + 0.02; th_ <- 2.50
  for (k in 1:2) {
    m <- c("k2016", "k2020")[k]; x <- GEO$ML + wd[1] + (k - 1) * sum(wd[2:4])
    deck_text(tx(sprintf("A3.label.%s", m)), c(x, ly, sum(wd[2:4]), lh), size = 16, label = sprintf("label_%s", m), gap_pt = 0)
  }
  deck_table(df, box = c(GEO$ML, ty, GEO$CW, th_), widths = wd, size = 13, label = "table_params")
  # 두 모델 열 묶음 사이 세로 구분선(1 pt): 모델 이름표 위에서 표 아래 테두리까지. deck_table은 세로 테두리를 두지 않으므로 글자 없는 가는 도형으로 그린다.
  # 표 행 높이는 내용으로 정해진다(실측: 한 줄 행 23.5 pt, 두 줄 행 39 pt; 머리글 + 한 줄 행 5개 + 두 줄 행 1개 = 180 pt = 2.50 in). 선 끝은 그 실측 아래 테두리에 맞춘다
  xd <- GEO$ML + sum(wd[1:4]); lw <- 1 / 72; tb <- ty + 180 / 72
  deck_text(" ", c(xd - lw / 2, ly + 0.06, lw, tb - ly - 0.06), size = 16, bg = PAL$muted, geom = "rect", label = "rule_models", gap_pt = 0)

  # ---- 아래: 75 kg 환산(CL, Q, 말초 용적), 2020 흡수 구조, 잔차 ----
  cv <- function(y, file, m) {
    Vc <- y$theta$Vc$value; ke <- y$theta$ke$value; k12 <- y$theta$k12$value; k21 <- y$theta$k21$value
    list(cl = dderived(sprintf("%s, CL = ke x Vc at 75 kg (L/day), derived", m), file.path("config", file), "theta.ke.value x theta.Vc.value", ke * Vc, fnum(ke * Vc, 3)),
         q = dderived(sprintf("%s, Q = k12 x Vc at 75 kg (L/day), derived", m), file.path("config", file), "theta.k12.value x theta.Vc.value", k12 * Vc, fnum(k12 * Vc, 3)),
         vp = dderived(sprintf("%s, peripheral volume = Vc x k12 / k21 at 75 kg (L), derived", m), file.path("config", file), "theta.Vc.value x theta.k12.value / theta.k21.value", Vc * k12 / k21, fnum(Vc * k12 / k21, 2)))
  }
  c16 <- cv(pt, PT, "2016 model"); c20 <- cv(p20, P20, "2020 Model 1")
  premise(pv$residual$sigma_add$value == p20$residual$sigma_add$value && grepl("고정", pv$residual$sigma_add$source) && grepl("고정", p20$residual$sigma_add$source),
          "additive residual SD equal and fixed in both models")
  f <- list(cl16 = c16$cl, cl20 = c20$cl, q16 = c16$q, q20 = c20$q, vp16 = c16$vp, vp20 = c20$vp, wr = wr,
            ntr = a3_v(P20, c("theta", "n_transit", "value"), "2020 Model 1, number of transit compartments", 0),
            mtt = a3_v(P20, c("theta", "MTT", "value"), "2020 Model 1, mean transit time (day)", 3),
            mttw = a3_v(P20, c("iiv", "sd", "MTT", "omega2"), "2020 Model 1, IIV variance MTT", 3),
            sp16 = dderived("2016 model, proportional residual (%)", "config/params_variability.yaml", "residual.sigma_prop.value x 100", 100 * pv$residual$sigma_prop$value, paste0(fnum(100 * pv$residual$sigma_prop$value, 1), "%")),
            sp20 = dderived("2020 Model 1, proportional residual (%)", "config/params_k2020_model1.yaml", "residual.sigma_prop.value x 100", 100 * p20$residual$sigma_prop$value, paste0(fnum(100 * p20$residual$sigma_prop$value, 1), "%")),
            sa = a3_v(PV, c("residual", "sigma_add", "value"), "additive residual SD (mg/L), fixed, both models", 2))
  yb <- ty + th_ + 0.08
  deck_bullets(tx("A3.bullets", f), box = c(GEO$ML, yb, GEO$CW, GEO$BODY_BOTTOM - yb), size = 16, gap_pt = 4)

  deck_notes(tx("A3.notes", list(
    wr = wr, lloq = f_lloq(),
    sd_vc = a3_v(P20, c("iiv", "sd", "Vc", "sd"), "2020 Model 1, IIV SD Vc (Supplementary Table 2)", 3),
    sd_ke = a3_v(P20, c("iiv", "sd", "ke", "sd"), "2020 Model 1, IIV SD ke (Supplementary Table 2)", 3),
    sd_ka = a3_v(P20, c("iiv", "sd", "ka", "sd"), "2020 Model 1, IIV SD ka (Supplementary Table 2)", 3),
    sd_vm = a3_v(P20, c("iiv", "sd", "Vmax", "sd"), "2020 Model 1, IIV SD Vmax (Supplementary Table 2)", 3),
    sd_mtt = a3_v(P20, c("iiv", "sd", "MTT", "sd"), "2020 Model 1, IIV SD MTT (Supplementary Table 2)", 3),
    om_mult = a3_v(PV, c("sensitivity", "omega2_multiplier", "value"), "variability sensitivity: multiplier on all IIV variances", 1),
    sp_alt = dderived("residual sensitivity: alternative proportional residual (%)", "config/params_variability.yaml", "sensitivity.sigma_prop_alt.value x 100",
                      100 * pv$sensitivity$sigma_prop_alt$value, paste0(fnum(100 * pv$sensitivity$sigma_prop_alt$value, 0), "%")))))
  deck_end()
}
