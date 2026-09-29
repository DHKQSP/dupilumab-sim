# S3 질문과 방법 흐름도(지시 §2): 모델 2개 → 계획 시험 모의(60~90 kg, Syneos 채혈 13회) → Phoenix 규칙 NCA → 모델 참값과 비교.
# 모델 적격성 한 줄(연구 제형 300 mg 관측 AUClast ±15%, 별첨 A1): step1b_quant_gate.csv의 gate 역할 행 비가 두 모델 모두 허용 범위 안인지 전제로 확인한다.
slide_S3 <- function() {
  G16 <- "step1/step1b_quant_gate.csv"; G20 <- "step1_k2020/step1b_quant_gate.csv"
  deck_slide("S3", tag = "sim")
  tol <- dcfg("design_clot2021.yaml", c("gate", "auclast_mean_tol_pct"), "gate tolerance for the simulated AUClast mean (%)", num_fmt(0))
  tolv <- .read("config/design_clot2021.yaml")$gate$auclast_mean_tol_pct
  for (g in c(G16, G20)) { r <- rows(g, "gate_role=='gate' & dose_mg==300"); premise(nrow(r) >= 1 && all(abs(r$AUClast_ratio - 1) * 100 <= tolv), sprintf("%s: every 300 mg gate row within the tolerance (text)", g)) }
  rng <- dspan(G16, "gate_role=='gate' & dose_mg==300", "AUClast_ratio", "AUClast_ratio", 2, "", "simulated / observed AUClast mean, 300 mg gate rows, 2016 model")
  rng20 <- dspan(G20, "gate_role=='gate' & dose_mg==300", "AUClast_ratio", "AUClast_ratio", 2, "", "simulated / observed AUClast mean, 300 mg gate rows, 2020 model")
  y0 <- core_title(tx("S3.title"), tx("S3.kicker"))

  # 질문 띠
  deck_text(tx("S3.question"), c(GEO$ML, y0, GEO$CW, 0.56), size = 20, bold = TRUE, label = "body_question", bg = PAL$tint_grey)
  # ---- 흐름도(주 시각 요소): 단계 4개를 위에서 아래로(왼쪽 머리말, 오른쪽 내용), 단계 사이 아래 화살표 ----
  premise(identical(.read("config/literature_precedents.yaml")$ema_2012_mab$primary_endpoint_single_dose, "AUC0-inf") && identical(.read("config/literature_precedents.yaml")$ema_2012_mab$subcutaneous_co_primary, "Cmax"),
          "EMA 2012 mAb guideline excerpt: AUC0-inf primary, Cmax co-primary for SC (body)")
  body <- tx("S3.body", list(tol = tol, dose = f_dose(), yr = dcfg("literature_precedents.yaml", c("ema_2012_mab", "year"), "EMA mAb biosimilar guideline year", num_fmt(0))))
  B <- DK$txt$S3$boxes; top <- y0 + 0.56 + 0.20; bottom <- GEO$BODY_BOTTOM - core_body_h(body) - 0.10; ag <- 0.16
  rh <- (bottom - top - 3 * ag) / 4; hw <- 3.05; cx <- GEO$ML + hw + 0.12; cwid <- GEO$CW - hw - 0.12
  f <- list(wt = f_wt_range(), dose = f_dose(), nB0 = f_study_days("B0", "n"), last = f_study_days("B0", "last"), n_arm = f_n_arm(), lloq = f_lloq(),
            minpts = dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0)))
  for (k in 1:4) {
    y <- top + (k - 1) * (rh + ag); b <- B[[k]]
    deck_text(b$head, c(GEO$ML, y, hw, rh), size = 20, bold = TRUE, color = if (k == 4) PAL$blue else PAL$ink, label = "cmid_diag_head", bg = if (k == 4) PAL$tint_blue else PAL$tint_grey)
    deck_text(fill(b$text, f), c(cx, y, cwid, rh), size = SZ$body, label = "cmid_diag_text", bg = if (k == 4) PAL$tint_blue else PAL$tint_grey)
    if (k < 4) deck_box(c(GEO$ML + hw / 2 - 0.16, y + rh + 0.01, 0.32, ag - 0.02), fill = PAL$muted, geom = "downArrow", label = "arrow")
  }
  deck_visual(c(GEO$ML, top, GEO$CW, bottom - top))
  # ---- 본문 한 줄: 모델 적격성 ----
  core_body(body, GEO$BODY_BOTTOM)
  deck_notes(tx("S3.notes", c(f, list(tol = tol, rng = rng, rng20 = rng20))))
  deck_end()
}
