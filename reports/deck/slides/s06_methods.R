# S06 방법 개요: 모델 -> 대상자·채혈(B0) -> 비구획 분석(NCA) -> 가상 시험(TOST, M1) -> 경계 1종 오류(사전 등록) 흐름도,
# NCA 엔진 검증 표(results/nca_engine/engine_validation_summary.csv), 사전 등록 상자. 결과 수치는 다른 슬라이드에서 다룬다.
# B0 채혈일: config는 투여 후 일, 연구일 = 투여 후 일 + 1. (공용 f_study_days는 YAML 목록이 정수·실수 혼합이라 d + 1에서 멈추므로 여기서 unlist)
s06_days <- function(which = c("n", "last")) {
  which <- match.arg(which); d <- unlist(.read("config/trial_design.yaml")$schedules$B0$days); premise(length(d) > 0, "schedule B0")
  sd <- d + 1; p <- switch(which, n = as.character(length(d)), last = fnum(max(sd), 0))
  dderived(sprintf("schedule B0, %s (study day = days after dose + 1)", which), "config/trial_design.yaml", sprintf("schedules.B0.days :: %s", which), sd, p)
}
# 사전 등록 커밋(regulatory/tables/prespecification_register.csv의 Date (evidence) 열에서 읽음)
s06_commit <- function(item_rx, item) {
  PR <- "regulatory/tables/prespecification_register.csv"; w <- sprintf("grepl('%s', Item)", item_rx)
  r <- row1(PR, w); m <- regmatches(r[["Date (evidence)"]], regexec("commit ([0-9a-f]{7})", r[["Date (evidence)"]]))[[1]]
  premise(length(m) == 2, sprintf("commit hash in register row [%s]", w))
  dderived(item, PR, sprintf("%s :: Date (evidence), regex commit ([0-9a-f]{7})", w), m[2], m[2])
}
slide_S06 <- function() {
  EV <- "nca_engine/engine_validation_summary.csv"; T1 <- "oc_models/type1_models.csv"; EX <- "oc_models/extension_decision_models.csv"
  deck_slide("S06", tag = "litsim")
  f <- list(n_arm = f_n_arm(), n_rand = f_n_rand(), wt = f_wt_range(), split = f_split(), lim = f_limits(), ci = f_ci_level(),
            nb0 = s06_days("n"), last = s06_days("last"),
            ntr = dcfg("params_k2020_model1.yaml", c("theta", "n_transit", "value"), "2020 model: transit compartments", num_fmt(0)),
            minpts = dcfg("nca_rules.yaml", c("standard", "lambda_z", "min_points"), "lambda-z minimum points", num_fmt(0)),
            b = dcfg("oc_design.yaml", "boundary_targets", "boundary true AUC0-inf ratios", function(x) paste(fnum(x, 2), collapse = ", ")),
            ncell = dcount(T1, "analysis_model=='M1' & config=='P2'", "boundary cells (two PK models)"),
            reps = f_reps("boundary"), ext = f_reps("ext"),
            next_ = dcount(EX, "analysis_model=='M0' & selected==TRUE", "boundary cells extended"))
  names(f)[names(f) == "next_"] <- "next"
  km_fixed <- vapply(c("config/params_typical.yaml", "config/params_k2020_model1.yaml"), function(p_) isTRUE(.read(p_)$theta$Km$fixed), TRUE)
  premise(all(km_fixed), "Km is fixed in both models (card 1)")
  deck_kicker(tx("S06.kicker")); deck_title(tx("S06.title", f))

  # ---- 단계 카드 5개(①~⑤): 제목 + 설명. 카드 너비는 글 양에 맞춰 나누고(가중치), 높이는 가장 긴 카드의 추정 높이에 맞춘다 ----
  # ③은 식이 "AUC0-inf = AUC0-last" / "+ Clast/λz" 두 줄로 나뉘는 너비(식 전체 한 줄은 카드 두 개 폭이 필요), ①·②는 가장 긴 낱말 덩어리가 들어가는 최소 폭.
  # 높이 추정은 렌더링의 실제 글줄 폭(카드 너비 - 0.27 in)에 맞춰 0.04 in 좁게 잰다(추정 0.24 in 여백보다 렌더링이 조금 좁다)
  S <- DK$txt$S06$steps; n <- length(S); gap <- 0.12; y0 <- GEO$BODY_TOP; cgap <- 6
  wt_ <- c(2.1, 2.0, 2.6, 2.5, 2.55); w <- wt_ / sum(wt_) * (GEO$CW - (n - 1) * gap)
  paras <- lapply(seq_len(n), function(i) c(sprintf("__%s__", S[[i]]$head), vapply(S[[i]]$body, function(s) fill(s, f, sprintf("S06.steps.%d", i)), "")))
  hc <- max(vapply(seq_len(n), function(i) est_height(vapply(paras[[i]], nobreak, ""), w[i] - 0.04, 16, gap_pt = cgap, card = TRUE), 0)) + 0.08
  for (i in seq_len(n)) {
    x <- GEO$ML + sum(w[seq_len(i - 1)]) + (i - 1) * gap
    deck_text(paras[[i]], c(x, y0, w[i], hc), size = 16, label = sprintf("step%d", i), gap_pt = cgap,
              bg = if (i == n) PAL$tint_orange else PAL$tint_blue, geom = "roundRect")
  }

  # ---- NCA 엔진 검증 표(이 엔진·NonCompart·PKNCA 쌍별 비교) ----
  ds <- c(theoph = "Theoph (12 profiles)", indometh = "Indometh (6 profiles, extravascular rules)", dupi = "simulated dupilumab (1,000 profiles: 500 per model)")
  ev <- rows(EV); premise(all(ev$pass) && setequal(unique(ev$dataset), ds), "every NCA engine comparison passed; three data sets")
  premise(all(ev$lz_points_identical + ev$lz_points_mismatch == ev$n_profiles), "lambda-z window counts are per profile and comparison (table note: matched windows / (profiles x pairs))")
  premise(setequal(unique(ev$comparison), c("this engine vs NonCompart", "PKNCA vs NonCompart", "this engine vs PKNCA")), "three pairwise comparisons: this engine, NonCompart, PKNCA (table header)")
  efmt <- function(x) formatC(x, format = "e", digits = 1)   # 같은 형식(가수 한 자리 소수)
  one <- function(k) {
    w_ <- sprintf("dataset=='%s'", ds[[k]]); r <- rows(EV, w_)
    np <- dint(EV, sprintf("%s & comparison=='this engine vs NonCompart'", w_), "n_profiles", sprintf("NCA engine validation: profiles, %s", k))
    a <- sum(r$lz_points_identical); b <- sum(r$lz_points_identical + r$lz_points_mismatch)
    lz <- dderived(sprintf("lambda-z windows identical, %s, three comparisons", k), EV, sprintf("%s :: sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", w_), c(a, b), sprintf("%s / %s", fint(a), fint(b)))
    mx <- dderived(sprintf("largest relative parameter difference, %s", k), EV, sprintf("%s :: max(max_rel_diff)", w_), max(r$max_rel_diff), efmt(max(r$max_rel_diff)))
    c(DK$txt$S06$table$rows[[k]], np, lz, mx, fill(DK$txt$S06$table$pass, list(k = dcount(EV, sprintf("%s & pass==TRUE", w_), sprintf("comparisons passed, %s", k)))))
  }
  m <- do.call(rbind, lapply(names(ds), one)); df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- unlist(DK$txt$S06$table$head)
  yb <- y0 + hc + 0.2; hb <- GEO$BODY_BOTTOM - yb; tw <- 7.35
  th <- 1.34   # 머리글 한 줄 + 3행(LibreOffice는 행 높이를 가장 높은 행에 맞추므로 머리글을 한 줄로 둔다)
  deck_table(df, box = c(GEO$ML, yb, tw, th), widths = c(2.6, 1.05, 1.3, 1.2, 1.2), size = 13)
  deck_text(unlist(DK$txt$S06$table$caption), c(GEO$ML, yb + th + 0.04, tw, 0.72), size = 16, color = PAL$ink2, label = "table_note", gap_pt = 2)

  # ---- 사전 등록 상자 ----
  c779 <- s06_commit("^Operating-characteristic design", "pre-registration commit, operating-characteristic design")
  c521 <- s06_commit("^Analysis-model re-judgement", "pre-registration commit, analysis models (v1.0.1)")
  c68b <- s06_commit("^Extension to 20,000 trials", "commit of the post hoc extension rule (v1.0)")
  pr <- rows("regulatory/tables/prespecification_register.csv")
  premise(grepl("^post hoc", pr[grepl("^Extension to 20,000 trials", Item), Status]) && grepl("^pre-registered", pr[grepl("^Analysis-model re-judgement", Item), Status]),
          "register: extension rule post hoc in v1.0, analysis models and per-model extension rule pre-registered in v1.0.1 (prereg box)")
  premise(startsWith(rows("oc/prereg.csv")$prereg_commit, c779) && !isTRUE(as.logical(rows("oc/prereg.csv")$changed_since)), "oc/prereg.csv: same commit, design unchanged since")
  xr <- GEO$ML + tw + 0.3
  deck_text(tx("S06.prereg", list(c1 = c779, c2 = c521)), c(xr, yb, GEO$W - GEO$MR - xr, hb), size = 16, label = "prereg", bg = PAL$tint_grey, geom = "roundRect")

  # ---- 노트 ----
  PV <- "params_variability.yaml"
  EXW <- "pk_model=='k2020' & scenario=='V2_up_080' & analysis_model=='M0'"
  premise(!any(rows(T1, "analysis_model=='M1' & config=='P2'")$mechanism == "Km"), "no Km boundary cell (notes)")
  premise(rows(EX, "analysis_model=='M0' & selected==TRUE")$scenario == "V2_up_080" && rows(EX, "analysis_model=='M0' & selected==TRUE")$pk_model == "k2020", "the extended cell is the 2020 model V2 up cell")
  ex0 <- row1(EX, EXW); ex1 <- row1(EX, "pk_model=='k2020' & scenario=='V2_up_080' & analysis_model=='M1'")
  premise(isTRUE(ex0$triggers) && !isTRUE(ex1$triggers) && isTRUE(ex1$selected) && ex1$lo > 5, "extension triggered by M0 only; M1 not triggered (lower bound above 5%) but extended with the cell (card 5 and notes)")
  v10 <- row1("oc/extension_decision.csv", "model=='k2020' & scenario=='V2_up_080' & config=='P2'")
  premise(isTRUE(v10$selected) && abs(v10$pass_pct - ex0$pass_pct) < 1e-9, "the M0 value of the extended cell equals the v1.0 value (known before the v1.0.1 registration; notes)")
  deck_notes(tx("S06.notes", list(
    km = dcfg("params_typical.yaml", c("theta", "Km", "value"), "Km, both models (mg/L)", num_fmt(2)),
    s16 = dcfg(PV, c("residual", "sigma_prop", "value"), "proportional residual, 2016 model (%)", function(x) paste0(fnum(100 * x, 1), "%")),
    s20 = dcfg("params_k2020_model1.yaml", c("residual", "sigma_prop", "value"), "proportional residual, 2020 model (%)", function(x) paste0(fnum(100 * x, 1), "%")),
    wm = dcfg("trial_design.yaml", c("weight", "base", "mean"), "body weight mean (kg)", num_fmt(0)),
    wsd = dcfg("trial_design.yaml", c("weight", "base", "sd"), "body weight SD (kg)", num_fmt(0)),
    wt = f$wt, sex = dcfg("trial_design.yaml", c("weight", "sex_ratio_male", "value"), "share male", num_fmt(1)),
    nind = dcfg("trial_design.yaml", c("mc", "n_individual"), "individual-level subjects per model", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    ntruth = dcfg("oc_design.yaml", c("estimand", "population", "n_subjects"), "common virtual subjects for true values", function(x) fnum(as.numeric(x), 0, big = TRUE)),
    tie = dcfg("nca_rules.yaml", c("standard", "lambda_z", "tie_tolerance"), "lambda-z tie tolerance (adjusted R2)", function(x) format(x, scientific = FALSE)),
    ncmp = dcount(EV, "pass==TRUE", "NCA engine comparisons passed"), ncmp_all = dcount(EV, "TRUE", "NCA engine comparisons"),
    lz = local({ ev <- rows(EV); a <- sum(ev$lz_points_identical); b <- sum(ev$lz_points_identical + ev$lz_points_mismatch)
      dderived("lambda-z windows identical over all comparisons", EV, "sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", c(a, b), sprintf("%s / %s", fint(a), fint(b))) }),
    mx = local({ x <- max(rows(EV)$max_rel_diff); dderived("largest relative parameter difference", EV, "max(max_rel_diff)", x, efmt(x)) }),
    n_arm = f$n_arm, n_rand = f$n_rand, split = f$split,
    df0 = local({ n_ <- as.numeric(.read("config/trial_design.yaml")$n_per_arm); dderived("degrees of freedom of the pooled t-test (2 x n_per_arm - 2)", "config/trial_design.yaml", "n_per_arm :: 2 x n - 2", 2 * n_ - 2, fnum(2 * n_ - 2, 0)) }),
    npm = dderived("boundary cells per PK model", T1, "analysis_model=='M1' & config=='P2' :: rows per pk_model", as.integer(table(rows(T1, "analysis_model=='M1' & config=='P2'")$pk_model)),
                   { tb <- unique(as.integer(table(rows(T1, "analysis_model=='M1' & config=='P2'")$pk_model))); premise(length(tb) == 1, "same number of boundary cells per PK model"); as.character(tb) }),
    ncell = f$ncell, reps = f$reps, ext = f$ext,
    thr = dcfg("prereg_20260926.yaml", c("section1", "extension", "threshold_pct"), "extension threshold (%)", num_fmt(0)),
    m0 = dci(EX, EXW, "pass_pct", "lo", "hi", 2, "%", "M0 P2 boundary type I error at 10,000 trials, extended cell"),
    m1 = dci(EX, "pk_model=='k2020' & scenario=='V2_up_080' & analysis_model=='M1'", "pass_pct", "lo", "hi", 2, "%", "M1 P2 boundary type I error at 10,000 trials, extended cell"),
    c1 = c779, c2 = c521, c3 = c68b)))
  deck_end()
}
