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
  deck_kicker(tx("S06.kicker")); deck_title(tx("S06.title", f))

  # ---- 단계 카드 5개(①~⑤): 제목 + 설명. 카드 사이 간격을 두어 단계를 구분한다 ----
  S <- DK$txt$S06$steps; n <- length(S); gap <- 0.14; w <- (GEO$CW - (n - 1) * gap) / n; y0 <- GEO$BODY_TOP + 0.05; hc <- 3.1
  for (i in seq_len(n)) {
    x <- GEO$ML + (i - 1) * (w + gap)
    body <- vapply(S[[i]]$body, function(s) fill(s, f, sprintf("S06.steps.%d", i)), "")
    deck_text(c(sprintf("__%s__", S[[i]]$head), body), c(x, y0, w, hc), size = 16, label = sprintf("step%d", i), gap_pt = 7,
              bg = if (i == n) PAL$tint_orange else PAL$tint_blue, geom = "roundRect")
  }

  # ---- NCA 엔진 검증 표 ----
  ds <- c(theoph = "Theoph (12 profiles)", indometh = "Indometh (6 profiles, extravascular rules)", dupi = "simulated dupilumab (1,000 profiles: 500 per model)")
  ev <- rows(EV); premise(all(ev$pass) && setequal(unique(ev$dataset), ds), "every NCA engine comparison passed; three data sets")
  one <- function(k) {
    w_ <- sprintf("dataset=='%s'", ds[[k]]); r <- rows(EV, w_)
    np <- dint(EV, sprintf("%s & comparison=='this engine vs NonCompart'", w_), "n_profiles", sprintf("NCA engine validation: profiles, %s", k))
    a <- sum(r$lz_points_identical); b <- sum(r$lz_points_identical + r$lz_points_mismatch)
    lz <- dderived(sprintf("lambda-z windows identical, %s, three comparisons", k), EV, sprintf("%s :: sum(lz_points_identical) / sum(lz_points_identical + lz_points_mismatch)", w_), c(a, b), sprintf("%s / %s", fint(a), fint(b)))
    mx <- dderived(sprintf("largest relative parameter difference, %s", k), EV, sprintf("%s :: max(max_rel_diff)", w_), max(r$max_rel_diff), format(signif(max(r$max_rel_diff), 2)))
    c(DK$txt$S06$table$rows[[k]], np, lz, mx, fill(DK$txt$S06$table$pass, list(k = dcount(EV, sprintf("%s & pass==TRUE", w_), sprintf("comparisons passed, %s", k)))))
  }
  m <- do.call(rbind, lapply(names(ds), one)); df <- as.data.frame(m, stringsAsFactors = FALSE); names(df) <- unlist(DK$txt$S06$table$head)
  yb <- y0 + hc + 0.3; hb <- 1.3
  deck_table(df, box = c(GEO$ML, yb, 7.45, hb), widths = c(2.45, 0.95, 1.45, 1.3, 1.3), size = 14)

  # ---- 사전 등록 상자 ----
  c779 <- s06_commit("^Operating-characteristic design", "pre-registration commit, operating-characteristic design")
  c521 <- s06_commit("^Analysis-model re-judgement", "pre-registration commit, analysis models (v1.0.1)")
  premise(startsWith(rows("oc/prereg.csv")$prereg_commit, c779) && !isTRUE(as.logical(rows("oc/prereg.csv")$changed_since)), "oc/prereg.csv: same commit, design unchanged since")
  xr <- GEO$ML + 7.45 + 0.3
  deck_text(tx("S06.prereg", list(c1 = c779, c2 = c521)), c(xr, yb, GEO$W - GEO$MR - xr, hb), size = 16, label = "prereg", bg = PAL$tint_grey, geom = "roundRect")

  # ---- 노트 ----
  PV <- "params_variability.yaml"
  EXW <- "pk_model=='k2020' & scenario=='V2_up_080' & analysis_model=='M0'"
  premise(!any(rows(T1, "analysis_model=='M1' & config=='P2'")$mechanism == "Km"), "no Km boundary cell (notes)")
  premise(rows(EX, "analysis_model=='M0' & selected==TRUE")$scenario == "V2_up_080" && rows(EX, "analysis_model=='M0' & selected==TRUE")$pk_model == "k2020", "the extended cell is the 2020 model V2 up cell")
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
    mx = local({ x <- max(rows(EV)$max_rel_diff); dderived("largest relative parameter difference", EV, "max(max_rel_diff)", x, format(signif(x, 2))) }),
    n_arm = f$n_arm, n_rand = f$n_rand, split = f$split,
    df0 = local({ n_ <- as.numeric(.read("config/trial_design.yaml")$n_per_arm); dderived("degrees of freedom of the pooled t-test (2 x n_per_arm - 2)", "config/trial_design.yaml", "n_per_arm :: 2 x n - 2", 2 * n_ - 2, fnum(2 * n_ - 2, 0)) }),
    npm = dderived("boundary cells per PK model", T1, "analysis_model=='M1' & config=='P2' :: rows per pk_model", as.integer(table(rows(T1, "analysis_model=='M1' & config=='P2'")$pk_model)),
                   { tb <- unique(as.integer(table(rows(T1, "analysis_model=='M1' & config=='P2'")$pk_model))); premise(length(tb) == 1, "same number of boundary cells per PK model"); as.character(tb) }),
    ncell = f$ncell, reps = f$reps, ext = f$ext,
    thr = dcfg("prereg_20260926.yaml", c("section1", "extension", "threshold_pct"), "extension threshold (%)", num_fmt(0)),
    m0 = dci(EX, EXW, "pass_pct", "lo", "hi", 2, "%", "M0 P2 boundary type I error at 10,000 trials, extended cell"),
    c1 = c779, c2 = c521)))
  deck_end()
}
