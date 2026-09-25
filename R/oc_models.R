# oc_models.R — 분석 모형 재판정(추가 지시 2026-09-26 §1, config/prereg_20260926.yaml, D-053)
#   M0: log 척도 pooled 두 표본 t (현행, be_pooled_t)
#   M1: log 척도 ANOVA — 처리 + 체중 층(2수준, 배정 층과 같은 경계). 주분석 후보
#   M2: log 척도 ANCOVA — 처리 + log 체중(연속). 민감도
# M1·M2는 시험 × 평가변수 묶음마다 lm을 부르지 않고 정규방정식(3 × 3)을 닫힌 꼴로 푼다(y는 묶음 안에서 중심화해 상쇄 오차를 줄임).
# tests/testthat/test-oc-models.R가 lm/confint와 대조한다.

# 대칭 3 × 3 행렬(a11 a12 a13 / a22 a23 / a33)의 역행렬 원소와 해. 벡터화(묶음마다 한 값).
.solve_sym3 <- function(a11, a12, a13, a22, a23, a33, b1, b2, b3) {
  c11 <- a22 * a33 - a23^2; c12 <- -(a12 * a33 - a13 * a23); c13 <- a12 * a23 - a13 * a22
  c22 <- a11 * a33 - a13^2; c23 <- -(a11 * a23 - a12 * a13); c33 <- a11 * a22 - a12^2
  det <- a11 * c11 + a12 * c12 + a13 * c13
  scale <- pmax(abs(a11 * a22 * a33), .Machine$double.xmin)
  bad <- !is.finite(det) | abs(det) <= 1e-10 * scale          # 특이(층이나 군이 비어 있음 등) → NA
  det[bad] <- NA_real_
  i11 <- c11 / det; i12 <- c12 / det; i13 <- c13 / det; i22 <- c22 / det; i23 <- c23 / det; i33 <- c33 / det
  list(beta1 = i11 * b1 + i12 * b2 + i13 * b3, beta2 = i12 * b1 + i22 * b2 + i23 * b3, beta3 = i13 * b1 + i23 * b2 + i33 * b3, inv22 = i22)
}

# d: data.table(<key cols>, arm ("R"/"T"), y (양수, 원척도), stratum (0/1: 두 번째 층이면 1), lwt (log 체중)).
# 반환: key × model(M1, M2)별 est(log GMR), se, df, GMR, CI_lower, CI_upper, pass, n_R, n_T. 유효하지 않은 y(NA, ≤ 0)는 제외.
be_models_fast <- function(d, keys, ci_level = 0.90, limits = c(0.80, 1.25)) {
  if (!all(d$stratum %in% c(0, 1))) stop("be_models_fast: stratum은 0/1 부호여야 합니다(M1의 a33 = sum(stratum)은 0/1에서만 제곱합과 같음)")
  x <- d[is.finite(y) & y > 0 & !is.na(arm)]
  x[, `:=`(t = as.numeric(arm == "T"), ly = log(y))]
  x[, `:=`(lyc = ly - mean(ly), xc = lwt - mean(lwt)), by = keys]          # 묶음 안 중심화(절편이 있으므로 추정량 불변)
  g <- x[, .(n = .N, nT = sum(t), ns = sum(stratum), nts = sum(t * stratum),
             sy = sum(lyc), sty = sum(t * lyc), ssy = sum(stratum * lyc), syy = sum(lyc^2),
             sx = sum(xc), stx = sum(t * xc), sxx = sum(xc^2), sxy = sum(xc * lyc)), by = keys]
  g[, n_R := n - nT]; g[, n_T := nT]
  q_ <- function(df) qt(1 - (1 - ci_level) / 2, df)
  fit <- function(a13, a23, a33, b3, model) {
    s <- .solve_sym3(g$n, g$nT, a13, g$nT, a23, a33, g$sy, g$sty, b3)
    df <- g$n - 3
    rss <- g$syy - (s$beta1 * g$sy + s$beta2 * g$sty + s$beta3 * b3)
    rss[rss < 0 & rss > -1e-9] <- 0
    se <- sqrt(rss / df * s$inv22)
    ok <- g$n_R >= 2 & g$n_T >= 2 & df >= 1 & is.finite(se)
    est <- ifelse(ok, s$beta2, NA_real_); se[!ok] <- NA_real_
    lo <- exp(est - q_(pmax(df, 1)) * se); hi <- exp(est + q_(pmax(df, 1)) * se)
    out <- g[, keys, with = FALSE]
    out[, `:=`(model = model, est = est, se = se, df = df, GMR = exp(est), CI_lower = lo, CI_upper = hi,
               pass = lo >= limits[1] & hi <= limits[2], n_R = g$n_R, n_T = g$n_T)]
    out[]
  }
  rbind(fit(g$ns, g$nts, g$ns, g$ssy, "M1"), fit(g$sx, g$stx, g$sxx, g$sxy, "M2"))
}

# 같은 자료의 M0(pooled t) — be_pooled_t와 같은 식(대조·교차 확인용; 저장본 대조는 be_pooled_t 결과로 한다)
be_m0_fast <- function(d, keys, ci_level = 0.90, limits = c(0.80, 1.25)) {
  x <- d[is.finite(y) & y > 0 & !is.na(arm)]
  x[, ly := log(y)]
  g <- x[, .(n_R = sum(arm == "R"), n_T = sum(arm == "T"), mR = mean(ly[arm == "R"]), mT = mean(ly[arm == "T"]),
             vR = var(ly[arm == "R"]), vT = var(ly[arm == "T"])), by = keys]
  g[, df := n_R + n_T - 2]
  g[, se := sqrt(((n_R - 1) * vR + (n_T - 1) * vT) / df * (1 / n_R + 1 / n_T))]
  g[, est := mT - mR]
  g[!(n_R >= 2 & n_T >= 2), `:=`(est = NA_real_, se = NA_real_)]
  q <- qt(1 - (1 - ci_level) / 2, pmax(g$df, 1))
  g[, `:=`(model = "M0", GMR = exp(est), CI_lower = exp(est - q * se), CI_upper = exp(est + q * se))]
  g[, pass := CI_lower >= limits[1] & CI_upper <= limits[2]]
  g[, c(keys, "model", "est", "se", "df", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T"), with = FALSE]
}

# --- 시험 재생성(개체 수준) ------------------------------------------------------------------------------------------
# 시험 j를 run_trial_oc_ext와 같은 simulate_trial_arms 호출·시드·NCA 입력(대조군 id 1..nmax, 시험군 시나리오 k는 id + k·nmax)으로 다시 만들고,
# 시나리오 × 평가변수마다 분석 모형을 적용한다.
#   M0 = be_pooled_t(run_trial_oc_ext와 같은 입력 순서 c(대조군, 시험군) → 저장본과 대조), M1·M2 = be_models_fast.
#   공변량은 배정 때 기록한 대상자 표(sa$subj)의 체중과 배정 층(1번째 층 = 60–75 kg, 2번째 = >75–90 kg)이다.
# lloqs: NULL이면 모의에 쓴 연구 LLOQ(sa$lloq = study_lloq()). 여러 값이면 같은 y_raw(같은 대상자·채혈 시각·잔차)를 LLOQ마다 simulate_observations와
#   같은 규칙(투여 전 또는 y_raw < LLOQ → BLQ, conc NA)으로 다시 검열한다(LLOQ 간 쌍대 비교).
# resid: "fixed" = 잔차 모형 그대로(가산 SD는 모델 개발 자료의 추정값), "scaled" = 가산 잔차를 LLOQ / p$lloq 배로(LLOQ에서의 상대 정밀도 유지; §2 민감도).
#   scaled는 시험군·대조군 잔차를 simulate_trial_arms와 같은 시드(master, j, arm, "eps")로 다시 뽑아 y_raw = C·(1 + eps_p) + eps_a·배율로 만든다.
#   배율 1에서 다시 만든 y_raw가 저장 y_raw와 비트 단위로 같은지 시험마다 검사한다.
#   (LLOQ, resid) 조합 l의 id 오프셋은 (l − 1)·(시나리오 수 + 1)·nmax.
# 반환 list(be = trial, lloq, resid, scenario, endpoint, model, est, se, df, GMR, CI_lower, CI_upper, pass, n_R, n_T;
#           drop = trial, lloq, resid, scenario(REF = 대조군), arm, n, n_lambda, n_reliable(세트 ii), n_reliable_i(세트 i), flag_rsq, flag_extrap, flag_span, tlast_median;
#           strata = trial, arm, stratum, n)
recensor_obs <- function(obs, lloq) {
  l_ <- lloq                                                           # 인자를 지역 변수로(data.table 스코프, D-022)
  o <- copy(obs)
  o[, blq := planned == 0 | y_raw < l_]
  o[, conc := fifelse(blq, NA_real_, y_raw)][]
}
# 잔차를 같은 시드로 다시 뽑아 가산 잔차 배율을 바꾼 y_raw. eps: draw_eps 결과(id, planned, eps_p, eps_a). 배율 1이면 원래 y_raw와 같아야 한다.
rescale_additive <- function(obs, eps, scale) {
  s_ <- scale
  o <- merge(copy(obs), eps, by = c("id", "planned"), all.x = TRUE, sort = TRUE)
  if (anyNA(o$eps_p)) stop("rescale_additive: 잔차가 없는 관측")
  chk <- o$C * (1 + o$eps_p) + o$eps_a
  if (!identical(chk, o$y_raw)) stop("rescale_additive: 다시 뽑은 잔차가 저장 y_raw를 재현하지 않습니다(시드·격자 확인)")
  o[, y_raw := C * (1 + eps_p) + eps_a * s_]
  o[, c("eps_p", "eps_a") := NULL]
  setorder(o, id, time)[]
}

run_trial_oc_models <- function(j, p, design, scenarios, master_seed, wt_spec, schedule = "B0", jitter = TRUE, model_id = NULL,
                                lloqs = NULL, resid = "fixed", models = c("M0", "M1", "M2"), endpoints = OC_ENDPOINTS_EXT) {
  stopifnot(all(models %in% c("M0", "M1", "M2")), all(endpoints %in% c(OC_ENDPOINTS_EXT, OC_ENDPOINTS_CRIT)), length(resid) >= 1, all(resid %in% c("fixed", "scaled")), !anyDuplicated(resid))
  sa <- simulate_trial_arms(j, p, design, scenarios, schedule, master_seed, wt_spec, jitter = jitter, model_id = model_id)
  sd_days <- get_schedule(design, schedule)
  scs <- names(scenarios); nmax <- max(sa$subj$id); K <- length(scs) + 1L
  if (is.null(lloqs)) lloqs <- sa$lloq
  if (anyDuplicated(lloqs) || any(!is.finite(lloqs) | lloqs <= 0)) stop("run_trial_oc_models: lloqs는 서로 다른 양수")
  if (nlevels(sa$subj$stratum) != 2L) stop("run_trial_oc_models: 배정 층이 2개가 아닙니다")
  cv <- sa$subj[, .(sid = id, WT, s2 = as.integer(as.integer(stratum) == 2L))]
  arms_l <- c(list(REF = sa$arms$R$S00), setNames(lapply(scs, function(sc) sa$arms$T[[sc]]), scs))
  arm_of <- c(REF = "R", setNames(rep("T", length(scs)), scs))
  eps_l <- NULL
  if ("scaled" %in% resid) {                                          # simulate_trial_arms와 같은 시드·격자로 잔차 재추출
    gp <- sort(unique(c(0, sa$grid)))
    eps_l <- lapply(c(R = "R", T = "T"), function(a) with_seed(derive_seed(master_seed, j, a, "eps"), draw_eps(sa$subj[arm == a, id], gp, p$sigma)))
  }
  combos <- CJ(resid = resid, lloq = lloqs, sorted = FALSE)
  obs_l <- vector("list", nrow(combos) * K); tr_l <- obs_l; map_l <- obs_l; i <- 0L
  for (ci in seq_len(nrow(combos))) for (k in seq_len(K)) {
    L <- combos$lloq[ci]; rv_ <- combos$resid[ci]
    x <- arms_l[[k]]; o <- ((ci - 1L) * K + (k - 1L)) * nmax; i <- i + 1L
    ob <- subset_schedule(x$obs, sd_days)
    if (rv_ == "scaled") ob <- rescale_additive(ob, eps_l[[arm_of[[names(arms_l)[k]]]]], L / p$lloq)
    if (rv_ == "scaled" || L != sa$lloq) ob <- recensor_obs(ob, L)  # 모의 LLOQ·고정 잔차면 저장 검열 그대로(규칙이 같아 결과도 같음)
    obs_l[[i]] <- ob[, .(id = id + o, time, conc)]
    tr_l[[i]] <- x$truth[, .(id = id + o, AUCinf_true)]
    map_l[[i]] <- data.table(id = x$subj$id + o, sid = x$subj$id, lloq = L, resid = rv_, scenario = names(arms_l)[k])
  }
  nca <- run_nca(rbindlist(obs_l))
  nca <- rbindlist(tr_l)[nca, on = "id"]
  nca <- rbindlist(map_l)[nca, on = "id"]
  nca[, AUCinf_C := fifelse(reliable %in% TRUE, AUCinf, AUClast)]
  nca[, rel_i := (lambda_ok & !flag_rsq & !flag_extrap) %in% TRUE]
  nca[, AUCinf_Ci := fifelse(rel_i, AUCinf, AUClast)]
  nca[, `:=`(rel_iii = crit_ok(.SD, "iii"), rel_iv = crit_ok(.SD, "iv"))]            # 기준 세트 (iii)·(iv) (R/criteria.R; 기존 평가변수 불변)
  nca <- cv[nca, on = "sid"]
  setorder(nca, id)                                                   # run_trial_oc_ext와 같은 대상자 순서(NCA 출력은 id 순)
  ev <- function(ep) switch(ep,
    Cmax = nca$Cmax, AUClast = nca$AUClast, AUCinf_A = fifelse(nca$reliable %in% TRUE, nca$AUCinf, NA_real_),
    AUCinf_B = fifelse(nca$lambda_ok %in% TRUE, nca$AUCinf, NA_real_), AUCinf_C = nca$AUCinf_C, AUCinf_true = nca$AUCinf_true,
    AUCinf_Ai = fifelse(nca$rel_i, nca$AUCinf, NA_real_), AUCinf_Ci = nca$AUCinf_Ci,
    AUCinf_Aiii = fifelse(nca$rel_iii, nca$AUCinf, NA_real_), AUCinf_Aiv = fifelse(nca$rel_iv, nca$AUCinf, NA_real_),
    AUCinf_Ciii = fifelse(nca$rel_iii, nca$AUCinf, nca$AUClast), AUCinf_Civ = fifelse(nca$rel_iv, nca$AUCinf, nca$AUClast))
  long <- rbindlist(lapply(endpoints, function(ep) data.table(id = nca$id, lloq = nca$lloq, resid = nca$resid, scenario = nca$scenario, endpoint = ep, y = ev(ep),
                                                              stratum = nca$s2, lwt = log(nca$WT))))
  ref <- long[scenario == "REF"]
  d <- rbindlist(lapply(scs, function(sc) rbind(ref[, .(lloq, resid, scenario = sc, endpoint, arm = "R", id, y, stratum, lwt)],
                                                long[scenario == sc, .(lloq, resid, scenario, endpoint, arm = "T", id, y, stratum, lwt)])))
  # 묶음 안 순서: 대조군(id 순) 다음 시험군(id 순) = run_trial_oc_ext의 c(yR, yT)
  keys <- c("lloq", "resid", "scenario", "endpoint")
  lim <- as.numeric(unlist(design$be$limits)); cl <- design$be$ci_level
  out <- list()
  if ("M0" %in% models) {
    out$M0 <- d[, { r <- be_pooled_t(y, arm, cl, lim)
                    .(model = "M0", est = log(r$GMR), se = r$se, df = r$df, GMR = r$GMR, CI_lower = r$CI_lower, CI_upper = r$CI_upper, pass = r$pass, n_R = r$n_R, n_T = r$n_T) },
                by = keys]
  }
  if (any(c("M1", "M2") %in% models)) out$M12 <- be_models_fast(d, keys, cl, lim)[model %in% models]
  be <- rbindlist(out, use.names = TRUE)
  be[, trial := j]
  setcolorder(be, c("trial", keys, "model"))
  setorderv(be, c("resid", "lloq", "scenario", "endpoint", "model"))
  drop <- nca[, .(n = .N, n_lambda = sum(lambda_ok %in% TRUE), n_reliable = sum(reliable %in% TRUE), n_reliable_i = sum(rel_i), n_reliable_iii = sum(rel_iii), n_reliable_iv = sum(rel_iv),
                  flag_rsq = sum(flag_rsq %in% TRUE), flag_extrap = sum(flag_extrap %in% TRUE), flag_span = sum(flag_span %in% TRUE),
                  tlast_median = as.numeric(median(tlast, na.rm = TRUE))), by = .(lloq, resid, scenario)]
  drop[, `:=`(trial = j, arm = fifelse(scenario == "REF", "R", "T"))]
  strata <- sa$subj[, .(n = .N), by = .(arm, stratum = as.character(stratum))][, trial := j]
  list(be = be[], drop = drop[], strata = strata[])
}
