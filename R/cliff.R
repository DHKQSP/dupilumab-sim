# cliff.R — 절벽 채혈 분석 계산(scripts/34_cliff_analysis.R, scripts/46_lloq_sensitivity.R; config/oc_design.yaml cliff, D-041, D-054).
# 절벽 = tmax 이후 순간 반감기 ln2/(−d ln C/dt)가 처음 d일(defs) 미만이 되는 시점부터 참 농도가 LLOQ에 닿는 시점까지.
# d ln C/dt는 ODE 상태값으로 정확히 계산: dC/dt = (ka·흡수구획 − (ke + k12)·central + k21·periph − Vmax·C/(Km + C)·Vc)/Vc. 교차 시점은 격자 사이 보간.
# LLOQ는 벡터로 받는다: 같은 대상자·같은 참 농도 곡선에서 LLOQ별 도달 시각(쌍대). LLOQ가 하나면 scripts/34의 이전 인라인 계산과 같은 산술이다.

cliff_weight_spec <- function(design, w) if (identical(w, "base")) weight_spec_from_design(design, "base") else list(dist = w$dist, trunc = as.numeric(unlist(w$trunc)))

# 한 모델·체중군 N명(시드 derive_seed(seed, model, wname, "subj")). 반환: 대상자 × LLOQ 행(id, lloq, tmax, Cmax, t_lloq, start1, start2, c_start1, c_start2,
#   len1, len2, model, weight, WT). 격자 seq(step, 100, by = step), 1,000명 묶음 풀이(묶음·격자를 바꾸면 수치가 미세하게 달라진다).
cliff_subjects <- function(model, wname, N, seed, step, lloqs, defs, weight_def, design) {
  stopifnot(length(defs) == 2, length(lloqs) >= 1, all(is.finite(lloqs) & lloqs > 0))
  p <- load_params(model); mod <- get_model(p$model_id)
  subj <- with_seed(derive_seed(seed, model, wname, "subj"), make_subjects(N, p, cliff_weight_spec(design, weight_def), 0.5, 0))
  ip <- individual_params(p, subj)
  times <- seq(step, 100, by = step)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
  Ls <- lloqs
  res <- list(); chunk <- 1000L
  for (b in seq_len(ceiling(N / chunk))) {
    rows <- ((b - 1) * chunk + 1):min(b * chunk, N)
    prm <- as.data.frame(ip[rows, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    ev <- rxode2::add.sampling(rxode2::et(amt = design$dose_mg, cmt = "depot", time = 0), times)
    s <- rxode2::rxSolve(mod, params = prm, events = ev, atol = 1e-10, rtol = 1e-8, maxsteps = 2000000L, cores = 1L, returnType = "data.table", addDosing = FALSE)
    s[, id := prm$id[sim.id]]
    s <- merge(s, ip[rows, .(id, Vc_i, ke_i = ke, k12_i = k12, k21_i = k21, Vmax_i = Vmax, Km_i = Km, ka_i = ka)], by = "id")
    absorb <- if ("absc" %in% names(s)) s$absc else s$depot
    s[, dCdt := (ka_i * absorb - (ke_i + k12_i) * central + k21_i * periph - Vmax_i * C / (Km_i + C) * Vc_i) / Vc_i]
    s[, r := -dCdt / C]                                                   # −d ln C/dt
    setorder(s, id, time)
    res[[b]] <- s[, {
      i <- which.max(C); tt <- time; cc <- C; rr <- r; n <- .N
      post <- (i + 1):n
      tl <- vapply(Ls, function(L) {
        kl <- post[which(cc[post] < L)[1]]
        if (is.na(kl)) NA_real_ else tt[kl - 1] + (log(L) - log(cc[kl - 1])) / (log(cc[kl]) - log(cc[kl - 1])) * (tt[kl] - tt[kl - 1]) }, numeric(1))
      st <- vapply(defs, function(dd) {
        thr <- log(2) / dd; k <- post[which(rr[post] >= thr)[1]]
        if (is.na(k)) NA_real_ else if (k == i + 1) tt[k] else tt[k - 1] + (thr - rr[k - 1]) / (rr[k] - rr[k - 1]) * (tt[k] - tt[k - 1]) }, numeric(1))
      cst <- vapply(st, function(x) if (is.na(x)) NA_real_ else exp(approx(tt, log(cc), x)$y), numeric(1))
      .(lloq = Ls, tmax = tt[i], Cmax = cc[i], t_lloq = tl, start1 = st[1], start2 = st[2], c_start1 = cst[1], c_start2 = cst[2])
    }, by = id]
  }
  out <- rbindlist(res)
  out[, `:=`(len1 = pmax(t_lloq - start1, 0), len2 = pmax(t_lloq - start2, 0))]
  out[, `:=`(model = model, weight = wname, WT = subj$WT[match(id, subj$id)])][]
}

# 일정별 절벽 구간 채혈점 수(1, 2, ≥3점)와 Wilson 구간. cs: cliff_subjects의 한 LLOQ 부분(대상자당 1행). sched: 이름 있는 일정 목록(경과일),
# sched_label·min_interval: 일정별 표기·최소 간격. 채혈 시각 시드 derive_seed(seed, model, wname, 일정, "jitter")는 LLOQ와 무관(쌍대).
cliff_count_points <- function(cs, timing, model, wname, sched, sched_label, min_interval, defs, seed, design) {
  if (anyDuplicated(cs$id)) stop("cliff_count_points: 대상자당 1행(한 LLOQ)이어야 합니다")
  rbindlist(lapply(names(sched), function(sn) {
    days <- sched[[sn]]
    ob <- with_seed(derive_seed(seed, model, wname, sn, "jitter"), make_obs_times(cs$id, days, design, jitter = timing == "windowed"))[planned > 0]
    m <- merge(ob, cs[, .(id, start1, start2, t_lloq)], by = "id")
    rbindlist(lapply(seq_along(defs), function(k) {
      stc <- if (k == 1) "start1" else "start2"
      cnt <- m[, .(n = sum(!is.na(t_lloq) & time >= get(stc) & time < t_lloq), st = get(stc)[1], tl = t_lloq[1]), by = id]
      w1 <- wilson_ci(sum(cnt$n >= 1), nrow(cnt)); w2 <- wilson_ci(sum(cnt$n >= 2), nrow(cnt)); w3 <- wilson_ci(sum(cnt$n >= 3), nrow(cnt))
      inw <- cnt[!is.na(tl) & st >= 28 & tl <= 56]                     # 절벽이 Day 29–57(경과일 28–56) 안에 있는 대상자(매일 채혈 창과 비교용)
      data.table(model = model, weight = wname, timing = timing, definition_day = defs[k], schedule = sn, schedule_label = sched_label[[sn]], n_subjects = nrow(cnt),
                 min_interval_day = min_interval[[sn]], pct_ge1 = w1$est, ge1_lo = w1$lo, ge1_hi = w1$hi, pct_ge2 = w2$est, ge2_lo = w2$lo, ge2_hi = w2$hi,
                 pct_ge3 = w3$est, ge3_lo = w3$lo, ge3_hi = w3$hi, n_cliff_in_day29_57 = nrow(inw),
                 pct_ge1_cliff_in_day29_57 = 100 * mean(inw$n >= 1), pct_ge2_cliff_in_day29_57 = 100 * mean(inw$n >= 2), pct_ge3_cliff_in_day29_57 = 100 * mean(inw$n >= 3))
    }))
  }))
}
