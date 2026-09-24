#!/usr/bin/env Rscript
# §2 절벽 채혈 무의미성 분석 (config/oc_design.yaml cliff, D-041).
# 절벽 = tmax 이후 순간 반감기 ln2/(−d ln C/dt)가 처음 1일(민감도 2일) 미만이 되는 시점부터 참 농도가 LLOQ(0.078 mg/L)에 닿는 시점까지.
# d ln C/dt는 ODE 상태값으로 정확히 계산: dC/dt = (ka·흡수구획 − (ke + k12)·central + k21·periph − Vmax·C/(Km + C)·Vc)/Vc. 교차 시점은 0.05일 격자 사이 선형 보간.
# 조건: 두 모델 × 체중(60–90 kg 기본, 90–110, 110–130, 130–150 kg 균등) × 20,000명 × 채혈 시각(명목, 방문 허용창).
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
oc <- read_cfg("oc_design.yaml"); cf <- oc$cliff; design <- read_cfg("trial_design.yaml")
out_dir <- proj_path("results", if (length(commandArgs(trailingOnly = TRUE))) "cliff_test" else "cliff"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args)) as.integer(args[1]) else as.integer(cf$n_subjects); LLOQ <- as.numeric(cf$lloq); STEP <- as.numeric(cf$grid_step_day); SEED <- as.integer(cf$seed)
defs <- as.numeric(unlist(cf$definition_days))
logfile <- start_run_log("cliff_analysis", master_seed = SEED, run_mode = "final", extra = list(n = N, step = STEP))
b0 <- get_schedule(design, "B0")
SCHED <- list(current = b0, plus_39_46 = sort(c(b0, 38, 45)), plus_39_46_53 = sort(c(b0, 38, 45, 52)),
              plus_32_39_46_53 = sort(c(b0, 31, 38, 45, 52)), plus_40_47 = sort(c(b0, 39, 46)), daily_29_57 = sort(unique(c(b0, 28:56))))
SCHED_LABEL <- c(current = "현행", plus_39_46 = "+Day 39·46", plus_39_46_53 = "+Day 39·46·53", plus_32_39_46_53 = "+Day 32·39·46·53",
                 plus_40_47 = "+Day 40·47", daily_29_57 = "Day 29–57 매일")
min_interval <- vapply(SCHED, function(s) { x <- s[s >= 21]; min(diff(x)) }, numeric(1))    # 경과일 21(Day 22) 이후 최소 방문 간격

wt_spec_for <- function(w) if (identical(w, "base")) weight_spec_from_design(design, "base") else list(dist = w$dist, trunc = as.numeric(unlist(w$trunc)))

# 한 모델·체중군: 20,000명 → 대상자별 tmax, 절벽 시작(정의별), LLOQ 도달 시각, 절벽 시작 농도
cliff_subjects <- function(model, wname) {
  p <- load_params(model); mod <- get_model(p$model_id)
  subj <- with_seed(derive_seed(SEED, model, wname, "subj"), make_subjects(N, p, wt_spec_for(cf$weights[[wname]]), 0.5, 0))
  ip <- individual_params(p, subj)
  times <- seq(STEP, 100, by = STEP)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr")
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
      kl <- post[which(cc[post] < LLOQ)[1]]
      t_lloq <- if (is.na(kl)) NA_real_ else tt[kl - 1] + (log(LLOQ) - log(cc[kl - 1])) / (log(cc[kl]) - log(cc[kl - 1])) * (tt[kl] - tt[kl - 1])
      st <- vapply(defs, function(dd) {
        thr <- log(2) / dd; k <- post[which(rr[post] >= thr)[1]]
        if (is.na(k)) NA_real_ else if (k == i + 1) tt[k] else tt[k - 1] + (thr - rr[k - 1]) / (rr[k] - rr[k - 1]) * (tt[k] - tt[k - 1]) }, numeric(1))
      cst <- vapply(st, function(x) if (is.na(x)) NA_real_ else exp(approx(tt, log(cc), x)$y), numeric(1))
      .(tmax = tt[i], Cmax = cc[i], t_lloq = t_lloq, start1 = st[1], start2 = st[2], c_start1 = cst[1], c_start2 = cst[2])
    }, by = id]
  }
  out <- rbindlist(res)
  out[, `:=`(len1 = pmax(t_lloq - start1, 0), len2 = pmax(t_lloq - start2, 0))]
  out[, `:=`(model = model, weight = wname, WT = subj$WT[match(id, subj$id)])][]
}

count_points <- function(cs, timing, model, wname) {
  rbindlist(lapply(names(SCHED), function(sn) {
    days <- SCHED[[sn]]
    ob <- with_seed(derive_seed(SEED, model, wname, sn, "jitter"), make_obs_times(cs$id, days, design, jitter = timing == "windowed"))[planned > 0]
    m <- merge(ob, cs[, .(id, start1, start2, t_lloq)], by = "id")
    rbindlist(lapply(seq_along(defs), function(k) {
      stc <- if (k == 1) "start1" else "start2"
      cnt <- m[, .(n = sum(!is.na(t_lloq) & time >= get(stc) & time < t_lloq), st = get(stc)[1], tl = t_lloq[1]), by = id]
      w1 <- wilson_ci(sum(cnt$n >= 1), nrow(cnt)); w2 <- wilson_ci(sum(cnt$n >= 2), nrow(cnt)); w3 <- wilson_ci(sum(cnt$n >= 3), nrow(cnt))
      inw <- cnt[!is.na(tl) & st >= 28 & tl <= 56]                     # 절벽이 Day 29–57(경과일 28–56) 안에 있는 대상자(매일 채혈 창과 비교용)
      data.table(model = model, weight = wname, timing = timing, definition_day = defs[k], schedule = sn, schedule_label = SCHED_LABEL[[sn]], n_subjects = nrow(cnt),
                 min_interval_day = min_interval[[sn]], pct_ge1 = w1$est, ge1_lo = w1$lo, ge1_hi = w1$hi, pct_ge2 = w2$est, ge2_lo = w2$lo, ge2_hi = w2$hi,
                 pct_ge3 = w3$est, ge3_lo = w3$lo, ge3_hi = w3$hi, n_cliff_in_day29_57 = nrow(inw),
                 pct_ge1_cliff_in_day29_57 = 100 * mean(inw$n >= 1), pct_ge2_cliff_in_day29_57 = 100 * mean(inw$n >= 2), pct_ge3_cliff_in_day29_57 = 100 * mean(inw$n >= 3))
    }))
  }))
}

all_cs <- list(); pts <- list()
for (model in c("k2016", "k2020")) for (wname in names(cf$weights)) {
  t0 <- Sys.time(); cs <- cliff_subjects(model, wname)
  all_cs[[paste(model, wname)]] <- cs
  for (timing in unlist(cf$timing)) pts[[paste(model, wname, timing)]] <- count_points(cs, timing, model, wname)
  cat(sprintf("[%s] %s %s done in %s\n", format(Sys.time(), "%H:%M:%S"), model, wname, format(Sys.time() - t0)))
}
cs_all <- rbindlist(all_cs); pts <- rbindlist(pts)
saveRDS(cs_all, file.path(out_dir, "cliff_subjects.rds"))
qq <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- cs_all[, .(n = .N, lloq_never_pct = 100 * mean(is.na(t_lloq)),
                   lloq_studyday_median = median(t_lloq + 1, na.rm = TRUE), lloq_studyday_p05 = qq(t_lloq + 1, 0.05), lloq_studyday_p95 = qq(t_lloq + 1, 0.95),
                   lloq_after_day58_pct = 100 * mean(t_lloq + 1 > 58, na.rm = TRUE),
                   c_start1_median = median(c_start1, na.rm = TRUE), c_start2_median = median(c_start2, na.rm = TRUE),
                   len1_median = median(len1, na.rm = TRUE), len1_p05 = qq(len1, 0.05), len1_p95 = qq(len1, 0.95),
                   len2_median = median(len2, na.rm = TRUE), len2_p05 = qq(len2, 0.05), len2_p95 = qq(len2, 0.95),
                   tmax_median = median(tmax)), by = .(model, weight)]
fwrite(summ, file.path(out_dir, "cliff_summary.csv")); fwrite(pts, file.path(out_dir, "cliff_points.csv"))
fwrite(data.table(schedule = names(SCHED), schedule_label = unname(SCHED_LABEL[names(SCHED)]), min_interval_day_from_day22 = unname(min_interval),
                  added_study_days = vapply(SCHED, function(s) paste(setdiff(s, b0) + 1, collapse = " "), character(1))), file.path(out_dir, "cliff_schedules.csv"))
print(summ, digits = 3); print(pts[weight == "base" & timing == "nominal"], digits = 3)

source(proj_path("scripts", "34b_cliff_figures.R"))   # 그림 2-1 ~ 2-4 (저장된 결과에서 다시 그릴 수 있게 분리)
append_run_log(logfile, "done")
