# 단계 1: Clot 2021 설계 재현 (SPEC §4.2–4.3, DECISIONS D-003)
skip_if_no_rxode2()
p  <- load_params("dev")
dz <- read_cfg("design_clot2021.yaml")
days <- as.numeric(dz$sample_days_post_dose)
lloq <- dz$lloq_mg_L
targets <- rbindlist(lapply(dz$targets, function(x) data.table(dose = x$dose_mg, ratio = x$auc_ratio_pct, tlast = x$tlast_median_day)))

typical_profile <- function(dose, ka_val = NULL) {
  ip <- typical_subject(p); if (!is.null(ka_val)) ip[, ka := ka_val]   # data.table 스코프: 인자명을 열이름과 다르게
  tr <- true_auc(ip, dose, t_grid = days)
  prof <- tr$profile[time > 0]
  tl <- true_auc_to_tlast(prof, lloq)
  list(prof = prof, tlast = tl$tlast_true, AUClast = tl$AUClast_true, AUCinf = tr$inf$AUCinf_true,
       ratio = tl$AUClast_true / tr$inf$AUCinf_true)
}

test_that("(a) 대표 개체의 마지막 정량 시점이 각 용량의 목표 중앙값과 일치 (현재 ka 기준)", {
  for (k in seq_len(nrow(targets))) {
    r <- typical_profile(targets$dose[k])
    expect_equal(r$tlast, targets$tlast[k], info = sprintf("dose=%d mg, ka=%.3g: tlast=%s (목표 %s)", targets$dose[k], p$theta[["ka"]], r$tlast, targets$tlast[k]))
  }
})

test_that("(a') 진단: ka 격자별 tlast 표를 기록한다(ka 보정 전 민감도, 단정 없음 — DECISIONS D-008)", {
  grid <- as.numeric(read_cfg("calibration_targets.yaml")$ka_grid_1_day)
  tab <- rbindlist(lapply(seq_len(nrow(targets)), function(k) rbindlist(lapply(grid, function(ka) {
    r <- typical_profile(targets$dose[k], ka_val = ka)
    data.table(dose = targets$dose[k], ka = ka, tlast = r$tlast, target = targets$tlast[k], match = isTRUE(all.equal(r$tlast, targets$tlast[k])))
  }))))
  message("\nka 격자별 대표 개체 tlast:\n", paste(capture.output(print(dcast(tab, ka ~ dose, value.var = "tlast"))), collapse = "\n"))
  expect_true(nrow(tab) == length(grid) * nrow(targets))
})

test_that("(b) 모델 기반 참값 AUC(0-tlast)/AUC(0-inf) 진단 — 정상 범위 확인과 기록만 (NCA 목표와 직접 비교 불가, DECISIONS D-007)", {
  tab <- rbindlist(lapply(seq_len(nrow(targets)), function(k) {
    r <- typical_profile(targets$dose[k])
    data.table(dose = targets$dose[k], ratio_true_pct = 100 * r$ratio, nca_target_pct = targets$ratio[k], tlast = r$tlast)
  }))
  message("\n모델 기반 참값 비율(대표 개체; NCA 기반 목표는 참고용):\n", paste(capture.output(print(tab)), collapse = "\n"))
  expect_true(all(tab$ratio_true_pct > 90 & tab$ratio_true_pct <= 100))
})

test_that("(c) NCA 기반 AUClast/AUCinf 비 — BEmaster 필요", {
  skip_if_no_bemaster()
  skip_if_variability_pending(p)
  skip_if_ka_pending(p)
  design <- read_cfg("trial_design.yaml")
  res <- lapply(1:200, function(r) {
    s <- with_seed(derive_seed(1, "clot", r), make_subjects(dz$n_per_arm, p, design))
    obs <- make_obs_times(s$ipar$id, days, design, jitter = FALSE)
    sim <- simulate_observations(s$ipar, obs, 300, p$sigma, lloq)
    nca <- run_nca_bemaster(sim, 300, lloq)
    mean(nca$AUClast / nca$AUCinf)
  })
  expect_equal(100 * median(unlist(res)), targets[dose == 300, ratio], tolerance = dz$tolerances$auc_ratio_abs_pct / 100)
})

test_that("(d) IIV 반영 tlast 분포(n=8 반복): 중앙값·범위 — FDA 표 필요", {
  skip_if_variability_pending(p)
  skip_if_ka_pending(p)
  design <- read_cfg("trial_design.yaml")
  med <- range_lo <- range_hi <- numeric(0)
  for (r in 1:200) {
    s <- with_seed(derive_seed(2, "clot", r), make_subjects(dz$n_per_arm, p, design))
    obs <- make_obs_times(s$ipar$id, days, design, jitter = FALSE)
    sim <- with_seed(derive_seed(3, "clot", r), simulate_observations(s$ipar, obs, 300, p$sigma, lloq))
    tl <- observed_tlast(sim)$tlast_planned
    med <- c(med, median(tl)); range_lo <- c(range_lo, min(tl)); range_hi <- c(range_hi, max(tl))
  }
  expect_equal(median(med), targets[dose == 300, tlast])
  rg <- dz$targets[[1]]$tlast_range_day
  expect_lte(median(range_lo), rg[1]); expect_gte(median(range_hi), rg[2])
})

test_that("변동성: 300 mg AUClast CV가 Li 2020 범위(35–51%) 안 — FDA 표·BEmaster 필요", {
  skip_if_no_bemaster(); skip_if_variability_pending(p); skip_if_ka_pending(p)
  succeed("구현 예정: BEmaster NCA로 AUClast CV 산출 후 범위 비교")
})

test_that("변동성: 200 mg AUClast log SD ≈ 0.49 (Cohen 2022) — FDA 표·BEmaster 필요", {
  skip_if_no_bemaster(); skip_if_variability_pending(p); skip_if_ka_pending(p)
  succeed("구현 예정: BEmaster NCA로 log(AUClast) SD 산출 후 비교")
})
