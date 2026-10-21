# ICH M15 분석 함수(R/m15.R, D-068)의 검증: 집단 샘플러, 배정, 허용창 편차, NCA 보조, 중첩 표본 수 분석 모형, 부분 AUC, AUROC, 참값
source(proj_path("R", "m15.R"))

POP_T <- list(weight = list(dist = "range_truncnorm", lo = 58.0, hi = 95.1), weight_sensitivity = list(dist = "uniform", lo = 58.0, hi = 95.1),
              height = list(by_sex = FALSE, mean_cm = 172.5, sd_cm = 10.4))
POP_K <- list(weight = list(dist = "range_truncnorm", lo = 52.1, hi = 78.3), height = list(by_sex = FALSE, mean_cm = 167.0, sd_cm = 8.68))
INCL <- list(male_kg = c(60, 90), female_kg = c(50, 90))

test_that("집단 샘플러: 성별 선정 범위, 범위 기반 절단 정규, 체중-키 상관, 집단 간 공통 난수", {
  s1 <- m15_subjects(40000, POP_T, INCL, 0.513, 0.5, 11)
  s2 <- m15_subjects(40000, POP_K, INCL, 0.513, 0.5, 11)
  expect_true(all(s1[sex == "M", WT >= 60 & WT <= 90])); expect_true(all(s1[sex == "F", WT >= 58 & WT <= 90]))
  expect_true(all(s2[sex == "M", WT >= 60 & WT <= 78.3])); expect_true(all(s2[sex == "F", WT >= 52.1 & WT <= 78.3]))
  expect_equal(mean(s1$sex == "M"), 0.5, tolerance = 0.02)
  expect_equal(s1$z_Vc, s2$z_Vc); expect_equal(s1$sex, s2$sex)                 # η 점수와 성별이 집단과 무관
  # 역누적분포가 절단 정규와 일치(여성, 58-90 범위의 범위 기반 분포)
  sp <- m15_weight_spec(POP_T$weight); f <- s1[sex == "F", WT]
  ks <- suppressWarnings(ks.test(f, function(x) (pnorm((x - sp$mean) / sp$sd) - pnorm((58 - sp$mean) / sp$sd)) / (pnorm((90 - sp$mean) / sp$sd) - pnorm((58 - sp$mean) / sp$sd))))
  expect_gt(ks$p.value, 0.001)
  zw <- qnorm(with_seed(11, list(runif(40000), runif(40000)))[[2]])          # 같은 시드의 두 번째 흐름 = 체중 균등 난수
  zh <- (s1$HT - 172.5) / 10.4
  expect_equal(cor(zw, zh), 0.513, tolerance = 0.02)
  expect_equal(s1$BMI, s1$WT / (s1$HT / 100)^2)
})

test_that("η 표: 모델 IIV로 척도화, IIV 0인 파라미터는 0, 두 모델이 같은 점수", {
  s <- m15_subjects(500, POP_T, INCL, 0.513, 0.5, 3)
  p20 <- load_params("k2020"); p16 <- load_params("k2016")
  e20 <- m15_eta_table(s, p20); e16 <- m15_eta_table(s, p16)
  expect_equal(e20$eta_Vc, s$z_Vc * p20$omega[["Vc"]]); expect_equal(e16$eta_Vc, s$z_Vc * p16$omega[["Vc"]])
  expect_true(all(e20$eta_F == 0)); expect_true(all(e16$eta_k12 == 0))
  ip <- individual_params(p20, e20); expect_equal(nrow(ip), 500)
})

test_that("층 안 크기 2 블록 배정: 층마다 균형(차이 1 이하), 군 안 등록 순서 연속, 앞 n명 중첩", {
  set.seed(5); wt <- runif(330, 55, 92); u <- runif(330)
  a <- m15_allocate(wt, 75, u)
  for (k in 0:1) { x <- a[stratum == k]; expect_lte(abs(sum(x$arm == "R") - sum(x$arm == "T")), 1) }
  for (arm_ in c("R", "T")) expect_equal(sort(a[arm == arm_, arm_order]), seq_len(sum(a$arm == arm_)))
  expect_gte(min(table(a$arm)), 160)
  # 등록 순서상 앞 블록의 배정은 뒤 대상자에 의존하지 않는다(중첩)
  b <- m15_allocate(wt[1:200], 75, u[1:200]); expect_equal(b$arm, a$arm[1:200])
})

test_that("벡터화 허용창 편차는 R/sampling.R jitter_times와 같은 규칙", {
  design <- read_cfg("trial_design.yaml"); pl <- c(0, 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56, 63, 70)
  set.seed(9); U <- matrix(runif(50 * length(pl)), 50)
  ob <- m15_obs_times(1:50, pl, design, U)
  w <- window_for_time(pl, design)
  for (i in c(1, 17, 50)) {
    a <- pl + (2 * U[i, ] - 1) * w; a[pl == 0] <- 0; a <- pmax(a, 0)
    for (k in 2:length(a)) if (a[k] < a[k - 1] + 1 / 24) a[k] <- a[k - 1] + 1 / 24
    expect_equal(ob[id == i, time], a)
  }
  expect_true(all(ob[planned == 0, time] == 0))
})

test_that("m15_nca_basic은 run_nca의 Cmax·AUClast·Tlast와 같다", {
  set.seed(21)
  sim <- rbindlist(lapply(1:60, function(i) { t <- c(0, 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56)
    c_ <- 30 * (exp(-0.06 * t) - exp(-0.5 * t)) * exp(rnorm(length(t), 0, 0.2)) - runif(1, 0, 8) * (t > 30)
    data.table(id = i, time = t, conc = ifelse(t == 0 | c_ < 0.078, NA_real_, c_)) }))
  sim <- rbind(sim, data.table(id = 61L, time = c(0, 1, 3), conc = NA_real_))
  a <- run_nca(sim); b <- m15_nca_basic(sim)
  expect_equal(b$id, a$id); expect_equal(b$Cmax, a$Cmax); expect_equal(b$AUClast[!is.na(a$Cmax)], a$AUClast[!is.na(a$Cmax)]); expect_equal(b$Tlast, a$Tlast)
})

test_that("중첩 표본 수 분석: 군당 n마다 M0는 be_pooled_t, M1은 lm/confint와 같다(빈 층이면 층 항 제외)", {
  set.seed(7); N <- 330; wt <- runif(N, 55, 92); al <- m15_allocate(wt, 75, runif(N))
  d <- data.table(kk = 1L, arm = al$arm, arm_order = al$arm_order, s = al$stratum, ok = runif(N) > 0.1,
                  y = exp(6 + 0.05 * (al$arm == "T") - 0.3 * al$stratum + rnorm(N, 0, 0.4)))
  d <- d[arm_order <= 160]
  r <- m15_be_nested(d, "kk", c(80, 120, 160))
  for (n in c(80, 120, 160)) {
    z <- d[arm_order <= n & ok]
    b0 <- be_pooled_t(z$y, z$arm); r0 <- r[model == "M0" & n_rand == n]
    expect_equal(exp(r0$est), b0$GMR, tolerance = 1e-10); expect_equal(r0$CI_lower, b0$CI_lower, tolerance = 1e-10); expect_equal(r0$pass, b0$pass)
    f1 <- lm(log(y) ~ factor(arm, c("R", "T")) + factor(s), data = z); ci <- exp(confint(f1, level = 0.90)[2, ]); r1 <- r[model == "M1" & n_rand == n]
    expect_equal(exp(r1$est), unname(exp(coef(f1)[2])), tolerance = 1e-10); expect_equal(c(r1$CI_lower, r1$CI_upper), unname(ci), tolerance = 1e-10)
    expect_equal(r1$df, f1$df.residual)
  }
  d2 <- copy(d)[, s := 0L]; r2 <- m15_be_nested(d2, "kk", 160)
  expect_true(r2[model == "M1", m1_dropped]); expect_equal(r2[model == "M1", est], r2[model == "M0", est])
})

test_that("부분 AUC는 같은 정의에서 PKNCA pk.calc.aucint(lin up/log down)와 같다", {
  skip_if_not_installed("PKNCA")
  set.seed(13)
  sim <- rbindlist(lapply(1:200, function(i) { t <- sort(c(0, 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56) + c(0, runif(13, -0.2, 0.2)))
    c_ <- 30 * (exp(-0.05 * t) - exp(-0.5 * t)) * exp(rnorm(length(t), 0, 0.15)); data.table(id = i, time = t, conc = ifelse(t == 0, NA_real_, c_)) }))
  iv <- list(p0_14 = c(0, 14), p14_28 = c(14, 28), p28_40 = c(28, 40))
  a <- m15_pauc(sim, iv)
  pre <- preprocess_blq(sim)
  for (i in c(1, 50, 200)) {
    z <- pre[id == i]
    for (nm in names(iv)) {
      b <- PKNCA::pk.calc.aucint(conc = z$conc, time = z$time, start = iv[[nm]][1], end = iv[[nm]][2], method = "lin up/log down", auc.type = "AUClast")
      expect_equal(a[id == i, get(nm)], b, tolerance = 1e-8)
    }
  }
})

test_that("AUROC: Mann-Whitney(동점 1/2)와 DeLong 분산이 직접 계산과 같다", {
  set.seed(3); pos <- round(rnorm(80, 0.5), 1); neg <- round(rnorm(90, 0), 1)
  r <- m15_auroc(pos, neg)
  psi <- outer(pos, neg, function(x, y) (x > y) + 0.5 * (x == y))
  expect_equal(r$auroc, mean(psi))
  v10 <- rowMeans(psi); v01 <- colMeans(psi)
  expect_equal(r$se, sqrt(var(v10) / 80 + var(v01) / 90))
})


test_that("연속 시간 참 Tlast와 AUC: 0.1일 격자 보간이 0.001일 격자 직접 계산과 맞는다", {
  skip_if_no_rxode2()
  p <- load_params("k2020"); s <- m15_subjects(6, POP_T, INCL, 0.513, 0.5, 8); ip <- individual_params(p, m15_eta_table(s, p))
  a <- m15_true_tlast_auc(ip, 300, p$model_id, 0.078)
  mod <- get_model(p$model_id)
  prm <- as.data.frame(ip[, c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", "ktr", "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
  d <- .solve_param_table(mod, prm, seq(0.001, 100, by = 0.001), 300)
  ref <- d[, { i <- which.max(C); k <- which(seq_len(.N) > i & C < 0.078)[1]; .(t_ref = time[k - 1], a_ref = auc[k - 1]) }, by = id]
  m <- merge(a, ref, by = "id")
  expect_true(all(!m$censored)); expect_lt(max(abs(m$t_true_last - m$t_ref)), 0.01)   # 절벽 구간의 로그 선형 보간 오차(약 7분)
  expect_lt(max(abs(m$AUClast_true_cont / m$a_ref - 1)), 1e-4)
})

test_that("m15_truth 조밀 채혈 AUC0-last 참값: 0.25일 격자 보간이 0.01일 격자 기준과 맞고, min(참 Tlast, 마지막 계획 시각) 정의를 따른다", {
  skip_if_no_rxode2()
  for (m in c("k2020", "k2016")) for (dose in c(200, 300)) {
    p <- load_params(m); s <- m15_subjects(40, POP_T, INCL, 0.513, 0.5, 31); ip <- individual_params(p, m15_eta_table(s, p))
    sch <- list(B0 = c(0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56), A70 = c(0.25, 1, 28, 56, 63, 70))
    tr <- m15_truth(ip, dose, p$model_id, sch, 0.078, list())
    ref <- m15_true_tlast_auc(ip, dose, p$model_id, 0.078, step = 0.01, t_end = 70)
    mod <- get_model(p$model_id)
    prm <- as.data.frame(ip[, c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ip)) "ktr", "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
    a_end <- .solve_param_table(mod, prm, c(56, 70), dose)[, .(id, time, auc)]
    z <- merge(merge(tr, ref, by = "id"), dcast(a_end, id ~ time, value.var = "auc"), by = "id")
    setnames(z, c("56", "70"), c("a56", "a70"))
    z[, `:=`(exp_B0 = fifelse(!censored & t_true_last.y < 56, AUClast_true_cont, a56), exp_A70 = fifelse(!censored & t_true_last.y < 70, AUClast_true_cont, a70))]
    expect_lt(max(abs(z$AUClast_dense_B0 / z$exp_B0 - 1)), 1e-4); expect_lt(max(abs(z$AUClast_dense_A70 / z$exp_A70 - 1)), 1e-4)
    fin <- is.finite(z$t_true_last.x) & !z$censored
    expect_lt(max(abs(z$t_true_last.x[fin] - z$t_true_last.y[fin])), 0.1)
    expect_true(all(z$AUClast_dense_B0 <= z$AUClast_dense_A70 + 1e-9))
  }
})

test_that("농도 배율 g: AUC0-inf·Cmax·부분 AUC 참값은 정확히 g배, 조밀 AUC0-last 참값은 LLOQ 교차 이동만큼(|log 비 - log g| < 0.002)", {
  skip_if_no_rxode2()
  p <- load_params("k2020"); s <- m15_subjects(200, POP_T, INCL, 0.513, 0.5, 4); ip <- individual_params(p, m15_eta_table(s, p))
  mod <- get_model(p$model_id)
  prm <- as.data.frame(ip[, c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", "ktr", "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
  d <- .solve_param_table(mod, prm, c(seq(0.25, 30, 0.25), seq(30.5, 70, 0.5)), 300)[, .(id, time, C, auc)]; setorder(d, id, time)
  dense <- function(x) { tl <- .m15_tlast_grid(x, 0.078); a <- x[time == 56, .(id, a56 = auc)]; merge(tl, a, by = "id")[, fifelse(t_true_last >= 56, a56, auc_true_last)] }
  r <- dense(d)
  for (g in c(0.8, 1.25)) {
    dg <- copy(d)[, `:=`(C = g * C, auc = g * auc)]
    expect_equal(dg[, max(C), by = id]$V1 / d[, max(C), by = id]$V1, rep(g, 200), tolerance = 1e-12)
    expect_equal(dg[time == 14, auc] - dg[time == 0.25, auc], g * (d[time == 14, auc] - d[time == 0.25, auc]), tolerance = 1e-12)
    expect_lt(abs(log(m15_paired_gmr(dense(dg), r)[["gmr"]]) - log(g)), 0.002)
  }
})

test_that("Tango 점수 신뢰구간: 경계에서 수치 제한 최대우도 점수가 ±z, 불일치 0이면 z^2/(n + z^2)", {
  z <- qnorm(0.975)
  num_score <- function(b, c, n, D) {
    nll <- function(q) { p12 <- q + D; if (q <= 0 || p12 <= 0 || 1 - 2 * q - D <= 0) return(1e10); -((n - b - c) * log(1 - 2 * q - D) + b * log(p12) + c * log(q)) }
    lo_q <- max(0, -D) + 1e-12; hi_q <- (1 - D) / 2 - 1e-12
    q <- optimize(nll, c(lo_q, hi_q), tol = 1e-14)$minimum
    (b - c - n * D) / sqrt(n * (2 * q + D * (1 - D)))
  }
  for (bc in list(c(12, 30, 1000), c(5, 1, 200), c(40, 41, 1000))) {
    r <- m15_tango_ci(bc[1], bc[2], bc[3])
    expect_equal(num_score(bc[1], bc[2], bc[3], r[["lo"]]), z, tolerance = 1e-5); expect_equal(num_score(bc[1], bc[2], bc[3], r[["hi"]]), -z, tolerance = 1e-5)
    expect_equal(r[["diff"]], (bc[1] - bc[2]) / bc[3])
  }
  r0 <- m15_tango_ci(0, 0, 1000); expect_equal(r0[["hi"]], z^2 / (1000 + z^2), tolerance = 1e-8); expect_equal(r0[["lo"]], -z^2 / (1000 + z^2), tolerance = 1e-8)
})

test_that("순위 기반 AUROC와 재표집 값: m15_auroc와 같고, 항등 색인이면 점추정", {
  set.seed(6); neg <- round(abs(rnorm(200, 0, 0.05)), 2); pos <- cbind(round(abs(rnorm(200, 0.04, 0.05)), 2), round(abs(rnorm(200, -0.04, 0.05)), 2))
  expect_equal(.m15_mw(as.vector(pos), neg), m15_auroc(as.vector(pos), neg)$auroc)
  expect_equal(m15_auroc_draws(pos, neg, list(seq_len(200))), m15_auroc(as.vector(pos), neg)$auroc)
  idx <- m15_boot_index(200, 50, 3); expect_length(idx, 50); expect_equal(idx, m15_boot_index(200, 50, 3))
})

test_that("λz 외삽 부분 AUC는 지수 꼬리 적분과 같고, Tlast가 구간 끝 뒤이거나 λz가 없으면 그대로", {
  lz <- 0.12; cl <- 0.5
  f <- function(tl, a, b) integrate(function(t) cl * exp(-lz * (t - tl)), max(a, tl), b)$value
  expect_equal(m15_pauc_ext(3, 30, cl, lz, 28, 56), 3 + f(30, 28, 56), tolerance = 1e-8)
  expect_equal(m15_pauc_ext(0, 20, cl, lz, 28, 56), f(20, 28, 56), tolerance = 1e-8)
  expect_equal(m15_pauc_ext(c(2, 2), c(57, 40), c(cl, cl), c(lz, NA), 28, 56), c(2, 2))
})

test_that("시험 단위 짝 부트스트랩 AUROC: 점추정은 m15_auroc와 같고 구간이 점추정을 포함한다", {
  set.seed(4); neg <- abs(rnorm(300, 0, 0.05)); pos <- cbind(abs(rnorm(300, 0.05, 0.05)), abs(rnorm(300, -0.05, 0.05)))
  b <- m15_auroc_boot(pos, neg, B = 300, seed = 2)
  expect_equal(b$auroc, m15_auroc(as.vector(pos), neg)$auroc); expect_lt(b$lo, b$auroc); expect_gt(b$hi, b$auroc)
})

test_that("분석계획서(영문): ASCII만, 한글·en/em dash·U+2212 없음, 등록 수치가 설정 파일과 일치", {
  f <- proj_path("plans", "M15_analysis_plan_20261020.md"); skip_if_not(file.exists(f))
  b <- readBin(f, "raw", file.info(f)$size); expect_false(any(as.integer(b) > 127L))
  source(proj_path("regulatory", "src", "reg_helpers.R"), local = TRUE); expect_true(check_english(f))
  x <- paste(readLines(f, warn = FALSE), collapse = " "); cfg <- read_cfg("prereg_20261020_m15.yaml")
  expect_true(grepl(sprintf("%s replicate studies", format(cfg$validation$replicates, big.mark = ",")), x))
  expect_true(grepl("95% PI", x)); expect_true(grepl("99% PI", x)); expect_equal(cfg$validation$failure_rule$R4$pi_level, 0.99)
  expect_true(grepl("mean 167.0 cm, SD 8.68 cm", x)); expect_equal(cfg$populations$S2$height$mean_cm, 167.0)
  expect_true(grepl("so 21 values", x)); expect_equal(as.integer(cfg$oc$gmr_grid$n_total), 21L)
  expect_true(grepl("studies 1 to 500", x)); expect_true(grepl("c = 0.078 mg x day/L", x)); expect_equal(cfg$auroc$pauc_offset_mg_day_L, 0.078)
  expect_true(grepl("2,000 resamples", x)); expect_true(grepl("20,000 log-normal replicate sets", x))
})

test_that("시험 시드 조회: 시험·흐름마다 스칼라 하나, 서로 다르고, 충돌 없는 칸은 derive_seed와 같다(D-068 이탈: 모든 시험이 같은 시드였던 버그)", {
  tab <- m15_trial_seeds(20261020L, 300L); sm <- m15_seed_matrix(tab)
  expect_equal(dim(sm), c(300L, 5L)); expect_false(anyDuplicated(as.vector(sm)) > 0)
  for (jj in c(1L, 2L, 150L, 300L)) for (st in c("subj", "eps")) {
    expect_length(m15_seed_of(sm, jj, st), 1L)
    if (tab$resolved[tab$j == jj & tab$stream == st] == 0L) expect_equal(m15_seed_of(sm, jj, st), derive_seed(20261020L, "trial", jj, st))
  }
  expect_false(m15_seed_of(sm, 1L, "subj") == m15_seed_of(sm, 2L, "subj"))
  big <- m15_seed_matrix(m15_trial_seeds(20261020L, 10000L)); expect_false(anyDuplicated(as.vector(big)) > 0)
  expect_equal(big[1:300, ], sm[1:300, ][, colnames(big)])
})

test_that("가중 Mann-Whitney 재표집 AUROC는 재표집 자료에 직접 계산한 값과 같다(동점 포함)", {
  set.seed(8); nt <- 300; neg <- round(abs(rnorm(nt, 0, 0.05)), 3); pos <- cbind(round(abs(rnorm(nt, 0.03, 0.05)), 3), round(abs(rnorm(nt, -0.02, 0.05)), 3))
  idx <- m15_boot_index(nt, 40, 5)
  expect_equal(m15_auroc_draws_w(pos, neg, idx), m15_auroc_draws(pos, neg, idx), tolerance = 1e-12)
  expect_equal(m15_auroc_draws_w(pos[, 1, drop = FALSE], neg, list(seq_len(nt))), m15_auroc(pos[, 1], neg)$auroc, tolerance = 1e-12)
})

test_that("conditional label and writer: label only after a primary failure, caller table unchanged", {
  td <- tempfile("m15cond"); dir.create(file.path(td, "validation"), recursive = TRUE)
  expect_identical(m15_conditional_label(td, "k2020"), "")
  data.table::fwrite(data.table::data.table(pk_model = c("k2020", "k2020", "k2016"), wmode = c("primary", "wunc", "primary"), failure = c(TRUE, FALSE, FALSE)),
                     file.path(td, "validation", "validation_rule.csv"))
  lab <- m15_conditional_label(td, "k2020"); expect_match(lab, "^Conditional on a primary model")
  expect_identical(m15_conditional_label(td, "k2016"), "")
  x <- data.table::data.table(a = 1:3); f <- file.path(td, "x.csv")
  m15_fwrite_cond(x, f, lab); expect_identical(names(x), "a")
  y <- data.table::fread(f); expect_identical(names(y), c("a", "conditionality")); expect_true(all(y$conditionality == lab))
  m15_fwrite_cond(x, f, ""); expect_identical(names(data.table::fread(f)), "a")
})

test_that("conditional writer adds the 200 mg limitation column only where 200 mg rows exist", {
  f <- tempfile(fileext = ".csv")
  m15_fwrite_cond(data.table::data.table(dose_mg = c(300, 200, NA)), f, "")
  y <- data.table::fread(f); expect_identical(y$limitation[2], M15_LIM200); expect_true(all(is.na(y$limitation[-2]) | y$limitation[-2] == ""))
  m15_fwrite_cond(data.table::data.table(contrast = c("dose", "schedule"), design_a = c("200", "RM36_50")), f, "")
  y <- data.table::fread(f); expect_identical(y$limitation[1], M15_LIM200)
  m15_fwrite_cond(data.table::data.table(contrast = "dose 200 - 300", a = "k2020|S1|200|B0|F|down", b = "k2020|S1|300|B0|F|down"), f, "")
  expect_identical(data.table::fread(f)$limitation, M15_LIM200)
  m15_fwrite_cond(data.table::data.table(dose_mg = 300), f, ""); expect_false("limitation" %in% names(data.table::fread(f)))
})

test_that("체중 효과(개정 1, 계획서 16.4절 1): 모델마다 자기 지수를 쓰고, Vmax가 0에 가까우면 (w/75)^-theta_WT, 두 모델의 비는 다르다", {
  ratio <- function(p, w, dose = 300) {
    ip <- individual_params(p, data.table(id = 1:2, WT = c(w, 75)))
    s <- solve_model(ip, CJ(id = 1:2, time = c(seq(0, 400, by = 0.5), 2000, 4000)), dose, model_id = p$model_id)
    a <- s[time == 4000, auc]; a[1] / a[2] }
  p16 <- load_params("k2016"); p20 <- load_params("k2020")
  expect_equal(p16$cov$theta_WT, 0.705); expect_equal(p20$cov$theta_WT, 0.711)
  lin <- list()
  for (m in c("k2016", "k2020")) { q <- load_params(m); q$theta["Vmax"] <- q$theta[["Vmax"]] * 1e-9          # (a) 선형 소실이면 해석해
    lin[[m]] <- ratio(q, 63.2); expect_equal(lin[[m]], (63.2 / 75)^(-q$cov$theta_WT), tolerance = 1e-6) }
  expect_gt(abs(lin$k2020 - lin$k2016), 5e-4)                                            # (b) 선형 소실에서 두 모델의 비는 5e-4보다 더 다르다
  r16 <- ratio(p16, 63.2); r20 <- ratio(p20, 63.2); expect_gt(abs(r20 - r16), 1e-5)       # 전체 모델: 다른 비(1.206 근처로 가까운 것은 부분 상쇄)
  for (p in list(p16, p20)) { q <- p; q$cov$theta_WT <- q$cov$theta_WT + 0.01; expect_gt(abs(ratio(q, 63.2) - ratio(p, 63.2)), 1e-3) }   # (c) 지수 +0.01이면 비가 1e-3보다 더 바뀐다
  wr <- fread(proj_path("results", "m15", "weight", "weight_ratio.csv"))
  skip_if(!nrow(wr))
  for (m in c("k2016", "k2020")) for (pv in c("S1", "S2")) for (pr in c(0.05, 0.95)) {                    # (d) 저장된 비 = 독립 재계산(1e-6)
    z <- wr[pk_model == m & dose_mg == 300 & population == pv & abs(prob - pr) < 1e-9 & metric == "AUCinf"]
    expect_equal(nrow(z), 1L); expect_equal(ratio(load_params(m), z$WT), z$typical_ratio, tolerance = 1e-6) }
})

test_that("진단 모델(개정 1): k2020cov·k2020base는 Model 1 구조, 공변량·IIV 차용, ke_cov 열이 ke에 곱해지고 기존 모델은 그대로", {
  pc <- load_params("k2020cov"); pb <- load_params("k2020base"); p1 <- load_params("k2020")
  expect_equal(pc$model_id, p1$model_id); expect_equal(pc$cov$theta_WT, 0.817); expect_equal(pb$cov$theta_WT, 0.919)
  expect_equal(pc$omega2, p1$omega2); expect_equal(pc$sigma, p1$sigma); expect_true(all(pc$status$status == "assumption"))
  e <- data.table(id = 1:2, WT = 75, ke_cov = c(1, 1.123))
  ip <- individual_params(pc, e); expect_equal(ip$ke[2] / ip$ke[1], 1.123); expect_true("ktr" %in% names(ip))
  ip1 <- individual_params(p1, data.table(id = 1:2, WT = 75)); expect_equal(ip1$ke, rep(0.0534, 2))
})

test_that("등록 모델 불변(개정 1): k2016·k2020의 파라미터와 개체 파라미터가 커밋 4a948df와 같다(digest)", {
  CFG <- m15_cfg()
  pin <- list(k2016 = c("f39187e7ea23af7997915d219e88a3ac", "9f3ec6f4beb733831f55b8a20a2c83b2"), k2020 = c("5a0593014e273f8a16367a84bdaea398", "2787ed7aba27826aef0c3586dd094cff"))
  for (m in names(pin)) { p <- load_params(m)
    s <- m15_subjects(1000L, CFG$populations$S1, CFG$inclusion, 0.513, 0.5, 20261021L)
    ip <- individual_params(p, m15_eta_table(s, p))
    expect_identical(digest::digest(unclass(p)[c("theta", "cov", "omega2", "omega", "sigma", "model_id")]), pin[[m]][1], info = m)
    expect_identical(digest::digest(as.data.frame(ip)), pin[[m]][2], info = m) }
})

test_that("요청자 NCA 규칙(개정 1, m15_nca_blq0): 투여 뒤 BLQ = 0, 마지막 양의 농도까지, 잘라내기 없음", {
  d <- data.table(id = c(1, 1, 1, 1, 1, 2, 2), time = c(0, 1, 2, 3, 4, 0, 1), conc = c(NA, 4, NA, 2, NA, NA, NA))
  r <- m15_nca_blq0(d)
  expect_equal(r[id == 1, Cmax], 4); expect_equal(r[id == 1, Tlast], 3)
  expect_equal(r[id == 1, AUClast], 0.5 * 4 + 0.5 * 4 + 0.5 * 2)                       # 0->4 선형, 4->0 선형(0 끝점), 0->2 선형
  expect_true(is.na(r[id == 2, AUClast]))
  b <- m15_nca_basic(d); expect_equal(b[id == 1, AUClast], 0.5 * 4 + (4 - 2) / log(2) * 2)   # 등록 규칙: 중간 BLQ 제외, 4->2 로그 하강
})

test_that("파라미터 불확실성(개정 1, R/pu.R): SE 규칙, 차용, 세트 0 = config, 같은 이름은 같은 z", {
  source(proj_path("R", "pu.R"), local = TRUE)
  expect_error(pu_param_set(load_params("k2016"), transform(pu_entries_k2016(pu_cfg()), estimate = ifelse(group == "omega2", -1, estimate)), 1, 20261020L), "10")   # 등록 다시 뽑기 최대 10회(r1~r10)
  PU <- pu_cfg(); E16 <- pu_entries_k2016(PU)
  v <- E16[group == "theta" & name == "Vc"]; expect_equal(v$sd_log, (log(2.97) - log(2.61)) / 3.29)
  w <- E16[group == "theta" & name == "theta_WT"]; expect_equal(w$se_lin, (0.840 - 0.576) / 3.29)
  o <- E16[group == "omega2" & name == "ke"]; expect_equal(o$se_lin, (0.191 - 0.0738) / 3.29)
  E20b <- pu_entries_k2020(PU, "borrowed_k2016_relative"); expect_setequal(E20b$name, c("Vc", "ke", "k12", "k21", "ka", "Vmax", "F", "theta_WT"))   # 정정: Kovalenko 2016 상대 폭 차용
  for (nm in c("Vc", "ke", "k12", "k21", "ka", "Vmax", "F")) expect_equal(E20b[name == nm, sd_log], E16[group == "theta" & name == nm, sd_log])
  expect_equal(E20b[name == "theta_WT", se_lin], (0.840 - 0.576) / 3.29 / 0.705 * 0.711)
  expect_false(any(c("omega2", "sigma") %in% E20b$group))
  for (m in c("k2016", "k2020")) { p <- load_params(m); E <- if (m == "k2016") E16 else E20b; pu_check_estimates(p, E)
    s0 <- pu_param_set(p, E, 0, 20261020L); expect_identical(unclass(s0$p)[c("theta", "cov", "omega2", "omega", "sigma")], unclass(p)[c("theta", "cov", "omega2", "omega", "sigma")]) }
  a <- pu_param_set(load_params("k2016"), E16, 7, 20261020L); b <- pu_param_set(load_params("k2020"), E20b, 7, 20261020L)
  z16 <- log(a$values[["theta.Vc"]] / 2.74) / E16[group == "theta" & name == "Vc", sd_log]; z20 <- log(b$values[["theta.Vc"]] / 2.48) / E20b[name == "Vc", sd_log]
  expect_equal(z16, z20, tolerance = 1e-10)                                            # 같은 이름(Vc)은 모델 간에 같은 점수
  expect_equal(a$p$omega[["ke"]], sqrt(a$p$omega2[["ke"]]))
})

test_that("분석계획서 16절(개정 1)의 등록 수치가 설정 블록 amendment_20261021과 일치", {
  f <- proj_path("plans", "M15_analysis_plan_20261020.md"); skip_if_not(file.exists(f))
  x <- paste(readLines(f, warn = FALSE), collapse = " "); A <- read_cfg("prereg_20261020_m15.yaml")$amendment_20261021
  expect_true(grepl("## 16. Amendment 1 (2026-10-21)", x, fixed = TRUE))
  expect_equal(as.numeric(unlist(A$mechsuppx$multipliers)), c(0.60, 0.70, 1.30, 1.40, 1.50)); expect_true(grepl("multipliers 0.60, 0.70, 1.30, 1.40 and 1.50", x))
  expect_equal(as.numeric(unlist(A$mechsuppx$appended[[1]]$multipliers)), c(1.75, 2.00)); expect_true(grepl("ka also at 1.75 and 2.00", x))
  expect_equal(as.numeric(unlist(A$ada$elapsed_days)), c(14, 28, 56, 63, 70, 84)); expect_true(grepl("elapsed days 14, 28, 56, 63, 70 and 84", x))
  expect_equal(A$ppc$between_study_factor$cv_pct, 10); expect_true(grepl("CV 10% (median 1", x, fixed = TRUE))
  expect_equal(A$ppc$conclusion_template$threshold_fp10_pct, 10); expect_true(grepl("rule met and FP at most 10%", x, fixed = TRUE))
  expect_equal(A$param_uncertainty$n_sets, 1000L); expect_equal(A$param_uncertainty$n_subj_per_set, 2000L); expect_equal(as.integer(unlist(A$param_uncertainty$n_randomized_per_arm)), c(107L, 130L))
  expect_true(grepl("at 107 and 130 randomized per arm (96.3 and 117 evaluable)", x, fixed = TRUE))
  expect_equal(A$robustness$thresholds$pp, 2); expect_equal(A$robustness$thresholds$auroc, 0.02); expect_true(grepl("at least 2 percentage points", x)); expect_true(grepl("at least 0.02 (AUROC, RMSE of the log ratio)", x, fixed = TRUE))
  labs <- vapply(A$ppc$mappings, `[[`, "", "label"); expect_setequal(labs, c("a_model1", "a_k2016")); for (l in labs) expect_true(grepl(l, x, fixed = TRUE), info = l)   # 정정: 등록 두 모델만
  expect_identical(A$labels$in_sample$k2020[["PKM12350"]], "out-of-sample"); expect_true(grepl("PKM12350, PKM14271 and R668-HV-1108 out-of-sample (FDA 2017 review, section 4.6)", x, fixed = TRUE))
  expect_identical(A$ppc$between_arm_check$present_if, "Cochran Q p-value < 0.10"); expect_true(grepl("Cochran Q p-value is below 0.10", x, fixed = TRUE))
  expect_setequal(names(A$param_uncertainty$sets), c("k2016", "k2016_thetas", "k2020_b16")); expect_null(A$applicability$request_range_kg); expect_true(grepl("belongs to the phase 3 data of Models 3 and 4 and is withdrawn (correction of 2026-10-21)", x, fixed = TRUE))
  expect_identical(A$ppc$between_arm_check$present_p_below, 0.10); expect_true(grepl("below 0.10 (fixed here, before Q is computed)", x, fixed = TRUE))
  expect_equal(A$param_uncertainty$redraw$max, 10L)
})

# ---------------------------------------------------------------- 추록 1(계획서 17절, D-071): target-mediated 소실 비중 ----------------------------------------------------------------
test_that("f_TM(추록 1): MM 속도 몫과 같고, C_x에서 정확히 x, C*에서 0.5, 농도에 따라 감소", {
  Vm <- c(1.07, 0.968, 0.5); ke <- c(0.0534, 0.0459, 0.08); Km <- 0.01; C <- c(0.3, 5, 40)
  expect_equal(m15_ftm(C, Vm, ke, Km), (Vm * C / (Km + C)) / (Vm * C / (Km + C) + ke * C), tolerance = 1e-14)
  for (x in c(0.5, 0.8, 0.95)) expect_equal(m15_ftm(m15_tm_conc(x, Vm, ke, Km), Vm, ke, Km), rep(x, 3), tolerance = 1e-13)
  expect_equal(m15_tm_conc(0.5, 1.07, 0.0534, 0.01), 1.07 / 0.0534 - 0.01)
  expect_equal(m15_tm_conc(0.8, 1.07, 0.0534, 0.01), 1.07 / (4 * 0.0534) - 0.01)
  expect_true(all(diff(m15_ftm(seq(0, 100, by = 0.5), 1.07, 0.0534, 0.01)) < 0))
  expect_true(m15_tm_conc(0.999, 0.001, 1, 0.01) < 0)                  # f_TM(0) < x이면 C_x <= 0(도달 불가)
})

test_that("문턱 교차(추록 1): Bateman 곡선에서 하강 교차 시각과 누적 AUC가 해석해와 같고, 상태 분류가 맞다", {
  k1 <- 0.06; k2 <- 0.5; A <- 40; Cf <- function(t) A * (exp(-k1 * t) - exp(-k2 * t)); Af <- function(t) A * ((1 - exp(-k1 * t)) / k1 - (1 - exp(-k2 * t)) / k2)
  tt <- seq(0, 140, by = 0.01); C <- Cf(tt); au <- Af(tt); c0 <- 12
  tm <- optimize(Cf, c(0, 30), maximum = TRUE)$maximum
  t_dn <- uniroot(function(t) Cf(t) - c0, c(tm, 140), tol = 1e-12)$root; t_up <- uniroot(function(t) Cf(t) - c0, c(0, tm), tol = 1e-12)$root
  z <- .m15_tm_cross(tt, C, au, c0)
  expect_identical(z$status, "crossed"); expect_equal(z$n_up, 1L); expect_equal(z$n_dn, 1L)
  expect_equal(z$t_dn, t_dn, tolerance = 1e-6); expect_equal(z$t_up, t_up, tolerance = 1e-4)
  expect_equal(z$a_dn, Af(t_dn), tolerance = 1e-7); expect_equal(z$a_up, Af(t_up), tolerance = 1e-4)   # 상승 교차는 직선 보간(0.01일 격자에서 시각 오차 약 1e-5일)
  expect_identical(.m15_tm_cross(tt, C, au, 100)$status, "throughout"); expect_equal(.m15_tm_cross(tt, C, au, 100)$t_dn, 0)
  expect_identical(.m15_tm_cross(tt, C, au, -1)$status, "never")
  expect_identical(.m15_tm_cross(tt[tt <= 20], C[tt <= 20], au[tt <= 20], c0)$status, "censored")
})

test_that("단회 투여 지표(추록 1, m15_tm_single): AUC0-last 몫이 수치 적분과 같고, Tlast 잘라내기와 throughout 처리가 맞다", {
  k1 <- 0.08; k2 <- 0.6; A <- 35; Cf <- function(t) A * (exp(-k1 * t) - exp(-k2 * t)); Af <- function(t) A * ((1 - exp(-k1 * t)) / k1 - (1 - exp(-k2 * t)) / k2)
  tt <- round(seq(0.01, 140, by = 0.01), 8)
  s <- rbind(data.table(id = 1L, time = tt, C = Cf(tt), auc = Af(tt)), data.table(id = 2L, time = tt, C = Cf(tt) / 10, auc = Af(tt) / 10))
  pr <- data.table(id = 1:2, Vm = 1.07, ke = 0.0534, Km = 0.01); lloq <- 0.078
  o <- m15_tm_single(s, pr, c(0.5, 0.8, 0.95), lloq, t_end = 56, t_max = 140)
  cs <- m15_tm_conc(0.5, 1.07, 0.0534, 0.01)
  tl <- uniroot(function(t) Cf(t) - lloq, c(20, 140), tol = 1e-12)$root; Te <- min(tl, 56)
  num <- integrate(function(t) Cf(t) * (Cf(t) <= cs), 0, Te, subdivisions = 5000, rel.tol = 1e-10)$value
  expect_equal(o[id == 1, share_50], num / Af(Te), tolerance = 1e-6)
  expect_equal(o[id == 1, AUClast], Af(Te), tolerance = 1e-6)
  expect_identical(o[id == 2, status_50], "throughout"); expect_equal(o[id == 2, share_50], 1); expect_equal(o[id == 2, t_50], 0)
  tmx <- optimize(Cf, c(0, 30), maximum = TRUE)$maximum                                 # Cmax 뒤 몫(보충): crossed = t_50 뒤 AUC, throughout = Tmax 뒤 AUC
  t50_ <- uniroot(function(t) Cf(t) - cs, c(tmx, 140), tol = 1e-12)$root
  expect_equal(o[id == 1, share_post_50], (Af(Te) - Af(t50_)) / Af(Te), tolerance = 1e-6)
  Te2 <- min(uniroot(function(t) Cf(t) / 10 - lloq, c(20, 140), tol = 1e-12)$root, 56)   # 대상자 2(농도 1/10)는 참 Tlast가 더 이르다
  expect_equal(o[id == 2, AUClast], Af(Te2) / 10, tolerance = 1e-6)
  expect_equal(o[id == 2, share_post_50], (Af(Te2) - Af(tmx)) / Af(Te2), tolerance = 2e-3)   # Tmax는 0.01일 격자 점(포물선 보정 전)
  t50 <- uniroot(function(t) Cf(t) - cs, c(optimize(Cf, c(0, 30), maximum = TRUE)$maximum, 140), tol = 1e-12)$root
  expect_equal(o[id == 1, t_50], t50, tolerance = 1e-6)
  expect_equal(o[id == 1, ftm_cmax], m15_ftm(optimize(Cf, c(0, 30), maximum = TRUE)$objective, 1.07, 0.0534, 0.01), tolerance = 1e-6)
  rt <- function(t) 1.07 * Cf(t) / (0.01 + Cf(t)); rl <- function(t) 0.0534 * Cf(t)
  expect_equal(o[id == 1, amt_tm_share], integrate(rt, 0, 56, rel.tol = 1e-10)$value / (integrate(rt, 0, 56, rel.tol = 1e-10)$value + integrate(rl, 0, 56, rel.tol = 1e-10)$value), tolerance = 1e-4)   # 보충 지표: 첫 격자 구간(C < Km)의 사다리꼴 오차
})

test_that("정상상태 지표(추록 1, m15_tm_ss): 일정한 농도에서 f_TM, Cavg, 정상상태 비가 해석해와 같다", {
  c0 <- 15; tt <- round(c(98, seq(112, 126, by = 0.01)), 8)
  s <- data.table(id = 1L, time = tt, C = c0, auc = c0 * tt, resid = 1)
  o <- m15_tm_ss(s, data.table(id = 1L, Vm = 1.07, ke = 0.0534, Km = 0.01), c(112, 126), c(98, 112))
  f0 <- m15_ftm(c0, 1.07, 0.0534, 0.01)
  expect_equal(o$Cavg_ss, c0); expect_equal(o$AUCtau_ss, 14 * c0); expect_equal(o$ss_ratio_auc, 1); expect_equal(o$ss_ratio_trough, 1)
  expect_equal(c(o$ftm_trough, o$ftm_cavg, o$ftm_cmax, o$ftm_cmin, o$ftm_tmean), rep(f0, 5), tolerance = 1e-12)
  expect_equal(o$amt_tm_share, f0, tolerance = 1e-12)
})

test_that("투여 일정 풀이(추록 1, m15_solve_regimen): 단회 투여는 등록 단회 풀이(.solve_param_table)와 같고, 투여 시각 전에는 다음 투여가 영향이 없다", {
  skip_if_no_rxode2()
  p <- load_params("k2020"); ip <- individual_params(p, data.table(id = 1:3, WT = c(60, 75, 90), ada = 0, eta_Vmax = c(-0.2, 0, 0.3), eta_ke = c(0.1, 0, -0.2)))
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", "ktr")
  prm <- as.data.frame(ip[, c(pcols, "F"), with = FALSE]); names(prm)[names(prm) == "F"] <- "Fbio"
  tt <- c(0.25, 1, 3, 7, 14, 28, 56, 400)
  a <- m15_solve_regimen(ip, p$model_id, 0, 300, tt); b <- .solve_param_table(get_model(p$model_id), prm, tt, 300)
  expect_equal(a$C, b$C, tolerance = 1e-12); expect_equal(a$auc, b$auc, tolerance = 1e-12); expect_equal(a$id, as.integer(b$id))
  m <- m15_solve_regimen(ip, p$model_id, c(0, 14), c(600, 300), c(1, 7, 13.99, 14 - 1e-6, 14 + 1e-6, 20))
  s6 <- m15_solve_regimen(ip, p$model_id, 0, 600, c(1, 7, 13.99, 20))
  expect_equal(m[time < 13.995, C], s6[time < 13.995, C], tolerance = 1e-10)           # 두 번째 투여 전에는 첫 투여만
  jm <- m[time > 13.995 & time < 14.005][order(id, time)][, .(jump = diff(resid)), by = id]
  expect_equal(jm$jump, 300 * ip$F, tolerance = 1e-4)                                     # 두 번째 투여: 몸 안 총량이 300 x F만큼 는다
  expect_true(all(m[time == 20, C] > s6[time == 20, C] + 1))                              # 두 번째 투여 뒤 농도가 첫 투여만보다 높다
})

test_that("추록 1 설정(addendum_tm_20261021)과 계획서 17절이 일치하고 채혈일 환산이 등록 일정과 같다", {
  f <- proj_path("plans", "M15_analysis_plan_20261020.md"); skip_if_not(file.exists(f))
  x <- paste(readLines(f, warn = FALSE), collapse = " "); C <- read_cfg("prereg_20261020_m15.yaml"); A <- C$addendum_tm_20261021
  expect_true(grepl("## 17. Addendum 1 (2026-10-21)", x, fixed = TRUE)); expect_identical(A$decision, "D-071"); expect_identical(A$plan_section, "17")
  expect_equal(A$n_subjects, 1000L); expect_equal(as.numeric(unlist(A$weights_kg)), c(60, 75, 90)); expect_equal(as.numeric(unlist(A$thresholds)), c(0.50, 0.80, 0.95))
  expect_equal(as.numeric(unlist(A$sampling_days$elapsed)), c(0, as.numeric(unlist(C$schedules$B0))))
  expect_equal(as.integer(unlist(A$sampling_days$study_day)), c(0L, 1L, 2L, 4L, 6L, 8L, 11L, 15L, 22L, 29L, 36L, 43L, 50L, 57L))
  expect_equal(as.numeric(unlist(A$scenarios$B$dose_times)), c(0, seq(14, 112, by = 14))); expect_equal(as.numeric(unlist(A$scenarios$B$dose_mg)), c(600, rep(300, 8)))
  expect_equal(as.numeric(unlist(A$scenarios$B$ss_interval)), c(112, 126)); expect_true(grepl("days 112 to 126 after the first dose", x, fixed = TRUE))
  expect_equal(as.numeric(unlist(A$scenarios$A$grid$fine)), c(0.01, 100, 0.01)); expect_true(grepl("Every 0.01 day to day 100, every 0.05 day to day 140", x, fixed = TRUE))
  expect_equal(A$mc_error$n_boot, 1000L); expect_true(grepl("1,000 resamples, derive_seed(20261020, \"tm_boot\")", x, fixed = TRUE))
  expect_setequal(names(A$checks), letters[1:7])
  expect_equal(A$post_hoc_true_ss$dose_times_last, 2240L); expect_equal(as.numeric(unlist(A$post_hoc_true_ss$interval)), c(2240, 2254))   # 17.11절 사후 점검(첫 시도 1,120일 수렴 실패 뒤 연장, 17.11절 3)
  expect_true(grepl("### 17.11 Changes after the first run (2026-10-21)", x, fixed = TRUE)); expect_true(grepl("continued to day 2,240 (160 repeat doses)", x, fixed = TRUE))
  expect_equal(as.numeric(unlist(A$elasticity$multipliers)), c(0.8, 1.2)); expect_identical(as.character(unlist(A$elasticity$parameters)), c("Vmax", "ke", "F"))
  expect_true(grepl("n = 1,000 virtual subjects", x, fixed = TRUE)); expect_true(grepl("derive_seed(20261020, \"tm\")", x, fixed = TRUE))
  expect_true(grepl("60, 75 and 90 kg", x, fixed = TRUE)); expect_equal(as.numeric(A$lloq_mg_L), study_lloq())
})
