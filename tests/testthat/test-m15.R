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
