# pu.R — 개정 1(계획서 16.8절, D-069) 파라미터 불확실성: 입력 해석(P5/P95 또는 표준오차), SE 규칙, 파라미터 세트 추출.
# 입력 config/param_uncertainty_m15.yaml. scripts/76_m15_param_assurance.R와 tests/testthat/test-m15.R가 직접 source한다(source_project()에 넣지 않는다).
# 기존 scripts/67(Day 57 지시, RSE 형식)은 이 파일을 쓰지 않는다.

# z 점수의 고정 이름 목록(모델·세트와 무관): 같은 이름의 파라미터는 모델 간에 같은 z를 쓴다(공통 난수)
PU_NAMES <- c("theta.Vc", "theta.ke", "theta.k12", "theta.k21", "theta.ka", "theta.Vmax", "theta.F", "theta.theta_WT",
              "omega2.Vc", "omega2.ke", "omega2.ka", "omega2.Vmax", "sigma.prop")

pu_cfg <- function(file = "param_uncertainty_m15.yaml") read_cfg(file)

# k2016: P5/P95 → 로그 척도 sd_log 또는 선형 척도 se_lin. 반환: group, name, estimate, scale, sd_log, se_lin, p05, p95
pu_entries_k2016 <- function(PU) {
  m <- PU$models$k2016; div <- as.numeric(PU$se_rule$divisor)
  E <- rbindlist(lapply(c("theta", "omega2", "sigma"), function(g) rbindlist(lapply(names(m[[g]]), function(nm) { e <- m[[g]][[nm]]
    if (isTRUE(e$fixed)) return(NULL)
    p05 <- as.numeric(e$p05); p95 <- as.numeric(e$p95); est <- as.numeric(e$estimate)
    if (!(p05 <= est && est <= p95)) stop(sprintf("k2016 %s.%s: 추정값 %g가 [P5, P95] = [%g, %g] 밖", g, nm, est, p05, p95))
    data.table(group = g, name = nm, estimate = est, scale = e$scale, p05 = p05, p95 = p95,
               sd_log = if (e$scale == "log") (log(p95) - log(p05)) / div else NA_real_, se_lin = if (e$scale == "linear") (p95 - p05) / div else NA_real_) }))))
  stopifnot(all(E$scale %in% c("log", "linear")), all(paste(E$group, E$name, sep = ".") %in% PU_NAMES))
  E[]
}

# k2020(Model 1): 세트 이름(borrowed_covariate_model 또는 borrowed_covariate_model_and_k2016) → 같은 형식
pu_entries_k2020 <- function(PU, set) {
  S <- PU$models$k2020$sets; s <- S[[set]]; if (is.null(s)) stop("k2020 세트 없음: ", set)
  base <- if (!is.null(s$include)) pu_entries_k2020(PU, s$include) else NULL
  D <- pu_entries_k2016(PU)
  E <- rbindlist(lapply(names(s$theta), function(nm) { e <- s$theta[[nm]]; est <- as.numeric(e$estimate)
    sdl <- NA_real_; sel <- NA_real_
    if (!is.null(e$sd_log_from)) sdl <- as.numeric(e$sd_log_from$se) / as.numeric(e$sd_log_from$estimate)          # 상대 SE를 로그 척도 SD로
    if (!is.null(e$se)) sel <- as.numeric(e$se)                                                                       # 절대 SE(선형 척도)
    if (!is.null(e$sd_log_borrow_k2016)) { d <- D[group == "theta" & name == e$sd_log_borrow_k2016]; stopifnot(nrow(d) == 1, d$scale == "log"); sdl <- d$sd_log }
    data.table(group = "theta", name = nm, estimate = est, scale = e$scale, p05 = NA_real_, p95 = NA_real_, sd_log = sdl, se_lin = sel) }))
  if (!is.null(base)) E <- rbind(base[!base$name %in% E$name], E)
  stopifnot(all((E$scale == "log" & is.finite(E$sd_log) & E$sd_log > 0) | (E$scale == "linear" & is.finite(E$se_lin) & E$se_lin > 0)),
            all(paste(E$group, E$name, sep = ".") %in% PU_NAMES))
  E[]
}

# 추정값이 모델 config와 같은지(상대 1e-9; 입력 오타·단위 실수 검출)
pu_check_estimates <- function(p, E) {
  cur <- function(g, nm) switch(g, theta = if (nm == "theta_WT") p$cov$theta_WT else p$theta[[nm]], omega2 = p$omega2[[nm]], sigma = p$sigma[[nm]])
  for (i in seq_len(nrow(E))) { v <- cur(E$group[i], E$name[i])
    if (abs(v / E$estimate[i] - 1) > 1e-9) stop(sprintf("추정값 불일치 %s.%s: config %g, 입력 %g", E$group[i], E$name[i], v, E$estimate[i])) }
  invisible(TRUE)
}

# 세트 s의 파라미터: z는 derive_seed(ms, "pset", s)에서 PU_NAMES 순서로(모델·세트와 무관), 로그 척도는 estimate x exp(z sd_log), 선형은 estimate + z se_lin.
# 분산 또는 잔차가 0 이하이거나 F >= 1이면 세트 전체를 하위 흐름 "r1", "r2", ...로 다시 뽑고 횟수를 센다. s = 0은 추정값 그대로(config와 같음).
pu_param_set <- function(p, E, s, ms) {
  key <- paste(E$group, E$name, sep = "."); r <- 0L
  repeat {
    z <- if (s == 0) setNames(rep(0, length(PU_NAMES)), PU_NAMES) else
      with_seed(if (r == 0) derive_seed(ms, "pset", s) else derive_seed(ms, "pset", s, paste0("r", r)), setNames(rnorm(length(PU_NAMES)), PU_NAMES))
    val <- ifelse(E$scale == "log", E$estimate * exp(z[key] * E$sd_log), E$estimate + z[key] * E$se_lin)
    Fv <- if (any(key == "theta.F")) val[key == "theta.F"] else p$theta[["F"]]
    if (all(val[E$group %in% c("omega2", "sigma")] > 0) && Fv > 0 && Fv < 1) break
    r <- r + 1L; if (r > 20) stop("세트 ", s, ": 다시 뽑기 20회 초과")
  }
  ps <- p
  for (i in seq_len(nrow(E))) { g <- E$group[i]; nm <- E$name[i]; v <- val[i]
    if (g == "theta") { if (nm == "theta_WT") ps$cov$theta_WT <- v else ps$theta[[nm]] <- v }
    if (g == "omega2") { ps$omega2[[nm]] <- v; ps$omega[[nm]] <- sqrt(v) }
    if (g == "sigma") ps$sigma[[nm]] <- v }
  list(p = ps, values = as.list(setNames(val, key)), redraws = r)
}
