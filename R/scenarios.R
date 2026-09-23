# scenarios.R — 시나리오 정의 로드 (SPEC §5.5)
load_scenarios <- function() {
  s <- read_cfg("scenarios.yaml")
  main <- lapply(names(s$scenarios), function(k) {
    x <- s$scenarios[[k]]
    list(code = k, label = x$label, T_multipliers = if (length(x$T_multipliers)) x$T_multipliers else list(),
         both_arms = list(), both_arms_iiv = list(), status = "confirmed")
  })
  names(main) <- names(s$scenarios)
  sens <- lapply(names(s$sensitivity), function(k) {
    x <- s$sensitivity[[k]]
    list(code = k, label = x$label, T_multipliers = list(),
         both_arms = if (length(x$both_arms)) x$both_arms else list(),
         both_arms_iiv = if (length(x$both_arms_iiv)) x$both_arms_iiv else list(),
         status = if (!is.null(x$status)) x$status else "confirmed")
  })
  names(sens) <- names(s$sensitivity)
  list(main = main, sensitivity = sens)
}

# 시나리오를 파라미터 객체 쌍(R, T)에 적용
apply_scenario <- function(p, scen) {
  pR <- p; pT <- p
  for (nm in names(scen$both_arms)) {
    v <- scen$both_arms[[nm]]
    if (is.null(v)) stop("시나리오 ", scen$code, ": ", nm, " 값이 PENDING 입니다 (SPEC Q8)")
    pR$theta[nm] <- v; pT$theta[nm] <- v
  }
  for (nm in names(scen$both_arms_iiv)) {
    v <- scen$both_arms_iiv[[nm]]
    if (is.null(v)) stop("시나리오 ", scen$code, ": ", nm, " 값이 PENDING 입니다 (SPEC Q8)")
    par <- sub("_omega$", "", nm)
    pR$omega[par] <- v; pT$omega[par] <- v
  }
  for (nm in names(scen$T_multipliers)) pT$theta[nm] <- pT$theta[nm] * scen$T_multipliers[[nm]]
  list(R = pR, T = pT)
}
