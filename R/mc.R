# mc.R — 2층 몬테카를로 (SPEC §5.2). 바깥층: 파라미터 불확실성, 안쪽층: 시험 반복.
# 시드: master → (scenario, outer i) → (scenario, outer i, inner j)  (R/seeds.R)

# 바깥층 theta 추출: 로그 척도 정규. covariance 있으면 사용, 없으면 RSE 대각.
draw_outer_theta <- function(p, seed) {
  th <- p$theta
  nm <- names(p$rse_pct)
  with_seed(seed, {
    if (!is.null(p$covariance)) {
      stopifnot(all(rownames(p$covariance) %in% names(th)))
      nm <- rownames(p$covariance)
      z <- MASS::mvrnorm(1, mu = log(th[nm]), Sigma = p$covariance)
      th[nm] <- exp(z)
    } else {
      rse <- p$rse_pct[nm] / 100; rse[is.na(rse)] <- 0
      sdlog <- sqrt(log(1 + rse^2))
      th[nm] <- th[nm] * exp(rnorm(length(nm), 0, sdlog))
    }
  })
  # F는 (0,1) 범위 유지: logit 척도로 재투영하지 않고 상한 절단(문서화)
  th["F"] <- min(th[["F"]], 0.999)
  th
}

run_mc <- function(scenario, p_base, design, n_outer, n_inner, master_seed,
                   schedule_name = NULL, n_per_arm = NULL, endpoints = c("AUClast", "AUCinf", "Cmax"),
                   nca_fn = run_nca_bemaster, be_fn = run_be_bemaster, do_be = TRUE,
                   progress = TRUE, parallel = FALSE) {
  pp <- apply_scenario(p_base, scenario)
  one_outer <- function(i) {
    seed_i <- derive_seed(master_seed, scenario$code, "outer", i)
    th_i <- draw_outer_theta(p_base, seed_i)
    pR <- pp$R; pT <- pp$T
    # 시나리오 승수를 유지한 채 바깥층 theta 반영
    mult <- pp$T$theta / pp$R$theta
    pR$theta <- th_i; pT$theta <- th_i * mult
    res <- lapply(seq_len(n_inner), function(j) {
      seed_ij <- derive_seed(master_seed, scenario$code, "outer", i, "inner", j)
      tr <- run_trial(pR, pT, design, seed = seed_ij, schedule_name = schedule_name, n_per_arm = n_per_arm,
                      nca_fn = nca_fn, be_fn = be_fn, endpoints = endpoints, do_be = do_be)
      list(outer = i, inner = j, seed = seed_ij, theta = th_i,
           subjects = tr$subjects[, `:=`(outer = i, inner = j)],
           nca = if (!is.null(tr$nca)) tr$nca[, `:=`(outer = i, inner = j)] else NULL,
           be  = if (!is.null(tr$be))  tr$be[,  `:=`(outer = i, inner = j)] else NULL)
    })
    if (progress) cat(sprintf("  [%s] outer %d/%d done\n", scenario$code, i, n_outer))
    res
  }
  outer_res <- if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(n_outer), one_outer, future.seed = NULL)
  } else lapply(seq_len(n_outer), one_outer)
  flat <- unlist(outer_res, recursive = FALSE)
  list(
    scenario = scenario$code, label = scenario$label, master_seed = master_seed,
    n_outer = n_outer, n_inner = n_inner,
    theta_outer = rbindlist(lapply(outer_res, function(o) as.data.table(as.list(o[[1]]$theta))[, outer := o[[1]]$outer])),
    subjects = rbindlist(lapply(flat, `[[`, "subjects")),
    nca = rbindlist(lapply(flat, `[[`, "nca")),
    be = rbindlist(lapply(flat, `[[`, "be"))
  )
}
