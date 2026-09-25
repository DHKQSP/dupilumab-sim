# 연구 LLOQ 민감도(추가 지시 2026-09-26 §2, R/lloq.R, R/cliff.R, R/oc_models.R recensor_obs·rescale_additive): 재검열·잔차 배율·절벽 LLOQ 벡터의 항등성

test_that("재검열은 모의 LLOQ에서 원래 관측과 같고, 잔차 배율 1은 y_raw를 그대로 재현하며, y_raw는 BLQ 행에도 남는다", {
  skip_if_no_rxode2()
  design <- read_cfg("trial_design.yaml"); rv <- resolve_variant("base", design); p <- rv$p
  pop <- run_individual_population(40, p, design, "B0", 11L, rv$wt_spec, jitter = TRUE, model_id = p$model_id, tag = "t_lloq")
  expect_equal(pop$lloq, study_lloq())
  expect_false(anyNA(pop$obs$y_raw))
  expect_identical(recensor_obs(pop$obs, pop$lloq), pop$obs)
  eps <- with_seed(derive_seed(11L, "t_lloq", "eps"), draw_eps(pop$subj$id, sort(unique(c(0, pop$grid))), p$sigma))
  r1 <- rescale_additive(pop$obs, eps, 1)
  expect_identical(r1$y_raw, pop$obs$y_raw)
  r2 <- rescale_additive(pop$obs, eps, 0.5)
  expect_identical(r2$C, pop$obs$C); expect_false(identical(r2$y_raw, pop$obs$y_raw))
  bad <- copy(eps)[, eps_a := eps_a + 1e-9]
  expect_error(rescale_additive(pop$obs, bad, 1), "재현")
  li <- lloq_individual(pop, design, p, c(0.02, study_lloq()), c("fixed", "scaled"), 11L, "t_lloq")
  s <- li$summary
  expect_equal(s[resid == "scaled" & lloq == study_lloq(), reliable_pct], s[resid == "fixed" & lloq == study_lloq(), reliable_pct])
  expect_true(all(li$nca[, !is.na(tlast_planned) | is.na(tlast)]))
  expect_true(all(li$nca[!is.na(tlast), tlast_planned %in% c(0, get_schedule(design, "B0"))]))
})

test_that("절벽 도달 시각은 LLOQ를 벡터로 줘도 하나씩 줄 때와 같고 LLOQ가 낮을수록 늦다", {
  skip_if_no_rxode2()
  skip_if_fast_scope("절벽 곡선 모의 20명, 약 5초")
  design <- read_cfg("trial_design.yaml"); cf <- read_cfg("oc_design.yaml")$cliff
  a <- cliff_subjects("k2016", "base", 20L, 5L, 0.05, c(0.02, 0.078, 0.5), c(1, 2), cf$weights$base, design)
  b <- cliff_subjects("k2016", "base", 20L, 5L, 0.05, 0.078, c(1, 2), cf$weights$base, design)
  expect_identical(a[lloq == 0.078, t_lloq], b$t_lloq)
  w <- dcast(a, id ~ lloq, value.var = "t_lloq")
  expect_true(all(w[["0.02"]] > w[["0.078"]] & w[["0.078"]] > w[["0.5"]]))
})
