test_that("dev 모드는 PENDING 값을 DEV-가정으로 채우고 상태를 표시한다", {
  p <- load_params("dev")
  expect_s3_class(p, "dupi_params")
  expect_equal(unname(p$theta[c("Vc", "ke", "k12", "k21", "F", "Vmax", "Km")]),
               c(2.74, 0.0459, 0.0652, 0.129, 0.607, 0.968, 0.01))
  expect_true(is.finite(p$theta[["ka"]]))
  expect_equal(p$status[item == "theta.ka", status], "dev_assumption")
  expect_true(all(p$status[grepl("^theta\\.(Vc|ke|k12|k21|F|Vmax|Km)$", item), status] == "confirmed"))
})

test_that("final 모드는 출처 미확정 값이 있으면 중단한다 (DECISIONS D-005)", {
  p <- load_params("final")
  expect_true(is.na(p$theta[["ka"]]))
  expect_error(assert_final_ready(p), "출처 미확정")
  # 대표값만 필요한 경우(theta 중 ka 제외)는 통과해야 하므로 need를 좁혀 확인
  p2 <- p; p2$status <- p2$status[!item %in% c("theta.ka")]
  expect_error(assert_final_ready(p2, need = "theta"), NA)
})

test_that("시나리오 승수는 지정 파라미터에만 곱해진다", {
  p0 <- load_params("dev"); p1 <- load_params("dev", scenario_multipliers = list(F = 0.9, Km = 2))
  expect_equal(p1$theta[["F"]], p0$theta[["F"]] * 0.9)
  expect_equal(p1$theta[["Km"]], 0.02)
  expect_equal(p1$theta[["ke"]], p0$theta[["ke"]])
})

test_that("CV <-> omega 변환", {
  expect_equal(omega_to_cv(cv_to_omega(30)), 30)
  expect_equal(cv_to_omega(0), 0)
})
