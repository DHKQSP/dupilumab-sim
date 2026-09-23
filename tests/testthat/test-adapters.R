test_that("BEmaster 부재 시 어댑터는 명확한 오류를 내고 자체 계산을 하지 않는다", {
  if (bemaster_available()) skip("BEmaster가 있어 부재 시 동작은 검사하지 않음")
  sim <- data.table(id = 1L, planned = c(0, 1), time = c(0, 1), C = c(0, 1), auc = c(0, 0.5), y_raw = c(0, 1), blq = c(TRUE, FALSE), conc = c(NA, 1))
  expect_error(run_nca_bemaster(sim, 300, 0.078), "BEmaster")
  expect_error(run_be_bemaster(data.table(), ), "BEmaster")
})

test_that("반환 열 계약이 정의되어 있다", {
  expect_true(all(c("AUClast", "AUCinf", "pct_extrap", "lambda_ok") %in% NCA_COLUMNS))
  expect_true(all(c("GMR", "CI_lower", "CI_upper", "pass") %in% BE_COLUMNS))
})
