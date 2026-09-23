test_that("시드 파생은 결정적이고 인자에 민감하다", {
  a <- derive_seed(20260923, "S0", "outer", 1, "inner", 2)
  b <- derive_seed(20260923, "S0", "outer", 1, "inner", 2)
  c <- derive_seed(20260923, "S0", "outer", 1, "inner", 3)
  expect_identical(a, b); expect_false(a == c)
  expect_true(is.integer(a) && a >= 0)
})

test_that("with_seed는 바깥 RNG 상태를 복원한다", {
  set.seed(1); x1 <- runif(1)
  set.seed(1); y <- with_seed(99, runif(3)); x2 <- runif(1)
  expect_equal(x1, x2)
  expect_equal(y, with_seed(99, runif(3)))
})
