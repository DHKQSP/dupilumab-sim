design <- read_cfg("trial_design.yaml")

test_that("채혈 편차는 허용창 안이고 시점 0은 고정, 순서 유지", {
  planned <- get_schedule(design, "d57_base")
  for (k in 1:50) {
    a <- with_seed(k, jitter_times(planned, design))
    w <- window_for_time(planned, design)
    expect_true(all(abs(a - planned) <= w + 1e-9))
    expect_equal(a[planned == 0], 0)
    expect_true(all(diff(a[order(planned)]) > 0))
  }
})

test_that("일정 이름별 채혈 일수", {
  expect_equal(max(get_schedule(design, "d57_base")), 56)
  expect_equal(max(get_schedule(design, "d85_ext")), 84)
  expect_gt(length(get_schedule(design, "d57_dense")), length(get_schedule(design, "d57_base")))
  expect_error(get_schedule(design, "nope"))
})
