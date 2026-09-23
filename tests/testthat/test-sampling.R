design <- read_cfg("trial_design.yaml")

test_that("채혈 일정 정의(지시서 §5)", {
  expect_equal(get_schedule(design, "B0"), c(0.25, 1, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49, 56))
  expect_equal(setdiff(get_schedule(design, "B0"), get_schedule(design, "Bminus")), 49)
  expect_equal(setdiff(get_schedule(design, "D1"), get_schedule(design, "B0")), c(38, 45))
  expect_equal(setdiff(get_schedule(design, "D2"), get_schedule(design, "B0")), c(38, 45, 52))
  expect_equal(setdiff(get_schedule(design, "D3"), get_schedule(design, "B0")), c(31, 38, 45, 52))
  expect_equal(setdiff(get_schedule(design, "D4"), get_schedule(design, "B0")), c(39, 46))
  expect_false("d85_ext" %in% names(design$schedules))
  g <- union_grid(design, design$schedule_analysis)
  expect_equal(g, sort(unique(c(0, 0.25, 1, 3, 5, 7, 10, 14, 21, 28, 31, 35, 38, 39, 42, 45, 46, 49, 52, 56))))
})

test_that("채혈 편차는 허용창 안, 시점 0 고정, 순서 유지", {
  planned <- union_grid(design, "B0")
  for (k in 1:30) {
    a <- with_seed(k, jitter_times(planned, design)); w <- window_for_time(planned, design)
    expect_true(all(abs(a - planned) <= w + 1e-9)); expect_equal(a[planned == 0], 0); expect_true(all(diff(a[order(planned)]) > 0))
  }
  o <- make_obs_times(1:3, get_schedule(design, "B0"), design, jitter = FALSE)
  expect_equal(nrow(o), 3 * 14); expect_true(all(o$time == o$planned))
})
