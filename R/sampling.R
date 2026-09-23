# sampling.R — 채혈 일정과 허용창 내 무작위 편차 (SPEC §5.1)

get_schedule <- function(design, name = NULL) {
  if (is.null(name)) name <- design$default_schedule
  s <- design$schedules[[name]]
  if (is.null(s)) stop("알 수 없는 채혈 일정: ", name)
  as.numeric(s$days)
}

window_for_time <- function(t, design) {
  rules <- design$sampling_window_days$rule
  w <- vapply(t, function(tt) {
    for (r in rules) if (tt <= r$max_day) return(as.numeric(r$window))
    as.numeric(rules[[length(rules)]]$window)
  }, numeric(1))
  w
}

# 계획 시각 → 실제 시각(계획 + U(-w, w)); 순서가 뒤바뀌지 않도록 정렬 후 최소 간격 보장
jitter_times <- function(planned, design, min_gap = 1 / 24) {
  w <- window_for_time(planned, design)
  actual <- planned + runif(length(planned), -w, w)
  actual[planned == 0] <- 0
  actual <- pmax(actual, 0)
  o <- order(planned)
  a <- actual[o]
  for (k in seq_along(a)[-1]) if (a[k] < a[k - 1] + min_gap) a[k] <- a[k - 1] + min_gap
  actual[o] <- a
  actual
}

# 모든 피험자의 관측 시각 표: data.table(id, planned, time)
make_obs_times <- function(ids, planned, design, jitter = TRUE) {
  rbindlist(lapply(ids, function(i) {
    data.table(id = i, planned = planned, time = if (jitter) jitter_times(planned, design) else planned)
  }))
}
