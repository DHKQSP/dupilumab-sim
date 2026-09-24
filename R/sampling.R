# sampling.R — 채혈 일정(투여 후 경과일)과 허용창 편차 (지시서 §4–5)

get_schedule <- function(design, name = NULL) {
  if (is.null(name)) name <- design$default_schedule
  s <- design$schedules[[name]]
  if (is.null(s)) stop("알 수 없는 채혈 일정: ", name)
  sort(as.numeric(s$days))
}

# 여러 일정의 합집합 격자(투여 전 0 포함)
union_grid <- function(design, names) sort(unique(c(0, unlist(lapply(names, function(n) get_schedule(design, n))))))

window_for_time <- function(t, design) {
  rules <- design$sampling_window_days$rule
  vapply(t, function(tt) { for (r in rules) if (tt <= r$max_day) return(as.numeric(r$window)); as.numeric(rules[[length(rules)]]$window) }, numeric(1))
}

# 계획 → 실제 시각. 시점 0은 고정, 순서 유지(최소 간격 1시간)
jitter_times <- function(planned, design, min_gap = 1 / 24) {
  pl <- planned                                        # 인자명을 부분집합 조건에 직접 쓰지 않는다(D-025)
  w <- window_for_time(pl, design)
  actual <- pl + runif(length(pl), -w, w)
  actual[pl == 0] <- 0
  actual <- pmax(actual, 0)
  o <- order(pl); a <- actual[o]
  for (k in seq_along(a)[-1]) if (a[k] < a[k - 1] + min_gap) a[k] <- a[k - 1] + min_gap
  actual[o] <- a
  actual
}

# 모든 피험자의 관측 시각 표: data.table(id, planned, time)
make_obs_times <- function(ids, planned, design, jitter = TRUE) {
  pl <- sort(unique(c(0, planned)))
  if (!jitter) { g <- CJ(id = ids, planned = pl); g[, time := g[["planned"]]]; return(g[]) }
  rbindlist(lapply(ids, function(i) data.table(id = i, planned = pl, time = jitter_times(pl, design))))
}
