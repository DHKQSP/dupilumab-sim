# seeds.R — 결정적 시드 파생 (DECISIONS D-006)
# master seed → outer i → inner j 를 digest::digest 기반으로 파생한다.
# 규칙 버전은 로그에 기록되며, 바뀌면 SEED_RULE_VERSION을 올린다.
SEED_RULE_VERSION <- "1"

derive_seed <- function(master, ...) {
  key <- paste(c("seedrule", SEED_RULE_VERSION, master, ...), collapse = "|")
  h <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
  as.integer(strtoi(substr(h, 1, 7), 16L))   # 0 .. 2^28-1, R 정수 범위 안
}

with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}
