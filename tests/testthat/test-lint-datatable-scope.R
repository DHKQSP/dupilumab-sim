# 정적 검사: data.table `[` 안에서 함수 인자와 같은 이름의 열을 쓰지 않는다 (검토 의견 3차 §4, D-022 재발 방지).
# 함수 인자 이름이 알려진 열 이름과 같고, 그 이름이 함수 본문의 `[` 호출 인자(첫 번째 대상 제외) 안에 기호로 나타나면 위반.
# `..var`, `x[["col"]]`, 문자열은 허용된다. 열 이름 목록은 이 저장소 표들의 열 이름.
KNOWN_COLUMNS <- c("id", "time", "planned", "conc", "C", "auc", "blq", "y_raw", "arm", "WT", "sex", "stratum", "ada",
                   "scenario", "schedule", "endpoint", "method", "trial", "variant", "dose", "dose_mg", "tlast", "tmax", "Cmax",
                   "AUClast", "AUCinf", "ka", "ke", "k12", "k21", "Vmax", "Km", "F", "Vc", "Vc_i", "pass", "width", "GMR",
                   "reliable", "lambda_ok", "lambda_z", "n_lambda", "pct_extrap", "n", "N", "cohort", "outer", "inner", "weight_shape", "ref")

collect_symbols <- function(e) {
  if (is.symbol(e)) return(as.character(e))
  if (is.call(e)) {
    fn <- e[[1]]
    if (identical(fn, as.name("[[")) || identical(fn, as.name("$"))) return(unlist(lapply(as.list(e)[2], collect_symbols)))   # x[["col"]], x$col 의 필드는 제외
    if (identical(fn, as.name("function"))) return(character(0))
    return(unlist(lapply(as.list(e)[-1], collect_symbols)))
  }
  character(0)
}
find_violations <- function(body, formals_) {
  bad <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("["))) {
        syms <- unique(unlist(lapply(as.list(e)[-(1:2)], collect_symbols)))
        v <- intersect(syms, intersect(formals_, KNOWN_COLUMNS))
        if (length(v)) bad <<- c(bad, paste0(paste(v, collapse = ","), " in ", paste(deparse(e, width.cutoff = 60)[1], collapse = "")))
      }
      for (a in as.list(e)[-1]) if (!missing(a)) walk(a)
    }
  }
  walk(body); bad
}

# 검토한 알려진 예외(정확한 위반 문자열로만 허용). 대상 표에 그 이름의 열이 없어 함수 인자로 해석되는 것이 확실한 경우만 둔다.
# 해당 줄이 고쳐지거나 바뀌면 이 목록은 아무것도 허용하지 않게 되므로(문자열 불일치) 다시 검사 대상이 된다. 고친 뒤 지운다.
LINT_ALLOW <- character(0)   # 허용 목록(정확한 문자열). 비워 두는 것이 원칙

test_that("R/와 scripts/의 함수에서 인자-열 이름 충돌이 없다", {
  files <- c(list.files(proj_path("R"), pattern = "\\.R$", full.names = TRUE), list.files(proj_path("scripts"), pattern = "\\.R$", full.names = TRUE))
  viol <- character(0)
  for (f in files) {
    exprs <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    visit <- function(e) {
      if (is.call(e)) {
        if (identical(e[[1]], as.name("function"))) {
          fm <- names(e[[2]]); v <- find_violations(e[[3]], fm)
          if (length(v)) viol <<- c(viol, paste0(basename(f), ": ", v))
        }
        for (a in as.list(e)[-1]) if (!missing(a)) visit(a)
      }
    }
    for (e in exprs) visit(e)
  }
  viol <- setdiff(viol, LINT_ALLOW)
  if (length(viol)) message(paste(viol, collapse = "\n"))
  expect_length(viol, 0)
})

test_that("검사기가 실제로 위반을 잡는다(자기 검증)", {
  bad_fn <- quote(function(be, method) be[be$method == method])
  expect_gt(length(find_violations(bad_fn[[3]], names(bad_fn[[2]]))), 0)
  ok_fn <- quote(function(be, method) { m_ <- method; be[be[["method"]] == m_] })
  expect_length(find_violations(ok_fn[[3]], names(ok_fn[[2]])), 0)
})
