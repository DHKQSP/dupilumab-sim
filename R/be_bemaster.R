# be_bemaster.R — 시험 수준 BE 판정 어댑터(로그 변환 ANOVA + 체중 층, GMR 90% CI, 80-125%).
# 규칙과 계산은 BEmaster가 수행한다. 여기서는 데이터 형식 변환만 한다 (SPEC §5.4, §6).
#
# 반환 계약: data.table(endpoint, GMR, CI_lower, CI_upper, pass, n_R, n_T, method)

BE_COLUMNS <- c("endpoint", "GMR", "CI_lower", "CI_upper", "pass", "n_R", "n_T", "method")

.bemaster_call_be <- function(dt, endpoint, ci_level, limits, strata_var) {
  stop("BEmaster BE 인터페이스 바인딩이 아직 구현되지 않았습니다. R/be_bemaster.R::.bemaster_call_be 를 BEmaster 함수 시그니처에 맞춰 작성하세요 (SPEC §6).", call. = FALSE)
}

# nca: run_nca_bemaster() 출력 + arm, WT_stratum 열. endpoints: 분석할 열 이름들.
run_be_bemaster <- function(nca, endpoints = c("AUClast", "AUCinf", "Cmax"),
                            ci_level = 0.90, limits = c(0.80, 1.25), strata_var = "WT_stratum") {
  if (!bemaster_available()) stop_bemaster_missing()
  stopifnot(all(c("arm", strata_var) %in% names(nca)))
  res <- lapply(endpoints, function(ep) {
    d <- nca[!is.na(get(ep)) & get(ep) > 0, c("id", "arm", strata_var, ep), with = FALSE]
    r <- .bemaster_call_be(d, endpoint = ep, ci_level = ci_level, limits = limits, strata_var = strata_var)
    r <- as.data.table(as.list(r)); r[, endpoint := ep]
    r[, `:=`(n_R = sum(d$arm == "R"), n_T = sum(d$arm == "T"))]
    r
  })
  out <- rbindlist(res, fill = TRUE)
  miss <- setdiff(BE_COLUMNS, names(out))
  if (length(miss)) stop("BEmaster BE 반환에 필요한 열이 없습니다: ", paste(miss, collapse = ", "))
  setcolorder(out, BE_COLUMNS)[]
}
