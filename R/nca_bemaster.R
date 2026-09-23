# nca_bemaster.R — 비구획분석 어댑터. 이 저장소에서 NCA(AUClast, AUCinf, λz 등)를 계산하는 유일한 진입점.
# 규칙(λz 선택, BLQ 취급, 사다리꼴 방식)은 BEmaster가 정의한다. 여기서는 데이터 형식 변환만 한다 (SPEC §6, DECISIONS D-004).
#
# [PENDING] BEmaster 입수 후 확인할 것:
#   - 패키지/모듈 이름과 버전(커밋 SHA 또는 태그) → renv.lock 고정
#   - 개체별 NCA 함수 시그니처(입력: time, conc, dose, lloq; 출력 열 이름)
#   - BLQ 입력 규약(NA vs 0 vs 문자 "BLQ")
# 이 파일의 .bemaster_call_nca() 본문만 바꾸고, 반환 열 이름 계약(NCA_COLUMNS)은 유지한다.

NCA_COLUMNS <- c("id", "AUClast", "AUCinf", "Cmax", "tmax", "tlast", "Clast",
                 "lambda_z", "r2_adj", "n_lambda", "pct_extrap", "lambda_ok")

bemaster_available <- function() {
  nzchar(Sys.getenv("BEMASTER_PATH")) || requireNamespace("BEmaster", quietly = TRUE)
}

bemaster_version <- function() {
  if (requireNamespace("BEmaster", quietly = TRUE)) return(as.character(utils::packageVersion("BEmaster")))
  if (nzchar(Sys.getenv("BEMASTER_PATH"))) {
    p <- Sys.getenv("BEMASTER_PATH")
    sha <- tryCatch(system(sprintf("git -C %s rev-parse HEAD", shQuote(p)), intern = TRUE, ignore.stderr = TRUE), error = function(e) NA)
    return(paste0("path:", p, " sha:", sha))
  }
  NA_character_
}

stop_bemaster_missing <- function() {
  stop("BEmaster 모듈을 찾을 수 없습니다. 패키지 설치 또는 BEMASTER_PATH 환경변수 설정이 필요합니다 (SPEC Q1).", call. = FALSE)
}

# 실제 호출부 — BEmaster API 확인 후 여기만 구현한다.
.bemaster_call_nca <- function(time, conc, dose_mg, lloq, opts) {
  stop("BEmaster NCA 인터페이스 바인딩이 아직 구현되지 않았습니다. R/nca_bemaster.R::.bemaster_call_nca 를 BEmaster 함수 시그니처에 맞춰 작성하세요 (SPEC §6).", call. = FALSE)
}

# 개체별 NCA 실행. sim: simulate_observations() 출력. 반환: NCA_COLUMNS 열을 가진 data.table
run_nca_bemaster <- function(sim, dose_mg, lloq, opts = list()) {
  if (!bemaster_available()) stop_bemaster_missing()
  ids <- sort(unique(sim$id))
  if (length(dose_mg) == 1) dose_mg <- setNames(rep(dose_mg, length(ids)), ids)
  res <- lapply(ids, function(i) {
    d <- sim[id == i][order(time)]
    r <- .bemaster_call_nca(time = d$time, conc = d$conc, dose_mg = dose_mg[[as.character(i)]], lloq = lloq, opts = opts)
    r <- as.data.table(as.list(r))
    r[, id := i]
    r
  })
  out <- rbindlist(res, fill = TRUE)
  miss <- setdiff(NCA_COLUMNS, names(out))
  if (length(miss)) stop("BEmaster NCA 반환에 필요한 열이 없습니다: ", paste(miss, collapse = ", "))
  setcolorder(out, NCA_COLUMNS)[]
}
