# model.R — rxode2 모델 정의 (SPEC §3.1, DECISIONS D-002)
#
# 상태(mg): depot, central, periph. 농도 C = central / Vc_i (mg/L). 시간 day.
# auc 상태는 모델 기반 참값 AUC(0-t) 적분용(진단 전용, NCA 대체 아님 — SPEC §6).
# 개체별 파라미터(Vc_i, ke, k12, k21, Vmax, Km, ka, F→Fbio)는 R 쪽에서 미리 계산해 params로 전달한다.
# (모델 코드 안에서는 F 대신 Fbio를 쓴다: F/T는 R에서 FALSE/TRUE 별칭이라 파서 충돌을 피한다)
# rxode2 내부 난수는 사용하지 않는다(재현성: 시드는 R/seeds.R).

DUPI_MODEL_CODE <- "
  C           = central / Vc_i
  MMrate      = Vmax * C / (Km + C) * Vc_i
  f(depot)    = Fbio
  d/dt(depot)   = -ka * depot
  d/dt(central) =  ka * depot - (ke + k12) * central + k21 * periph - MMrate
  d/dt(periph)  =  k12 * central - k21 * periph
  d/dt(auc)     =  C
"

.dupi_model_env <- new.env(parent = emptyenv())

get_dupi_model <- function() {
  if (is.null(.dupi_model_env$model)) {
    if (!requireNamespace("rxode2", quietly = TRUE)) stop("rxode2가 설치되어 있지 않습니다.")
    wd <- file.path(tempdir(), "rxode2_models"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)
    .dupi_model_env$model <- rxode2::rxode2(DUPI_MODEL_CODE, modName = "dupi_2cmt_linMM", wd = wd)
  }
  .dupi_model_env$model
}

# 개체별 Vc: 체중 공변량 (SPEC §3.3)
individual_Vc <- function(Vc_pop, WT, theta_WT, WT_ref, eta_Vc = 0) {
  Vc_pop * (WT / WT_ref)^theta_WT * exp(eta_Vc)
}

# 개체별 파라미터 표 생성: theta(대표값) + eta(개체 편차) → data.table(id, Vc_i, ke, k12, k21, Vmax, Km, ka, F)
# eta: data.table(id, WT, eta_Vc, eta_ke, eta_k12, eta_k21, eta_F, eta_Vmax, eta_ka, eta_Km) (없는 열은 0)
individual_params <- function(p, eta) {
  stopifnot(is.data.table(eta), "id" %in% names(eta))
  g <- function(nm) if (nm %in% names(eta)) eta[[nm]] else rep(0, nrow(eta))
  WT <- if ("WT" %in% names(eta)) eta$WT else rep(p$cov$WT_ref, nrow(eta))
  th <- p$theta
  F_logit <- qlogis(th[["F"]]) + g("eta_F")      # logit-normal (SPEC §3.5)
  data.table(
    id   = eta$id,
    WT   = WT,
    Vc_i = individual_Vc(th[["Vc"]], WT, p$cov$theta_WT, p$cov$WT_ref, g("eta_Vc")),
    ke   = th[["ke"]]   * exp(g("eta_ke")),
    k12  = th[["k12"]]  * exp(g("eta_k12")),
    k21  = th[["k21"]]  * exp(g("eta_k21")),
    Vmax = th[["Vmax"]] * exp(g("eta_Vmax")),
    Km   = th[["Km"]]   * exp(g("eta_Km")),
    ka   = th[["ka"]]   * exp(g("eta_ka")),
    F    = plogis(F_logit)
  )
}

# 대표 개체(모든 eta = 0, 기준 체중) 1명
typical_subject <- function(p, id = 1L, WT = NULL) {
  if (is.null(WT)) WT <- p$cov$WT_ref
  individual_params(p, data.table(id = id, WT = WT))
}

# ODE 풀기. ipar: individual_params() 출력. obs: data.table(id, time) 관측 시각(투여 후 일, 0 포함 가능).
# dose_mg: 스칼라 또는 id별 벡터. 반환: data.table(id, time, C, auc, depot, central, periph)
solve_model <- function(ipar, obs, dose_mg, atol = 1e-10, rtol = 1e-8, maxsteps = 500000L) {
  stopifnot(all(c("id", "time") %in% names(obs)))
  mod <- get_dupi_model()
  ids <- ipar$id
  if (length(dose_mg) == 1) dose_mg <- rep(dose_mg, length(ids))
  dose_ev <- data.table(id = ids, time = 0, cmt = "depot", amt = dose_mg, evid = 1L)
  obs_ev  <- data.table(id = obs$id, time = obs$time, cmt = NA_character_, amt = NA_real_, evid = 0L)
  ev <- rbind(dose_ev, obs_ev)
  setorder(ev, id, time, -evid)
  params <- as.data.frame(ipar[, .(id, Vc_i, ke, k12, k21, Vmax, Km, ka, Fbio = F)])
  sol <- rxode2::rxSolve(mod, params = params, events = as.data.frame(ev),
                         atol = atol, rtol = rtol, maxsteps = maxsteps,
                         returnType = "data.table", cores = 1L, addDosing = FALSE)
  out <- sol[, .(id = as.integer(as.character(id)), time, C, auc, depot, central, periph)]
  out
}
