# model.R — rxode2 모델 (SPEC §3). 기본: Kovalenko 2016 2구획 + 선형/MM 병렬 소실 + 1차 흡수.
# 구조 민감도: Kovalenko 2020 Model 1 — Savic transit: 투여 구획 → transit 3개 → 흡수 구획(모두 ktr) → 중심(ka), ktr = (n+1)/MTT. ADA 유사 소집단: onset 이후 ke × 배율(지시서 §5 선택 민감도).
# 상태(mg): depot [, tr1..tr3], central, periph; auc = ∫C dt (참값 진단 전용). 농도 C = central/Vc_i.
# 개체별 파라미터는 R에서 계산해 params로 넘긴다. rxode2 내부 난수는 쓰지 않는다. F는 파서 충돌을 피해 Fbio로 넘긴다.

MODEL_CODE <- list(
  k2016_2cmt_linMM_ka = "
    C      = central / Vc_i
    keEff  = ke
    if (ada > 0.5 && t >= t_ada) { keEff = ke * ada_mult }
    MMrate = Vmax * C / (Km + C) * Vc_i
    f(depot) = Fbio
    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot - (keEff + k12) * central + k21 * periph - MMrate
    d/dt(periph)  =  k12 * central - k21 * periph
    d/dt(auc)     =  C
  ",
  k2020_model1_transit = "
    C      = central / Vc_i
    keEff  = ke
    if (ada > 0.5 && t >= t_ada) { keEff = ke * ada_mult }
    MMrate = Vmax * C / (Km + C) * Vc_i
    f(depot) = Fbio
    d/dt(depot)   = -ktr * depot
    d/dt(tr1)     =  ktr * depot - ktr * tr1
    d/dt(tr2)     =  ktr * tr1   - ktr * tr2
    d/dt(tr3)     =  ktr * tr2   - ktr * tr3
    d/dt(absc)    =  ktr * tr3   - ka * absc
    d/dt(central) =  ka * absc - (keEff + k12) * central + k21 * periph - MMrate
    d/dt(periph)  =  k12 * central - k21 * periph
    d/dt(auc)     =  C
  ")


.model_env <- new.env(parent = emptyenv())
get_model <- function(model_id = "k2016_2cmt_linMM_ka") {
  if (is.null(.model_env[[model_id]])) {
    if (!requireNamespace("rxode2", quietly = TRUE)) stop("rxode2가 설치되어 있지 않습니다.")
    wd <- file.path(tempdir(), "rxode2_models"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)
    .model_env[[model_id]] <- rxode2::rxode2(MODEL_CODE[[model_id]], modName = model_id, wd = wd)
  }
  .model_env[[model_id]]
}

individual_Vc <- function(Vc_pop, WT, theta_WT, WT_ref, eta_Vc = 0) Vc_pop * (WT / WT_ref)^theta_WT * exp(eta_Vc)

# eta: data.table(id, WT[, eta_Vc, eta_ke, eta_k12, eta_k21, eta_ka, eta_Vmax, eta_Km, eta_F, ada])
individual_params <- function(p, eta) {
  stopifnot(is.data.table(eta), "id" %in% names(eta))
  g <- function(nm) if (nm %in% names(eta)) eta[[nm]] else rep(0, nrow(eta))
  WT <- if ("WT" %in% names(eta)) eta$WT else rep(p$cov$WT_ref, nrow(eta))
  th <- p$theta
  ip <- data.table(
    id = eta$id, WT = WT,
    Vc_i = individual_Vc(th[["Vc"]], WT, p$cov$theta_WT, p$cov$WT_ref, g("eta_Vc")),
    ke = th[["ke"]] * exp(g("eta_ke")) * (if (!is.null(p$cov$ke_bmi_exp) && p$cov$ke_bmi_exp != 0) { if (!"BMI" %in% names(eta)) stop("ke~BMI 공변량에 BMI 열이 필요합니다"); (eta$BMI / p$cov$BMI_ref)^p$cov$ke_bmi_exp } else 1), k12 = th[["k12"]] * exp(g("eta_k12")), k21 = th[["k21"]] * exp(g("eta_k21")),
    Vmax = th[["Vmax"]] * exp(g("eta_Vmax")), Km = th[["Km"]] * exp(g("eta_Km")), ka = th[["ka"]] * exp(g("eta_ka")),
    F = plogis(qlogis(th[["F"]]) + g("eta_F")),
    ada = if ("ada" %in% names(eta)) as.numeric(eta$ada) else 0,
    t_ada = p$ada$onset_day, ada_mult = p$ada$ke_multiplier)
  if (p$model == "k2020") { mtt_i <- th[["MTT"]] * exp(g("eta_MTT")); ip[, ktr := (th[["n_transit"]] + 1) / mtt_i] }   # MTT IIV(D-029)
  ip[]
}

typical_subject <- function(p, id = 1L, WT = NULL) {
  if (is.null(WT)) WT <- p$cov$WT_ref
  individual_params(p, data.table(id = id, WT = WT))
}

# obs: data.table(id, time). dose_mg: 스칼라 또는 id별. 반환 data.table(id, time, C, auc, depot, central, periph)
solve_model <- function(ipar, obs, dose_mg, model_id = "k2016_2cmt_linMM_ka", atol = 1e-10, rtol = 1e-8, maxsteps = 500000L) {
  stopifnot(all(c("id", "time") %in% names(obs)))
  mod <- get_model(model_id)
  ids <- ipar$id
  if (length(dose_mg) == 1) dose_mg <- rep(dose_mg, length(ids))
  ev <- rbind(data.table(id = ids, time = 0, cmt = "depot", amt = dose_mg, evid = 1L),
              data.table(id = obs$id, time = obs$time, cmt = NA_character_, amt = NA_real_, evid = 0L))
  setorder(ev, id, time, -evid)
  pcols <- c("id", "Vc_i", "ke", "k12", "k21", "Vmax", "Km", "ka", "ada", "t_ada", "ada_mult", if ("ktr" %in% names(ipar)) "ktr")
  params <- as.data.frame(ipar[, c(pcols, "F"), with = FALSE]); names(params)[names(params) == "F"] <- "Fbio"
  sol <- rxode2::rxSolve(mod, params = params, events = as.data.frame(ev), atol = atol, rtol = rtol, maxsteps = maxsteps,
                         returnType = "data.table", cores = 1L, addDosing = FALSE)
  if (!"id" %in% names(sol)) sol[, id := ids[1]]
  sol[, .(id = as.integer(as.character(id)), time, C, auc, depot, central, periph,
          absorb_chain = depot + (if ("tr1" %in% names(sol)) tr1 + tr2 + tr3 + absc else 0))]
}
