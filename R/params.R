# params.R — config 읽기와 파라미터 객체 (지시서 2026-09-23 §2)
# load_params(model, variant) 반환:
#   $theta   : 명명 벡터 (Vc, ke, k12, k21, ka, Vmax, Km, F [, MTT, n_transit])
#   $cov     : list(theta_WT, WT_ref)
#   $omega   : 명명 벡터, 로그 척도 SD (= sqrt(ω²))
#   $omega2  : 명명 벡터, 분산(원자료)
#   $sigma   : c(prop, add)
#   $lloq    : 정량한계
#   $model_id: "k2016_2cmt_linMM_ka" | "k2020_model1_transit"
#   $ada     : list(fraction, onset_day, ke_multiplier) — ADA 유사 소집단(기본 fraction 0)
#   $status  : data.table(item, status, source)

read_cfg <- function(name) yaml::read_yaml(proj_path("config", name), fileEncoding = "UTF-8", readLines.warn = FALSE)

.num <- function(x) if (is.list(x) && !is.null(x$value)) as.numeric(x$value) else as.numeric(x)

load_params <- function(model = c("k2016", "k2020"), variant = "base") {
  model <- match.arg(model)
  typ <- if (model == "k2016") read_cfg("params_typical.yaml") else read_cfg("params_k2020_model1.yaml")
  var <- read_cfg("params_variability.yaml")
  st <- list()
  th_names <- c("Vc", "ke", "k12", "k21", "ka", "Vmax", "Km", "F", if (model == "k2020") c("MTT", "n_transit"))
  theta <- setNames(numeric(length(th_names)), th_names)
  for (nm in th_names) {
    e <- typ$theta[[nm]]; theta[nm] <- .num(e)
    st[[length(st) + 1]] <- data.table(item = paste0("theta.", nm), status = e$status, source = e$source)
  }
  cov <- list(theta_WT = .num(typ$covariates$WT_on_Vc$theta_WT), WT_ref = .num(typ$covariates$WT_on_Vc$WT_ref), ke_bmi_exp = 0, BMI_ref = 26)
  st[[length(st) + 1]] <- data.table(item = "cov.theta_WT", status = typ$covariates$WT_on_Vc$theta_WT$status, source = typ$covariates$WT_on_Vc$theta_WT$source)
  st[[length(st) + 1]] <- data.table(item = "cov.WT_ref", status = typ$covariates$WT_on_Vc$WT_ref$status, source = typ$covariates$WT_on_Vc$WT_ref$source)
  # k2016: params_variability.yaml(Kovalenko 2016). k2020: 자체 IIV·잔차(Supplementary Table 2, D-029), MTT IIV는 마지막(k2016 난수 순서 불변)
  om_names <- c("Vc", "ke", "k12", "k21", "ka", "Vmax", "Km", "F", if (model == "k2020") "MTT")
  omega2 <- setNames(numeric(length(om_names)), om_names)
  if (model == "k2020") var <- list(iiv_omega2 = lapply(typ$iiv$sd, function(e) list(omega2 = e$omega2, status = "confirmed", source = typ$iiv$source)), residual = typ$residual)
  for (nm in om_names) {
    e <- var$iiv_omega2[[nm]]; omega2[nm] <- as.numeric(e$omega2)
    st[[length(st) + 1]] <- data.table(item = paste0("omega2.", nm), status = e$status, source = e$source)
  }
  sigma <- c(prop = as.numeric(var$residual$sigma_prop$value), add = as.numeric(var$residual$sigma_add$value))
  st[[length(st) + 1]] <- data.table(item = "sigma.prop", status = var$residual$sigma_prop$status, source = var$residual$sigma_prop$source)
  st[[length(st) + 1]] <- data.table(item = "sigma.add", status = var$residual$sigma_add$status, source = var$residual$sigma_add$source)
  lloq <- if (!is.null(typ$lloq_mg_L)) .num(typ$lloq_mg_L) else 0.078
  p <- list(theta = theta, cov = cov, omega2 = omega2, omega = sqrt(omega2), sigma = sigma, lloq = lloq,
            model_id = typ$model_id, model = model, variant = variant,
            ada = list(fraction = 0, onset_day = 14, ke_multiplier = 1),
            status = rbindlist(st))
  class(p) <- c("dupi_params", "list")
  p
}

# 시나리오 배율(시험군) 적용 — 고정효과에 곱한다
apply_multipliers <- function(p, mult) {
  if (length(mult) == 0) return(p)
  for (nm in names(mult)) { stopifnot(nm %in% names(p$theta)); p$theta[nm] <- p$theta[nm] * as.numeric(mult[[nm]]) }
  p
}

# 민감도 변형(지시서 §2, config/scenarios.yaml sensitivity_variants)
apply_variant <- function(p, variant_spec, design = NULL) {
  if (!is.null(variant_spec$theta_multipliers)) p <- apply_multipliers(p, variant_spec$theta_multipliers)   # 양 군 공통(곡률 민감도)
  if (!is.null(variant_spec$omega2_multiplier)) { p$omega2 <- p$omega2 * variant_spec$omega2_multiplier; p$omega <- sqrt(p$omega2) }
  if (!is.null(variant_spec$sigma_prop)) p$sigma["prop"] <- variant_spec$sigma_prop
  if (!is.null(variant_spec$sigma_add)) p$sigma["add"] <- variant_spec$sigma_add
  if (!is.null(variant_spec$ke_bmi)) { p$cov$ke_bmi_exp <- as.numeric(variant_spec$ke_bmi$exponent); p$cov$BMI_ref <- as.numeric(variant_spec$ke_bmi$bmi_ref) }   # §6 (c)(d)
  if (!is.null(variant_spec$theta_WT_override)) p$cov$theta_WT <- as.numeric(variant_spec$theta_WT_override)                                                   # §6 (d)
  if (isTRUE(variant_spec$ada) && !is.null(design)) p$ada <- list(fraction = design$ada_sensitivity$fraction, onset_day = design$ada_sensitivity$onset_day, ke_multiplier = design$ada_sensitivity$ke_multiplier)
  p
}

# 확정 여부 검사: 모델 파라미터에 pending이 남아 있으면 중단
assert_params_confirmed <- function(p) {
  bad <- p$status[!status %in% c("confirmed", "assumption")]
  if (nrow(bad)) stop("출처 미확정 파라미터:\n", paste(sprintf("  %s (%s)", bad$item, bad$status), collapse = "\n"))
  invisible(TRUE)
}

print.dupi_params <- function(x, ...) {
  cat("dupilumab params [", x$model_id, " / variant=", x$variant, "]\n", sep = "")
  print(round(x$theta, 4)); cat("omega2:\n"); print(x$omega2); cat("sigma:\n"); print(x$sigma)
  invisible(x)
}

cv_to_omega <- function(cv_pct) sqrt(log(1 + (cv_pct / 100)^2))
omega_to_cv <- function(omega) 100 * sqrt(exp(omega^2) - 1)
log_cv_pct <- function(x) { x <- x[is.finite(x) & x > 0]; 100 * sqrt(exp(sd(log(x))^2) - 1) }
geo_mean <- function(x) { x <- x[is.finite(x) & x > 0]; exp(mean(log(x))) }
