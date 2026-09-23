# params.R — config 읽기, 출처 상태 검사, dev/final 모드 (DECISIONS D-005)
#
# load_params(run_mode) 는 다음을 반환한다.
#   $theta      : 명명 수치 벡터 (Vc, ke, k12, k21, F, Vmax, Km, ka)
#   $cov        : list(theta_WT, WT_ref)
#   $omega      : 명명 수치 벡터(로그/logit 척도 SD): Vc, ke, k12, k21, F, Vmax, ka, Km
#   $sigma      : c(prop, add)
#   $rse_pct    : 명명 수치 벡터 (바깥층용) 또는 NULL
#   $covariance : 로그 척도 분산-공분산 행렬 또는 NULL
#   $status     : data.table(item, status, source) — 각 값의 출처 상태
#   $run_mode   : "dev" | "final"

read_cfg <- function(name) yaml::read_yaml(proj_path("config", name), fileEncoding = "UTF-8", readLines.warn = FALSE)

.val_or_pending <- function(entry, item, dev_entry = NULL, status_tbl, field = "value") {
  v <- entry[[field]]
  st <- if (!is.null(entry$status)) entry$status else if (is.null(v)) "pending" else "confirmed"
  src <- if (!is.null(entry$source)) entry$source else NA_character_
  if (is.null(v) || (length(v) == 1 && is.na(v))) {
    if (!is.null(dev_entry) && !is.null(dev_entry[[field]])) {
      v <- dev_entry[[field]]; st <- "dev_assumption"
      src <- paste0("[DEV-가정] ", if (!is.null(dev_entry$note)) dev_entry$note else "")
    } else {
      v <- NA_real_; st <- "pending"
    }
  }
  status_tbl[[length(status_tbl) + 1]] <- data.table(item = item, status = st, value = as.numeric(v), source = src)
  list(value = as.numeric(v), status = status_tbl)
}

load_params <- function(run_mode = c("dev", "final"), scenario_multipliers = NULL,
                        theta_override = NULL) {
  run_mode <- match.arg(run_mode)
  typ <- read_cfg("params_typical.yaml")
  var <- read_cfg("params_variability.yaml")
  dev <- if (run_mode == "dev") read_cfg("dev_assumptions.yaml") else list()
  status_tbl <- list()

  # --- theta ---
  th_names <- c("Vc", "ke", "k12", "k21", "F", "Vmax", "Km", "ka")
  theta <- setNames(numeric(length(th_names)), th_names)
  for (nm in th_names) {
    r <- .val_or_pending(typ$theta[[nm]], paste0("theta.", nm), dev$theta[[nm]], status_tbl)
    theta[nm] <- r$value; status_tbl <- r$status
  }
  # --- covariate ---
  cov <- list()
  for (nm in c("theta_WT", "WT_ref")) {
    r <- .val_or_pending(typ$covariates$WT_on_Vc[[nm]], paste0("cov.", nm), dev$covariates$WT_on_Vc[[nm]], status_tbl)
    cov[[nm]] <- r$value; status_tbl <- r$status
  }
  # --- omega ---
  om_names <- c("Vc", "ke", "k12", "k21", "F", "Vmax", "ka", "Km")
  omega <- setNames(numeric(length(om_names)), om_names)
  for (nm in om_names) {
    r <- .val_or_pending(var$iiv[[nm]], paste0("omega.", nm), dev$iiv[[nm]], status_tbl, field = "omega")
    omega[nm] <- r$value; status_tbl <- r$status
  }
  # --- sigma ---
  sig <- c(prop = NA_real_, add = NA_real_)
  r <- .val_or_pending(var$residual$sigma_prop, "sigma.prop", dev$residual$sigma_prop, status_tbl); sig["prop"] <- r$value; status_tbl <- r$status
  r <- .val_or_pending(var$residual$sigma_add,  "sigma.add",  dev$residual$sigma_add,  status_tbl); sig["add"]  <- r$value; status_tbl <- r$status
  # --- uncertainty ---
  covariance <- var$uncertainty$covariance
  if (!is.null(covariance)) covariance <- as.matrix(do.call(rbind, covariance))
  rse_names <- c("Vc", "ke", "k12", "k21", "F", "Vmax", "ka")
  rse <- setNames(numeric(length(rse_names)), rse_names)
  for (nm in rse_names) {
    dv <- if (!is.null(dev$uncertainty$rse_pct[[nm]])) list(value = dev$uncertainty$rse_pct[[nm]], note = "자리표시자 RSE") else NULL
    r <- .val_or_pending(var$uncertainty$rse_pct[[nm]], paste0("rse.", nm), dv, status_tbl)
    rse[nm] <- r$value; status_tbl <- r$status
  }
  status <- rbindlist(status_tbl)

  # --- 시나리오 승수(시험약 arm) / 직접 override ---
  if (!is.null(scenario_multipliers)) for (nm in names(scenario_multipliers)) {
    stopifnot(nm %in% names(theta)); theta[nm] <- theta[nm] * scenario_multipliers[[nm]]
  }
  if (!is.null(theta_override)) for (nm in names(theta_override)) {
    stopifnot(nm %in% names(theta)); theta[nm] <- theta_override[[nm]]
  }

  out <- list(theta = theta, cov = cov, omega = omega, sigma = sig, rse_pct = rse,
              covariance = covariance, status = status, run_mode = run_mode,
              units = typ$units, model_id = typ$model_id)
  class(out) <- c("dupi_params", "list")
  out
}

# final 모드 게이트: PENDING/DEV 값이 남아 있으면 중단. 필요한 항목만 검사(need)
assert_final_ready <- function(p, need = c("theta", "cov", "omega", "sigma", "rse")) {
  bad <- p$status[status %in% c("pending", "dev_assumption") &
                  sub("\\..*$", "", item) %in% need]
  if (nrow(bad) > 0) {
    stop("run_mode = 'final' 인데 출처 미확정 값이 있습니다:\n",
         paste(sprintf("  %-14s %-15s %s", bad$item, bad$status, bad$source), collapse = "\n"))
  }
  invisible(TRUE)
}

print.dupi_params <- function(x, ...) {
  cat("dupilumab params [", x$run_mode, "]\n", sep = "")
  print(x$status)
  invisible(x)
}

# %CV → 로그 척도 SD
cv_to_omega <- function(cv_pct) sqrt(log(1 + (cv_pct / 100)^2))
omega_to_cv <- function(omega) 100 * sqrt(exp(omega^2) - 1)
