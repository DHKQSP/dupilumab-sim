# scenarios.R — 제품 차이 시나리오와 민감도 변형 (config/scenarios.yaml)
load_scenarios <- function() {
  s <- read_cfg("scenarios.yaml")
  sc <- lapply(names(s$scenarios), function(k) { x <- s$scenarios[[k]]; list(code = k, label = x$label, T_multipliers = if (length(x$T_multipliers)) x$T_multipliers else list()) })
  names(sc) <- names(s$scenarios)
  list(scenarios = sc, schedule_analysis = s$schedule_analysis_scenarios, variants = s$sensitivity_variants, km_caveat = s$km_caveat)
}

# 변형 이름 → (params, wt_spec, model_id)
resolve_variant <- function(variant_name, design, sc = load_scenarios()) {
  v <- sc$variants[[variant_name]]; if (is.null(v)) stop("알 수 없는 변형: ", variant_name)
  model <- if (!is.null(v$model) && v$model == "k2020_model1") "k2020" else "k2016"
  p <- load_params(model, variant = variant_name)
  p <- apply_variant(p, v, design)
  wt <- weight_spec_from_design(design, if (!is.null(v$weight)) v$weight else "base")
  list(p = p, wt_spec = wt, model_id = p$model_id, label = v$label)
}
