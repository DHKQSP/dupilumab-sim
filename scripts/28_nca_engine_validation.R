#!/usr/bin/env Rscript
# NCA 엔진 검증 (검토 의견 2026-09-24 §1-2, D-039): 자체 엔진(run_nca) 대 NonCompart 0.8.4(참조) 대 PKNCA 0.12.1(보조).
# (i) NonCompart 문서 예제 Theoph, Indometh (ii) 두필루맙 모의 프로필 1,000개(2016 모델 500, Model 1 500; B0, 채혈 허용창·잔차·BLQ 포함).
# 세 엔진에 같은 전처리 자료(BLQ 규칙 적용 후)를 넣는다. 통과 기준: λz 선택 점(하한·상한·점 수) 전부 일치, 파라미터 상대 차이 ≤ 1e-6.
# Indometh는 정맥 투여 자료이나 λz 알고리즘 비교를 위해 세 엔진 모두 혈관외(Cmax 제외) 규칙으로 계산한다.
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages({ library(NonCompart); library(PKNCA) })
out_dir <- proj_path("results", "nca_engine"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
MASTER_SEED <- 20260924L
logfile <- start_run_log("nca_engine_validation", master_seed = MASTER_SEED, run_mode = "final",
                         extra = list(NonCompart = as.character(packageVersion("NonCompart")), PKNCA = as.character(packageVersion("PKNCA"))))
PARS <- c("Cmax", "Tmax", "Tlast", "Clast", "AUClast", "Lambda_z", "HL_Lambda_z", "Rsq", "Rsq_adjusted", "No_points_lambda_z",
          "Lambda_z_lower", "Lambda_z_upper", "Clast_pred", "AUCINF_obs", "AUCINF_pred", "AUC_%Extrap_obs", "AUC_%Extrap_pred", "Span_ratio")

# 자체 엔진: 이미 전처리된 자료(0 = 첫 정량 이전 BLQ)를 받는다 → conc 0은 NA(BLQ)로 되돌려 넣는다(전처리 재적용은 멱등)
own_engine <- function(d) { x <- copy(d)[, conc := fifelse(conc > 0, conc, NA_real_)]; r <- run_nca(x); r[, c("id", PARS), with = FALSE] }

nc_engine <- function(d) rbindlist(lapply(split(d, by = "id"), function(x) {
  r <- NonCompart::sNCA(x$time, x$conc, dose = 1, adm = "Extravascular", down = "Log", R2ADJ = 0, excludeDelta = 1)
  data.table(id = x$id[1], Cmax = r[["CMAX"]], Tmax = r[["TMAX"]], Tlast = r[["TLST"]], Clast = r[["CLST"]], AUClast = r[["AUCLST"]],
             Lambda_z = r[["LAMZ"]], HL_Lambda_z = r[["LAMZHL"]], Rsq = r[["R2"]], Rsq_adjusted = r[["R2ADJ"]], No_points_lambda_z = as.integer(r[["LAMZNPT"]]),
             Lambda_z_lower = r[["LAMZLL"]], Lambda_z_upper = r[["LAMZUL"]], Clast_pred = r[["CLSTP"]], AUCINF_obs = r[["AUCIFO"]], AUCINF_pred = r[["AUCIFP"]],
             `AUC_%Extrap_obs` = r[["AUCPEO"]], `AUC_%Extrap_pred` = r[["AUCPEP"]], Span_ratio = (r[["LAMZUL"]] - r[["LAMZLL"]]) / r[["LAMZHL"]])
}))

pk_engine <- function(d) {
  o <- list(auc.method = "lin up/log down", adj.r.squared.factor = 1e-4, min.hl.points = 3, allow.tmax.in.half.life = FALSE,
            conc.blq = list(first = "keep", middle = "drop", last = "drop"), first.tmax = TRUE)
  miss0 <- setdiff(unique(d$id), d[time == 0, id])                     # NonCompart·자체 엔진처럼 투여 전 (0, 0) 앵커(Indometh에 0시점 없음)
  dd <- as.data.frame(setorder(rbind(d, data.table(id = miss0, time = 0, conc = 0)), id, time))
  cobj <- PKNCA::PKNCAconc(dd, conc ~ time | id)
  dobj <- PKNCA::PKNCAdose(data.frame(id = unique(d$id), time = 0, dose = 1), dose ~ time | id, route = "extravascular")
  iv <- data.frame(start = 0, end = Inf, cmax = TRUE, tmax = TRUE, tlast = TRUE, clast.obs = TRUE, auclast = TRUE, half.life = TRUE,
                   lambda.z = TRUE, r.squared = TRUE, adj.r.squared = TRUE, lambda.z.n.points = TRUE, lambda.z.time.first = TRUE,
                   lambda.z.time.last = TRUE, clast.pred = TRUE, aucinf.obs = TRUE, aucinf.pred = TRUE, aucpext.obs = TRUE, aucpext.pred = TRUE, span.ratio = TRUE)
  res <- suppressWarnings(PKNCA::pk.nca(PKNCA::PKNCAdata(cobj, dobj, intervals = iv, options = o), verbose = FALSE))
  w <- dcast(as.data.table(as.data.frame(res))[!is.na(id), .(id, PPTESTCD, PPORRES)], id ~ PPTESTCD, value.var = "PPORRES")
  g <- function(nm) if (nm %in% names(w)) as.numeric(w[[nm]]) else NA_real_
  data.table(id = w$id, Cmax = g("cmax"), Tmax = g("tmax"), Tlast = g("tlast"), Clast = g("clast.obs"), AUClast = g("auclast"),
             Lambda_z = g("lambda.z"), HL_Lambda_z = g("half.life"), Rsq = g("r.squared"), Rsq_adjusted = g("adj.r.squared"),
             No_points_lambda_z = as.integer(g("lambda.z.n.points")), Lambda_z_lower = g("lambda.z.time.first"), Lambda_z_upper = g("lambda.z.time.last"),
             Clast_pred = g("clast.pred"), AUCINF_obs = g("aucinf.obs"), AUCINF_pred = g("aucinf.pred"),
             `AUC_%Extrap_obs` = g("aucpext.obs"), `AUC_%Extrap_pred` = g("aucpext.pred"), Span_ratio = g("span.ratio"))
}

compare <- function(a, b, label_a, label_b, dataset) {
  m <- merge(a, b, by = "id", suffixes = c(".a", ".b"))
  pts <- m[, .(id, lz_a = !is.na(Lambda_z.a), lz_b = !is.na(Lambda_z.b),
               same_points = (is.na(Lambda_z.a) & is.na(Lambda_z.b)) |
                 (!is.na(Lambda_z.a) & !is.na(Lambda_z.b) & No_points_lambda_z.a == No_points_lambda_z.b & abs(Lambda_z_lower.a - Lambda_z_lower.b) < 1e-9 & abs(Lambda_z_upper.a - Lambda_z_upper.b) < 1e-9))]
  rd <- rbindlist(lapply(setdiff(PARS, "No_points_lambda_z"), function(p) {
    x <- m[[paste0(p, ".a")]]; y <- m[[paste0(p, ".b")]]; both <- !is.na(x) & !is.na(y)
    rel <- abs(x[both] - y[both]) / pmax(abs(y[both]), 1e-300)
    data.table(parameter = p, n_both = sum(both), n_na_mismatch = sum(xor(is.na(x), is.na(y))), max_rel_diff = if (any(both)) max(rel) else NA_real_)
  }))
  list(summary = data.table(dataset = dataset, comparison = paste(label_a, "vs", label_b), n_profiles = nrow(m), lz_points_identical = sum(pts$same_points),
                            lz_points_mismatch = sum(!pts$same_points), max_rel_diff = max(rd$max_rel_diff, na.rm = TRUE),
                            na_mismatch = sum(rd$n_na_mismatch), pass = all(pts$same_points) & max(rd$max_rel_diff, na.rm = TRUE) <= 1e-6 & sum(rd$n_na_mismatch) == 0),
       params = rd[, `:=`(dataset = dataset, comparison = paste(label_a, "vs", label_b))], mismatch_ids = pts[same_points == FALSE, id])
}

run_all_engines <- function(d, dataset) {
  own <- own_engine(d); nc <- nc_engine(d); pk <- pk_engine(d)
  c1 <- compare(own, nc, "자체", "NonCompart", dataset); c2 <- compare(pk, nc, "PKNCA", "NonCompart", dataset); c3 <- compare(own, pk, "자체", "PKNCA", dataset)
  list(summary = rbind(c1$summary, c2$summary, c3$summary), params = rbind(c1$params, c2$params, c3$params),
       mismatch = list(own_nc = c1$mismatch_ids, pk_nc = c2$mismatch_ids), own = own, nc = nc, pk = pk)
}

# (i) 문서 예제
th <- as.data.table(datasets::Theoph)[, .(id = as.integer(as.character(Subject)), time = Time, conc = conc)]; setorder(th, id, time)
im <- as.data.table(datasets::Indometh)[, .(id = as.integer(as.character(Subject)), time = time, conc = conc)]; setorder(im, id, time)
r_th <- run_all_engines(th, "Theoph (12명)"); r_im <- run_all_engines(im, "Indometh (6명, 혈관외 규칙)")

# (ii) 두필루맙 모의 프로필 1,000개: B0, 채혈 허용창·잔차·BLQ. 자체 전처리(preprocess_blq) 후 세 엔진에 같은 자료
design <- read_cfg("trial_design.yaml")
dupi <- rbindlist(lapply(c(k2016 = "base", k2020 = "struct2020"), function(v) {
  rv <- resolve_variant(v, design)
  pop <- run_individual_population(500, rv$p, design, "B0", MASTER_SEED, rv$wt_spec, jitter = TRUE, model_id = rv$p$model_id, tag = paste0("ncaval_", v))
  ob <- subset_schedule(pop$obs, get_schedule(design, "B0"))
  pp <- preprocess_blq(ob[, .(id, time, conc)])
  pp[, id := id + if (v == "base") 0L else 100000L]
}))
r_du <- run_all_engines(dupi, "두필루맙 모의 1,000개 (2016 500 + Model 1 500)")

# (iii) 규칙 차이 확인용 합성 사례(통과 판정 대상 아님): 말단에 상승 구간이 있어 adj R² 최대 창의 기울기가 양수인 경우,
#       NonCompart·자체는 양의 기울기 창을 먼저 제외하고 고르며, PKNCA는 전체 창 최댓값 기준이라 λz가 산출되지 않을 수 있다.
edge <- data.table(id = 1L, time = c(0, 1, 2, 4, 8, 12, 16, 24), conc = c(0, 10, 20, 12, 7, 6, 6.6, 7.3))
r_edge <- run_all_engines(edge, "합성 사례: 말단 상승(adj R² 최대 창의 기울기 양수)")
edge_tab <- unique(rbind(copy(r_edge$own)[, engine := "자체"], copy(r_edge$nc)[, engine := "NonCompart"], copy(r_edge$pk)[, engine := "PKNCA"]))
fwrite(edge_tab, file.path(out_dir, "engine_validation_edge_case.csv"))
print(edge_tab[, .(engine, Lambda_z, No_points_lambda_z, Lambda_z_lower, Lambda_z_upper, Rsq_adjusted)])

summ <- rbind(r_th$summary, r_im$summary, r_du$summary)
pars <- rbind(r_th$params, r_im$params, r_du$params)
fwrite(summ, file.path(out_dir, "engine_validation_summary.csv")); fwrite(pars, file.path(out_dir, "engine_validation_parameters.csv"))
# 불일치 사례 상세(있으면): PKNCA 대 NonCompart
mm <- rbindlist(lapply(list(r_th, r_im, r_du), function(r) {
  ids <- r$mismatch$pk_nc; if (!length(ids)) return(NULL)
  merge(r$pk[id %in% ids, .(id, pk_n = No_points_lambda_z, pk_lo = Lambda_z_lower, pk_up = Lambda_z_upper, pk_lz = Lambda_z, pk_r2adj = Rsq_adjusted)],
        r$nc[id %in% ids, .(id, nc_n = No_points_lambda_z, nc_lo = Lambda_z_lower, nc_up = Lambda_z_upper, nc_lz = Lambda_z, nc_r2adj = Rsq_adjusted)], by = "id")
}))
fwrite(mm, file.path(out_dir, "engine_validation_pknca_mismatch.csv"))
fwrite(r_du$own[, dataset := "dupilumab"], file.path(out_dir, "engine_validation_dupilumab_own.csv"))
print(summ); print(pars[max_rel_diff > 1e-9 | n_na_mismatch > 0]); if (nrow(mm)) print(mm)
append_run_log(logfile, "done")
