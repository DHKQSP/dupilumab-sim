#!/usr/bin/env Rscript
# 기존 엔진(D-010) 대비 새 엔진(Phoenix 호환, D-039)의 산출물 차이표 — 시험·개체 수준 핵심 지표 (검토 의견 2026-09-24 §1-2).
# 이전 값: results/nca_engine/legacy_snapshot/ (엔진 교체 직전 커밋의 요약 CSV), 새 값: 재실행 결과.
source("R/00_setup.R"); source_project()
L <- proj_path("results", "nca_engine", "legacy_snapshot"); out <- proj_path("results", "nca_engine")
rd_old <- function(...) { f <- file.path(L, ...); if (file.exists(f)) fread(f) else NULL }
rd_new <- function(...) { f <- proj_path("results", ...); if (file.exists(f)) fread(f) else NULL }
rows <- list()
add <- function(area, item, old, new, unit = "") rows[[length(rows) + 1]] <<- data.table(area = area, item = item, old = old, new = new, diff = new - old, unit = unit)
# 1) 개인 수준(B0, 20,000명): individual_<variant>.csv
for (v in c("base", "struct2020")) {
  o <- rd_old("individual", sprintf("individual_%s.csv", v)); n <- rd_new("individual", sprintf("individual_%s.csv", v))
  if (is.null(o) || is.null(n)) next
  o <- o[schedule == "B0"]; n <- n[schedule == "B0"]
  for (m in c("reliable_pct", "lambda_ok_pct", "extrap_median", "extrap_gt20_pct", "tlast_median", "AUClast_geo", "extrap_true_median")) add(sprintf("개인 B0 (%s)", v), m, o[[m]], n[[m]])
}
# 2) 채혈 일정 판정(D3 대 B0)
for (v in c("base", "struct2020", "vmax080_both", "vmax125_both")) {
  o <- rd_old("trials", sprintf("schedule_decision_%s.csv", v)); n <- rd_new("trials", sprintf("schedule_decision_%s.csv", v))
  if (is.null(o) || is.null(n)) next
  for (sh in c("D3", "Bminus")) { oo <- o[schedule == sh]; nn <- n[schedule == sh]; if (!nrow(oo) || !nrow(nn)) next
    add(sprintf("일정 판정 %s (%s)", sh, v), "(a) CI 폭 상대 감소 %", 100 * oo$a_mean_width_rel_decrease, 100 * nn$a_mean_width_rel_decrease)
    add(sprintf("일정 판정 %s (%s)", sh, v), "(c) 신뢰 충족 증가 %p", oo$c_reliable_gain_pp, nn$c_reliable_gain_pp)
    add(sprintf("일정 판정 %s (%s)", sh, v), "(d) 외삽 20% 초과 비", oo$d_extrap20_ratio, nn$d_extrap20_ratio)
    add(sprintf("일정 판정 %s (%s)", sh, v), "권고(1 = 추가 채혈)", as.numeric(any(oo$recommend)), as.numeric(any(nn$recommend))) }
}
# 3) 제품 시나리오(500회): 통과율
for (v in c("", "_struct2020")) {
  o <- rd_old("trials", sprintf("products_per_endpoint%s.csv", v)); n <- rd_new("trials", sprintf("products_per_endpoint%s.csv", v))
  if (is.null(o) || is.null(n)) next
  for (sc in c("S00", "F090", "VM125", "KE120")) for (ep in c("AUClast", "AUCinf_reliable", "AUCinf_all", "Cmax")) {
    oo <- o[scenario == sc & endpoint == ep]; nn <- n[scenario == sc & endpoint == ep]; if (!nrow(oo) || !nrow(nn)) next
    add(sprintf("제품 500회 (%s)", if (v == "") "base" else "struct2020"), sprintf("%s %s 통과율 %%", sc, ep), oo$pass_rate, nn$pass_rate)
  }
}
# 4) 5,000회 이상 인용 비율
o <- rd_old("trials5000", "products5000_props_base.csv"); n <- rd_new("trials5000", "products5000_props_base.csv")
if (!is.null(o) && !is.null(n)) {
  m <- merge(o[, .(scenario, metric, old = est)], n[, .(scenario, metric, new = est)], by = c("scenario", "metric"))
  for (i in seq_len(nrow(m))) if (m$metric[i] %in% c("AUCinf_reliable", "joint_3rel", "discord_last_pass_inf_fail", "AUClast")) add("제품 5,000회 이상 (base)", paste(m$scenario[i], m$metric[i], "%"), m$old[i], m$new[i])
}
# 5) 기준 (d) 20만 명
for (v in c("base", "struct2020", "vmax080_both")) {
  o <- rd_old("individual200k", sprintf("criterion_d_200k_%s.csv", v)); n <- rd_new("individual200k", sprintf("criterion_d_200k_%s.csv", v))
  if (!is.null(o) && !is.null(n)) { add(sprintf("기준 (d) 20만 명 (%s)", v), "D3 외삽 20% 초과 비", o$extrap_gt20_ratio, n$extrap_gt20_ratio); add(sprintf("기준 (d) 20만 명 (%s)", v), "D3 신뢰 충족 증가 %p", o$reliable_gain_pp, n$reliable_gain_pp) }
}
# 6) 단계 1 gate
for (m_ in c("step1", "step1_k2020")) {
  o <- rd_old(m_, "step1b_quant_gate.csv"); n <- rd_new(m_, "step1b_quant_gate.csv"); if (is.null(o) || is.null(n)) next
  mm <- merge(o[gate_role == "gate", .(id, old = AUClast_ratio)], n[, .(id, new = AUClast_ratio)], by = "id")
  add(sprintf("단계 1 gate (%s)", m_), "AUClast 모의/관측 비 최대 절대 차이", 0, max(abs(mm$new - mm$old)))
}
res <- rbindlist(rows)
fwrite(res, file.path(out, "engine_difference_outputs.csv")); print(res, digits = 4)
