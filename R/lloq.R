# lloq.R — 연구 LLOQ 민감도 (추가 지시 2026-09-26 §2, config/assay.yaml, config/prereg_20260926.yaml section2, D-054).
# 같은 대상자·같은 채혈 시각·같은 잔차(y_raw)를 LLOQ마다 다시 검열한다(쌍대). 잔차 변형:
#   fixed  = 잔차 모형 그대로(가산 SD는 모델 개발 자료, LLOQ p$lloq에서 추정한 값)
#   scaled = 가산 잔차를 LLOQ / p$lloq 배(LLOQ에서의 상대 정밀도 유지). 잔차는 같은 시드로 다시 뽑고 배율 1에서 저장 y_raw 재현을 검사(rescale_additive).

# 개인 수준: run_individual_population 결과(pop: subj, obs(y_raw 포함), truth, lloq, grid)에서 일정 하나의 NCA를 LLOQ × 잔차 변형마다.
# 반환 list(summary = summarize_individual + lloq, resid; nca = 대상자별(id, lloq, resid, 지표))
lloq_individual <- function(pop, design, p, lloqs, resid, master_seed, tag, schedule = "B0") {
  stopifnot(all(resid %in% c("fixed", "scaled")), !anyDuplicated(lloqs), all(lloqs > 0))
  sd_days <- get_schedule(design, schedule)
  ob0 <- subset_schedule(pop$obs, sd_days)
  eps <- if ("scaled" %in% resid) with_seed(derive_seed(master_seed, tag, "eps"), draw_eps(pop$subj$id, sort(unique(c(0, pop$grid))), p$sigma)) else NULL
  summ <- list(); ncas <- list()
  for (rv_ in resid) for (L in lloqs) {
    ob <- ob0
    if (rv_ == "scaled") ob <- rescale_additive(ob, eps, L / p$lloq)
    if (rv_ == "scaled" || L != pop$lloq) ob <- recensor_obs(ob, L)
    nca <- attach_truth(run_nca(ob), ob, pop$truth)
    nca <- merge(nca, pop$subj[, intersect(c("id", "WT", "sex", "ada", "HT", "BMI"), names(pop$subj)), with = FALSE], by = "id")
    k <- paste(rv_, L)
    summ[[k]] <- summarize_individual(nca, last_planned = max(sd_days))[, `:=`(lloq = L, resid = rv_)]
    tp <- ob[, .(id, tlast = time, tlast_planned = planned)]           # tlast(실제 시각)의 계획 채혈일
    nn <- tp[nca, on = c("id", "tlast")]
    ncas[[k]] <- nn[, .(id, lloq = L, resid = rv_, tlast, tlast_planned, lambda_ok, reliable, rel_i = (lambda_ok & !(flag_rsq %in% TRUE) & !(flag_extrap %in% TRUE)) %in% TRUE,
                        pct_extrap, pct_extrap_true, coverage_true, AUClast, Cmax, n_quant)]
  }
  list(summary = rbindlist(summ), nca = rbindlist(ncas))
}
