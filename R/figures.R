# figures.R — 보고서 그림(ggplot2, 정적 HTML). 팔레트: dataviz 참조 팔레트 1–3번(검증 통과, 전쌍),
# aqua는 대비 3:1 미만 → 모든 그림에 같은 값의 표를 보고서에 함께 둔다(relief rule). 텍스트는 잉크색, 계열 색은 표식에만.
VIZ <- list(s1 = "#2a78d6", s2 = "#eb6834", s3 = "#1baf7a", ink = "#0b0b0b", ink2 = "#52514e", muted = "#8a8984", grid = "#e6e5e0", surface = "#fcfcfb")

theme_dupi <- function() {
  ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = VIZ$surface, colour = NA),
    panel.grid.major = ggplot2::element_line(colour = VIZ$grid, linewidth = 0.3), panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(colour = VIZ$ink2), axis.title = ggplot2::element_text(colour = VIZ$ink2),
    plot.title = ggplot2::element_text(colour = VIZ$ink, face = "bold"), plot.subtitle = ggplot2::element_text(colour = VIZ$ink2),
    legend.position = "top", legend.text = ggplot2::element_text(colour = VIZ$ink2), strip.text = ggplot2::element_text(colour = VIZ$ink, face = "bold"))
}

# 단계 1: 코호트 중앙값 tlast 분포(용량별 소형 다중), 관측 중앙값은 세로선
fig_cohort_tlast <- function(dist_list, obs) {
  d <- rbindlist(lapply(names(dist_list), function(k) dist_list[[k]][, dose := paste0(k, " mg")]))
  o <- data.table(dose = paste0(obs$dose_mg, " mg"), obs = obs$obs_median)
  ggplot2::ggplot(d, ggplot2::aes(tlast_median, pct)) +
    ggplot2::geom_col(fill = VIZ$s1, width = 2.8) +
    ggplot2::geom_vline(data = o, ggplot2::aes(xintercept = obs), colour = VIZ$ink, linewidth = 0.6, linetype = "22") +
    ggplot2::geom_text(data = o, ggplot2::aes(x = obs, y = Inf, label = paste0("관측 ", obs, "일")), vjust = 1.4, hjust = 1.08, colour = VIZ$ink, size = 3.2) +
    ggplot2::facet_wrap(~dose, nrow = 1) + ggplot2::labs(x = "8명 코호트의 중앙값 tlast (투여 후 일)", y = "코호트 비율 (%)",
    title = "Clot 2021 재현: 관측 중앙값은 모의 분포 안에 있다", subtitle = "IIV·잔차 포함 8명 코호트 2,000회, 용량군별 실제 체중") + theme_dupi()
}

# 단계 1: 정량 gate — 데이터셋별 모의/관측 AUClast 비, 두 모델
fig_gate <- function(gate_base, gate_alt) {
  d <- rbind(gate_base[, .(id, ratio = AUClast_ratio, model = "Kovalenko 2016 (기본)")], gate_alt[, .(id, ratio = AUClast_ratio, model = "Kovalenko 2020 Model 1")])
  ord <- gate_base[order(dose_mg, AUClast_ratio), id]; d[, id := factor(id, levels = ord)]
  ggplot2::ggplot(d, ggplot2::aes(ratio, id, colour = model, shape = model)) +
    ggplot2::annotate("rect", xmin = 0.85, xmax = 1.15, ymin = -Inf, ymax = Inf, fill = VIZ$grid, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 1, colour = VIZ$muted, linewidth = 0.4) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = c(VIZ$s1, VIZ$s2), name = NULL) + ggplot2::scale_shape_manual(values = c(16, 17), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(min(0.78, min(d$ratio) - 0.02), max(1.22, max(d$ratio) + 0.02))) +
    ggplot2::labs(x = "AUClast 평균 모의/관측 (회색 띠 = ±15%)", y = NULL, title = "정량 gate: 200 mg ≈80 kg에서 기본 모델이 낮다") + theme_dupi()
}

# 단계 2: 일정별 외삽률(비구획 vs 실제) 중앙값·95백분위
fig_extrap_by_schedule <- function(ind) {
  d <- rbind(ind[, .(schedule, what = "비구획 추정", med = extrap_median, p95 = extrap_p95)], ind[, .(schedule, what = "실제(진적분)", med = extrap_true_median, p95 = extrap_true_p95)])
  d[, schedule := factor(schedule, levels = c("Bminus", "B0", "D4", "D1", "D2", "D3"))]
  ggplot2::ggplot(d, ggplot2::aes(schedule, med, colour = what, shape = what)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = med, ymax = p95), position = ggplot2::position_dodge(width = 0.5), linewidth = 0.8) +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = c(VIZ$s1, VIZ$s2), name = NULL) + ggplot2::scale_shape_manual(values = c(16, 17), name = NULL) +
    ggplot2::labs(x = "채혈 일정", y = "외삽률 (%) — 점 중앙값, 선 끝 95백분위", title = "일정별 외삽률: 비구획 추정과 실제의 괴리") + theme_dupi()
}

# 단계 2: 일정별 AUClast CI 평균 폭 상대 감소(B0 대비), 모델별 소형 다중. 점선 = 판정 임계 2%
fig_width_by_schedule <- function(dec_list) {
  d <- rbindlist(lapply(names(dec_list), function(m) dec_list[[m]][, .(schedule, model = m, v = 100 * a_mean_width_rel_decrease, lo = 100 * a_paired_lo, hi = 100 * a_paired_hi)]))
  d[, schedule := factor(schedule, levels = c("Bminus", "D4", "D1", "D2", "D3"))]
  ggplot2::ggplot(d, ggplot2::aes(v, schedule)) +
    ggplot2::geom_vline(xintercept = 0, colour = VIZ$muted, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 2, colour = VIZ$ink2, linewidth = 0.4, linetype = "22") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2, colour = VIZ$s1) +
    ggplot2::geom_point(size = 3, colour = VIZ$s1) + ggplot2::facet_wrap(~model, nrow = 1) +
    ggplot2::labs(x = "AUClast 90% CI 평균 폭의 B0 대비 상대 감소 (%) — 점선 = 권고 임계 2%", y = NULL,
                  title = "추가 채혈이 AUClast 정밀도에 주는 효과", subtitle = "동일 제품 500회, 같은 시험의 대응 비교 95% CI") + theme_dupi()
}

# 단계 3: 제품 시나리오별 GMR — AUClast vs AUCinf(신뢰) vs 진적분 AUCinf
fig_products_gmr <- function(pe) {
  d <- pe[endpoint %in% c("AUClast", "AUCinf_reliable", "AUCinf_true"), .(scenario, endpoint, GMR_mean)]
  d[, endpoint := factor(endpoint, levels = c("AUClast", "AUCinf_reliable", "AUCinf_true"), labels = c("AUClast", "AUCinf (신뢰군)", "AUCinf 진적분"))]
  ord <- d[endpoint == "AUClast"][order(GMR_mean), scenario]; d[, scenario := factor(scenario, levels = ord)]
  ggplot2::ggplot(d, ggplot2::aes(GMR_mean, scenario, colour = endpoint, shape = endpoint)) +
    ggplot2::annotate("rect", xmin = 0.8, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = VIZ$grid, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 1, colour = VIZ$muted, linewidth = 0.4) +
    ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = c(VIZ$s1, VIZ$s2, VIZ$s3), name = NULL) + ggplot2::scale_shape_manual(values = c(16, 17, 15), name = NULL) +
    ggplot2::scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.4), labels = c("0.70", "0.80", "0.90", "1.00", "1.10", "1.25", "1.40")) +
    ggplot2::coord_cartesian(xlim = c(min(0.75, min(d$GMR_mean) * 0.98), max(1.3, max(d$GMR_mean) * 1.02))) +
    ggplot2::labs(x = "GMR 평균 (시험/대조, 로그 척도; 회색 띠 = 80–125%)", y = NULL, title = "제품 차이 시나리오: 세 AUC 지표의 GMR") + theme_dupi()
}
