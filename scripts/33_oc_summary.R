#!/usr/bin/env Rscript
# §3–§4 요약: 운용특성 곡선, 경계 1종 오류, 검정력, 무작위 제품 공간의 소비자·생산자 위험, 세 구성 비교, 그림 3-A~3-D, 자동 문구.
# 입력: results/oc/ (30·31·32 산출물), config/oc_design.yaml. 사전 고정 커밋 해시는 git에서 읽는다.
source("R/00_setup.R"); source_project()
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
oc <- read_cfg("oc_design.yaml"); out_dir <- if (length(args)) args[1] else proj_path("results", "oc"); fig_dir <- out_dir   # 인수: 시험용 디렉터리
models <- intersect(c("k2016", "k2020"), unique(sub("^oc_trials_be_(.*)\\.csv\\.gz$", "\\1", list.files(out_dir, pattern = "^oc_trials_be_.*\\.csv\\.gz$"))))
MODEL_LABEL <- c(k2016 = "Kovalenko 2016 (주)", k2020 = "Kovalenko 2020 Model 1")
CFG <- c("P2", "F3A", "F3B", "F3C", "G2", "AUClast_only", "AUCinf_only"); CFG_LABEL <- c(P2 = "P2", F3A = "F3-A", F3B = "F3-B", F3C = "F3-C", G2 = "G2", AUClast_only = "AUClast 단독", AUCinf_only = "AUCinf 단독")
MECH <- c("F", "ka", "ke", "Vmax", "Km", "V2")
bnd <- as.numeric(unlist(oc$boundary_targets))

# 사전 고정 해시: config/oc_design.yaml을 처음 추가한 커밋, 그 뒤 변경 여부
git <- function(...) tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE), error = function(e) character(0))
prereg_hash <- tail(git("log", "--diff-filter=A", "--format=%H", "--", "config/oc_design.yaml"), 1)
prereg_changed <- length(prereg_hash) && length(git("diff", "--name-only", prereg_hash, "--", "config/oc_design.yaml")) > 0
fwrite(data.table(prereg_commit = if (length(prereg_hash)) prereg_hash else NA_character_, changed_since = prereg_changed), file.path(out_dir, "prereg.csv"))

inv_all <- rbindlist(lapply(list.files(out_dir, pattern = "^inversion_k20(16|20)_.*\\.csv$", full.names = TRUE), fread))
fwrite(inv_all[order(model, mechanism, direction, target)], file.path(out_dir, "inversion_all.csv"))

cfg_pass <- function(be) { w <- config_pass(be, oc); w }
tab_all <- list(); bnd_all <- list(); pow_all <- list(); comp_all <- list(); gmr_all <- list()
for (mdl in models) {
  be <- fread(file.path(out_dir, sprintf("oc_trials_be_%s.csv.gz", mdl)))
  sc <- fread(file.path(out_dir, sprintf("oc_scenarios_%s.csv", mdl)))
  w <- cfg_pass(be)
  long <- melt(w, id.vars = c("trial", "scenario"), measure.vars = paste0("cfg_", CFG), variable.name = "config", value.name = "pass")
  long[, config := sub("^cfg_", "", config)]
  tab <- long[, { wc <- wilson_ci(sum(pass), .N); .(n_trials = .N, pass_pct = wc$est, lo = wc$lo, hi = wc$hi) }, by = .(scenario, config)]
  info <- rbind(sc[, .(scenario = code, mechanism, direction, target, multiplier, auc_ratio, cmax_ratio)],
                data.table(scenario = "S00", mechanism = "동일", direction = "-", target = 1, multiplier = 1, auc_ratio = 1, cmax_ratio = 1))
  tab <- merge(tab, info, by = "scenario"); tab[, model := mdl]
  tab_all[[mdl]] <- tab
  # (b) 경계 1종 오류
  bt <- tab[abs(target - bnd[1]) < 1e-9 | abs(target - bnd[2]) < 1e-9]
  bnd_all[[mdl]] <- bt
  # (c) 검정력: 1.00(동일 제품), 0.95, 1.05
  pow_all[[mdl]] <- tab[abs(target - 1) < 1e-9 | abs(target - 0.95) < 1e-9 | abs(target - 1.05) < 1e-9]
  # 4절 구성 비교: 같은 시험의 쌍대 차이 (P2 − F3x: 추가 보호/추가 탈락, G2 − P2)
  pr <- function(a, b) paired_prop_diff_ci(a, b)
  comp <- w[, {
    d1 <- pr(cfg_P2, cfg_F3A); d2 <- pr(cfg_P2, cfg_F3C); d3 <- pr(cfg_G2, cfg_P2); d4 <- pr(cfg_P2, cfg_F3B)
    .(n_trials = .N, P2_minus_F3A = d1$est, P2_minus_F3A_lo = d1$lo, P2_minus_F3A_hi = d1$hi, P2_minus_F3B = d4$est, P2_minus_F3B_lo = d4$lo, P2_minus_F3B_hi = d4$hi,
      P2_minus_F3C = d2$est, P2_minus_F3C_lo = d2$lo, P2_minus_F3C_hi = d2$hi, G2_minus_P2 = d3$est, G2_minus_P2_lo = d3$lo, G2_minus_P2_hi = d3$hi,
      P2_pass_AUCinfA_fail = 100 * mean(cfg_P2 & !(AUCinf_A %in% TRUE)))
  }, by = scenario]
  comp <- merge(comp, info, by = "scenario"); comp[, `:=`(model = mdl, region = fifelse(abs(auc_ratio - bnd[1]) < 0.002 | abs(auc_ratio - bnd[2]) < 0.002, "경계",
                                                                               fifelse(auc_ratio < bnd[1] | auc_ratio > bnd[2], "범위 밖", "범위 안")))]
  comp_all[[mdl]] <- comp
  # 평가변수별 시험 GMR의 기하평균(시험 간) 대 참값: 비구획 AUCinf의 치우침(규칙 A = 추정 + 선택, 규칙 B = 추정)을 분해
  g <- be[, .(gmr_geo = exp(mean(log(GMR), na.rm = TRUE)), n_R_mean = mean(n_R), n_T_mean = mean(n_T), n_trials = .N), by = .(scenario, endpoint)]
  g <- merge(g, info, by = "scenario")
  g[, truth := fifelse(endpoint == "Cmax", cmax_ratio, auc_ratio)]
  g[, rel_bias_pct := 100 * (gmr_geo / truth - 1)][, model := mdl]
  gmr_all[[mdl]] <- g
}
tab_all <- rbindlist(tab_all); bnd_all <- rbindlist(bnd_all); pow_all <- rbindlist(pow_all); comp_all <- rbindlist(comp_all)
gmr_all <- rbindlist(gmr_all); if (nrow(gmr_all)) fwrite(gmr_all[order(model, mechanism, direction, target, endpoint)], file.path(out_dir, "gmr_by_endpoint.csv"))
fwrite(tab_all, file.path(out_dir, "oc_curves.csv")); fwrite(bnd_all, file.path(out_dir, "boundary_type1.csv"))
fwrite(pow_all, file.path(out_dir, "power.csv")); fwrite(comp_all, file.path(out_dir, "config_comparison.csv"))

# (d) 무작위 제품 공간
rs_models <- intersect(models, unique(sub("^random_trials_be_(.*)\\.csv\\.gz$", "\\1", list.files(out_dir, pattern = "^random_trials_be_.*\\.csv\\.gz$"))))
risk_all <- list(); bins_all <- list(); smooth_all <- list(); rs_prod_all <- list()
nb <- oc$random_space$near_boundary
SPLINE_DF <- as.integer(sub(".*df = ([0-9]+).*", "\\1", oc$random_space$smoothing)); stopifnot(!is.na(SPLINE_DF))   # 사전 고정 문구의 df
for (mdl in rs_models) {
  tr <- fread(file.path(out_dir, sprintf("random_truth_%s.csv", mdl)))
  be <- fread(file.path(out_dir, sprintf("random_trials_be_%s.csv.gz", mdl)))
  w <- config_pass(be, oc); w[, product := trial]
  d <- merge(w, tr, by = "product")
  d[, inside := true_auc_ratio >= bnd[1] & true_auc_ratio <= bnd[2]]
  d[, inside_both := inside & true_cmax_ratio >= bnd[1] & true_cmax_ratio <= bnd[2]]
  d[, near_out := (true_auc_ratio >= nb$low[[1]] & true_auc_ratio < bnd[1]) | (true_auc_ratio > bnd[2] & true_auc_ratio <= nb$high[[2]])]
  d[, near_in := (true_auc_ratio >= bnd[1] & true_auc_ratio <= nb$low[[2]]) | (true_auc_ratio >= nb$high[[1]] & true_auc_ratio <= bnd[2])]
  rs_prod_all[[mdl]] <- d[, .(model = mdl, product, true_auc_ratio, true_cmax_ratio, auc_se_log, inside, inside_both)]
  risk_all[[mdl]] <- rbindlist(lapply(CFG, function(cf) { pcol <- paste0("cfg_", cf)
    f <- function(sel, fail = FALSE) { x <- d[sel][[pcol]]; if (fail) x <- !x; w_ <- wilson_ci(sum(x), length(x)); c(w_$est, w_$lo, w_$hi, length(x)) }
    rbind(data.table(model = mdl, config = cf, truth = "AUC0-inf", metric = "소비자 위험(범위 밖 통과)", scope = "전체", t(f(!d$inside))),
          data.table(model = mdl, config = cf, truth = "AUC0-inf", metric = "소비자 위험(범위 밖 통과)", scope = "경계 근처", t(f(d$near_out))),
          data.table(model = mdl, config = cf, truth = "AUC0-inf", metric = "생산자 위험(범위 안 불통과)", scope = "전체", t(f(d$inside, TRUE))),
          data.table(model = mdl, config = cf, truth = "AUC0-inf", metric = "생산자 위험(범위 안 불통과)", scope = "경계 근처", t(f(d$near_in, TRUE))),
          data.table(model = mdl, config = cf, truth = "AUC0-inf·Cmax 모두", metric = "소비자 위험(범위 밖 통과)", scope = "전체", t(f(!d$inside_both))),
          data.table(model = mdl, config = cf, truth = "AUC0-inf·Cmax 모두", metric = "생산자 위험(범위 안 불통과)", scope = "전체", t(f(d$inside_both, TRUE)))) }))
  setnames(risk_all[[mdl]], c("V1", "V2", "V3", "V4"), c("pct", "lo", "hi", "n_products"))
  bw <- as.numeric(oc$random_space$bins_log_width)
  d[, bin := round(floor(log(true_auc_ratio) / bw) * bw + bw / 2, 6)]
  bins_all[[mdl]] <- rbindlist(lapply(CFG, function(cf) d[, { x <- get(paste0("cfg_", cf)); wc <- wilson_ci(sum(x), .N); .(n = .N, pass_pct = wc$est, lo = wc$lo, hi = wc$hi) }, by = bin][, `:=`(config = cf, model = mdl, ratio = exp(bin))]))
  grid <- data.table(lr = seq(min(log(d$true_auc_ratio)), max(log(d$true_auc_ratio)), length.out = 200))   # 자료 범위 안에서만(외삽 없음)
  smooth_all[[mdl]] <- rbindlist(lapply(CFG, function(cf) {
    fit <- glm(as.formula(sprintf("cfg_%s ~ splines::ns(log(true_auc_ratio), df = %d)", cf, SPLINE_DF)), family = binomial(), data = d)
    pr <- predict(fit, newdata = data.table(true_auc_ratio = exp(grid$lr)), type = "response")
    data.table(model = mdl, config = cf, ratio = exp(grid$lr), pass_pct = 100 * pr) }))
}
if (length(rs_models)) {
  fwrite(rbindlist(risk_all), file.path(out_dir, "random_space_risks.csv")); fwrite(rbindlist(bins_all), file.path(out_dir, "random_space_bins.csv"))
  fwrite(rbindlist(smooth_all), file.path(out_dir, "random_space_smooth.csv"))
  rs_prod <- rbindlist(rs_prod_all)
  fwrite(rs_prod[, .(n = .N, inside_pct = 100 * mean(inside), lt080_pct = 100 * mean(true_auc_ratio < bnd[1]), gt125_pct = 100 * mean(true_auc_ratio > bnd[2]),
                     near_low_pct = 100 * mean(true_auc_ratio >= 0.75 & true_auc_ratio <= 0.85), near_high_pct = 100 * mean(true_auc_ratio >= 1.18 & true_auc_ratio <= 1.33),
                     truth_se_log_median = median(auc_se_log), truth_se_log_p95 = quantile(auc_se_log, 0.95)), by = model], file.path(out_dir, "random_space_truth_distribution.csv"))
}

# ----- 그림 3-A~3-D: 한국어(보고서)와 영문(summary_en.md, 파일명 _en) 두 벌. 색 + 선 모양 + 점 모양으로 구성 식별 ------------------------------
FX <- list(
  ko = list(model = MODEL_LABEL, cfg = CFG_LABEL, cap_n = "시험 반복: 경계·동일 제품 %s회, 나머지 %s회, arm당 117명, B0, 60–90 kg", xA = "참 AUC0-inf 비 (로그 척도, 200,000명 공통 난수)", yA = "통과 확률 (%)",
            tA = "그림 3-A. 기전별 운용특성 곡선 — %s", sA = ".\n세로선 0.80·1.25, 가로선 5%", same = "동일", panel = "참값 %.2f · %s %s (×%.3g)", yB = "경계 1종 오류 (%, Wilson 95% 구간)",
            tB = "그림 3-B. 경계 1종 오류 — %s", sB = "기전 × 구성, 경계 시나리오 각 %s회. 점선 = 5%%. 막대 아래 구성 이름으로 식별", yC1 = "제품 수", tC = "그림 3-C. 무작위 제품 공간 — %s",
            sC = "라틴 하이퍼큐브 %s개 제품(기전별 로그 균등), 제품당 시험 1회, 제품별 참값 1,000명 공통 난수.\n위: 참 AUC0-inf 비 분포, 아래: 구간별 통과율(점, n ≥ 20)과 로지스틱 평활(선)",
            xC = "참 AUC0-inf 비 (로그)", yC = "통과율 (%)", xD = "목표 참 AUC0-inf 비 (로그)", yD = "필요한 시험군 배율 (로그)", tD = "그림 3-D. 목표 참값 비에 필요한 기전별 배율",
            sD = "×: 탐색 범위 끝에서도 도달 불가(표시 위치 = 범위 끝 배율). 200,000명 공통 난수, 이분법 ±0.1%"),
  en = list(model = c(k2016 = "Kovalenko 2016 (primary)", k2020 = "Kovalenko 2020 Model 1"), cfg = c(CFG_LABEL[1:5], AUClast_only = "AUClast only", AUCinf_only = "AUCinf only"),
            cap_n = "Trials: boundary and identical-product scenarios %s each, others %s each; 117 per arm, B0, 60 to 90 kg", xA = "True AUC0-inf ratio (log scale, 200,000 CRN subjects)", yA = "Pass probability (%)",
            tA = "Figure 3-A. Operating characteristic curves by mechanism, %s", sA = ".\nVertical lines 0.80 and 1.25, horizontal line 5%", same = "identical", panel = "True %.2f · %s %s (x%.3g)",
            yB = "Boundary type I error (%, Wilson 95% CI)", tB = "Figure 3-B. Boundary type I error, %s", sB = "Mechanism by configuration, %s trials per boundary scenario. Dashed line = 5%%. Configurations named under the bars",
            yC1 = "Products", tC = "Figure 3-C. Random product space, %s",
            sC = "Latin hypercube, %s products (log-uniform per mechanism), one trial each, truth per product from 1,000 CRN subjects.\nTop: true AUC0-inf ratio; bottom: pass rate per bin (points, n >= 20) and logistic smooth (lines)",
            xC = "True AUC0-inf ratio (log)", yC = "Pass rate (%)", xD = "Target true AUC0-inf ratio (log)", yD = "Required test-arm multiplier (log)", tD = "Figure 3-D. Multiplier required for each target true ratio, by mechanism",
            sD = "x: not reachable even at the end of the search range (plotted at the range-end multiplier).\n200,000 CRN subjects, bisection to within 0.1%"))
cols4 <- c(P2 = VIZ$s1, F3A = VIZ$s2, F3C = VIZ$s3, G2 = "#eda100"); lt4 <- c(P2 = "solid", F3A = "22", F3C = "42", G2 = "12"); sh4 <- c(P2 = 16, F3A = 17, F3C = 15, G2 = 18)
cfg5 <- c("P2", "F3A", "F3B", "F3C", "G2"); fill5 <- c(VIZ$s1, VIZ$s2, "#e87ba4", VIZ$s3, "#eda100")
for (lg in names(FX)) {
  T_ <- FX[[lg]]; sfx <- if (lg == "en") "_en" else ""
  cL <- setNames(cols4, T_$cfg[names(cols4)]); lL <- setNames(lt4, T_$cfg[names(lt4)]); sL <- setNames(sh4, T_$cfg[names(sh4)])
  for (mdl in models) {
    x <- tab_all[model == mdl & config %in% names(cols4)]
    isb <- x$scenario == "S00" | abs(x$target - bnd[1]) < 1e-9 | abs(x$target - bnd[2]) < 1e-9      # 캡션의 반복 수는 실제 시험 수
    cap_n <- sprintf(T_$cap_n, format(max(x$n_trials[isb]), big.mark = ","), format(if (any(!isb)) max(x$n_trials[!isb]) else 0, big.mark = ","))
    s0 <- x[scenario == "S00"]
    xx <- rbind(x[mechanism %in% MECH], rbindlist(lapply(MECH, function(mc) copy(s0)[, mechanism := mc])))
    xx[, config := factor(config, levels = names(cols4), labels = T_$cfg[names(cols4)])]
    xx[, mechanism := factor(mechanism, levels = MECH)]
    setorder(xx, mechanism, config, auc_ratio)
    g <- ggplot(xx, aes(auc_ratio, pass_pct, colour = config, linetype = config, shape = config)) +
      geom_vline(xintercept = bnd, colour = VIZ$muted, linewidth = 0.4) + geom_hline(yintercept = 5, colour = VIZ$muted, linewidth = 0.4, linetype = "22") +
      geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
      scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.11, 1.25, 1.43)) + scale_colour_manual(values = cL, name = NULL) + scale_linetype_manual(values = lL, name = NULL) +
      scale_shape_manual(values = sL, name = NULL) + facet_wrap(~mechanism, ncol = 3) +
      labs(x = T_$xA, y = T_$yA, title = sprintf(T_$tA, T_$model[[mdl]]), subtitle = paste0(cap_n, T_$sA)) + theme_dupi()
    ggsave(file.path(fig_dir, sprintf("fig3A_oc_curves_%s%s.png", mdl, sfx)), g, width = 11, height = 7, dpi = 120)
    b <- bnd_all[model == mdl & config %in% cfg5]
    if (nrow(b)) {
      b[, config := factor(config, levels = cfg5, labels = T_$cfg[cfg5])]
      # 패널 = 참값 경계 × 기전(도달 방향 표시). 빈 패널 없이 참값 0.80 줄, 1.25 줄
      b[, panel := sprintf(T_$panel, target, mechanism, fifelse(direction == "down", "↓", "↑"), multiplier)]
      b[, panel := factor(panel, levels = unique(b[order(target, match(mechanism, MECH), direction), panel]))]
      g <- ggplot(b, aes(config, pass_pct, fill = config)) + geom_col(width = 0.7) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = VIZ$ink2) +
        geom_hline(yintercept = 5, colour = VIZ$ink, linetype = "22") + geom_text(aes(label = sprintf("%.1f", pass_pct), y = hi), vjust = -0.4, size = 2.5, colour = VIZ$ink2) +
        scale_fill_manual(values = fill5, guide = "none") + facet_wrap(~panel, ncol = max(b[, uniqueN(panel), by = target]$V1)) +   # 한 줄 = 한 경계(0.80이 더 많을 때) scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(x = NULL, y = T_$yB, title = sprintf(T_$tB, T_$model[[mdl]]), subtitle = sprintf(T_$sB, format(max(b$n_trials), big.mark = ","))) +
        theme_dupi() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      ggsave(file.path(fig_dir, sprintf("fig3B_boundary_type1_%s%s.png", mdl, sfx)), g, width = 12, height = 6, dpi = 120)
    }
  }
  if (length(rs_models)) for (mdl in rs_models) {
    c4n <- c("P2", "F3A", "F3C", "G2")
    bb <- rbindlist(bins_all)[model == mdl & config %in% c4n & n >= 20]
    sm <- rbindlist(smooth_all)[model == mdl & config %in% c4n]
    lab <- T_$cfg[c4n]; c4 <- setNames(unname(cols4[c4n]), lab); l4 <- setNames(unname(lt4[c4n]), lab); s4 <- setNames(unname(sh4[c4n]), lab)
    bb[, config := factor(config, levels = c4n, labels = lab)]; sm[, config := factor(config, levels = c4n, labels = lab)]
    pd <- rbindlist(rs_prod_all)[model == mdl]
    xl <- range(pd$true_auc_ratio) * c(0.98, 1.02)   # 모든 제품이 보이도록 자료 범위로
    g1 <- ggplot(pd, aes(true_auc_ratio)) + geom_histogram(bins = 80, fill = VIZ$s1) + geom_vline(xintercept = bnd, colour = VIZ$ink, linetype = "22") +
      scale_x_log10(limits = xl) + labs(x = NULL, y = T_$yC1, title = sprintf(T_$tC, T_$model[[mdl]]), subtitle = sprintf(T_$sC, format(nrow(pd), big.mark = ","))) + theme_dupi()
    g2 <- ggplot() + geom_vline(xintercept = bnd, colour = VIZ$muted) + geom_hline(yintercept = 5, colour = VIZ$muted, linetype = "22") +
      geom_point(data = bb, aes(ratio, pass_pct, colour = config, shape = config), size = 1.6, alpha = 0.8) +
      geom_line(data = sm, aes(ratio, pass_pct, colour = config, linetype = config), linewidth = 0.7) +
      scale_x_log10(limits = xl) + scale_colour_manual(values = c4, name = NULL) + scale_linetype_manual(values = l4, name = NULL) + scale_shape_manual(values = s4, name = NULL) +
      labs(x = T_$xC, y = T_$yC) + theme_dupi()
    png(file.path(fig_dir, sprintf("fig3C_random_space_%s%s.png", mdl, sfx)), width = 10 * 120, height = 8 * 120, res = 120)
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 1)))
    print(g1, vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1)); print(g2, vp = grid::viewport(layout.pos.row = 3:5, layout.pos.col = 1)); dev.off()
  }
  if (nrow(inv_all)) {
    iv <- copy(inv_all)
    iv[, mechanism := factor(mechanism, levels = MECH)]; iv[, model_l := factor(T_$model[model], levels = unname(T_$model))]
    g <- ggplot(iv[reachable == TRUE], aes(target, multiplier, colour = model_l, shape = model_l)) +
      geom_hline(yintercept = 1, colour = VIZ$muted) + geom_vline(xintercept = bnd, colour = VIZ$muted, linetype = "22") +
      geom_line(aes(group = interaction(model_l, direction)), linewidth = 0.5) + geom_point(size = 2) +
      geom_point(data = iv[reachable == FALSE], aes(target, end_multiplier), shape = 4, size = 2.4, stroke = 0.9) +
      scale_x_log10(breaks = c(0.7, 0.8, 0.9, 1, 1.11, 1.25, 1.43)) + scale_y_log10(labels = function(v) format(v, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)) +
      scale_colour_manual(values = c(VIZ$s1, VIZ$s2), name = NULL, drop = FALSE) + scale_shape_manual(values = c(16, 17), name = NULL, drop = FALSE) + facet_wrap(~mechanism, ncol = 3, scales = "free_y") +
      labs(x = T_$xD, y = T_$yD, title = T_$tD, subtitle = T_$sD) + theme_dupi()
    ggsave(file.path(fig_dir, sprintf("fig3D_inversion_multipliers%s.png", sfx)), g, width = 11, height = 7, dpi = 120)
  }
}
cat("summary written\n"); print(bnd_all[config %in% c("P2", "F3A", "F3C", "G2"), .(model, mechanism, direction, target, config, pass_pct, lo, hi)], digits = 3)

# ----- 자동 결론 문구(한국어·영문). 고정 문구의 전제는 결과로 검사하고 어긋나면 중단한다(§4) ---------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)
MECH_KO <- c(F = "F(흡수량)", ka = "ka(흡수 속도)", ke = "ke(선형 소실)", Vmax = "Vmax(표적 매개 소실)", Km = "Km(결합)", V2 = "V2(분포)")
MECH_EN <- c(F = "bioavailability (F)", ka = "absorption rate (ka)", ke = "linear elimination (ke)", Vmax = "target-mediated elimination capacity (Vmax)", Km = "binding constant (Km)", V2 = "peripheral volume (V2)")
ko <- c(); en <- c()
if (nrow(bnd_all)) {
  p2 <- bnd_all[config == "P2"]
  ex <- p2[pass_pct > 5][order(-pass_pct)]
  if (nrow(ex)) {
    stopifnot(all(ex$pass_pct > 5))
    ko <- c(ko, sprintf("**P2(AUClast + Cmax)의 경계 1종 오류가 5%%를 넘는 경우가 있다.** %s. 해당 기전·크기·구간은 3절 표와 그림 3-B에 있다.",
                        paste(sprintf("%s %s, 참값 %.2f(배율 %.3g): %s%% [%s, %s]%s", MODEL_LABEL[ex$model], MECH_KO[ex$mechanism], ex$target, ex$multiplier, f1(ex$pass_pct), f1(ex$lo), f1(ex$hi),
                                      fifelse(ex$lo > 5, " — 구간 하한도 5% 초과", "")), collapse = "; ")))
    en <- c(en, sprintf("**The boundary type I error of P2 (AUClast + Cmax) exceeds 5%% in some cases.** %s.",
                        paste(sprintf("%s, %s, true ratio %.2f (multiplier %.3g): %s%% (95%% CI %s to %s)%s", c(k2016 = "2016 model", k2020 = "Model 1")[ex$model], MECH_EN[ex$mechanism], ex$target, ex$multiplier,
                                      f1(ex$pass_pct), f1(ex$lo), f1(ex$hi), fifelse(ex$lo > 5, ", lower bound also above 5%", "")), collapse = "; ")))
  } else {
    mx <- p2[which.max(pass_pct)]
    stopifnot(max(p2$pass_pct) <= 5)
    ko <- c(ko, sprintf("P2(AUClast + Cmax)의 경계 1종 오류는 모든 기전·방향·모델에서 5%% 이하다(%d개 경계 시나리오, 각 %s회). 최대 %s%% [%s, %s]: %s %s, 참값 %.2f.",
                        nrow(p2), format(min(p2$n_trials), big.mark = ","), f1(mx$pass_pct), f1(mx$lo), f1(mx$hi), MODEL_LABEL[[mx$model]], MECH_KO[[mx$mechanism]], mx$target))
    en <- c(en, sprintf("The boundary type I error of P2 (AUClast + Cmax) is at most 5%% for every mechanism, direction and model (%d boundary scenarios, %s trials each); the largest is %s%% (95%% CI %s to %s) for %s, %s, true ratio %.2f.",
                        nrow(p2), format(min(p2$n_trials), big.mark = ","), f1(mx$pass_pct), f1(mx$lo), f1(mx$hi), c(k2016 = "the 2016 model", k2020 = "Model 1")[[mx$model]], MECH_EN[[mx$mechanism]], mx$target))
  }
  g2 <- bnd_all[config == "G2"]; exg <- g2[pass_pct > 5]
  ko <- c(ko, sprintf("G2(AUCinf 규칙 A + Cmax)의 경계 1종 오류: 범위 %s–%s%%%s.", f1(min(g2$pass_pct)), f1(max(g2$pass_pct)),
                      if (nrow(exg)) sprintf(", 5%% 초과 %d건(%s)", nrow(exg), paste(sprintf("%s %s %.2f: %s%%", c(k2016 = "2016", k2020 = "Model 1")[exg$model], exg$mechanism, exg$target, f1(exg$pass_pct)), collapse = "; ")) else ", 모두 5% 이하"))
  en <- c(en, sprintf("G2 (AUCinf rule A + Cmax) boundary type I error ranges from %s%% to %s%%%s.", f1(min(g2$pass_pct)), f1(max(g2$pass_pct)),
                      if (nrow(exg)) sprintf("; above 5%% in %d cases (%s)", nrow(exg), paste(sprintf("%s %s %.2f: %s%%", c(k2016 = "2016", k2020 = "Model 1")[exg$model], exg$mechanism, exg$target, f1(exg$pass_pct)), collapse = "; ")) else "; all at most 5%"))
  # 원인: G2가 가장 높은 경계 시나리오(모델별)의 평가변수별 평균 GMR 대 참값. 전제(비구획 AUCinf가 AUClast보다 참값에서 멀다)를 검사해 문구를 고른다
  for (mdl in unique(g2$model)) {
    x <- g2[model == mdl][which.max(pass_pct)]; gg <- gmr_all[model == mdl & scenario == x$scenario]
    gv <- setNames(gg$gmr_geo, gg$endpoint); tr <- x$auc_ratio
    d_last <- abs(log(gv[["AUClast"]] / tr)); d_A <- abs(log(gv[["AUCinf_A"]] / tr)); d_B <- abs(log(gv[["AUCinf_B"]] / tr))
    toward1 <- abs(log(gv[["AUCinf_A"]])) < abs(log(tr))
    ko <- c(ko, sprintf("%s G2 최대 시나리오(%s %s, 참 AUC0-inf 비 %s): 시험 GMR 기하평균은 AUCinf 규칙 A %s, 규칙 B %s, AUClast %s, Cmax %s(참값 %s). %s",
                        MODEL_LABEL[[mdl]], MECH_KO[[x$mechanism]], x$direction, formatC(tr, format = "f", digits = 3), formatC(gv[["AUCinf_A"]], format = "f", digits = 3), formatC(gv[["AUCinf_B"]], format = "f", digits = 3),
                        formatC(gv[["AUClast"]], format = "f", digits = 3), formatC(gv[["Cmax"]], format = "f", digits = 3), formatC(x$cmax_ratio, format = "f", digits = 3),
                        if (d_A > d_last && toward1) sprintf("비구획 AUCinf가 참값보다 1 쪽으로 치우쳐(규칙 A %s%%, 규칙 B %s%%; AUClast %s%%) 범위 밖 제품이 통과하기 쉽다.", formatC(100 * (gv[["AUCinf_A"]] / tr - 1), format = "f", digits = 1, flag = "+"),
                                                            formatC(100 * (gv[["AUCinf_B"]] / tr - 1), format = "f", digits = 1, flag = "+"), formatC(100 * (gv[["AUClast"]] / tr - 1), format = "f", digits = 1, flag = "+"))
                        else "비구획 AUCinf의 치우침이 AUClast보다 크지 않다(차이는 표본 수·변동 때문)."))
    en <- c(en, sprintf("%s, scenario with the largest G2 error (%s, %s, true AUC0-inf ratio %s): geometric mean of trial GMRs AUCinf rule A %s, rule B %s, AUClast %s, Cmax %s (true %s). %s",
                        c(k2016 = "2016 model", k2020 = "Model 1")[[mdl]], MECH_EN[[x$mechanism]], x$direction, formatC(tr, format = "f", digits = 3), formatC(gv[["AUCinf_A"]], format = "f", digits = 3), formatC(gv[["AUCinf_B"]], format = "f", digits = 3),
                        formatC(gv[["AUClast"]], format = "f", digits = 3), formatC(gv[["Cmax"]], format = "f", digits = 3), formatC(x$cmax_ratio, format = "f", digits = 3),
                        if (d_A > d_last && toward1) sprintf("Non-compartmental AUCinf is biased toward 1 relative to the truth (rule A %s%%, rule B %s%%; AUClast %s%%), so products outside the limits pass more often.", formatC(100 * (gv[["AUCinf_A"]] / tr - 1), format = "f", digits = 1, flag = "+"),
                                                            formatC(100 * (gv[["AUCinf_B"]] / tr - 1), format = "f", digits = 1, flag = "+"), formatC(100 * (gv[["AUClast"]] / tr - 1), format = "f", digits = 1, flag = "+"))
                        else "The bias of non-compartmental AUCinf is not larger than that of AUClast (differences reflect sample size and variability)."))
  }
}
if (nrow(comp_all)) {
  cb <- comp_all[region == "경계"]; ci_ <- comp_all[region == "범위 안" & mechanism != "동일"]; co <- comp_all[region == "범위 밖"]
  rng <- function(x, col) sprintf("%s–%s%%p", f2(min(x[[col]])), f2(max(x[[col]])))
  rng_en <- function(x, col) sprintf("%s to %s percentage points", f2(min(x[[col]])), f2(max(x[[col]])))
  ko <- c(ko, sprintf("F3 대 P2: 경계 시나리오에서 F3-A가 추가로 막는 비율(P2 통과·F3-A 불통과, 같은 시험 쌍대) %s, 범위 밖 %s; 범위 안(참값 0.85–1.18)에서 추가 탈락 %s. G2 대 P2: 경계에서 G2 − P2 통과율 차이 %s.",
                      rng(cb, "P2_minus_F3A"), if (nrow(co)) rng(co, "P2_minus_F3A") else "해당 없음", rng(ci_, "P2_minus_F3A"), rng(cb, "G2_minus_P2")))
  en <- c(en, sprintf("F3 versus P2: at the boundaries the additional protection of F3-A (P2 passes, F3-A fails; paired within trials) is %s, outside the limits %s; inside the limits (true ratio 0.85 to 1.18) the additional failure is %s. G2 versus P2 at the boundaries: pass-rate difference G2 minus P2 %s.",
                      rng_en(cb, "P2_minus_F3A"), if (nrow(co)) rng_en(co, "P2_minus_F3A") else "not applicable",
                      rng_en(ci_, "P2_minus_F3A"), rng_en(cb, "G2_minus_P2")))
}
if (nrow(pow_all)) {
  s0 <- pow_all[scenario == "S00" & config %in% c("P2", "F3A", "F3C", "G2")]
  ko <- c(ko, sprintf("검정력(참값 1.00, 동일 제품, 시험 %s회): %s.", format(max(s0$n_trials), big.mark = ","), paste(sprintf("%s %s: %s%% [%s, %s]", c(k2016 = "2016", k2020 = "Model 1")[s0$model], CFG_LABEL[s0$config], f1(s0$pass_pct), f1(s0$lo), f1(s0$hi)), collapse = "; ")))
  en <- c(en, sprintf("Power at true ratio 1.00 (identical products, %s trials): %s.", format(max(s0$n_trials), big.mark = ","), paste(sprintf("%s %s: %s%% (%s to %s)", c(k2016 = "2016", k2020 = "Model 1")[s0$model], sub("단독", "only", CFG_LABEL[s0$config]), f1(s0$pass_pct), f1(s0$lo), f1(s0$hi)), collapse = "; ")))
}
if (length(rs_models)) {
  rk <- rbindlist(risk_all)[config %in% c("P2", "F3A", "F3C", "G2")]
  lst <- function(x, d = 2, en = FALSE) paste(x[, sprintf(if (en) "%s %s%% (%s to %s)" else "%s %s%% [%s, %s]", CFG_LABEL[config], formatC(pct, format = "f", digits = d), formatC(lo, format = "f", digits = d), formatC(hi, format = "f", digits = d))], collapse = ", ")
  for (mdl in rs_models) { x <- rk[model == mdl]; pd <- rbindlist(rs_prod_all)[model == mdl]
    cA <- x[truth == "AUC0-inf" & startsWith(metric, "소비자") & scope == "전체"]; cN <- x[truth == "AUC0-inf" & startsWith(metric, "소비자") & scope == "경계 근처"]
    pA <- x[truth == "AUC0-inf" & startsWith(metric, "생산자") & scope == "전체"]; pN <- x[truth == "AUC0-inf" & startsWith(metric, "생산자") & scope == "경계 근처"]
    pB <- x[truth != "AUC0-inf" & startsWith(metric, "생산자")]; cB <- x[truth != "AUC0-inf" & startsWith(metric, "소비자")]
    ko <- c(ko, sprintf("무작위 제품 공간(%s, %s개 제품, 제품당 시험 1회; 참 AUC0-inf 비 범위 안 %s%%, 참 AUC0-inf·Cmax 모두 범위 안 %s%%): 소비자 위험(참 AUC0-inf 범위 밖 제품 중 통과, %s개) %s; 경계 근처(0.75–0.80, 1.25–1.33, %s개) %s. 생산자 위험(참 AUC0-inf 범위 안 제품 중 불통과, %s개) %s; 경계 근처(0.80–0.85, 1.18–1.25, %s개) %s. 참 AUC0-inf·Cmax 모두 기준: 소비자 위험 %s, 생산자 위험 %s.",
      MODEL_LABEL[[mdl]], format(nrow(pd), big.mark = ","), f1(100 * mean(pd$inside)), f1(100 * mean(pd$inside_both)),
      format(cA$n_products[1], big.mark = ","), lst(cA), format(cN$n_products[1], big.mark = ","), lst(cN), format(pA$n_products[1], big.mark = ","), lst(pA, 1), format(pN$n_products[1], big.mark = ","), lst(pN, 1), lst(cB), lst(pB, 1)))
    en <- c(en, sprintf("Random product space (%s, %s products, one trial each; true AUC0-inf ratio inside the limits for %s%%, true AUC0-inf and Cmax both inside for %s%%): consumer risk (pass among %s products with true AUC0-inf outside) %s; near the boundaries (0.75 to 0.80 and 1.25 to 1.33, %s products) %s. Producer risk (failure among %s products with true AUC0-inf inside) %s; near the boundaries (0.80 to 0.85 and 1.18 to 1.25, %s products) %s. Against the joint truth (AUC0-inf and Cmax): consumer risk %s, producer risk %s.",
      c(k2016 = "2016 model", k2020 = "Model 1")[[mdl]], format(nrow(pd), big.mark = ","), f1(100 * mean(pd$inside)), f1(100 * mean(pd$inside_both)),
      format(cA$n_products[1], big.mark = ","), lst(cA, en = TRUE), format(cN$n_products[1], big.mark = ","), lst(cN, en = TRUE), format(pA$n_products[1], big.mark = ","), lst(pA, 1, TRUE),
      format(pN$n_products[1], big.mark = ","), lst(pN, 1, TRUE), lst(cB, en = TRUE), lst(pB, 1, TRUE))) }
}
writeLines(ko, file.path(out_dir, "oc_conclusion_ko.md")); writeLines(en, file.path(out_dir, "oc_conclusion_en.md"))
cat(ko, sep = "\n\n")
