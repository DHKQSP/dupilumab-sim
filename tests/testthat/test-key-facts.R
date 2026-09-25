# 보고 우선순위와 핵심 문장(지시 2026-09-25 §4): key_facts의 값·전제 검사, 연장 문구, 영문 표의 문자 규칙
kf_copy <- function(root, rel) { for (r in rel) { f <- proj_path("results", r); dir.create(file.path(root, dirname(r)), recursive = TRUE, showWarnings = FALSE); file.copy(f, file.path(root, r), overwrite = TRUE) } }
KF_FILES <- c("oc/inversion_all.csv", "oc/inversion_scan_k2016_Km.csv", "oc/inversion_scan_k2020_Km.csv", "oc/boundary_type1.csv",
              "trials5000/products5000_props_base.csv", "fallback/consumer_risk.csv", "fallback/discordance_classification.csv", "rationale/pillar2_products_B0.csv")

test_that("key_facts: 결과가 없으면 원소가 NULL", {
  kf <- key_facts(file.path(tempdir(), "no_such_results"))
  expect_null(kf$km); expect_null(kf$inv_bnd); expect_null(kf$vm150); expect_null(kf$ke120); expect_null(kf$p2max)
  expect_equal(nrow(key_summary_table(kf, "en")), 0)
})

test_that("key_facts: 저장 결과에서 Km 범위·도달 목표, 역산 경계표, VM150·KE120", {
  skip_if_not(all(file.exists(proj_path("results", KF_FILES))), "OC·예비 제품 결과 없음")
  root <- file.path(tempdir(), "kf_ok"); unlink(root, recursive = TRUE); kf_copy(root, KF_FILES)
  kf <- key_facts(root)
  expect_true(kf$km$ratio[1] > 0.8 && kf$km$ratio[2] < 1.25)
  expect_length(kf$km$target, 1); expect_named(kf$km$mult, c("k2016", "k2020"))
  expect_equal(kf$km$range_mult, as.numeric(unlist(read_cfg("oc_design.yaml")$mechanisms$Km$range)))
  expect_equal(nrow(kf$inv_bnd), uniqueN(kf$inv_bnd[, .(model, mechanism, target)]))
  expect_true(all(kf$inv_bnd[reachable == FALSE, abs(log(end_auc_ratio)) < abs(log(target)) & sign(log(end_auc_ratio)) == sign(log(target))]))
  expect_true(kf$vm150$rel$est > 5 && kf$vm150$true_ratio < 0.8)
  expect_true(kf$ke120$true_ratio >= 0.8 && kf$ke120$true_ratio <= 1.25)
  en <- c(unlist(key_summary_table(kf, "en")), unlist(key_inv_bnd_table(kf$inv_bnd, "en")), key_km_sentence(kf$km, "en"))
  expect_false(any(grepl("[가-힣]|—|–|−", en)))                      # 영문: 한글·em dash·en dash·U+2212 없음
})

test_that("key_facts: 문구의 전제와 어긋나면 중단", {
  skip_if_not(all(file.exists(proj_path("results", KF_FILES))), "OC·예비 제품 결과 없음")
  root <- file.path(tempdir(), "kf_bad"); unlink(root, recursive = TRUE)
  kf_copy(root, setdiff(KF_FILES, "rationale/pillar2_products_B0.csv"))
  # VM150 참값이 범위 안이면 "참값이 한계 밖인데 5% 초과" 문구의 전제가 깨진다
  cr <- fread(file.path(root, "fallback", "consumer_risk.csv")); cr[scenario == "VM150", true_ratio := 0.85]; fwrite(cr, file.path(root, "fallback", "consumer_risk.csv"))
  expect_error(key_facts(root), "VM150 참 AUC0-inf 비 < 0.80", fixed = TRUE)
  kf_copy(root, "fallback/consumer_risk.csv")
  # VM150 AUCinf 신뢰군 통과율이 5% 이하면 "명목 5% 초과" 문구의 전제가 깨진다
  pr <- fread(file.path(root, "trials5000", "products5000_props_base.csv")); pr[scenario == "VM150" & metric == "AUCinf_reliable", est := 4.9]
  fwrite(pr, file.path(root, "trials5000", "products5000_props_base.csv"))
  expect_error(key_facts(root), "명목 5%", fixed = TRUE)
  kf_copy(root, "trials5000/products5000_props_base.csv")
  # KE120 참값이 범위 밖이면 "AUCinf 위음성" 분류의 전제가 깨진다
  dc <- fread(file.path(root, "fallback", "discordance_classification.csv")); dc[scenario == "KE120", true_ratio := 0.75]
  fwrite(dc, file.path(root, "fallback", "discordance_classification.csv"))
  expect_error(key_facts(root), "KE120 참 AUC0-inf 비가 80–125% 안", fixed = TRUE)
  kf_copy(root, "fallback/discordance_classification.csv")
  # Km으로 닿는 목표가 둘이면 "도달 가능한 목표는 하나" 문구의 전제가 깨진다
  iv <- fread(file.path(root, "oc", "inversion_all.csv")); i <- iv[, which(model == "k2016" & mechanism == "Km" & direction == "up" & abs(target - 1.11) < 1e-9)]
  iv[i, `:=`(reachable = TRUE, multiplier = 99, auc_ratio = 1.11, cmax_ratio = 1.01, within_tol = TRUE)]; fwrite(iv, file.path(root, "oc", "inversion_all.csv"))
  expect_error(key_facts(root), "도달 가능한 사전 고정 목표는 하나", fixed = TRUE)
  kf_copy(root, "oc/inversion_all.csv")
  expect_silent(key_facts(root))
})

test_that("key_facts: 연장 문구의 시험 수는 boundary_type1.csv가 연장 파일을 모두 반영했을 때만(아니면 중단)", {
  skip_if_not(all(file.exists(proj_path("results", KF_FILES))), "OC·예비 제품 결과 없음")
  root <- file.path(tempdir(), "kf_ext"); unlink(root, recursive = TRUE); kf_copy(root, KF_FILES)
  bt <- fread(file.path(root, "oc", "boundary_type1.csv")); r <- bt[config == "P2"][which.max(pass_pct)]
  fwrite(data.table(model = r$model, scenario = r$scenario, config = "P2", n_before = r$n_trials, pass_pct = r$pass_pct, lo = r$lo, hi = r$hi, threshold = 5,
                    rule = "test", selected = TRUE, n_after = 2L * r$n_trials), file.path(root, "oc", "extension_decision.csv"))
  kf <- key_facts(root)                                              # 판정만 있고 연장 파일 없음: "미실행"
  expect_equal(kf$p2max$ext$n_after, 2L * r$n_trials)
  expect_match(key_trials_text(kf$p2max$row, kf$p2max$ext, "en"), "not yet run", fixed = TRUE)
  ext <- data.table(trial = r$n_trials + rep(1:2, each = 2), scenario = r$scenario, endpoint = rep(c("Cmax", "AUClast"), 2), GMR = 0.8, CI_lower = 0.7, CI_upper = 0.9,
                    pass = FALSE, n_R = 117L, n_T = 117L)
  fwrite(ext, file.path(root, "oc", sprintf("oc_trials_ext_be_%s.csv.gz", r$model)))
  expect_error(key_facts(root), "scripts/33이 연장 파일을 반영", fixed = TRUE)   # 연장 2회가 있는데 표는 사전 고정 수: "미실행" 문구가 틀린다
  bt[model == r$model & scenario == r$scenario, n_trials := n_trials + 2L]; fwrite(bt, file.path(root, "oc", "boundary_type1.csv"))
  kf <- key_facts(root)
  expect_match(key_trials_text(kf$p2max$row, kf$p2max$ext, "en"), "adaptive extension in progress", fixed = TRUE)
})

test_that("key_trials_text: 연장 없음, 미실행, 진행 중, 완료", {
  r <- data.table(model = "k2020", scenario = "V2_up_080", n_trials = 10000L)
  e <- data.table(n_before = 10000L, n_after = 20000L, pass_pct = 5.32, lo = 4.9, hi = 5.78)
  expect_equal(key_trials_text(r, NULL, "en"), "10,000 trials")
  expect_equal(key_trials_text(r, NULL, "ko"), "시험 10,000회")
  expect_match(key_trials_text(r, e, "en"), "selected for adaptive extension to 20,000 trials, not yet run", fixed = TRUE)
  expect_match(key_trials_text(copy(r)[, n_trials := 12000L], e, "en"), "12,000 trials (adaptive extension in progress, target 20,000; pre-registered 10,000 trials: 5.32% (95% CI 4.90 to 5.78))", fixed = TRUE)
  expect_equal(key_trials_text(copy(r)[, n_trials := 20000L], e, "ko"), "시험 20,000회(적응적 연장; 사전 고정 10,000회 5.32% [4.90, 5.78])")
})

test_that("random_space_ranges_text: config 범위 문구", {
  oc <- list(random_space = list(ranges = list(F = c(0.8, 1.25), Km = c(0.2, 5))))
  expect_equal(random_space_ranges_text(oc, " to "), "F 0.80 to 1.25, Km 0.20 to 5")
  expect_equal(random_space_ranges_text(oc), "F 0.80–1.25, Km 0.20–5")
})
