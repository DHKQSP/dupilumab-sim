# 신뢰 충족률·탈락률 문구(검토 의견 W2 §1): fmt_flag_pair 형식과 reliability_facts 전제 검사
test_that("fmt_flag_pair: (i) 주값, (ii) 병기, 범위·원소별·부호·반올림", {
  expect_equal(fmt_flag_pair(86.565, 80.735), "(i) 86.6% [(ii) 80.7%]")
  expect_equal(fmt_flag_pair(c(8.165, 13.435), c(16.085, 19.265)), "(i) 8.2–13.4% [(ii) 16.1–19.3%]")
  expect_equal(fmt_flag_pair(c(8.165, 13.435), c(16.085, 19.265), sep = " to "), "(i) 8.2 to 13.4% [(ii) 16.1 to 19.3%]")
  expect_equal(fmt_flag_pair(c(2.755, 1.235), c(1.69, -1.935), 2, "%p", signed = TRUE, each = TRUE), c("(i) +2.76%p [(ii) +1.69%p]", "(i) +1.24%p [(ii) -1.94%p]"))
  expect_equal(fmt_num(-1.015, 2), "-1.02")      # CSV에서 읽은 0.005 배수는 0에서 먼 쪽으로
  expect_equal(fmt_num(-0.001, 2), "0.00")       # 음의 0 없음
})

test_that("reliability_facts: 39 결과가 없으면 NULL, 전제와 어긋나면 중단", {
  expect_null(reliability_facts(file.path(tempdir(), "no_such_results")))
  need <- proj_path("results", "reliability", c("reliability_two_flag_sets_summary.csv", "reliability_paired_vs_B0.csv", "dropout_reasons_by_schedule.csv"))
  skip_if_not(all(file.exists(need)), "scripts/39_reliability_flags.R 결과 없음")
  root <- file.path(tempdir(), "relfacts"); dir.create(file.path(root, "reliability"), recursive = TRUE, showWarnings = FALSE)
  file.copy(need, file.path(root, "reliability"), overwrite = TRUE)
  rf <- reliability_facts(root)
  expect_false(is.null(rf))
  expect_true(all(rf$dense$d3$gain_ii_pp < rf$dense$d3$gain_i_pp))
  expect_equal(rf$rng$fail_i, 100 - rev(rf$rng$rel_i))
  # (ii) 증감이 (i)보다 큰 칸을 만들면 "span 플래그가 증감을 낮춘다" 문구의 전제가 깨진다
  pr <- fread(need[2]); pr[variant == "base" & schedule == "D1", gain_ii_pp := gain_i_pp + 1]
  fwrite(pr, file.path(root, "reliability", basename(need[2])))
  expect_error(reliability_facts(root), "(ii) 증감 < (i) 증감", fixed = TRUE)
  # 기준 (c)를 충족한 칸이 있으면 "두 세트 모두 미충족" 문구의 전제가 깨진다
  pr <- fread(need[2]); pr[variant == "base" & schedule == "D3", crit_c_i := TRUE]
  fwrite(pr, file.path(root, "reliability", basename(need[2])))
  expect_error(reliability_facts(root), "기준 (c)", fixed = TRUE)
})
