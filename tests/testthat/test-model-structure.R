# 단계 1 구조 테스트 (SPEC §4.1) — BEmaster·IIV 불필요
skip_if_no_rxode2()
p <- load_params("dev")

# 해석해: 선형 2구획 + 1차 흡수 (행렬지수)
analytic_linear <- function(t, ka, ke, k12, k21, Vc, F, dose) {
  A <- matrix(c(-ka, 0, 0,
                ka, -(ke + k12), k21,
                0, k12, -k21), 3, 3, byrow = TRUE)
  x0 <- c(F * dose, 0, 0)
  sapply(t, function(tt) (as.matrix(Matrix::expm(A * tt)) %*% x0)[2] / Vc)
}

test_that("Vmax=0이면 해석적 선형 2구획 해와 일치", {
  ip <- typical_subject(p); ip[, Vmax := 0]
  t <- c(0.5, 1, 3, 7, 14, 28, 56)
  sol <- solve_model(ip, data.table(id = 1L, time = t), dose_mg = 300)
  ref <- analytic_linear(t, ip$ka, ip$ke, ip$k12, ip$k21, ip$Vc_i, ip$F, 300)
  expect_equal(sol$C, ref, tolerance = 1e-6)
})

test_that("소실이 없으면 총량 = F*Dose (질량 보존)", {
  ip <- typical_subject(p); ip[, `:=`(Vmax = 0, ke = 0)]
  sol <- solve_model(ip, data.table(id = 1L, time = c(1, 10, 100)), dose_mg = 300)
  expect_equal(sol$depot + sol$central + sol$periph, rep(ip$F * 300, 3), tolerance = 1e-8)
})

test_that("F는 depot 초기량에 반영된다", {
  ip <- typical_subject(p); ip[, `:=`(F = 0.5, ka = 1e-6, Vmax = 0, ke = 0)]
  sol <- solve_model(ip, data.table(id = 1L, time = 1e-6), dose_mg = 300)
  expect_equal(sol$depot, 150, tolerance = 1e-5)
})

test_that("Km << C 에서 MM 항은 0차 소실(Vmax*Vc mg/day)로 수렴", {
  ip <- typical_subject(p); ip[, `:=`(ke = 0, k12 = 0, ka = 1e3, Km = 1e-6)]
  sol <- solve_model(ip, data.table(id = 1L, time = c(1, 10)), dose_mg = 300)
  slope <- (sol$central[1] - sol$central[2]) / 9
  expect_equal(slope, ip$Vmax * ip$Vc_i, tolerance = 1e-3)
})

test_that("auc 상태는 농도 적분과 일치(선형 케이스, 사다리꼴 근사와 비교)", {
  ip <- typical_subject(p); ip[, Vmax := 0]
  t <- seq(0, 56, by = 0.05)
  sol <- solve_model(ip, data.table(id = 1L, time = t), dose_mg = 300)
  trap <- sum(diff(t) * (head(sol$C, -1) + tail(sol$C, -1)) / 2)
  expect_equal(tail(sol$auc, 1), trap, tolerance = 1e-4)
})

test_that("MM 소실로 용량 비례성이 깨진다: 용량당 AUCinf가 용량과 함께 증가", {
  ip <- rbind(typical_subject(p, id = 1L), typical_subject(p, id = 2L), typical_subject(p, id = 3L))
  tr <- true_auc(ip, dose_mg = c(200, 300, 600))
  per_dose <- tr$inf$AUCinf_true / c(200, 300, 600)
  expect_true(all(diff(per_dose) > 0))
})

test_that("여러 개체를 한 번에 풀어도 개별 풀이와 같다", {
  ip <- rbind(typical_subject(p, id = 1L), typical_subject(p, id = 2L, WT = 100))
  t <- c(1, 7, 28)
  both <- solve_model(ip, CJ(id = 1:2, time = t), 300)
  one2 <- solve_model(ip[id == 2L], data.table(id = 2L, time = t), 300)
  expect_equal(both[id == 2L, C], one2$C, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(both[id == 1L, C], both[id == 2L, C])))
})
