skip_if_no_rxode2()
p <- load_params("k2016")

analytic_linear <- function(t, ka, ke, k12, k21, Vc, F, dose) {
  A <- matrix(c(-ka, 0, 0, ka, -(ke + k12), k21, 0, k12, -k21), 3, 3, byrow = TRUE); x0 <- c(F * dose, 0, 0)
  sapply(t, function(tt) (as.matrix(Matrix::expm(A * tt)) %*% x0)[2] / Vc)
}

test_that("Vmax=0이면 해석적 선형 2구획 해와 일치", {
  ip <- typical_subject(p); ip[, Vmax := 0]; t <- c(0.5, 1, 3, 7, 14, 28, 56)
  sol <- solve_model(ip, data.table(id = 1L, time = t), 300)
  expect_equal(sol$C, analytic_linear(t, ip$ka, ip$ke, ip$k12, ip$k21, ip$Vc_i, ip$F, 300), tolerance = 1e-6)
})

test_that("질량 보존, F 반영, MM 0차 극한, auc 상태", {
  ip <- typical_subject(p); ip[, `:=`(Vmax = 0, ke = 0)]
  sol <- solve_model(ip, data.table(id = 1L, time = c(1, 10, 100)), 300)
  expect_equal(sol$depot + sol$central + sol$periph, rep(ip$F * 300, 3), tolerance = 1e-8)
  ip2 <- typical_subject(p); ip2[, `:=`(F = 0.5, ka = 1e-6, Vmax = 0, ke = 0)]
  expect_equal(solve_model(ip2, data.table(id = 1L, time = 1e-6), 300)$depot, 150, tolerance = 1e-5)
  ip3 <- typical_subject(p); ip3[, `:=`(ke = 0, k12 = 0, ka = 1e3, Km = 1e-6)]
  s3 <- solve_model(ip3, data.table(id = 1L, time = c(1, 10)), 300)
  expect_equal((s3$central[1] - s3$central[2]) / 9, ip3$Vmax * ip3$Vc_i, tolerance = 1e-3)
  ip4 <- typical_subject(p); ip4[, Vmax := 0]; t <- seq(0, 56, by = 0.05)
  s4 <- solve_model(ip4, data.table(id = 1L, time = t), 300)
  expect_equal(tail(s4$auc, 1), sum(diff(t) * (head(s4$C, -1) + tail(s4$C, -1)) / 2), tolerance = 1e-4)
})

test_that("MM 소실로 용량당 AUCinf가 용량과 함께 증가; 다개체 일관성", {
  ip <- rbind(typical_subject(p, 1L), typical_subject(p, 2L), typical_subject(p, 3L))
  tr <- true_auc(ip, dose_mg = c(200, 300, 600)); expect_true(all(diff(tr$inf$AUCinf_true / c(200, 300, 600)) > 0))
  ip2 <- rbind(typical_subject(p, 1L), typical_subject(p, 2L, WT = 90)); t <- c(1, 7, 28)
  both <- solve_model(ip2, CJ(id = 1:2, time = t), 300); one <- solve_model(ip2[id == 2L], data.table(id = 2L, time = t), 300)
  expect_equal(both[id == 2L, C], one$C, tolerance = 1e-8)
})

test_that("ADA 스위치: onset 이전 동일, 이후 ke 배율만큼 소실 증가", {
  a0 <- typical_subject(p); a1 <- typical_subject(p); a1[, `:=`(ada = 1, t_ada = 14, ada_mult = 2)]
  s0 <- solve_model(a0, data.table(id = 1L, time = c(7, 13.9, 21, 28)), 300); s1 <- solve_model(a1, data.table(id = 1L, time = c(7, 13.9, 21, 28)), 300)
  expect_equal(s0[time < 14, C], s1[time < 14, C], tolerance = 1e-8); expect_true(all(s1[time > 14, C] < s0[time > 14, C]))
})

test_that("Kovalenko 2020 transit 모델: 질량 보존, 평균 통과 시간 = MTT, 대표 프로파일이 2016 모델과 같은 자릿수", {
  p2 <- load_params("k2020"); mid <- "k2020_model1_transit"
  ip0 <- typical_subject(p2); ip0[, `:=`(ke = 0, Vmax = 0)]
  s0 <- solve_model(ip0, data.table(id = 1L, time = c(0.1, 1, 10, 100)), 300, model_id = mid)
  expect_equal(s0$absorb_chain + s0$central + s0$periph, rep(0.643 * 300, 4), tolerance = 1e-8)
  ipm <- typical_subject(p2); ipm[, `:=`(ke = 0, Vmax = 0, k12 = 0, ka = 1e4)]; tt <- seq(0, 1, by = 0.005)
  sm <- solve_model(ipm, data.table(id = 1L, time = tt), 300, model_id = mid)
  cum <- sm$central / (0.643 * 300); dens <- diff(cum) / diff(tt)
  mtt <- sum(dens * (head(tt, -1) + diff(tt) / 2) * diff(tt)) / sum(dens * diff(tt))
  expect_equal(mtt, 0.105, tolerance = 0.02)
  s16 <- solve_model(typical_subject(p), data.table(id = 1L, time = c(7, 28)), 300)
  s20 <- solve_model(typical_subject(p2), data.table(id = 1L, time = c(7, 28)), 300, model_id = mid)
  expect_true(all(abs(log(s20$C / s16$C)) < 0.5))
})
