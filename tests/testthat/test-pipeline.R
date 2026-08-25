test_that("pipeline is deterministic and internally coherent", {
  sim <- simulate_lid_data(N = 400, p = 10, rho = .5, n_groups = 2, seed = 3)
  r1 <- bias_aware_lid(sim$data, sim$group, B = 200, seed = 11)
  r2 <- bias_aware_lid(sim$data, sim$group, B = 200, seed = 11)
  expect_identical(r1$edges, r2$edges)
  expect_true(all(!r1$edges$final | r1$edges$E0))          # final subset of E0
  expect_true(all(r1$edges$wTO >= 0 & r1$edges$wTO <= 1))  # implementation range
  expect_true(all(is.finite(r1$edges$Q3_star_signed)))
  expect_true(r1$beta_star %in% c(0,.1,.25,.5,.75,1,1.25,1.5,2,3,5))
})
