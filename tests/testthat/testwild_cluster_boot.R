library(testthat)
library(wildclusterboot)


test_that('wild_cluster_boot returns correct value',{

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  boot_reps <- 4
  clusterby <- 'clusterby'
  bootby <- clusterby
  x_interest <- 'X'
  boot_dist <- c(-1, 1)
  H0 <- 1
  enum <- F
  bound <- c('upper','lower')


  set.seed(42)
  test_p <- wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                              boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0, bound = bound)

  p_vals <- c(upper = 1, lower = 0)

  expect_equal(test_p, p_vals)

})
