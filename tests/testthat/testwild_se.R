###### TESTS wild_se FUNCTION ######

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
  data_wild <- wild_data(data, model, x_interest, H0)
  y_wild <- wild_y(data_wild = data_wild, bootby = bootby, boot_dist = boot_dist, boot_reps = boot_reps, enum = enum)
  X <- model.matrix(model)
  x_ind <- grep(x_interest, colnames(X))
  bread <- bread_cpp(X)
  B <- beta_cpp(X, bread, y_wild)
  E <- y_wild - X %*% B

  test_se <- wild_se(data_wild, E, X, bread, clusterby, x_ind)

  se_wild <- c(X1 = 0.08499428, X2 = 0.37472011, X3 = 0.08499428, X4 = 0.37472011)

  expect_equal(test_se, se_wild, tolerance = 0.0000001)

})
