library(testthat)
library(wildclusterboot)

test_that('wild_cluster_boot returns correct value',{

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  boot_reps <- 50
  clusterby <- 'clusterby'
  x_interest <- 'X'
  boot_dist <- c(-1, 1)
  H0 <- 1

  set.seed(42)
  test_p <- t_wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.52))

  set.seed(36)
  test_p <- t_wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.48))


  set.seed(128)
  test_p <- t_wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.56))

  data <- read.csv('test_files/multiway_problems.csv')
  model <- lm(data = data, formula = Y ~ W)
  boot_reps <- 50
  clusterby <- ~ G + H
  x_interest <- 'W'
  boot_dist <- 'six_pt'
  H0 <- 1

  set.seed(42)
  test_p <- t_wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.08))

})

test_that('multiply_cluster_meat calculates correct tensor', {

  x <- matrix(c(1, 1, 1, 2, 2, 2), ncol = 2, nrow = 3)
  e <- matrix(1:6, ncol = 2, nrow = 3, byrow = TRUE)
  k <- 2
  boot_reps <- 2
  xe <- cbind(x, e)

  cluster_meat <- array(c(81, 162, 162, 324, 144, 288, 288, 576), c(2, 2, 2))

  test_cluster_meat <- suppressWarnings(multiply_cluster_meat(xe = xe, boot_reps = boot_reps, k = k))

  expect_equal(test_cluster_meat, cluster_meat, tolerance = 0.0001)

})

test_that('gen_boot_weights calculates correct matrix', {

  set.seed(42)
  boot_dist <- c(-1, 1)
  boot_reps <- 4
  boot_unique <- matrix(c(1,2))
  enum <- FALSE

  test_weight_mat <- gen_boot_weights(boot_dist = boot_dist, boot_unique = boot_unique, boot_reps = boot_reps, enum = enum)
  weight_mat <- data.frame(matrix(c(1, 1, -1, 1, 1, 1, 1, -1), ncol = 4))

  expect_equal(test_weight_mat, weight_mat)

  enum <- TRUE
  boot_unique <- matrix(c(1,2))
  test_weight_mat <- gen_boot_weights(boot_dist = boot_dist, boot_unique = boot_unique, boot_reps = boot_reps, enum = enum)
  weight_mat <- data.frame(matrix(c(-1, -1, 1, -1, -1, 1, 1, 1), ncol = 4))
  row.names(weight_mat) <- c('Var1', 'Var2')

  expect_equal(test_weight_mat, weight_mat)


})

test_that('t_boot_p_val returns correct p-values', {

  H0 <- 1
  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(formula = Y ~ X, data = data)
  clusterby <- 'clusterby'
  x_interest <- 'X'
  beta0 <- coef(model)[x_interest]
  se0 <- clustered_se(data = data, model = model, clusterby = clusterby)[x_interest]
  beta <- c(0:10/5, beta0)
  se <- c(rep(0.1,11), se0)
  boot_reps <- 12
  absval <- FALSE
  bound <- c('upper', 'lower', 'mid', 'uniform', 'density')

  set.seed(42)
  test_p_val <- t_boot_p_val(se = se, beta = beta, H0 = H0, data = data, model = model,
                             clusterby = clusterby, x_interest = x_interest,
                             boot_reps = boot_reps, absval = absval, bound = bound)

  p_val <- c('upper' = 8/12, 'lower' = 6/12, 'mid' = 7/12, 'uniform' = (6+0.914806*2)/12, 'density' = 2*0.2666798)

  expect_equal(test_p_val, p_val, tolerance = 0.00001)

})

test_that('density_p_val returns correct p-values', {

  t <- c(-10.000000, -8.000000, -6.000000, -4.000000, -2.000000, 0.000000, 2.000000, 4.000000, 6.000000, 8.000000, 10.000000, -5.941251)
  t0 <- -5.941251
  boot_reps <- 12

  test_p_val <- density_p_val(t = t, t0 = t0, boot_reps = boot_reps)

  p_val <- 0.2666798

  expect_equal(test_p_val, p_val, tolerance = 0.00001)

})


test_that('cluster_meat calculates correct tensor', {

  x <- matrix(c(1, 1, 1, 2, 2, 2), ncol = 2, nrow = 3)
  X <- rbind(x, x)
  e <- matrix(1:6, ncol = 2, nrow = 3, byrow = TRUE)
  E <- rbind(e, e)
  clustervars <- c(1,1,1,2,2,2)
  k <- 2
  boot_reps <- 2

  cluster_meat <- 2*array(c(81, 162, 162, 324, 144, 288, 288, 576), c(2, 2, 2))

  test_cluster_meat <- suppressWarnings(cluster_meat(clustervars = clustervars, X = X, E = E, boot_reps = boot_reps, k = k))

  expect_equal(test_cluster_meat, cluster_meat, tolerance = 0.0001)

})

test_that('cluster_sandwich calculates correct tensor', {

  X <- cbind(1, 1:6)
  bread <- solve(crossprod(X))
  e <- matrix(1:6, ncol = 2, nrow = 3, byrow = TRUE)
  E <- rbind(e, e)
  clustervars <- c(1,1,1,2,2,2)
  k <- 2
  boot_reps <- 2

  cluster_sandwich <- c(3.236736, 5.551019)

  test_cluster_sandwich <- suppressWarnings(cluster_sandwich(clustervars = clustervars, X = X, bread = bread, x_ind = 2,
                                                         E = E, boot_reps = boot_reps, k = k, n = nrow(X), comb_n = 1))

  expect_equal(as.numeric(test_cluster_sandwich[2, ,2]), cluster_sandwich, tolerance = 0.0001)

})
