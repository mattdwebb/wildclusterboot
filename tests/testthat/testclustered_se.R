library(testthat)
library(wildclusterboot)

test_that('calc_cluster_se calculates proper values',{

  #Set out X matrix, Y vector, uhat vector and clusterby vector
  X <- c(1.9157903, 1.9157903, 1.4053975, 1.4053975, 0.3096093, 0.3096093, -0.2007836, -0.2007836)
  X1 <- cbind('(Intercept)' = 1, X)
  Y <- c(5.3051725, 3.2240347, 0.3796652, 1.4897371, 2.4207392, 3.4025150, 0.6618861, 1.4311579)
  uhat <- c(2.1504555, 0.0693176, -2.3577072, -1.2476353, 0.5793847, 1.5611605, -0.7621238, 0.0071479)
  clusterby <- c(1, 1, 2, 2, 1, 1, 2, 2)

  #Create vector of what result should be
  test_se <- c('(Intercept)' = 0.7938156, 'X' = 0.4473510)

  #Calculate se's and test if results correc
  se <- calc_cluster_se(X = X1, uhat = uhat, clusterby = clusterby, cluster_ids = c(1,2), G = 2)
  expect_equal(se, test_se, tolerance = 0.0001)
  expect_equal(length(se),2)

  #Create extra column and check to make sure dimensions correct
  X2 <- cbind(X1, x_2 = rnorm(length(X)))
  se <- calc_cluster_se(X = X2, uhat = uhat, clusterby = clusterby, cluster_ids = c(1,2), G = 2)
  expect_equal(length(se),3)

})

test_that('clustered_se calculates proper values',{

  #Set out X and Y vectors
  X <- c(1.9157903, 1.9157903, 1.4053975, 1.4053975, 0.3096093, 0.3096093, -0.2007836, -0.2007836)
  Y <- c(5.3051725, 3.2240347, 0.3796652, 1.4897371, 2.4207392, 3.4025150, 0.6618861, 1.4311579)

  #Set out clusterby vector
  G <- c(1, 1, 2, 2, 1, 1, 2, 2)

  dat <- data.frame(X = X, Y = Y, G = G)

  #Result of what se vector should be
  test_se <- c('(Intercept)' = 0.7938156, 'X' = 0.4473510)

  #Generate model and test if result is expected value and size
  model <- lm(Y ~ X)
  se <- clustered_se(data = dat, model = model, clusterby = 'G')
  expect_equal(se, test_se, tolerance = 0.0001)
  expect_equal(length(se),2)

  #Create second X value and check that resulting vector correct size and names
  X2 <- rnorm(length(X))
  model <- lm(Y ~ X + X2, x = TRUE)
  se <- clustered_se(model = model, clusterby = 'G')
  expect_equal(length(se), 3)
  expect_equal(names(se), c('(Intercept)', 'X', 'X2'))

  #Make sure throws correct error on model not lm class
  expect_error(clustered_se(model = '', clusterby = 'G'),
               'Model variable must be lm class')

})

test_that('clustered_se calculates proper values',{

  #Set out X and Y vectors
  X <- c(1.9157903, 1.9157903, 1.4053975, 1.4053975, 0.3096093, 0.3096093, -0.2007836, -0.2007836)
  Y <- c(5.3051725, 3.2240347, 0.3796652, 1.4897371, 2.4207392, 3.4025150, 0.6618861, 1.4311579)

  #Set out clusterby vector
  G <- c(1, 1, 2, 2, 1, 1, 2, 2)
  H <- c(1, 2, 1, 2, 1, 2, 1, 2)

  dat <- data.frame(X = X, Y = Y, G = G, H = H)

  #Result of what se vector should be
  test_se <- c('(Intercept)' = 0.8054447, 'X' = 0.5150004)

  #Generate model and test if result is expected value and size
  model <- lm(Y ~ X)
  se <- clustered_se(data = dat, model = model, clusterby = ~ G + H)
  expect_equal(se, test_se, tolerance = 0.0001)
  expect_equal(length(se), 2)

})
