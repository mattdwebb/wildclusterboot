library(testthat)
library(wildclusterboot)

test_that('crve_sandwich works with 2 group dimensions',{

  #Set out X and Y vectors
  data <- data.frame(X = c(1.9157903, 1.9157903, 1.4053975, 1.4053975, 0.3096093, 0.3096093, -0.2007836, -0.2007836),
                     Y = c(5.3051725, 3.2240347, 0.3796652, 1.4897371, 2.4207392, 3.4025150, 0.6618861, 1.4311579),
                     G = c(1, 1, 2, 2, 1, 1, 2, 2),
                     H = c(1, 2, 1, 2, 1, 2, 1, 2),
                     GH = c(0, 1, 2, 3, 0, 1, 2, 3))

  cluster_ind <- as.matrix(data.frame(G = c(0, 0, 1, 1, 0, 0, 1, 1),
                                      H = c(0, 1, 0, 1, 0, 1, 0, 1),
                                      GH = c(0, 1, 2, 3, 0, 1, 2, 3)))


  #Result of what se vector should be
  test_var <- c(0.6487412, 0.2652254)

  #Generate model and test if result is expected value and size
  model <- lm(Y ~ X, data = data)
  uhat <- resid(model)
  X <- model.matrix(model)
  bread <- solve(crossprod(X))
  comb <- c(1, 1, -1)
  G <- c(2, 2, 4)
  var_mat <- crve_sandwich(X, bread, uhat, cluster_ind, comb, G)
  expect_equal(diag(var_mat), test_var, tolerance = 0.0001)

})

test_that('crve_sandwich works with 1 group dimension',{

  #Set out X and Y vectors
  data <- data.frame(X = c(1.9157903, 1.9157903, 1.4053975, 1.4053975, 0.3096093, 0.3096093, -0.2007836, -0.2007836),
                     Y = c(5.3051725, 3.2240347, 0.3796652, 1.4897371, 2.4207392, 3.4025150, 0.6618861, 1.4311579),
                     G = c(1, 1, 2, 2, 1, 1, 2, 2))

  cluster_ind <- as.matrix(data.frame(G = c(0, 0, 1, 1, 0, 0, 1, 1)))


  #Result of what se vector should be
  test_var <- c(0.6301431, 0.2001230)

  #Generate model and test if result is expected value and size
  model <- lm(Y ~ X, data = data)
  uhat <- resid(model)
  X <- model.matrix(model)
  bread <- solve(crossprod(X))
  comb <- 1
  G <- 2
  var_mat <- crve_sandwich(X, bread, uhat, cluster_ind, comb, G)
  expect_equal(diag(var_mat), test_var, tolerance = 0.0001)

})
