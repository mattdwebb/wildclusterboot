###### TESTS wild_data FUNCTION ######

test_that('wild_data works with H0=0 and no covariates',{

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  x_interest <- 'X'
  H0 <- 0

  fitted_data <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  short_model <- lm(data = data, formula = Y ~ 1)

  data['fitted_data'] <- predict(short_model)
  data['uhat'] <- resid(short_model)


  expect_equal(data, fitted_data)

})

test_that('wild_data works with H0=1 and no covariates',{

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  x_interest <- 'X'
  H0 <- 1

  fitted_data <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  data['Y'] <- data['Y'] - data['X']*H0
  short_model <- lm(data = data, formula = Y ~ 1)
  data['fitted_data'] <- predict(short_model) + data['X']*H0
  data['uhat'] <- resid(short_model)
  data['Y'] <- data['Y'] + data['X']*H0

  expect_equal(data, fitted_data)

})

test_that('wild_data works with H0=0 and one covariate',{

  data <- read.csv('test_files/test_boot_p.csv')
  data['W'] <- rnorm(10)
  model <- lm(data = data, formula = Y ~ X + W)
  x_interest <- 'X'
  H0 <- 0

  fitted_data <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  short_model <- lm(data = data, formula = Y ~ W)
  data['fitted_data'] <- predict(short_model) + data['X']*H0
  data['uhat'] <- resid(short_model)

  expect_equal(data, fitted_data)

})

test_that('wild_data works with H0=1 and one covariate',{

  data <- read.csv('test_files/test_boot_p.csv')
  data['W'] <- rnorm(10)
  model <- lm(data = data, formula = Y ~ X + W)
  x_interest <- 'X'
  H0 <- 1

  fitted_data <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  data['Y'] <- data['Y'] - data['X']*H0
  short_model <- lm(data = data, formula = Y ~ W)
  data['fitted_data'] <- predict(short_model) + data['X']*H0
  data['uhat'] <- resid(short_model)
  data['Y'] <- data['Y'] + data['X']*H0

  expect_equal(data, fitted_data)

})
