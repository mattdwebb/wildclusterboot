library(testthat)
library(wildclusterboot)

test_data <- read.csv(file = 'test_files/test_data.csv', stringsAsFactors = FALSE)
test_ywild <- read.csv('test_files/test_ywild.csv', stringsAsFactors = FALSE)
x_interest <- 'X'
bootby <- 'clusterby'
clusterby <- 'clusterby'
H0 <- 1
boot_dist <- 'six_pt'
model <- lm(data = test_data, formula = 'Y ~ X')

fitted_uhat_data <- read.csv(file = 'test_files/fitted_uhat.csv', stringsAsFactors = FALSE)

test_mle <- wild_clust_mle(model = model,
                           x_interest = x_interest,
                           bootby = bootby,
                           boot_dist = boot_dist,
                           H0 = H0)

test_that('wild_clust_mle returns proper list',{

  expect_equal(length(test_mle), 5)
  expect_equal(as.numeric(test_mle[['uhat']]), as.numeric(fitted_uhat_data[,'er']), tolerance = 0.0001)
  expect_equal(as.numeric(test_mle[['fitted_data']]), as.numeric(fitted_uhat_data[,'xbr']), tolerance = 0.0001)
  expect_equal(test_mle[['y_name']], 'Y')
  expect_equal(test_mle[['bootby']], bootby)
  expect_equal(test_mle[['boot_dist']], c(-sqrt(3/2),-sqrt(2/2),-sqrt(1/2),sqrt(1/2),sqrt(2/2),sqrt(3/2)))

})

set.seed(42)
ran_data <- wild_clust_ran(data = test_data, mle = test_mle)

mult_test_data <- read.csv(file = 'test_files/mult_test_data.csv', stringsAsFactors = FALSE)
mult_test_ywild <- read.csv(file = 'test_files/mult_test_ywild.csv', stringsAsFactors = FALSE)
mult_test_mle <- wild_clust_mle(model = model,
                           x_interest = x_interest,
                           bootby = ~ clusterby + clusterby2,
                           boot_dist = boot_dist,
                           H0 = H0)

set.seed(42)
mult_ran_data <- wild_clust_ran(data = mult_test_data, mle = mult_test_mle)

test_that('wild_clust_ran returns returns the correct dataframe',{

  expect_equal(nrow(ran_data), 100)
  expect_equal(ncol(ran_data), 3)
  expect_equal(as.numeric(ran_data[, 'X']), as.numeric(test_data[, 'X']))
  expect_equal(as.numeric(ran_data[, 'Y']), as.numeric(test_ywild[, 'Y']), tolerance = 0.0001)
  expect_equal(as.numeric(ran_data[, 'clusterby']), as.numeric(test_ywild[, 'clusterby']))

  expect_equal(nrow(mult_ran_data), 100)
  expect_equal(ncol(mult_ran_data), 4)
  expect_equal(as.numeric(mult_ran_data[, 'X']), as.numeric(mult_test_data[, 'X']))
  expect_equal(as.numeric(mult_ran_data[, 'Y']), as.numeric(mult_test_ywild[, 'Y']), tolerance = 0.0001)
  expect_equal(as.numeric(mult_ran_data[, 'clusterby']), as.numeric(mult_test_ywild[, 'clusterby']))

})

model <- lm(data = test_ywild, formula = 'Y ~ X')
test_se <- wild_clust_statistic(data = test_ywild, model = model, x_interest = x_interest, clusterby = clusterby, H0 = H0)

test_that('wild_clust_statistic returns correct se',{

  expect_equal(length(test_se), 1)
  expect_equal(test_se, -1.5917145, tolerance = 0.000001)

})

t <- read.csv(file = 'test_files/test_boot.csv', stringsAsFactors = FALSE)
t0 <- c(-0.4443353, -0.3948289)
R <- 399
boot.out <- list(t = t, t0 = t0, R = R)
test_t <- boot_p_val(boot.out)

test_that('boot_p_vals returns the correct values', {

  expect_equal(test_t, c(0.6165414, 0.601504), tolerance = 0.0001)
  expect_equal(length(test_t), 2)

})
