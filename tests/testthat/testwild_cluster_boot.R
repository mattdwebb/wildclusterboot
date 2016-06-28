library(testthat)
library(wildclusterboot)

test_data <- read.csv(file = 'test_files/test_data.csv', stringsAsFactors = FALSE)
y <- test_data[, 'Y']
X <- test_data[, 'X']
bootby <- test_data[, 'clusterby']
clusterby <- rep(c(1,2), each = 50)
param <- 'X'
H0 <- 1
boot_dist <- c(-sqrt(3/2),-sqrt(2/2),-sqrt(1/2),sqrt(1/2),sqrt(2/2),sqrt(3/2))
model <- lm(y ~ X)

fitted_uhat_data <- read.csv(file = 'test_files/fitted_uhat.csv', stringsAsFactors = FALSE)
test_mle <- wild_clust_mle(model = model,
                         param = param,
                         clusterby = clusterby,
                         bootby = bootby,
                         boot_dist = boot_dist,
                         H0 = H0)

test_that('wild_clust_mle returns proper list',{

  expect_equal(length(test_mle), 9)
  expect_equal(as.numeric(test_mle[['uhat']]), as.numeric(fitted_uhat_data[,'er']), tolerance = 0.0001)
  expect_equal(as.numeric(test_mle[['fitted_data']]), as.numeric(fitted_uhat_data[,'xbr']), tolerance = 0.0001)
  expect_equal(test_mle[['y_name']], 'y')
  expect_equal(test_mle[['param']], param)
  expect_equal(test_mle[['clusterby']], clusterby)
  expect_equal(test_mle[['bootby']], bootby)
  expect_equal(test_mle[['boot_dist']], boot_dist)
  expect_equal(test_mle[['bG']], length(unique(bootby)))
  expect_equal(test_mle[['boot_ids']], unique(bootby))

})

set.seed(42)
rand_data <- wild_clust_ran(data = cbind(y, X, clusterby), mle = test_mle)
test_ywild <- read.csv('test_files/test_ywild.csv', stringsAsFactors = FALSE)

test_that('wild_clust_ran returns the correct values',{

  expect_equal(nrow(rand_data), 100)
  expect_equal(ncol(rand_data), 3)
  expect_equal(as.numeric(rand_data[, 2]), as.numeric(test_data[, 'X']))
  expect_equal(as.numeric(rand_data[, 1]), as.numeric(test_ywild[, 'y_wild']), tolerance = 0.0001)
  expect_equal(as.numeric(rand_data[, 3]), as.numeric(clusterby))

})

test_se <- wild_clust_statistic(as.matrix(test_ywild))

test_that('wild_clust_statistic returns correct se',{

  expect_equal(length(test_se), 2)
  expect_equal(as.numeric(test_se['se']), 0.4034162, tolerance = 0.000001)
  expect_equal(as.numeric(test_se['beta']), 0.3578766, tolerance = 0.000001)

})
