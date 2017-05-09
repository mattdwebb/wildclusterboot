###### TESTS wild_y FUNCTION ######

library(testthat)
library(wildclusterboot)

test_that('wild_y calculates correct matrix', {

  set.seed(42)

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  x_interest <- 'X'
  H0 <- 0

  data_wild <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)
  bootby <- 'clusterby'
  boot_dist <- 'two_pt'
  boot_reps <- 4
  enum <- F

  test_y_wild <- wild_y(data_wild = data_wild, bootby = bootby, boot_dist = boot_dist,
                        boot_reps = boot_reps, enum = enum)

  y_wild <- matrix(c(2.197272152,
                     1.243391271,
                     -1.50428824,
                     -0.124482112,
                     -0.207608963,
                     -0.036807853,
                     0.660635341,
                     -3.317748194,
                     -0.98867695,
                     0.690765512,
                     -2.474781758,
                     -1.520900878,
                     1.226778633,
                     -0.153027495,
                     -0.069900644,
                     -0.036807853,
                     0.660635341,
                     -3.317748194,
                     -0.98867695,
                     0.690765512,
                     2.197272152,
                     1.243391271,
                     -1.50428824,
                     -0.124482112,
                     -0.207608963,
                     -0.036807853,
                     0.660635341,
                     -3.317748194,
                     -0.98867695,
                     0.690765512,
                     2.197272152,
                     1.243391271,
                     -1.50428824,
                     -0.124482112,
                     -0.207608963,
                     -0.240701754,
                     -0.938144948,
                     3.040238587,
                     0.711167343,
                     -0.968275119),
                   ncol = 4)

  expect_equal(y_wild, test_y_wild)

})
