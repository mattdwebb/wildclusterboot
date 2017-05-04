test_that('wild_cluster_boot returns correct value',{

  data <- read.csv('test_files/test_boot_p.csv')
  model <- lm(data = data, formula = Y ~ X)
  boot_reps <- 50
  clusterby <- 'clusterby'
  x_interest <- 'X'
  boot_dist <- c(-1, 1)
  H0 <- 1

  set.seed(42)
  test_p <- wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.52))

  set.seed(36)
  test_p <- wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.48))


  set.seed(128)
  test_p <- wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
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
  test_p <- wild_cluster_boot(data = data, model = model, x_interest = x_interest, clusterby = clusterby,
                                boot_dist = boot_dist, boot_reps = boot_reps, H0 = H0)

  expect_equal(test_p, c(upper = 0.08))

})
