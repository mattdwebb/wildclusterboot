###### TESTS eigen_fix_cpp FUNCTION ######

test_that('eigen_fix_cpp returns correct sandiwhich matrix',{

  sandwich <- matrix(c(-2.021514e-05, -1.109287e-05, -1.109287e-05, 1.916697e-05), nrow = 2)

  test_sandwich <- matrix(c(1.421063e-06, -5.417819e-06, -5.417819e-06, 2.065550e-05), nrow = 2)
  fixed_sandwich <- eigen_fix_cpp(sandwich)
  eigen_fix_cpp(sandwich)

  expect_equal(test_sandwich, fixed_sandwich)

})
