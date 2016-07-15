
<!-- README.md is generated from README.Rmd. Please edit that file -->
wildclusterboot
===============

wildclusterboot is a package designed to implement clustered standard errors and the wild cluster bootstrap based on Cameron, Gelbach, and Miller (2008) while adding the ability to specify any bootstrap distribution and to cluster and bootstrap on different variables.

Installation
------------

The easiest way to install this package is via devtools. To install devtools, run the following code:

    install.packages('devtools')

Once devtools is installed, go get a GitHub [personal API key](https://github.com/settings/tokens) (you only need to give it repo access). Make sure to copy the generated key and then run the following code:

    devtools::install_github(repo = 'scottmcneil/wildclusterboot', auth_token = [auth token you just generated])

Usage
-----

The library provides three main functions. The first is `clustered_se`, which can be used to get the clustered standard errors out of an `lm` object:

    #Create dataframe of random normal variables
    test_data <- data.frame(Y = rnorm(10), X = rnorm(10), G = c(1,1,1,2,2,2,3,3,3,3))

    #Fit model for test data
    test_model <- lm(data = test_data, formula = 'Y ~ X')

    #Get a vector fo clustered standard errors for model
    SEs <- clustered_se(model = test_model, clusterby = test_data$G)

The second is `wild_cluster_boot`, which can be used to perform a wild clustered bootstrap on an lm object and a set of data for a given independent variable of interest:

    test.boot.out <- wild_cluster_boot(data = test_data,
                                       model = test_model,
                                       x_interest = 'X',
                                       clusterby = 'G',
                                       bootby = 'G',
                                       boot_dist = 'six_pt',
                                       boot_reps = 399,
                                       cores = 3)

The `clusterby` and `bootby` arguments can be different. The `clusterby` variable can also be a vector of strings, allowing you to calculate standard errors on multiple group variables at once.

Finally, the `boot_p_val` function takes a `boot.out` object from `wild_cluster_boot` and returns a vector of p-values for each `clusterby` variable:

    p_vec <- boot_p_val(test.boot.out)
