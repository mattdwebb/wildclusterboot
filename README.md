
<!-- README.md is generated from README.Rmd. Please edit that file -->
wildclusterboot
===============

wildclusterboot is a package designed to implement clustered standard errors and the wild cluster bootstrap based on Cameron, Gelbach, and Miller (2008) while adding the ability to specify any bootstrap distribution and to cluster and bootstrap on different variables.

Installation
------------

The easiest way to install this package is via devtools. To install devtools, run the following code:

    install.libraries(devtools)

Once devtools is installed, go get a GitHub [personal API key](https://github.com/settings/tokens) (you only need to give it repo access). Make sure to copy the generated key and then run the following code:

    devtools::install_github(repo = 'scottmcneil/wildclusterboot', auth_token = [auth token you just generated])
