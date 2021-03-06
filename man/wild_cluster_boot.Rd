% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wild_cluster_boot.R
\name{wild_cluster_boot}
\alias{wild_cluster_boot}
\title{Calculate wild cluster bootstrap}
\usage{
wild_cluster_boot(data, model, x_interest, clusterby, boot_dist, boot_reps,
  bootby = clusterby, H0 = 0, enum = FALSE, absval = FALSE,
  bound = c("upper", "lower", "mid", "uniform", "density"))
}
\arguments{
\item{data}{Dataframe with all data, including group indices}

\item{model}{lm object to perform the boostrap with}

\item{x_interest}{String indicating the name of the paramater of interest}

\item{clusterby}{String or formula with name of clusterby variable in data.
Use a string for one-way clustering and a formula to indicate
 multiway clustering (e.g. use ~ G + H to a clustering
 using the G and H dimensions}

\item{boot_dist}{Specifies the distribution from which to draw the wild boostrap weights.
This can be one of three things: a numeric vector,
a string specifying default distribution ("rad" or
"six" for now)vor a function that takes only the argument
"n" (e.g. rnorm)}

\item{boot_reps}{Integer indicating number of bootstrap repetitions}

\item{bootby}{String with name of bootby variable.
The default for this is same as clusterby, however
a single dimension must be specified if the clusterby
variable is a multi-way formula}

\item{H0}{Float or integer inticating the null hypothesis, default is 0}

\item{enum}{Boolean indicating whether to calculate all possible wild bootstrap combinations.
Only valid if using a vector or default distribution.
If this is set to TRUE, then the boot_reps variable
is ignored and a new enumerated total is calculated}

\item{absval}{Boolean indicating whether or not to use absolute valued t-statistics}

\item{bound}{String or vector of string indicating which bootstrap "bound" to use when tie-breaking.
upper and lower indicate using the highest or lowest
value when actual value is tied with bootstrapped values.
mid takes the half-way point between the two. uniform
calcualtes a random interval betwen the two and density
uses a kernel smoothing correction due to Racine-MacKinnon
(2007)}
}
\value{
p-value or vector of p-values corresponding to bootstrap result
}
\description{
This function provides a flexible framework for
calculating the wild cluster bootsrap, this function
implements the original wild cluster boostrap, due
to Cameron, Gelbach and Miller (2008), and a number of
modifications, including flexible boostrap distributions
and implementing a "multi-way" version.
}
