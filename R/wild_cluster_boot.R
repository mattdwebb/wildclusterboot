#' Performs wild clustered bootstrap
#'
#' @param data Dataframe with all data, including group indices
#' @param model lm object of interest
#' @param x_interest X paramater of interest
#' @param clusterby String or formila with name of clusterby variable in data
#' @param boot_dist Vector of weights for wild bootstrap, string specifying default distribution or a function that only takes the argument "n"
#' @param boot_reps Number of repititions for resampling data
#' @param bootby String with name of bootby variable in data, default is same as clusterby
#' @param H0 Float of integer inticating the null hypothesis, default is 0
#' @param enum Boolean indicating whether to calculate all possible wild bootstrap combinations, will override boot_reps and report upper and lower bounds
#' @param absval Boolean indicating whether or not to use absolute valued t-statistics
#' @param bound Boolean indicating whether or not to use bound-MacKinnon correction
#' @return p-value corresponding to bootstrap result
#' @importFrom magrittr "%>%"
#' @export
wild_cluster_boot <- function(data, model, x_interest, clusterby, boot_dist, boot_reps, bootby = clusterby, H0 = 0, enum = FALSE, absval = FALSE, bound = c('upper', 'lower', 'mid', 'uniform', 'density')){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  #Check to ensure no incorrect bound values supplied
  bound <- if(missing(bound)) 'upper' else match.arg(bound, several.ok = TRUE)

  #Get dataframe from model
  data <- eval(model$call$data)

  #Get wild bootstrap fitted data
  data_wild <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  #Generate matrix of y-wild values
  y_wild <- wild_y(data_wild = data_wild, bootby = bootby, boot_dist = boot_dist, boot_reps = boot_reps, enum = enum)

  #Get model matrix from model
  X <- model.matrix(model)

  #Set x_ind as index from model matrix
  x_ind <- grep(x_interest, colnames(X))

  #Create bread, B and E matrices
  bread <- bread_cpp(X)
  B <- beta_cpp(X, bread, y_wild)
  E <- y_wild - X %*% B

  se <- wild_se(data_wild, E, X, bread, clusterby, x_ind)

  #Get beta for x of interest
  beta <- B[x_ind, ]

  p_value <- boot_p_val(se = se, beta = beta, H0 = H0, data = data, model = model,
                          clusterby = clusterby, x_interest = x_interest,
                          boot_reps = boot_reps, bound = bound, absval = absval)

  return(p_value)

}


#' Return wild bootstrap residuals and fitted values
#'
#' @param model lm object of interest
#' @param x_interest X paramater of interest
#' @param H0 Float of integer inticating the null hypothesis, default is 0
#' @return dataframe with original data and added fitted values and residuals
wild_data <- function(data, model, x_interest, H0){

  #Set out model variable names
  model_vars <- all.vars(formula(model))
  y_name <- model_vars[1]
  short_vars <- model_vars[model_vars != x_interest]
  short_x <- short_vars[short_vars != y_name]

  #H0 transformation
  h0_mutate <- lazyeval::interp(~ a - c*b, a = as.name(y_name), b = as.name(x_interest), c = H0)

  #Imposes null hypothesis
  short_data <- data %>%
    dplyr::mutate_(.dots = setNames(list(h0_mutate), y_name))

  #Check if short model is just intercept, otherwise create short formula
  short_formula <- if(length(short_x) == 0){
    reformulate(termlabels = c('1'), response = y_name)
  } else{
    reformulate(termlabels = short_x, response = y_name)
  }

  #Estimate short model
  short_model <- lm(data = short_data, formula = short_formula)

  reverse_h0_mutate <- lazyeval::interp(~ a + c*b, a = as.name(y_name), b = as.name(x_interest), c = H0)
  fitted_mutate <- lazyeval::interp(~ fitted_data + c*b, b = as.name(x_interest), c = H0)

  #Get residuals and fitted data from short model
  fitted_data <- short_data %>%
    modelr::add_predictions(short_model) %>%
    modelr::add_residuals(short_model) %>%
    dplyr::rename(uhat = resid,
                  fitted_data = pred) %>%
    dplyr::mutate_(.dots = setNames(list(fitted_mutate), 'fitted_data')) %>%
    dplyr::mutate_(.dots = setNames(list(reverse_h0_mutate), y_name))

  return(fitted_data)

}

#' Return matrix of wild response variables
#'
#' @param data_wild Dataframe of original data with wild fitted data and residuals
#' @param bootby string or formula indicating variables in data_wild to group bootstrap weights by
#' @param boot_dist Vector of weights for wild bootstrap, string specifying default distribution or a function that only takes the argument "n"
#' @param boot_reps Number of repititions for resampling data
#' @param enum Boolean indicating whether to calculate all possible wild bootstrap combinations, will override boot_reps
#' @return Matrix of bootstrap weights
#'
wild_y <- function(data_wild, bootby, boot_dist, boot_reps, enum){

  #TODO: Split this out as a separate function
  #If boot_dist is character, resolve to default distribution
  if(class(boot_dist) == 'character'){

    default_dist <- list(six_pt = c(-sqrt(3/2),-sqrt(2/2),-sqrt(1/2),sqrt(1/2),sqrt(2/2),sqrt(3/2)),
                         two_pt = c(-1,1))

    if(!boot_dist %in% names(default_dist)){

      error_message <- paste('Only the following default distributions supported:', paste(names(default_dist), collapse = ', '))
      stop(error_message)

    }

    boot_dist <- default_dist[[boot_dist]]

  }

  bootby <- if(class(bootby) == 'formula') all.vars(bootby) else bootby

  #Create unique vector of group ids
  boot_unique <- unique(data_wild[bootby])
  boot_paste <- Reduce(paste0, data_wild[bootby])
  boot_ind <- as.numeric(factor(boot_paste, levels = unique(boot_paste)))-1

  boot_reps <- if(enum) 2^nrow(boot_unique) else boot_reps

  # #Add weights for each group
  weights <- gen_boot_weights(boot_dist = boot_dist, boot_unique = boot_unique, boot_reps = boot_reps, enum = enum) %>%
    setNames(paste0('wild_boot_weight_', 1:boot_reps))

  weight_names <- names(weights)
  boot_weights <- cbind(boot_unique, weights)
  expanded_weights <- merge(x = data_wild, y = boot_weights)

  #Generate matrix of y-wild values
  y_wild <- as.matrix(expanded_weights[,weight_names]) * data_wild[,'uhat'] + data_wild[,'fitted_data']

  return(y_wild)

}
#' Calculate vector of clustered standard errors
#'
#' @param data_wild Dataframe of original data with wild fitted data and residuals
#' @param E Matrix of residuals, which each column representing a bootstrap replication
#' @param X Model matrix
#' @param bread X*X' matrix
#' @param clusterby String or formila with name of clusterby variable in data
#' @param x_ind integer indicating the index of x of interest in model matrix
#' @return Vector of standard errors
#'
wild_se <- function(data_wild, E, X, bread, clusterby, x_ind){

  clusterby_list <- if(class(clusterby) == 'formula') all.vars(clusterby) else clusterby
  comb_vars <- unlist(lapply(X = 1:length(clusterby_list), FUN = combn, x = clusterby_list, simplify = FALSE), recursive = FALSE)

  cluster_ind <- sapply(X = comb_vars, FUN = crve_ind, data_wild = data_wild)
  comb <- (-1)^(sapply(X = comb_vars, length)+1)
  G <- apply(cluster_ind, 2, function(x) length(unique(x)))

  sandwich_array_cpp <- mapply(uhat = data.frame(E), FUN = crve_sandwich, MoreArgs = list(X = X, bread = bread, clusterby = cluster_ind, comb, G))
  se_cpp <- sqrt(sandwich_array_cpp[x_ind^2,])

  return(se_cpp)

}


#' Generate bootstrap weight matrix
#'
#' @param boot_dist Vector of weights for wild bootstrap, string specifying default distribution or a function that only takes the argument "n"
#' @param boot_unique Matrix of unique group IDs
#' @param boot_reps Number of repititions for resampling data
#' @param enum Boolean indicating whether to calculate all possible wild bootstrap combinations, will override boot_reps and report upper and lower bounds
#' @return Matrix of bootstrap weights
#'
crve_ind <- function(data_wild, clusterby){

  clustervars <- Reduce(data_wild[clusterby], f = paste0)
  cluster_ind <- as.numeric(factor(clustervars, levels = unique(clustervars)))-1

}


#' Calculate bootstrap p-values
#'
#' @param se Vector of standard errors
#' @param beta Vector of coefficients for x of interest
#' @param H0 Float of integer inticating the null hypothesis, default is 0
#' @param data Dataframe of initial data
#' @param model lm object of interest
#' @param clusterby String or formula indicating variable in data for clustering
#' @param x_interest X paramater of interest
#' @param boot_reps Integer for number of replications
#' @param absval Boolean indicating whether or not to use absolute valued t-statistics
#' @param bound Character or character vector indicating which bootstrap p-values to generate
#' @return Bootstrap p-value
#'
boot_p_val <- function(se, beta, H0, data, model, clusterby, x_interest, boot_reps, absval, bound){

  #Calculate t vector
  t <- if(absval) abs((beta - H0)/se) else (beta - H0)/se

  #Calculate initial t
  se0 <- clustered_se(data = data, model = model, clusterby = clusterby)
  beta0 <- coef(model)
  t0 <- if(absval) abs((beta0[x_interest] - H0)/se0[x_interest]) else (beta0[x_interest] - H0)/se0[x_interest]

  p <- sapply(bound, p_eval, t = t, t0 = t0,  boot_reps = boot_reps, absval = absval)

  return(p)

}

#' Calculate bootstrap p-value
#'
#' @param t Vector of t-values
#' @param t0 Original t-value
#' @param boot_reps Integer for number of replications
#' @param bound Boolean indicating whether or not to use bound-MacKinnon correction
#' @param absval Boolean indicating whether or not to use absolute valued t-statistics
#' @return Bootstrap p-value
#'
p_eval <- function(t, t0, boot_reps, bound, absval){

  #Get proportion below and equal, while accounting for rounding errors
  prop_less <- sum((t - t0) < -1e-12)/boot_reps
  prop_equal <- sum(abs(t0 - t) < 1e-12)/boot_reps

  prop <- switch (bound,
                  'upper' = prop_less + prop_equal,
                  'lower' = prop_less,
                  'mid' = prop_less + prop_equal/2,
                  'uniform' = prop_less + runif(1)*prop_equal,
                  'density' = density_p_val(t = t, t0 = t0, boot_reps = boot_reps)
  )

  p <- if(absval) prop else 2*min(prop, 1 - prop)

  return(p)

}

#' Calculate density p-value
#'
#' @param t Vector of t-values
#' @param t0 Original t-value
#' @param boot_reps Integer for number of replications
#' @return Bootstrap p-value
#'
density_p_val <- function(t, t0, boot_reps){

  h <- mlcv(t)

  gauss <- (t0-t)/h

  p <- sum(pnorm(q = gauss))/boot_reps

  return(p)

}

#' Calculate best kernel for density function
#'
#' @param t Vector of t-values
#'

mlcv <- function(t) {

  n <- length(t)

  kdenest.mlcv <- function(h) {

    p_val <- sapply(X = t, FUN = function(x, t) sum(dnorm((x-t)/h), -dnorm(0))/((n-1)*h), t = t)

    p <- if(h > 0) {
      -sum(log(ifelse(p_val > 0, p_val,.Machine$double.xmin)))/n
    } else {
      .Machine$double.xmax
    }

    return(p)
  }

  return(nlm(kdenest.mlcv, 1.06*sd(t)*n^{-1/5})$estimate*n^{-2/15})

}
