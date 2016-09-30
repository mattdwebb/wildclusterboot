#' Performs wild clustered bootstrap
#'
#' @param data Dataframe with all data, including group indices
#' @param model lm object of interest
#' @param x_interest X paramater of interest
#' @param clusterby String or formila with name of clusterby variable in data
#' @param boot_dist Vector of weights for wild bootstrap or string specifying default distribution
#' @param boot_reps Number of repititions for resampling data
#' @param bootby String with name of bootby variable in data, default is same as clusterby
#' @param H0 Float of integer inticating the null hypothesis, default is 0
#' @return p-value corresponding to bootstrap result
#' @export
t_wild_cluster_boot <- function(data, model, x_interest, clusterby, boot_dist, boot_reps, bootby = clusterby, H0 = 0){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

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

  #Set out model variable names
  model_vars <- all.vars(formula(model))
  y_name <- model_vars[1]
  short_vars <- model_vars[!model_vars %in% c(x_interest, y_name)]

  #This will impose the H0
  short_data <- data.frame(Y = data[,y_name] - H0* data[,x_interest], X = data[, x_interest])

  #Check if short model is just intercept, otherwise create short formula
  if(length(short_vars) == 0){
    form <- reformulate(termlabels = c('1'), response = y_name)
  } else{
    form <- reformulate(termlabels = short_vars, response = y_name)
  }

  #Create short data frame and estimate short model
  short_model <- lm(data = short_data, formula = form)

  #Get residuals and fitted data from short model
  uhat <- resid(short_model)
  fitted_data <- fitted(short_model) + H0 * data[, x_interest]

  if(class(bootby) == 'formula'){

    bootby <- all.vars(bootby)

  }

  #Create unique vector of group ids
  boot_unique <- unique(data[bootby])

  #Add weights for each group
  weights <- data.frame(matrix(sample(x = boot_dist, size = nrow(boot_unique)*boot_reps, replace = TRUE),
                               nrow = nrow(boot_unique), ncol = boot_reps, byrow = FALSE))

  weight_names <- names(weights)
  boot_weights <- cbind(boot_unique, weights)
  expanded_weights <- merge(x = data, y = boot_weights)

  #Generate matrix of y-wild values
  y_wild <- as.matrix(expanded_weights[,weight_names]) * uhat + fitted_data

  #Get model matrix from model
  X <- model.matrix(model)

  #Set x_ind as index from model matrix
  x_ind <- grep(x_interest, colnames(X))

  #Create bread, B and E matrices
  bread <- solve(crossprod(X))
  B <- bread %*% t(X) %*% y_wild
  E <- y_wild - X %*% B

  #Check if supplied clusterby is a formula
  if(class(clusterby) == 'formula'){

    #extract variables from formula
    form_vars <- all.vars(clusterby)

    #Create combinations of group dimensions as well as concatenated names
    comb_list <- lapply(X = 1:length(form_vars), FUN = combn, x = form_vars, simplify = FALSE)
    comb_vars <- unlist(comb_list, recursive = FALSE)
    comb_names <- lapply(X = comb_vars, FUN = paste0, collapse = '')
    comb_n <- lapply(X = comb_vars, length)

    #Create sub-dataframes for each combination of group dimensions
    data_combs <- lapply(X = comb_vars, function(comb) data[comb])

    #Create matrix of all group dimension combinations
    clustervars <- lapply(X = data_combs, FUN = Reduce, f = paste0)
    names(clustervars) <- comb_names

  } else{

    clustervars <- data[clusterby]
    comb_n <- 1

  }

  #Set out k and n variables from model matrix
  n <- nrow(X)
  k <- ncol(X)

  #Calculate sandwitch matrix for each cluster combinations
  sandwich_list <- suppressWarnings(mapply(FUN = cluster_sandwich, clustervars = clustervars, comb_n = comb_n,
                                           MoreArgs = list(X = X, bread = bread, x_ind = x_ind, E = E, k = k, boot_reps = boot_reps, n = n), SIMPLIFY = FALSE))

  #Reduce by addtion and convert to vector
  se <- sqrt(Reduce('+', sandwich_list))
  beta <- B[x_ind, ]

  p_value <- t_boot_p_val(se = se, beta = beta, H0 = H0, data = data, model = model,
                          clusterby = clusterby, x_interest = x_interest, boot_reps = boot_reps)

}

#TODO Create bootstrap p-value function

#' Calculate bootstrap p-value
#'
#' @param se Vector of standard errors
#' @param beta Vector of coefficients for x of interest
#' @param H0 Float of integer inticating the null hypothesis, default is 0
#' @param data Dataframe of initial data
#' @param model lm object of interest
#' @param clusterby String or formula indicating variable in data for clustering
#' @param x_interest X paramater of interest
#' @param boot_reps Integer for number of replications
#' @return Bootstrap p-value
#'
t_boot_p_val <- function(se, beta, H0, data, model, clusterby, x_interest, boot_reps){

  #Calculate t vector
  t <- (beta - H0)/se

  #Calculate initial t
  se0 <- clustered_se(data = data, model = model, clusterby = clusterby)
  beta0 <- coef(model)
  t0 <- (beta0[x_interest] - H0)/se0[x_interest]

  #Calculate p-value
  prop <- sum(t0 < t)/boot_reps
  p <- 2*min(prop, 1 - prop)

  return(p)

}


#' Calculate sandiwch tensor
#'
#' @param clustervars Vector of cluster indices
#' @param X model matrix
#' @param bread Bread matrix
#' @param x_ind Integer for index of regressor of interest
#' @param E Matrix of residuals
#' @param k Number of regressors
#' @param boot_reps Integer for number of replications
#' @param n Number of observations
#' @param comb_n Number of combined cluster dimensions
#' @return Sandwich tensor
#'
cluster_sandwich <- function(clustervars, X, bread, x_ind, E, k, boot_reps, n, comb_n){

  #Create tensor that is one meat matrix for each bootstrap
  meat_tensor <- cluster_meat(clustervars, X, E, k, boot_reps)

  #Get constant for correction
  G <- length(unique(clustervars))
  const <- G/(G-1) * (n-1)/(n-k)

  #Create tensor of half sandwiches
  half_sandwich <- tensorA::mul.tensor(X = bread, i = 2, Y = meat_tensor, j = 1)

  #Create full tensor of sandwiches
  sandwich_tensor <- tensorA::mul.tensor(X = half_sandwich, i = 2, Y = t(bread), j = 1)

  #Return element for x of interest from each sandwich matrix in tensor
  return((-1)^(comb_n - 1) * const * as.vector(sandwich_tensor[x_ind, , x_ind]))

}

#' Calculate meat tensor
#'
#' @param clustervars Vector of cluster indices
#' @param X model matrix
#' @param k Number of regressors
#' @param boot_reps Integer for number of replications
#' @return Meat tensor for all clusters
#'
cluster_meat <- function(clustervars, X, E, k, boot_reps){

  #Split X and E matrices by cluster indices
  clusters <- split(data.frame(X, E), clustervars)

  #Calculate list of meat tensors, reduce by addition and return
  meat_list <- lapply(clusters, multiply_cluster_meat, boot_reps = boot_reps, k = k)
  return(Reduce('+', meat_list))

}

#' Calculate meat tensor for one cluster
#'
#' @param xe Combined X and residual matrix
#' @param boot_reps Integer for number of replications
#' @param k Number of regressors
#' @return Meat tensor for one cluster
#'
multiply_cluster_meat <- function(xe, boot_reps, k){

  #Split X and residual matrices
  x <- as.matrix(xe[,1:k])
  e <- as.matrix(xe[,1:boot_reps+k])

  #Get number of obs
  n <- nrow(x)

  #Calculate meat tensor and return
  half_meat <- tensorA::mul.tensor(X = x,i = 1, Y = array(e, c(n, 1, boot_reps)), j = 1)
  return(array(tensorA::mul.tensor(X = half_meat, Y = half_meat, by = 3), c(k, k, boot_reps)))

}
