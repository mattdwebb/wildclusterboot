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

  #Check to ensure no incorrect bound values supplied
  bound <- if(missing(bound)) 'upper' else match.arg(bound, several.ok = TRUE)

  #Get dataframe from model
  data <- eval(model$call$data)

  #Get wild bootstrap fitted data
  wild_data <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  #Generate matrix of y-wild values
  y_wild <- wild_y(wild_data = wild_data, bootby = bootby, boot_dist = boot_dist, boot_reps = boot_reps, enum = enum)

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

  sandwich <- Reduce('+', sandwich_list)

  if(class(clusterby) == 'formula' & any(sandwich < 0)){

    #Run correction on individual sandwich matrices
    eigen_sandwich <- apply(X = sandwich, MARGIN = 2, FUN = eigen_fix)

    #Reform array and extract SE
    sandwich_array <- array(eigen_sandwich, dim = c(k, k, boot_reps))
    se <- sqrt(sandwich_array[x_ind, x_ind, ])

  } else {

    #Extract SE
    se <- sqrt(as.vector(sandwich[x_ind, ,x_ind]))
  }

  #Get beta for x of interest
  beta <- B[x_ind, ]

  p_value <- t_boot_p_val(se = se, beta = beta, H0 = H0, data = data, model = model,
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
#' @param wild_data Dataframe of original data with wild fitted data and residuals
#' @param bootby string or formula indicating variables in wild_data to group bootstrap weights by
#' @param boot_dist Vector of weights for wild bootstrap, string specifying default distribution or a function that only takes the argument "n"
#' @param boot_reps Number of repititions for resampling data
#' @param enum Boolean indicating whether to calculate all possible wild bootstrap combinations, will override boot_reps
#' @return Matrix of bootstrap weights
#'
wild_y <- function(wild_data, bootby, boot_dist, boot_reps, enum){

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
  boot_unique <- unique(wild_data[bootby])

  boot_reps <- if(enum) 2^nrow(boot_unique) else boot_reps

  # #Add weights for each group
  weights <- gen_boot_weights(boot_dist = boot_dist, boot_unique = boot_unique, boot_reps = boot_reps, enum = enum) %>%
    setNames(paste0('wild_boot_weight_', 1:boot_reps))

  weight_names <- names(weights)
  boot_weights <- cbind(boot_unique, weights)
  expanded_weights <- merge(x = wild_data, y = boot_weights)

  #Generate matrix of y-wild values
  y_wild <- as.matrix(expanded_weights[,weight_names]) * wild_data[,'uhat'] + wild_data[,'fitted_data']

  return(y_wild)

}
