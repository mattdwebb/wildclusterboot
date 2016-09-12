#' Takes an lm model and returns named vector of clustered standard errors
#'
#' @param data Dataframe containing both model data and clusterby variables
#' @param model lm object for which clustered errors sought
#' @param clusterby Vector of cluster indices corresponding to X variables
#' @param G Total number of groups, should be specified for performance
#' @param cluster_ids Vector of unique cluster_ids, should be specified for performance
#' @return Named vector of clustered standard errors
#' @export
clustered_se <- function(data, model, clusterby, G = NULL, cluster_ids = NULL){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

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

  #If unique clusterids not specified, generate them
  if(is.null(cluster_ids)){
    cluster_ids <- lapply(X = clustervars, FUN = unique)
  }

  #If G is not specified, calculate it
  if(is.null(G)){
    G <- lapply(cluster_ids, length)
  }

  #Recover model matrix from lm model
  X <- model.matrix(model)

  #Calculate 'bread' matrix
  bread <- solve(crossprod(X))

  #Get residuals from lm model
  uhat <- resid(model)

  #TODO: generate sandwich mat for each comb of cluster dimensions, then add subtract and extract SEs

  sandwich_list <- mapply(FUN = calc_sandwich_mat, clustervars = clustervars, cluster_ids = cluster_ids, G = G, comb_n = comb_n,
                          MoreArgs = list(X = X, bread = bread, uhat = uhat), SIMPLIFY = FALSE, USE.NAMES = FALSE)

  #Calculate SE and return
  sandwich <- Reduce('+', sandwich_list)

  #Take diagonal of sandwich to get SE vector and return it
  se <- sqrt(diag(sandwich))

  return(se)

}

#' Calculate clustered SEs by sandwich matrix method
#'
#' @param bread Bread matrix as X'X
#' @param uhat Vector of residual values from model
#' @param clusterby Vector of cluster indices
#' @param cluster_ids Vector of unique cluster indices
#' @param G Integer for the number of groups
#' @return Named vector of clustered standard errors
calc_sandwich_mat <- function(X, bread, uhat, clustervars, cluster_ids, G, comb_n){

  #Calculate equation constant
  n <- dim(X)[1]
  k <- dim(X)[2]
  const <- G/(G-1) * (n-1)/(n-k)

  #Calculate 'meat' matrix
  meat_list <- lapply(cluster_ids, function(x) crossprod(t(crossprod(X[clustervars==x, ,drop = FALSE], uhat[clustervars==x]))))
  meat <- Reduce('+', meat_list)

  # Calculate sandwich matrix
  sandwich <-  (-1)^(comb_n - 1) * const * bread %*% meat %*% bread

  return(sandwich)

}
