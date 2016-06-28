#' Takes an lm model and returns named vector of clustered standard errors
#'
#' @param model lm object for which clustered errors sought
#' @param clusterby Vector of cluster indices corresponding to X variables
#' @param G Total number of groups, should be specified for performance
#' @param cluster_ids Vector of unique cluster_ids, should be specified for performance
#' @return Named vector of clustered standard errors
#' @export
clustered_se <- function(model, clusterby, G = NULL, cluster_ids = NULL){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  #If unique clusterids not specified, generate them
  if(is.null(cluster_ids)){
    cluster_ids <- unique(clusterby)
  }

  #If G is not specified, calculate it
  if(is.null(G)){
    G <- length(cluster_ids)
  }

  #Recover model matrix from lm model
  X <- model.matrix(model)

  #Get residuals from lm model
  uhat <- resid(model)

  #Calculate SE and return
  se <- calc_cluster_se(X = X, uhat = uhat, clusterby = clusterby, cluster_ids = cluster_ids, G = G)
  return(se)

}


#' Calculate clustered SEs by sandwich matrix method
#'
#' @param X Vector or matrix of x variables
#' @param uhat Vector of residual values from model
#' @param clusterby Vector of cluster indices
#' @param cluster_ids Vecotr of unique cluster indices
#' @param G Integer for the number of groups
#' @return Named vector of clustered standard errors
#' @examples
#' cluster_se(X = x, uhat = uhat, clusterby = clusterby)
calc_cluster_se <- function(X, uhat, clusterby, cluster_ids, G){

  #Calculate equation constant
  n <- dim(X)[1]
  k <- dim(X)[2]
  const <- G/(G-1) * (n-1)/(n-k)

  #Calculate 'bread' matrix
  bread <- solve(crossprod(X))

  #Calculate 'meat' matrix
  meat_list <- lapply(cluster_ids,function(x) crossprod(t(crossprod(X[clusterby==x,],uhat[clusterby==x]))))
  meat <- meat <- Reduce('+', meat_list)

  # Calculate sandwich matrix
  sandwich <-  const * bread %*% meat %*% bread

  #Take diagonal of sandwich to get SE vector and return it
  se <- sqrt(diag(sandwich))
  return(se)

}
