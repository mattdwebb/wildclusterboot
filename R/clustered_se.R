#' Calculate clustered standard errors for lm model
#'
#' This is a flexible function that will calculate cluster-robust
#' standard errors for an `lm` object. This includes the standard,
#' one-way clustered standard error (similar to Stata's `cluster` option)
#' and a "multi-way" clustered standard error, due to Cameron, Gelbach
#' and Miller (2011)
#'
#' @param data Dataframe object containing both model data and clusterby variables
#' @param model lm object to calculate clustered standard errors for
#' @param clusterby String or fomula indicating the column in `data` to cluster on
#'   Use a string for one-way clustering. a formula to indicate
#'   multiway clustering (e.g. use `~ G + H` to a clustering
#'   using the G and H dimensions
#' @return Named vector of clustered standard errors
#' @export
clustered_se <- function(data, model, clusterby){

  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  if(class(clusterby) == 'formula'){

    form_vars <- all.vars(clusterby)

    comb_list <- lapply(X = 1:length(form_vars), FUN = combn, x = form_vars, simplify = FALSE)
    comb_vars <- unlist(comb_list, recursive = FALSE)
    comb_names <- lapply(X = comb_vars, FUN = paste0, collapse = '')
    comb_n <- lapply(X = comb_vars, length)

    data_combs <- lapply(X = comb_vars, function(comb) data[comb])

    clustervars <- lapply(X = data_combs, FUN = Reduce, f = paste0)
    names(clustervars) <- comb_names

  } else{

    clustervars <- data[clusterby]
    comb_n <- 1

  }

  cluster_ids <- lapply(X = clustervars, FUN = unique)

  G <- lapply(cluster_ids, length)

  X <- model.matrix(model)


  bread <- solve(crossprod(X))

  uhat <- resid(model)

  sandwich_list <- mapply(FUN = calc_sandwich_mat, clustervars = clustervars, cluster_ids = cluster_ids, G = G, comb_n = comb_n,
                          MoreArgs = list(X = X, bread = bread, uhat = uhat), SIMPLIFY = FALSE, USE.NAMES = FALSE)

  sandwich <- Reduce('+', sandwich_list)

  #Spectral decomposition correction
  if(class(clusterby) == 'formula' | any(sandwich < 0)){

    cnames <- colnames(sandwich)
    rnames <- rownames(sandwich)

    eig <- eigen(sandwich)

    eig_vals <- ifelse(eig$values < 0, 0, eig$values)

    sandwich <- eig$vectors %*% diag(eig_vals) %*% t(eig$vectors)

    rownames(sandwich) <- rnames
    colnames(sandwich) <- cnames

  }

  se <- sqrt(diag(sandwich))

  return(se)

}

# Helper functions -------------------------------------------------------------

calc_sandwich_mat <- function(X, bread, uhat, clustervars, cluster_ids, G, comb_n){

  n <- dim(X)[1]
  k <- dim(X)[2]
  const <- G/(G-1) * (n-1)/(n-k)

  meat_list <- mapply(FUN = x_u_cross, id = cluster_ids, MoreArgs = list(X = X, uhat = uhat, clustervars = clustervars), SIMPLIFY = FALSE)
  meat <- Reduce('+', meat_list)

  sandwich <-  (-1)^(comb_n - 1) * const * bread %*% meat %*% bread

  return(sandwich)

}

x_u_cross <- function(id, X, uhat, clustervars){

  half_mat <- crossprod(X[clustervars==id, ,drop = FALSE], uhat[clustervars==id])

  full_mat <- crossprod(t(half_mat))
  return(full_mat)

}
