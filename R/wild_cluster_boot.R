#' Function to curry function for ran.gen parameter of boot function
#'
#' @param model lm model of interest
#' @param x_interest String with name of independent variable of interest
#' @param bootby String with name of bootby variable in data
#' @param boot_dist Vector of bootstrap distribution or string indicating built-in distribution
#' @param H0 Number indicating the null hypothesis, default is 0
#' @return Function to pass to ran.gen parameter of boot function
#' @export
wild_clust_ran <- function(model, x_interest, bootby, boot_dist, H0 = 0){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  #If boot_dist is character, resolve to default distribution
  if(class(boot_dist) == 'character'){

    default_dist <- list(six_pt = c(-sqrt(3/2),-sqrt(2/2),-sqrt(1/2),sqrt(1/2),sqrt(2/2),sqrt(3/2)),
                         two_pt = c(-1,1),
                         norm = rnorm)

    if(!boot_dist %in% names(default_dist)){
      error_message <- paste('Only the following default distributions supported:', paste(names(default_dist), collapse = ', '))
      stop(error_message)
    }

    boot_dist <- default_dist[[boot_dist]]

  }

  #Extract data from the model
  data <- model.frame(model)

  #Set out model variable names
  model_vars <- all.vars(formula(model))
  y_name <- model_vars[1]
  short_vars <- model_vars[!model_vars %in% c(x_interest, y_name)]

  #This will impose the H0
  data[y_name] <- data[y_name] - H0 * data[x_interest]

  #Check if short model is just intercept, otherwise create short formula
  if(length(short_vars) == 0){
    form <- reformulate(termlabels = c('1'), response = y_name)
  } else{
    form <- reformulate(termlabels = short_vars, response = y_name)
  }

  #Create short data frame and estimate short model
  short_model <- lm(data = data, formula = form)

  #Get residuals and fitted data from short model
  uhat <- resid(short_model)
  fitted_data <- fitted(short_model) + H0 * data[, x_interest]

  #Returns a curried version of curry_wild_clust_ran with variables filled in
  return(functional::Curry(FUN = curry_wild_clust_ran,
                           y_name = y_name,
                           uhat = uhat,
                           fitted_data = fitted_data,
                           bootby = bootby,
                           boot_dist = boot_dist))

}

#' Function to pass to ran.gen parameter of boot function
#'
#' @param data Dataframe, which will be the initial data passed to boot data parameter
#' @param y_name String with name of response variable in data
#' @param uhat Vector of residuals from model
#' @param fitted_data Vector of fitted data from model
#' @param bootby String with name of bootby variable in data
#' @param boot_dist Vector of bootstrap distribution
#' @param mle Placeholder variable from boot function
#' @return Matrix of y_wild, X variables and clusterby, to be passed to boot statistic function
curry_wild_clust_ran <- function(data, y_name, uhat, fitted_data, bootby, boot_dist, mle){

  #Create dataframe of sampled weights, distributed by boot group
  boot_unique <- unique(data[bootby])
  boot_unique['weight'] <- sample(x = boot_dist, size = nrow(boot_unique), replace = TRUE)
  data[y_name] <- fitted_data + uhat * merge(x = data, y = boot_unique)['weight']

  return(data)

}

#' Function to curry function to pass to statistic parameter of boot function
#'
#' @param model lm model of interest
#' @param x_interest String with name of the independent variable of interest
#' @param clusterby String or list of strings with clusterby variable name or names
#' @return List containing clustered standard error and coefficients for wild bootstrap
#' @param H0 Number indicating the null hypothesis, default is 0
#' @export
wild_clust_statistic <- function(model, x_interest, clusterby, H0 = 0){

  #Get names for all variables from df
  form <- formula(model)

  return(functional::Curry(FUN = curry_wild_clust_statistic,
                           form = form,
                           x_interest = x_interest,
                           clusterby = clusterby,
                           H0 = H0))
}

#' Function for wild clustered bootstrap to pass to statistic parameter of boot function
#'
#' @param data Dataframe, which will be the data passed from boot ran.gen function at each iteration
#' @param formula Formula object to be used on data
#' @param x_interest String with name of the independent variable of interest
#' @param clusterby Vector of group ids for clustering
#' @param H0 Number indicating the null hypothesis, default is 0
#' @return List containing clustered standard error and coefficients for wild bootstrap

curry_wild_clust_statistic <- function(data, form, x_interest, clusterby, H0){

  #Estimate model from data
  model <- lm(data = data, formula = form)

  #Get clustered SE and beta for variable of interest
  se_list <- lapply(X = clusterby, FUN = function(x) clustered_se(model = model, clusterby = data[,x])[x_interest])
  # ses <- clustered_se(model = model, clusterby = clusterby)

  se <- as.numeric(se_list)
  beta <- as.numeric(coef(model)[x_interest])

  return((beta - H0)/se)

}

#' Function for getting p-values out of boot.out object
#'
#' @param boot.out Bootstrap object returned from the boot function
#' @return Named vector of p-values corresponding to each column of boot.out t-matrix
#' @export
boot_p_val <- function(boot.out){

  t0 <- boot.out$t0
  t_df <- data.frame(boot.out$t)
  boot_reps <- boot.out$R

  prop_vec <- mapply(FUN = function(t, t_vec) sum(t < t_vec)/boot_reps, t = t0, t_vec = t_df)
  p_vec <- sapply(prop_vec, function(d) 2*min(d, 1-d))

  return(p_vec)

}
