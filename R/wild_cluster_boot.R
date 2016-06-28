#' Generates wild clustered bootstrap list to pass to mle parameter of boot function
#'
#' @param model lm object of interest
#' @param param X paramater of interest
#' @param clusterby Vector of cluster indices corresponding to X variables
#' @param boot_dist Vector of weights for wild bootstrap
#' @param bootby Vector of bootstrap group indices corresponding to X variables, default is same as clusterby
#' @param H0 Integer inticating the alternative hypothesis, default is 0
#' @return List of specifications to pass to mle option of boot function
#' @export
wild_clust_mle <- function(model, param, clusterby, boot_dist, bootby = clusterby, H0 = 0){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  #Generate unique bootstrap group ids and calculate total groups
  boot_ids <- unique(bootby)
  bG <- length(boot_ids)

  #Extract data from the model
  X <- data.frame(model.matrix(model))
  Y <- model.frame(model)[!names(model.frame(model)) %in% names(X)]

  #Separating the X variable of interest and other variable
  X_short <- X[names(X) != param]
  X_interest <- as.numeric(X[,param])

  #This will impose the H0
  Y_temp <- Y - H0 * X_interest

  #Generating formula for short regression without X of interest
  part_form <- paste(names(X_short), collapse = ' +')
  full_form <- paste(names(Y_temp), part_form, sep = ' ~ ')

  #Create short data frame and estimate short model
  short_data <- cbind(Y_temp, X_short)
  short_model <- lm(data = short_data, formula = as.formula(full_form))

  #Get residuals from short model
  uhat <- resid(short_model)

  #Get fitted data from model and add in X of interest * H_0
  fitted_data <- fitted(short_model) + H0 * X_interest

  #Set out list of mle options
  mle_list <- list(uhat = uhat,
                   fitted_data = fitted_data,
                   y_name = names(Y),
                   param = param,
                   boot_dist = boot_dist,
                   clusterby = clusterby,
                   bootby = bootby,
                   bG = bG,
                   boot_ids = boot_ids)

  return(mle_list)

}

#' Function for wild clustered bootstrap randomization to pass to ran.gen parameter of boot function
#'
#' @param data Dataframe, which will be the initial data passed to boot data parameter
#' @param mle List of wild bootstrap options, generated from gen_wild_mle
#' @return Matrix of y_wild, X variables and clusterby, to be passed to boot statistic function
#' @export
wild_clust_ran <- function(data, mle){

  #Unpack mle list
  uhat <- mle[['uhat']]
  fitted_data <- mle[['fitted_data']]
  y_name <- mle[['y_name']]
  param <- mle[['param']]
  boot_dist <- mle[['boot_dist']]
  clusterby <- mle[['clusterby']]
  bootby <- mle[['bootby']]
  bG <- mle[['bG']]
  boot_ids <- mle[['boot_ids']]

  #Create dataframe of sampled weights, distributed by boot group
  bootby <- data.frame(bootby)
  g_weights <- cbind(bootby = boot_ids, weight = sample(x = boot_dist, size = bG, replace = TRUE))
  u_weights <- merge(x = bootby, y = g_weights, by = 'bootby')

  #Create y_wild as sum of fitted plus residual times sampled weight
  y_wild <- fitted_data + uhat * u_weights[, 'weight']

  #Seperate X data from Y and clusterby
  X <- data[, -c(1, ncol(data))]

  #Check if only one X
  if(is.null(dim(X))){
    X_final <- X
  }
  #Otherwise re-order X so X of interest first
  else{
    X_short <- X[names(X) != param]
    X_interest <- as.numeric(X[, param])
    X_final <- cbind(X_interest, X_short)
  }

  #Return matrix of y_wild, X_data,
  boot_data <- cbind(y_wild, X_final, clusterby)
  return(boot_data)

}

#' Function for wild clustered bootstrap to pass to statistic parameter of boot function
#'
#' @param data Dataframe, which will be the data passed from boot ran.gen function at each iteration
#' @return List containing clustered standard error and coefficients for wild bootstrap
#' @note Dataframe must be in order dependent-independent of interest-other independent-cluster
#' @export
wild_clust_statistic <- function(data){

  #Get names for all variables from df
  y_param <- colnames(data)[1]
  x_interest <- colnames(data)[2]
  x_names <- colnames(data)[-c(1,ncol(data))]
  model_data <- data.frame(data[, -ncol(data)])
  clusterby <- data[, ncol(data)]

  #Create formula for regression
  part_form <- paste(x_names, collapse = ' +')
  full_form <- paste(y_param, part_form, sep = ' ~ ')

  #Estimate model from data
  model <- lm(data = model_data, formula = as.formula(full_form))

  #Get clustered SE and beta for variable of interest
  se <- as.numeric(clustered_se(model = model, clusterby = clusterby)[x_interest])
  beta <- as.numeric(coef(model)[x_interest])

  return(c(se = se, beta = beta))

}
