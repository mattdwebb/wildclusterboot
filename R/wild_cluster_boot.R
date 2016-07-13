#' @param data Dataframe with all data, including group indices
#' @param model lm object of interest
#' @param x_interest X paramater of interest
#' @param clusterby String with name of clusterby variable in data
#' @param boot_dist Vector of weights for wild bootstrap or string specifying default distribution
#' @param boot_reps Number of repititions for resampling data
#' @param bootby String with name of bootby variable in data, default is same as clusterby
#' @param H0 Integer inticating the alternative hypothesis, default is 0
#' @param cores Interger indicating cores for multicore processing
#' @return boot.out object from built-in boot function
#' @export
wildclusterboot <- function(data, model, x_interest, clusterby, boot_dist, boot_reps, bootby = clusterby, H0 = 0, cores = 1){

  #Check if model is class lm
  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  #Set parallel option based on number of cores
  if(cores == 1){
    parallel <- 'no'
  } else {
    parallel <- 'snow'
  }

  mle <- wild_clust_mle(model = model, x_interest = x_interest, clusterby = clusterby, boot_dist = boot_dist, bootby = bootby, H0 = H0)

  boot.out <- boot::boot(data = data,
                         statistic = wild_clust_statistic,
                         R = boot_reps,
                         sim = 'parametric',
                         ran.gen = wild_clust_ran,
                         mle = mle,
                         parallel = parallel,
                         ncpus = cores,
                         model = model,
                         x_interest = x_interest,
                         clusterby = clusterby,
                         H0 = H0)

  return(boot.out)

}


#' Generates wild clustered bootstrap list to pass to mle parameter of boot function
#'
#' @param model lm object of interest
#' @param x_interest X paramater of interest
#' @param clusterby Vector of cluster indices corresponding to X variables
#' @param boot_dist Vector of weights for wild bootstrap or string specifying default distribution
#' @param bootby String with name of bootby variable in data
#' @param H0 Integer inticating the alternative hypothesis, default is 0
#' @return List of specifications to pass to mle option of boot function
#'
wild_clust_mle <- function(model, x_interest, boot_dist, bootby, H0 = 0){

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

  #Set out list of mle options
  mle_list <- list(uhat = uhat,
                   fitted_data = fitted_data,
                   y_name = y_name,
                   boot_dist = boot_dist,
                   bootby = bootby)

  return(mle_list)

}

#' Function to curry function for ran.gen parameter of boot function
#'
#' @param data Data of interest, passed by boot function
#' @param mle List of optional parameters, passed by boot function
#' @return Randomized data using wild bootstrap proccess
#'
wild_clust_ran <- function(data, mle){

  mle$data <- data

  return(do.call(what = wild_clust_ran_, args = mle))

}

#' Function to pass to ran.gen parameter of boot function
#'
#' @param data Dataframe, which will be the initial data passed to boot data parameter
#' @param y_name String with name of response variable in data
#' @param uhat Vector of residuals from model
#' @param fitted_data Vector of fitted data from model
#' @param bootby String with name of bootby variable in data
#' @param boot_dist Vector of bootstrap distribution
#' @return Randomized data using wild bootstrap proccess
wild_clust_ran_ <- function(data, y_name, uhat, fitted_data, bootby, boot_dist){

  #Create unique vector of group ids
  boot_unique <- unique(data[bootby])

  #Add weights for each group
  boot_weights <- cbind(boot_unique, weight = sample(x = boot_dist, size = nrow(boot_unique), replace = TRUE))
  expanded_weights <- merge(x = data, y = boot_weights)

  #Combine weights, fitted data and residuals to create new y values
  data[y_name] <- fitted_data + uhat * expanded_weights[, 'weight']

  return(data)

}

#' Function for wild clustered bootstrap to pass to statistic parameter of boot function
#'
#' @param data Dataframe, which will be the data passed from boot ran.gen function at each iteration
#' @param model lm model of interest
#' @param x_interest String with name of the independent variable of interest
#' @param clusterby String with name of clusterby variable in data
#' @param H0 Number indicating the null hypothesis, default is 0
#' @return List containing clustered standard error and coefficients for wild bootstrap

wild_clust_statistic <- function(data, model, x_interest, clusterby, H0){

  form <- formula(model)

  #Estimate model from data
  new_model <- lm(data = data, formula = form)

  #Get clustered SE and beta for variable of interest
  se_list <- lapply(X = clusterby, FUN = function(clust) clustered_se(model = new_model, clusterby = data[,clust])[x_interest])
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
