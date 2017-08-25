#' Calculate wild cluster bootstrap
#'
#' This function provides a flexible framework for
#' calculating the wild cluster bootsrap, this function
#' implements the original wild cluster boostrap, due
#' to Cameron, Gelbach and Miller (2008), and a number of
#' modifications, including flexible boostrap distributions
#' and implementing a "multi-way" version.
#'
#' @param data Dataframe with all data, including group indices
#' @param model lm object to perform the boostrap with
#' @param x_interest String indicating the name of the paramater of interest
#' @param clusterby String or formula with name of clusterby variable in data.
#'  Use a string for one-way clustering and a formula to indicate
#'   multiway clustering (e.g. use ~ G + H to a clustering
#'   using the G and H dimensions
#' @param boot_reps Integer indicating number of bootstrap repetitions
#' @param boot_dist Specifies the distribution from which to draw the wild boostrap weights.
#'   This can be one of three things: a numeric vector,
#'   a string specifying default distribution ("rad" or
#'   "six" for now)vor a function that takes only the argument
#'   "n" (e.g. rnorm)
#' @param bootby String with name of bootby variable.
#'   The default for this is same as clusterby, however
#'   a single dimension must be specified if the clusterby
#'   variable is a multi-way formula
#' @param H0 Float or integer inticating the null hypothesis, default is 0
#' @param enum Boolean indicating whether to calculate all possible wild bootstrap combinations.
#'  Only valid if using a vector or default distribution.
#'  If this is set to TRUE, then the boot_reps variable
#'  is ignored and a new enumerated total is calculated
#' @param absval Boolean indicating whether or not to use absolute valued t-statistics
#' @param bound String or vector of string indicating which bootstrap "bound" to use when tie-breaking.
#'   upper and lower indicate using the highest or lowest
#'   value when actual value is tied with bootstrapped values.
#'   mid takes the half-way point between the two. uniform
#'   calcualtes a random interval betwen the two and density
#'   uses a kernel smoothing correction due to Racine-MacKinnon
#'   (2007)
#' @return p-value or vector of p-values corresponding to bootstrap result
#' @importFrom magrittr "%>%"
#' @export
wild_cluster_boot <- function(data, model, x_interest, clusterby, boot_reps, boot_dist = 'rad', bootby = clusterby, H0 = 0, enum = FALSE, absval = FALSE, bound = c('upper', 'lower', 'mid', 'uniform', 'density')){

  if(class(model) != 'lm'){
    stop('Model variable must be lm class')
  }

  bound <- if(missing(bound)) 'upper' else match.arg(bound, several.ok = TRUE)

  data <- if(missing(data)) eval(model$call$data) else data

  data_wild <- wild_data(data = data, model = model, x_interest = x_interest, H0 = H0)

  y_wild <- wild_y(data_wild = data_wild, bootby = bootby, boot_dist = boot_dist, boot_reps = boot_reps, enum = enum)

  X <- model.matrix(model)

  x_ind <- grep(x_interest, colnames(X))

  bread <- bread_cpp(X)
  B <- beta_cpp(X, bread, y_wild)
  E <- y_wild - X %*% B

  se <- wild_se(data_wild, E, X, bread, clusterby, x_ind)

  beta <- B[x_ind, ]

  p_value <- boot_p_val(se = se, beta = beta, H0 = H0, data = data, model = model,
                          clusterby = clusterby, x_interest = x_interest,
                          boot_reps = boot_reps, bound = bound, absval = absval)

  return(p_value)

}

# Helper functions -------------------------------------------------------------

wild_data <- function(data, model, x_interest, H0){

  model_vars <- all.vars(formula(model))
  y_name <- model_vars[1]
  short_vars <- model_vars[model_vars != x_interest]
  short_x <- short_vars[short_vars != y_name]

  h0_mutate <- lazyeval::interp(~ a - c*b, a = as.name(y_name), b = as.name(x_interest), c = H0)

  short_data <- data %>%
    dplyr::mutate_(.dots = setNames(list(h0_mutate), y_name))

  short_formula <- if(length(short_x) == 0){
    reformulate(termlabels = c('1'), response = y_name)
  } else{
    reformulate(termlabels = short_x, response = y_name)
  }

  short_model <- lm(data = short_data, formula = short_formula)

  reverse_h0_mutate <- lazyeval::interp(~ a + c*b, a = as.name(y_name), b = as.name(x_interest), c = H0)
  fitted_mutate <- lazyeval::interp(~ fitted_data + c*b, b = as.name(x_interest), c = H0)

  fitted_data <- short_data %>%
    modelr::add_predictions(short_model) %>%
    modelr::add_residuals(short_model) %>%
    dplyr::rename(uhat = resid,
                  fitted_data = pred) %>%
    dplyr::mutate_(.dots = setNames(list(fitted_mutate), 'fitted_data')) %>%
    dplyr::mutate_(.dots = setNames(list(reverse_h0_mutate), y_name))

  return(fitted_data)

}

wild_y <- function(data_wild, bootby, boot_dist, boot_reps, enum){

  #TODO: Split this out as a separate function
  if(class(boot_dist) == 'character'){

    default_dist <- list(sixpt = c(-sqrt(3/2),-sqrt(2/2),-sqrt(1/2),sqrt(1/2),sqrt(2/2),sqrt(3/2)),
                         rad = c(-1,1))

    if(!boot_dist %in% names(default_dist)){

      error_message <- paste('Only the following default distributions supported:', paste(names(default_dist), collapse = ', '))
      stop(error_message)

    }

    boot_dist <- default_dist[[boot_dist]]

  }

  bootby <- if(class(bootby) == 'formula') all.vars(bootby) else bootby

  boot_unique <- unique(data_wild[bootby])
  boot_paste <- Reduce(paste0, data_wild[bootby])
  boot_ind <- as.numeric(factor(boot_paste, levels = unique(boot_paste)))-1

  boot_reps <- if(enum) 2^nrow(boot_unique) else boot_reps

  weights <- gen_boot_mat(boot_dist = boot_dist, boot_unique = boot_unique, boot_reps = boot_reps, enum = enum)

  y_wild <- y_weights(data_wild[,'fitted_data'], data_wild[,'uhat'], weights, boot_ind)

  return(y_wild)

}

gen_boot_mat <- function(boot_dist, boot_unique, boot_reps, enum){

  weights <- if(class(boot_dist) != 'function' & !enum){

    matrix(sample(x = boot_dist, size = nrow(boot_unique)*boot_reps, replace = TRUE),
                      nrow = nrow(boot_unique), ncol = boot_reps, byrow = FALSE)

  } else if(class(boot_dist) != 'function' & enum){

    t(expand.grid(rep(x = list(boot_dist), times = nrow(boot_unique))))

  } else{

    matrix(boot_dist(n = nrow(boot_unique)*boot_reps),
           nrow = nrow(boot_unique), ncol = boot_reps, byrow = FALSE)

  }

  return(weights)

}

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

crve_ind <- function(data_wild, clusterby){

  clustervars <- Reduce(data_wild[clusterby], f = paste0)
  cluster_ind <- as.numeric(factor(clustervars, levels = unique(clustervars)))-1

}

boot_p_val <- function(se, beta, H0, data, model, clusterby, x_interest, boot_reps, absval, bound){

  t <- if(absval) abs((beta - H0)/se) else (beta - H0)/se

  se0 <- clustered_se(data = data, model = model, clusterby = clusterby)
  beta0 <- coef(model)
  t0 <- if(absval) abs((beta0[x_interest] - H0)/se0[x_interest]) else (beta0[x_interest] - H0)/se0[x_interest]

  p <- sapply(bound, p_eval, t = t, t0 = t0,  boot_reps = boot_reps, absval = absval)

  return(p)

}

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

density_p_val <- function(t, t0, boot_reps){

  h <- mlcv(t)

  gauss <- (t0-t)/h

  p <- sum(pnorm(q = gauss))/boot_reps

  return(p)

}

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
