.onAttach <- function(...) {
  ver <- utils::packageVersion("gplite")
  packageStartupMessage("This is gplite version ", ver)
}



prepare_inputmat <- function(x, vars=NULL) {
  # ensure x is a matrix, and pick the correct columns
  if (is.null(vars))
    return(as.matrix(x))
  else
    return(as.matrix(x)[,vars,drop=F])
}







is_fixed <- function(object, param_names) {
  # identify fixed hyperparameter(s) of a given object (cf or lik)
  fixed <- sapply(param_names, function(name) {
    class(object$priors[[name]]) == 'prior_fixed'
  })
  return(fixed)
}

filter_fixed <- function(object, param_names) {
  # take only parameter names that are not fixed
  param_names[!is_fixed(object, param_names)]
}

add_obj_name <- function(object, param_names) {
  # attach the object name (cf of lik) to the parameter names
  paste0(get_name(object), '.', param_names)
}

rm_obj_name <- function(object, param_names) {
  # remove the object name from parameter names
  sapply(param_names, function(name) {
    unlist(strsplit(name,'.', fixed=T))[2]
  })
}



