#' The 'gplite' package.
#'
#'
#' @docType package
#' @name gplite-package
#' @aliases gplite
#' @useDynLib gplite, .registration = TRUE
#' @import methods
#' @import Rcpp 
#' @importFrom rstan sampling optimizing
#' 
#' @description 
#' 
#' \pkg{gplite} provides a convenient interface for fitting some of the most common 
#' Gaussian process (GP) models. The package offers tools for integrating out
#' the latent values analytically using Laplace approximation and then estimating the 
#' hyperparameters based on maximizing the (approximate) marginal likelihood. 
#' It is also possible to run MCMC for the latent values conditioned on the estimated 
#' values for the hyperparameters using Stan. The package also implements some common
#' sparse approximations for larger datasets.
#' 
#' 
#'
#' @references
#' Stan Development Team. RStan: the R interface to Stan. \url{https://mc-stan.org}
#'
NULL
