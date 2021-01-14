#' The 'gplite' package.
#'
#'
#' @docType package
#' @name gplite-package
#' @aliases gplite
#' @useDynLib gplite, .registration = TRUE
#' @import methods
#' @import Rcpp
#'
#' @description
#'
#' \pkg{gplite} implements some of the most common Gaussian process (GP) models. 
#' The package offers tools for integrating out the latent values analytically 
#' using Laplace or expectation propagation (EP) approximation and for estimating the
#' hyperparameters based on maximizing the (approximate) marginal likelihood or posterior.
#' The package also implements some common sparse approximations for larger datasets.
#'
#' @section Functions:
#'
#' Here's a list of the most important functions:
#' \describe{
#'  \item{\link{gp_init}}{
#'  Set up the GP model.
#'  }
#'  \item{\link{cf}, \link{lik}, \link{method}, \link{approx}}{
#'  Choose the covariance functions, likelihood (observation model), type of the GP
#'  (full or some sparse approximation) and the latent function approximation method
#'  (Laplace, EP).
#'  }
#'  \item{\link{gp_optim}, \link{gp_fit}}{
#'  Optimize the model hyperparameters, or just fit the model with the current
#'   hyperparameter values.
#'  }
#'  \item{\link{gp_pred}, \link{gp_draw}}{
#'  Make predictions with the fitted model. Can also be used before fitting to obtain
#'  prior predictive distribution or draws.
#'  }
#'  \item{\link{gp_loo}, \link{gp_compare}}{
#'  Model assessment and comparison using leave-one-out (LOO) cross-validation.
#'  }
#' }
#'
#'
#'
NULL
