#' Save and load a GP model
#'
#' Convenience functions for saving and loading GP models.
#' 
#' @name gp_saveload
#' 
#' @param gp The gp model object to be saved.
#' @param filename Where to save or load from. 
#'
#' @return \code{gp_load} returns the loaded GP model object.
#'  
#' @examples
#' \donttest{
#' gp_save(gp, 'gp.rda')
#' gp <- gp_load('gp.rda')
#' 
#' }
#'
NULL

#' @rdname gp_saveload
#' @export
gp_save <- function(gp, filename) {
  save(gp, file=filename)
}

#' @rdname gp_saveload
#' @export
gp_load <- function(filename) {
  load(filename)
  gp
}