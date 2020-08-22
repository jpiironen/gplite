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
#' 
#' # init the model
#' gp <- gp_init()
#' 
#' # fit the model (skipped here)
#' 
#' # save the model
#' gp_save(gp, 'gp.rda')
#' 
#' # load the model and remove the file
#' gp <- gp_load('gp.rda')
#' file.remove('gp.rda')
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
  model <- load(filename)
  get(model)
}