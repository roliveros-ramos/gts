#' Regridding
#'
#' @param object The object
#' @param grid The grid
#' @inheritParams interpolate
#' @return A regrided object
#' @export
#'
#' @examples
#' regrid(object, grid)
regrid = function(object, grid, ...) {
  UseMethod("regrid")
}

if(!isGeneric("regrid")) {
  setGeneric("regrid", function(object, grid, ...)
    standardGeneric("regrid"))
}


# Internal functions ------------------------------------------------------

regrid_gts = function(object, grid, method="bilinear", extrap=FALSE, control=list(), ...) {

  if(inherits(grid, "gts")) grid = grid$grid

  MARGIN = seq_along(dim(object$x))[-c(1,2)]
  xx = apply(object$x, MARGIN, .interp, x=object$longitude, y=object$latitude,
             xout=grid$longitude, yout=grid$latitude, method=method, extrap=extrap,
             control=control, ...)
  dim(xx) = c(dim(grid$LAT), dim(object$x)[-c(1,2)])

  # mask correction
  xx = .mask_correction(x=xx, mask=grid$mask)

  object$grid = grid
  object$longitude = object$grid$longitude
  object$latitude = object$grid$latitude

  if(!is.matrix(object$longitude) & !is.matrix(object$latitude)) {
    object$breaks[[1]] = .getBreaks(object$longitude)
    object$breaks[[2]] = .getBreaks(object$latitude)
    object$info$dim[[1]] = object$longitude
    object$info$dim[[2]] = object$latitude
    names(object$info$dim)[1:2] = c("longitude", "latitude")
    object$info$var = "x"
    object$info$dim.units[1:2] = c("degrees East", "degrees North")
    names(object$info$dim.units)[1:2] = c("longitude", "latitude")
    object$info$units = tail(object$info$units, 1)
    object$info$ovarid = tail(object$info$ovarid, 1)
  } else {
    object$breaks[[1]] = NA
    object$breaks[[2]] = NA
    object$info$dim[[1]] = seq_len(nrow(object$longitude))
    object$info$dim[[2]] = seq_len(ncol(object$latitude))
    names(object$info$dim)[1:2] = c("i", "j")
    object$info$var = c("longitude", "latitude", "x")
    object$info$dim.units[1:2] = c("", "")
    names(object$info$dim.units)[1:2] = c("i", "j")
    object$info$units = c("degrees East", "degrees North", tail(object$info$units, 1))
    object$info$ovarid = c("longitude", "latitude", tail(object$info$ovarid, 1))
  }

  object$x = xx

  return(object)

}


# Methods -----------------------------------------------------------------


setMethod('regrid', signature(object='gts', grid='gts'), regrid_gts)
setMethod('regrid', signature(object='gts', grid='grid'), regrid_gts)

