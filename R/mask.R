#' Extract, create, or assign spatial masks
#'
#' `mask()` extracts or derives a spatial mask from gridded objects and arrays.
#'
#' For `grid`, `gts`, and `static` objects, the function first returns an
#' existing mask when one is already stored in the object. Otherwise, it tries to
#' derive a mask from the available data or, for `grid` objects, to create one
#' from coastline geometry.
#'
#' `mask<-` assigns a new mask to an object. In the current implementation, a
#' replacement method is provided for `gts` objects.
#'
#' @param x An object from which to extract or derive a mask.
#' @param n Number of internal points per grid cell used when creating a
#'   coastline-based mask for a `grid` object.
#' @param thr Threshold between `0` and `1` used to convert the ocean- or
#'   land-proportion surface to a binary mask when creating a mask from a
#'   `grid` object.
#' @param hires Logical; if `TRUE`, use the high-resolution coastline when
#'   creating a mask from a `grid` object.
#' @param ocean Logical; if `TRUE`, valid cells are interpreted as ocean cells.
#'   If `FALSE`, valid cells are interpreted as land cells.
#' @param value A mask to assign to an object.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' For `grid` objects, `mask.grid()` returns `x$mask` when available. Otherwise,
#' it creates a coastline-based mask using the grid geometry. The returned mask
#' is typically a matrix with `1` for valid cells and `NA` for excluded cells.
#'
#' For `gts` objects, `mask.gts()` returns `x$grid$mask` when available.
#' Otherwise, it attempts to derive a mask from the underlying data array
#' `x$x`. This automatic derivation currently relies on the array method and is
#' therefore intended for three-dimensional arrays.
#'
#' For `static` objects, `mask.static()` returns `x$grid$mask` when available.
#' Otherwise, it drops redundant dimensions and attempts to derive a mask from
#' the data. Automatic derivation is currently supported only for two-dimensional
#' static matrices.
#'
#' For arrays, `mask.array()` drops redundant dimensions and derives a mask by
#' checking whether each spatial cell contains at least one non-missing value
#' across the third dimension. This is currently supported only for
#' three-dimensional arrays.
#'
#' For matrices, `mask.matrix()` returns `1` for non-missing cells and `0` for
#' missing cells.
#'
#' `mask<-.gts()` replaces the mask stored in `x$grid$mask`.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`mask.grid()`}{A spatial mask matrix, usually with `1` for valid
#'   cells and `NA` for excluded cells.}
#'   \item{`mask.gts()`}{The stored grid mask when available, otherwise a mask
#'   derived from the data array, or `NULL` when automatic derivation is not
#'   supported.}
#'   \item{`mask.static()`}{The stored grid mask when available, otherwise a
#'   mask derived from the data, or `NULL` when automatic derivation is not
#'   supported.}
#'   \item{`mask.array()`}{A numeric matrix with `1` for cells containing at
#'   least one non-missing value across the third dimension and `0` otherwise.}
#'   \item{`mask.matrix()`}{A numeric matrix with `1` for non-missing cells and
#'   `0` for missing cells.}
#'   \item{`mask<-.gts()`}{The modified `gts` object.}
#' }
#'
#' @seealso [grid-class], [gts-class], [static-class], [make_grid()],
#'   [update_grid()]
#'
#' @examples
#' \dontrun{
#' m1 <- mask(grd)
#' m2 <- mask(x)
#' m3 <- mask(bathy)
#'
#' mask(x) <- m1
#' }
#' @name mask
NULL

#' @rdname mask
#' @export
mask = function(x, ...) {
  UseMethod("mask")
}

#' @describeIn mask Extract or create a mask from a `grid` object.
#' @export
mask.grid = function(x, n=2, thr=0.8, hires=FALSE, ocean=TRUE, ...) {
  if(!is.null(x$mask)) return(x$mask)
  m = .create_mask(grid=x, n=n, thr=thr, hires=hires, ocean=ocean, ...)
  return(m$mask)
}

#' @describeIn mask Extract or derive a mask from a `gts` object.
#' @export
mask.gts = function(x, n=2, thr=0.8, hires=FALSE, ocean=TRUE, ...) {
  if(!is.null(x$grid$mask)) return(x$grid$mask)
  return(mask(x$x, ...))
}

#' @describeIn mask Extract or derive a mask from a `static` object.
#' @export
mask.static = function(x, n=2, thr=0.8, hires=FALSE, ocean=TRUE, ...) {
  if(!is.null(x$grid$mask)) return(x$grid$mask)
  x = drop(x)
  if(length(dim(x))!=2) {
    warning("Automatic mask extraction only for 2D static matrices.")
    return(NULL)
  }
  return(mask(x$x, ...))
}

#' @describeIn mask Derive a mask from a three-dimensional array.
#' @export
mask.array = function(x, ...) {
  x = drop(x)
  if(length(dim(x))!=3) {
    warning("Automatic mask extraction only for 3D arrays.")
    return(NULL)
  }
  out = apply(x, 1:2, FUN = function(x) !all(is.na(x)))
  out = 0 + out
  return(out)
}

#' @describeIn mask Derive a mask from a matrix.
#' @export
mask.matrix = function(x, ...) {
  x = drop(x)
  out = !(is.na(x))
  out = 0 + out
  return(out)
}

#' @rdname mask
#' @export
'mask<-' = function(x, value) {
  UseMethod('mask<-')
}

#' @describeIn mask Replace the mask stored in a `gts` object.
#' @export
'mask<-.gts' = function(x, value) {
  x$grid$mask = value
  x
}

# if(!isGeneric('mask<-')) {
#   setGeneric('mask<-', function(object, mask, ...)
#     standardGeneric('mask<-'))
# }


# put_mask = function(object, mask, n=2, thr=0.8, hires=FALSE, ...) {
#
#   grid = NULL
#   if(inherits(mask, "grid")) grid = mask
#   if(inherits(mask, "gts"))  grid = mask$grid
#
#   if(!is.null(grid)) {
#     grid = update_grid(grid=grid, thr=thr)
#     mask = grid$mask
#   }
#
#   if(!identical(dim(object$LAT), dim(mask)))
#     stop("Dimension of mask is not compatible.")
#
#   if(!is.null(grid)) object$grid = grid
#
#   if(is.null(grid)) {
#     object$grid$mask = mask
#     object$grid$prob  = NULL
#     object$grid$n     = NA
#     object$grid$hires = FALSE
#   }
#
#   return(object)
#
# }
#
#
# setMethod('mask<-', signature(object='gts', mask='matrix'), put_mask)
# setMethod('mask<-', signature(object='gts', mask='array'), put_mask)
# setMethod('mask<-', signature(object='gts', mask='gts'), put_mask)
# setMethod('mask<-', signature(object='gts', mask='grid'), put_mask)


# Internal ----------------------------------------------------------------

.create_mask = function(lon, lat, dx, dy, n=NULL, thr=NULL, hires=NULL, grid=NULL, ocean=TRUE) {

  if(is.null(n))     n     = 2
  if(is.null(thr))   thr   = 0.8
  if(is.null(hires)) hires = FALSE

  if(!is.null(grid)) {
    dx = median(diff(grid$LON))
    dy = median(t(diff(t(grid$LAT))))
  }

  if(is.null(grid))
    grid = .create_grid(lon=lon, lat=lat, dx=dx, dy=dy)

  off  = .grid_offset(n=n, dx=dx, dy=dy)

  if(is.null(off)) {
    out = grid$df
  } else {
    out = list()
    for(i in seq_len(nrow(off))) {
      tmp = grid$df
      tmp$lon = tmp$lon + off$x[i]
      tmp$lat = tmp$lat + off$y[i]
      out[[i]] = tmp
    }
    out = do.call(rbind, out)
  }

  ind = if(isTRUE(ocean)) {
    is_ocean(lon=out$lon, lat=out$lat, hires=hires)
  } else {
    is_land(lon=out$lon, lat=out$lat, hires=hires)
  }

  gg = array(0 + ind, dim=c(dim(grid$LAT), nrow(off)))
  gg = apply(gg, 1:2, mean)

  mask = 0 + (gg >= thr)
  mask[mask==0] = NA # remove 'land'

  output = list(mask=mask, ocean=gg, n=nrow(off)+1)

  return(output)
}

