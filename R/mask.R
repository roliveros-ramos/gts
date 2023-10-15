
#' Extract, create or assign ocean/land masks.
#'
#' @param x An object containing a mask.
#' @param ... Additional arguments, currently not used.
#'
#' @return The mask.
#' @export
#'
mask = function(x, ...) {
  UseMethod("mask")
}

#' @param n Number of internal points used to calculate the mask.
#' @param thr Threshold to assign a point to the ocean.
#' @param hires Boolean, use high-resolution coastline?
#' @param ocean Boolean, are valid points in the ocean? Default is TRUE. FALSE for land.
#'
#' @rdname mask
#' @export
mask.grid = function(x, n=2, thr=0.8, hires=FALSE, ocean=TRUE, ...) {
  if(!is.null(x$mask)) return(x$mask)
  m = .create_mask(grid=x, n=n, thr=thr, hires=hires, ocean=ocean, ...)
  return(m$mask)
}

#' @rdname mask
#' @export
mask.gts = function(x, n=2, thr=0.8, hires=FALSE, ocean=TRUE, ...) {
  if(!is.null(x$grid$mask)) return(x$grid$mask)
  return(mask(x$grid, n=n, thr=thr, hires=hires, ocean=ocean, ...))
}

#' @rdname mask
#' @export
mask.array = function(x, ...) {
  x = drop(x)
  if(length(dim(x))!=3) stop("Automatic mask extraction only for 3D arrays.")
  out = apply(x, 1:2, FUN = function(x) !all(is.na(x)))
  out = 0 + out
  return(out)
}

#' @param object The object to add or replace the mask.
#' @param mask A mask.
#'
#' @export
#'
#' @rdname mask
'mask<-' = function(x, value) {
  UseMethod('mask<-')
}

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

