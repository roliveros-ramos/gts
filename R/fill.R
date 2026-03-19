#' Fill missing values in gridded objects
#'
#' `fill()` replaces missing values in matrices, arrays, and supported gridded
#' package classes by interpolating from neighbouring values.
#'
#' For plain matrices, interpolation is performed on the matrix index space using
#' row and column positions. For higher-dimensional arrays, filling is applied
#' independently to each layer beyond the first two dimensions. For `grid` and
#' `gts` objects, the method also updates associated spatial grid components.
#'
#' @param x An object to fill. Supported inputs include:
#'   \itemize{
#'   \item a numeric matrix;
#'   \item a numeric array whose first two dimensions are spatial;
#'   \item a `grid` object;
#'   \item a `gts` object.
#'   }
#' @param method Interpolation method used to fill missing values. In the current
#'   implementation, only `"akima"` should be considered supported.
#' @param control A named list of control options. Supported entries include:
#'   \describe{
#'   \item{`mask`}{Optional mask used to identify cells that should not trigger
#'   filling. If all non-masked cells are already complete, the input is returned
#'   unchanged.}
#'   \item{`keep_range`}{Logical; if `TRUE`, interpolated values outside the
#'   observed range of the original data are discarded and remain `NA`.}
#'   \item{`create_psi`}{Logical; for `grid` objects, controls whether a missing
#'   `psi` grid should be created. By default, this is enabled for regular grids.}
#'   \item{`padding`}{Currently accepted but not used by the shared
#'   implementation.}
#'   }
#' @param ... Additional arguments passed to the underlying interpolation
#'   machinery.
#'
#' @details
#' The default method fills missing values in a matrix by interpolating over the
#' matrix row and column indices. Existing non-missing values are left unchanged.
#'
#' For arrays, `fill()` is applied independently over all trailing dimensions,
#' assuming that the first two dimensions represent space.
#'
#' For `grid` objects, the function first calls `interp_grid()` and then ensures
#' that the `psi` grid is available. For regular grids, `psi` is constructed by
#' extending the longitude and latitude boundaries and averaging adjacent cell
#' corners. For irregular grids, missing `psi` coordinates are padded, filled,
#' and then smoothed with `roll()`.
#'
#' For `gts` objects, the function first fills the associated grid, then fills
#' the data array stored in `x$x`. If a grid mask is available, it is used during
#' filling and re-applied afterwards with `.mask_correction()`.
#'
#' @return
#' An object of the same general type as `x`:
#' \describe{
#'   \item{matrix input}{A matrix with missing values replaced where possible.}
#'   \item{array input}{An array of the same dimensions as `x`, with filling
#'   applied independently to each spatial layer.}
#'   \item{`grid` input}{A `grid` object, potentially with a completed `psi`
#'   component.}
#'   \item{`gts` input}{A `gts` object with its grid and data component filled,
#'   subject to the grid mask when present.}
#' }
#'
#' @seealso [interpolate()]
#'
#' @examples
#' \dontrun{
#' x <- volcano
#' x[10:15, 12:18] <- NA
#'
#' xf <- fill(x)
#' }
#' @export
fill = function(x, ...) {
  UseMethod("fill")
}

#' @rdname fill
#' @export
fill.default = function(x, method="akima", control=list(), ...) {

  if(all(!is.na(x))) return(x)
  if(!is.null(control$mask)) {
    if(all(!is.na(x[!is.na(control$mask)]))) return(x)
  }

  if(is.null(control$keep_range)) control$keep_range = FALSE
  if(is.null(control$padding)) control$padding=0

  xlim = range(x, na.rm=TRUE)

  x0 = row(x)
  y0 = col(x)

  # if(!is.null(control$link)) {
  #   trans = gaussian(link=control$link)
  #   xi = suppressWarnings(.interp(x=x0, y=y0, z=trans$linkfun(x), xout=x0, yout=y0, method=method, extrap=TRUE, control=control, ...))
  #   xi = trans$linkinv(xi)
  # } else {
  #   xi = suppressWarnings(.interp(x=x0, y=y0, z=x, xout=x0, yout=y0, method=method, extrap=TRUE, control=control, ...))
  # }

  xi = suppressWarnings(.interp(x=x0, y=y0, z=x, xout=x0, yout=y0, method=method, extrap=TRUE, control=control, ...))

  if(isTRUE(control$keep_range)) {
    xi[(xi < xlim[1]) | (xi > xlim[2])] = NA
  }

  x[is.na(x)] = xi[is.na(x)]

  return(x)
}

#' @rdname fill
#' @export
fill.array = function(x, method="akima", control=list(), ...) {
  if(inherits(x, "matrix")) return(fill.default(x=x, method=method, control=control, ...))
  out = apply(x, MARGIN=-(1:2), FUN=fill, method=method, control=control, ...)
  dim(out) = dim(x)
  return(out)
}

#' @rdname fill
#' @export
fill.grid = function(x, method="akima", control=list(), ...) {
  # 1. check if LAT, LON have NAs, use inter_grid to use monotonicity.
  x = interp_grid(x)

  is_regular = .is_regular_grid(x)
  if(is.null(control$create_psi)) control$create_psi = is_regular

  is_null_psi = is.null(x$psi) | (is.null(x$psi$LAT) & is.null(x$psi$LON))

  if(!isTRUE(control$create_psi) & is_null_psi) return(x)

  # 2. check if psi is NULL, create psi using c(NA, x, NA) and fill.
  if(is.null(x$psi) | is.null(x$psi$LAT) | is.null(x$psi$LON)) {
    if(is_regular) {

      xlat = x$LAT[1, ]
      xlon = x$LON[, 1]

      nlat = length(xlat)
      nlon = length(xlon)

      xlat = c(2*xlat[1] - xlat[2], xlat, 2*xlat[nlat] - xlat[nlat-1])
      xlon = c(2*xlon[1] - xlon[2], xlon, 2*xlon[nlon] - xlon[nlon-1])

      xlat = matrix(xlat, nrow=nrow(x$LAT)+2, ncol=ncol(x$LAT)+2, byrow=TRUE)
      xlon = matrix(xlon, nrow=nrow(x$LON)+2, ncol=ncol(x$LON)+2, byrow=FALSE)

      pgrid = list(LON=xlon, LAT=xlat)
      pgrid$LON = roll(roll(pgrid$LON, 2, mean), 1, mean)
      pgrid$LAT = roll(roll(pgrid$LAT, 1, mean), 2, mean)
      x$psi = pgrid

    } else {

      if(is.null(x$psi$LAT)) x$psi$LAT = cbind(NA, rbind(NA, x$LAT, NA), NA)
      if(is.null(x$psi$LON)) x$psi$LON = cbind(NA, rbind(NA, x$LON, NA), NA)
      pgrid = x$psi
      pgrid$LON = fill(pgrid$LON, method=method, control=control, ...)
      pgrid$LAT = fill(pgrid$LAT, method=method, control=control, ...)
      pgrid$LON = roll(roll(pgrid$LON, 2, mean), 1, mean)
      pgrid$LAT = roll(roll(pgrid$LAT, 1, mean), 2, mean)
      x$psi = pgrid

    }
  }
  return(x)
}

#' @rdname fill
#' @export
fill.gts = function(x, method="akima", control=list(), ...) {
  # check if grid need to be filled.
  x$grid = fill(x$grid, method=method, control=control, ...)
  # use interpolate to fill x
  if(is.null(control$mask)) control$mask = x$grid$mask
  x$x = fill(x$x, method=method, control=control, ...)
  if(!is.null(x$grid$mask))
    x$x = .mask_correction(x=x$x, mask=x$grid$mask)
  return(x)
}

