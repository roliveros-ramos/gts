#' Compute the area of grid cells
#'
#' `area()` computes the horizontal area of each cell in a spatial grid.
#'
#' For `grid` objects, the method returns a matrix with the area of each grid
#' cell. For `gts` and `static` objects, the calculation is delegated to the
#' associated grid.
#'
#' When an area matrix is already stored in a `grid` object and no output units
#' are requested, that cached matrix is returned directly.
#'
#' @param x A `grid`, `gts`, or `static` object.
#' @param units Character string giving the output area units. Supported values
#'   are `"km2"`, `"m2"`, and `"cm2"`. If omitted, the current implementation
#'   defaults to `"km2"` when the area is computed.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' `area.grid()` computes cell area from the `psi` grid coordinates stored in the
#' object. The calculation uses the corner geometry defined by `x$psi$LON` and
#' `x$psi$LAT`, together with the latitude-dependent cosine correction used in
#' the current implementation.
#'
#' If `x$area` is already available and `units` is `NULL`, the stored matrix is
#' returned unchanged.
#'
#' If the `psi` grid is missing, the low-level area calculation currently returns
#' `NULL`.
#'
#' `area.gts()` and `area.static()` are convenience methods that apply `area()`
#' to `x$grid`.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`area.grid()`}{A numeric matrix with the area of each grid cell, or
#'   `NULL` when the required `psi` coordinates are not available.}
#'   \item{`area.gts()`}{A numeric matrix with the area of each cell of the
#'   associated grid, or `NULL` when the required `psi` coordinates are not
#'   available.}
#'   \item{`area.static()`}{A numeric matrix with the area of each cell of the
#'   associated grid, or `NULL` when the required `psi` coordinates are not
#'   available.}
#' }
#'
#' @seealso [grid-class], [gts-class], [static-class], [make_grid()], [fill()]
#'
#' @examples
#' \dontrun{
#' a1 <- area(grd)
#' a2 <- area(grd, units = "m2")
#' a3 <- area(x)
#' a4 <- area(bathy)
#' }
#' @name area
NULL

#' @rdname area
#' @export
area = function(x, ...) {
  UseMethod("area")
}

#' @describeIn area Compute the area of each cell of a `grid` object.
#' @export
area.grid = function(x, units=NULL, ...) {
  if(!is.null(x$area) & is.null(units)) return(x$area)
  return(calculate_area(x, units=units))
}

#' @describeIn area Compute the area of the cells of the grid associated with a
#'   `gts` object.
#' @export
area.gts = function(x, units=NULL, ...) {
  return(area(x$grid, units=units))
}

#' @describeIn area Compute the area of the cells of the grid associated with a
#'   `static` object.
#' @export
area.static = area.grid

# Internal functions ------------------------------------------------------

calculate_area = function(x, units=NULL) {

  if(is.null(units)) {
    units = "km2"
    message(sprintf("Area units are %s.", units))
  }

  punits = c(km2=1, m2=1e6, cm2=1e10)
  if(!(units %in% names(punits))) stop(sprintf("Unit '%s' not supported.", units))
  fac = punits[units]

  is_null_psi = is.null(x$psi) | (is.null(x$psi$LAT) & is.null(x$psi$LON))
  if(is_null_psi) return(NULL)

  pLON = x$psi$LON
  pLAT = x$psi$LAT
  dyp = t(diff(t(pLAT)))
  dxp = diff(pLON)
  LAT = x$LAT

  xm = roll(pLON, 2, FUN=min, na.rm=TRUE)
  xM = roll(pLON, 2, FUN=max, na.rm=TRUE)

  ym = roll(pLAT, 1, FUN=min, na.rm=TRUE)
  yM = roll(pLAT, 1, FUN=max, na.rm=TRUE)

  dx = xM[,-1] - xm[, -ncol(xm)]
  dx1 = dxp[, -ncol(dxp)]
  dx2 = dxp[, -1]
  dx1p = dx - dx1
  dx2p = dx - dx2

  dy = yM[-1, ] - ym[-nrow(ym), ]
  dy1 = dyp[-nrow(dyp), ]
  dy2 = dyp[-1, ]
  dy1p = dy - dy1
  dy2p = dy - dy2

  A0 = dx*dy - (dx1p*dy2p + dx2p*dy1p + 0.5*dx2p*dy1 + 0.5*dx1p*dy2)
  A = A0*(111^2)*cospi(LAT/180)

  return(A*fac)

}

roll = function(X, MARGIN=NULL, FUN, lag=1, ...) {

  FUN = match.fun(FUN)

  .roll = function(x, .FUN, lag=1, ...) {
    .FUN = match.fun(.FUN)
    n = length(x) - lag
    out = numeric(n)
    for(i in seq_len(n)) {
      out[i] = .FUN(x[i + seq_len(lag+1) - 1], ...)
    }
    # ..FUN = function(i) .FUN(x[i + seq_len(lag+1) - 1])
    # out = sapply(1:n, FUN=..FUN)
    return(out)
  }

  if(is.null(MARGIN)) return(.roll(x=X, .FUN=FUN, lag=lag))

  out = apply(X, MARGIN, FUN=.roll, .FUN=FUN, lag=lag, ...)
  if(MARGIN==1) out = t(out)
  return(out)
}

