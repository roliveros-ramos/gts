
#' Calculate the area of every square of a grid
#'
#' @param x A grid, linear o curvilinear.
#' @param ... Additional parameters for grid computation.
#'
#' @return A matrix, with the area of each square of the grid.
#' @export
#'
area = function(x, ...) {
  UseMethod("area")
}


# S3 methods --------------------------------------------------------------

#' @export
area.grid = function(x, units=NULL, ...) {
  if(!is.null(x$area) & is.null(units)) return(x$area)
  return(calculate_area(x, units=units))
}

#' @export
area.gts = function(x, units=NULL, ...) {
  return(area(x$grid, units=units))
}


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
    return(out)
  }

  if(is.null(MARGIN)) return(.roll(x=X, .FUN=FUN, lag=lag))

  out = apply(X, MARGIN, FUN=.roll, .FUN=FUN, lag=lag, ...)
  if(MARGIN==1) out = t(out)
  return(out)
}

