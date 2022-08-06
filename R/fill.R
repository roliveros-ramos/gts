
#' Fill missing values in an array
#'
#' @param x A matrix or a supported array (e.g. a gts::gts object).
#' @return The array without missing values.
#' @export
#'
fill = function(x, ...) {
  UseMethod("fill")
}

#' @param method The interpolation method used, currently only 'akima' available.
#' @param control Additional control argument passed to the interpolation method.
#' @param ... Additional arguments, currently ignored.
#' @rdname fill
#' @export
fill.default = function(x, method, control=list(), ...) {
  if(all(!is.na(x))) return(x)
  x0 = row(x)
  y0 = col(x)
  xi = suppressWarnings(.interp(x=x0, y=y0, z=x, xout=x0, yout=y0, method="akima", extrap=TRUE, control=control, ...))
  x[is.na(x)] = xi[is.na(x)]
  return(x)
}


#' @rdname fill
#' @export
fill.grid = function(x, method, control=list(), ...) {
  # 1. check if LAT, LON have NAs, use inter_grid to use monotonicity.
  x = interp_grid(x)
  # 2. check if psi is NULL, create psi using c(NA, x, NA) and fill.
  if(is.null(x$psi)) {
    pgrid = list(LAT = cbind(NA, rbind(NA, x$LAT, NA), NA),
                 LON = cbind(NA, rbind(NA, x$LON, NA), NA))
    pgrid$LON = fill(pgrid$LON, method=method, control=control, ...)
    pgrid$LAT = fill(pgrid$LAT, method=method, control=control, ...)
    pgrid$LON = roll(roll(pgrid$LON, 2, mean), 1, mean)
    pgrid$LAT = roll(roll(pgrid$LAT, 1, mean), 2, mean)
    x$psi = pgrid
  }
  return(x)
}

#' @rdname fill
#' @export
fill.gts = function(x, method, control=list(), ...) {
  # check if grid need to be filled.
  x$grid = fill(x$grid, method=method, control=control, ...)
  # use interpolate to fill x
  x$x = fill(x$x, method=method, control=control, ...)
  if(!is.null(x$grid$mask))
    x$x = .mask_correction(x=x$x, mask=x$grid$mask)
  return(x)
}

