
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

  if(is.null(control$create_psi)) control$create_psi = TRUE

  is_null_psi = is.null(x$psi) | (is.null(x$psi$LAT) & is.null(x$psi$LON))

  if(!isTRUE(control$create_psi) & is_null_psi) return(x)

  # 2. check if psi is NULL, create psi using c(NA, x, NA) and fill.
  if(is.null(x$psi) | is.null(x$psi$LAT) | is.null(x$psi$LON)) {
    if(is.null(x$psi$LAT)) x$psi$LAT = cbind(NA, rbind(NA, x$LAT, NA), NA)
    if(is.null(x$psi$LON)) x$psi$LON = cbind(NA, rbind(NA, x$LON, NA), NA)
    pgrid = x$psi
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

