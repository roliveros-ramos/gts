
#' @export
integrate = function(f, ...) {
  UseMethod("integrate")
}

#' @export
integrate.default = stats::integrate

#' Vertical integration of a gridded time series.
#'
#' @param f An object of class gts.
#' @param lower Lower limit for vertical integration.
#' @param upper Upper limit for vertical integration.
#' @param by Width of the integration interval, 1 by default.
#' @param ... Additional arguments, currently unused.
#'
#' @return The original object after vertical integration
#' @export
#'
#' @examples integrate()
integrate.gts = function(f, lower, upper=0, by=NULL, subdivisions = NULL, ...) {

  out = vertical_integration(f$x, lower=lower, upper=upper, dim=f$depth,
                             by=by, subdivisions=subdivisions, ...)
  f$x = out
  f$depth = NULL

  return(f)

}


# Auxiliar functions ------------------------------------------------------

#' @export
vertical_integration = function(f, lower, upper=0, dim, by=NULL, subdivisions = NULL, ...) {

  if(length(dim(f))<3) stop("At least three dimensions are needed: third dimension assumed as depth.")

  if(!is.array(lower) & !is.array(upper)) {
    .internal = FALSE
    new_f = f
  } else {
    .internal = TRUE
    lower_array = is.array(lower)
    upper_array = is.array(upper)
    if(lower_array & !upper_array) upper = array(upper, dim=dim(lower))
    if(!lower_array & upper_array) lower = array(lower, dim=dim(upper))
    if(!identical(dim(lower)[1:2], dim(upper)[1:2])) stop("Spatial dimensions of lower and upper do not match.")
    if(!identical(dim(lower), dim(upper))) {
      if(length(dim(lower))<length(dim(upper))) lower = array(lower, dim=dim(upper))
      if(length(dim(lower))>length(dim(upper))) upper = array(upper, dim=dim(lower))
      if(prod(dim(lower))<prod(dim(upper)))  lower = array(lower, dim=dim(upper))
      if(prod(dim(upper))<prod(dim(lower)))  upper = array(upper, dim=dim(lower))
    }
    if(!identical(dim(lower), dim(upper))) stop("Dimensions of lower and upper do not match.")
    ref_dim = dim(f)[-3]
    if(!identical(dim(lower)[1:2], ref_dim[1:2])) stop("Spatial dimensions do not match, check upper and lower.")
    if(length(dim(lower))>length(ref_dim)) stop("Dimensions of lower and upper incompatible with array f.")
    lower = array(lower, dim=ref_dim)
    upper = array(upper, dim=ref_dim)
    new_dim = rep(0, length(dim(f)))
    new_dim[3] = 2
    new_dim = new_dim + dim(f)
    new_f = array(dim=new_dim)
    if(length(new_dim)==3) {
      new_f[,,1] = lower
      new_f[,,2] = upper
      new_f[,,-c(1,2)] = f
    } else {
      new_f[,,1,] = lower
      new_f[,,2,] = upper
      new_f[,,-c(1,2),] = f
    }
  }

  MARGIN = seq_along(dim(f))[-3]
  out = apply(new_f, MARGIN = MARGIN, FUN=.vertical_integration,
              x=dim, lower=lower, upper=upper, by=by, subdivisions=subdivisions, .internal=.internal)
  return(out)

}

#' @export
locate = function(y, x, ref, allover=NA, allbelow=NA) {
  if(all(is.na(y))) return(NA)
  yM = y > ref
  if(all(yM)) return(allover)
  if(all(y<ref)) return(allbelow)
  ym = rle(yM)
  ym$lengths = cumsum(ym$lengths)
  if(!tail(ym$values, 1)) warning("non-monotonicity detected close to the surface.")
  ind = max(which(ym$values)) - 1
  if(ind==0) return(NA) # start over zero, returning NA
  val = ym$lengths[ind]
  xx = x[c(val, val+1)]
  yy = y[c(val, val+1)]
  # fun = splinefun(x=x, y=y)
  # xsim = seq(from=rr[1], to=rr[2], length.out=100)
  # ysim = fun(xsim)
  # linear interpolation
  comb = (ref - yy[1])/(yy[2]-yy[1])
  xref = xx[1]*(1-comb) + xx[2]*comb
  return(xref)
}



# Internal functions ------------------------------------------------------

# .vertical_integration = function(f, lower, upper=0, dim, by=NULL, subdivisions=NULL, .internal=FALSE, ...) {
#
#   if(length(dim(f))<3) stop("At least three dimensions are needed: third dimension assumed as depth.")
#   MARGIN = seq_along(dim(f))[-3]
#   out = apply(f, MARGIN = MARGIN, FUN=..vertical_integration,
#               x=dim, lower=lower, upper=upper, by=by, subdivisions=subdivisions, .internal=.internal)
#   return(out)
#
# }

.vertical_integration = function(y, x, lower, upper, by=NULL, subdivisions=NULL, .internal=FALSE) {

  if(.internal) {
    lower = y[1]
    upper = y[2]
    y = y[-c(1,2)]
  }

  if(any(is.na(lower), is.na(upper))) return(NA)

  xr = range(lower, upper)
  lower = max(xr[1], min(x, na.rm=TRUE))
  upper = min(xr[2], max(x, na.rm=TRUE))

  if(!is.null(subdivisions) & !is.null(by))
    stop("Either 'by' or 'subdivisions' must be specified, too many arguments.")
  if(is.null(subdivisions) & is.null(by)) by = 1
  if(!is.null(subdivisions)) by = ((upper - lower)/(subdivisions - 1))

  if(all(is.na(y))) return(NA)
  xfun = splinefun(x=x, y=y)
  xsim = seq(from=lower+by/2, to=upper, by=by)
  ysim = xfun(xsim)
  out = by*sum(ysim)
  return(out)

}

