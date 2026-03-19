#' Integrate functions and gridded time-series objects
#'
#' `integrate()` is an S3 generic.
#'
#' For ordinary functions, the default method delegates to [stats::integrate()].
#' For `gts` objects, `integrate.gts()` performs vertical integration along the
#' third dimension of the data array, which is assumed to represent depth.
#'
#' @param f For `integrate.default()`, a function to integrate. For
#'   `integrate.gts()`, a `gts` object with a depth dimension.
#' @param lower Lower integration bound. For `integrate.gts()`, this may be a
#'   scalar or an array defining spatially varying lower bounds.
#' @param upper Upper integration bound. For `integrate.gts()`, this may be a
#'   scalar or an array defining spatially varying upper bounds. The default is
#'   `0`.
#' @param by Width of the integration interval used by the numerical
#'   approximation in `integrate.gts()`. If both `by` and `subdivisions` are
#'   omitted, `by = 1` is used.
#' @param subdivisions Number of subdivisions used to define the integration
#'   step size in `integrate.gts()`. This is an alternative to `by`.
#' @param ... Additional arguments passed to the method. For
#'   `integrate.default()`, these are passed to [stats::integrate()]. For
#'   `integrate.gts()`, they are passed to [vertical_integration()].
#'
#' @details
#' `integrate.gts()` applies [vertical_integration()] to `f$x` using the depth
#' coordinates stored in `f$depth`. The third dimension of the data array is
#' interpreted as depth and removed by the integration. The returned object
#' keeps the original `gts` structure but replaces `f$x` by the vertically
#' integrated array and sets `f$depth` to `NULL`.
#'
#' Integration is performed numerically after spline interpolation along the
#' depth coordinate.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`integrate.default()`}{The result returned by [stats::integrate()].}
#'   \item{`integrate.gts()`}{A `gts` object after vertical integration, with
#'   the depth dimension removed from `f$x` and `f$depth` set to `NULL`.}
#' }
#'
#' @seealso [vertical_integration()], [locate()], [stats::integrate()]
#'
#' @examples
#' \dontrun{
#' out <- integrate(x, lower = -100, upper = 0)
#' }
#' @name integrate
NULL

#' @rdname integrate
#' @export
integrate = function(f, ...) {
  UseMethod("integrate")
}

#' @describeIn integrate Default integration method for ordinary functions.
#' @export
integrate.default = stats::integrate

#' @describeIn integrate Vertically integrate a `gts` object over its depth
#'   dimension.
#' @export
integrate.gts = function(f, lower, upper=0, by=NULL, subdivisions = NULL, ...) {

  out = vertical_integration(f$x, lower=lower, upper=upper, dim=f$depth,
                             by=by, subdivisions=subdivisions, ...)
  f$x = out
  f$depth = NULL

  return(f)

}

# Auxiliar functions ------------------------------------------------------

#' Vertical integration and level location for gridded profiles
#'
#' `vertical_integration()` performs numerical integration along the third
#' dimension of an array, which is assumed to represent depth or another
#' vertical coordinate.
#'
#' `locate()` finds the position along a reference coordinate where a profile
#' crosses a target value.
#'
#' @param f Numeric array to integrate. It must have at least three dimensions,
#'   with the third dimension interpreted as the vertical coordinate.
#' @param lower Lower integration bound. This may be a scalar or an array of
#'   spatially varying bounds.
#' @param upper Upper integration bound. This may be a scalar or an array of
#'   spatially varying bounds. The default is `0`.
#' @param dim Numeric vector giving the vertical coordinate associated with the
#'   third dimension of `f`.
#' @param by Width of the integration interval. If both `by` and `subdivisions`
#'   are omitted, `by = 1` is used.
#' @param subdivisions Number of subdivisions used to define the integration
#'   step size. This is an alternative to `by`.
#' @param y Numeric profile values for [locate()].
#' @param x Numeric coordinate associated with `y` in [locate()].
#' @param ref Reference value to locate within the profile.
#' @param allover Value returned by [locate()] when all profile values are above
#'   `ref`.
#' @param allbelow Value returned by [locate()] when all profile values are
#'   below `ref`.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' `vertical_integration()` removes the third dimension of `f` by integrating
#' each vertical profile independently. When `lower` and `upper` are scalars,
#' the same bounds are used for every profile. When one or both are arrays,
#' integration bounds may vary spatially.
#'
#' Numerical integration is performed after spline interpolation of each profile
#' over the vertical coordinate `dim`. The integral is approximated from the
#' interpolated profile evaluated at regularly spaced points between the bounds.
#'
#' `locate()` assumes that the profile is sufficiently monotonic near the
#' crossing of interest. It identifies the interval where `y` crosses `ref` and
#' returns the corresponding coordinate by linear interpolation.
#'
#' @return
#' Depending on the function:
#' \describe{
#'   \item{`vertical_integration()`}{An array with the same dimensions as `f`
#'   except that the third dimension is removed.}
#'   \item{`locate()`}{A numeric coordinate giving the interpolated position
#'   where `y` crosses `ref`, or one of `allover`, `allbelow`, or `NA` when no
#'   crossing can be identified.}
#' }
#'
#' @seealso [integrate.gts()], [stats::splinefun()]
#'
#' @examples
#' \dontrun{
#' zint <- vertical_integration(arr, lower = -200, upper = 0, dim = depth)
#' zref <- locate(y = temp_profile, x = depth, ref = 15)
#' }
#' @name vertical_integration
NULL

#' @rdname vertical_integration
#' @export
vertical_integration = function(f, lower, upper=0, dim, by=NULL, subdivisions = NULL, ...) {

  # TODO: check if all metadata for the depth dimension is removed.

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

#' @describeIn vertical_integration Locate where a profile crosses a reference
#'   value.
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

.vertical_integration = function(y, x, lower, upper, by=NULL, subdivisions=NULL, .internal=FALSE, xxx=TRUE) {

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
  out = if(xxx) by*sum(ysim) else mean(ysim)
  return(out)

}

.vertical_interpolation = function(y, x, xout, .internal=FALSE) {

  if(.internal) {
    y = y[-c(1,2)]
  }

  if(all(is.na(y))) return(NA)

  xfun = splinefun(x=x, y=y)
  out = xfun(xout)
  return(out)

}
