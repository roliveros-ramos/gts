#' Interpolate bivariate data on grids or point locations
#'
#' `interpolate()` provides a unified interface for interpolating two-dimensional
#' data from regular grids, irregular grids, or point locations to new spatial
#' locations. The default method is `"bilinear"`.
#'
#' The function accepts a spatial field as a matrix, or a higher-dimensional array
#' whose first two dimensions represent space. When `z` has additional trailing
#' dimensions, interpolation is applied independently to each layer and the extra
#' dimensions are preserved in the output.
#'
#' @param x,y Coordinates describing the input data. Supported input formats are:
#'   \itemize{
#'   \item a regular grid, given as increasing numeric vectors, with
#'   `dim(z)[1:2] == c(length(x), length(y))`;
#'   \item an irregular grid, given as numeric matrices with the same dimensions
#'   as `z`;
#'   \item point data, given as numeric vectors of the same length as `z`.
#'   }
#' @param z Numeric values to interpolate. `z` must be a matrix for a single
#'   spatial field, or an array whose first two dimensions match the spatial
#'   dimensions implied by `x` and `y`.
#' @param xout,yout Coordinates of the target locations. Supported output formats
#'   are:
#'   \itemize{
#'   \item a regular grid, given as increasing numeric vectors;
#'   \item an irregular grid, given as numeric matrices of identical dimensions;
#'   \item paired point locations, given as numeric vectors of equal length.
#'   }
#' @param method Interpolation method. Currently supported methods in the shared
#'   implementation are `"bilinear"`, `"nearest"`, and `"akima"`. Not all methods
#'   support all combinations of input and output geometry; see \strong{Details}.
#' @param extrap Logical; if `TRUE`, allow extrapolation outside the input domain
#'   when supported by the selected method. If `FALSE`, unsupported target
#'   locations are returned as `NA`.
#' @param control A named list of method-specific control options. Common entries
#'   include:
#'   \describe{
#'   \item{`mask`}{Optional mask used by bilinear interpolation. Values equal to
#'   `0` are treated as missing.}
#'   \item{`link`}{Optional link function name used to transform `z` before
#'   interpolation and back-transform it afterwards via [stats::gaussian()].}
#'   \item{`linear`}{Logical; for `"akima"` interpolation, request linear rather
#'   than Akima interpolation.}
#'   }
#'   Additional named entries are passed to the underlying Akima routines; see
#'   \strong{Details}.
#' @param ... Additional arguments passed to the underlying interpolation helper.
#'
#' @details
#' Inputs are classified internally as one of three spatial data layouts:
#' regular grids, irregular grids, or point data.
#'
#' Missing values in gridded inputs are handled differently depending on the
#' interpolation method. For Akima-based interpolation, complete cases are
#' extracted and treated as point data. For bilinear and nearest-neighbour
#' interpolation, support is currently limited to regular-grid origins.
#'
#' In the current implementation:
#' \itemize{
#' \item `"bilinear"` is implemented for regular-grid origins and regular-grid or
#' irregular-grid targets;
#' \item `"nearest"` is implemented for regular-grid origins and regular-grid or
#' irregular-grid targets;
#' \item `"akima"` supports regular grids, irregular grids, and point data, and is
#' the method that most clearly supports paired point output;
#' \item bilinear and nearest-neighbour interpolation are not implemented for
#' irregular-grid or point-data origins.
#' }
#'
#' For Akima-based interpolation, `control` may include arguments used by the
#' underlying Akima routines, such as `duplicate`, `dupfun`, `nx`, `ny`,
#' `jitter`, `jitter.iter`, `jitter.random`, and `remove`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`x`, `y`}{The requested output coordinates, returned in the same form as
#'   supplied by `xout` and `yout`.}
#'   \item{`z`}{Interpolated values. For regular-grid output this is typically a
#'   matrix of dimension `length(xout) x length(yout)`. For irregular-grid output,
#'   it has `dim(xout)`. If the input `z` has additional trailing dimensions, they
#'   are preserved after the two spatial dimensions.}
#' }
#'
#' @seealso [akima::interp()], [akima::interpp()]
#'
#' @examples
#' \dontrun{
#' data(volcano)
#' x <- seq(0, 1, length.out = nrow(volcano))
#' y <- seq(0, 1, length.out = ncol(volcano))
#' xout <- seq(0, 1, length.out = 3 * nrow(volcano))
#' yout <- seq(0, 1, length.out = 3 * ncol(volcano))
#'
#' out <- interpolate(x, y, volcano, xout, yout, method = "bilinear")
#' str(out)
#' }
#' @export
interpolate = function(x, y, z, ...) {
  UseMethod("interpolate", z)
}


# Auxiliar functions ------------------------------------------------------

#' Internal helpers for `interpolate()`
#'
#' Developer-oriented helper functions used to apply interpolation methods over
#' one or more layers and to manage optional transformations and error handling.
#'
#' @name interpolate-internal
#' @keywords internal
NULL

.interp = function(z, x, y, xout, yout, iFUN, extrap=FALSE, control=list(), ...) {

  iFUN = match.fun(iFUN)

  use_link = !is.null(control$link)
  if(use_link) {
    trans = gaussian(link=control$link)
    z = suppressWarnings(trans$linkfun(z))
    control$link = NULL
  }

  hasNA = control$hasNA

  out = try(iFUN(x, y, z, xout, yout, extrap, hasNA, control, ...))

  if(inherits(out, "try-error")) {
    warning(out)
    return(NA*xout)
  }

  if(use_link) out$z = trans$linkinv(out$z)

  return(out$z)

}

#' Vectorise interpolation over trailing dimensions of an array
#'
#' `interp()` applies `.interp()` over all layers of `z` beyond the first two
#' spatial dimensions and reconstructs the output array.
#'
#' @param z Input matrix or array.
#' @param x,y Input coordinates.
#' @param xout,yout Output coordinates.
#' @param FUN Interpolation helper function.
#' @param extrap Logical; passed to `.interp()`.
#' @param control Named list of internal control options.
#' @param ... Additional arguments passed to `FUN`.
#'
#' @return An array of interpolated values whose first dimensions correspond to
#'   the requested output geometry and whose remaining dimensions match the
#'   trailing dimensions of `z`.
#' @rdname interpolate-internal
interp = function(z, x, y, xout, yout, FUN, extrap=FALSE, control=list(), ...) {

  ndim = dim(xout)
  if(is.null(ndim) | length(ndim)==1) ndim = c(length(xout), length(yout))
  MARGIN = seq_along(dim(z))[-c(1,2)] # all dimensions, including depth
  if(length(MARGIN)>0) {
    xx = apply(z, MARGIN, FUN=.interp, x=x, y=y, xout=xout, yout=yout, iFUN=FUN,
               extrap=extrap, control=control, ...)
  } else {
    xx = .interp(z, x=x, y=y, xout=xout, yout=yout, iFUN=FUN,
                 extrap=extrap, control=control, ...)
  }
  dim(xx) = c(ndim, dim(z)[-c(1,2)])
  return(xx)
}

# Methods -----------------------------------------------------------------

#' @rdname interpolate
#' @method interpolate default
#' @export
interpolate.default = function(x, y, z, xout, yout, method="bilinear", extrap=FALSE, control=list(), ...) {

  # check for arguments: x, y, z, xout, yout.
  # decide input and output types

  input  = .check_input(x=x, y=y, z=z)
  output = .check_output(x=xout, y=yout)

  case = sprintf("%s%s", c("r", "i", "i")[input$case], c("r", "i", "i")[output$case])
  FUN = get(sprintf(".%s_%s", method, case), mode="function")
  control$hasNA = input$hasNA

  if(output$case["is_gridR"]) {
    dat = expand.grid(xout=xout, yout=yout)
  } else {
    dat = list(xout=xout, yout=yout)
  }

  if(method=="bilinear") {
    control$W = get_bilinear_weights(x=x, y=y, xout=dat$xout, yout=dat$yout, mask=control$mask, extrap=extrap)
  }

  if(method=="nearest") {
    control$ind = get_nearest_index(x=x, y=y, xout=dat$xout, yout=dat$yout)
  }

  if(method=="fill") {
    control$ind = get_nearest_index(x=x, y=y, xout=dat$xout, yout=dat$yout)
    control$W = get_bilinear_weights(x=x, y=y, xout=dat$xout, yout=dat$yout, mask=control$mask, extrap=extrap)
  }

  if(input$case["is_point"]) {
    # we use input$data because it has NAs removed
    if(nrow(input$data) < 3) {
      out = list(x=xout, y=yout, z=NA*xout)
      # if(use_link) out$z = trans$linkinv(out$z)
      warning("Not enough data points to interpolate, returning NAs.")
      return(out)
    }
    x = input$data$x
    y = input$data$y
    z = input$data$z
  }

  # barebone interpolation, plus transformation and error handling
  out = list(x=xout, y=yout)
  out$z = .interp(z=z, x=x, y=y, xout=xout, yout=yout, iFUN=FUN, extrap=extrap, control=control, ...)

  if(output$case["is_gridI"]) {
    dim(out$z) = dim(xout)
  }

  return(out)

}


# Bilinear interpolation --------------------------------------------------

.bilinear_rr = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  W = control$W
  invalid = attr(W, "invalid")

  out = list(x=xout, y=yout,
             z=matrix(as.vector(W %*% as.vector(z)), nrow=length(xout), ncol=length(yout)))

  out$z[invalid] = NA

  return(out)

}

# .bilinear_rr = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
#   .akima_rr(x=x, y=y, z=z, xout=xout, yout=yout,
#             method="linear", extrap=extrap, control=control, ...)
# }

.bilinear_ri = .bilinear_rr

.bilinear_ir = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'bilinear' is not yet implemented for irregular grids (origin).")
  # .akima_ir(x=x, y=y, z=z, xout=xout, yout=yout,
  #           method="linear", extrap=extrap, control=control, ...)
}

.bilinear_ii = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'bilinear' is not yet implemented for irregular grids (origin).")
  # .akima_ii(x=x, y=y, z=z, xout=xout, yout=yout,
  #           method="linear", extrap=extrap, control=control, ...)
}


# Akima bicubic interpolation ---------------------------------------------

# splines
.akima_rr = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  # method = c("linear", "akima")
  # bilinear does not extrapolate (put 0s)
  # bicubic does extrapolate. Could be deactivated (add it).

  if(!is.null(control$linear)) control$linear = FALSE
  method = if(isTRUE(control$linear)) "linear" else "akima"

  if(!hasNA) {

    dat = expand.grid(xout=xout, yout=yout)
    out = .akima_ri(x=x, y=y, z=z, xout=dat$xout, yout=dat$yout,
                    method=method, extrap=extrap, control=control, ...)
    out = list(x=xout, y=yout,
               z=matrix(out$z, nrow=length(xout), ncol=length(yout)))
    return(out)

  } else {

    if(!is.matrix(x)) x = matrix(x, ncol=ncol(z), nrow=nrow(z))
    if(!is.matrix(y)) y = matrix(y, ncol=ncol(z), nrow=nrow(z), byrow=TRUE)
    dat = data.frame(x=as.numeric(x), y=as.numeric(y), z=as.numeric(z))
    dat = dat[complete.cases(dat), ] # remove NAs
    .akima_ir(x=dat$x, y=dat$y, z=dat$z, xout=xout, yout=yout,
              method=method, extrap=extrap, control=control, ...)

  }

}

.akima_ri = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  # bilinear does not extrapolate (put 0s)
  # bicubic does extrapolate. Could be desactivated (add it).

  if(!is.null(control$linear)) control$linear = FALSE
  method = if(isTRUE(control$linear)) "linear" else "akima"

  if(hasNA) {
    linear = (method == "linear")
    if(isTRUE(linear)) {
      out = bilinear(x=x, y=y, z=z, x0=xout, y0=yout)
    } else {
      out = bicubic(x=x, y=y, z=z, x0=xout, y0=yout)
    }
    return(out)

  } else {

    if(!is.matrix(x)) x = matrix(x, ncol=ncol(z), nrow=nrow(z))
    if(!is.matrix(y)) y = matrix(y, ncol=ncol(z), nrow=nrow(z), byrow=TRUE)
    dat = data.frame(x=as.numeric(x), y=as.numeric(y), z=as.numeric(z))
    dat = dat[complete.cases(dat), ] # remove NAs
    .akima_ii(x=dat$x, y=dat$y, z=dat$z, xout=xout, yout=yout,
              method=method, extrap=extrap, control=control, ...)

  }

}

.akima_ir = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  if(!is.null(control$linear)) control$linear = FALSE
  method = if(isTRUE(control$linear)) "linear" else "akima"

  linear = (method == "linear")

  # options for akima::interp
  con = list(duplicate = "error", dupfun = NULL,
             nx = 40, ny = 40, jitter = 10^-12, jitter.iter = 6,
             jitter.random = FALSE, remove = !linear)

  con[names(control)]  = control

  out = akima::interp(x=x, y=y, z=z, xo=xout, yo=yout, linear = linear,
                      extrap=extrap, duplicate=con$duplicate, dupfun=con$dupfun,
                      nx=con$nx, ny=con$ny, jitter=con$jitter, jitter.iter = con$jitter.iter,
                      jitter.random = con$jitter.random, remove = con$remove)

  # options for interp::interp
  # con = list(duplicate = "error", dupfun = NULL,
  #            nx = 40, ny = 40, deltri = "shull", h=0,
  #            kernel="gaussian", solver="QR", degree=3,
  #            baryweight=TRUE, autodegree=FALSE, adtol=0.1,
  #            smoothpde=FALSE, akimaweight=TRUE, nweight=25)
  #
  # con[names(control)]  = control
  # out = interp::interp(x=x, y = y, z=z, xo=xout, yo=yout,
  #              linear = (method == "bilinear"), extrap = extrap,
  #              input="points", output = "grid",
  #              method = method, duplicate = con$duplicate, dupfun = con$dupfun,
  #              nx = con$nx, ny = con$ny, deltri = con$deltri, h=con$h,
  #              kernel=con$kernel, solver=con$solver, degree=con$degree,
  #              baryweight=con$baryweight, autodegree=con$autodegree,
  #              adtol=con$adtol, smoothpde=con$smoothpde, akimaweight=con$akimaweight,
  #              nweight=con$nweight)

  return(out)

}

.akima_ii = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  if(!is.null(control$linear)) control$linear = FALSE
  method = if(isTRUE(control$linear)) "linear" else "akima"

  linear = (method == "bilinear")

  # options for akima::interpp
  con = list(duplicate = "error", dupfun = NULL,
             nx = 40, ny = 40, jitter = 10^-12, jitter.iter = 6,
             jitter.random = FALSE, remove = !linear)

  con[names(control)]  = control

  out = akima::interpp(x=x, y=y, z=z, xo=xout, yo=yout, linear = linear,
                       extrap=extrap, duplicate=con$duplicate, dupfun=con$dupfun,
                       jitter=con$jitter, jitter.iter = con$jitter.iter,
                       jitter.random = con$jitter.random, remove = con$remove)

  return(out)

}


# Nearest neighbour interpolation -----------------------------------------

.nearest_rr = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  if(!is.null(control$ind)) {
    ind = control$ind
  } else {
    x = .getBreaks(x)
    y = .getBreaks(y)

    dat = expand.grid(xout=xout, yout=yout)

    indx = cut(dat$xout, breaks=x, labels=FALSE)
    indy = cut(dat$yout, breaks=y, labels=FALSE)
    ind = cbind(indx, indy)
  }

  out = list(x=xout, y=yout,
             z=matrix(z[ind], nrow=length(xout), ncol=length(yout)))

  return(out)

}

.nearest_ri = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  if(!is.null(control$ind)) {
    ind = control$ind
  } else {
    x = .getBreaks(x)
    y = .getBreaks(y)

    indx = cut(as.numeric(xout), breaks=x, labels=FALSE)
    indy = cut(as.numeric(yout), breaks=y, labels=FALSE)
    ind = cbind(indx, indy)
  }

  out = list(x=xout, y=yout,
             z=matrix(z[ind], nrow=nrow(xout), ncol=ncol(yout)))

  return(out)

}

.nearest_ir = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'nearest' is not yet implemented for irregular grids (origin).")
}

.nearest_ii = function(x, y, z, xout, yout, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'nearest' is not yet implemented for irregular grids (origin).")
}

