

#' Interpolation of bivariate data
#'
#' @param x x input
#' @param y y input
#' @param z z input
#' @param xout x output
#' @param yout y output
#' @param method The method
#' @param extrap Do you want to extrapolate outside the convex-hull of your data?
#' @param control Control arguments
#' @param ... Additional arguments
#'
#' @return A list
#' @export
#'
#' @examples
#'data(volcano)
#'x = seq(from=0, to=1, length.out = nrow(volcano))
#'y = seq(from=0, to=1, length.out = ncol(volcano))
#'xout = seq(from=0, to=1, length.out = 3*nrow(volcano))
#'yout = seq(from=0, to=1, length.out = 3*ncol(volcano))
#'zout = interpolate(x, y, volcano, xout, yout)
interpolate = function(x, y, z, ...) {
  UseMethod("interpolate", z)
}


# Auxiliar functions ------------------------------------------------------

# this only return the array z, and has it as first argument to be used with apply.
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

interp = function(z, x, y, xout, yout, FUN, extrap=FALSE, control=list(), ...) {

  ndim = dim(xout)
  if(is.null(ndim) | length(ndim)==1) ndim = c(length(xout), length(yout))

  MARGIN = seq_along(dim(z))[-c(1,2)] # all dimensions, including depth
  if(length(MARGIN)>0) {
    xx = apply(z, MARGIN, FUN=.interp, x=x, y=y, xout=xout, yout=yout, iFUN=FUN,
               extrap=extrap, control=control, ...)
  } else {
    xx = .interp(z, x=x, y=y, xout=xout, yout=xout, iFUN=FUN,
                 extrap=extrap, control=control, ...)
  }
  dim(xx) = c(ndim, dim(z)[-c(1,2)])
  return(xx)
}

# Methods -----------------------------------------------------------------

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

  out = interp(x=x, y=y, z=z, xo=xout, yo=yout, linear = linear,
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

  out = interpp(x=x, y=y, z=z, xo=xout, yo=yout, linear = linear,
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

