

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


.interpolate = function(x, y, z, xout, yout, method="bilinear", extrap=FALSE, control=list(), ...) {

  # check for arguments: x, y, z, xout, yout.
  # decide input and output types

  use_link = !is.null(control$link)
  if(use_link) {
    trans = gaussian(link=control$link)
    z = suppressWarnings(trans$linkfun(z))
    control$link = NULL
  }

  input  = .check_input(x=x, y=y, z=z)
  output = .check_output(x=xout, y=yout)

  case = sprintf("%s%s", c("r", "i", "i")[input$case], c("r", "i", "i")[output$case])
  FUN = get(sprintf(".%s_%s", method, case), mode="function")
  hasNA = input$hasNA

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

  out = try(FUN(x, y, z, xout, yout, method, extrap, hasNA, control, ...))

  if(inherits(out, "try-error")) {
    out = list(x=xout, y=yout, z=NA*xout)
    warning(out)
    return(out)
  }

  if(output$case["is_gridI"]) {
    dim(out$z) = dim(xout)
    out$x = xout
    out$y = yout
  }

  if(use_link) out$z = trans$linkinv(out$z)

  return(out)

}

# this only return the array z, and has it as first argument to be used with apply.
.interp = function(z, x, y, xout, yout, method="bilinear", extrap=FALSE, control=list(), ...) {
  # to be optimize to reduce duplication of check_input and check_output
  .interpolate(x=x, y=y, z=z, xout=xout, yout=yout, method=method, extrap=extrap, control=control, ...)$z
}

# Methods -----------------------------------------------------------------

#' @export
interpolate.default = .interpolate

# Internal functions: akima -----------------------------------------------

# bilinear interpolation
.bilinear_rr = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  .akima_rr(x=x, y=y, z=z, xout=xout, yout=yout,
            method="linear", extrap=extrap, control=control, ...)
}

.bilinear_ri = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  .akima_ri(x=x, y=y, z=z, xout=xout, yout=yout,
            method="linear", extrap=extrap, control=control, ...)
}

.bilinear_ir = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  .akima_ir(x=x, y=y, z=z, xout=xout, yout=yout,
            method="linear", extrap=extrap, control=control, ...)
}

.bilinear_ii = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  .akima_ii(x=x, y=y, z=z, xout=xout, yout=yout,
            method="linear", extrap=extrap, control=control, ...)
}


# splines
.akima_rr = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  # method = c("linear", "akima")
  # bilinear does not extrapolate (put 0s)
  # bicubic does extrapolate. Could be deactivated (add it).

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

.akima_ri = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  # bilinear does not extrapolate (put 0s)
  # bicubic does extrapolate. Could be desactivated (add it).

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

.akima_ir = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

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
  #              linear = (method == "linear"), extrap = extrap,
  #              input="points", output = "grid",
  #              method = method, duplicate = con$duplicate, dupfun = con$dupfun,
  #              nx = con$nx, ny = con$ny, deltri = con$deltri, h=con$h,
  #              kernel=con$kernel, solver=con$solver, degree=con$degree,
  #              baryweight=con$baryweight, autodegree=con$autodegree,
  #              adtol=con$adtol, smoothpde=con$smoothpde, akimaweight=con$akimaweight,
  #              nweight=con$nweight)

  return(out)

}

.akima_ii = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  linear = (method == "linear")

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

# nearest (regular) interpolation
.nearest_rr = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  x = .getBreaks(x)
  y = .getBreaks(y)

  dat = expand.grid(xout=xout, yout=yout)

  indx = cut(dat$xout, breaks=x, labels=FALSE)
  indy = cut(dat$yout, breaks=y, labels=FALSE)
  ind = cbind(indx, indy)

  out = list(x=xout, y=yout,
             z=matrix(z[ind], nrow=length(xout), ncol=length(yout)))

  return(out)

}

.nearest_ri = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {

  x = .getBreaks(x)
  y = .getBreaks(y)

  indx = cut(as.numeric(xout), breaks=x, labels=FALSE)
  indy = cut(as.numeric(yout), breaks=y, labels=FALSE)
  ind = cbind(indx, indy)

  out = list(x=xout, y=yout,
             z=matrix(z[ind], nrow=nrow(xout), ncol=ncol(yout)))

  return(out)

}

.nearest_ir = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'nearest' is not yet implemented for irregular grids (origin).")
}

.nearest_ii = function(x, y, z, xout, yout, method, extrap=FALSE, hasNA=FALSE, control=list(), ...) {
  stop("Method 'nearest' is not yet implemented for irregular grids (origin).")
}


# Internal functions: interpolate -----------------------------------------

.check_input = function(x, y, z) {

  x = drop(x)
  y = drop(y)

  if(length(x)==1) stop("Degenerate dimension x.")
  if(length(y)==1) stop("Degenerate dimension y.")

  z = drop(z)

  isMatrix = is.matrix(x) & is.matrix(y)

  # Either both inputs are vectors or matrices.
  if(!isMatrix & is.matrix(x)) stop("Input y must be a matrix too.")
  if(!isMatrix & is.matrix(y)) stop("Input x must be a matrix too.")

  # check if matrices can be reduced to vectors!

  isMonot  = monot(x, side=2) & monot(y, side=1)
  sameDim  = identical(dim(x), dim(y)) & identical(dim(x), dim(z))
  sameLen  = (length(x)==length(y)) & (length(x)==length(z))
  consDim  = identical(c(length(x), length(y)), dim(z))

  is_gridR = isTRUE(consDim & isMonot)
  is_gridI = isTRUE(sameDim & isMonot & isMatrix)
  is_point = isTRUE(sameLen & !isMatrix)

  if(!any(is_gridR, is_gridI, is_point)) {
    if(!is.matrix(x) & any(is.na(x))) stop("Input vector 'x' cannot contain NAs.")
    if(!is.matrix(y) & any(is.na(y))) stop("Input vector 'y' cannot contain NAs.")
    # pending validation
    stop("Inputs x, y and z are not consistent.")
  }

  hasNA = any(is.na(x)) | any(is.na(y)) | any(is.na(z))

  if(isTRUE(hasNA)) {

    if(!isTRUE(is_point)) {
      if(!is.matrix(x)) x = matrix(x, ncol=ncol(z), nrow=nrow(z))
      if(!is.matrix(y)) y = matrix(y, ncol=ncol(z), nrow=nrow(z), byrow=TRUE)
    }

    dat = data.frame(x=as.numeric(x), y=as.numeric(y), z=as.numeric(z))
    dat = dat[complete.cases(dat), ]

    out = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point),
               hasNA=TRUE, data=dat)

    return(out)

  }

  out = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point),
             hasNA=FALSE, data=NULL)
  return(out)

}

.check_output = function(x, y) {

  x = drop(x)
  y = drop(y)

  isMatrix = is.matrix(x) & is.matrix(y)

  if(!isMatrix & is.matrix(x)) stop("Input y must be a matrix too.")
  if(!isMatrix & is.matrix(y)) stop("Input x must be a matrix too.")

  isMonot  = monot(x, side=2) & monot(y, side=1)
  sameDim  = identical(dim(x), dim(y))
  sameLen  = (length(x)==length(y))

  is_gridR = isTRUE(isMonot & !isMatrix)
  is_gridI = isTRUE(isMonot &  isMatrix)
  is_point = isTRUE(!isMonot & sameLen & !isMatrix)

  if(!any(is_gridR, is_gridI, is_point)) {
    # pending validation
    stop("Arguments 'xout' and 'yout' are not consistent.")
  }

  output = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point))

  return(output)

}

monot = function(x, side) {
  .monot  = function(x) all(diff(x) > 0)
  .monotM = function(x, side) all(apply(x, side, .monot))
  if(is.matrix(x)) return(.monotM(x, side))
  return(.monot(x))
}

