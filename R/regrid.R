#' Regrid a gridded object to a new spatial grid
#'
#' `regrid()` interpolates the spatial dimensions of a `gts` or `static` object
#' to a target grid while preserving any trailing non-spatial dimensions. The
#' target grid can be supplied directly as a `grid` object, or extracted from a
#' `gts` or `static` object.
#'
#' This is a higher-level wrapper around the package interpolation machinery for
#' classed gridded objects.
#'
#' @param object A `gts` or `static` object containing gridded data to be
#'   regridded.
#' @param grid Target grid definition. This can be a `grid` object, or a `gts`
#'   or `static` object from which the grid definition is extracted.
#' @param method Interpolation method. Supported methods in the current
#'   implementation are `"bilinear"`, `"nearest"`, and `"fill"`. If `NULL`,
#'   `"bilinear"` is used.
#' @param extrap Logical; if `TRUE`, allow extrapolation outside the source grid
#'   when supported by the selected method.
#' @param control A named list of control options. Common entries include:
#'   \describe{
#'   \item{`use_mask`}{Logical; if `TRUE` (default), use the source mask when
#'   available.}
#'   \item{`mask`}{Optional mask to use during interpolation. If omitted and
#'   `use_mask = TRUE`, the mask of `object` is used.}
#'   \item{`verbose`}{Logical; if `TRUE`, report progress messages.}
#'   \item{`maxit`, `reltol`, `abstol`}{Control parameters used by
#'   `method = "fill"` to iteratively reduce masked gaps before the final
#'   bilinear interpolation.}
#'   }
#' @param ... Additional arguments passed to the underlying interpolation
#'   helpers.
#'
#' @details
#' `regrid()` uses the spatial coordinates stored in `object` and interpolates
#' the data component `object$x` to the spatial geometry defined by `grid`.
#'
#' When `grid` is itself a `gts` or `static` object, only its grid definition is
#' used. If needed, longitudes are first converted to the prime meridian
#' convention of the target grid.
#'
#' `method = "fill"` performs an iterative pre-filling step on the source grid by
#' repeatedly applying bilinear interpolation with extrapolation enabled. This is
#' used to reduce masked gaps before carrying out the final bilinear regridding
#' to the target grid.
#'
#' If the target grid has a mask, it is applied to the interpolated result after
#' regridding.
#'
#' Only regular and irregular grids are supported by `regrid()`. Point-data
#' inputs are not supported.
#'
#' The returned object has the same class as `object`, with updated grid,
#' longitude, latitude, breaks, and metadata fields consistent with the target
#' grid.
#'
#' @return An object of the same class as `object`, regridded to the spatial
#'   grid defined by `grid`. The data component `x` is interpolated to the new
#'   grid, and the corresponding spatial metadata are updated accordingly.
#'
#' @seealso [interpolate()]
#' @name regrid
#' @export
regrid = function(object, grid, ...) {
  UseMethod("regrid")
}

if(!isGeneric("regrid")) {
  setGeneric("regrid", function(object, grid, ...)
    standardGeneric("regrid"))
}


# Internal functions ------------------------------------------------------

#' Internal implementation for `regrid()` methods
#'
#' Developer-oriented helper used by the S4 `regrid()` methods for `gts` and
#' `static` objects.
#'
#' @param object A `gts` or `static` object.
#' @param grid A `grid`, `gts`, or `static` object defining the target grid.
#' @param method,extrap,control,... See [regrid()].
#'
#' @return An object of the same class as `object`, regridded to the target
#'   spatial grid.
#' @keywords internal
regrid_gts = function(object, grid, method="bilinear", extrap=FALSE, control=list(), ...) {

  on.exit(gc(verbose=FALSE))

  if(method=="fill") {
    ctrl = control
    if(is.null(ctrl$maxit)) ctrl$maxit = 100
    if(is.null(ctrl$reltol)) ctrl$reltol = 0.05
    if(is.null(ctrl$abstol)) ctrl$abstol = 0.01
  }

  if(is.null(control$verbose)) control$verbose = FALSE

  pm = attr(longitude(grid), "pm")
  if(attr(longitude(object), "pm") != pm) {
    longitude(object) = longitude(object, pm)
  }

  if(inherits(grid, "gts") | inherits(grid, "static")) {
    if(method=="fill")
      mask(grid) = mask(grid) # needed to fill the grid mask
    grid = grid$grid
    gc(verbose=FALSE)
  }

  if(is.null(control$use_mask)) control$use_mask = TRUE

  if(isTRUE(control$use_mask)) {
    if(is.null(control$mask)) control$mask = mask(object)
  }

  ndim = if(!is.null(grid$mask)) dim(grid$mask) else dim(grid$LAT)
  if(is.null(ndim)) ndim = dim(grid$longitude)
  if(is.null(ndim) | length(ndim)==1) ndim = c(length(grid$longitude), length(grid$latitude))
  if(is.null(ndim)) stop("The grid does not contain the required spatial information.")

  if(is.null(method)) method = "bilinear"

  # validation
  input  = .check_input(x=object$longitude, y=object$latitude, z=object$x)
  output = .check_output(x=grid$longitude, y=grid$latitude)

  case = sprintf("%s%s", c("r", "i", "i")[input$case], c("r", "i", "i")[output$case])
  control$hasNA = input$hasNA

  if(output$case["is_gridR"]) {
    dat = expand.grid(xout=grid$longitude, yout=grid$latitude)
  } else {
    dat = list(xout=grid$longitude, yout=grid$latitude)
  }

  if(method=="fill") {

    ctrl$verbose = FALSE
    xgrid = object$grid
    ind = get_nearest_index(x=object$longitude, y=object$latitude,
                            xout=dat$xout, yout=dat$yout)[which(mask(grid)==1), ]
    maxit = ctrl$maxit
    reltol = ctrl$reltol
    om = mask(object)
    n0 = numeric(maxit+1)
    n0[1] = sum(om[ind]==0, na.rm=TRUE)
    atol = n0[1]*ctrl$abstol

    for(i in seq_len(maxit)) {
      if(n0[i]==0) break
      if(isTRUE(control$verbose)) message(n0[i])
      if(isTRUE(control$verbose)) message("Filling Iteration", i)
      object = regrid_gts(object, grid=xgrid, method="bilinear", extrap=TRUE, control=ctrl, ...)
      om = mask(object)
      n0[i+1] = sum(om[ind]==0, na.rm=TRUE)
      if(n0[i+1] > atol) next
      rtol = (n0[i]-n0[i+1])/n0[i]
      if(rtol < reltol) break
    }

    method = "bilinear"
    control$mask = mask(object) # always use mask (update)

  }

  if(method=="bilinear") {
    if(isTRUE(control$verbose)) message("Computing Bilinear weights")
    control$W = get_bilinear_weights(x=object$longitude, y=object$latitude,
                                     xout=dat$xout, yout=dat$yout, mask=control$mask, extrap=extrap)
  }

  if(method=="nearest") {
    control$ind = get_nearest_index(x=object$longitude, y=object$latitude,
                                    xout=dat$xout, yout=dat$yout)
  }

  if(input$case["is_point"]) {
    stop("Only regular and irregular grids are supported by regrid.")
  }
  # validation

  if(isTRUE(control$verbose)) message("Performing Interpolation")
  FUN = get(sprintf(".%s_%s", method, case), mode="function")
  xx = interp(object$x, x=object$longitude, y=object$latitude,
              xout=grid$longitude, yout=grid$latitude, FUN=FUN,
              extrap=extrap, control=control, ...)

  # MARGIN = seq_along(dim(object$x))[-c(1,2)] # all dimensions, including depth
  # if(length(MARGIN)>0) {
  #   xx = apply(object$x, MARGIN, FUN=.interp, x=object$longitude, y=object$latitude,
  #              xout=grid$longitude, yout=grid$latitude, iFUN=FUN, extrap=extrap,
  #              control=control, ...)
  # } else {
  #   xx = .interp(object$x, x=object$longitude, y=object$latitude,
  #                xout=grid$longitude, yout=grid$latitude, iFUN=FUN, extrap=extrap,
  #                control=control, ...)
  # }
  #
  # dim(xx) = c(ndim, dim(object$x)[-c(1,2)])

  # mask correction
  xx = .mask_correction(x=xx, mask=grid$mask)

  object$grid = grid
  object$longitude = object$grid$longitude
  object$latitude = object$grid$latitude

  if(!is.matrix(object$longitude) & !is.matrix(object$latitude)) {
    object$breaks[[1]] = .getBreaks(object$longitude)
    object$breaks[[2]] = .getBreaks(object$latitude)
    object$info$dim[[1]] = object$longitude
    object$info$dim[[2]] = object$latitude
    names(object$info$dim)[1:2] = c("longitude", "latitude")
    object$info$var = "x"
    object$info$dim.units[1:2] = c("degrees East", "degrees North")
    names(object$info$dim.units)[1:2] = c("longitude", "latitude")
    object$info$units = tail(object$info$units, 1)
    object$info$ovarid = tail(object$info$ovarid, 1)
  } else {
    object$breaks[[1]] = NA
    object$breaks[[2]] = NA
    object$info$dim[[1]] = seq_len(nrow(object$longitude))
    object$info$dim[[2]] = seq_len(ncol(object$latitude))
    names(object$info$dim)[1:2] = c("i", "j")
    object$info$var = c("longitude", "latitude", "x")
    object$info$dim.units[1:2] = c("", "")
    names(object$info$dim.units)[1:2] = c("i", "j")
    object$info$units = c("degrees East", "degrees North", tail(object$info$units, 1))
    object$info$ovarid = c("longitude", "latitude", tail(object$info$ovarid, 1))
  }

  object$x = xx

  return(object)

}


# Methods -----------------------------------------------------------------

#' @describeIn regrid Regrid a `gts` object to the grid carried by another
#'   `gts` object.
#' @aliases regrid,gts,gts-method
#' @export
setMethod('regrid', signature(object='gts', grid='gts'), regrid_gts)

#' @describeIn regrid Regrid a `gts` object to the grid carried by a `static`
#'   object.
#' @aliases regrid,gts,static-method
#' @export
setMethod('regrid', signature(object='gts', grid='static'), regrid_gts)

#' @describeIn regrid Regrid a `gts` object to a `grid` object.
#' @aliases regrid,gts,grid-method
#' @export
setMethod('regrid', signature(object='gts', grid='grid'), regrid_gts)

#' @describeIn regrid Regrid a `static` object to the grid carried by a `gts`
#'   object.
#' @aliases regrid,static,gts-method
#' @export
setMethod('regrid', signature(object='static', grid='gts'), regrid_gts)

#' @describeIn regrid Regrid a `static` object to the grid carried by another
#'   `static` object.
#' @aliases regrid,static,static-method
#' @export
setMethod('regrid', signature(object='static', grid='static'), regrid_gts)

#' @describeIn regrid Regrid a `static` object to a `grid` object.
#' @aliases regrid,static,grid-method
#' @export
setMethod('regrid', signature(object='static', grid='grid'), regrid_gts)

