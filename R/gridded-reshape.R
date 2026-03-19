#' Combine and reshape gridded objects
#'
#' These methods combine `gts` objects along their time dimension, convert
#' gridded objects to a long tabular representation, or drop redundant array
#' dimensions.
#'
#' The page covers two related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension.
#' }
#'
#' `c.gts()` concatenates several compatible `gts` objects along their time
#' dimension and rebuilds a common time axis covering the full time span of the
#' inputs.
#'
#' `melt.gts()` and `melt.static()` convert gridded objects to long
#' `data.frame`s, with one row per spatial or spatiotemporal observation.
#'
#' `drop.gts()` and `drop.static()` simplify the underlying data array by
#' removing redundant dimensions.
#'
#' @param ... For `c.gts()`, `gts` objects to combine. For `melt.gts()`,
#'   additional variables to append after being stripped to match the target
#'   object. Additional arguments are currently ignored by `melt.static()` and
#'   the `drop()` methods.
#' @param data A `gts` or `static` object to reshape to long format.
#' @param na.rm Currently accepted by `melt.gts()` and `melt.static()`, but not
#'   used by the current implementation.
#' @param value.name Optional name for the value column created by
#'   `melt.gts()`. For `melt.static()`, this argument is currently accepted but
#'   ignored; the value column is named from `names(data)[1]`.
#' @param short Logical; if `TRUE`, rename `longitude` and `latitude` columns to
#'   `lon` and `lat` in the output of `melt.gts()`. Not used by
#'   `melt.static()`.
#' @param x An object.
#'
#' @details
#' `c.gts()` currently supports only three-dimensional `gts` objects. All inputs
#' must have:
#' \itemize{
#'   \item identical horizontal dimensions;
#'   \item identical grid coordinates;
#'   \item identical time frequency.
#' }
#'
#' The method determines the earliest start and latest end among the inputs,
#' constructs a complete time index covering that range, and inserts each input
#' into the corresponding time window. Time slots not covered by any input are
#' left as `NA`.
#'
#' `melt.gts()` uses the flattened grid stored in `data$grid$df` and expands it
#' over the time dimension. For objects with a depth dimension, depth is also
#' expanded explicitly before time is replicated. The output time column is named
#' according to `data$info$time_var` when available, and `"time"` otherwise.
#'
#' `melt.static()` produces a long `data.frame` over space only. If a depth
#' component is present, depth is included as an additional column.
#'
#' Additional arguments supplied to `melt.gts()` are passed through an internal
#' stripping helper and appended as extra columns. This is intended for companion
#' variables aligned with the target object.
#'
#' `drop.gts()` is intentionally narrow in scope in the current implementation.
#' It returns the object unchanged when it already has three dimensions. For
#' higher-dimensional objects, it only drops the third dimension when its length
#' is exactly one, and then removes the `depth` component.
#'
#' `drop.static()` drops redundant dimensions from the underlying data array, but
#' does not currently update other metadata or grid components.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`c.gts()`}{A `gts` object with a reconstructed time axis and data
#'   array spanning the union of the input time ranges.}
#'   \item{`melt.gts()`}{A long `data.frame` containing longitude, latitude,
#'   optional depth, a time-like variable, the data values, and any appended
#'   stripped variables.}
#'   \item{`melt.static()`}{A long `data.frame` containing longitude, latitude,
#'   optional depth, and the data values.}
#'   \item{`drop.default()`}{The result of [base::drop()].}
#'   \item{`drop.gts()`}{A `gts` object with a singleton third dimension removed
#'   when applicable; otherwise the original object is returned.}
#'   \item{`drop.static()`}{A `static` object with `x$x` simplified using
#'   [base::drop()].}
#' }
#'
#' @seealso [gridded-summary], [quantile.gts()], [gridded-plot], [gts()], [read_gts()], [read_static()]
#'
#' @examples
#' \dontrun{
#' xy <- c(x1, x2)
#' tab1 <- melt(x)
#' tab2 <- melt(bathy)
#' x2 <- drop(x)
#' s2 <- drop(bathy)
#' }
#' @name gridded-reshape
NULL

#' @describeIn gridded-reshape Concatenate compatible `gts` objects along the
#'   time dimension.
#' @export
c.gts = function(...) {

  out = list(...)

  .identical = function(x) all(sapply(x, identical, y=x[[1]]))

  # check time dimension, only 3D for now, extend.
  time_check = all(sapply(out, FUN=function(x) length(dim(x))) == 3)
  if(!time_check) stop("Only 3D gts objects are supported.")
  # check horizontal dimensions are compatible
  dims = lapply(lapply(out, dim), FUN="[", i=1:2)
  dim_check = .identical(dims)
  if(!dim_check) stop("Spatial dimensions are not compatible.")
  # check for identical grids.
  gLON = lapply(out, FUN=function(x) c(x$grid$LON))
  gLAT = lapply(out, FUN=function(x) c(x$grid$LAT))
  coord_check = .identical(gLON) & .identical(gLAT)
  if(!coord_check) stop("Grid coordinates do not match.")
  # check for time frequency
  freqs = sapply(out, frequency)
  freq_check = all(sapply(freqs, identical, freqs[1]))
  if(!freq_check) stop("Frequencies are not compatible.")
  frequency = freqs[1]

  starts = sapply(out, start.gts)
  starts[2,] = starts[2,]/frequency
  starts = colSums(starts)
  istart = which.min(starts)

  ends = sapply(out, end.gts)
  ends[2,] = ends[2,]/frequency
  ends = colSums(ends)
  iend = which.max(ends)

  start = start(out[[istart]])
  end   = end(out[[iend]])

  the_ts = ts(NA, start=start, end=end, frequency = frequency)
  the_ts = ts(seq_along(the_ts), start=start, end=end, frequency = frequency)

  the_array = array(NA, dim=c(dims[[1]], length(the_ts)))
  the_time  = array(NA, dim=c(length(the_ts)))
  info_time  = array(NA, dim=c(length(the_ts)))

  starts = sapply(out, start.gts)
  ends = sapply(out, end.gts)

  for(i in seq_along(out)) {
    ind = as.numeric(window(the_ts, start=starts[,i], end=ends[,i]))
    the_array[,, ind] = out[[i]]$x
    the_time[ind] = out[[i]]$time
    info_time[ind] = out[[i]]$info$time$time
  }
  class(the_time) = class(out[[1]]$time)

  x = out[[1]]
  x$x = the_array
  x$time = the_time
  x$info$time$time = info_time
  x$breaks$time = .getBreaks(x$info$time$time)
  x$info$dim$time = info_time
  x$info$ts = ts(seq_along(x$time), start=start(the_ts), frequency=frequency(the_ts))

  # class(x) = c("gts", "ts")
  class(x) = "gts"
  return(x)

}

#' @describeIn gridded-reshape Convert a `gts` object to long format.
#' @export
melt.gts = function(data, ..., na.rm=FALSE, value.name=NULL, short=FALSE) {

  time_var = data$info$time_var
  if(is.null(time_var)) time_var = "time"
  if(is.null(value.name)) value.name = names(data)[1]

  if(is.null(data$depth)) {
    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     time=rep(data$time, each=nrow(data$grid$df)),
                     value=as.numeric(data$x))

    names(out)[3:4] = c(time_var, value.name)
    if(isTRUE(short)) names(out)[1:2] = c("lon", "lat")

  } else {

    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     depth = rep(data$depth, each=nrow(data$grid$df)))
    out = data.frame(longitude=out$longitude, latitude=out$latitude,
                     depth = out$depth,
                     time=rep(data$time, each=nrow(out)),
                     value=as.numeric(data$x))

    names(out)[4:5] = c(time_var, value.name)
    if(isTRUE(short)) names(out)[1:2] = c("lon", "lat")

  }

  other = list(...)
  if(length(other)>0) {
    other = lapply(other, FUN=.strip, target=data)
    other = as.data.frame(other)
    out = cbind(out, other)
  }
  return(out)

}

#' @describeIn gridded-reshape Convert a `static` object to long format.
#' @export
melt.static = function(data, ..., na.rm=FALSE, value.name=NULL) {

  if(is.null(value.name)) value.name = names(data)[1]

  if(is.null(data$depth)) {
    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     value=as.numeric(data$x))

    names(out)[3] = value.name

  } else {

    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     depth = rep(data$depth, each=nrow(data$grid$df)))
    out$value = as.numeric(data$x)

    names(out)[4] = value.name

  }

  return(out)

}

#' Drop redundant dimensions
#'
#' `drop()` is an S3 generic in this package so that redundant dimensions can be
#' simplified for gridded objects.
#'
#' @param x An object.
#' @param ... Additional arguments passed to the method.
#'
#' @return The result of the selected method.
#' @rdname gridded-reshape
#' @export
drop = function(x, ...) {
  UseMethod("drop")
}

#' @describeIn gridded-reshape Default `drop()` method.
#' @export
drop.default = function(x, ...) {
  base::drop(x=x)
}

#' @describeIn gridded-reshape Drop a singleton third dimension from a `gts`
#'   object.
#' @export
drop.gts = function(x, ...) {
  ndim = dim(x)
  if(length(ndim)==3) return(x)
  if(ndim[3]==1) {
    dim(x$x) = ndim[-3]
    x$depth = NULL
    return(x)
  }
}

#' @describeIn gridded-reshape Simplify the data array of a `static` object by
#'   dropping redundant dimensions.
#' @export
drop.static = function(x, ...) {
  ndim = dim(x)
  if(length(ndim)==2) return(x)
  x$x = drop(x$x)
  return(x)
}
