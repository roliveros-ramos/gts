
if(!isGeneric("merge")) {
  setGeneric("merge", function(x, y, ...)
    standardGeneric("merge"))
}

#' Merge point records with a gridded time-series object
#'
#' This page documents the behaviour of merging a `data.frame` with a `gts`
#' object through `merge(x, y)`.
#'
#' The package provides an S4 method for `merge(x = "data.frame", y = "gts")`
#' that matches rows of the data frame to the cells of a `gts` object using the
#' spatial and temporal breaks stored in the `gts` object, then appends the
#' corresponding gridded values to the data frame.
#'
#' `merge_gts()` is the exported implementation helper used by that S4 method.
#' It is exported for method dispatch and package interoperability, but typical
#' user code should call `merge(data.frame, gts)` directly.
#'
#' @param x A `data.frame` containing point records to be matched to the gridded
#'   object.
#' @param y A `gts` object providing the gridded variable and the dimension
#'   breaks used for matching.
#' @param dt Time unit used when `y` contains only a single time value. In that
#'   case, a time interval is constructed with [lubridate::floor_date()] and
#'   [lubridate::ceiling_date()] using `dt`. The default is `"month"`.
#' @param dimnames Named character vector giving the column names in `x`
#'   corresponding to the dimensions of `y`. By default,
#'   `c(longitude = "lon", latitude = "lat", depth = "depth", time = "time")`.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' Matching is performed by cutting each relevant column of `x` into the
#' corresponding breaks stored in `y$breaks`. Only dimensions present in
#' `y$breaks` are used.
#'
#' Rows of `x` are matched only when all required dimension indices are
#' available. Unmatched rows are kept in the output and receive `NA` for the
#' merged variable.
#'
#' The merged values are written to a column named `y$info$varid`. If that
#' column does not already exist in `x`, it is created. If it already exists, it
#' is overwritten for matched rows.
#'
#' When `length(y$time) == 1`, the time breaks stored in `y` are replaced
#' temporarily by an interval spanning the requested unit `dt`, so that point
#' records can be matched against a single gridded time slice.
#'
#' @return A `data.frame` containing the original columns of `x` plus one column
#'   named after `y$info$varid`, filled with values extracted from `y`.
#'
#' @seealso [base::merge()], [subset.gts()], [melt.gts()], [gts-class]
#'
#' @examples
#' \dontrun{
#' pts <- data.frame(
#'   lon = c(-75, -74.5),
#'   lat = c(-12, -11.5),
#'   time = as.Date(c("2000-01-15", "2000-02-15"))
#' )
#'
#' out <- merge(pts, x)
#' }
#' @name merge-gts
#' @aliases merge,data.frame,gts-method merge_gts
NULL

#' @describeIn merge-gts Implementation helper for merging a `data.frame` with a
#'   `gts` object.
#' @export
merge_gts = function(x, y, dt="month", dimnames=NULL, ...) {

  if(is.null(dimnames))
    dimnames = c(longitude="lon", latitude="lat", depth="depth", time="time")

  z = y$x
  breaks = y$breaks
  if(length(y$time) == 1) {
    breaks$time = as.Date(c(floor_date(y$time, unit = dt),
                            ceiling_date(y$time, unit=dt)))
  }

  i_lon = if(!is.null(breaks$longitude)) {
    cut(x[[dimnames["longitude"]]], breaks=breaks$longitude, labels = FALSE)
  } else NULL

  i_lat = if(!is.null(breaks$latitude)) {
    cut(x[[dimnames["latitude"]]], breaks=breaks$latitude, labels = FALSE)
  } else NULL

  i_dep = if(!is.null(breaks$depth)) {
    cut(x[[dimnames["depth"]]], breaks=breaks$depth, labels = FALSE)
  } else NULL

  i_tim = if(!is.null(breaks$time)) {
    cut(x[[dimnames["time"]]], breaks=breaks$time, labels = FALSE)
  } else NULL

  ind = cbind(i_lon, i_lat, i_dep, i_tim)
  ii = which(complete.cases(ind))

  if(is.null(x[[y$info$varid]])) x[[y$info$varid]] = NA
  x[[y$info$varid]][ii] = z[ind][ii]

  return(x)
}

#' @rdname merge-gts
setMethod('merge', signature(x='data.frame', y='gts'), merge_gts)
