#' Summarise and extract values from gridded objects
#'
#' These methods summarise the data stored in gridded objects or extract their
#' values in a simpler form.
#'
#' The page covers two related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension.
#' }
#'
#' `mean()`, `sum()`, and `sd()` compute scalar summaries for both classes.
#' For `gts` objects, these methods can also aggregate by a retained dimension.
#'
#' `as.double()` extracts the underlying data values as a numeric vector.
#'
#' The package also defines an `sd()` S3 generic so that `sd.gts()` and
#' `sd.static()` can coexist with a default method that delegates to
#' [stats::sd()].
#'
#' @param x A `gts` or `static` object.
#' @param by Optional aggregation target for `gts` objects. Supported values are
#'   `"time"`, `"space"`, `"latitude"`, and `"longitude"`.
#' @param trim Proportion to trim from each tail when computing the mean.
#' @param na.rm Logical; whether missing values should be removed.
#' @param ... Additional arguments passed to the underlying summary function.
#'
#' @details
#' For both `gts` and `static` objects, `as.double()` extracts values from the
#' underlying data array `x$x`. When a grid mask is present, only values at
#' non-missing mask cells are returned.
#'
#' For `static` objects, `mean.static()`, `sum.static()`, and `sd.static()` are
#' simple scalar summaries computed from `as.double(x)`.
#'
#' For `gts` objects:
#' \itemize{
#'   \item when `by = NULL`, `mean.gts()`, `sum.gts()`, and `sd.gts()` first
#'   coerce the object to a numeric vector with `as.double()` and then apply the
#'   corresponding summary function;
#'   \item when `by` is supplied, aggregation is performed with an internal
#'   helper that retains the requested dimension.
#' }
#'
#' Supported `by` values for `gts` objects are:
#' \describe{
#'   \item{`by = "time"`}{Summarise over the horizontal spatial dimensions and
#'   return a time series. For four-dimensional inputs, the result may be a
#'   multivariate `ts` object indexed by time.}
#'   \item{`by = "space"`}{Summarise over the last dimension and return a
#'   spatial array.}
#'   \item{`by = "latitude"`}{Retain the latitude dimension and summarise over
#'   longitude. Other trailing dimensions, such as depth or time, are preserved.}
#'   \item{`by = "longitude"`}{Retain the longitude dimension and summarise over
#'   latitude. Other trailing dimensions, such as depth or time, are preserved.}
#' }
#'
#' In the current implementation, grouped `gts` summaries always pass
#' `na.rm = TRUE` internally, regardless of the value supplied by the user.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`mean.gts()`, `sum.gts()`, `sd.gts()` with `by = NULL`}{A numeric
#'   scalar.}
#'   \item{`mean.gts()`, `sum.gts()`, `sd.gts()` with `by = "time"`}{A `ts`
#'   object.}
#'   \item{`mean.gts()`, `sum.gts()`, `sd.gts()` with another `by` value}{An
#'   array or vector retaining the requested dimension(s).}
#'   \item{`mean.static()`, `sum.static()`, `sd.static()`}{A numeric scalar.}
#'   \item{`as.double.gts()`, `as.double.static()`}{A numeric vector of extracted
#'   values.}
#'   \item{`sd.default()`}{A numeric result returned by [stats::sd()].}
#' }
#'
#' @seealso [quantile.gts()], [Summary.gts()], [Summary.static()], [stats::ts()]
#'
#' @examples
#' \dontrun{
#' mean(x)
#' mean(x, by = "time")
#' sd(x, by = "space")
#' mean(bathy)
#' as.double(bathy)
#' }
#' @name gridded-summary
NULL

#' @describeIn gridded-summary Compute the mean of a `gts` object, optionally by
#'   a retained dimension.
#' @export
mean.gts = function(x, by=NULL, trim=0, na.rm=FALSE, ...) {
  if(is.null(by)) return(mean(as.double(x), trim=trim, na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=mean, type=by, trim=trim, na.rm=TRUE, ...)
  return(out)
}

#' @describeIn gridded-summary Compute the mean of a `static` object.
#' @export
mean.static = function(x, trim=0, na.rm=FALSE, ...) {
  return(mean(as.double(x), trim=trim, na.rm=na.rm, ...))
}

#' @describeIn gridded-summary Compute the sum of a `gts` object, optionally by
#'   a retained dimension.
#' @export
sum.gts = function(x, by=NULL, na.rm=FALSE, ...) {
  if(is.null(by)) return(sum(as.double(x), na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=sum, type=by, na.rm=TRUE, ...)
  return(out)
}

#' @describeIn gridded-summary Compute the sum of a `static` object.
#' @export
sum.static = function(x, na.rm=FALSE, ...) {
  return(sum(as.double(x), na.rm=na.rm, ...))
}

#' @describeIn gridded-summary Extract the data values of a `gts` object as a
#'   numeric vector.
#' @export
as.double.gts = function(x, ...) {

  if(!is.null(x$grid$mask)) {
    return(as.double(x$x[!is.na(x$grid$mask)]))
  }
  return(as.double(x$x))
}

#' @describeIn gridded-summary Extract the data values of a `static` object as a
#'   numeric vector.
#' @export
as.double.static = as.double.gts

#' Standard deviation generic with methods for gridded objects
#'
#' `sd()` is an S3 generic in this package so that standard deviations can be
#' computed directly for `gts` and `static` objects.
#'
#' @param x An object.
#' @param na.rm Logical; whether missing values should be removed.
#' @param ... Additional arguments passed to the method.
#'
#' @return The result of the selected method.
#' @rdname gridded-summary
#' @seealso [quantile.gts()], [Summary.gts()], [Summary.static()], [stats::ts()]
#' @export
sd = function(x, ...) {
  UseMethod("sd")
}

#' @describeIn gridded-summary Default `sd()` method.
#' @export
sd.default = function(x, na.rm=FALSE, ...) {
  stats::sd(x=x, na.rm=na.rm)
}

#' @describeIn gridded-summary Compute the standard deviation of a `gts` object,
#'   optionally by a retained dimension.
#' @export
sd.gts = function(x, by=NULL, na.rm=FALSE, ...) {
  if(is.null(by)) return(sd(as.double(x), na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=sd, type=by, na.rm=TRUE, ...)
  return(out)
}

#' @describeIn gridded-summary Compute the standard deviation of a `static`
#'   object.
#' @export
sd.static = function(x, na.rm=FALSE, ...) {
  return(sd(as.double(x), na.rm=na.rm, ...))
}
