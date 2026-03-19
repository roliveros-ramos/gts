#' Inspect and modify gridded objects
#'
#' These methods provide compact inspection tools for gridded objects and access
#' to selected metadata fields stored in `gts` and `static` objects.
#'
#' The page covers three related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension;
#'   \item `grid` objects, which store spatial geometry and associated metadata.
#' }
#'
#' `dim.gts()` and `dim.static()` return the dimensions of the data array stored
#' in `x$x`. `dim.grid()` returns the dimensions of the horizontal grid stored in
#' `x$LON`.
#'
#' `is.na.gts()` and `is.na.static()` return a copy of the input object with
#' `x$x` replaced by a logical array indicating missing values.
#'
#' `print.gts()`, `print.static()`, and `print.grid()` print a compact summary
#' of the object, including variable identifier, units, spatial extent, and grid
#' resolution. For `gts` objects, the printed summary also includes time
#' coverage and frequency.
#'
#' `names.gts()` and `names.static()` return the variable identifier and long
#' name stored in the metadata. The corresponding replacement methods update
#' those fields.
#'
#' `units.gts()` and `units.static()` return the units of the data variable,
#' that is, the last entry of `x$info$units`. The corresponding replacement
#' methods update that value.
#'
#' `str.gts()` and `str.static()` display the underlying list structure of the
#' object.
#'
#' `is_climatology()` tests whether a `gts` object is marked as a climatology.
#'
#' @param x A `gts`, `static`, or `grid` object, depending on the method.
#' @param object A `gts` or `static` object.
#' @param value For `names<-` methods, a character vector giving the variable
#'   identifier and long name. For `units<-` methods, a character scalar giving
#'   the data-variable units.
#' @param ... Additional arguments passed to the underlying generic. These are
#'   mainly relevant for `print.*()` and `str.*()`.
#'
#' @details
#' `names.gts()` and `names.static()` expose selected variable metadata rather
#' than dimension names. They return the variable identifier and long name stored
#' in the object's `info` component.
#'
#' `units.gts()` and `units.static()` expose only the units of the data variable,
#' not the units of coordinate dimensions.
#'
#' `str.gts()` and `str.static()` temporarily coerce the object to a plain list
#' before calling [utils::str()], so they display the underlying structure rather
#' than a shortened class-specific summary.
#'
#' `is_climatology()` is only defined for `gts` objects, because climatology is a
#' time-series concept in this package.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`is.na.gts()`, `is.na.static()`}{An object of the same class, with
#'   `x$x` replaced by a logical array of missing-value indicators.}
#'   \item{`dim.gts()`, `dim.static()`, `dim.grid()`}{An integer vector of
#'   dimensions.}
#'   \item{`print.gts()`, `print.static()`, `print.grid()`, `str.gts()`,
#'   `str.static()`}{The object is printed for inspection and the function
#'   returns `invisible(NULL)`.}
#'   \item{`names.gts()`, `names.static()`}{A character vector of length two
#'   containing the variable identifier and long name.}
#'   \item{`names<-.gts()`, `names<-.static()`}{The modified object.}
#'   \item{`units.gts()`, `units.static()`}{A character scalar giving the
#'   data-variable units.}
#'   \item{`units<-.gts()`, `units<-.static()`}{The modified object.}
#'   \item{`is_climatology()`}{A logical scalar.}
#' }
#'
#' @seealso [gts-class], [static-class], [grid-class], [gts()], [read_gts()], [read_static()], [make_grid()]
#'
#' @examples
#' \dontrun{
#' dim(x)
#' names(x)
#' units(x)
#' print(x)
#' str(x)
#' }
#' @name gridded-inspect
NULL

#' @describeIn gridded-inspect Identify missing values in a `gts` object.
#' @exportS3Method is.na gts
is.na.gts = function(x) {
  x$x = is.na(x$x)
  return(x)
}

#' @describeIn gridded-inspect Identify missing values in a `static` object.
#' @exportS3Method is.na static
is.na.static = is.na.gts

#' @describeIn gridded-inspect Return the dimensions of the data array stored in
#'   a `gts` object.
#' @export
dim.gts = function(x) dim(x$x)

#' @describeIn gridded-inspect Return the dimensions of the data array stored in
#'   a `static` object.
#' @export
dim.static = dim.gts

#' @describeIn gridded-inspect Return the dimensions of the horizontal grid
#'   stored in a `grid` object.
#' @export
dim.grid = function(x) dim(x$LON)

#' @describeIn gridded-inspect Print a compact summary of a `gts` object.
#' @export
print.gts = function(x, ...) {
  cat(sprintf("Gridded Time Series: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  cat(sprintf("Start     = %s\n", min(x$time)))
  cat(sprintf("End       = %s\n", max(x$time)))
  cat(sprintf("Frequency = %s\n", frequency(x)))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  print(resolution(x))
  return(invisible(NULL))
}

#' @describeIn gridded-inspect Print a compact summary of a `static` object.
#' @export
print.static = function(x, ...) {
  cat(sprintf("Static gridded object: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  print(resolution(x))
  return(invisible(NULL))
}

#' @describeIn gridded-inspect Print a compact summary of a `grid` object.
#' @export
print.grid = function(x, ...) {
  cat(sprintf("Grid object: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  print(resolution(x))
  return(invisible(NULL))
}

#' @describeIn gridded-inspect Return the variable identifier and long name
#'   stored in a `gts` object.
#' @export
names.gts = function(x) {
  return(c(x$info$varid, x$info$long_name))
}

#' @describeIn gridded-inspect Return the variable identifier and long name
#'   stored in a `static` object.
#' @export
names.static = names.gts

#' Replace variable identifier and long name
#'
#' Internal helper used by the `names<-` methods for gridded objects.
#'
#' @param x A `gts` or `static` object.
#' @param value A character vector containing the variable identifier and long
#'   name.
#'
#' @return The modified object.
#' @keywords internal
setnames_gts = function(x, value) {

  if(any(is.na(value))) stop("NAs are not allowed in names for a gridded object.")
  if(length(value) > 2) stop("A maximum of two values (varid, long_name) must be provided.")

  # TODO: decide behaviour for value[2] = NA

  x$info$varid = value[1]
  x$info$ovarid[length(x$info$ovarid)] = value[1]
  x$info$long_name = value[2]

  x

}

#' @describeIn gridded-inspect Replace the variable identifier and long name
#'   stored in a `gts` object.
#' @export
'names<-.gts' = setnames_gts

#' @describeIn gridded-inspect Replace the variable identifier and long name
#'   stored in a `static` object.
#' @export
'names<-.static' = setnames_gts

#' @describeIn gridded-inspect Return the data-variable units stored in a `gts`
#'   object.
#' @export
units.gts = function(x) {
  return(tail(x$info$units, 1))
}

#' @describeIn gridded-inspect Return the data-variable units stored in a
#'   `static` object.
#' @export
units.static = units.gts

#' Replace data-variable units
#'
#' Internal helper used by the `units<-` methods for gridded objects.
#'
#' @param x A `gts` or `static` object.
#' @param value A character scalar giving the units of the data variable.
#'
#' @return The modified object.
#' @keywords internal
setunits_gts = function(x, value) {

  if(length(value) != 1) stop("Only one value for units must be provided.")

  if(is.na(value)) value = ""

  x$info$units[length(x$info$units)] = value

  x

}

#' @describeIn gridded-inspect Replace the data-variable units stored in a `gts`
#'   object.
#' @export
'units<-.gts' = setunits_gts

#' @describeIn gridded-inspect Replace the data-variable units stored in a
#'   `static` object.
#' @export
'units<-.static' = setunits_gts

#' @describeIn gridded-inspect Display the underlying list structure of a `gts`
#'   object.
#' @export
str.gts = function(object, ...) {
  class(object) = "list"
  str(object, ...)
  return(invisible())
}

#' @describeIn gridded-inspect Display the underlying list structure of a
#'   `static` object.
#' @export
str.static = str.gts

#' @describeIn gridded-inspect Test whether a `gts` object is marked as a
#'   climatology.
#' @export
is_climatology = function(x, ...) {
  UseMethod("is_climatology")
}

#' @describeIn gridded-inspect Test whether a `gts` object is marked as a
#'   climatology.
#' @export
is_climatology.gts = function(x, ...) {
  return(x$info$climatology)
}
