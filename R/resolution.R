#' Compute the spatial resolution of a gridded object
#'
#' `resolution()` computes the horizontal spatial resolution of a gridded object
#' from its longitude and latitude coordinates.
#'
#' The default method computes successive differences in longitude and latitude,
#' rounds them to 12 decimal places, and returns the unique values found in each
#' direction. This allows both regular grids, with a single resolution value per
#' axis, and irregular grids, with varying spacing, to be summarised.
#'
#' Objects returned by `resolution()` have class `"gts.resolution"` and can be
#' printed with a dedicated method.
#'
#' @param x A gridded object providing `longitude()` and `latitude()` methods,
#'   such as a `gts`, `static`, or `grid` object.
#' @param ... Additional arguments passed to the method.
#'
#' @return
#' `resolution()` returns an object of class `"gts.resolution"`, implemented as
#' a list with components:
#' \describe{
#'   \item{`dx`}{Unique longitude spacings.}
#'   \item{`dy`}{Unique latitude spacings.}
#' }
#'
#' For regular grids, `dx` and `dy` typically each contain a single value. For
#' grids with varying spacing, they may contain several values.
#'
#' @seealso [print.gts.resolution()], [longitude()], [latitude()]
#'
#' @examples
#' \dontrun{
#' r <- resolution(x)
#' print(r)
#' }
#' @name resolution
#' @export
resolution = function(x, ...) {
  UseMethod("resolution")
}

#' @describeIn resolution Compute the spatial resolution from longitude and
#'   latitude coordinates.
#' @export
resolution.default = function(x, ...) {
  # TODO: check this for irregular grids
  dx = unique(round(diff(longitude(x)), 12))
  dy = unique(round(diff(latitude(x)), 12))
  out = list(dx=dx, dy=dy)
  class(out) = "gts.resolution"
  return(out)
}

#' Print a spatial resolution summary
#'
#' `print.gts.resolution()` prints a compact summary of the horizontal spatial
#' resolution returned by [resolution()].
#'
#' If a single spacing is present, that value is printed directly. If several
#' spacings are present, the method prints the mean spacing together with the
#' minimum and maximum values.
#'
#' @param x An object of class `"gts.resolution"`.
#' @param ... Additional unused arguments.
#'
#' @return The input object, invisibly.
#'
#' @rdname print.gts.resolution
#' @export
print.gts.resolution = function(x, ...) {
  degree = "\U00B0"
  msg0 = "%0.2g\U00B0"
  msg1 = "%0.2g\U00B0 (average), [%0.3g-%0.3g]"
  mlon = if(length(x$dx)==1) sprintf(msg0, x$dx) else sprintf(msg1, mean(x$dx, na.rm=TRUE), min(x$dx), max(x$dx))
  mlat = if(length(x$dy)==1) sprintf(msg0, x$dy) else sprintf(msg1, mean(x$dy, na.rm=TRUE), min(x$dy), max(x$dy))
  cat("Spatial resolution:\n")
  cat("\tLongitude :", mlon, "\n")
  cat("\tLatitude  :", mlat)
  invisible(x)
}
