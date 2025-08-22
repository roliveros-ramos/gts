
#' Compute the spatial resolution of a grid
#'
#' @param x A gridded object, gts, static or grid.
#' @param ... Additional parameters for computation.
#'
#' @return The spatial resolution.
#' @export
#'
resolution = function(x, ...) {
  UseMethod("resolution")
}


# S3 methods --------------------------------------------------------------

#' @export
resolution.default = function(x, ...) {
  dx = unique(round(diff(longitude(x)), 12))
  dy = unique(round(diff(latitude(x)), 12))
  out = list(dx=dx, dy=dy)
  class(out) = "gts.resolution"
  return(out)
}

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
