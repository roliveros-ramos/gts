

# merge -------------------------------------------------------------------

#' Merge a data.frame and a GTS object.
#'
#' @param x A data frame with dates and coordinates.
#' @param y A GTS object.
#' @param dt Time step when only 1 time value is used (the default is "month").
#' @param dimnames The names of the dimensions in the data.frame.
#' @param ... Additional arguments.
#'
#' @return A data.frame with the corresponding values of the variable from the GTS object.
#' @export
#'
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


