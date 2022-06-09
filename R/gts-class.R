
#' Read a netCDF variable into a Gridded Time-Series.
#'
#' @param filename The filename of the netCDF file.
#' @param varid Name of the variable in the file.
#' @param control Additional options, currently not used.
#' @param ... Additional arguments.
#'
#' @return
#' @export
#'
read_gts = function(filename, varid, control=list(), ...) {

  nc = ncdf4::nc_open(filename = filename)
  on.exit(nc_close(nc))

  x = ncdf4::ncvar_get(nc, varid=varid, collapse_degen = FALSE)

  hasdepth = ifelse(length(dim(x))==4, TRUE, FALSE)

  time_unit = control$time_unit
  origin = control$origin

  dimx = ncvar_dim(nc, varid=varid, value=TRUE)
  dims = sapply(dimx, length)

  if(length(dimx)<3) stop("At leat (lon, lat, time) array is needed.")
  if(length(dimx)>4) stop("More than 4-dimensions are not supported.")

  time = tail(names(dimx), 1)

  if(is.null(time_unit)) {
    tmp = ncatt_get(nc, varid=time, attname = "units")
    tunit = ifelse(tmp$hasatt, tmp$value, NULL)
    punits = c("second", "minute", "hour", "day", "week", "month", "year")
    time_unit = punits[sapply(punits, FUN=grepl, x=tunit)]
    origin = parse_date_time(tunit, orders = c("Ymd", "YmdHMS", "dmY", "dmYHMS"))
  }

  time_conf = list(time_unit=time_unit, origin=origin)

  .getBreaks = function(x, conf) {
    out = c(x[1] - 0.5*(diff(x[1:2])),
            head(x, -1) + 0.5*diff(x),
            tail(x, 1) + 0.5*(diff(tail(x, 2))))
    return(out)
  }

  breaks = lapply(dimx, FUN=.getBreaks, conf=time_conf)

  pp = setNames(list(dimx[[time]]), nm=time_conf$time_unit)
  dimx[[time]] = time_conf$origin + do.call(period, args=pp)

  output = list(x = x,
                longitude = dimx[[1]],
                latitude = dimx[[2]],
                depth = if(hasdepth) dimx[[3]] else NULL,
                time = if(hasdepth) dimx[[4]] else dimx[[3]],
                breaks = breaks,
                info = list(varid=varid, time=time_conf))

  class(output) = c("gts", class(output))

  return(output)

}


# merge -------------------------------------------------------------------

#' Merge a data.frame and a GTS object.
#'
#' @param x A data frame with dates and coordinates.
#' @param y A GTS object.
#' @param dt Time step when only 1 time value is used (the default is "month").
#' @param dimnames The names of the dimensions in the data.frame.
#' @param ... Additional arguments.
#'
#' @return
#' @export
#'
merge_gts = function(x, y, dt="month", dimnames=NULL, ...) {

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


