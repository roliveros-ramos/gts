
#' @export
plot.gts = function(x, time=1, depth=1, ...) {
  args = list(...)
  if(mode(x$x)=="logical") x$x = 0L + x$x
  time_var = "time"
  if(!is.null(x$info$time_var)) time_var = x$info$time_var
  if(length(dim(x$x))==4) {
    image.map(x$grid$LON, x$grid$LAT, x$x[,,depth,], slice=time, ...)
    lab = sprintf("%s (depth=%s, %s=%s)", x$info$varid, x$depth[depth],
                  time_var, x$time[time])
  } else {
    image.map(x$grid$LON, x$grid$LAT, x$x, slice=time, ...)
    lab = sprintf("%s (%s=%s)", x$info$varid, time_var, x$time[time])
  }
  main = if(is.null(args$main)) lab else NULL
  title(main=main)
  return(invisible(NULL))
}

#' @export
plot.static = function(x, depth=1, ...) {
  args = list(...)
  if(mode(x$x)=="logical") x$x = 0L + x$x
  lon_check = !is.null(x$grid$longitude)
  lat_check = !is.null(x$grid$latitude)
  if(lon_check & lat_check) {
    ilon = x$grid$longitude
    ilat = x$grid$latitude
  } else {
    ilon = x$grid$LON
    ilat = x$grid$LAT
  }
  if(length(dim(x$x))==3) {
    image.map(ilon, ilat, x$x, slice=depth, ...)
    lab = sprintf("%s (depth=%s)", x$info$varid, x$depth[depth])
  } else {
    image.map(ilon, ilat, x$x, ...)
    lab = sprintf("%s", x$info$varid)
  }
  main = if(is.null(args$main)) lab else NULL
  title(main=main)
  return(invisible(NULL))
}

