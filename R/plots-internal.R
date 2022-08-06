
#' @export
plot.gts = function(x, slice=1, ...) {
  if(mode(x$x)=="logical") x$x = 0L + x$x
  image.map(x$grid$LON, x$grid$LAT, x$x, slice=slice, ...)
  lab = sprintf("%s (%s)", x$info$varid, x$time[slice])
  title(main=lab)
  return(invisible(NULL))
}
