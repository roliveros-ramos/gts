
#' @export
plot.gts = function(x, time=1, depth=1, ...) {
  args = list(...)
  if(mode(x$x)=="logical") x$x = 0L + x$x
  if(length(dim(x$x))==4) {
    image.map(x$grid$LON, x$grid$LAT, x$x[,,depth,], slice=time, ...)
    lab = sprintf("%s (depth=%s, time=%s)", x$info$varid, x$depth[depth], x$time[time])
  } else {
    image.map(x$grid$LON, x$grid$LAT, x$x, slice=time, ...)
    lab = sprintf("%s (time=%s)", x$info$varid, x$time[time])
  }
  main = if(is.null(args$main)) lab else NULL
  title(main=main)
  return(invisible(NULL))
}
