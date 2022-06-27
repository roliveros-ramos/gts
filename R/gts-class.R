
#' Read a netCDF variable into a Gridded Time-Series.
#'
#' @param filename The filename of the netCDF file.
#' @param varid Name of the variable in the file.
#' @param control Additional options, currently not used.
#' @param ... Additional arguments.
#'
#' @return A gts object
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

