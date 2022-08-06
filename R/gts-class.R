
#' Read a netCDF variable into a Gridded Time-Series.
#'
#' @param filename The filename of the netCDF file.
#' @param varid String with  the name of the variable in the file to read the
#' data from. If left unspecified, the function will determine if there is only
#' one variable in the file and, if so, read from that. If left unspecified and
#' there are multiple variables in the file, an error is generated.
#' @param control A list with additional options, particularly the time units and origin when not specified in the file.
#' @param ... Additional arguments.
#'
#' @return A gts object
#' @export
#'
read_gts = function(filename, varid=NULL, control=list(), ...) {

  nc = nc_open(filename = filename)
  on.exit(nc_close(nc))

  output = gts(x=nc, varid=varid, control=control, ...)

  return(output)

}


# gts constructor ---------------------------------------------------------

#' Create a 'gts' object.
#'
#' @param x An array-like object to be a gridded time-series.
#' @param varid Name of the variable (e.g. from a ncdf file).
#' @param control A list with control arguments for the creation of a gts object. See details.
#' @param ... Additional arguments to be passed to different methods.
#'
#' @return A 'gts' object.
#' @export
#'
gts = function(x, ...) {
  UseMethod("gts")
}

#' @rdname gts
#' @export
gts.ncdf4 = function(x, varid=NULL, control=list(), ...) {

  nc = x

  if(is.null(varid)) {
    varid = names(nc$var)
    ndims = sapply(nc$var, FUN="[[", i="ndims")
    ind   = which(ndims>2)
    if(length(ind)==0) stop("No suitable variables found in file.")
    if(length(ind)!=1) {
      vs = paste(varid[ind], sep="", collapse=", ")
      message("Variables found in file:", vs)
      stop("More than one variable in the file, must specify 'varid'.")
    }
    varid = varid[ind]
    message(sprintf("Reading variable '%s'.", varid))
  }

  x = ncvar_get(nc, varid=varid, collapse_degen = FALSE)

  gatt = ncatt_get(nc, varid=0)

  tmp = ncatt_get(nc, varid=varid, attname = "units")
  var_unit = if(tmp$hasatt) tmp$value else ""

  hasdepth = if(length(dim(x))==4) TRUE else FALSE

  time_unit = control$time_unit
  origin = control$origin

  dimx = ncvar_dim(nc, varid=varid, value=TRUE)
  dims = sapply(dimx, length)

  if(length(dimx)<3) stop("At leat (lon, lat, time) array is needed.")
  if(length(dimx)>4) stop("More than 4-dimensions are not supported.")

  time = tail(names(dimx), 1)

  if(is.null(time_unit)) {
    tmp = ncatt_get(nc, varid=time, attname = "units")
    tunit = if(tmp$hasatt) tmp$value else NULL
    if(is.null(tunit)) stop("You must especify time unit.")
    punits = c("second", "minute", "hour", "day", "week", "month", "year")
    time_unit = punits[sapply(punits, FUN=grepl, x=tunit)]
    if(length(time_unit)==0) {
      message(sprintf("Time units '%s' are not recognized.", tunit))
      tab = guess_origin(dimx[[time]])
      message("Here are some guesses from time values, please check carefully.")
      print(tab)
      message(sprintf("Try read_gts(..., control=list(time_unit='%s', origin='%s'))",
                      tab$units[1], tab$origin[1]))
      stop("Must specify time unit and origin manually.")
    }
    origin = parse_date_time(tunit, orders = c("Ymd", "YmdHMS", "dmY", "dmYHMS"))
  }

  depth_conf = list(depth_unit = NULL, depth=NULL)
  depth_name = NULL
  if(hasdepth) {
    tmp = ncatt_get(nc, varid=names(dimx)[3], attname = "units")
    dunit = if(tmp$hasatt) tmp$value else NULL
    depth_conf = list(depth_unit=dunit, depth=dimx[[3]])
    depth_name = "depth"
  }

  if(is.null(origin)) stop("You must especify time origin.")

  time_conf = list(time_unit=time_unit, origin=origin, time=dimx[[time]],
                   units = sprintf("%s since %s", time_unit, origin))

  breaks = lapply(dimx, FUN=.getBreaks)

  new_time = time2date(dimx[[time]], units=time_conf$time_unit,
                       origin=time_conf$origin)

  tt = get_time(new_time[1])
  ff = get_freq(new_time)
  myts = ts(seq_along(new_time), start=c(year(new_time[1]),
                                         ceiling((tt%%1)*ff)), freq=ff)

  ilat = grep(x=tolower(names(nc$var)), pattern="lat")
  ilon = grep(x=tolower(names(nc$var)), pattern="lon")

  if(length(ilat)==1 & length(ilon)==1) {
    longitude = ncvar_get(nc, names(nc$var)[ilon])
    latitude  = ncvar_get(nc, names(nc$var)[ilat])
  } else {
    longitude = dimx[[1]]
    latitude  = dimx[[2]]
  }

  if(length(dim(longitude))<2) {
    LON = matrix(longitude, nrow=dims[1], ncol=dims[2])
    nlon = longitude
    lon_name = "longitude"
    lon_var = NULL
    dlon_unit = "degrees East"
    lon_unit = NULL
  } else {
    LON = longitude
    nlon = seq_len(dims[1])
    breaks[[1]] = NA
    lon_name = "i"
    lon_var = "longitude"
    dlon_unit = ""
    lon_unit = "degrees East"
  }

  if(length(dim(latitude))<2) {
    LAT = matrix(latitude, nrow=dims[1], ncol=dims[2], byrow=TRUE)
    nlat = latitude
    lat_name = "latitude"
    lat_var = NULL
    dlat_unit = "degrees North"
    lat_unit = NULL
  } else {
    LAT = latitude
    nlat = seq_len(dims[2])
    breaks[[2]] = NA
    lat_name = "j"
    lat_var = "latitude"
    dlat_unit = ""
    lat_unit = "degrees North"
  }

  dimnames(LON) = dimnames(LAT) = list(x=nlon, y=nlat)

  odim = if(hasdepth) list(nlon, nlat, dimx[[3]], dimx[[time]]) else
    list(nlon, nlat, dimx[[time]])
  names(odim) = c(lon_name, lat_name, depth_name, "time")

  ovar = c(lon_var, lat_var, "x")
  ovarid = c(lon_var, lat_var, varid)

  dim.units = c(dlon_unit, dlat_unit, depth_conf$depth_unit, time_conf$units)
  names(dim.units) = c(lon_name, lat_name, depth_name, "time")

  units = c(lon_unit, lat_unit, var_unit)

  grid = list(longitude=longitude, latitude=latitude, rho=list(LON=LON, LAT=LAT),
              psi=NULL, LON=LON, LAT=LAT, area=NULL, mask=NULL,
              df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))
  class(grid) =  c("grid", class(grid))

  grid = fill(grid)
  grid$area = suppressMessages(area(grid))

  output = list(x = x,
                longitude = longitude,
                latitude  = latitude,
                depth     = if(hasdepth) dimx[[3]] else NULL,
                time      = new_time,
                breaks    = breaks,
                grid      = grid,
                info      = list(varid=varid, time=time_conf, depth=depth_conf,
                                 dim = odim, var=ovar, dim.units=dim.units,
                                 units=units, ovarid=ovarid, global=gatt, ts=myts))

  class(output) = c("gts", "ts")

  return(output)

}
