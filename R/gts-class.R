
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
gts.ncdf4 = function(x, varid=NULL, climatology=FALSE, control=list(), ...) {

  nc = x

  # find the variable to read
  if(is.null(varid)) {
    varid = names(nc$var)
    ndims = sapply(nc$var, FUN="[[", i="ndims")
    ind   = which(ndims>2)
    if(length(ind)==0) stop("No suitable variables found in file.")
    if(length(ind)!=1) {
        vs = paste(varid[ind], sep="", collapse=", ")
      if(!is.null(control$use_first)) {
        ind = ind[1]
        if(!isTRUE(control$use_first)) {
          message("Variables found in file: ", vs)
          stop("More than one variable in the file, must specify 'varid'.")
        }
      } else {
        message("Variables found in file: ", vs)
        stop("More than one variable in the file, must specify 'varid'.")
      }
    }
    varid = varid[ind]
    message(sprintf("Reading variable '%s'.", varid))
  }

  # get global attributes, units and long_name
  gatt = ncatt_get(nc, varid=0)

  tmp = ncatt_get(nc, varid=varid, attname = "units")
  var_unit = if(tmp$hasatt) tmp$value else ""

  tmp = ncatt_get(nc, varid=varid, attname = "long_name")
  long_name = if(tmp$hasatt) tmp$value else varid

  # get dimensions, including depth
  ndimx = nc$var[[varid]]$varsize
  hasdepth = if(length(ndimx)==4) TRUE else FALSE

  # get time information
  time_unit = control$time_unit
  origin = control$origin
  tcal = control$calendar

  # get dimensions information
  dimx = ncvar_dim(nc, varid=varid, value=TRUE)
  dims = sapply(dimx, length)
  time = tail(names(dimx), 1)

  if(length(dimx)<3) stop("At least (lon, lat, time) array is needed.")
  if(length(dimx)>4) stop("More than 4-dimensions are not supported.")

  # check for monotonicity in dimensions
  is_monot = all(sapply(dimx, FUN=is_monotonic))
  if(!is_monot) stop("Variable dimensions have non-monotonic values.")

  # get the data
  x = ncvar_get(nc, varid=varid, collapse_degen = FALSE)

  ## validation of dimensions
  if(is_decreasing(dimx[[1]])) {
    warning("First dimension has decreasing values: fixing it.")
    dimx[[1]] = rev(dimx[[1]])
    ind = rev(seq_along(dimx[[1]]))
    x = if(hasdepth) x[ind, , , , drop=FALSE] else x[ind, , , drop=FALSE]
  }

  if(is_decreasing(dimx[[2]])) {
    warning("Second dimension has decreasing values: fixing it.")
    dimx[[2]] = rev(dimx[[2]])
    ind = rev(seq_along(dimx[[2]]))
    x = if(hasdepth) x[, ind, , , drop=FALSE] else x[, ind, , drop=FALSE]
  }

  if(hasdepth) {
    if(is_decreasing(dimx[[3]])) {
      warning("Third dimension has decreasing values: fixing it.")
      dimx[[3]] = rev(dimx[[3]])
      ind = rev(seq_along(dimx[[3]]))
      x = x[, , ind, , drop=FALSE]
    }
  }

  if(is_decreasing(dimx[[time]])) {
    warning("Time dimension has decreasing values: fixing it.")
    dimx[[time]] = rev(dimx[[time]])
    ind = rev(seq_along(dimx[[time]]))
    x = if(hasdepth) x[, , , ind, drop=FALSE] else x[, , ind, drop=FALSE]
  }
  ## END validation of dimensions

  # try to guess time unit and origin
  if(is.null(time_unit)) {
    tmp = ncatt_get(nc, varid=time, attname = "units")
    tunit = if(tmp$hasatt) tmp$value else NULL
    if(is.null(tunit)) stop("You must especify time unit.")
    punits = c("second", "minute", "hour", "day", "week", "fortnight", "month", "year")
    time_unit = punits[sapply(punits, FUN=grepl, x=tunit)]
    if(length(time_unit)==0) {
      message(sprintf("Time units '%s' are not recognized.", tunit))
      tab = guess_origin(dimx[[time]])
      message("Here are some guesses from time values, please check carefully.")
      print(tab)
      message(sprintf("Try 'read_gts(..., control=list(time_unit='%s', origin='%s'))'",
                      tab$time_unit[1], tab$origin[1]))
      stop("Must specify time unit and origin manually.")
    }
    origin = parse_date_time(tunit, orders = c("Ymd", "YmdHMS", "dmY", "dmYHMS"))
    tmp = ncatt_get(nc, varid=time, attname = "calendar")
    tcal = if(tmp$hasatt) tmp$value else NULL
    if(!is.null(tcal)) {
      if(tcal=="gregorian") tcal = "365.2425"
      tcal = as.numeric(gsub(x=tcal, pattern="[^0-9].*", replacement = ""))
      if(!(tcal %in% c(360, 365, 364, 364.25, 365.2425, 365.25)))
        warning(sprintf("Using calendar year of %s days.", tcal))
    } else {
      tcal = 365
    }
  }

  # get depth information
  depth_conf = list(depth_unit = NULL, depth=NULL)
  depth_name = NULL
  if(hasdepth) {
    tmp = ncatt_get(nc, varid=names(dimx)[3], attname = "units")
    dunit = if(tmp$hasatt) tmp$value else NULL
    depth_conf = list(depth_unit=dunit, depth=dimx[[3]])
    depth_name = "depth"
  }

  # if no origin, ask and throw an error
  if(is.null(origin)) stop("You must especify time origin.")

  time_conf = list(time_unit=time_unit, origin=origin, time=dimx[[time]],
                   units = sprintf("%s since %s", time_unit, origin),
                   calendar=tcal)

  # compute breaks
  breaks = lapply(dimx, FUN=.getBreaks)

  # create time
  new_time = time2date(dimx[[time]], units=time_conf$time_unit,
                       origin=time_conf$origin, calendar=time_conf$calendar)
  tt = get_time(new_time[1])
  ff = get_freq(new_time, freq=control$frequency)
  myts = ts(seq_along(new_time), start=c(year(new_time[1]),
                                         floor((tt%%1)*ff) + 1), frequency=ff)

  if(isTRUE(climatology) & (length(new_time)!=ff)) {
    warning("Data is not climatological.")
    climatology = FALSE
  }

  # get longitude and latitude
  ilat = grep(x=tolower(names(nc$var)), pattern="^lat")
  ilon = grep(x=tolower(names(nc$var)), pattern="^lon")

  if(length(ilat)==1 & length(ilon)==1) {
    longitude = ncvar_get(nc, names(nc$var)[ilon])
    latitude  = ncvar_get(nc, names(nc$var)[ilat])
  } else {
    longitude = dimx[[1]]
    latitude  = dimx[[2]]
  }

  if(length(dim(longitude))<2) {
    LON = matrix(longitude, nrow=dims[1], ncol=dims[2])
    pLON = matrix(breaks[[1]], nrow=dims[1]+1, ncol=dims[2]+1)
    nlon = longitude
    lon_name = "longitude"
    lon_var = NULL
    dlon_unit = "degrees East"
    lon_unit = NULL
  } else {
    LON = longitude
    pLON = NULL
    nlon = seq_len(dims[1])
    breaks[[1]] = NA
    lon_name = "i"
    lon_var = "longitude"
    dlon_unit = ""
    lon_unit = "degrees East"
  }

  if(length(dim(latitude))<2) {
    LAT = matrix(latitude, nrow=dims[1], ncol=dims[2], byrow=TRUE)
    pLAT = matrix(breaks[[2]], nrow=dims[1]+1, ncol=dims[2]+1)
    nlat = latitude
    lat_name = "latitude"
    lat_var = NULL
    dlat_unit = "degrees North"
    lat_unit = NULL
  } else {
    LAT = latitude
    pLAT = NULL
    nlat = seq_len(dims[2])
    breaks[[2]] = NA
    lat_name = "j"
    lat_var = "latitude"
    dlat_unit = ""
    lat_unit = "degrees North"
  }

  dimnames(LON) = dimnames(LAT) = list(x=nlon, y=nlat)
  # END get longitude and latitude

  # variables and units
  odim = if(hasdepth) list(nlon, nlat, dimx[[3]], dimx[[time]]) else
    list(nlon, nlat, dimx[[time]])
  names(odim) = c(lon_name, lat_name, depth_name, "time")

  ovar = c(lon_var, lat_var, "x")
  ovarid = c(lon_var, lat_var, varid)

  dim.units = c(dlon_unit, dlat_unit, depth_conf$depth_unit, time_conf$units)
  names(dim.units) = c(lon_name, lat_name, depth_name, "time")

  units = c(lon_unit, lat_unit, var_unit)

  # create grid
  grid = list(longitude=longitude, latitude=latitude, rho=list(LON=LON, LAT=LAT),
              psi=list(LON=pLON, LAT=pLAT), LON=LON, LAT=LAT, area=NULL, mask=NULL,
              df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))
  class(grid) =  c("grid", class(grid))

  grid = fill(grid, control=control)
  grid$area = suppressMessages(area(grid))
  # END: create grid

  # constructor
  output = list(x = x,
                longitude = longitude,
                latitude  = latitude,
                depth     = if(hasdepth) dimx[[3]] else NULL,
                time      = new_time,
                breaks    = breaks,
                grid      = grid,
                info      = list(varid=varid, long_name=long_name,
                                 time=time_conf, depth=depth_conf,
                                 dim = odim, var=ovar, dim.units=dim.units,
                                 units=units, ovarid=ovarid, global=gatt, ts=myts,
                                 climatology=climatology))

  class(output) = "gts"

  return(output)

}



#' @rdname gts
#' @export
gts.array = function(x, varid=NULL, grid, data = NA, start = 1, end=numeric(),
                     frequency = 1, deltat = 1, climatology=FALSE, long_name=NULL, unit=NULL, control=list(), ...) {

  if(is.null(varid)) varid = ""
  if(is.null(long_name)) long_name = ""
  var_unit = if(is.null(unit))  "" else unit

  if(length(dim(x))>3) x = drop(x)

  ntime = ndata = tail(dim(x), 1)

  if(isTRUE(climatology)) {
    frequency = tail(dim(x), 1)
    start = 0
  }

  if(missing(frequency)) {
    frequency = 1/deltat
  } else {
    if(missing(deltat)) deltat = 1/frequency
  }

  if(!missing(start)) start = as.numeric(start)
  if (length(start) > 1L) {
    start = start[1L] + (start[2L] - 1)/frequency
  }
  if(!missing(end)) end = as.numeric(end)
  if(length(end) > 1L) {
    end = end[1L] + (end[2L] - 1)/frequency
  }
  if(missing(end)) {
    end = start + (ndata - 1)/frequency
  } else {
    if(missing(start)) start = end - (ndata - 1)/frequency
  }

  if(inherits(grid, "gts")) grid = grid$grid
  if(!inherits(grid, "grid")) stop("grid must be of class 'grid'.")

  ndimx = dims = dim(x)

  if(length(ndimx)<3) stop("At least (lon, lat, time) array is needed.")
  if(length(ndimx)>4) stop("More than 4-dimensions are not supported.")

  dimx = lapply(dimnames(x), as.numeric)
  if(length(dimx)==0) dimx = lapply(dim(x), seq_len)

  hasdepth = if(length(ndimx)==4) TRUE else FALSE

  # get depth information
  depth_conf = list(depth_unit = NULL, depth=NULL)
  depth_name = NULL
  if(hasdepth) {
    dunit = ifelse(is.null(control$depth_unit), "", control$depth_unit)
    depth_conf = list(depth_unit=dunit, depth=dimx[[3]])
    depth_name = "depth"
  }

  # check grid is appropiate

  # compute breaks
  breaks = lapply(dimx, FUN=.getBreaks)

  # create time
  origin = "0000-01-01 12:00:00.5"
  time_unit = "year"
  myts = ts(seq_len(ntime), start=start, end=end, frequency=frequency, deltat=deltat, ...)
  tt = time(myts)
  new_time = time2date(tt, units=time_unit, origin=origin)

  time_conf = list(time_unit="year", origin=origin, time=tt,
                   units = sprintf("%s since %s", time_unit, origin),
                   calendar=NULL)

  time = length(ndimx)

  longitude = grid$longitude
  latitude  = grid$latitude

  if(length(dim(longitude))<2) {
    nlon = longitude
    lon_name = "longitude"
    lon_var = NULL
    dlon_unit = "degrees East"
    lon_unit = NULL
  } else {
    nlon = seq_len(dims[1])
    lon_name = "i"
    lon_var = "longitude"
    dlon_unit = ""
    lon_unit = "degrees East"
  }

  if(length(dim(latitude))<2) {
    nlat = latitude
    lat_name = "latitude"
    lat_var = NULL
    dlat_unit = "degrees North"
    lat_unit = NULL
  } else {
    nlat = seq_len(dims[2])
    lat_name = "j"
    lat_var = "latitude"
    dlat_unit = ""
    lat_unit = "degrees North"
  }

  # compute breaks
  dimx[[1]] = nlon
  dimx[[2]] = nlat
  dimx[[time]] = tt

  breaks = lapply(dimx, FUN=.getBreaks)

  # variables and units
  odim = if(hasdepth) list(nlon, nlat, dimx[[3]], dimx[[time]]) else
    list(nlon, nlat, dimx[[time]])
  names(odim) = c(lon_name, lat_name, depth_name, "time")

  ovar = c(lon_var, lat_var, "x")
  ovarid = c(lon_var, lat_var, varid)

  dim.units = c(dlon_unit, dlat_unit, depth_conf$depth_unit, time_conf$units)
  names(dim.units) = c(lon_name, lat_name, depth_name, "time")

  units = c(lon_unit, lat_unit, var_unit)

  # constructor
  output = list(x = x,
                longitude = grid$longitude,
                latitude  = grid$latitude,
                depth     = if(hasdepth) dimx[[3]] else NULL,
                time      = new_time,
                breaks    = breaks,
                grid      = grid,
                info      = list(varid=varid, long_name=long_name,
                                 time=time_conf, depth=depth_conf,
                                 dim = odim, var=ovar, dim.units=dim.units,
                                 units=units, ovarid=ovarid, global=NULL, ts=myts,
                                 climatology=climatology))

  class(output) = "gts"

  return(output)

}








