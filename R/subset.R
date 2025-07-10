

#' Subsetting a spatial time series object
#'
#' @param x the object to be subsetted.
#' @param longitude The range of longitude to subset the grid.
#' @param latitude The range of latitude to subset the grid.
#' @param grid A grid object, or a 'gts' object.
#' @param expand Number of units to expand the grid, generally degrees.
#' It can be one number or two (longitude, latitude).
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return A 'gts' object after subsetting.
#' @export
#'
subset.gts = function(x, longitude=NULL, latitude=NULL, grid=NULL, frequency=NULL, expand=0, ...) {

  if(!is.null(frequency)) {
    frequency = as.integer(frequency)
    if(any(is.na(frequency))) stop("Frequency argument must specify integers with the ts cycle.")
    if(any(frequency<1)) stop("Frecuency must be positive.")
    if(any(frequency > frequency(x))) stop("Frecuency values must be lower than the ts frequency.")
  }

  if(!is.null(longitude) | !is.null(latitude) | !is.null(grid)) {

    x$grid = subset(x$grid, longitude=longitude, latitude=latitude, index.return=TRUE, grid=grid, expand=expand, ...)
    x$longitude = x$grid$longitude
    x$latitude = x$grid$latitude

    if(!is.matrix(x$longitude) & !is.matrix(x$latitude)) {
      x$breaks[[1]] = .getBreaks(x$longitude)
      x$breaks[[2]] = .getBreaks(x$latitude)
    }

    if(length(dim(x$x))==3) x$x = x$x[x$grid$index$lon, x$grid$index$lat, , drop=FALSE]
    if(length(dim(x$x))==4) x$x = x$x[x$grid$index$lon, x$grid$index$lat, , , drop=FALSE]

    x$info$dim[[1]] = x$info$dim[[1]][x$grid$index$lon]
    x$info$dim[[2]] = x$info$dim[[2]][x$grid$index$lat]

    if("i" %in% names(x$info$dim)) x$info$dim$i = seq_along(x$info$dim$i)
    if("j" %in% names(x$info$dim)) x$info$dim$j = seq_along(x$info$dim$j)

    x$grid$index = NULL

  }

  if(!is.null(frequency)) {

    if(any(is.na(frequency))) stop("Frequency argument must specify integers with the ts cycle.")

    ind = which(cycle(x) %in% frequency)
    if(length(dim(x$x))==3) x$x = x$x[, , ind, drop=FALSE]
    if(length(dim(x$x))==4) x$x = x$x[, , , ind, drop=FALSE]

    x$time = time(obs)[ind]
    x$breaks$time = .getBreaks(x$time)
    x$info$time$time = x$time
    x$info$dim$time = x$time
    x$info$ts = NULL

  }

  return(x)

}


#' @rdname subset.gts
#' @export
subset.grid = function(x, longitude=NULL, latitude=NULL, index.return=FALSE, grid=NULL, expand=0, ...) {

  if(is.null(longitude)&is.null(latitude)&is.null(grid)) return(x)

  if(length(expand)==1) expand = rep(expand, 2)
  if(length(expand)>2) stop("Argument 'expand' must be of length 1 or 2.")

  if(!is.null(grid)) {
    if(inherits(grid, "gts")) grid = grid$grid
    if(!inherits(grid, "grid")) stop("Argument 'grid' must be a grid object.")
    nlongitude = range(grid$longitude, na.rm=TRUE)
    nlatitude  = range(grid$latitude, na.rm=TRUE)
    if(is.null(longitude)) longitude = nlongitude + expand[1]*c(-1, 1)
    if(is.null(latitude))  latitude  = nlatitude + expand[2]*c(-1, 1)
  }

  lon = if(!is.null(longitude)) range(longitude, na.rm=TRUE) else NULL
  lat = if(!is.null(latitude))  range(latitude, na.rm=TRUE)  else NULL

  ilon = TRUE
  ilat = TRUE
  if(!is.null(lon)) ilon = x$LON >= lon[1] & x$LON <= lon[2]
  if(!is.null(lat)) ilat = x$LAT >= lat[1] & x$LAT <= lat[2]
  ind = ilat & ilon
  i_lon = setNames(apply(ind, 1, any), NULL)
  i_lat = setNames(apply(ind, 2, any), NULL)

  x$longitude = if(length(dim(x$longitude))==2) x$longitude[i_lon, i_lat, drop=FALSE] else x$longitude[i_lon]
  x$latitude = if(length(dim(x$latitude))==2) x$latitude[i_lon, i_lat, drop=FALSE] else x$latitude[i_lat]

  x$LON = x$LON[i_lon, i_lat, drop=FALSE]
  x$LAT = x$LAT[i_lon, i_lat, drop=FALSE]

  if(!is.null(x$area)) x$area = x$area[i_lon, i_lat, drop=FALSE]
  if(!is.null(x$mask)) x$mask = x$mask[i_lon, i_lat, drop=FALSE]

  x$df = data.frame(lon=as.numeric(x$LON), lat=as.numeric(x$LAT))

  if(isTRUE(index.return)) x$index = list(lon=i_lon, lat=i_lat)

  if(!inherits(x, "grid")) class(x) = c("grid", class(x))

  return(x)

}


#' @rdname subset.gts
#' @export
subset.static = function(x, longitude=NULL, latitude=NULL, grid=NULL, expand=0, ...) {

  x$grid = subset(x$grid, longitude=longitude, latitude=latitude, index.return=TRUE, grid=grid, expand=expand, ...)
  x$longitude = x$grid$longitude
  x$latitude = x$grid$latitude

  if(!is.matrix(x$longitude) & !is.matrix(x$latitude)) {
    x$breaks[[1]] = .getBreaks(x$longitude)
    x$breaks[[2]] = .getBreaks(x$latitude)
  }

  if(length(dim(x$x))==2) x$x = x$x[x$grid$index$lon, x$grid$index$lat, drop=FALSE]
  if(length(dim(x$x))==3) x$x = x$x[x$grid$index$lon, x$grid$index$lat, , drop=FALSE]

  x$info$dim[[1]] = x$info$dim[[1]][x$grid$index$lon]
  x$info$dim[[2]] = x$info$dim[[2]][x$grid$index$lat]

  if("i" %in% names(x$info$dim)) x$info$dim$i = seq_along(x$info$dim$i)
  if("j" %in% names(x$info$dim)) x$info$dim$j = seq_along(x$info$dim$j)

  x$grid$index = NULL

  return(x)

}
