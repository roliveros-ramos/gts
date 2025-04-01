
#' @export
longitude = function(x, ...) {
  UseMethod("longitude")
}

#' @export
longitude.gts = function(x, prime_meridian=NULL) {
  lon = x$longitude
  if(!is.null(prime_meridian)) {
    pm = match.arg(prime_meridian, c("center", "left"))
    lon = checkLongitude(lon, primeMeridian=pm)
    attr(lon, "pm") = pm
  }
  return(lon)
}

#' @export
latitude = function(x, ...) {
  UseMethod("longitude")
}

#' @export
latitude.gts = function(x, prime_meridian=NULL) {
  lat = x$latitude
  return(lat)
}


#' @export
'longitude<-' = function(x, value) {
  UseMethod('longitude<-')
}

#' @export
'longitude<-.gts' = function(x, value) {

  if(is.matrix(x$longitude)) stop("longitude modification is only implemented for regular grids.")

  nlon = length(x$longitude)

  if(length(value)!=nlon) stop("Longitud replacement does not have the right length.")

  if(any(is.na(value))) stop("NA values are not allowed in longitude replacement.")

  i_lon = sort(value, index.return=TRUE)$ix

  x$longitude = value
  x$x = if(length(dim(x$x))==3) x$x[i_lon, , , drop=FALSE] else x$x[i_lon, , , , drop=FALSE]
  x$longitude = if(length(dim(x$longitude))==2) x$longitude[i_lon, , drop=FALSE] else x$longitude[i_lon]

  x$breaks$lon = .getBreaks(x$longitude)

  grid = x$grid

  grid$longitude = x$longitude
  if(!is.null(grid$area)) grid$area = grid$area[i_lon, , drop=FALSE]
  if(!is.null(grid$mask)) grid$mask = grid$mask[i_lon, , drop=FALSE]

  grid$LON[] = value
  grid$LON = grid$LON[i_lon, , drop=FALSE]
  grid$df = data.frame(lon=as.numeric(grid$LON), lat=as.numeric(grid$LAT))

  if(!is.null(attr(value, "pm"))) {
    pm = attr(value, "pm")
    if(!is.null(grid$rho$LON))
      grid$rho$LON = checkLongitude(grid$rho$LON, primeMeridian=pm)
    if(!is.null(grid$psi$LON))
      grid$rho$LON = checkLongitude(grid$psi$LON, primeMeridian=pm)
  }

  x$grid = grid

  x

}


