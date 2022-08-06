
#' @export
subset.gts = function(x, longitude=NULL, latitude=NULL, ...) {

  x$grid = subset(x$grid, longitude=longitude, latitude=latitude, index.return=TRUE, ...)
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

  return(x)

}

#' @export
subset.grid = function(x, longitude=NULL, latitude=NULL, index.return=FALSE, ...) {

  if(is.null(longitude)&is.null(latitude)) return(x)

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

