#' Extract or modify horizontal coordinates of gridded objects
#'
#' These methods provide access to the horizontal coordinates stored in gridded
#' objects.
#'
#' `longitude()` and `latitude()` extract the horizontal coordinate definitions
#' from `gts`, `static`, and `grid` objects.
#'
#' `longitude<-` replaces the longitude coordinate of regular-grid `gts` and
#' `static` objects and reorders the associated data and grid components
#' accordingly.
#'
#' @param x A `gts`, `static`, or `grid` object.
#' @param prime_meridian Optional prime meridian convention used by
#'   `longitude()`. Supported values are `"center"` for longitudes centred on
#'   `[-180, 180]` and `"left"` for longitudes on `[0, 360]`. If omitted, the
#'   convention is inferred from the stored longitude values and recorded in the
#'   `"pm"` attribute of the returned object.
#' @param value Replacement longitude values for `longitude<-`.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' `longitude()` returns the stored longitude coordinates. When
#' `prime_meridian` is supplied, the returned values are converted to the
#' requested convention using the package longitude helpers.
#'
#' `latitude()` returns the stored latitude coordinates unchanged.
#'
#' `longitude<-.gts()` and `longitude<-.static()` currently support only regular
#' grids, that is, objects whose longitude coordinate is stored as a vector
#' rather than as a matrix. The replacement values must have the same length as
#' the original longitude coordinate and must not contain missing values.
#'
#' When longitude values are replaced, the methods reorder the data array and
#' associated grid components to keep the object aligned with the new longitude
#' ordering. They also update the longitude breaks and selected grid metadata.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`longitude()`}{The longitude coordinate vector or matrix. The
#'   returned object has a `"pm"` attribute indicating the prime meridian
#'   convention.}
#'   \item{`latitude()`}{The latitude coordinate vector or matrix.}
#'   \item{`longitude<-`}{The modified object.}
#' }
#'
#' @seealso [grid-class], [gts-class], [static-class], [transform_coordinates()]
#'
#' @examples
#' \dontrun{
#' lon <- longitude(x)
#' lon360 <- longitude(x, prime_meridian = "left")
#' lat <- latitude(x)
#'
#' longitude(x) <- seq(-80, -70, length.out = length(longitude(x)))
#' }
#' @name coordinates
NULL

#' @rdname coordinates
#' @export
longitude = function(x, prime_meridian = NULL, ...) {
  UseMethod("longitude")
}

#' @describeIn coordinates Extract the longitude coordinates of a `gts` object.
#' @export
longitude.gts = function(x, prime_meridian=NULL, ...) {
  lon = x$longitude
  if(!is.null(prime_meridian)) {
    pm = match.arg(prime_meridian, c("center", "left"))
    lon = checkLongitude(lon, primeMeridian=pm)
    attr(lon, "pm") = pm
  } else {
    pm = .findPrimeMeridian(lon)
    attr(lon, "pm") = pm
  }
  return(lon)
}

#' @describeIn coordinates Extract the longitude coordinates of a `static`
#'   object.
#' @export
longitude.static = longitude.gts

#' @describeIn coordinates Extract the longitude coordinates of a `grid` object.
#' @export
longitude.grid = longitude.gts

#' @rdname coordinates
#' @export
latitude = function(x, ...) {
  UseMethod("latitude")
}

#' @describeIn coordinates Extract the latitude coordinates of a `gts` object.
#' @export
latitude.gts = function(x, ...) {
  lat = x$latitude
  return(lat)
}

#' @describeIn coordinates Extract the latitude coordinates of a `static`
#'   object.
#' @export
latitude.static = latitude.gts

#' @describeIn coordinates Extract the latitude coordinates of a `grid` object.
#' @export
latitude.grid = latitude.gts

#' @rdname coordinates
#' @export
'longitude<-' = function(x, value) {
  UseMethod('longitude<-')
}

#' @describeIn coordinates Replace the longitude coordinates of a regular-grid
#'   `gts` object.
#' @export
'longitude<-.gts' = function(x, value) {

  if(identical(as.numeric(x$longitude), as.numeric(value))) return(x)

  nlon = length(x$longitude)
  if(length(value)!=nlon) stop("Longitud replacement does not have the right length.")
  if(any(is.na(value))) stop("NA values are not allowed in longitude replacement.")

  if(is.matrix(x$longitude)) stop("longitude modification is only implemented for regular grids.")

  i_lon = sort(value, index.return=TRUE)$ix
  value = value[i_lon]

  x$longitude = value
  x$x = if(length(dim(x$x))==3) x$x[i_lon, , , drop=FALSE] else x$x[i_lon, , , , drop=FALSE]
  # x$longitude = if(length(dim(x$longitude))==2) x$longitude[i_lon, , drop=FALSE] else x$longitude[i_lon]

  x$breaks$longitude = .getBreaks(x$longitude)

  grid = x$grid

  grid$longitude = x$longitude
  if(!is.null(grid$area)) grid$area = grid$area[i_lon, , drop=FALSE]
  if(!is.null(grid$mask)) grid$mask = grid$mask[i_lon, , drop=FALSE]

  grid$LON[] = value
  # grid$LON = grid$LON[i_lon, , drop=FALSE]
  grid$df = data.frame(lon=as.numeric(grid$LON), lat=as.numeric(grid$LAT))

  if(!is.null(grid$rho$LON)) grid$rho$LON = grid$LON
  if(!is.null(grid$psi$LON)) grid$psi$LON[] = x$breaks$longitude

  x$grid = grid

  x

}

#' @describeIn coordinates Replace the longitude coordinates of a regular-grid
#'   `static` object.
#' @export
'longitude<-.static' = function(x, value) {

  if(identical(as.numeric(x$longitude), as.numeric(value))) return(x)

  nlon = length(x$longitude)
  if(length(value)!=nlon) stop("Longitud replacement does not have the right length.")
  if(any(is.na(value))) stop("NA values are not allowed in longitude replacement.")

  if(is.matrix(x$longitude)) stop("longitude modification is only implemented for regular grids.")

  i_lon = sort(value, index.return=TRUE)$ix

  x$longitude = value
  x$x = if(length(dim(x$x))==2) x$x[i_lon, , drop=FALSE] else x$x[i_lon, , , drop=FALSE]
  x$longitude = if(length(dim(x$longitude))==2) x$longitude[i_lon, , drop=FALSE] else x$longitude[i_lon]

  x$breaks$longitude = .getBreaks(x$longitude)

  grid = x$grid

  grid$longitude = x$longitude
  if(!is.null(grid$area)) grid$area = grid$area[i_lon, , drop=FALSE]
  if(!is.null(grid$mask)) grid$mask = grid$mask[i_lon, , drop=FALSE]

  grid$LON[] = value
  grid$LON = grid$LON[i_lon, , drop=FALSE]
  grid$df = data.frame(lon=as.numeric(grid$LON), lat=as.numeric(grid$LAT))

  if(!is.null(grid$rho$LON)) grid$rho$LON = grid$LON
  if(!is.null(grid$psi$LON)) grid$psi$LON[] = x$breaks$longitude

  x$grid = grid

  x

}


# Auxiliar functions ------------------------------------------------------

#' Transform coordinates between coordinate reference systems
#'
#' `transform_coordinates()` transforms a set of two-dimensional coordinates
#' from one coordinate reference system (CRS) to another using `sf`.
#'
#' @param object A two-column matrix- or data-frame-like object containing
#'   coordinates. The first column is interpreted as `x` and the second as `y`.
#' @param crs_old Original CRS of the coordinates. This is passed to `sf` when
#'   constructing the input geometry.
#' @param crs_new Target CRS of the transformed coordinates. The default is
#'   `"EPSG:4326"`.
#'
#' @details
#' The function converts each coordinate pair to an `sf` point geometry,
#' constructs a simple-feature geometry set with the supplied input CRS, and
#' transforms it to `crs_new` with [sf::st_transform()].
#'
#' The returned coordinates are provided as a plain `data.frame` with columns
#' `x` and `y`.
#'
#' @return A `data.frame` with columns `x` and `y` giving the transformed
#'   coordinates.
#'
#' @seealso [longitude()], [latitude()], [sf::st_transform()]
#'
#' @examples
#' \dontrun{
#' pts <- data.frame(x = c(-75, -74), y = c(-12, -11))
#' transform_coordinates(pts, crs_old = "EPSG:4326", crs_new = "EPSG:3857")
#' }
#' @name transform_coordinates
#' @export
transform_coordinates = function(object, crs_old, crs_new="EPSG:4326") {
  x = object[, 1]
  y = object[, 2]
  .createPoint = function(i) st_point(c(x[i], y[i]))
  points = lapply(seq_along(x), FUN=.createPoint)
  points$crs = crs_old
  coords = do.call(st_sfc, points)
  dcoords = st_transform(coords, crs_new)
  newcoords = data.frame(x=sapply(dcoords, FUN="[", i=1),
                         y=sapply(dcoords, FUN="[", i=2))
  return(newcoords)
}


# Internal functions ------------------------------------------------------

.longitude2Center = function(x, ...) {
  if (!any(x > 180, na.rm=TRUE))
    return(x)
  x[which(x > 180)] = x[which(x > 180)] - 360
  return(x)
}

.longitude2Left = function(x, ...) {
  if (!any(x < 0, na.rm=TRUE))
    return(x)
  x[which(x < 0)] = x[which(x < 0)] + 360
  return(x)
}

checkLongitude = function(x, primeMeridian="center", ...) {
  if(all(is.na(x))) return(x)
  out = switch(primeMeridian,
               center = .longitude2Center(x, ...),
               left = .longitude2Left(x, ...))
  return(out)
}

.findPrimeMeridian = function(x) {
  if(all(is.na(x))) return("center")
  if(any(x<0)) return("center")
  if(any(x>180)) return("left")
  return("center")
}

