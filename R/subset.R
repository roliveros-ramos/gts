#' Subset gridded objects
#'
#' These are methods for the base [base::subset()] generic applied to gridded
#' objects.
#'
#' The page covers three related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension;
#'   \item `grid` objects, which store spatial geometry and associated metadata.
#' }
#'
#' Spatial subsetting is performed by longitude and latitude ranges, or by using
#' the extent of another grid-like object. For `gts` objects, an additional
#' temporal subset can be requested by cycle position within the regular time
#' series.
#'
#' @param x A `gts`, `grid`, or `static` object.
#' @param longitude,latitude Numeric vectors defining the horizontal range to
#'   retain. Only the range of the supplied values is used.
#' @param grid A `grid` object, or a `gts` object from which the grid is
#'   extracted. When supplied, its spatial extent is used as the target subset
#'   extent.
#' @param frequency Integer vector of cycle values to retain when subsetting a
#'   `gts` object. Values must be positive integers not exceeding
#'   `frequency(x)`.
#' @param expand Numeric expansion applied to the target grid extent, typically
#'   in coordinate units. This may be a scalar, used for both longitude and
#'   latitude, or a vector of length two.
#' @param index.return Logical; if `TRUE`, `subset.grid()` includes the logical
#'   longitude and latitude indices used for subsetting in the returned object.
#' @param ... Additional arguments passed to the method.
#'
#' @details
#' `subset.grid()` is the low-level spatial subsetting method. It subsets the
#' horizontal coordinates, coordinate matrices, optional `area`, optional
#' `mask`, and the flattened coordinate table `df`. When `index.return = TRUE`,
#' it also stores the logical indices used for longitude and latitude in
#' `x$index`.
#'
#' `subset.gts()` first applies `subset.grid()` to the object's `grid`
#' component, then subsets the data array `x$x` along the first two spatial
#' dimensions. It updates the top-level longitude and latitude components,
#' recomputes spatial breaks for regular grids, updates the first two entries of
#' `x$info$dim`, and removes the temporary grid indices afterwards.
#'
#' If `frequency` is supplied, `subset.gts()` also subsets the last dimension by
#' cycle position using `cycle(x) %in% frequency`. After this operation, it
#' updates the stored time vector, time breaks, and time metadata, and sets
#' `x$info$ts` to `NULL`.
#'
#' `subset.static()` performs the same spatial subsetting workflow as
#' `subset.gts()`, but without any time-related operations.
#'
#' In the current implementation, `subset.grid()` does not update all possible
#' auxiliary grid components, such as `psi` or `prob`, and `subset.gts()` uses
#' cycle-based rather than date-based temporal filtering when `frequency` is
#' supplied.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`subset.gts()`}{A subsetted `gts` object.}
#'   \item{`subset.grid()`}{A subsetted `grid` object. If `index.return = TRUE`,
#'   the returned object also includes an `index` component with logical
#'   longitude and latitude indices.}
#'   \item{`subset.static()`}{A subsetted `static` object.}
#' }
#'
#' @seealso [base::subset()], [window.gts()], [grid-class], [gts-class],
#'   [static-class]
#'
#' @examples
#' \dontrun{
#' x1 <- subset(x, longitude = c(-80, -70), latitude = c(-20, -10))
#' x2 <- subset(x, grid = grd, expand = 1)
#' x3 <- subset(x, frequency = c(1, 2, 3))
#' grd2 <- subset(grd, longitude = c(2, 8), latitude = c(40, 45))
#' s2 <- subset(bathy, grid = grd2)
#' }
#' @name gridded-subset
NULL


#' @describeIn gridded-subset Subset a `gts` object in space and, optionally, by
#'   cycle position.
#' @export
subset.gts = function(x, longitude=NULL, latitude=NULL, grid=NULL, frequency=NULL, expand=0, ...) {

  if(!is.null(frequency)) {
    frequency = as.integer(frequency)
    if(any(is.na(frequency))) stop("Frequency argument must specify integers with the ts cycle.")
    if(any(frequency<1)) stop("Frequency must be positive.")
    if(any(frequency > frequency(x))) stop("Frequency values must be lower than the ts frequency.")
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

    x$time = time(x)[ind]
    x$breaks$time = .getBreaks(x$time)
    x$info$time$time = x$time
    x$info$dim$time = x$time
    x$info$ts = NULL

  }

  return(x)

}

#' @describeIn gridded-subset Subset a `grid` object by spatial extent.
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


#' @describeIn gridded-subset Subset a `static` object in space.
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
