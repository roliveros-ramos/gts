#' Classify points relative to land, ocean, or a polygon layer
#'
#' These helpers classify point locations from longitude and latitude
#' coordinates.
#'
#' `is_ocean()` returns `TRUE` for points that are not covered by the land
#' polygon layer used by the package.
#'
#' `is_land()` is the logical complement of `is_ocean()`.
#'
#' `is_over()` returns `TRUE` for points that are covered by a user-supplied
#' polygon layer.
#'
#' @param lon,lat Numeric vectors of longitude and latitude coordinates.
#' @param hires Logical; if `TRUE`, use the high-resolution land polygon layer in
#'   [is_ocean()] and [is_land()].
#' @param layer A sf polygon.
#'
#' @details
#' All functions preserve missing coordinate pairs as `NA` in the output.
#'
#' `is_ocean()` tests whether each complete point falls outside the package land
#' layer (`land_lr` or `land_hr`).
#'
#' `is_land()` is implemented as `!is_ocean(...)`, so `TRUE` indicates land,
#' `FALSE` indicates ocean, and missing values are preserved.
#'
#' `is_over()` uses [sp::over()] with the supplied polygon layer and returns
#' `TRUE` for points that intersect that layer.
#'
#' @return A logical vector of the same length as `lon` and `lat`.
#'
#' @seealso [make_grid()], [mask()]
#'
#' @examples
#' n <- 100
#' 
#' exDF <- data.frame(
#'   lon = runif(n = n, min = -90, max = -70),
#'   lat = runif(n = n, min = -20, max = -2)
#' )
#' 
#' is_ocean(lon = exDF$lon, lat = exDF$lat)
#' is_land(lon = exDF$lon, lat = exDF$lat)
#' @name is_ocean
NULL

#' @describeIn is_ocean Test whether points fall in the ocean.
#' @export
is_ocean = function(lon, lat, hires=FALSE) !is_land(lon, lat, hires)

#' @describeIn is_ocean Test whether points fall on land.
#' @export
is_land = function(lon, lat, hires=FALSE){
  is_over(lon = lon, lat = lat, layer = if(isTRUE(hires)) land_hr else land_lr)
} 

#' @describeIn is_ocean Test whether points are covered by a polygon layer.
#' @export
#' @inheritParams is_ocean
is_over = function(lon, lat, layer) {
  
  coords = data.frame(lon = lon, lat = lat)
  coords = st_as_sf(x = coords, coords = c("lon", "lat"))
  
  st_crs(coords) = st_crs(layer)
  
  out = st_join(x = coords, y = cbind(layer, xind = TRUE))
  
  return(!is.na(out$xind))
}

# Internal functions ------------------------------------------------------

.create_grid = function(lon, lat, dx, dy, center=FALSE) {

  # Create a rectangular grid given lat, lon and dxy.
  # No correction by Earth curvature
  if(dx <= 0 || dy <= 0) stop("dx and dy must be positive.")
  if(diff(lon) <= 0) stop("Longitude values must not be decreasing.")
  if(diff(lat) <= 0) stop("Latitude values must not be decreasing.")

  nx = floor(diff(lon)/dx)
  ny = floor(diff(lat)/dy)

  dlon = dx*nx - diff(lon)
  dlat = dy*ny - diff(lat)

  if(abs(dlon) > 1e-6*dx) {
    nx  = nx + 1
    dlon = dx*nx - diff(lon)
    lon = lon + 0.5*dlon*c(-1, +1)
  }
  if(abs(dlat) > 1e-6*dy) {
    ny  = ny + 1
    dlat = dy*ny - diff(lat)
    lat = lat + 0.5*dlat*c(-1, +1)
  }

  if(isTRUE(center)) {
    lat[which.min(lat)] = lat[which.min(lat)] + 0.5*dy
    lat[which.max(lat)] = lat[which.max(lat)] - 0.5*dy
    lon[which.min(lon)] = lon[which.min(lon)] + 0.5*dx
    lon[which.max(lon)] = lon[which.max(lon)] - 0.5*dx
    nx = nx + 1
    ny = ny + 1
  }

  lats.psi = seq(from=min(lat),to=max(lat), length=ny + 1)
  lons.psi = seq(from=min(lon),to=max(lon), length=nx + 1)

  lats.rho = seq(from=min(lat) + 0.5*dy, to=max(lat) - 0.5*dy, length=ny)
  lons.rho = seq(from=min(lon) + 0.5*dx, to=max(lon) - 0.5*dx, length=nx)

  rho = list(lat=lats.rho, lon=lons.rho)
  psi = list(lat=lats.psi, lon=lons.psi)

  nlat = length(rho$lat)
  nlon = length(rho$lon)

  LAT = matrix(rho$lat, ncol=nlat, nrow=nlon, byrow=TRUE)
  LON = matrix(rho$lon, ncol=nlat, nrow=nlon)

  psi$LAT = matrix(psi$lat, ncol=nlat+1, nrow=nlon+1, byrow=TRUE)
  psi$LON = matrix(psi$lon, ncol=nlat+1, nrow=nlon+1)

  area = (111*dy)*(111*dx*cospi(LAT/180))

  output = list(longitude=lons.rho, latitude=lats.rho, psi=psi,
                LON=LON, LAT=LAT, area=area,
                df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))

  return(output)

}

.grid_offset = function(n=0, dx, dy) {
  if(n==0) return(NULL)
  n = 2*n + 1
  x = seq(from=0, to=dx, length=n) - 0.5*dx
  y = seq(from=0, to=dy, length=n) - 0.5*dy
  return(expand.grid(y=y, x=x)[, 2:1])
}

.getBreaks = function(x) {
  out = c(x[1] - 0.5*(diff(x[1:2])),
          head(x, -1) + 0.5*diff(x),
          tail(x, 1) + 0.5*(diff(tail(x, 2))))
  return(out)
}

interp_grid = function(grid) {

  .interpx = function(y) {
    if(all(!is.na(y))) return(y)
    if(sum(!is.na(y)) < 4) {
      warning("Not enough points to do something.")
      return(y)
    }
    dy = diff(y)
    x = seq_along(y)
    dat = data.frame(x=x, y=y, dy=c(NA,dy))
    ff = gam(dy ~ s(x), data=dat, family=gaussian(link="log"))
    dat$pred = predict(ff, type="response", newdata = dat)
    dat$use = dat$dy
    dat$use[is.na(dat$use)] = dat$pred[is.na(dat$use)]
    dat$use[1] = 0
    dat$cuse = cumsum(dat$use)
    ind = which(!is.na(dat$y))[1]
    dat$cuse = dat$cuse - dat$cuse[ind]
    dat$yx = dat$cuse + dat$y[ind]
    dat$yf = dat$y
    dat$yf[is.na(dat$yf)] = dat$yx[is.na(dat$yf)]
    return(dat$yf)
  }

  if(any(is.na(grid$LON)))
    if(length(dim(grid$LON))==2) grid$LON = apply(grid$LON, 2, .interpx)
  if(any(is.na(grid$LAT)))
    if(length(dim(grid$LAT))==2) grid$LAT = t(apply(grid$LAT, 1, .interpx))

  return(grid)

}

.mask_correction = function(x, mask) {
  if(is.null(mask)) return(x)
  mask[mask==0] = NA
  x = x*as.numeric(mask)
  return(x)
}

.grid_info = function(grid) {
  dim = list()
  units = NULL
  dim.units = NULL
  ovarid = NULL
  if(length(dim(grid$longitude))==2) {
    dim$i = seq_len(nrow(grid$longitude))
    ovarid = c(ovarid, "longitude")
    units = c(units, longitude="degrees East")
    dim.units = c(dim.units, i="index")
  } else {
    dim$longitude = grid$longitude
    dim.units = c(dim.units, longitude="degrees East")
  }

  if(length(dim(grid$latitude))==2) {
    dim$j = seq_len(ncol(grid$latitude))
    ovarid = c(ovarid, "latitude")
    units = c(units, latitude="degrees North")
    dim.units = c(dim.units, j="index")
  } else {
    dim$latitude = grid$latitude
    dim.units = c(dim.units, latitude="degrees North")
  }

  ovarid = c(ovarid, "mask", "area")
  # TO_DO: store area units in the grid file!
  units = c(units, mask="0=land/1=ocean", area="km2")
  # end of info to write_ncdf.grid

  out = list(varid="mask", units=units, dim = dim, dim.units=dim.units,
             ovarid=ovarid, var=ovarid)

  return(out)
}

.is_regular_grid = function(x) {

  if(inherits(x, "gts") | inherits(x, "static")) x = x$grid
  almost.zero = function(x) (all(diff(x) < 1e-6))

  reg_lon = all(apply(x$LON, 1, almost.zero))
  reg_lat = all(apply(x$LAT, 2, almost.zero))

  is_regular =  reg_lon & reg_lat

  return(is_regular)

}

.compare_grids = function(e1, e2, strict=TRUE) {

  if(inherits(e1, "gts") | inherits(e1, "static")) e1=e1$grid
  if(inherits(e2, "gts") | inherits(e2, "static")) e2=e2$grid

  if(isTRUE(strict)) return(identical(e1, e2))
  e1$mask = e2$mask = NULL
  e1$area = e2$area = NULL
  e1$rho = e2$rho = NULL
  e1$psi = e2$psi = NULL
  e1$df = e2$df = NULL
  e1$info = e2$info = NULL
  return(identical(e1, e2))
}
