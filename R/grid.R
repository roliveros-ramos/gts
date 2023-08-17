
#' Make a rectangular grid
#'
#' @param lon Longitude of the left and right corners of the grid.
#' @param lat Latitude of the upper and lower corners of the grid.
#' @param dx Resolution in the longitude axis (in degrees).
#' @param dy Resolution in the latitude axis (in degrees).
#' @param n Number of points within a square to test ocean/land mask. See details.
#' @param hires Use high resolution coast line? Default to FALSE, TRUE will take considerly longer time to compute.
#'
#' @return A list, including the grid and land/ocean mask.
#' @export
#'
#' @examples
#' grid = make_grid(lon=c(2, 9), lat=c(40,45), dx=1/12, dy=1/12, n=1)
#' plot(grid)
make_grid = function(lon, lat, dx, dy=dx, n=1, thr=0.8, hires=FALSE, mask=TRUE) {

  lon = range(lon, na.rm=TRUE)
  lat = range(lat, na.rm=TRUE)

  grid = .create_grid(lon=lon, lat=lat, dx=dx, dy=dy)
  if(isTRUE(mask)) {
    mm = .create_mask(lon=lon, lat=lat, dx=dx, dy=dy, n=n, thr=thr, hires=hires)
    grid$mask  = mm$mask
    grid$prob  = mm$ocean
    grid$n     = mm$n
    grid$hires = hires
  } else {
    grid$mask  = NULL
    grid$prob  = NULL
    grid$n     = 0
    grid$hires = FALSE
  }

  # info to write_ncdf.grid
  grid$info = .grid_info(grid)

  class(grid) = c("grid", class(grid))
  return(grid)

}

#' Update grid
#'
#' @param grid An object of class 'grid'
#' @param thr A value between 0 and 1, used as threshold to define an ocean cell.
#'
#' @return A grid object.
#' @export
#'
update_grid = function(grid, thr) {
  if(!inherits(grid, "grid")) stop("This is not a valid grid object.")
  grid$mask = 0 + (grid$prob >= thr)
  return(grid)
}


#' Read the grid from a file
#'
#' @param mask Variable with an actual mask for the grid, must be of dimension 2.
#' @param varid Variable used to create the mask, see details.
#' @param create_mask Should we create a new mask using the same algorithm as \link{make_grid}?
#' @param ... Additional arguments passed to the \link{make_grid} function, mainly 'n', 'thr' and 'hires'.
#'
#' @inheritParams read_gts
#' @return Return the grid associated to the variables in the file.
#' @export
#'
#' @examples read_grid()
read_grid = function(filename, varid=NULL, mask=NULL, create_mask=FALSE, ...) {

  control = list(...)

  nc = nc_open(filename = filename)
  on.exit(nc_close(nc))

  if(!is.null(mask) & !is.null(varid))
    warning("Ignoring 'varid', using 'mask' for calculations.")
  if(!is.null(mask)) varid = mask

  if(is.null(varid)) {
    varid = names(nc$var)
    ndims = sapply(nc$var, FUN="[[", i="ndims")
    ind   = which(ndims>=2)
    if(length(ind)==0) stop("No suitable variables found in file.")
    if(length(ind)!=1) {
      ind = ind[1]
    }
    varid = varid[ind]
  }

  ndimx = nc$var[[varid]]$varsize
  hasdepth = if(length(ndimx)==4) TRUE else FALSE

  dimx = ncvar_dim(nc, varid=varid, value=TRUE)
  dims = sapply(dimx, length)

  is_monot = all(sapply(dimx, FUN=is_monotonic))
  if(!is_monot) stop("Variable dimensions have non-monotonic values.")

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

  depth_conf = list(depth_unit = NULL, depth=NULL)
  depth_name = NULL
  if(hasdepth) {
    tmp = ncatt_get(nc, varid=names(dimx)[3], attname = "units")
    dunit = if(tmp$hasatt) tmp$value else NULL
    depth_conf = list(depth_unit=dunit, depth=dimx[[3]])
    depth_name = "depth"
  }

  breaks = lapply(dimx, FUN=.getBreaks)

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

  odim = if(hasdepth) list(nlon, nlat, dimx[[3]]) else list(nlon, nlat)
  names(odim) = c(lon_name, lat_name, depth_name)

  dim.units = c(dlon_unit, dlat_unit, depth_conf$depth_unit)
  names(dim.units) = c(lon_name, lat_name, depth_name)

  grid = list(longitude=longitude, latitude=latitude, rho=list(LON=LON, LAT=LAT),
              psi=list(LON=pLON, LAT=pLAT), LON=LON, LAT=LAT, area=NULL, mask=NULL,
              df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))
  class(grid) =  c("grid", class(grid))


  grid = fill(grid, control=list(create_psi=.is_regular_grid(grid)))
  grid$area = suppressMessages(area(grid))

  # creating mask
  if(isTRUE(create_mask)) {
    if(is.null(control$hires)) control$hires = FALSE
    mm = .create_mask(grid=grid, n=control$n, thr=control$thr, hires=control$hires)
    grid$mask  = mm$mask
    grid$prob  = mm$ocean
    grid$n     = mm$n
    grid$hires = control$hires
  } else {
    # we will modify the data so it ends being a matrix with 1s and NAs.
    .mask = function(x) 0 + !all(is.na(x))
    x = drop(x)
    if(!is.null(mask)) {
      if(length(dim(x))!=2) stop("'mask' must be a variable of dimension 2.")
    } else {
      if(length(dim(x))>4) stop("Variables of dimensions greater than 4 are not supported.")
      if(length(dim(x))>=3) x = apply(x, head(seq_along(dim(x)), -1), FUN=.mask)
      if(length(dim(x))==3) x = x[,,1] # after sorting depth, this is 'surface'
      if(length(dim(x))<2) stop(sprintf("Variable '%s' has dimension lower than 2.", varid))
    }
    x = x/x # 0s become NaN, so NA.
    x[is.na(x)] = NA

    grid$mask  = x
    grid$prob  = NULL
    grid$n     = NULL
    grid$hires = NULL

  }

  # info to write_ncdf.grid
  grid$info = .grid_info(grid)

  return(grid)

}

# Methods -----------------------------------------------------------------

#' @exportS3Method plot grid
plot.grid = function(x, land.col="darkolivegreen4", sea.col="aliceblue", prob=FALSE,
                     boundaries.col="black", grid=TRUE, grid.lwd=1, grid.col="lightgray", lwd=2, ...) {

  if(is.null(x$mask)) {
    z = x$LAT
    col = sea.col
    zlim = range(z, na.rm=TRUE)
    } else {
    if(isTRUE(prob)) {
      z = x$prob
      col = divergencePalette(n=x$n, zlim=c(0,1), col=c(land.col, sea.col), center = 0.5, p=0.4)
      zlim = c(0, 1)
    } else {
      z = x$mask
      col = c(land.col, sea.col)
      zlim = c(0,1)
      }
  }

  hires = if(!is.null(x$hires)) x$hires else FALSE

  image.map(x$LON, x$LAT, z, land=FALSE, legend=prob, zlim=zlim,
            col=col, grid=grid, grid.lwd=grid.lwd, grid.col=grid.col)
  map_details(fill=FALSE, col=boundaries.col, lwd=lwd, hires=hires, ...)
  return(invisible(NULL))

}

# Auxiliar functions ------------------------------------------------------

#' Are the points in the ocean (or land)?
#'
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes.
#' @param hires Boolean. Use high resolution mask to compute ocean/land?
#'
#' @return A logical vector. For is_ocean, TRUE when the point is on water, FALSE for land or NA. Similiar for is_land.
#' @export
#'
is_ocean = function(lon, lat, hires=FALSE) {

  coords = cbind(lon=lon, lat=lat)
  out = rep(FALSE, length=nrow(coords))
  out[which(!complete.cases(coords))] = NA
  xind = which(complete.cases(coords))
  coords = coords[xind, ]

  layer = if(isTRUE(hires)) land_hr else land_lr

  crs = suppressWarnings(CRS(proj4string(layer)))
  sets = SpatialPoints(cbind(lon=coords[, "lon"], lat=coords[, "lat"]),
                       proj4string=crs)
  ind = xind[which(is.na(over(sets, layer)))]
  out[ind] = TRUE
  return(out)
}

#' @rdname is_ocean
#' @export
is_land = function(lon, lat, hires=FALSE) return(!is_ocean(lon, lat, hires))


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
    ff = gam(dy ~ s(x), dat=dat, family=gaussian(link="log"))
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

  if(inherits(x, "gts")) x = x$grid
  almost.zero = function(x) (all(diff(x) < 1e-6))

  reg_lon = all(apply(x$LON, 1, almost.zero))
  reg_lat = all(apply(x$LAT, 2, almost.zero))

  is_regular =  reg_lon & reg_lat

  return(is_regular)

}
