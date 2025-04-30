
#' Read an static variable from a netCDF file
#'
#' @param varid Variable to read static variable (time-invariant).
#' @param ... Additional arguments, currently unused.
#'
#' @returns
#' @export
#' @inheritParams read_grid
#' @inheritParams read_gts
#' @examples
read_static = function(filename, varid=NULL, create_mask=FALSE, ...) {

  control = list(...)

  nc = nc_open(filename = filename)
  on.exit(nc_close(nc))

  # if(!is.null(mask) & !is.null(varid))
  #   warning("Ignoring 'varid', using 'mask' for calculations.")
  # if(!is.null(mask)) varid = mask

  # find the variable to read
  if(is.null(varid)) {
    varid = names(nc$var)
    ndims = sapply(nc$var, FUN="[[", i="ndims")
    ind   = which(ndims>1)
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

  ndimx = nc$var[[varid]]$varsize
  hasdepth = if(length(ndimx)==3) TRUE else FALSE

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
    x = if(hasdepth) x[ind, , , drop=FALSE] else x[ind, , drop=FALSE]
  }

  if(is_decreasing(dimx[[2]])) {
    warning("Second dimension has decreasing values: fixing it.")
    dimx[[2]] = rev(dimx[[2]])
    ind = rev(seq_along(dimx[[2]]))
    x = if(hasdepth) x[, ind, , drop=FALSE] else x[, ind, drop=FALSE]
  }

  if(hasdepth) {
    if(is_decreasing(dimx[[3]])) {
      warning("Third dimension has decreasing values: fixing it.")
      dimx[[3]] = rev(dimx[[3]])
      ind = rev(seq_along(dimx[[3]]))
      x = x[, , ind, drop=FALSE]
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
    xx = drop(x)
    if(!is.null(mask)) {
      if(length(dim(xx))!=2) stop("'mask' must be a variable of dimension 2.")
    } else {
      if(length(dim(xx))>3) stop("Variables of dimensions greater than 4 are not supported.")
      if(length(dim(xx))==3) xx = xx[,,1] # after sorting depth, this is 'surface'
      if(length(dim(xx))<2) stop(sprintf("Variable '%s' has dimension lower than 2.", varid))
    }
    xx = xx/xx # 0s become NaN, so NA.
    xx[is.na(xx)] = NA

    grid$mask  = xx
    grid$prob  = NULL
    grid$n     = NULL
    grid$hires = NULL

  }

  # info to write_ncdf.grid
  grid$info = .grid_info(grid)

  out = list(x=x, grid=grid)
  class(out) = c("static")

  return(out)

}

