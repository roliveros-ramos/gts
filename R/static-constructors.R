#' Read a static gridded object from a netCDF file
#'
#' `read_static()` reads a time-invariant gridded variable from a netCDF file
#' and returns it as a `static` object.
#'
#' A `static` object stores a spatially structured field without a time axis,
#' together with its horizontal coordinates, grid geometry, dimension metadata,
#' and variable metadata.
#'
#' @param filename Path to a netCDF file.
#' @param varid Variable identifier. This is the variable name in the netCDF
#'   file. If omitted, the function tries to identify a single suitable variable
#'   with more than one dimension.
#' @param create_mask Logical; if `TRUE`, create a coastline-based spatial mask
#'   using the same masking workflow as [make_grid()]. If `FALSE`, derive the
#'   mask from the availability of the selected variable.
#' @param ... Additional control arguments used during variable selection and
#'   optional mask creation.
#'
#'   Supported entries include:
#'   \describe{
#'   \item{`use_first`}{Logical; if `TRUE`, use the first suitable netCDF
#'   variable when several candidates are present.}
#'   \item{`n`}{Number of test points per grid cell used when creating a
#'   coastline-based mask.}
#'   \item{`thr`}{Threshold used to convert the ocean-proportion surface to a
#'   binary mask when creating a coastline-based mask.}
#'   \item{`hires`}{Logical; if `TRUE`, use the high-resolution coastline when
#'   creating a mask.}
#'   }
#'
#' @details
#' If `varid` is not supplied, `read_static()` searches the netCDF file for a
#' suitable variable with more than one dimension. If no such variable exists,
#' it errors. If several candidates exist, it errors unless `use_first = TRUE`,
#' in which case the first candidate is used.
#'
#' The current implementation is intended for two-dimensional or
#' three-dimensional variables. The first two dimensions are interpreted as
#' horizontal space. When a third dimension is present, it is treated as depth
#' and its metadata are stored in `info$depth`.
#'
#' The function checks that all coordinate dimensions are monotonic. If longitude,
#' latitude, or depth are stored in decreasing order, both the coordinates and
#' the data array are reversed to restore increasing order.
#'
#' Longitude and latitude are read from variables whose names begin with `lon`
#' and `lat` when available. Otherwise, the first two dimensions of the selected
#' variable are used as horizontal coordinates.
#'
#' For regular grids, longitude and latitude are stored as numeric vectors and
#' the associated grid geometry is built from those vectors. For irregular grids,
#' coordinate geometry is stored as matrices and the first two dimensions are
#' indexed as `i` and `j`.
#'
#' When `create_mask = TRUE`, a coastline-based mask is generated with the same
#' helper used by [make_grid()]. When `create_mask = FALSE`, the mask is derived
#' from the availability of the selected variable, so cells containing only
#' missing values are excluded from the spatial mask.
#'
#' The returned object includes a `grid` component with horizontal geometry and
#' mask metadata, plus an `info` list describing the variable, dimensions,
#' units, and global attributes.
#'
#' @return A `static` object. This is a list-like object with components:
#' \describe{
#'   \item{`x`}{Numeric matrix or array containing the spatial field.}
#'   \item{`grid`}{Associated `grid` object describing the horizontal geometry
#'   and mask.}
#'   \item{`longitude`, `latitude`}{Horizontal coordinates as vectors or
#'   matrices.}
#'   \item{`breaks`}{List of cell boundaries for each stored dimension.}
#'   \item{`info`}{Metadata list containing the variable identifier, long name,
#'   depth metadata when present, dimension values, dimension units, variable
#'   units, original variable names, and netCDF global attributes.}
#' }
#'
#' @seealso [static-class], [read_gts()], [make_grid()], [read_grid()],
#'   [grid-class]
#'
#' @examples
#' \dontrun{
#' bathy <- read_static("bathymetry.nc", varid = "depth")
#' dist  <- read_static("distance_to_coast.nc")
#' }
#' @name static
#' @export
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

  # get global attributes, units and long_name
  gatt = ncatt_get(nc, varid=0)

  tmp = ncatt_get(nc, varid=varid, attname = "units")
  var_unit = if(tmp$hasatt) tmp$value else ""

  tmp = ncatt_get(nc, varid=varid, attname = "long_name")
  long_name = if(tmp$hasatt) tmp$value else varid

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

  ovar = c(lon_var, lat_var, "x")
  ovarid = c(lon_var, lat_var, varid)

  dim.units = c(dlon_unit, dlat_unit, depth_conf$depth_unit)
  names(dim.units) = c(lon_name, lat_name, depth_name)

  units = c(lon_unit, lat_unit, var_unit)

  # create grid
  grid = list(longitude=longitude, latitude=latitude, rho=list(LON=LON, LAT=LAT),
              psi=list(LON=pLON, LAT=pLAT), LON=LON, LAT=LAT, area=NULL, mask=NULL,
              df=data.frame(lon=as.numeric(LON), lat=as.numeric(LAT)))
  class(grid) =  c("grid", class(grid))

  #grid = fill(grid, control=list(create_psi=.is_regular_grid(grid)))
  #grid$area = suppressMessages(area(grid))

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
  out$longitude = longitude
  out$latitude = latitude
  out$breaks = breaks

  out$info  = list(varid=varid, long_name=long_name,
                   depth=depth_conf,
                   dim = odim, var=ovar, dim.units=dim.units,
                   units=units, ovarid=ovarid, global=gatt)

  class(out) = c("static")

  return(out)

}

