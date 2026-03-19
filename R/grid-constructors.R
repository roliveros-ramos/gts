#' Create or read a spatial grid object
#'
#' `make_grid()` creates a regular rectangular grid from horizontal coordinate
#' ranges or vectors.
#'
#' `read_grid()` reads the horizontal geometry associated with a variable in a
#' netCDF file and returns it as a `grid` object.
#'
#' `update_grid()` rebuilds the binary spatial mask of an existing `grid` object
#' from its probability surface.
#'
#' Grid objects store horizontal coordinates, coordinate matrices, optional cell
#' corner coordinates, optional cell area, and optional mask information used by
#' downstream gridded objects.
#'
#' @param lon,lat Horizontal coordinate ranges or vectors used to define a
#'   regular grid in `make_grid()`. If vectors longer than two are supplied and
#'   `dx` or `dy` is omitted, the corresponding resolution is inferred from the
#'   median spacing.
#' @param dx,dy Horizontal grid resolution in degrees for longitude and latitude.
#' @param n Number of test points per grid cell used when computing a
#'   coastline-based mask.
#' @param thr Threshold between `0` and `1` used to convert the ocean-proportion
#'   surface to a binary spatial mask.
#' @param hires Logical; if `TRUE`, use the high-resolution coastline when
#'   computing a coastline-based mask.
#' @param mask Logical in `make_grid()`, indicating whether a coastline-based
#'   mask should be created, or a character variable name in `read_grid()`,
#'   indicating which netCDF variable should be used to derive the mask.
#' @param grid A `grid` object to update.
#' @param filename Path to a netCDF file.
#' @param varid Variable identifier. In `read_grid()`, this is the variable name
#'   used to infer the grid and, when `create_mask = FALSE`, to derive a spatial
#'   mask from data availability.
#' @param create_mask Logical; if `TRUE`, create a coastline-based mask using the
#'   same workflow as [make_grid()]. If `FALSE`, derive the mask from the
#'   selected variable or from `mask` when supplied.
#' @param value Logical; if `TRUE`, `read_grid()` returns a list with both the
#'   derived mask-like array (`x`) and the `grid` object. If `FALSE`, it returns
#'   only the `grid` object.
#' @param ... Additional control arguments used by `read_grid()` during
#'   coastline-based mask creation, notably `n`, `thr`, and `hires`.
#'
#' @details
#' `make_grid()` creates a regular rectangular grid from the range of `lon` and
#' `lat`. If `mask = TRUE`, it also creates a coastline-based spatial mask and
#' stores the underlying ocean-proportion surface in `prob`.
#'
#' `update_grid()` updates `grid$mask` by thresholding `grid$prob` at the value
#' given by `thr`. It therefore requires a valid `grid` object with a probability
#' surface already present.
#'
#' If `varid` is not supplied, `read_grid()` searches the netCDF file for a
#' suitable variable, excluding variables whose names begin with `lon` or `lat`.
#' It first prefers variables with more than two dimensions, then falls back to
#' two-dimensional variables when needed.
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
#' the grid also stores corresponding coordinate matrices and cell-corner
#' coordinates. For irregular grids, the geometry is represented by coordinate
#' matrices and the first two dimensions are indexed as `i` and `j`.
#'
#' `read_grid()` then fills missing auxiliary grid information, computes cell
#' area, and creates the spatial mask in one of two ways:
#' \itemize{
#'   \item if `create_mask = TRUE`, it creates a coastline-based mask using the
#'   same helper as [make_grid()];
#'   \item otherwise, it derives the mask from the availability of the selected
#'   variable, or from the variable named by `mask` when supplied.
#' }
#'
#' @return
#' Depending on the function:
#' \describe{
#'   \item{`make_grid()`}{A `grid` object containing horizontal coordinates,
#'   coordinate matrices, optional cell-corner coordinates, optional spatial
#'   mask, optional probability surface, and grid metadata.}
#'   \item{`update_grid()`}{The updated `grid` object.}
#'   \item{`read_grid()`}{By default, a `grid` object containing the horizontal
#'   geometry associated with the selected netCDF variable, plus area, mask, and
#'   metadata. If `value = TRUE`, a list with components `x` and `grid`, where
#'   `x` is the intermediate mask-like array derived from the selected variable.}
#' }
#'
#' @seealso [grid-class], [make_grid()], [read_static()], [read_gts()]
#'
#' @examples
#' \dontrun{
#' grd <- make_grid(lon = c(2, 9), lat = c(40, 45), dx = 1/12, dy = 1/12)
#' grd2 <- read_grid("temperature.nc", varid = "thetao")
#' grd3 <- update_grid(grd2, thr = 0.8)
#' }
#' @name grid
#' @export
make_grid = function(lon, lat, dx, dy=dx, n=1, thr=0.8, hires=FALSE, mask=TRUE) {

  if(length(lon)>2 & missing(dx)) {
    dlon = diff(lon)
    check = all((unique(dlon) - median(dlon)) < 1e-6)
    dx = median(dlon)
  }

  if(length(lat)>2 & missing(dy)) {
    dlat = diff(lat)
    check = all((unique(dlat) - median(dlat)) < 1e-6)
    dy = median(dlat)
  }

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

#' @rdname grid
#' @export
update_grid = function(grid, thr) {
  if(!inherits(grid, "grid")) stop("This is not a valid grid object.")
  grid$mask = 0 + (grid$prob >= thr)
  return(grid)
}


#' @rdname grid
#' @export
read_grid = function(filename, varid=NULL, mask=NULL, create_mask=FALSE,
                     value=FALSE, ...) {

  control = list(...)

  nc = nc_open(filename = filename)
  on.exit(nc_close(nc))

  if(!is.null(mask) & !is.null(varid))
    warning("Ignoring 'varid', using 'mask' for calculations.")
  if(!is.null(mask)) varid = mask

  if(is.null(varid)) {
    varid = names(nc$var)
    ndims = sapply(nc$var, FUN="[[", i="ndims")
    # exclude latitude and longitude
    ind_lon = grep(x=varid, pattern="^lon", ignore.case=TRUE)
    ind_lat = grep(x=varid, pattern="^lat", ignore.case=TRUE)
    ind_exc = c(ind_lon, ind_lat)
    varid = varid[-ind_exc]
    ndims = ndims[-ind_exc]
    # first checks for arrays of dim 3 or more.
    ind   = which(ndims>2)
    # If none is found, check for dim 2. Sometimes file has just one layer...
    if(length(ind)==0) ind = which(ndims>=2)
    # If nothing found, throw an error
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

  if(isTRUE(value)) return(list(x=x, grid=grid))

  return(grid)

}
