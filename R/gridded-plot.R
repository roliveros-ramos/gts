#' Plot gridded objects
#'
#' These methods provide quick visualisation of gridded objects.
#'
#' The page covers three related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension;
#'   \item `grid` objects, which store spatial geometry and optional mask
#'   information.
#' }
#'
#' `plot.gts()` displays a spatial map for a selected time slice, and for
#' four-dimensional objects also a selected depth slice.
#'
#' `plot.static()` displays a spatial map of a static field, optionally selecting
#' a depth slice when the underlying data array has three dimensions.
#'
#' `plot.grid()` displays the grid geometry itself, optionally showing the grid
#' mask or the mask probability surface when available.
#'
#' @param x A `gts`, `static`, or `grid` object, depending on the method.
#' @param time Integer index of the time slice to plot for a `gts` object.
#' @param depth Integer index of the depth slice to plot for a `gts` or `static`
#'   object when the underlying data array has a depth dimension.
#' @param land.col Colour used for land cells or low probabilities in
#'   `plot.grid()`.
#' @param sea.col Colour used for sea cells or high probabilities in
#'   `plot.grid()`.
#' @param prob Logical; if `TRUE` and probability information is available,
#'   `plot.grid()` displays `x$prob` instead of the binary mask.
#' @param boundaries.col Colour used to draw coastline or polygon boundaries in
#'   `plot.grid()`.
#' @param grid Logical or `NULL`; whether to draw grid lines in `plot.grid()`.
#'   If `NULL`, the current implementation chooses a default based on whether a
#'   mask or probability surface is being displayed.
#' @param grid.lwd Line width of grid lines in `plot.grid()`.
#' @param grid.col Colour of grid lines in `plot.grid()`.
#' @param hires Logical; request high-resolution coastlines in `plot.grid()`.
#'   In the current implementation, this value is overridden by `x$hires` when
#'   that component is present.
#' @param lwd Line width used for coastline or boundary overlays in `plot.grid()`.
#' @param ... Additional graphical arguments. For `plot.gts()` and
#'   `plot.static()`, these are passed to `image_map()`. For `plot.grid()`, they
#'   are passed to `map_details()`.
#'
#' @details
#' `plot.gts()` converts logical data arrays to integer before plotting. It then
#' calls `image_map()` using the grid geometry stored in `x$grid$LON` and
#' `x$grid$LAT`.
#'
#' For four-dimensional `gts` objects, the selected `depth` and `time` slices
#' are plotted. For lower-dimensional `gts` objects, only the selected `time`
#' slice is used.
#'
#' Unless a `main` argument is supplied in `...`, `plot.gts()` constructs a
#' default title from the variable identifier and the selected slice labels. The
#' time label uses `x$info$time_var` when available, and `"time"` otherwise.
#'
#' `plot.static()` also converts logical arrays to integer before plotting. When
#' `x$grid$longitude` and `x$grid$latitude` are available they are used as the
#' plotting coordinates; otherwise `x$grid$LON` and `x$grid$LAT` are used.
#'
#' For three-dimensional `static` objects, the selected `depth` slice is
#' plotted. Otherwise the full two-dimensional field is plotted. Unless a `main`
#' argument is supplied, the default title is built from `x$info$varid` and, if
#' relevant, the selected depth.
#'
#' `plot.grid()` visualises the spatial grid itself. Its behaviour depends on
#' the components available in the object:
#' \itemize{
#'   \item if `x$mask` is `NULL`, a background field based on `x$LAT` is used to
#'   display the grid extent;
#'   \item if `prob = TRUE`, `x$prob` is displayed using a diverging palette and
#'   a legend is requested;
#'   \item otherwise the binary mask is displayed, with land and sea coloured by
#'   `land.col` and `sea.col`.
#' }
#'
#' For regular grids, `plot.grid()` passes longitude and latitude vectors to
#' `image_map()`. For irregular grids, it passes the coordinate matrices.
#' Coastline or boundary outlines are then added with `map_details()`.
#'
#' @return These methods are called for their side effect of drawing a plot and
#'   return `invisible(NULL)`.
#'
#' @seealso [gridded-reshape], [gts()], [read_gts()], [read_static()], [make_grid()]
#'
#' @examples
#' \dontrun{
#' plot(x, time = 1)
#' plot(bathy)
#' plot(grd, prob = TRUE)
#' }
#' @name gridded-plot
NULL

#' @describeIn gridded-plot Plot a spatial slice of a `gts` object.
#' @export
plot.gts = function(x, time=1, depth=1, ...) {
  args = list(...)
  if(mode(x$x)=="logical") x$x = 0L + x$x
  time_var = "time"
  if(!is.null(x$info$time_var)) time_var = x$info$time_var
  if(length(dim(x$x))==4) {
    image_map(x$grid$LON, x$grid$LAT, x$x[,,depth,], slice=time, ...)
    lab = sprintf("%s (depth=%s, %s=%s)", x$info$varid, x$depth[depth],
                  time_var, x$time[time])
  } else {
    image_map(x$grid$LON, x$grid$LAT, x$x, slice=time, ...)
    lab = sprintf("%s (%s=%s)", x$info$varid, time_var, x$time[time])
  }
  main = if(is.null(args$main)) lab else NULL
  title(main=main)
  return(invisible(NULL))
}

#' @describeIn gridded-plot Plot a spatial slice of a `static` object.
#' @export
plot.static = function(x, depth=1, ...) {
  args = list(...)
  if(mode(x$x)=="logical") x$x = 0L + x$x
  lon_check = !is.null(x$grid$longitude)
  lat_check = !is.null(x$grid$latitude)
  if(lon_check & lat_check) {
    ilon = x$grid$longitude
    ilat = x$grid$latitude
  } else {
    ilon = x$grid$LON
    ilat = x$grid$LAT
  }
  if(length(dim(x$x))==3) {
    image_map(ilon, ilat, x$x, slice=depth, ...)
    lab = sprintf("%s (depth=%s)", x$info$varid, x$depth[depth])
  } else {
    image_map(ilon, ilat, x$x, ...)
    lab = sprintf("%s", x$info$varid)
  }
  main = if(is.null(args$main)) lab else NULL
  title(main=main)
  return(invisible(NULL))
}


#' @describeIn gridded-plot Plot a `grid` object, optionally showing its mask or
#'   probability surface.
#' @exportS3Method plot grid
plot.grid = function(x, land.col="darkolivegreen4", sea.col="aliceblue", prob=FALSE,
                     boundaries.col="black", grid=NULL, grid.lwd=1, grid.col="lightgray",
                     hires=FALSE, lwd=2, ...) {

  is_regular = .is_regular_grid(x)

  if(is.null(x$mask)) {
    z = x$LAT
    col = sea.col
    zlim = range(z, na.rm=TRUE)
  } else {
    if(isTRUE(prob)) {
      z = x$prob
      col = divergencePalette(n=x$n, zlim=c(0,1), col=c(land.col, sea.col), center = 0.5, p=0.4)
      zlim = c(0, 1)
      if(is.null(grid)) grid=TRUE
    } else {
      z = x$mask
      z[is.na(z)] = 0
      col = c(land.col, sea.col)
      zlim = c(0, 1)
      if(is.null(grid)) grid=FALSE
    }
  }

  hires = if(!is.null(x$hires)) x$hires else FALSE

  if(is_regular) {
    lon = x$LON[, 1]
    lat = x$LAT[1, ]
  } else {
    lon = x$LON
    lat = x$LAT
  }

  image_map(lon, lat, z, land=FALSE, legend=prob, zlim=zlim, hires=hires,
            col=col, grid=grid, grid.lwd=grid.lwd, grid.col=grid.col)
  map_details(fill=FALSE, col=boundaries.col, lwd=lwd, hires=hires, ...)
  return(invisible(NULL))

}
