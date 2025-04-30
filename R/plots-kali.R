
#' Draw an image plot with geographical coordinates and coastlines.
#'
#' @param lon x-axis coordinates, assumed to be longitude (degrees)
#' @param lat y-axis coordinates, assumed to be latitude (degrees)
#' @param z A matrix of dimension [lon x lat]
#' @param zlim numeric vector of length 2, giving the z matrix plotted range.
#' @param legend Boolean, add color legend? FALSE can be used to construct multiple plots with a common legend.
#' @param hires Boolean, use high resolution coastline?
#' @param slice If z is an 3D array, slice indicates which layer will be used.
#' @param col Color table to use for image, by default rev(rainbow(nlevel/10, start = 0/6, end = 4/6))
#' @param axes Boolean, add axes to the plot? Default is TRUE.
#' @param land Boolean, add a polygon over the land? If FALSE, only coastline is added.
#' @param land.col Color of the land, if land is TRUE.
#' @param labels Boolean, add LATITUDE and LONGITUDE labels?
#' @param sea.col Color of the water bodies.
#' @param boundaries.col Color of the country boundaries.
#' @param grid Boolean, add grid lines over the map?
#' @param grid.col Color of the grid.
#' @param grid.lwd Width of the grid lines.
#' @param ... Additional parameters passed to plot methods.
#'
#' @returns
#' @export
#' @inheritParams fields::image.plot
#' @examples
image.map = function (lon, lat, z, zlim=NULL, legend=TRUE, hires=FALSE, add = FALSE, nlevel = 1000, horizontal = FALSE,
                      legend.shrink = 0.9, legend.width = 1.2, slice=NULL,
                      legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL, graphics.reset = FALSE,
                      bigplot = NULL, smallplot = NULL, legend.only = FALSE,
                      col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)),
                      lab.breaks = NULL, axis.args = NULL, legend.args = NULL, axes=TRUE,
                      midpoint = FALSE, border = TRUE, lwd = 1, land=TRUE, land.col="black", labels = TRUE,
                      sea.col="white", boundaries.col="grey", grid=FALSE, grid.col="white", grid.lwd=0.5, ...) {

  if(!is.null(attr(z, "longitude"))) lon = attr(z, "lon")
  if(!is.null(attr(z, "latitude"))) lat = attr(z, "lat")

  if(is.null(zlim)) zlim = range(as.numeric(z), na.rm=TRUE)

  if(length(dim(z))==3) {
    if(!is.null(slice)) {
      z = z[, , slice]
    }
    if(is.null(slice)) stop("Must specify an slice for plotting arrays.")
  }

  if(length(dim(z))>3) stop("Only arrays of 3-dimensions are supported.")

  nx = nrow(z)
  ny = ncol(z)

  pm = .findPrimeMeridian(lon)

  if(!isTRUE(legend)) {
    .image.mapnl(lon=lon, lat=lat, z=z, zlim=zlim, hires=hires, add=add, nlevel=nlevel,
                 col=col, land=land, land.col=land.col, sea.col=sea.col, boundaries.col=boundaries.col,
                 grid.col=grid.col, grid=grid, axes=axes, border=border, labels=labels, grid.lwd=grid.lwd, ...)
    return(invisible())
  }

  old.par = par(no.readonly = TRUE)
  info = imageplot.info(x=lon, y=lat, z=z, zlim=zlim, ...)
  if (add) {
    big.plot = old.par$plt
  }
  if (legend.only) {
    graphics.reset = TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar = ifelse(horizontal, 3.1, 5.1)
  }
  temp = imageplot.setup(add = add, legend.shrink = legend.shrink,
                         legend.width = legend.width, legend.mar = legend.mar,
                         horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot = temp$smallplot
  bigplot = temp$bigplot
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(x=lon, y=lat, z=z, add = add, col = col, axes=FALSE,
            xlab="", ylab="", zlim=zlim, xaxs="i", yaxs="i", ...)
      map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE,
                  axes=axes, border=border, boundaries.col=boundaries.col, nx=nx, ny=ny,
                  grid=grid, grid.col=grid.col, water=sea.col, land=land, grid.lwd=grid.lwd, labels=labels)
    }
    else {
      poly.image(x=lon, y=lat, z=z, add = add, col = col, midpoint = midpoint,
                 border = NA, lwd.poly = lwd, axes=FALSE, zlim=zlim, xaxs="i", yaxs="i",
                 xlab="", ylab="",...)
      map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE,
                  axes=axes, border=border, boundaries.col=boundaries.col, nx=nx, ny=ny,
                  grid=grid, grid.col=grid.col, water=sea.col, land=land, grid.lwd=grid.lwd, labels=labels)
    }
    big.par = par(no.readonly = TRUE)
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  ix = 1
  minz = info$zlim[1]
  maxz = info$zlim[2]
  binwidth = (maxz - minz)/nlevel
  midpoints = seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy = midpoints
  iz = matrix(iy, nrow = 1, ncol = length(iy))
  breaks = list(...)$breaks
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  if (!is.null(breaks) & !is.null(lab.breaks)) {
    axis.args = c(list(side = ifelse(horizontal, 1, 4),
                       mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
                       at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args = c(list(side = ifelse(horizontal, 1, 4),
                       mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                  axis.args)
  }
  if (!horizontal) {
    if (is.null(breaks)) {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }
  else {
    if (is.null(breaks)) {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col)
    }
    else {
      image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
            ylab = "", col = col, breaks = breaks)
    }
  }
  do.call("axis", axis.args)
  box()
  if (!is.null(legend.lab)) {
    legend.args = list(text = legend.lab, side = ifelse(horizontal,
                                                        1, 4), line = legend.mar - 2)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  mfg.save = par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    par(mar = par("mar"))
    invisible()
  }

  return(invisible())
}


#' Add details to a map
#'
#' @param primeMeridian Prime Meridian used in the map, possible values are "center" and "left", the former is the default.
#' @param water Color of the water masses.
#' @param countries Boolean, do add country boundaries? Default is FALSE.
#' @param ...
#'
#' @returns
#' @export
#' @inheritParams image.map
#' @inheritParams map_details
#' @inheritParams graphics::grid
#' @inheritParams maps::map
#' @examples
map_details = function(primeMeridian="center", hires=FALSE, col="darkolivegreen4",
                       interior=FALSE, axes=FALSE, border=FALSE, boundaries.col="black",
                       grid=FALSE, grid.col="white", grid.lwd=0.5, cex.axis=0.75, fill=TRUE,
                       boundary = TRUE, water=NULL, land=TRUE, countries=FALSE, nx=NULL, ny=nx, labels=TRUE, ...) {

  primeMeridian = match.arg(primeMeridian, choices=c("center", "left"))

  if(hires) {
    if(!requireNamespace("mapdata", quietly = TRUE)) {
      warning("You need to install the 'mapdata' package, using hires=FALSE.")
      hires = FALSE
    }
  }

  mapa =  if(hires) {
    if(primeMeridian=="center") "mapdata::worldHires" else "mapdata::world2Hires"
  } else {
    if(primeMeridian=="center") "world" else "world2"
  }

  # mapa = if(primeMeridian=="center") "world" else "world2"

  wrap = ifelse(primeMeridian=="center", c(-180,180), c(0,360))

  if(isTRUE(grid)) grid(nx=nx, ny=ny, col=grid.col, lty=1, lwd=grid.lwd)

  if(isTRUE(land)) {
    if(isTRUE(fill)) {
      map(database=mapa, fill = TRUE, col = col, add = TRUE, interior=FALSE,
          border=col, boundary = boundary, ...)
    } else {
      map(database=mapa, fill = FALSE, col = col, add = TRUE, interior=FALSE,
          boundary = boundary, ...)
    }

    if(!is.null(water)) {
      if(primeMeridian=="center") {
        map("lakes", fill = TRUE, col = water, add = TRUE, border=water, ...)
      } else {
        lakes = map("lakes", plot=FALSE)
        lakes$x = checkLongitude(lakes$x, primeMeridian = "left")
        map(lakes, fill = TRUE, col = water, add = TRUE, border=water, ...)
      }
    }
    if(isTRUE(countries))
      map(database=mapa, fill = FALSE, col = boundaries.col, add = TRUE)
  }

  if(axes) {
    map.axes2(cex.axis=cex.axis)
    if(isTRUE(labels)) {
      mtext("LONGITUDE", 1, line = 1.8*cex.axis/0.75, cex = 0.9*par("cex"))
      mtext("LATITUDE", 2, line = 2.4*cex.axis/0.75, cex = 0.9*par("cex"))
    }
  } else {
    if(border) box()
  }

  return(invisible())
}

#' Add geographical axes to a map
#'
#' @param sides Side of the axis, 1 and 2 by default. See \code{graphics::axis}.
#' @param cex.axis Character expansion factor for axis labels.
#'
#' @returns
#' @export
#' @inheritParams graphics::axis
#' @examples
map.axes2 = function(sides=c(1,2), cex.axis=0.75, line=-0.4) {

  .axis.map(sides[1], "lon", las=1, cex.axis=cex.axis, line=line, tick=FALSE)
  .axis.map(sides[2], "lat", las=1, cex.axis=cex.axis, line=line, tick=FALSE)
  axis(1, labels=FALSE)
  axis(2, labels=FALSE)
  axis(3, labels=FALSE)
  axis(4, labels=FALSE)
  box()

  return(invisible(NULL))
}


# Auxiliar functions ------------------------------------------------------


.image.mapnl = function(lon, lat, z, zlim, hires=FALSE, add = FALSE, nlevel=1000,
                        col = rev(rainbow(nlevel/10, start = 0/6, end = 4/6)), land=TRUE,
                        land.col="darkolivegreen4", sea.col="aliceblue", boundaries.col = "black",
                        grid.col="white", grid.lwd=0.5, grid=FALSE, axes=TRUE, border=!axes, labels=TRUE, ...) {

  pm = .findPrimeMeridian(as.numeric(lon))

  if(is.matrix(lon) & is.matrix(lat)) {
    poly.image(x=lon, y=lat, z=z, zlim=zlim, col=col, axes=FALSE, add=add, xlab="", ylab="", xaxs="i", yaxs="i", ...)
    map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE,
                axes=axes, border=border, boundaries.col=boundaries.col, nx=nrow(z), ny=ncol(z),
                grid=FALSE, grid.col=grid.col, grid.lwd=grid.lwd, water=sea.col, land=land, labels=labels)
  } else {
    image(x=lon, y=lat, z=z, zlim=zlim, col=col, axes=FALSE, add=add, xlab="", ylab="", xaxs="i", yaxs="i", ...)
    map_details(primeMeridian=pm, hires=hires,col=land.col, interior=FALSE,
                axes=axes, border=border, boundaries.col=boundaries.col, nx=nrow(z), ny=ncol(z),
                grid=grid, grid.col=grid.col, grid.lwd=grid.lwd, water=sea.col, land=land, labels=labels)

  }

  return(invisible())
}


.axis.map = function(side, type, usr=NULL, n=5, ...) {

  if(is.na(side)) return(invisible())

  old_usr = par("usr")
  on.exit(par(usr=old_usr))
  if(is.null(usr)) usr = old_usr
  par(usr=usr)

  is.x <- side%%2 == 1

  x = pretty(usr[if(is.x) 1:2 else 3:4], n=n)
  xc = checkLongitude(x, "center")

  axis(side=side, at=x, labels=coord2text(coord=xc, type=type), ...)

  return(invisible(x))
}

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

.coord2text = function(coord, type) {

  # write nicely coordinates
  degree = "\U00B0"
  hemi = if(coord==0) {
    rep(degree, 2)
  } else {
    if(coord>0) paste0(degree, c("N","E")) else paste0(degree, c("S","W"))
  }

  out = switch(type,
               lat = paste0(abs(coord), hemi[1]),
               lon = paste0(abs(coord), hemi[2]),
               as.character(coord)
  )

  return(out)
}

coord2text = function(coord, type) {

  if(!is.character(type)) type=deparse(substitute(type))
  out = sapply(coord, FUN=.coord2text, type=type)

  return(out)
}


