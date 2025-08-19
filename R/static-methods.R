
# S3 methods --------------------------------------------------------------

#' @exportS3Method is.na static
is.na.static = is.na.gts

#' @export
dim.static = dim.gts

#' @export
Math.static = function(x, ...) {
  if(.Generic %in% c("cummax", "cummin", "cumprod", "cumsum"))
    stop(gettextf("'%s' not defined for \"static\" objects",
                  .Generic), domain = NA)
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}

#' @export
log.static = log.gts

#' @export
Math2.static = Math2.gts

#' @export
chooseOpsMethod.static = function(x, y, mx, my, cl, reverse) {
  if(inherits(y, "gts")) return(FALSE)
  return(TRUE)
}

#' @export
Ops.static = function(e1, e2) {

  # One argument Ops
  if(missing(e2)) {
    e1$x = get(.Generic, mode="function")(e1$x)
    return(e1)
  }

  # Ops.static only deals with static + static or static + default.
  # Operations involving a gts object are handled with Ops.gts

  # If e2 is not static, we deal with the case it is a matrix and move on.
  if(inherits(e1, "static") & !inherits(e2, "static")) {
    if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }

  # Now, if both e1 and e2 are static objects.

  # We must ensure:
  # - The dimensions are compatible
  ok0 = length(dim(e1)) == length(dim(e2))
  if(!ok0) stop(sprintf("Dimensions of both objects are not compatible (%d!=%d).",
                        length(dim(e1)), length(dim(e2))))
  # - The grids are the same
  ok1 = .compare_grids(e1$grid, e2$grid, strict=FALSE)
  if(!ok1) stop("Grids are not compatible.")

  # Finally, we make the operation over the arrays.
  e1$x = get(.Generic, mode="function")(e1$x, e2$x)
  # and we return e1.
  return(e1)

}

#' @export
Summary.static = Summary.gts

#' @export
print.static = function(x, ...) {
  cat(sprintf("Static gridded object: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  print(resolution(x))
  return(invisible(NULL))
}

#' @export
names.static = names.gts

#' @export
'names<-.static' = setnames_gts

#' @export
units.static = units.gts

# generic implemented in base.
#' @export
'units<-.static' = setunits_gts

#' @export
str.static = str.gts

# #' @exportS3Method reshape2::melt
#' @export
melt.static = function(data, ..., na.rm=FALSE, value.name=NULL) {

  if(is.null(data$depth)) {
    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     value=as.numeric(data$x))

    names(out)[3] = names(data)[1]

  } else {

    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     depth = rep(data$depth, each=nrow(data$grid$df)))
    out$value = as.numeric(data$x)

    names(out)[4] = names(data)[1]

  }

  return(out)

}

#' @export
mean.static = function(x, trim=0, na.rm=FALSE, ...) {
  return(mean(as.double(x), trim=trim, na.rm=na.rm, ...))
}

#' @export
sum.gts = function(x, na.rm=FALSE, ...) {
  return(sum(as.double(x), na.rm=na.rm, ...))
}

#' @export
as.double.static = as.double.gts

#' @export
sd.static = function(x, na.rm=FALSE, ...) {
  return(sd(as.double(x), na.rm=na.rm, ...))
}

#' @export
drop.static = function(x, ...) {
  ndim = dim(x)
  if(length(ndim)==2) return(x)
  x$x = drop(x$x)
  # check grid and attributes? do it at read time?
  return(x)
}

