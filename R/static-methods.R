
# S3 methods --------------------------------------------------------------

#' @exportS3Method is.na static
is.na.static = is.na.gts

#' @export
dim.static = dim.gts

#' @export
Math.gts = function(x, ...) {
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

# #' @export
# Ops.gts = function(e1, e2) {
#
#   # One argument Ops
#   if(missing(e2)) {
#     e1$x = get(.Generic, mode="function")(e1$x)
#     return(e1)
#   }
#
#   # For simplicity, just assume the e1 must be a gts object
#   if(!inherits(e1, "gts") & inherits(e2, "gts")) {
#     return(get(.Generic, mode="function")(e2, e1))
#   }
#
#   # If e2 is not gts, we deal with the case it is a matrix and move on.
#   if(inherits(e1, "gts") & !inherits(e2, "gts")) {
#     if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
#     e1$x = get(.Generic, mode="function")(e1$x, e2)
#     return(e1)
#   }
#
#   # Now, if both e1 and e2 are gts objects.
#
#   # We must ensure:
#   # - The dimensions are compatible
#   ok0 = length(dim(e1)) == length(dim(e2))
#   if(!ok0) stop(sprintf("Dimensions of both objects are not compatible (%d!=%d).",
#                         length(dim(e1)), length(dim(e2))))
#   # - The grids are the same
#   ok1 = .compare_grids(e1$grid, e2$grid, strict=FALSE)
#   if(!ok1) stop("Grids are not compatible.")
#   # - The time frequency is the same
#   ok2 = identical(frequency(e1),frequency(e2))
#   if(!ok2) stop("Time frequencies are not compatible.")
#   # TO DO: more checks for time consistency.
#
#   # If e1 is a climatology, move it to e2
#   if(is_climatology(e1) & !is_climatology(e2)) {
#     return(get(.Generic, mode="function")(e2, e1))
#   }
#   # Now, e1 can be a climatology too. But it doesn't matter, e1 leads.
#   # If e2 is a climatology, we make it match the e1 cycles.
#   if(is_climatology(e2)) {
#     if(length(dim(e2))==3) e2$x = e2$x[, , as.numeric(cycle(e1)), drop=FALSE]
#     if(length(dim(e2))==4) e2$x = e2$x[, , , as.numeric(cycle(e1)), drop=FALSE]
#   }
#   # After this, the dimensions of e2 must match e1.
#   ok3 = identical(dim(e1$x), dim(e2$x))
#   if(!ok3) stop(gettextf("Input 'gts' objects dimensions not compatible for '%s'.",
#                          .Generic), domain = NA)
#   # Finally, we make the operation over the arrays.
#   # Currently, we're not checking the time is identical, must do?
#   e1$x = get(.Generic, mode="function")(e1$x, e2$x)
#   # This is not needed, but just in case :)
#   e1$info$climatology = is_climatology(e1) & is_climatology(e2)
#   # and we return e1.
#   return(e1)
#
# }

#' @export
Summary.static = Summary.gts

#' @export
print.static = function(x, ...) {
  cat(sprintf("Static gridded object: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  return(invisible(NULL))
}

#' @export
names.static = names.gts

#' @export
'names<-.static' = function(x, value) {

  if(any(is.na(value))) stop("NAs are not allowed in names for a gts static object.")
  if(length(value)>2) stop("A maximum of two values (varid, long_name) must be provided.")

  x$info$varid = value[1]
  x$info$ovarid[length(x$info$ovarid)] = value[1]
  x$info$long_name = value[2]

  x

}

#' @export
units.static = units.gts

# generic implemented in base.
#' @export
'units<-.static' = function(x, value) {

  if(length(value)!=1) stop("Only one value for units must be provided.")

  if(is.na(value)) value = ""

  x$info$units[length(x$info$units)] = value

  x

}

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

