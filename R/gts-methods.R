
# S3 methods --------------------------------------------------------------

#' @exportS3Method is.na gts
is.na.gts = function(x) {
  x$x = is.na(x$x)
  return(x)
}

#' @export
dim.gts = function(x) dim(x$x)

#' @export
dim.grid = function(x) dim(x$LON)

#' @export
Math.gts = function(x, ...) {
  if(.Generic %in% c("cummax", "cummin", "cumprod", "cumsum"))
    stop(gettextf("'%s' not defined for \"gts\" objects",
                  .Generic), domain = NA)
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}

#' @export
log.gts = function(x, base = exp(1)) {
  x$x = log(x$x, base=base)
  return(x)
}

#' @export
Math2.gts = function(x, ...) {
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}


#' @export
Ops.gts = function(e1, e2) {
  if(inherits(e1, "gts") & !inherits(e2, "gts")) {
    if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }
  if(inherits(e2, "gts") & !inherits(e1, "gts")) {
    if(identical(dim(e2$x)[1:2], dim(e1))) e1 = as.numeric(e1)
    e2$x = get(.Generic, mode="function")(e1, e2$x)
    return(e2)
  }
  if(inherits(e1, "gts") & inherits(e2, "gts")) {
    ok1 = identical(e1$grid, e2$grid)
    ok2 = identical(dim(e1$x), dim(e2$x))
    if(ok1 & ok2) {
      e1$x = get(.Generic, mode="function")(e1$x, e2$x)
      return(e1)
    } else {
      stop(gettextf("Input 'gts' objects dimensions not compatible for '%s'.",
                    .Generic), domain = NA)
    }
  }
}

#' @export
Summary.gts = function(x, ..., na.rm = TRUE) {

  # is_gts = sapply(..., FUN=is, class2="gts")
  # if(any(is_gts)) {
  #   for(i in seq_along(...)) {
  #     if(is_gts[i]) ...[[i]] = ...[[i]]$x
  #   }
  # }
  get(.Generic, mode="function")(x$x, na.rm=na.rm)
}

#' @export
print.gts = function(x, ...) {
  cat(sprintf("Gridded Time Series: %s (%s)\n", x$info$varid, tail(x$info$units, 1)))
  cat(sprintf("Start     = %s\n", min(x$time)))
  cat(sprintf("End       = %s\n", max(x$time)))
  cat(sprintf("Frequency = %s\n", frequency(x)))
  # cat(sprintf("Spatial Dimension = [%s,%s]\n", dim(x)[1], dim(x)[2]))
  rlon = range(x$longitude)
  rlat = range(x$latitude)
  cat(sprintf("Longitude = [%0.2f, %0.2f]\n", rlon[1], rlon[2]))
  cat(sprintf("Latitude  = [%0.2f, %0.2f]\n", rlat[1], rlat[2]))
  return(invisible(NULL))
}
# S4 compatibility --------------------------------------------------------

setOldClass("gts")
setOldClass("grid")

# setMethod('Ops', signature(e1='gts', e2='ANY'), Ops.gts)
# setMethod('Arith', signature(e1='ANY', e2='gts'), Arith_gts)
# setMethod('Arith', signature(e1='gts', e2='ANY'), Ops.gts)

