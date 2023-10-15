
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

  # One argument Ops
  if(missing(e2)) {
    e1$x = get(.Generic, mode="function")(e1$x)
    return(e1)
  }

  # For simplicity, just assumme the e1 must be a gts object
  if(!inherits(e1, "gts") & inherits(e2, "gts")) {
    return(get(.Generic, mode="function")(e2, e1))
  }

  # If e2 is not gts, we deal with the case it is a matrix and move on.
  if(inherits(e1, "gts") & !inherits(e2, "gts")) {
    if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }

  # Now, if both e1 and e2 are gts objects.

  # We must ensure:
  # - The dimensions are compatible
  ok0 = length(dim(e1)) == length(dim(e2))
  if(!ok0) stop(sprintf("Dimensions of both objects are not compatible (%d!=%d).",
                        length(dim(e1)), length(dim(e2))))
  # - The grids are the same
  ok1 = .compare_grids(e1$grid, e2$grid, strict=FALSE)
  if(!ok1) stop("Grids are not compatible.")
  # - The time frequency is the same
  ok2 = identical(frequency(e1),frequency(e2))
  if(!ok2) stop("Time frequencies are not compatible.")
  # TO DO: more checks for time consistency.

  # If e1 is a climatology, move it to e2
  if(is_climatology(e1) & !is_climatology(e2)) {
    return(get(.Generic, mode="function")(e2, e1))
  }
  # Now, e1 can be a climatology too. But it doesn't matter, e1 leads.
  # If e2 is a climatology, we make it match the e1 cycles.
  if(is_climatology(e2)) {
    if(length(dim(e2))==3) e2$x = e2$x[, , as.numeric(cycle(e1))]
    if(length(dim(e2))==4) e2$x = e2$x[, , , as.numeric(cycle(e1))]
  }
  # After this, the dimensions of e2 must match e1.
  ok3 = identical(dim(e1$x), dim(e2$x))
  if(!ok3) stop(gettextf("Input 'gts' objects dimensions not compatible for '%s'.",
                         .Generic), domain = NA)
  # Finally, we make the operation over the arrays.
  # Currently, we're not checking the time is identical, must do?
  e1$x = get(.Generic, mode="function")(e1$x, e2$x)
  # This is not needed, but just in case :)
  e1$info$climatology = is_climatology(e1) & is_climatology(e2)
  # and we return e1.
  return(e1)

}

#' @export
Summary.gts = function(x, ..., na.rm=TRUE) {
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

#' @export
c.gts = function(...) {

  out = list(...)

  .identical = function(x) all(sapply(x, identical, y=x[[1]]))

  # check time dimension, only 3D for now, extend.
  time_check = all(sapply(out, FUN=function(x) length(dim(x))) == 3)
  if(!time_check) stop("Only 3D gts objects are supported.")
  # check horizontal dimensions are compatible
  dims = lapply(lapply(out, dim), FUN="[", i=1:2)
  dim_check = .identical(dims)
  if(!dim_check) stop("Spatial dimensions are not compatible.")
  # check for identical grids.
  gLON = lapply(out, FUN=function(x) c(x$grid$LON))
  gLAT = lapply(out, FUN=function(x) c(x$grid$LAT))
  coord_check = .identical(gLON) & .identical(gLAT)
  if(!coord_check) stop("Grid coordinates do not match.")
  # check for time frequency
  freqs = sapply(out, frequency)
  freq_check = all(sapply(freqs, identical, freqs[1]))
  if(!freq_check) stop("Frequencies are not compatible.")
  frequency = freqs[1]

  starts = sapply(out, start.gts)
  starts[2,] = starts[2,]/frequency
  starts = colSums(starts)
  istart = which.min(starts)

  ends = sapply(out, end.gts)
  ends[2,] = ends[2,]/frequency
  ends = colSums(ends)
  iend = which.max(ends)

  start = start(out[[istart]])
  end   = end(out[[iend]])

  the_ts = ts(NA, start=start, end=end, frequency = frequency)
  the_ts = ts(seq_along(the_ts), start=start, end=end, frequency = frequency)

  the_array = array(NA, dim=c(dims[[1]], length(the_ts)))
  the_time  = array(NA, dim=c(length(the_ts)))
  info_time  = array(NA, dim=c(length(the_ts)))

  starts = sapply(out, start.gts)
  ends = sapply(out, end.gts)

  for(i in seq_along(out)) {
    ind = as.numeric(window(the_ts, start=starts[,i], end=ends[,i]))
    the_array[,, ind] = out[[i]]$x
    the_time[ind] = out[[i]]$time
    info_time[ind] = out[[i]]$info$time$time
  }
  class(the_time) = class(out[[1]]$time)

  x = out[[1]]
  x$x = the_array
  x$time = the_time
  x$info$time$time = info_time
  x$breaks$time = .getBreaks(x$info$time$time)
  x$info$dim$time = info_time
  x$info$ts = ts(seq_along(x$time), start=start(the_ts), frequency=frequency(the_ts))

  # class(x) = c("gts", "ts")
  class(x) = "gts"
  return(x)

}

#' @export
names.gts = function(x) {
  return(c(x$info$varid, x$info$long_name))
}

#' @export
'names<-.gts' = function(x, value) {

  if(any(is.na(value))) stop("NAs are not allowed in names for a gts object.")
  if(length(value)>2) stop("A maximum of two values (varid, long_name) must be provided.")

  x$info$varid = value[1]
  x$info$ovarid[length(x$info$ovarid)] = value[1]
  x$info$long_name = value[2]

  x

}

#' @export
units.gts = function(x) {
  return(tail(x$info$units, 1))
}

# generic implemented in base.
#' @export
'units<-.gts' = function(x, value) {

  if(length(value)!=1) stop("Only one value for units must be provided.")

  if(is.na(value)) value = ""

  x$info$units[length(x$info$units)] = value

  x

}

#' @export
str.gts = function(object, ...) {
  class(object) = "list"
  str(object, ...)
  return(invisible())
}

#' @export
quantile.gts = function(x, probs = seq(0, 1, 0.25), na.rm = TRUE, names = TRUE,
                        type = 7, digits = 7, ...) {

  MARGIN = seq_along(dim(x))
  MARGIN = MARGIN[-length(MARGIN)]

  out = apply(x$x, MARGIN=MARGIN, FUN=quantile, probs=probs, na.rm=na.rm, names=names,
        type=type, digits=digits, ...)

  lans = length(out)/prod(dim(x)[MARGIN])
  if(lans>1) out = aperm(out, perm = c(seq_along(dim(out))[-1], 1))

  x$x = out
  x$time = probs
  x$breaks$time = .getBreaks(x$time)
  x$info$time$time = probs
  x$info$dim$time = probs
  x$info$ts = NULL
  x$info$time_var = "prob"

  return(x)

}

#' @export
melt = function (data, ..., na.rm = FALSE, value.name = "value") {
  UseMethod("melt", data)
}

#' @export
melt.gts = function(data, ..., na.rm=FALSE, value.name=NULL) {

  time_var = data$info$time_var
  if(is.null(time_var)) time_var = "time"

  if(is.null(data$depth)) {
    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     time=rep(data$time, each=nrow(data$grid$df)),
                     value=as.numeric(data$x))

    names(out)[3:4] = c(time_var, names(data)[1])

  } else {

    out = data.frame(longitude=data$grid$df$lon, latitude=data$grid$df$lat,
                     depth = rep(data$depth, each=nrow(data$grid$df)))
    out = data.frame(longitude=out$longitude, latitude=out$latitude,
                     depth = out$depth,
                     time=rep(data$time, each=nrow(out)),
                     value=as.numeric(data$x))

    names(out)[4:5] = c(time_var, names(data)[1])

  }

  return(out)

}

#' @export
mean.gts = function(x, by=NULL, trim=0, na.rm=FALSE, ...) {
  if(is.null(by)) return(mean(as.double(x), trim=trim, na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=mean, type=by, trim=trim, na.rm=TRUE, ...)
  return(out)
}

#' @export
sum.gts = function(x, by=NULL, na.rm=FALSE, ...) {
  if(is.null(by)) return(sum(as.double(x), na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=sum, type=by, na.rm=TRUE, ...)
  return(out)
}

#' @export
as.double.gts = function(x, ...) {

  if(!is.null(x$grid$mask)) {
    return(as.double(x$x[!is.na(x$grid$mask)]))
  }
  return(as.double(x$x))
}

#' @export
sd = function(x, na.rm, ...) {
  UseMethod("sd")
}

#' @export
sd.default = function(x, na.rm=FALSE, ...) {
  stats::sd(x=x, na.rm=na.rm)
}

#' @export
sd.gts = function(x, by=NULL, na.rm=FALSE, ...) {
  if(is.null(by)) return(sd(as.double(x), na.rm=na.rm, ...))
  out = apply_gts(x=x, FUN=sd, type=by, na.rm=TRUE, ...)
  return(out)
}

#' @export
is_climatology = function (x, ...) {
  UseMethod("is_climatology")
}

#' @export
is_climatology.gts = function(x, ...) {
  return(x$info$climatology)
}

# S4 compatibility --------------------------------------------------------

setOldClass("gts")
setOldClass("grid")


# setMethod('Ops', signature(e1='gts', e2='ANY'), Ops.gts)
# setMethod('Arith', signature(e1='ANY', e2='gts'), Arith_gts)
# setMethod('Arith', signature(e1='gts', e2='ANY'), Ops.gts)


# Auxiliar ----------------------------------------------------------------

apply_gts = function(x, FUN, type=c("time", "space"), ...) {

  type = match.arg(type)

  if(type=="time") MARGIN = seq_along(dim(x))[-c(1,2)]
  if(type=="space") MARGIN = head(seq_along(dim(x)), -1)

  out = apply(x$x, MARGIN=MARGIN, FUN=FUN, ...)

  if(type=="time") {
    out = ts(out, start=start(x), frequency=frequency(x))
  }
  return(out)
}


