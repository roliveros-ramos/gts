#' Time-series methods for `gts` objects
#'
#' These methods expose the regular time-series metadata stored in a `gts`
#' object through familiar time-series generics.
#'
#' `time.gts()`, `frequency.gts()`, `deltat.gts()`, `start.gts()`, `end.gts()`,
#' and `cycle.gts()` all delegate to the regular `ts` object stored in
#' `x$info$ts`.
#'
#' `window.gts()` subsets a `gts` object using [stats::window()] on that regular
#' time-series representation and updates the data array, time metadata, and
#' time breaks accordingly.
#'
#' @param x A `gts` object.
#' @param ... Additional arguments passed to the underlying method. For
#'   `time.gts()`, `frequency.gts()`, `deltat.gts()`, `start.gts()`, `end.gts()`,
#'   and `cycle.gts()`, these are passed to the corresponding generic. For
#'   `window.gts()`, they are passed to [stats::window()].
#'
#' @details
#' These methods operate on the regular time-series representation stored in
#' `x$info$ts`, not directly on `x$time`. In particular, `time.gts()` returns
#' the regular numeric time scale of the underlying `ts` object, whereas
#' `x$time` typically stores date-like values.
#'
#' `window.gts()` currently uses the standard `ts` indexing rules of
#' [stats::window()]. It therefore expects `start`, `end`, and related
#' arguments on the regular time scale of `x$info$ts`. Date or character
#' expansion is not currently implemented.
#'
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`time.gts()`}{A numeric vector giving the regular time scale of
#'   `x$info$ts`.}
#'   \item{`frequency.gts()`}{A numeric scalar giving the observation frequency
#'   of `x$info$ts`.}
#'   \item{`deltat.gts()`}{A numeric scalar giving the time step of
#'   `x$info$ts`.}
#'   \item{`start.gts()`, `end.gts()`}{Numeric vectors in the usual `ts`
#'   representation of the start or end of the series.}
#'   \item{`cycle.gts()`}{An integer vector giving the cycle position of each
#'   observation.}
#'   \item{`window.gts()`}{A subsetted `gts` object with updated data, time
#'   vector, time metadata, and time breaks.}
#' }
#'
#' @seealso [gts-class], [is_climatology()], [gridded-summary],
#'   [stats::window()], [stats::ts()]
#'
#' @examples
#' \dontrun{
#' frequency(x)
#' cycle(x)
#' x2 <- window(x, start = c(2000, 1), end = c(2005, 12))
#' clim <- climatology(x, FUN = "mean")
#' }
#' @name gts-time
NULL

#' @describeIn gts-time Return the regular time scale of a `gts` object.
#' @export
time.gts = function(x, ...) {
  time(x$info$ts)
}

#' @describeIn gts-time Return the observation frequency of a `gts` object.
#' @export
frequency.gts = function(x, ...) {
  frequency(x$info$ts)
}

#' @describeIn gts-time Return the time step of a `gts` object.
#' @export
deltat.gts = function(x, ...) {
  deltat(x$info$ts)
}

#' @describeIn gts-time Return the start of a `gts` time series.
#' @export
start.gts = function(x, ...) {
  start(x$info$ts)
}

#' @describeIn gts-time Return the end of a `gts` time series.
#' @export
end.gts = function(x, ...) {
  end(x$info$ts)
}

#' @describeIn gts-time Return the cycle positions of a `gts` time series.
#' @export
cycle.gts = function(x, ...) {
  cycle(x$info$ts)
}

#' @describeIn gts-time Subset a `gts` object using `ts`-style time indexing.
#' @export
window.gts = function(x, ...) {
  nts = window(x$info$ts, ...)
  ind = as.numeric(nts)
  hasdepth = !is.null(x$depth)
  x$x = if(hasdepth) x$x[, , , ind] else x$x[, , ind]
  x$time = x$time[ind]
  x$info$time$time = x$info$time$time[ind]
  x$breaks$time = .getBreaks(x$info$time$time)
  x$info$dim$time = x$info$dim$time[ind]
  x$info$ts = ts(seq_along(x$time), start=start(nts), frequency=frequency(nts))
  return(x)
}


#' Compute climatologies from gridded time-series objects
#'
#' `climatology()` is an S3 generic for computing climatological cycles from
#' time-indexed objects.
#'
#' `climatology.gts()` aggregates a `gts` object by `cycle(x)` and returns a new
#' `gts` object representing one climatological cycle.
#'
#' @param x A `gts` object.
#' @param FUN Aggregation function used to compute the climatology. This can be
#'   a function or the name of a function. The function is applied over all
#'   observations belonging to the same cycle.
#' @param ... Additional arguments passed to `FUN`.
#'
#' @details
#' `climatology.gts()` groups observations by `cycle(x)` using the regular
#' `ts` representation stored in `x$info$ts`.
#'
#' For three-dimensional objects, aggregation is performed for each horizontal
#' grid cell over the last dimension. For four-dimensional objects, aggregation
#' is performed for each horizontal grid cell and depth level over the last
#' dimension.
#'
#' The returned object preserves the `gts` structure but redefines the time axis
#' so that it represents cycle positions rather than the original sequence of
#' observations. Specifically:
#' \itemize{
#'   \item `x$x` is replaced by the aggregated climatology array;
#'   \item `x$time` is replaced by `seq_len(frequency(x))`;
#'   \item `x$breaks$time` is recomputed from those cycle values;
#'   \item `x$info$time$time` is replaced by the mean of the original time
#'   labels within each cycle;
#'   \item `x$info$dim$time` is replaced by the cycle values;
#'   \item `x$info$ts` is reset to a cycle-based [stats::ts()] object starting
#'   at `0`;
#'   \item `x$info$climatology` is set to `TRUE`.
#' }
#'
#' The internal helper currently supports only arrays of dimension 3 or 4.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`climatology()`}{The result of the selected method.}
#'   \item{`climatology.gts()`}{A `gts` object representing a climatological
#'   cycle.}
#' }
#'
#' @seealso [gts-time], [is_climatology()], [cycle.gts()], [stats::ts()]
#'
#' @examples
#' \dontrun{
#' clim <- climatology(x, FUN = "mean")
#' }
#' @name climatology-gts
NULL

#' @describeIn climatology-gts Compute a climatology from an object.
#' @export
climatology = function(x, ...) {
  UseMethod("climatology")
}

#' @describeIn climatology-gts Aggregate a `gts` object by cycle to create a
#'   climatology.
#' @export
climatology.gts = function(x, FUN="mean", ...) {

  x = drop(x)

  index = as.numeric(cycle(x))
  values = seq_len(frequency(x))

  clim = .climatology(x$x, index=index, values=values, FUN=FUN, ...)

  x$x = clim
  x$time = values
  x$breaks$time = .getBreaks(x$time)
  x$info$time$time = tapply(x$info$time$time, INDEX = cycle(x), FUN="mean")
  x$info$dim$time = values
  x$info$ts = ts(values, start=0, frequency=frequency(x))
  x$info$climatology = TRUE

  return(x)
}

# Auxiliar functions ------------------------------------------------------

get_time = function(x) {
  .this_time = function(x) {
    ini = ymd_hms(sprintf("%s-01-01 00:00:00", year(x)))
    end = ymd_hms(sprintf("%s-01-01 00:00:00", year(x)+1))
    out = year(x) + as.numeric(difftime(x, ini, units="secs"))/as.numeric(difftime(end, ini, units="secs"))
    return(out)
  }
  out = sapply(x, FUN=.this_time)
  return(out)
}

get_date = function(x) {
  year = floor(x)
  rest = x - year
  out = ymd(sprintf("%s-01-01", year))
  out = out + seconds(rest*(365 + leap_year(year))*24*60*60)
  return(out)
}

get_freq = function(x, tolerance=0, freq=NULL) {

  if(length(x)==1) {
    if(!is.null(freq)) return(freq)
    warning("Not enough values to compute frequency.")
    return(NA)
  }
  x = as.POSIXct(x)

  tol = max(sqrt(.Machine$double.eps), tolerance)

  dt = diff(as.numeric(x)) # difference in seconds
  dm = days_in_month(x)

  yfreq_d = 24*60*60/dt
  yfreq_m = head(dm, -1)*24*60*60/dt
  yfreq_y = (365+head(leap_year(x), -1))*24*60*60/dt

  freq_scale = c(day=count_integers(yfreq_d, tolerance=tol),
                 month=count_integers(yfreq_m, tolerance=tol),
                 year=count_integers(yfreq_y, tolerance=tol))

  all_freq = c(day=365*popular_integer(yfreq_d, tolerance=tol),
               month=12*popular_integer(yfreq_m, tolerance=tol),
               year=popular_integer(yfreq_y, tolerance=tol))

  freq = round(all_freq[which.max(freq_scale)],0)

  if(is.na(freq)) return(get_freq(x, tolerance=tolerance+0.0025))

  return(freq)
}

guess_origin = function(time, nyear_max=100) {

  useyear  = if(mean(diff(time))>nyear_max) FALSE else TRUE
  usemonth = if(mean(diff(time))>12*nyear_max) FALSE else TRUE
  useday   = if(mean(diff(time))>366*nyear_max) FALSE else TRUE

  guess = expand.grid(time_unit=c("second", "minute", "hour",
                              if(useday) "day" else NULL,
                              if(usemonth) "month" else NULL,
                              if(useyear) "year" else NULL),
                      origin=c("1970-01-01", "1950-01-01", "1900-01-01"), stringsAsFactors = FALSE)
  guess$freq = NA
  guess$start = NA
  guess$end = NA

  for(i in seq_len(nrow(guess))) {
    nt = time2date(time, units=guess$time_unit[i], origin=guess$origin[i])
    gf = get_freq(nt)
    if(gf==0) next
    guess$freq[i] = gf
    guess$start[i] = as.character(nt[1])
    guess$end[i] = as.character(tail(nt, 1))
  }
  guess = guess[complete.cases(guess), ]
  return(guess)
}

.climatology = function(object, index, values, FUN="mean", ...) {

  FUN = match.fun(FUN)
  dims = dim(object)
  if(!(length(dims) %in% c(3,4)))
    stop("climatology computation only supported for arrays of dimension 3 or 4.")
  ind  = values
  if(is.numeric(ind)) ind = sort(ind)
  out = array(NA, dim=c(dims[-length(dims)], length(ind)))
  for(i in seq_along(ind)) {
    if(length(dims)==3)
      out[,,i] = apply(object[,, which(index == ind[i])], 1:2, FUN, na.rm=TRUE, ...)
    if(length(dims)==4)
      out[,,,i] = apply(object[,,, which(index == ind[i])], 1:3, FUN, na.rm=TRUE, ...)
  }
  dimnames(out)[[length(dims)]] = ind
  return(out)

}

# Internal functions ------------------------------------------------------

count_integers = function(x, tolerance=sqrt(.Machine$double.eps)) {
  isint = sapply(x %% 1, FUN=all.equal, current=0, tolerance=tolerance)
  isint = as.logical(isint)
  isint[is.na(isint)] = FALSE
  return(sum(isint))
}

popular_integer = function(x, tolerance=sqrt(.Machine$double.eps)) {
  isint = sapply(x %% 1, FUN=all.equal, current=0, tolerance=tolerance)
  isint = as.logical(isint)
  isint[is.na(isint)] = FALSE
  ints = table(x[isint])
  out = as.numeric(names(ints)[which.max(ints)])
  if(length(out)==0) return(NA)
  return(out)
}


time2date = function(x, units="seconds", origin="1970-01-01", calendar=NULL) {
  if(inherits(origin, "character"))
    origin = parse_date_time(origin, orders = c("Ymd", "YmdHMS", "dmY", "dmYHMS"))
  if(is.na(origin)) stop("Origin must be a rightful date.")
  nt = try(do.call(period, args=setNames(list(x), nm=units)), silent = TRUE)
  if(inherits(nt, "try-error")) {
    if(units=="month") {
      xnt = do.call(period, args=setNames(list(floor(x)), nm=units))
      xtime = origin + xnt
      dim = if(calendar==365) days_in_month(xtime) else 30
      xx = (x%%1)*dim*24*60*60
      new_time = xtime + seconds(xx)
    }
    if(units=="year") {
      xnt = do.call(period, args=setNames(list(floor(x)), nm=units))
      xtime = origin + xnt
      if(is.null(calendar)) calendar = 365 + leap_year(xtime)
      xx = (x%%1)*calendar*24*60*60
      new_time = xtime + seconds(xx)
    }
    if(units=="week") {
      xnt = do.call(period, args=setNames(list(floor(x)), nm=units))
      xtime = origin + xnt
      xx = (x%%1)*7*24*60*60
      new_time = xtime + seconds(xx)
    }
    if(units=="fortnight") {
      xnt = do.call(period, args=setNames(list(floor(14*x)), nm="day"))
      xtime = origin + xnt
      xx = (x%%1)*24*60*60
      new_time = xtime + seconds(xx)
    }
    if(units=="day") {
      xnt = do.call(period, args=setNames(list(floor(x)), nm=units))
      xtime = origin + xnt
      xx = (x%%1)*24*60*60
      new_time = xtime + seconds(xx)
    }
  } else {
    new_time = origin + nt
  }
  return(new_time)
}

