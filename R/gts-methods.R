#' Quantiles of a gridded time-series object
#'
#' `quantile.gts()` computes empirical quantiles over the last dimension of a
#' `gts` object, which in standard `gts` objects corresponds to time.
#'
#' Quantiles are computed independently for each spatial cell, and for each depth
#' level when present. The returned object keeps the `gts` structure, but the
#' time axis is replaced by the requested probability values.
#'
#' @param x A `gts` object.
#' @param probs Numeric vector of probabilities passed to [stats::quantile()].
#' @param na.rm Logical; whether missing values should be removed before
#'   computing quantiles.
#' @param names Logical; passed to [stats::quantile()].
#' @param type Integer quantile algorithm passed to [stats::quantile()].
#' @param digits Integer used by [stats::quantile()] when `names = TRUE`.
#' @param ... Additional arguments passed to [stats::quantile()].
#'
#' @details
#' The function applies [stats::quantile()] over all dimensions except the last
#' one. For a three-dimensional `gts` object, quantiles are computed for each
#' horizontal grid cell over time. For a four-dimensional `gts` object, they are
#' computed for each horizontal grid cell and depth level over time.
#'
#' The returned object preserves the class and most metadata structure of the
#' input, but its time representation is redefined so that the last dimension
#' indexes probabilities rather than dates or periods. Specifically:
#' \itemize{
#'   \item `x$time` is replaced by `probs`;
#'   \item `x$breaks$time` is recomputed from `probs`;
#'   \item `x$info$time$time` and `x$info$dim$time` are replaced by `probs`;
#'   \item `x$info$ts` is set to `NULL`;
#'   \item `x$info$time_var` is set to `"prob"`.
#' }
#'
#' This means that downstream functions which use `time_var` for labelling, such
#' as `melt.gts()` or `plot.gts()`, will treat the last dimension as
#' probabilities rather than time.
#'
#' @return A modified `gts` object whose last dimension indexes probabilities
#'   rather than time.
#'
#' @seealso [gridded-summary], [stats::quantile()], [melt.gts()], [plot.gts()]
#'
#' @examples
#' \dontrun{
#' qx <- quantile(x, probs = c(0.1, 0.5, 0.9))
#' }
#' @name quantile-gts
NULL

#' @describeIn quantile-gts Compute empirical quantiles for each spatial cell of
#'   a `gts` object.
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


# S4 compatibility --------------------------------------------------------

setOldClass("gts")
setOldClass("grid")
setOldClass("static")
# setMethod('Ops', signature(e1='gts', e2='ANY'), Ops.gts)
# setMethod('Arith', signature(e1='ANY', e2='gts'), Arith_gts)
# setMethod('Arith', signature(e1='gts', e2='ANY'), Ops.gts)


# Auxiliar ----------------------------------------------------------------

apply_gts = function(x, FUN, type=c("time", "space", "latitude", "longitude"), ...) {

  type = match.arg(type)

  x = drop(x)

  if(type=="time") MARGIN = seq_along(dim(x))[-c(1,2)]
  if(type=="space") MARGIN = head(seq_along(dim(x)), -1)
  if(type=="latitude") MARGIN = seq_along(dim(x))[-1]
  if(type=="longitude") MARGIN = seq_along(dim(x))[-2]

  out = apply(x$x, MARGIN=MARGIN, FUN=FUN, ...)

  if(type=="time") {
    out = ts(out, start=start(x), frequency=frequency(x), names=names(x)[1])
  }
  if(type=="space") {
    x$x = out
    x$time = NULL
    x$breaks$time = NULL
    x$info$time$time = NULL
    x$info$dim$time = NULL
    x$info$ts = NULL
    x$info$climatology = NULL
    class(x) = "static"
    return(x)
  }
  return(out)
}


