

#' Seasonal Decomposition of Gridded Time Series
#'
#' @param x An array with a time dimension.
#' @param frequency The frequency in time-steps per year, default to 12.
#' @param s.window Only "periodic" is supported now.
#' @param ... Additional arguments controling the methods.
#'
#' @return A list, with the trend, seasonal and residual components.
#' @export
#'
stl_array = function(x, frequency=12, s.window="periodic", ...) {
  # optimize, do loop?

  .stl_array = function(x, frequency, s.window, ...) {
    xna = rep(NA, length=length(x))
    outna = list(data=data.frame(raw=xna, seasonal=xna, trend=xna,
                                 remainder=xna, weights=xna))
    if(all(is.na(x))) return(outna)
    out = try(stlplus(ts(x, frequency = frequency), s.window = s.window, ...))
    if(inherits(out, "try-error")) return(outna)
    return(out)
  }

  xx = apply(x, 1:2, FUN=.stl_array, frequency=frequency,
             s.window=s.window, ...)
  .xget = function(x, var) x[[1]]$data[[var]]

  z1 = aperm(apply(xx, 1:2, FUN=.xget, var="trend"),
             perm = c(2,3,1))
  z2 = aperm(apply(xx, 1:2, FUN=.xget, var="seasonal"),
             perm = c(2,3,1))
  z3 = aperm(apply(xx, 1:2, FUN=.xget, var="remainder"),
             perm = c(2,3,1))

  out = list(trend=z1, seasonal=z2, remainder=z3)

  return(out)

}


#' Seasonal Decomposition of Gridded Time Series
#'
#' @param x An array with a time dimension.
#' @param frequency The frequency in time-steps per year, default to 12.
#' @param s.window Only "periodic" is supported now.
#' @param ... Additional arguments controling the methods.
#'
#' @return A list, with the trend, seasonal and residual components.
#' @export
#'
stl_array2 = function(x, frequency=12, s.window="periodic", ...) {
  # optimize, do loop?

  DateStamp("Starting...")

  .stl_array = function(x, frequency, s.window, ...) {
    xna = rep(NA, length=length(x))
    outna = list(data=data.frame(raw=xna, seasonal=xna, trend=xna,
                                 remainder=xna, weights=xna))
    if(all(is.na(x))) return(outna)
    out = try(stlplus(ts(x, frequency = frequency), s.window = s.window, ...),
              silent=TRUE)
    if(inherits(out, "try-error")) return(outna)
    return(out)
  }

  z1 = array(NA_real_, dim=dim(x))
  z2 = array(NA_real_, dim=dim(x))
  z3 = array(NA_real_, dim=dim(x))

  for(i in seq_len(ncol(x))) {

    pb = txtProgressBar(style=3)
    setTxtProgressBar(pb, (i-1)/ncol(x))

    for(j in seq_len(nrow(x))) {


      tmp = .stl_array(x[j, i, ], frequency=frequency, s.window=s.window, ...)
      z1[j, i, ] = tmp$data$trend
      z2[j, i, ] = tmp$data$seasonal
      z3[j, i, ] = tmp$data$remainder

    }
  }

  out = list(trend=z1, seasonal=z2, remainder=z3)

  pb = txtProgressBar(style=3)
  setTxtProgressBar(pb, 1)

  DateStamp("Finishing at")

  return(out)

}
