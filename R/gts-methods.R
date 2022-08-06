
# S3 methods --------------------------------------------------------------

#' @exportS3Method is.na gts
is.na.gts = function(x) {
  x$x = is.na(x$x)
  return(x)
}

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


# S4 compatibility --------------------------------------------------------


setOldClass("gts")
setOldClass("grid")

# setMethod('Ops', signature(e1='gts', e2='ANY'), Ops.gts)
# setMethod('Arith', signature(e1='ANY', e2='gts'), Arith_gts)
# setMethod('Arith', signature(e1='gts', e2='ANY'), Ops.gts)

