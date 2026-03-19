#' Mathematical operations for gridded objects
#'
#' These methods apply mathematical transformations and operators to gridded
#' objects while preserving their object structure where possible.
#'
#' The page covers two related classes:
#' \itemize{
#'   \item `gts` objects, which store gridded time-series data;
#'   \item `static` objects, which store gridded spatial fields without a time
#'   dimension.
#' }
#'
#' Unary mathematical functions from the `Math` group, such as `abs()`,
#' `sqrt()`, `exp()`, `sin()`, and `cos()`, are applied element-wise to the
#' data array stored in `x$x`.
#'
#' Two-argument mathematical functions from the `Math2` group, such as
#' `round()` and `signif()`, are likewise applied element-wise.
#'
#' Arithmetic, comparison, and logical operators are handled through the `Ops`
#' group. These methods support operations between gridded objects and scalars,
#' arrays, matrices, and, in some cases, between different gridded classes.
#'
#' @param x A `gts` or `static` object.
#' @param e1,e2 Operands for unary or binary operators.
#' @param base Base of the logarithm for `log.gts()` and `log.static()`.
#' @param na.rm Logical; passed to the summary function.
#' @param ... Additional arguments passed to the underlying base generic.
#' @param y,mx,my,cl,reverse Formal arguments required by
#'   `chooseOpsMethod.gts()` and `chooseOpsMethod.static()`.
#'
#' @details
#' For `Math.gts()` and `Math.static()`, the operation is applied element-wise
#' to `x$x` and the modified object is returned.
#'
#' Cumulative `Math` operations are not supported for either class. In
#' particular, `cumsum()`, `cumprod()`, `cummin()`, and `cummax()` are rejected.
#'
#' `log.gts()` and `log.static()` are defined explicitly so that the `base`
#' argument is handled directly.
#'
#' `Math2.gts()` and `Math2.static()` apply the selected two-argument function
#' element-wise to the underlying data array.
#'
#' For binary operators:
#'
#' \strong{`gts` objects}
#' \itemize{
#'   \item Unary operators are applied directly to `e1$x`.
#'   \item If `e2` is not a `gts` object, the operation is applied between
#'   `e1$x` and `e2`.
#'   \item If `e2` inherits from `static`, its `x` component is used.
#'   \item If `e2` is a matrix, only the first two spatial dimensions are
#'   checked for compatibility in the current implementation.
#'   \item If both operands are `gts` objects, they must have the same number of
#'   dimensions, compatible grids, and identical frequencies.
#'   \item If one operand is a climatology and the other is not, the
#'   climatology is expanded to the cycles of the non-climatology object before
#'   the operation is applied.
#' }
#'
#' The current implementation does not fully verify that two `gts` objects have
#' identical time axes; it checks frequency, but not complete time identity.
#'
#' \strong{`static` objects}
#' \itemize{
#'   \item Unary operators are applied directly to `e1$x`.
#'   \item If `e2` is not a `static` object, the operation is applied between
#'   `e1$x` and `e2`.
#'   \item If both operands are `static` objects, they must have the same number
#'   of dimensions and compatible grids.
#' }
#'
#' \strong{Mixed `gts`/`static` operations}
#'
#' Mixed operations are intentionally handled by the `gts` operator method.
#' `chooseOpsMethod.static()` returns `FALSE` when the other operand is a `gts`
#' object, so dispatch falls through to `Ops.gts()`. As a result, operations
#' such as `gts + static` and `static + gts` return a `gts`-style result and
#' follow the `gts` compatibility rules.
#'
#' `Summary.gts()` and `Summary.static()` apply summary group generics such as
#' `min()`, `max()`, `range()`, `prod()`, `any()`, and `all()` to the full
#' underlying data array and do not preserve the original class.
#'
#' @return
#' Depending on the method:
#' \describe{
#'   \item{`Math.gts()`, `Math.static()`, `Math2.gts()`, `Math2.static()`,
#'   `log.gts()`, `log.static()`, `Ops.gts()`, `Ops.static()`}{A modified object
#'   of the same class as the leading operand, except for mixed `gts`/`static`
#'   operations, which are handled by `Ops.gts()` and therefore return a `gts`
#'   object.}
#'   \item{`Summary.gts()`, `Summary.static()`}{The result returned by the
#'   corresponding summary generic applied to the full data array.}
#'   \item{`chooseOpsMethod.gts()`}{`TRUE`, indicating that the `gts` method
#'   should be used for operator dispatch.}
#'   \item{`chooseOpsMethod.static()`}{`FALSE` when the other operand is a `gts`
#'   object, otherwise `TRUE`.}
#' }
#'
#' @seealso [gridded-summary], [gts()], [read_gts()], [read_static()]
#'
#' @examples
#' \dontrun{
#' y <- log(x)
#' z <- abs(x)
#' w <- x + 1
#' s <- bathy + temp
#' }
#' @name gridded-math
NULL

#' @describeIn gridded-math Apply unary mathematical functions element-wise to a
#'   `gts` object.
#' @export
Math.gts = function(x, ...) {
  if(.Generic %in% c("cummax", "cummin", "cumprod", "cumsum"))
    stop(gettextf("'%s' not defined for \"gts\" objects",
                  .Generic), domain = NA)
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}

#' @describeIn gridded-math Apply unary mathematical functions element-wise to a
#'   `static` object.
#' @export
Math.static = function(x, ...) {
  if(.Generic %in% c("cummax", "cummin", "cumprod", "cumsum"))
    stop(gettextf("'%s' not defined for \"static\" objects",
                  .Generic), domain = NA)
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}

#' @describeIn gridded-math Apply `log()` to a `gts` object.
#' @export
log.gts = function(x, base = exp(1)) {
  x$x = log(x$x, base=base)
  return(x)
}

#' @describeIn gridded-math Apply `log()` to a `static` object.
#' @export
log.static = log.gts

#' @describeIn gridded-math Apply two-argument mathematical functions
#'   element-wise to a `gts` object.
#' @export
Math2.gts = function(x, ...) {
  x$x = do.call(.Generic, list(x=x$x, ...))
  return(x)
}

#' @describeIn gridded-math Apply two-argument mathematical functions
#'   element-wise to a `static` object.
#' @export
Math2.static = Math2.gts

#' @describeIn gridded-math Internal dispatch helper for operators involving
#'   `gts` objects.
#' @export
chooseOpsMethod.gts = function(x, y, mx, my, cl, reverse) TRUE

#' @describeIn gridded-math Internal dispatch helper for operators involving
#'   `static` objects.
#' @export
chooseOpsMethod.static = function(x, y, mx, my, cl, reverse) {
  if(inherits(y, "gts")) return(FALSE)
  return(TRUE)
}

#' @describeIn gridded-math Apply arithmetic, comparison, or logical operators
#'   to `gts` objects.
#' @export
Ops.gts = function(e1, e2) {

  if(missing(e2)) {
    e1$x = get(.Generic, mode="function")(e1$x)
    return(e1)
  }

  if(!inherits(e1, "gts") & inherits(e2, "gts")) {
    return(get(.Generic, mode="function")(e2, e1))
  }

  if(inherits(e1, "gts") & !inherits(e2, "gts")) {
    if(inherits(e2, "static")) e2 = e2$x
    if(is.matrix(e2)) {
      if(!identical(dim(e1$x)[1:2], dim(e2))) {
        stop("Dimensions of both objects are not compatible.")
      }
      e2 = as.numeric(e2)
    }
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }

  ok0 = length(dim(e1)) == length(dim(e2))
  if(!ok0) stop(sprintf("Dimensions of both objects are not compatible (%d!=%d).",
                        length(dim(e1)), length(dim(e2))))

  ok1 = .compare_grids(e1$grid, e2$grid, strict=FALSE)
  if(!ok1) stop("Grids are not compatible.")

  ok2 = identical(frequency(e1), frequency(e2))
  if(!ok2) stop("Time frequencies are not compatible.")

  if(is_climatology(e1) & !is_climatology(e2)) {
    return(get(.Generic, mode="function")(e2, e1))
  }

  if(is_climatology(e2)) {
    if(length(dim(e2)) == 3) e2$x = e2$x[, , as.numeric(cycle(e1)), drop=FALSE]
    if(length(dim(e2)) == 4) e2$x = e2$x[, , , as.numeric(cycle(e1)), drop=FALSE]
  }

  ok3 = identical(dim(e1$x), dim(e2$x))
  if(!ok3) stop(gettextf("Input 'gts' objects dimensions not compatible for '%s'.",
                         .Generic), domain = NA)

  e1$x = get(.Generic, mode="function")(e1$x, e2$x)
  e1$info$climatology = is_climatology(e1) & is_climatology(e2)

  return(e1)
}

#' @describeIn gridded-math Apply arithmetic, comparison, or logical operators
#'   to `static` objects.
#' @export
Ops.static = function(e1, e2) {

  if(missing(e2)) {
    e1$x = get(.Generic, mode="function")(e1$x)
    return(e1)
  }

  if(inherits(e1, "static") & !inherits(e2, "static")) {
    if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }

  ok0 = length(dim(e1)) == length(dim(e2))
  if(!ok0) stop(sprintf("Dimensions of both objects are not compatible (%d!=%d).",
                        length(dim(e1)), length(dim(e2))))

  ok1 = .compare_grids(e1$grid, e2$grid, strict=FALSE)
  if(!ok1) stop("Grids are not compatible.")

  e1$x = get(.Generic, mode="function")(e1$x, e2$x)
  return(e1)
}

#' @describeIn gridded-math Apply summary group generics to the full data array
#'   of a `gts` object.
#' @export
Summary.gts = function(x, ..., na.rm=TRUE) {
  get(.Generic, mode="function")(x$x, na.rm=na.rm)
}

#' @describeIn gridded-math Apply summary group generics to the full data array
#'   of a `static` object.
#' @export
Summary.static = Summary.gts
