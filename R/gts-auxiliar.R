#' @export
slog = function(x) sign(x)*log(abs(x))

#' Reduce the object of a GAM object by eliminating caller environment references
#'
#' @param object A GAM or GLM object
#' @param ... Additional arguments
#'
#' @returns
#' @export
compact_gam = function(object, ...) {
  if(!inherits(object, "glm")) stop("object must at least inherit class 'glm'.")
  oclass = class(object)
  if(!is.null(object$dinfo)) {
    environment(object$dinfo$gp$pf) = baseenv()
    environment(object$dinfo$gp$fake.formula) = environment()
    environment(object$dinfo$gp$pred.formula) = baseenv()
  }
  object$Sl <- object$qrx <- object$R <- object$F <- NULL
  object$Ve <- object$Vc <- object$G <- object$residuals <- NULL
  object$fitted.values <- object$linear.predictors <- NULL
  object$na.action = NULL
  object$weights = NULL
  if(!is.null(object$model) & length(object$model!=0)) {
    object$model = object$model[FALSE, ]
    attr(object$model, "na.action") = NULL
  }
  if(!is.null(object$data) & !identical(object$data, NA))  object$data  = object$data[FALSE, ]
  object$prior.weights = NULL
  object$wt = NULL
  object$offset = NULL
  object$y = NULL
  object$call = NULL

  # # Remove environments from formula and terms
  if (!is.null(object$formula)) environment(object$formula) = baseenv()
  if (!is.null(object$pred.formula)) environment(object$pred.formula) = baseenv()
  if (!is.null(object$terms)) environment(object$terms) = baseenv()
  if (!is.null(object$pterms)) {
    environment(object$pterms) = environment()
    assign(x="offset", value=stats::offset, envir = environment(object$pterms))
  }
  # Remove environment from smooth terms (each smooth can hold a reference!)
  if (!is.null(object$smooth)) {
    for (i in seq_along(object$smooth)) {
      if (!is.null(object$smooth[[i]]$term)) environment(object$smooth[[i]]$term) = baseenv()
    }
  }
  class(object) = oclass
  return(object)
}
