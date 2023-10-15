#' @rdname gts
#' @export
as.gts = function(x, ...) {
  UseMethod("as.gts")
}

#' @export
as.gts.array = function(x, grid, ...) {

  if(inherits(grid, "gts")) grid = grid$grid
  if(!inherits(grid, "grid")) stop("grid must be of class 'grid'.")
  out = list(x=x,
             longitude=grid$longitude,
             latitude=grid$latitude,
             depth=NULL,
             time=NULL,
             breaks=NULL,
             grid=grid,
             )

  list(varid=, long_name=, time=, depth=, dim=, var=,
  dim.units=, units=, ovarid=, global=, ts=)


}

#' @export
is.gts = function(x) {
  inherits(x, "gts")
}

