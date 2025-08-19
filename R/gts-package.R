
#' Read a netCDF variable into a Gridded Time-Series.
#'
#' @param filename The filename of the netCDF file.
#' @param varid String with  the name of the variable in the file to read the
#' data from. If left unspecified, the function will determine if there is only
#' one variable in the file and, if so, read from that. If left unspecified and
#' there are multiple variables in the file, an error is generated.
#' @param control A list with additional options, particularly the time units and origin when not specified in the file.
#' @param ... Additional arguments.
#'
#' @return A gts object
#' @export
#'
read_gts = function(filename, varid=NULL, control=list(), ...) {

  nc = nc_open(filename = filename)
  on.exit(nc_close(nc))

  output = gts(x=nc, varid=varid, control=control, ...)

  return(output)

}


