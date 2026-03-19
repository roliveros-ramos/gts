#' Write gridded objects to a netCDF file
#'
#' These are S3 methods for the [nctools::write_ncdf()] generic.
#'
#' The methods write `gts`, `static`, and `grid` objects to a new netCDF file
#' using metadata already stored in the object. In particular, variable names,
#' dimensions, dimension units, and variable units are taken from the `info`
#' component of the object and passed to [nctools::write_ncdf()].
#'
#' For regular objects containing a single data variable, this typically writes
#' one variable to the output file. For objects whose metadata describe several
#' variables, such as irregular-grid objects with stored coordinate variables,
#' all variables listed in `x$info$var` are written.
#'
#' @param x A `gts`, `static`, or `grid` object.
#' @param filename Path to the output netCDF file.
#' @param varid Character string giving the name of the main output variable.
#'   For `gts` and `static` methods, if supplied, this replaces the last entry of
#'   `x$info$ovarid` before writing. In the current `grid` method
#'   implementation, variable names are taken from `x$info$ovarid`.
#' @param longname Optional long name or vector of long names passed to
#'   [nctools::write_ncdf()]. When several variables are written, this may be a
#'   length-one value to be recycled or one value per output variable.
#' @param prec Character string giving the storage precision passed to
#'   [ncdf4::ncvar_def()], for example `"short"`, `"integer"`, `"float"`, or
#'   `"double"`.
#' @param missval Missing value used in the netCDF variable definition.
#' @param compression Numeric compression level passed to
#'   [ncdf4::ncvar_def()]. Compression usually requires netCDF4 output.
#' @param chunksizes Optional chunk sizes for compressed variables.
#' @param verbose Logical; should `ncdf4` report progress while creating and
#'   filling the file?
#' @param dim.longname Character vector giving long names for the dimensions.
#'   If omitted, empty strings are used by [nctools::write_ncdf()].
#' @param unlim Character string naming the unlimited dimension. Use `FALSE`
#'   for no unlimited dimension.
#' @param global Named list of global attributes to add to the output file.
#'   A `history` attribute documenting the function call is added automatically
#'   by [nctools::write_ncdf()].
#' @param force_v4 Logical; should the output file be forced to netCDF4 format?
#' @param ... Additional arguments passed to [nctools::write_ncdf()].
#'
#' @details
#' These methods do not require the user to specify `dim`, `dim.units`, or
#' `units` explicitly. Instead, they derive those values from the object
#' metadata:
#' \describe{
#'   \item{`x$info$ovarid`}{Output variable names.}
#'   \item{`x$info$dim`}{Dimension coordinate values.}
#'   \item{`x$info$dim.units`}{Dimension units.}
#'   \item{`x$info$units`}{Variable units.}
#' }
#'
#' `write_ncdf.gts()` and `write_ncdf.static()` share the same implementation.
#' Both write the variables listed in `x$info$var`, using the stored metadata to
#' define the netCDF structure.
#'
#' `write_ncdf.grid()` behaves similarly for `grid` objects.
#'
#' @return These methods are called for their side effect of creating a netCDF
#'   file and return `invisible(NULL)`.
#'
#' @seealso [nctools::write_ncdf()], [gts-class], [static-class], [grid-class]
#'
#' @examples
#' \dontrun{
#' write_ncdf(x, "temperature.nc")
#' write_ncdf(bathy, "bathymetry.nc")
#' write_ncdf(grd, "grid.nc")
#' }
#' @name write_ncdf
NULL

#' @describeIn write_ncdf Write a `gts` object to a netCDF file.
#' @export
write_ncdf.gts = function(x, filename, varid, longname, prec="float",
                          missval=NA, compression=9, chunksizes=NA, verbose=FALSE,
                          dim.longname, unlim=FALSE, global=list(),
                          force_v4=FALSE, ...) {

  if(!missing(varid)) x$info$ovarid[length(x$info$ovarid)] = varid

  write_ncdf(x=x[x$info$var], filename=filename, varid = x$info$ovarid, longname=longname,
             dim=x$info$dim, dim.units=x$info$dim.units, units=x$info$units, prec=prec,
             missval=missval, compression=compression, chunksizes=chunksizes, verbose=verbose,
             dim.longname=dim.longname, unlim=unlim, global=global, force_v4=force_v4, ...)

  return(invisible(NULL))

}

#' @describeIn write_ncdf Write a `static` object to a netCDF file.
#' @export
write_ncdf.static = write_ncdf.gts

#' @describeIn write_ncdf Write a `grid` object to a netCDF file.
#' @export
write_ncdf.grid = function(x, filename, varid, longname, prec="float",
                           missval=NA, compression=9, chunksizes=NA, verbose=FALSE,
                           dim.longname, unlim=FALSE, global=list(),
                           force_v4=FALSE, ...) {

  if(!missing(varid)) x$info$ovarid[length(x$info$ovarid)] = varid

  write_ncdf(x=x[x$info$var], filename=filename, varid = x$info$ovarid, longname=longname,
             dim=x$info$dim, dim.units=x$info$dim.units, units=x$info$units, prec=prec,
             missval=missval, compression=compression, chunksizes=chunksizes, verbose=verbose,
             dim.longname=dim.longname, unlim=unlim, global=global, force_v4=force_v4, ...)

  return(invisible(NULL))

}
