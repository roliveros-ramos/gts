


#' @export
write_ncdf.gts = function(x, filename, longname, prec="float",
                          missval=NA, compression=9, chunksizes=NA, verbose=FALSE,
                          dim.longname, unlim=FALSE, global=list(),
                          force_v4=FALSE, ...) {

  write_ncdf(x=x[x$info$var], filename=filename, varid = x$info$ovarid, longname=longname,
             dim=x$info$dim, dim.units=x$info$dim.units, units=x$info$units, prec=prec,
             missval=missval, compression=compression, chunksizes=chunksizes, verbose=verbose,
             dim.longname=dim.longname, unlim=unlim, global=global, force_v4=force_v4, ...)

  return(invisible(NULL))

}

