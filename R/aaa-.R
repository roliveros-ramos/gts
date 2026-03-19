#' @importFrom ncdf4 nc_open nc_close ncvar_get ncatt_get
#' @importFrom nctools ncvar_dim write_ncdf
#' @export write_ncdf
#' @importFrom fields poly.image image.plot interp.surface.grid tim.colors
#' imageplot.info imageplot.setup
#' @importFrom lubridate yday<- year month ymd leap_year parse_date_time period
#' ceiling_date floor_date ymd_hms days_in_month parse_date_time
#' days_in_month seconds
#' @importFrom mgcv gam
#' @importFrom akima interp interpp bilinear bicubic
#' @importFrom sf st_sfc st_transform st_point st_as_sf st_crs 'st_crs<-' st_join
#' @importFrom colorful divergencePalette
#' @importFrom stats gaussian integrate quantile sd cycle frequency deltat time
#' start end window complete.cases median predict setNames splinefun ts
#' @importFrom maps map
#' @importFrom Matrix sparseMatrix
#' @importFrom methods Math2
#' @importFrom grDevices rainbow
#' @importFrom graphics axis box image mtext par title
#' @importFrom utils head setTxtProgressBar str tail txtProgressBar
#' @importFrom stlplus stlplus
NULL
