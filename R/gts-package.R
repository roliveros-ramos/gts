#' Gridded time-series objects
#'
#' A `gts` object stores a gridded environmental time series together with its
#' spatial coordinates, time axis, grid geometry, and metadata.
#'
#' In this package, the first two dimensions represent horizontal space, an
#' optional third dimension represents depth, and the last dimension represents
#' time.
#'
#' `gts` objects are designed to support a workflow in which gridded
#' environmental data are imported, transformed, aligned on a common grid,
#' summarised, reshaped, and combined with other gridded objects.
#'
#' @section Main components:
#' A `gts` object is a list-like object with the following main components:
#' \describe{
#'   \item{`x`}{Numeric array containing the data values. The first two
#'   dimensions are horizontal space, an optional third dimension is depth, and
#'   the last dimension is time.}
#'   \item{`longitude`, `latitude`}{Horizontal coordinates. For regular grids
#'   these are usually numeric vectors. For irregular grids they may be stored as
#'   matrices.}
#'   \item{`depth`}{Optional numeric vector of depth coordinates. Present when
#'   the data include an explicit depth dimension.}
#'   \item{`time`}{Time vector associated with the last dimension of `x`.}
#'   \item{`breaks`}{List of cell-boundary coordinates for each dimension.}
#'   \item{`grid`}{Associated [grid-class] object describing the horizontal
#'   geometry.}
#'   \item{`info`}{Metadata list containing variable identifiers, display names,
#'   original dimension values, dimension units, variable units, time metadata,
#'   optional global attributes, an optional `ts` representation of the time
#'   axis, and a logical `climatology` flag.}
#' }
#'
#' @section Time representation:
#' `gts` objects store both a time vector in `time` and a richer description of
#' the time axis in `info$time`. When possible, a regular [stats::ts()]
#' representation is also stored in `info$ts`.
#'
#' @section Climatologies:
#' A climatology is a `gts` object for which `info$climatology` is `TRUE`. Such
#' an object represents one full cycle at the object frequency rather than a
#' continuous historical time series.
#'
#' @section Regular and irregular grids:
#' The horizontal geometry of a `gts` object is defined by its associated
#' [grid-class] object. For regular grids, longitude and latitude are usually
#' stored as vectors. For irregular grids, the geometry is represented by
#' coordinate matrices.
#'
#' @seealso [gts()], [read_gts()], [grid-class], [static-class]
#' @name gts-class
#' @keywords classes
NULL


#' Spatial grid objects
#'
#' A `grid` object stores the horizontal geometry and associated spatial metadata
#' used by gridded objects in the package.
#'
#' `grid` objects are used by both [gts-class] and [static-class] objects to
#' define the spatial arrangement of cells, masks, areas, and plotting
#' coordinates.
#'
#' @section Main components:
#' A `grid` object is a list-like object with the following main components:
#' \describe{
#'   \item{`longitude`, `latitude`}{Horizontal coordinate definitions. For
#'   regular grids these are typically numeric vectors. For irregular grids they
#'   may be matrices.}
#'   \item{`LON`, `LAT`}{Matrix representations of the horizontal coordinates on
#'   the main grid.}
#'   \item{`rho`}{List containing `LON` and `LAT` on the main grid.}
#'   \item{`psi`}{Optional staggered or cell-corner coordinates, used notably by
#'   plotting and interpolation helpers.}
#'   \item{`area`}{Optional matrix of cell areas.}
#'   \item{`mask`}{Optional spatial mask indicating which cells are considered
#'   valid or active.}
#'   \item{`prob`}{Optional spatial probability or proportion surface used to
#'   derive a binary mask.}
#'   \item{`n`}{Optional integer controlling the number of colours or levels
#'   used when visualising the probability surface.}
#'   \item{`hires`}{Logical flag indicating whether high-resolution coastlines
#'   should be used by plotting helpers.}
#'   \item{`df`}{Flattened coordinate table, typically with longitude and
#'   latitude columns, used by reshaping methods such as `melt()`.}
#'   \item{`info`}{Metadata list describing the grid variables and dimensions.}
#' }
#'
#' @section Regular and irregular grids:
#' A grid is regular when longitude and latitude can be represented as monotonic
#' coordinate vectors. A grid is irregular when the horizontal geometry must be
#' represented by coordinate matrices.
#'
#' Many functions in the package work with both types, but some operations are
#' currently implemented only for regular grids.
#'
#' @section Masks and probabilities:
#' When present, `mask` is the main binary spatial mask used by downstream
#' methods. `prob` stores an underlying continuous surface, such as the
#' proportion of ocean coverage in each cell, from which a mask may be derived.
#'
#' @seealso [make_grid()], [read_grid()], [update_grid()], [gts-class],
#'   [static-class]
#' @name grid-class
#' @keywords classes
NULL


#' Static gridded objects
#'
#' A `static` object stores a spatially structured field that does not vary over
#' time.
#'
#' Typical examples include bathymetry, distance to coast, cell area, habitat
#' suitability masks, or other spatial covariates that are constant through
#' time.
#'
#' `static` objects are designed to share a grid-aware interface with
#' [gts-class] objects so they can be inspected, plotted, reshaped, summarised,
#' and combined arithmetically with time-varying gridded data.
#'
#' @section Main components:
#' A `static` object is a list-like object with the following main components:
#' \describe{
#'   \item{`x`}{Numeric matrix or array containing the spatial field. The first
#'   two dimensions represent horizontal space. An optional third dimension may
#'   represent depth.}
#'   \item{`grid`}{Associated [grid-class] object describing the horizontal
#'   geometry.}
#'   \item{`longitude`, `latitude`}{Horizontal coordinates associated with the
#'   field.}
#'   \item{`depth`}{Optional depth coordinate when the field has a third
#'   dimension.}
#'   \item{`breaks`}{List of cell-boundary coordinates for each dimension.}
#'   \item{`info`}{Metadata list containing variable identifiers, display names,
#'   original dimension values, dimension units, variable units, and optional
#'   global attributes.}
#' }
#'
#' @section Relation to `gts` objects:
#' `static` objects do not have a time axis, but otherwise aim to behave
#' similarly to `gts` objects where this is meaningful. In particular, they can
#' interact with `gts` objects in mathematical operations, allowing time-varying
#' fields to be combined with time-invariant spatial covariates.
#'
#' @section Regular and irregular grids:
#' As with `gts` objects, the horizontal geometry of a `static` object is
#' defined by its associated [grid-class] object and may correspond to either a
#' regular or an irregular grid.
#'
#' @seealso [read_static()], [gts-class], [grid-class]
#' @name static-class
#' @keywords classes
NULL



