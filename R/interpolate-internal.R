
# Internal functions: interpolate -----------------------------------------

.check_input = function(x, y, z) {

  x = drop(x)
  y = drop(y)

  if(length(x)==1) stop("Degenerate dimension x.")
  if(length(y)==1) stop("Degenerate dimension y.")

  if(length(dim(z))==4) z = z[, , 1, 1, drop=FALSE]
  if(length(dim(z))==3) z = z[, , 1, drop=FALSE]
  z = drop(z)

  isMatrix = is.matrix(x) & is.matrix(y)

  # Either both inputs are vectors or matrices.
  if(!isMatrix & is.matrix(x)) stop("Input y must be a matrix too.")
  if(!isMatrix & is.matrix(y)) stop("Input x must be a matrix too.")

  # check if matrices can be reduced to vectors!

  isMonot  = monot(x, side=2) & monot(y, side=1)
  sameDim  = identical(dim(x), dim(y)) & identical(dim(x), dim(z))
  sameLen  = (length(x)==length(y)) & (length(x)==length(z))
  consDim  = identical(c(length(x), length(y)), dim(z))

  is_gridR = isTRUE(consDim & isMonot)
  is_gridI = isTRUE(sameDim & isMonot & isMatrix)
  is_point = isTRUE(sameLen & !isMatrix)

  if(!any(is_gridR, is_gridI, is_point)) {
    if(!is.matrix(x) & any(is.na(x))) stop("Input vector 'x' cannot contain NAs.")
    if(!is.matrix(y) & any(is.na(y))) stop("Input vector 'y' cannot contain NAs.")
    # pending validation
    stop("Inputs x, y and z are not consistent.")
  }

  hasNA = any(is.na(x)) | any(is.na(y)) | any(is.na(z))

  if(isTRUE(hasNA)) {

    if(!isTRUE(is_point)) {
      if(!is.matrix(x)) x = matrix(x, ncol=ncol(z), nrow=nrow(z))
      if(!is.matrix(y)) y = matrix(y, ncol=ncol(z), nrow=nrow(z), byrow=TRUE)
    }

    dat = data.frame(x=as.numeric(x), y=as.numeric(y), z=as.numeric(z))
    dat = dat[complete.cases(dat), ]

    out = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point),
               hasNA=TRUE, data=dat)

    return(out)

  }

  out = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point),
             hasNA=FALSE, data=NULL)
  return(out)

}

.check_output = function(x, y) {

  x = drop(x)
  y = drop(y)

  isMatrix = is.matrix(x) & is.matrix(y)

  if(!isMatrix & is.matrix(x)) stop("Input y must be a matrix too.")
  if(!isMatrix & is.matrix(y)) stop("Input x must be a matrix too.")

  isMonot  = monot(x, side=2) & monot(y, side=1)
  sameDim  = identical(dim(x), dim(y))
  sameLen  = (length(x)==length(y))

  is_gridR = isTRUE(isMonot & !isMatrix)
  is_gridI = isTRUE(isMonot &  isMatrix)
  is_point = isTRUE(!isMonot & sameLen & !isMatrix)

  if(!any(is_gridR, is_gridI, is_point)) {
    # pending validation
    stop("Arguments 'xout' and 'yout' are not consistent.")
  }

  output = list(case=c(is_gridR=is_gridR, is_gridI=is_gridI, is_point=is_point))

  return(output)

}

monot = function(x, side) {
  .monot  = function(x) all(diff(x) > 0)
  .monotM = function(x, side) all(apply(x, side, .monot))
  if(is.matrix(x)) return(.monotM(x, side))
  return(.monot(x))
}

get_nearest_index = function(x, y, xout, yout, mask=NULL, extrap=FALSE) {
  indx = cut(as.numeric(xout), breaks=.getBreaks(x), labels=FALSE)
  indy = cut(as.numeric(yout), breaks=.getBreaks(y), labels=FALSE)
  return(cbind(indx, indy))
}

get_bilinear_weights = function(x, y, xout, yout, mask=NULL, extrap=FALSE) {

  nx = length(x)
  ny = length(y)
  n_out = length(xout)
  if(length(yout)!=n_out) stop("The lengths of 'xout' and 'yout' do not match.")

  if(!is.null(mask)) {
    mask[mask==0] = NA
  }

  dx = diff(x)
  dy = diff(y)

  # Preallocate vectors for sparseMatrix()
  row_idx = integer(n_out * 4)
  col_idx = integer(n_out * 4)
  weights = numeric(n_out * 4)
  invalid = logical(n_out)

  counter = 1

  for (i in seq_len(n_out)) {
    # Find grid cell (lower-left corner)
    ix = findInterval(xout[i], x, rightmost.closed=TRUE)
    iy = findInterval(yout[i], y, rightmost.closed=TRUE)

    # Skip if out of bounds
    if (ix <= 0 || ix >= nx || iy <= 0 || iy >= ny) {
      invalid[i] = TRUE
      next
    }

    # Get relative position inside the cell
    x0 = x[ix]
    y0 = y[iy]
    tx = (xout[i] - x0)/dx[ix]
    ty = (yout[i] - y0)/dy[iy]

    # Bilinear weights
    w00 = (1 - tx) * (1 - ty)
    w10 = tx * (1 - ty)
    w01 = (1 - tx) * ty
    w11 = tx * ty

    # Linear indices in vectorized input (column-major)
    idx00 = (ix - 1) + (iy - 1) * nx + 1
    idx10 = ix + (iy - 1) * nx + 1
    idx01 = (ix - 1) + iy * nx + 1
    idx11 = ix + iy * nx + 1

    ww  = c(w00, w10, w01, w11)
    idx = c(idx00, idx10, idx01, idx11)

    if(!is.null(mask)) {
      imask = mask[idx] # mask points, NA is land.
      no_na = sum(!is.na(imask))
      if(!isTRUE(extrap)) {
        if(no_na==0) {
          invalid[i] = TRUE
          next
        }
        if(no_na %in% c(1,2)) {
          iw = which(is.na(imask))
          nw = sum(ww[iw])
          if(nw > 0) {
            invalid[i] = TRUE
            next
          }
        }
        if(no_na==3) {
          iw = which(is.na(imask)) # only one
          ww[-iw] = ww[-iw] + ww[iw]/3
          ww[iw] = 0
        }
      } else {
        if (no_na < 4) {
          iw_na = which(is.na(imask))     # NAs in mask
          iw_ok = which(!is.na(imask))    # alid values
          w_na  = sum(ww[iw_na])          # Total "orphaned" weight
          if(length(iw_ok) == 0) {
            invalid[i] = TRUE
            next
          }
          ww[iw_ok] = ww[iw_ok] + w_na/length(iw_ok)  # Redistribute evenly
          ww[iw_na] = 0                   # Remove NAs
        }
      }
    }

    # Store in sparse matrix vectors
    row_idx[counter:(counter + 3)] = rep(i, 4)
    col_idx[counter:(counter + 3)] = idx
    weights[counter:(counter + 3)] = ww
    counter = counter + 4
  }

  xind = which(weights > 0)
  row_idx = row_idx[xind]
  col_idx = col_idx[xind]
  weights = weights[xind]

  # Build sparse matrix
  W = sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = weights,
    dims = c(n_out, nx * ny)
  )

  invalid = which(invalid)
  attr(W, "invalid") = invalid

  return(W)
}
