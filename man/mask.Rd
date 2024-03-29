% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mask.R
\name{mask}
\alias{mask}
\alias{mask.grid}
\alias{mask.gts}
\alias{mask.array}
\alias{mask<-}
\title{Extract, create or assign ocean/land masks.}
\usage{
mask(x, ...)

\method{mask}{grid}(x, n = 2, thr = 0.8, hires = FALSE, ocean = TRUE, ...)

\method{mask}{gts}(x, n = 2, thr = 0.8, hires = FALSE, ocean = TRUE, ...)

\method{mask}{array}(x, ...)

mask(x) <- value
}
\arguments{
\item{x}{An object containing a mask.}

\item{...}{Additional arguments, currently not used.}

\item{n}{Number of internal points used to calculate the mask.}

\item{thr}{Threshold to assign a point to the ocean.}

\item{hires}{Boolean, use high-resolution coastline?}

\item{ocean}{Boolean, are valid points in the ocean? Default is TRUE. FALSE for land.}

\item{object}{The object to add or replace the mask.}

\item{mask}{A mask.}
}
\value{
The mask.
}
\description{
Extract, create or assign ocean/land masks.
}
