\name{RleDataFrame-views}
\alias{RleDataFrame-views}

\alias{rangeSums}
\alias{rangeMeans}
\alias{rangeMins}
\alias{rangeMaxs}
\alias{rangeWhichMins}
\alias{rangeWhichMaxs}
\alias{rangeSums,RleDataFrame-method}
\alias{rangeMeans,RleDataFrame-method}
\alias{rangeMeans,numeric-method}
\alias{rangeMeans,matrix-method}
\alias{rangeMins,RleDataFrame-method}
\alias{rangeMaxs,RleDataFrame-method}
\alias{rangeWhichMins,RleDataFrame-method}
\alias{rangeWhichMaxs,RleDataFrame-method}
\alias{rangeColMeans}
\title{Calculate summary statistics on views of an RleDataFrame}

\description{
  These methods mirror the viewMeans type functions from IRanges for
  SimpleRleList. They differ in that they work on an RleDataFrame and an
  IRanges directly and also have a simplify argument. This works out to be
  faster (compute-wise) and also convenient.

  Still, an RleDataFrame inherits from SimpleRleList, so all of the
  views functions will work.
}

\details{
  The "range" name prefixes here serve to differentiate these functions
  from the "view" functions. This may change. I will be asking the
  IRanges team to add "..." and "simplify" to the "view" methods so that
  I can just make additional methods for \code{RleDataFrame}.
}

\usage{
rangeSums(x, bounds, na.rm=FALSE, simplify=TRUE)
rangeMeans(x, bounds, na.rm=FALSE, simplify=TRUE, ...)
rangeMins(x, bounds, na.rm=FALSE, simplify=TRUE)
rangeMaxs(x, bounds, na.rm=FALSE, simplify=TRUE)
rangeWhichMins(x, bounds, na.rm=FALSE, simplify=TRUE)
rangeWhichMaxs(x, bounds, na.rm=FALSE, simplify=TRUE)
}

\arguments{
  \item{x}{
    RleDataFrame
  }
  \item{bounds}{
    Matrix with two columns or IRanges representing ranges of rows of
  \code{x} to process. If \code{bounds} is a matrix, an IRanges is
  constructed assuming the first two columns represent the start and end
  of the ranges. The names for the IRanges is taken from the rownames of
  the matrix. Such a matrix can constructed with
  \code{boundingIndicesByChr} and is the preferred input.
  }
  \item{na.rm}{
    Scalar logical. Ignore NAs in calculations?
  }
  \item{simplify}{
    Scalar logical. Simplify result? If TRUE, the return value will be a
    vector or matrix. For a single view, a vector will be
    returned. Otherwise a matrix with one row per view and one column per
    column of x will be returned. If FALSE, the return value will be a
    list of length \code{ncol(x)} of vectors of length \code{nrow(bounds)}.
  }
  \item{...}{
    Additional arguments for other methods.
  }
}

\value{
   With \code{simplify == TRUE}, a vector for single view or a matrix
otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(bounds)}.
}

\examples{
df = RleDataFrame(list(a=Rle(1:5, rep(2, 5))), b=Rle(1:5, rep(2, 5)),
row.names=LETTERS[1:10])
mat = matrix(c(1,4,3,5),ncol=2,dimnames=list(c("Gene1","Gene2"),c("start","end")))
bounds = IRanges(start=c(1, 4), end=c(3, 5), names=c("Gene1","Gene2"))

rangeMeans(df,bounds,simplify=FALSE)
rangeMeans(df,bounds,simplify=TRUE)
rangeMeans(df,mat,simplify=TRUE)

rangeMeans(df,bounds)
rangeSums(df,bounds)
rangeMins(df,bounds)
rangeMaxs(df,bounds)
rangeWhichMins(df,bounds)
rangeWhichMaxs(df,bounds)

# RleDataFrame isa SimpleRleList, so all the IRanges view* methods work too:
viewMeans( Views( df, bounds) )
}

\seealso{
  \linkS4class{RleDataFrame} \link{boundingIndicesByChr}
}

\keyword{methods}
