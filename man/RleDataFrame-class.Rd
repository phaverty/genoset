\name{RleDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{RleDataFrame-class}
\alias{as.matrix,RleDataFrame-method}
\alias{coerce,RleDataFrame,matrix-method}
\alias{colMeans,RleDataFrame-method}
\alias{colSums,RleDataFrame-method}
\alias{colMeans,DataFrame-method}
\alias{rowMeans,RleDataFrame-method}
\alias{rowSums,RleDataFrame-method}
\alias{show,RleDataFrame-method}
\alias{RleDataFrame}
\title{Class \code{"RleDataFrame"}}
\description{
The \code{RleDataFrame} class serves to hold a collection of Run Length
Encoded vectors (Rle objects) of the same length. For example, it could
be used to hold information along the genome for a number of samples,
such as sequencing coverage, DNA copy number, or GC content. This class
inherits from both DataFrame and SimpleRleList (one of the AtomicVector
types). This means that all of the usual subsetting and applying
functions will work. Also, the AtomicList functions, like \code{mean} and
\code{sum}, that automatically apply over the list elements will
work. The scalar mathematical \code{AtomicList} methods can make this
class behave much like a matrix (see Examples).

New objects can be created with the \code{RleDataFrame} constructor:
\code{RleDataFrame(..., row.names=NULL)}, where \code{...} can be a list
of Rle objects, or one or more individual Rle objects.
}

\section{Use in Biobase eSet objects}{
  The \code{genoset} class defines an \code{annotatedDataFrameFrom}
  method for \code{DataFrame}, which makes it possible to include
  \code{DataFrames} as \code{assayData} elements. The column names for
  \code{DataFrame} cannot be \code{NULL}, which makes it impossible to
  use them as assays in \code{SummarizedExperiment} at this time.
}

\section{Row and Column Summaries}{
  These objects will sometimes be in place of a matrix, as in the eSet
  example above. It is convenient to have some of the summarization
  methods for matrices. Each of these methods takes an RleDataFrame and
  returns a single Rle. The time required is similar to that required
  for a matrix. For an RleDataFrame x,
  \describe{
    \item{}{\code{rowSums}:} Sum across 'rows'.
    \item{}{\code{rowMeans}:} Means across 'rows'.
    \item{}{\code{colSums}:} Sum each Rle. This is just the \code{sum}
  method for SimpleRleList.
    \item{}{\code{colSums}:} Mean of each Rle. This is just the
  \code{mean} method for SimpleRleList.
  }
}
\section{Slots}{
  \describe{
    \item{\code{rownames}:}{Object of class \code{"characterORNULL"}
      Names to describe each row of the DataFrame. These may end up
      taking more space than your collection of Rle objects, so consider
    leaving this NULL.}
    \item{\code{nrows}:}{Object of class \code{"integer"} Number of rows. }
    \item{\code{elementType}:}{Object of class \code{"character"} Notes
      that elements of the internal list are Rle objects. }
    \item{\code{elementMetadata}:}{Object of class
      \code{"DataTableORNULL"} Metadata on the elements, see DataFrame. }
    \item{\code{metadata}:}{Object of class \code{"list"} Metadata on
      the whole object, see DataFrame. }
    \item{\code{listData}:}{Object of class \code{"list"} Base list
      containing the Rle objects. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SimpleRleList}"}, directly.
Class \code{"\linkS4class{DataFrame}"}, directly.
}
\section{Methods}{
  \describe{
    \item{as.matrix}{\code{signature(x = "RleDataFrame")}: Convert to
  matrix. }
    \item{coerce}{\code{signature(x = "RleDataFrame")}: Convert to other
  classes. }
    \item{colMeans}{\code{signature(x = "RleDataFrame")}: Mean of each column. }
    \item{colSums}{\code{signature(x = "RleDataFrame")}: Sum of each column. }
    \item{rowMeans}{\code{signature(x = "RleDataFrame")}: Mean of each 'row'. }
    \item{rowSums}{\code{signature(x = "RleDataFrame")}: Sum of each 'row'. }
    \item{show}{\code{signature(object = "RleDataFrame")}: Short and
      pretty description of an object of this type. }
  }
}
\author{
Peter M. Haverty, design suggestion from Michael Lawrence.
}
\seealso{
  \code{\linkS4class{DataFrame}}
  \code{\linkS4class{AtomicList}}
  \code{\linkS4class{Rle}}
  \code{\linkS4class{RleList}}
  \code{\link{rowMeans}}
  \code{\link{colMeans}}
  \code{\link{rowSums}}
  \code{\link{colSums}}
  \code{\link{view-summarization-methods}}
}
\examples{
showClass("RleDataFrame")

## Constructors
df = new("RleDataFrame", listData=list(A=Rle(c(NA, 2:3, NA, 5), rep(2,
5)), B=Rle(c(6:7, NA, 8:10),c(3,2,1,2,1,1))), nrows=10L)

df2 = RleDataFrame(list(A=Rle(c(NA, 2:3, NA, 5), rep(2, 5)),
B=Rle(c(6:7, NA, 8:10),c(3,2,1,2,1,1))))

df3 = RleDataFrame(A=Rle(c(NA, 2:3, NA, 5), rep(2, 5)), B=Rle(c(6:7,
NA, 8:10),c(3,2,1,2,1,1)))

## AtomicList Methods
runValue(df)
runLength(df)
ranges(df)
mean(df)
sum(df)
df + 5
log2(df) - 1

## Row and Column Summaries
rowSums(df)
colSums(df)
rowMeans(df)
colMeans(df)

## Coercion
as(df, "matrix")
as(df, "list")
as(df, "RleList")
as(df, "DataFrame")
as(df, "data.frame")
}
\keyword{classes}
\keyword{methods}
