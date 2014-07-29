\name{GenoSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GenoSet-class}
\alias{GenoSetOrGenomicRanges-class}


\title{Class \code{"GenoSet"}}
\description{
\code{GenoSet} extends \code{eSet} by adding genome location information in the form of the locData slot.
GenoSet uses this location information to allow quick subsetting and summarization by a set of genome locations (GRanges). 
GenoSet implements and extends the GRanges API for access to the underlying location information.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GenoSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
However, as per BioConductor standard practice the object creation function \code{GenoSet} is recommended.
}
\section{Slots}{
  \describe{
    \item{\code{locData}:}{Object of class \code{"GenomicRanges"} Locations of features on the genome }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} From eSet }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} From eSet }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} From eSet }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} From eSet }
    \item{\code{annotation}:}{Object of class \code{"character"} From eSet }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} From eSet }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} From eSet }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "GenoSet", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "GenoSet", i = "character", j = "ANY", drop = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "GenoSet", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{chr}{\code{signature(object = "GenoSet")}: ... }
    \item{chrNames}{\code{signature(object = "GenoSet")}: ... }
    \item{elementLengths}{\code{signature(x = "GenoSet")}: ... }
    \item{featureNames}{\code{signature(object = "GenoSet")}: ... }
    \item{featureNames<-}{\code{signature(object = "GenoSet")}: ... }
    \item{sampleNames}{\code{signature(object = "GenoSet")}: ... }
    \item{dim}{\code{signature(object = "GenoSet")}: ... }
    \item{genoPlot}{\code{signature(x = "GenoSet", y = "ANY")}: ... }
    \item{locData}{\code{signature(object = "GenoSet")}: ... }
    \item{names}{\code{signature(x = "GenoSet")}: ... }
    \item{ranges}{\code{signature(x = "GenoSet")}: ... }
    \item{chrInfo}{\code{signature(x = "GenoSet")}: ... }
    \item{chrIndices}{\code{signature(x = "GenoSet")}: ... }
    \item{show}{\code{signature(object = "GenoSet")}: ... }
    \item{toGenomeOrder}{\code{signature(ds = "GenoSet")}: ... }
    \item{isGenomeOrder}{\code{signature(ds = "GenoSet")}: ... }
  }
}

\seealso{
\code{\link{GenoSet}}
}
\examples{
showClass("GenoSet")
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
gs = GenoSet(
   locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
   cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
   pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
   annotation="SNP6"
)

}
\keyword{classes}
