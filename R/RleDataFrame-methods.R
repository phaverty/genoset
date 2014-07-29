##' @include RleDataFrame-class.R
NULL

##' @export colMeans
setGeneric("colMeans", function(x, na.rm=TRUE, dims=1L) standardGeneric("colMeans") )

setMethod("colMeans", "RleDataFrame",
          function(x, na.rm=TRUE) {
            mean(x, na.rm=na.rm)
          })

setMethod("colMeans", signature(x="DataFrame"),
          function(x,na.rm=TRUE,dims=1L) {
            .Deprecated("colMeans", msg="colMeans on a DataFrame is Deprecated. It is kind of odd given that the column type are arbitrary. Try RleDataFrame, or another class that inherits from AtomicList and DataFrame. But, if you find this DataFrame version useful, let me know.")
            return( vapply(x,mean,na.rm=na.rm, FUN.VALUE=numeric(1), USE.NAMES=TRUE) ) } )

##' @export colSums
setGeneric("colSums", function(x, na.rm=TRUE, dims=1L) standardGeneric("colSums") )
setMethod("colSums", "RleDataFrame",
          function(x, na.rm=TRUE) {
            sum(x, na.rm=na.rm)
          })
##' @export rowMeans
setGeneric("rowMeans", function(x, na.rm=FALSE, dims=1L) standardGeneric("rowMeans") )
setMethod("rowMeans", signature(x="ANY"), base::rowMeans)
setMethod("rowMeans", signature(x="RleDataFrame"),
          function(x, na.rm=FALSE, dims=1L) {
            if (na.rm==TRUE) {
              sums = x[[1L]]
              na.sums = is.na(runValue(sums))
              runValue(sums)[na.sums] = 0
              na.count = Rle(na.sums, runLength(sums))             
              for (i in 2L:ncol(x)) {
                current = x[[i]]
                is.na.current = is.na(runValue(current))
                na.current = Rle(is.na.current, runLength(current))
                runValue(current)[is.na.current] = 0
                sums = sums + current
                na.count = na.count + na.current
              }
              means = sums / (ncol(x) - na.count)
            } else {
              sums = x[[1L]]
              for (i in 2L:ncol(x)) {
                sums = sums + x[[i]]
              }
              means = sums / ncol(x)
            }
            return(means)
          })

##' @export rowSums
setGeneric("rowSums", function(x, na.rm=FALSE, dims=1L) standardGeneric("rowSums") )
setMethod("rowSums", signature(x="ANY"), base::rowSums)
setMethod("rowSums", signature(x="RleDataFrame"),
          function(x, na.rm=FALSE, dims=1L) {
          if (na.rm==TRUE) {
              sums = x[[1L]]
              runValue(sums)[is.na(runValue(sums))] = 0
              for (i in 2L:ncol(x)) {
                current = x[[i]]
                runValue(current)[is.na(runValue(current))] = 0
                sums = sums + current
              }
            } else {
              sums = x[[1L]]
              for (i in 2L:ncol(x)) {
                sums = sums + x[[i]]
              }
            }
            return(sums)
        })

# Allow eSet constructor to make featureNames from a DataFrame as if it were a matrix
##' @export annotatedDataFrameFrom
setMethod("annotatedDataFrameFrom",
          signature(object="DataFrame"),
          Biobase:::annotatedDataFrameFromMatrix)
