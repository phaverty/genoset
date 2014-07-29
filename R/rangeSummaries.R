### Methods to do view operations on an RleDataFrame
###  viewMeans(Views(rledf,iranges)) works, but gives list.  Consistency good, but I really wanna simplify via vapply here ...

##' @include RleDataFrame-class.R
NULL

##' Calculate min/max/sum/mean/whichmin/whichmax over each view on each column of an RleDataFrame.
##'
##' Loop over the Rle objects in an RleDataFrame, calculate the appropriate statistic for each view. If simplify == FALSE,
##' this function returns a vector for each Rle. If simplify == TRUE, it returns a vector for the case of a single view, otherwise,
##' a matrix. Rownames for the matrix are taken from the names of the argument \code{bounds}.
##' @param x RleDataFrame
##' @param bounds IRanges or matrix, views on every Rle. If \code{bounds} is a matrix, it is converted to an IRanges using the first
##' two columns as the starts and stops. Names for the IRanges are taken from the rownames of the matrix. Such a matrix can be
##' constructed with \code{boundingIndicesByChr}.
##' @param na.rm scalar logical, ignore NAs in calculations?
##' @param simplify scalar logical, simplify result? For a single view, a vector, otherwise a matrix with one row per view.
##' @param RLEFUN function, internal rle view summary function like .rle_view_means
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @return With \code{simplify == TRUE}, a vector for single view or a matrix
##' otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(bounds)}.
##' @keywords internal
##' @rdname do_rledf_views
##' @seealso RleDataFrame boundingIndicesByChr
##' @family views
.do_rledf_views <- function(x, bounds, na.rm=FALSE, simplify=TRUE, RLEFUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
  # Make an IRanges from ranges matrix if necessary
  if (is.matrix(bounds)) {
    bounds = IRanges(start=bounds[, 1], end=bounds[, 2], names=rownames(bounds))
  }
  # Trim IRanges once if necessary
  start(bounds)[ start(bounds) < 1L ] = 1L
  end(bounds)[ end(bounds) > nrow(x) ] = nrow(x)
  # Hoist the Views dispatch
  myviewfun = getMethod("Views", "Rle", where="IRanges")
  # Calculate the view stats
  if (simplify == TRUE) {
    FUN.TYPE = match.arg(FUN.TYPE)
    nviews = length(bounds)
    val = vapply(x,
           FUN=function(rle) {
             RLEFUN(myviewfun(rle, bounds), na.rm=na.rm)
           }, USE.NAMES=TRUE,
           FUN.VALUE=structure( vector( FUN.TYPE, nviews ), names=names(bounds) ) )
  } else {
    val = lapply(x,
      function(rle) {
        RLEFUN(myviewfun(rle, bounds), na.rm=na.rm)
      })
    }
  return(val)
}

##' Calculate min/max/sum/mean/whichmin/whichmax over ranges on each column of an RleDataFrame.
##'
##' Loop over the Rle objects in an RleDataFrame, calculate the appropriate statistic for each view. If simplify == FALSE,
##' this function returns a vector for each Rle. If simplify == TRUE, it returns a vector for the case of a single range, otherwise,
##' a matrix. Rownames for the matrix are taken from the names of the argument \code{bounds}.
##' @param x RleDataFrame
##' @param bounds IRanges or matrix, views on every Rle. If \code{bounds} is a matrix, the first
##' two columns are used as as the starts and stops. Names for the ranges are taken from rownames of the matrix. Such a matrix can be
##' constructed with \code{boundingIndicesByChr}.
##' @param na.rm scalar logical, ignore NAs in calculations?
##' @param simplify scalar logical, simplify result? For a single view, a vector, otherwise a matrix with one row per view.
##' @param RLEFUN function, internal rle view summary function like .rle_range_means
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @return With \code{simplify == TRUE}, a vector for single view or a matrix
##' otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(bounds)}.
##' @keywords internal
##' @rdname do_rledf_range_summary
##' @seealso RleDataFrame boundingIndicesByChr
##' @family views
.do_rledf_range_summary <- function(x, bounds, na.rm=FALSE, simplify=TRUE, RLEFUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
  # Make an IRanges from ranges matrix if necessary
  if (is.matrix(bounds)) {
    start=bounds[, 1]
    end=bounds[, 2]
    names=rownames(bounds)
  } else if (is(bounds, "IRanges")) {
    start = start(bounds)
    end = end(bounds)
    names = names(bounds)
  } else {
    stop("x must be a two-column matrix or an IRanges.")
  }
  # Trim IRanges once if necessary
  start[ start < 1L ] = 1L
  end[ end > nrow(x) ] = nrow(x)

  # Calculate the view stats
  if (simplify == TRUE) {
    FUN.TYPE = match.arg(FUN.TYPE)
    val = vapply(x,
           FUN=function(rle) {
             RLEFUN(start, end, runValue(rle), runLength(rle), na.rm=na.rm)
           }, USE.NAMES=TRUE,
           FUN.VALUE=structure( vector( FUN.TYPE, length(start) ), names=names) )
  } else {
    val = lapply(x,
      function(rle) {
        structure( RLEFUN(start, end, runValue(rle), runLength(rle), na.rm=na.rm), names=names)
      })
    }
  return(val)
}

##' @export rangeSums
setGeneric("rangeSums", function(x, bounds, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeSums") })
setMethod("rangeSums", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_sums, FUN.TYPE="numeric")
          })

##' @export rangeMins
setGeneric("rangeMins", function(x, bounds, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeMins") })
setMethod("rangeMins", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_mins, FUN.TYPE="numeric")
          })

##' @export rangeMaxs
setGeneric("rangeMaxs", function(x, bounds, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeMaxs") })
setMethod("rangeMaxs", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_maxs, FUN.TYPE="numeric")
          })

##' @export rangeWhichMins
setGeneric("rangeWhichMins", function(x, bounds, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeWhichMins") })
setMethod("rangeWhichMins", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_mins, FUN.TYPE="integer")
          })

##' @export rangeWhichMaxs
setGeneric("rangeWhichMaxs", function(x, bounds, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeWhichMaxs") })
setMethod("rangeWhichMaxs", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_maxs, FUN.TYPE="integer")
          })


##' @export rangeMeans
setGeneric("rangeMeans", function(x, bounds, na.rm=FALSE, simplify=TRUE, ...) { standardGeneric("rangeMeans") })
setMethod("rangeMeans", signature=signature(x="RleDataFrame"),
          function(x, bounds, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_range_summary(x, bounds, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_range_means, FUN.TYPE="numeric")
          })

setMethod("rangeMeans", signature=signature(x="numeric"),
          function(x, bounds, na.rm=FALSE) {
              if (!is.double(x)) {
                  storage.mode(x) = "double"
              }
              ans = .Call("rangeMeans_numeric", bounds, x, na.rm)
              return(ans)
          })

setMethod("rangeMeans", signature=signature(x="matrix"), # S4 does not see a class relationship between numeric and matrix. Hulk smash S4.
          function(x, bounds, na.rm=FALSE) {
              if (!is.double(x)) {
                  storage.mode(x) = "double"
              }
              ans = .Call("rangeMeans_numeric", bounds, x, na.rm)
              return(ans)
          })

##' @export rangeColMeans
rangeColMeans <- function(x, all.indices) {
  .Deprecated("rangeMeans", msg="rangeColMeans has changed to rangeMeans. Please note that the order of arguments is different too.")
  rangeMeans(x, all.indices, na.rm=TRUE)
}

### Internal methods to get directly to summary functions, using Views, but skipping trim
.rle_view_sums <- function(x, na.rm) { .Call("RleViews_viewSums", x, na.rm, PACKAGE = "IRanges") }
.rle_view_means <- function(x, na.rm) { .Call("RleViews_viewMeans", x, na.rm, PACKAGE = "IRanges") }
.rle_view_mins <- function(x, na.rm) { .Call("RleViews_viewMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_maxs <- function(x, na.rm) { .Call("RleViews_viewMaxs", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_mins <- function(x, na.rm) { .Call("RleViews_viewWhichMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_maxs <- function(x, na.rm) { .Call("RleViews_viewWhichMaxs", x, na.rm, PACKAGE = "IRanges") }

### Internal methods to get directly to summary functions, skipping trim and Views
.rle_range_means <- function(start, end, values, lengths, na.rm) {
  .Call("rangeMeans_rle", as.integer(start), as.integer(end), as.numeric(values), lengths, na.rm=na.rm, PACKAGE = "genoset")
}

##' Count Rle positions >= min
##'
##' For Rle coverage vector, count number of positions where value >= min, think callable bases.
##' @param rle integer Rle, no NAs
##' @param bounds IRanges or matrix, positions in Rle to consider. If \code{bounds} is a matrix, the first
##' two columns are used as start and end.
##' @param min scalar integer, count Rle positions >= this value.
##' @return integer vector of length nrow(bounds)
##' @export
numCallable <- function(rle, bounds, min) {
    #### disjoin !
    if (is.matrix(bounds)) {
        if (storage.mode(bounds) != "integer") { storage.mode(bounds) = "integer" }
        start=bounds[, 1]
        end=bounds[, 2]
    } else if (is(bounds, "IRanges")) {
        start = start(bounds)
        end = end(bounds)
    } else {
        stop("x must be a two-column matrix or an IRanges.")
    }
    .Call("numCallable_rle", start, end, runValue(rle), runLength(rle), as.integer(min), PACKAGE = "genoset")
}
