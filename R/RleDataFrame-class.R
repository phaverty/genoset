##### Additional methods for DataFrame class to allow use of DataTable of Rle vectors as assayData element
##' @include genoset-class.R
NULL

##' @exportClass RleDataFrame
setClass("RleDataFrame",
         representation(
           rownames = "character_OR_NULL",
           nrows = "integer"
           ),
         prototype(rownames = NULL,
                   nrows = 0L,
                   listData = structure(list(),  names = character())),
         contains = c("SimpleRleList", "DataFrame")
)

##' @export RleDataFrame
RleDataFrame <- function(..., row.names=NULL, check.names=TRUE) {
  x = DataFrame(..., row.names=row.names, check.names=check.names)
  return( new("RleDataFrame", listData=x@listData, rownames=rownames(x), nrows=nrow(x)) )
}

##' @exportMethod show
setMethod("show", "RleDataFrame",
          function(object) {
            message(sprintf("RleDataFrame with %i rows and %i columns\n", nrow(object), ncol(object)))
            if (is.null(rownames(object))) {
              message("rownames: NULL\n")
            } else if (nrow(object) > 5){
              message(sprintf("rownames: %s, ...\n", paste(head(rownames(object)), collapse=", ")))
            } else {
              message(sprintf("rownames: %s\n", paste(rownames(object), collapse=", ")))
            }
            show(object@listData)
          })


### Coercion

##' @export
setAs("RleDataFrame", "matrix",
      function(from){
        mat = vapply(from, as.vector, vector(storage.mode(runValue(from[[1]])), nrow(from)), USE.NAMES=TRUE)
        dimnames(mat) = dimnames(from)
        return(mat)
      })

##' @exportMethod as.matrix
setMethod("as.matrix", "RleDataFrame", function(x) {as(x, "matrix")})

