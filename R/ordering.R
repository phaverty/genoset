#############################################################################
#####  Functions related to the order of chromosomes and range features #####
#############################################################################


##' Order chromosome names in proper genome order
##'
##' Chromosomes make the most sense orded by number, then by letter.
##' 
##' @param chr.names character, vector of unique chromosome names
##' @return character vector of chromosome names in proper order
##' @export chrOrder
##' @family "genome ordering"
##' @examples
##'    chrOrder(c("chr5","chrX","chr3","chr7","chrY"))  #  c("chr3","chr5","chr7","chrX","chrY")
chrOrder <- function(chr.names) {
  simple.names = gsub("^chr","",chr.names)
  name.is.numeric = grepl("^[0-9]+$",simple.names,perl=T)
  numeric.names = chr.names[name.is.numeric][ order(as.numeric(simple.names[name.is.numeric])) ]
  non.numeric.names = chr.names[! name.is.numeric][ order(chr.names[ !name.is.numeric]) ]
  all.names = c(numeric.names,non.numeric.names)
  return(all.names)
}

##' Check if a GRanges orGenoSet is in genome order
##'
##' Checks that rows in each chr are ordered by start.  If strict=TRUE, then chromosomes
##' must be in order specified by chrOrder. isGenomeOrder for GRanges differs from order
##' in that it orders by chromsome and start position only,
##' rather than chromsome, strand, start, and width.
##' 
##' @param ds GenoSet or GRanges
##' @param strict logical, should space/chromosome order be identical to that from chrOrder?
##' @return logical
##' @export isGenomeOrder
##' @family "genome ordering"
##' @examples
##'   data(genoset)
##'   isGenomeOrder( locData(genoset.ds) )
##' @rdname isGenomeOrder-methods
setGeneric("isGenomeOrder", function(ds, strict=TRUE) standardGeneric("isGenomeOrder"))

##' @rdname isGenomeOrder-methods
setMethod("isGenomeOrder",signature=signature(ds="GenoSet"),
          function(ds, strict=TRUE) {
            if (strict) {
              if ( ! all( chrNames(ds) == chrOrder( chrNames(ds) ) ) ) {
                return(FALSE)
              }
            }
            # Check each chr for ordered start
            chr.ind = chrIndices(ds)
            return(!any(unlist(tapply, start(ds), chr(ds), is.unsorted)))
          })

##' @rdname isGenomeOrder-methods
setMethod("isGenomeOrder",signature=signature(ds="GRanges"),
          function(ds, strict=TRUE) {
            if ( any(duplicated(runValue(seqnames(ds)))) ) { return(FALSE) }
            if (strict == TRUE) {
              if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
                return(FALSE)
              }
            }
            chr.ind = chrIndices(ds)
            return(!any(unlist(tapply( start(ds), chr(ds), FUN=is.unsorted))))
          })

##' Set a GRanges or GenoSet to genome order
##'
##' Returns a re-ordered object sorted by chromosome and start position. If strict=TRUE, then
##' chromosomes must be in order specified by chrOrder.
##' If ds is already ordered, no re-ordering is done. Therefore, checking order with isGenomeOrder,
##' is unnecessary if order will be corrected if isGenomeOrder is FALSE.
##'
##' toGenomeOrder for GRanges differs from sort in that it orders by chromsome and start position only,
##' rather than chromsome, strand, start, and width.
##' 
##' @param ds GenoSet or GRanges
##' @param strict logical, should chromosomes be in order specified by chrOrder?
##' @return re-ordered ds
##' @export toGenomeOrder
##' @examples
##'   data(genoset)
##'   toGenomeOrder( genoset.ds, strict=TRUE )
##'   toGenomeOrder( genoset.ds, strict=FALSE )
##'   toGenomeOrder( locData(genoset.ds) )
##' @docType methods
##' @family "genome ordering"
##' @rdname toGenomeOrder-methods
setGeneric("toGenomeOrder", function(ds, strict=TRUE) standardGeneric("toGenomeOrder"))

##' @rdname toGenomeOrder-methods
setMethod("toGenomeOrder",signature=signature(ds="GRanges"),
          function(ds, strict=TRUE) {
            if (strict == TRUE) {
              if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
                seqlevels(ds) = chrOrder(seqlevels(ds))
              }
            }
            row.order = order(as.integer(seqnames(ds)),start(ds))
            if (is.unsorted(row.order)) {
              ds = ds[row.order,,drop=FALSE]
            }
            return(ds)
          })

##' @rdname toGenomeOrder-methods
setMethod("toGenomeOrder", signature=signature(ds="GenoSet"),
          function(ds,strict=TRUE) {
            locData(ds) = toGenomeOrder(locData(ds),strict=strict) # locData<- fixes row ordering in ds
            return(ds)
          })
