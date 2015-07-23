
###########
## Plots ##
###########

##' Plot data along the genome
##'
##' Plot location data and chromosome boundaries from a GenoSet or GRanges object
##' against data from a numeric or Rle. Specifying a chromosome name and optionally a 'xlim'
##' will zoom into one chromosome region. If more than one chromosome is present, the
##' chromosome boundaries will be marked. Alternatively, for a numeric x and a
##' numeric or Rle y, data in y can be plotted at genome positions x. In this case,
##' chromosome boundaries can be taken from the argument locs. If data for y-axis comes
##' from a Rle lines are plotted representing segments. X-axis tickmarks will be labeled
##' with genome positions in the most appropriate units.
##'
##' @section Methods:
##' \describe{
##' 
##' \item{\code{signature(x = "GenoSetOrGenomicRanges", y = "ANY")}}{
##' Plot feature locations and data from one sample.
##' }
##' 
##' \item{\code{signature(x = "numeric", y = "numeric")}}{
##' Plot numeric location and a vector of numeric data.
##' }
##' 
##' \item{\code{signature(x = "numeric", y = "Rle")}}{
##' Plot numeric location and a vector of Rle data. Uses lines for Rle runs.
##' }
##' }
##' 
##' @param x GenoSet (or descendant) or GRanges
##' @param y numeric or Rle
##' @param locs GRanges, like locData slot of GenoSet
##' @param chr Chromosome to plot, NULL by default for full genome
##' @param add Add plot to existing plot
##' @param xlab character, label for x-axis of plot
##' @param ylab character, label for y-axis of plot
##' @param col character, color to plot lines or points
##' @param lwd numeric, line width for segment plots from an Rle
##' @param pch character or numeric, printing character, see points
##' @param xlim integer, length two, bounds for genome positions. Used in conjunction with "chr" to subset data for plotting.
##' @param ... Additional plotting args
##' @return TRUE
##' @export genoPlot
##' @family "genome plots"
##' @examples
##' data(genoset)
##' genoPlot( x=genoset.ds,y=genoset.ds[,1,"lrr"] )
##' genoPlot( genoPos(genoset.ds), genoset.ds[,1,"lrr"], locs=locData(genoset.ds) ) # The same
##' genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
##' @rdname genoPlot-methods
setGeneric("genoPlot", function(x,y,...) { standardGeneric("genoPlot") } )

##' @rdname genoPlot-methods
setMethod("genoPlot",c(x="numeric",y="numeric"),
          function(x, y, add=FALSE, xlab="", ylab="", col="black", locs=NULL, ...) {
            if (add == FALSE) {
              plot(x,y,axes=FALSE,xlab=xlab,ylab=ylab,xaxs="i",col=col,...)
              genomeAxis(locs=locs)
            } else {
              points(x,y,col=col,...)
            }
            return(invisible(TRUE))
          })

##' @rdname genoPlot-methods
setMethod("genoPlot", c(x="numeric",y="Rle"),
          function(x, y, add=FALSE, xlab="", ylab="", col="red", locs=NULL, lwd=2, xlim=NULL, ...) {
            if (add == FALSE) {
              if (is.null(xlim)) {
                xlim=range(x,na.rm=TRUE)
              }
              plot(NA,type="n",xlim=xlim,ylim=range(y,na.rm=TRUE),xlab=xlab,ylab=ylab,xaxs="i",axes=FALSE,...)
              genomeAxis(locs=locs)
            }
            num.mark = runLength(y)
            loc.end.indices = cumsum(num.mark)
            loc.end = x[loc.end.indices]
            loc.start.indices = (loc.end.indices - num.mark) + 1
            loc.start = x[loc.start.indices]
            seg.mean = runValue(y)
            segments(loc.start, seg.mean, loc.end, seg.mean, col=col, lwd=lwd)
            return(invisible(TRUE))
          })

##' @rdname genoPlot-methods
setMethod("genoPlot", signature(x="GenoSetOrGenomicRanges",y="ANY"), function(x, y, chr=NULL, add=FALSE, pch=".", xlab="", ylab="", ...) {
  ## Note: zoom in by subset is much faster (10X) than xlim, so implement a zoom in with subsetting
  # Get position info, subset by chr if necessary
  dot.args = list(...)
  if ( !is.null(chr) ) {
    if ( "xlim" %in% names(dot.args) ) {
      xlim = dot.args[["xlim"]]
      zoom.gr = GRanges(ranges=IRanges(start=xlim[1],end=xlim[2]),seqnames=chr)
      bounds = boundingIndicesByChr(zoom.gr, x)[1,]
      indices = bounds[1]:bounds[2]
      positions = start(x)[indices]
    } else {
      indices = chrIndices(x,chr)
      positions = start(x)[indices]
    }
    y = y[indices]
    locs = NULL
  } else {
    positions = genoPos(x)
    if (length(chrNames(x)) > 1) {
      locs = x
    } else {
      locs = NULL
    }
  }
  genoPlot(positions,y,locs=locs,add=add,xlab=xlab,ylab=ylab,pch=pch,...)
})

##' Label an axis with base positions
##'
##' Label a plot with Mb, kb, bp as appropriate, using tick locations from axTicks
##'
##' @title Label axis with base pair units
##' @param locs GenomicRanges to be used to draw chromosome boundaries, if necessary.  Usually locData slot from a GenoSet.
##' @param side integer side of plot to put axis
##' @param log logical Is axis logged?
##' @param do.other.side logical, label non-genome side with data values at tick marks?
##' @return nothing
##' @export genomeAxis
##' @family "genome plots"
##' @examples
##'   data(genoset)
##'   genoPlot(genoPos(genoset.ds), genoset.ds[,1, "baf"])
##'   genomeAxis( locs=locData(genoset.ds) )  # Add chromosome names and boundaries to a plot assuming genome along x-axis
##'   genomeAxis( locs=locData(genoset.ds), do.other.side=FALSE ) # As above, but do not label y-axis with data values at tickmarks
##'   genomeAxis()           # Add nucleotide position in sensible units assuming genome along x-axis
genomeAxis <- function(locs=NULL, side=1, log=FALSE, do.other.side=TRUE) {
  if (is.null(locs)) {
    label.positions = axTicks(side=side,log=log)
    if ( max(label.positions) > 1e9 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e9),"Gb",sep="")}))
    } else if ( max(label.positions) > 1e6 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e6),"Mb",sep="")}))
    } else if ( max(label.positions) > 1e3 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e3),"kb",sep="")}))
    } else {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.0f",x),"bp",sep="")}))
    }
  } else {
    chr.info = chrInfo(locs)
    abline(v=chr.info[-1,"start"])
    chr.labels = rownames(chr.info)
    mtext(side=rep(c(3,1),len=length(chr.labels)), text=chr.labels, line=0, at=rowMeans(chr.info[,c("start","stop"),drop=FALSE]))
  }
  box()
  if (do.other.side == TRUE) {
    if (side == 1) {
      axis(side=2)
    } else if (side == 2) {
      axis(side=1)
    }
  }
}
