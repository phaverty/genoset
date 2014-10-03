###########################################################
##### Functions for working with segmented data
#####   Segmenting, converting among various representation
#####     data.frames, GRanges, Rle
###########################################################

##' Make Rle from segments for one sample
##'
##' Take output of CBS, make Rle representing all features in 'locs' ranges. CBS output contains
##' run length and run values for genomic segmetns, which could very directly be converted into a Rle.
##' However, as NA values are often removed, especially for mBAF data, these run lengths do not
##' necessarily cover all features in every sample. Using the start and top positions of each segment
##' and the location of each feature, we can make a Rle that represents all features.
##'
##' @param segs data.frame of segments, formatted as output of segment function from DNAcopy package
##' @param locs GenomicRanges, like locData slot of a GenoSet
##' @return Rle with run lengths and run values covering all features in the data set.
##' @export segs2Rle
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   segs = runCBS( genoset.ds[, , "lrr"], locData(genoset.ds), return.segs=TRUE )
##'   segs2Rle( segs[[1]], locData(genoset.ds) )  # Take a data.frame of segments, say from DNAcopy's segment function, and make Rle's using probe locations in the locs
segs2Rle <- function(segs, locs) {
  if ("num.mark" %in% colnames(segs)) {
    if (sum(segs[,"num.mark"]) == nrow(locs)) {
      return(Rle( segs[,"seg.mean"], segs[,"num.mark"]))
    }
  }
  seg.gr = GRanges( ranges=IRanges(start=segs[,"loc.start"], end=segs[,"loc.end"]),
    seqnames=factor(segs[,"chrom"],levels=chrOrder(unique(as.character(segs$chrom)))), "Value"=segs[,"seg.mean"])
  seg.gr = toGenomeOrder(seg.gr, strict=TRUE)
  bounds = boundingIndicesByChr( seg.gr, locs )
  temp.rle = bounds2Rle( bounds, values(seg.gr)$Value, nrow(locs) )  # Breaks unit test for 1st being NA
  return(temp.rle)
}

##' Given segments, make an RleDataFrame of Rle objects for each sample
##'
##' Take table of segments from CBS, convert DataTable of Rle objects for each sample.
##' @title CBS segments to probe matrix
##' @param seg.list list, list of data frames, one per sample, each is result from CBS
##' @param locs locData from a GenoSet object
##' @return RleDataFrame with nrows same as locs and one column for each sample
##' @export segs2RleDataFrame
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( genoset.ds[, , "lrr"], locData(genoset.ds), return.segs=TRUE )
##'   segs2RleDataFrame( seg.list, locData(genoset.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##' @family segments
segs2RleDataFrame <- function(seg.list, locs) {
  rle.list = lapply(seg.list, segs2Rle, locs)
  rle.data.frame = RleDataFrame(rle.list, row.names=rownames(locs))
  return(rle.data.frame)
}

##' GRanges from segment table
##'
##' GenoSet contains a number of functions that work on segments. Many work on a data.frame of
##' segments, like segTable and runCBS. This function converts one of these tables in a GRanges.
##' The three columns specifying the ranges become the GRanges and all other columns go into the 'mcols'
##' portion of the GRanges object.
##' @param segs data.frame with loc.start, loc.end, and chrom columns, like from segTable or runCBS
##' @export  segs2Granges
##' @return GRanges
##' @family "segmented data"
segs2Granges <- function(segs) {
  other.cols = setdiff(colnames(segs), c("loc.start", "loc.end", "chrom"))
  segs.gr = GRanges(IRanges(start=segs$loc.start, end=segs$loc.end), seqnames=segs$chrom)
  mcols(segs.gr) = segs[, other.cols]
  return(segs.gr)
}

##' Convert Rle objects to tables of segments
##'
##' Like the inverse of segs2Rle and segs2RleDataFrame. Takes a
##' Rle or a RleDataFrame and the locData both from a GenoSet object
##' and makes a list of data.frames each like the result of CBS's
##' segment.  Note the loc.start and loc.stop will correspond
##' exactly to probe locations in locData and the input to
##' segs2RleDataFrame are not necessarily so. For a DataFrame, the
##' argument \code{stack} combines all of the individual data.frames
##' into one large data.frame and adds a "Sample" column of sample ids.
##'
##' For a Rle, the user can provide \code{locs} or \code{chr.ind},
##' \code{start} and \code{stop}.  The latter is surprisingly much faster
##' and this is used in the DataFrame version.
##'
##' @param object Rle or RleDataFrame
##' @param locs GenomicRanges with rows corresponding to rows of df
##' @param chr.ind matrix, like from chrIndices method
##' @param start integer, vector of feature start positions
##' @param end integer, vector of feature end positions
##' @param factor.chr scalar logical, make 'chrom' column a factor?
##' @param ... in generic, for extra args in methods
##' @return one or a list of data.frames with columns chrom, loc.start, loc.end, num.mark, seg.mean
##' @export segTable
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( genoset.ds[, , "lrr"], locData(genoset.ds), return.segs=TRUE )
##'   df = segs2RleDataFrame( seg.list, locData(genoset.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##'   assayDataElement( genoset.ds, "lrr.segs" ) = df
##'   segTable( df, locData(genoset.ds) )
##'   segTable( genoset.ds[ , , "lrr.segs"], locData(genoset.ds) )
##'   segTable( genoset.ds[ , 1, "lrr.segs"], locData(genoset.ds), colnames(genoset.ds)[1] )
##' @docType methods
##' @rdname segTable-methods
setGeneric("segTable", function(object, ...) standardGeneric("segTable"))

##' @rdname segTable-methods
setMethod("segTable", signature(object="Rle"), function(object,locs=NULL,chr.ind=NULL,start=NULL,end=NULL,factor.chr=TRUE) {

  if (!is.null(locs)) {
    chr.ind = chrIndices(locs)
    start = start(locs)
    end = end(locs)
  } else {
    if (is.null(chr.ind) || is.null(start) || is.null(end)) {
      stop("If locs arg is not provided then chr.ind, start, and end must be provided.")
    }
  }

  # Get union of all breakpoints in Rle and chromosomes
  object.ends = cumsum(runLength(object))

  all.ends = sort.int(unique(c(chr.ind[,2],object.ends)))
  all.starts = c(1L,all.ends[-length(all.ends)]+1L)
  num.mark = (all.ends - all.starts) + 1L

  # Look up runValue with binary search on cumsum runValue starts. Starts rather than ends because findInterval is < rather than <=.
  object.starts = c(1L,object.ends[-length(object.ends)]+1L)
  object.vals = runValue(object)[ findInterval( all.ends, object.starts ) ]

  # Assign chrom,start,stop to each segment
  if (factor.chr == TRUE) {
    chrom = factor(rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ],levels=rownames(chr.ind))
  } else {
    chrom = rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ]
  }
  loc.end = end[all.ends]
  loc.start = start[all.starts]

#  sample.seg = data.frame(chrom = chrom, loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, seg.mean = object.vals, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE, check.rows=FALSE)
  sample.seg = list(chrom = chrom, loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, seg.mean = object.vals)
  class(sample.seg) = "data.frame"
  attr(sample.seg, "row.names") = .set_row_names(length(chrom))
  return(sample.seg)
})

##' @rdname segTable-methods
##' @param stack logical, rbind list of segment tables for each sample and add "Sample" column?
setMethod("segTable", signature(object="DataFrame"), function(object,locs,factor.chr=TRUE,stack=FALSE) {
  internal.factor.chr = ifelse(factor.chr == TRUE && stack == FALSE,TRUE,FALSE)
  chr.ind = chrIndices(locs)
  start = start(locs)
  end = end(locs)

  segTable.method = getMethod("segTable", "Rle")
  segs = lapply( object,
    function(x) {
      return(segTable.method(x,chr.ind=chr.ind, start=start, end=end,factor.chr=internal.factor.chr))
    })
  if (stack == FALSE) {
    return(segs)
  } else {
    segs.df = .simple_rbind_dataframe(segs, "Sample")
    if (factor.chr == TRUE) {
      chr.names = chrNames(locs)
      segs.df$chrom = factor(segs.df$chrom,levels=chr.names)
    }
    return(segs.df)
  }
})


##' Convert Rle objects to tables of segments
##'
##' Like segTable, but for two Rle objects. Takes a
##' pair of Rle or DataFrames with Rle
##' columns and makes one or more data.frames with bounds of each new
##' segment.  Rle objects are broken up so that each resulting segment
##' has one value from each Rle. For a DataFrame, the
##' argument \code{stack} combines all of the individual data.frames
##' into one large data.frame and adds a "Sample" column of sample ids.
##'
##' For a Rle, the user can provide \code{locs} or \code{chr.ind},
##' \code{start} and \code{stop}.  The latter is surprisingly much faster
##' and this is used in the DataFrame version.
##'
##' @param x Rle or list/DataFrame of Rle vectors
##' @param y Rle or list/DataFrame of Rle vectors
##' @param ... in generic, extra arguments for methods
##' @param locs GenomicRanges with rows corresponding to rows of df
##' @param chr.ind matrix, like from chrIndices method
##' @param start integer, vector of feature start positions
##' @param end integer, vector of feature end positions
##' @param factor.chr scalar logical, make 'chrom' column a factor?
##' @return one or a list of data.frames with columns chrom, loc.start, loc.end, num.mark, seg.mean
##' @export segPairTable
##' @family "segmented data"
##' @examples
##'   cn = Rle(c(3,4,5,6),rep(3,4))
##'   loh = Rle(c(2,4,6,8,10,12),rep(2,6))
##'   start = c(9:11,4:9,15:17)
##'   end = start
##'   locs = GRanges(IRanges(start=start,end=end),seqnames=c(rep("chr1",3),rep("chr2",6),rep("chr3",3)))
##'   segPairTable(cn,loh,locs)
##' @docType methods
##' @rdname segPairTable-methods
setGeneric("segPairTable", function(x, y, ...) standardGeneric("segPairTable"))

##' @rdname segPairTable-methods
setMethod("segPairTable", signature(x="Rle",y="Rle"), function(x,y,locs=NULL,chr.ind=NULL,start=NULL,end=NULL,factor.chr=TRUE) {
  # Fill in missing args if locs given
  # Maybe use ... rather than x and y and get names from that to use in colnames
  if (!is.null(locs)) {
    chr.ind = chrIndices(locs)
    start = start(locs)
    end = end(locs)
  } else {
    if (is.null(chr.ind) || is.null(start) || is.null(end)) {
      stop("If locs arg is not provided then chr.ind, start, and end must be provided.")
    }
  }

  # Get union of all breakpoints in two Rles and chromosomes
  x.ends = cumsum(runLength(x))
  y.ends = cumsum(runLength(y))

  all.ends = sort.int(unique(c(chr.ind[,2],x.ends,y.ends)))
  all.starts = c(1L,all.ends[-length(all.ends)]+1L)
  num.mark = (all.ends - all.starts) + 1L

  # Look up runValue with binary search on cumsum runValue starts. Starts rather than ends because findInterval is < rather than <=.
  x.starts = c(1L,x.ends[-length(x.ends)]+1L)
  x.vals = runValue(x)[ findInterval( all.ends, x.starts ) ]
  y.starts = c(1L,y.ends[-length(y.ends)]+1L)
  y.vals = runValue(y)[ findInterval( all.ends, y.starts ) ]

  # Assign chrom,start,stop to each segment
  if (factor.chr == TRUE) {
    chrom = factor(rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ],levels=rownames(chr.ind))
  } else {
    chrom = rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ]
  }
  loc.end = end[all.ends]
  loc.start = start[all.starts]

  # Make output object
  sample.seg = list(chrom=chrom,loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, x=x.vals, y=y.vals)
  class(sample.seg) = "data.frame"
  attr(sample.seg, "row.names") = .set_row_names(length(chrom))
  return(sample.seg)
})

##' @rdname segPairTable-methods
##' @param stack logical, rbind list of segment tables for each sample and add "Sample" column?
setMethod("segPairTable", signature(x="DataFrame",y="DataFrame"), function(x,y,locs,stack=FALSE,factor.chr=TRUE) {
  internal.factor.chr = ifelse(factor.chr == TRUE && stack == FALSE,TRUE,FALSE)
  chr.ind = chrIndices(locs)
  start = start(locs)
  end = end(locs)
  segs = mapply(
    function(one,two) {
      return(segPairTable(one,two,chr.ind=chr.ind, start=start, end=end,factor.chr=internal.factor.chr))
    },
    x,y,
    SIMPLIFY=FALSE
    )
  if (stack == FALSE) {
    return(segs)
  } else {
    segs.df = .simple_rbind_dataframe(segs, "Sample")
    if (factor.chr == TRUE) {
      chr.names = chrNames(locs)
      segs.df$chrom = factor(segs.df$chrom,levels=chr.names)
    }
    return(segs.df)
  }
})

##' Fix NA runs in a Rle
##'
##' Fix NA runs in a Rle when the adjacent runs have equal values
##' @param x Rle to be fixed
##' @param max.na.run integer, longest run of NAs that will be fixed
##' @return Rle
##' @export fixSegNAs
fixSegNAs <- function(x,max.na.run=3) {
  if (is.na(runValue(x)[1]) & runLength(x)[1] <= max.na.run) {
    runValue(x)[1] = runValue(x)[2]
  }
  if (is.na(runValue(x)[nrun(x)]) & runLength(x)[nrun(x)] <= max.na.run) {
    runValue(x)[nrun(x)] = runValue(x)[nrun(x)-1]
  }
  bad = which(is.na(runValue(x)) & runLength(x) <= max.na.run)
  bad = bad[ runValue(x)[bad-1] == runValue(x)[bad+1] ]
  runValue(x)[bad] = runValue(x)[bad+1]
  return(x)
}

##' Utility function to run CBS's three functions on one or more samples
##'
##' Takes care of running CBS segmentation on one or more samples. Makes appropriate
##' input, smooths outliers, and segment
##'
##' @title Run CBS Segmentation
##' @aliases runCBS
##' @param data numeric matrix with continuous data in one or more columns
##' @param locs GenomicRanges, like locData slot of GenoSet
##' @param return.segs logical, if true list of segment data.frames return, otherwise a DataFrame of Rle vectors. One Rle per sample.
##' @param n.cores numeric, number of cores to ask mclapply to use
##' @param smooth.region number of positions to left and right of individual positions to consider when smoothing single point outliers
##' @param outlier.SD.scale number of SD single points must exceed smooth.region to be considered an outlier
##' @param smooth.SD.scale floor used to reset single point outliers
##' @param trim fraction of sample to smooth
##' @param alpha pvalue cutoff for calling a breakpoint
##' @return data frame of segments from CBS
##' @family "segmented data"
##' @export runCBS
##' @examples
##'     sample.names = paste("a",1:2,sep="")
##'     probe.names =  paste("p",1:30,sep="")
##'     ds = matrix(c(c(rep(5,20),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
##'     locs = GRanges(ranges=IRanges(start=c(1:20,1:10),width=1,names=probe.names),seqnames=paste("chr",c(rep(1,20),rep(2,10)),sep=""))
##'
##'     seg.rle.result = RleDataFrame( a1 = Rle(c(rep(5,20),rep(3,10))), a2 = Rle(c(rep(2,10),rep(7,10),rep(9,10))), row.names=probe.names )
##'     seg.list.result = list(
##'       a1 = data.frame( ID=rep("a1",2), chrom=factor(c("chr1","chr2")), loc.start=c(1,1), loc.end=c(20,10), num.mark=c(20,10), seg.mean=c(5,3), stringsAsFactors=FALSE),
##'       a2 = data.frame( ID=rep("a2",3), chrom=factor(c("chr1","chr1","chr2")), loc.start=c(1,11,1), loc.end=c(10,20,10), num.mark=c(10,10,10), seg.mean=c(2,7,9), stringsAsFactors=FALSE)
##'       )
##'
##'     runCBS(ds,locs)  # Should give seg.rle.result
##'     runCBS(ds,locs,return.segs=TRUE) # Should give seg.list.result
runCBS <- function(data, locs, return.segs=FALSE, n.cores=1, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, alpha=0.001) {
  if (!requireNamespace("DNAcopy",quietly=TRUE)) {
    stop("Failed to require DNAcopy package.\n")
  }
  sample.name.list = colnames(data)
  names(sample.name.list) = sample.name.list
  loc.start = as.numeric(start(locs))
  loc.chr = chr(locs)
  presorted = isGenomeOrder(locs,strict=TRUE)

  # mclapply over samples. cbs can loop over the columns of data, but want to use multiple forks
  if (n.cores > 1 && is.loaded("mc_fork", PACKAGE="parallel")) {
    mcLapply <- get('mclapply', envir=getNamespace('parallel'))
    loopFunc = function(...) { mcLapply(...,mc.cores=n.cores, mc.preschedule=FALSE) }
    cat("Using mclapply for segmentation ...\n")
  } else {
    loopFunc = lapply
  }
  segs = loopFunc(sample.name.list,
    function(sample.name) {
      writeLines(paste("Working on segmentation for sample number",match(sample.name,sample.name.list),":",sample.name))
      temp.data = as.numeric(data[,sample.name,drop=TRUE])
      ok.indices = !is.na(temp.data)
      CNA.object <- DNAcopy::CNA(temp.data[ok.indices], loc.chr[ok.indices], loc.start[ok.indices], data.type = "logratio", sampleid = sample.name, presorted=presorted)
      smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object, smooth.region=smooth.region, outlier.SD.scale=outlier.SD.scale, smooth.SD.scale=smooth.SD.scale, trim=trim)
      segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, verbose=0, alpha=alpha)$output
      if (return.segs == TRUE) {
        segment.smoothed.CNA.object$chrom = factor(as.character(segment.smoothed.CNA.object$chrom),levels=chrNames(locs))
        return(segment.smoothed.CNA.object)
      } else {
        return(segs2Rle(segment.smoothed.CNA.object,locs))
      }
    })

  if (return.segs == TRUE) {
    return(segs)
  } else {
    return( RleDataFrame(segs, row.names=rownames(locs) ) )
  }
}

##' Get segment widths
##'
##' The width of a genomic segment helps inform us about the importance of a copy number value. Focal amplifications
##' are more interesting than broad gains, for example. Given a range of interesting regions (i.e. genes) this
##' function determines all genomics segments covered by each gene and returns the average length of the
##' segments covered by each gene in each sample. Often only a single segment covers a given gene in
##' a given sample.
##' @param range.gr GRanges, genome regions of interest, usually genes
##' @param segs data.frame of segments, like from segTable, or a list of these
##' @export rangeSegMeanLength
##' @return named vector of lengths, one per item in range.gr, or a range x length(segs) of these if segs is also list-like.
##' @family "segmented data"
##' @rdname rangeSegMeanLength-methods
setGeneric("rangeSegMeanLength", function(range.gr,segs) standardGeneric("rangeSegMeanLength"))

##' @rdname rangeSegMeanLength-methods
setMethod("rangeSegMeanLength", signature=signature(range.gr="GRanges", segs="list"),
  function(range.gr, segs) {
    vapply(segs, function(x) { .rangeSegMeanLength(range.gr, x) }, numeric(length(range.gr)) )
  })

##' @rdname rangeSegMeanLength-methods
setMethod("rangeSegMeanLength", signature=signature(range.gr="GRanges", segs="data.frame"),
  function(range.gr, segs) {
    .rangeSegMeanLength(range.gr, segs)
})

# Internal function for rangeSegMeanLength methods
.rangeSegMeanLength <- function(range.gr, segs) {
  segs.gr = segs2Granges(segs)
  bounds = boundingIndicesByChr(range.gr, segs.gr)
  rangeMeans(width(segs.gr), bounds)
}
