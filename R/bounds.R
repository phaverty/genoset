##' Convert bounding indices into a Rle
##'
##' Given a matrix of first/last indices, like from boundingIndicesByChr, and values for
##' each range, convert to a Rle.  This function takes the expected length of the Rle, n,
##' so that any portion of the full length not covered by a first/last range will be a
##' run with the value NA.  This is typical in the case where data is segmented with CBS
##' and some of the data to be segmented is NA.
##' @export bounds2Rle
##' @param bounds matrix, two columns, with first and last index, like from boundingIndicesByChr
##' @param values ANY, some value to be associated with each range, like segmented copy number.
##' @param n integer, the expected length of the Rle, i.e. the number of features in the
##' genome/target ranges processed by boundingIndicesByChr.
##' @return Rle
##' @family "segmented data"
bounds2Rle <- function( bounds, values, n ) {
  if ( length(values) != nrow(bounds) ) {
    stop("must have one value for each bound")
  }
  if (n < length(values)) {
    stop("n must be >= length(values)")
  }
  run.length = integer((2*nrow(bounds))+1)
  run.length[1] = bounds[1] - 1
  run.value = rep(NA_real_, (2*nrow(bounds))+1)
  data.indices = seq.int(from=2, by=2, length.out=nrow(bounds))
  widths = (bounds[, 2] - bounds[, 1]) + 1
  run.length[data.indices] = widths
  run.value[data.indices] = values
  run.length[data.indices+1] = diff(c(bounds[, 2], n)) - c(widths[-1], 0)

  if (sum(run.length) != n) {
    stop("Rle is the wrong length. Look for double counted features in your bounds table.")
  }
  return( Rle( run.value, run.length ) )
}

##' Find indices of features bounding a set of chromosome ranges/genes
##'
##' This function is similar to findOverlaps but it guarantees at least two features will be
##' covered. This is useful in the case of finding features corresponding to a set of genes.
##' Some genes will fall entirely between two features and thus would not return any ranges
##' with findOverlaps. Specifically, this function will find the indices of the features
##' (first and last) bounding the ends of a range/gene (start and stop) such that
##' first <= start < stop <= last. Equality is necessary so that multiple conversions between
##' indices and genomic positions will not expand with each conversion. Ranges/genes that are
##' outside the range of feature positions will be given the indices of the corresponding
##' first or last index rather than 0 or n + 1 so that genes can always be connected to some data.
##'
##' This function uses some tricks from findIntervals, where is for k queries and n features it
##' is O(k * log(n)) generally and ~O(k) for sorted queries. Therefore will be dramatically
##' faster for sets of query genes that are sorted by start position within each chromosome.
##' The index of the stop position for each gene is found using the left bound from the start
##' of the gene reducing the search space for the stop position somewhat. boundingIndices does not
##' check for NAs or unsorted data in the subject positions. These assumptions are safe for
##' position info coming from a GenoSet or GRanges.
##'
##' @param starts integer vector of first base position of each query range
##' @param stops integer vector of last base position of each query range
##' @param positions Base positions in which to search
##' @param all.indices logical, return a list containing full sequence of indices for each query
##' @return integer matrix of 2 columms for start and stop index of range in data or a list of full sequences of indices for each query (see all.indices argument)
##' @family "range summaries"
##' @export boundingIndices
##' @examples
##'   starts = seq(10,100,10)
##'   boundingIndices( starts=starts, stops=starts+5, positions = 1:100 )
boundingIndices <- function(starts, stops, positions, all.indices=FALSE) {
  if (length(starts) != length(stops)) {
    stop("starts and stops must be the same length.")
  }
  bounds = .Call("binary_bound", as.integer(starts), as.integer(stops), as.integer(positions), PACKAGE="genoset")

  if (all.indices == TRUE) { # Return all covered and bounding indices
    return( apply( bounds, 1, function(x) { seq(from=x[1], to=x[2]) }) )
  } else {  # Just return left and right indices
    return(bounds)
  }

}

##' Find indices of features bounding a set of chromosome ranges/genes, across chromosomes
##'
##' Finds subject ranges corresponding to a set of genes (query ranges), taking chromosome
##' into account. Specifically, this function will find the indices of the features
##' (first and last) bounding the ends of a range/gene (start and stop) such that
##' first <= start < stop <= last. Equality is necessary so that multiple conversions between
##' indices and genomic positions will not expand with each conversion. Ranges/genes that are
##' outside the range of feature positions will be given the indices of the corresponding
##' first or last index on that chromosome, rather than 0 or n + 1 so that genes can always be
##' connected to some data. Checking the left and right bound for equality will tell you when
##' a query is off the end of a chromosome.
##'
##' This function uses some tricks from findIntervals, where is for k queries and n features it
##' is O(k * log(n)) generally and ~O(k) for sorted queries. Therefore will be dramatically
##' faster for sets of query genes that are sorted by start position within each chromosome.
##' The index of the stop position for each gene is found using the left bound from the start
##' of the gene reducing the search space for the stop position somewhat.
##'
##' This function differs from boundingIndices in that 1. it uses both start and end positions for
##' the subject, and 2. query and subject start and end positions are processed in blocks corresponding
##' to chromosomes.
##'
##' Both query and subject must be in at least weak genome order (sorted by start within chromosome blocks).
##'
##' @param query GRanges or something coercible to GRanges
##' @param subject GenomicRanges
##' @return integer matrix with two columns corresponding to indices on left and right bound of queries in subject
##' @export boundingIndicesByChr
##' @family "range summaries"
boundingIndicesByChr <-function(query, subject) {
  if (!is(query,"GRanges")) {
    tryCatch({ query = as(query,"GRanges"); }, error=function(e) { stop("Could not convert query into GRanges.\n") })
  }

  # Subject must have features ordered by start within chromosome. Query need not really, but it's faster.
  # Just checking query genome order to assure data are in blocks by chromosome in a GRanges. Chromosome order doesn't matter.
  if (! isGenomeOrder(subject,strict=FALSE) ) {
    stop("subject must be in genome order.\n")
  }
  if (! isGenomeOrder(query,strict=FALSE) ) {
    stop("query must be in genome order.\n")
  }
  query.chr.indices = chrIndices(query)
  subject.chr.indices = chrIndices(subject)
  if (! all(rownames(query.chr.indices) %in% rownames(subject.chr.indices))) {
    stop("Some query chromosomes not represented in subject. Try query = keepSeqlevels( query, value=chrNames(subject) )")
  }
  subject.chr.indices = subject.chr.indices[rownames(query.chr.indices),,drop=FALSE]
  nquery = as.integer(sum(query.chr.indices[,2] - query.chr.indices[,3])) # !!!
  query.start = start(query)
  query.end = end(query)
  query.names = names(query)
  if (is.null(query.names)) { query.names = as.character(seq.int(from=1,to=nquery)) }
  subject.start = start(subject)
  subject.end = end(subject)
  return(.Call("binary_bound_by_chr", nquery, query.chr.indices, query.start, query.end, query.names, subject.chr.indices, subject.start, subject.end, PACKAGE="genoset"))
}

##' Average features in ranges per sample
##'
##' This function takes per-feature genomic data and returns averages for each of a set of genomic ranges.
##' The most obvious application is determining the copy number of a set of genes. The features
##' corresponding to each gene are determined with boundingIndices such that all features with the bounds
##' of a gene (overlaps). The features on either side of the gene unless those positions
##' exactly match the first or last base covered by the gene.  Therefore, genes falling between two features
##' will at least cover two features. Range bounding is performed by the boundingIndices function.
##'
##' @param query GRanges object representing genomic regions (genes) to be averaged.
##' @param subject A GenoSet object or derivative
##' @param assay.element character, name of element in assayData to use to extract data
##' @param na.rm scalar logical, ignore NAs?
##' @return numeric matrix of features in each range averaged by sample
##' @family "range summaries"
##' @export rangeSampleMeans
##' @examples
##'   data(genoset)
##'   my.genes = GRanges( ranges=IRanges(start=c(35e6,128e6),end=c(37e6,129e6),names=c("HER2","CMYC")), seqnames=c("chr17","chr8") )
##'   rangeSampleMeans( my.genes, genoset.ds, "lrr" )
rangeSampleMeans <- function(query, subject, assay.element, na.rm=FALSE) {
  ## Find feature bounds of each query in subject genoset, get feature data average for each sample
  all.indices = boundingIndicesByChr(query, subject)
  data.matrix = assayDataElement(subject,assay.element)
  range.means = rangeMeans(data.matrix, all.indices, na.rm=na.rm)
  return(range.means)
}
