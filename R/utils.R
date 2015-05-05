##' Take vector or matrix of log2 ratios, convert to copynumber
##'
##' Utility function for converting log2ratio units (zero is normal) to copynumber units (two is normal)
##' @param x numeric data in log2ratio values
##' @return data of same type as "x" transformed into copynumber units
##' @export cn2lr
##' @seealso cn2lr
lr2cn <- function(x) {
  return( 2 ^ (x + 1) )
}

##' Take vector or matrix of copynumber values, convert to log2ratios
##'
##' Utility function for converting copynumber units (2 is normal) to log2ratio units (two is normal). If ploidy
##' is provided lr is log2(cn/ploidy), otherwise log2(cn/2).
##' @param x numeric vector or matrix, or DataFrame with numeric-like columns (Rle typicaly). Assumed to be in copynumber units.
##' @param ploidy numeric, of length ncol(x). Ploidy of each sample.
##' @return data of same type as "x" transformed into log2ratio units
##' @export lr2cn
##' @seealso lr2cn
##' @rdname cn2lr-methods
setGeneric("cn2lr", function(x, ploidy) standardGeneric("cn2lr"))

##' @rdname cn2lr-methods
setMethod("cn2lr", signature(x="numeric"),
          function(x, ploidy) {
            x[x <= 1e-3 & !is.na(x)] = 1e-3
            if (missing(ploidy)){
              new.x = log2(x) - 1
            } else {
              if (length(ploidy) != 1) { stop("ploidy must be of length 1") }
              new.x = log2( x / ploidy)
            }
            return(new.x)
          })

##' @rdname cn2lr-methods
setMethod("cn2lr", signature(x="matrix"),
          function(x, ploidy) {
            x[x <= 1e-3 & !is.na(x)] = 1e-3
            if (missing(ploidy)){
              new.x = log2(x) - 1
            } else {
              if ( ncol(x) != length(ploidy) ) { stop("ploidy must have the length of ncol(x)") }
              new.x = log2( sweep(x, MARGIN=2, STATS=ploidy, FUN="/"))
            }
            return(new.x)
          })

##' @rdname cn2lr-methods
setMethod("cn2lr", signature(x="DataFrame"),
          function(x, ploidy) {
            if (missing(ploidy)){
              res.list = lapply( x, function(y) {
                y[y <= 1e-3 & !is.na(y)] = 1e-3
                log2(y) - 1
              } )
            } else {
              if ( ncol(x) != length(ploidy) ) { stop("ploidy must have the length of ncol(x)") }
              res.list = mapply( FUN=function(y, p) {
                y[y <= 1e-3 & !is.na(y)] = 1e-3
                log2(y/p)
              }, x, ploidy, SIMPLIFY=FALSE )
            }
            new.x = DataFrame( res.list, row.names=row.names(x), check.names=FALSE)
            return(new.x)
          })

##' Correct copy number for GC content
##'
##' Copy number estimates from various platforms show "Genomic Waves" (Diskin et al.,
##' Nucleic Acids Research, 2008, PMID: 18784189) where copy number trends with local GC content.
##' This function regresses copy number on GC percentage and removes the effect
##' (returns residuals). GC content should be smoothed along the genome in wide
##' windows >= 100kb.
##'
##' @param ds numeric matrix of copynumber or log2ratio values, samples in columns
##' @param gc numeric vector, GC percentage for each row of ds, must not have NAs
##' @param retain.mean logical, center on zero or keep same mean?
##' @return numeric matrix, residuals of ds regressed on gc
##' @export gcCorrect
##' @family "gc content"
##' @examples
##'   gc = runif(n=100, min=1, max=100)
##'   ds = rnorm(100) + (0.1 * gc)
##'   gcCorrect(ds, gc)
gcCorrect <- function(ds, gc, retain.mean=TRUE) {
  if (!requireNamespace("stats",quietly=TRUE)) {
    stop("Failed to require stats package.\n")
  }
  ds = na.exclude(ds)
  if ("na.action" %in% names(attributes(ds))) {
    gc = gc[ - attr(ds, "na.action") ]
  }
  mm = cbind(rep.int(1, length(gc)), gc)
  fit = stats::lm.fit(mm, ds)
  fit$na.action = attr(ds, "na.action")
  ds.fixed = stats::residuals(fit)
  if (retain.mean == TRUE) {
    if (is.null(dim(ds))) {
      ds.fixed = ds.fixed + (sum(ds,na.rm=TRUE)/length(ds))
    } else {
      ds.fixed = sweep(ds.fixed,2,colMeans(ds,na.rm=TRUE),FUN="+")
    }
  }
  return(ds.fixed)
}

##' Calculate mBAF from BAF
##'
##' Calculate Mirrored B-Allele Frequence (mBAF) from B-Allele Frequency (BAF) as in
##' Staaf et al., Genome Biology, 2008.  BAF is converted to mBAF by folding around 0.5 so
##' that is then between 0.5 and 1. HOM value are then made NA to leave only HET values that
##' can be easily segmented. Values > hom.cutoff are made NA. Then, if genotypes (usually from
##' a matched normal) are provided as the matrix 'calls' additional HOMs can be set to NA. The
##' argument 'call.pairs' is used to match columns in 'calls' to columns in 'baf'.
##'
##' @param baf numeric matrix of BAF values
##' @param hom.cutoff numeric, values above this cutoff to be made NA (considered HOM)
##' @param calls matrix of NA, CT, AG, etc. genotypes to select HETs (in normals). Dimnames must match baf matrix.
##' @param call.pairs list, names represent target samples for HOMs to set to NA. Values represent columns in "calls" matrix.
##' @return numeric matix of mBAF values
##' @examples
##'    data(genoset)
##'    mbaf = baf2mbaf( genoset.ds[, , "baf"], hom.cutoff=0.9 )
##'    calls = matrix(sample(c("AT","AA","CG","GC","AT","GG"),(nrow(genoset.ds) * 2),replace=TRUE),ncol=2,dimnames=list(rownames(genoset.ds),c("K","L")))
##'    mbaf = baf2mbaf( genoset.ds[, , "baf"], hom.cutoff=0.9, calls = calls, call.pairs = list(K="L",L="L") ) # Sample L is matched normal for tumor sample K, M only uses hom.cutoff
##'    genoset.ds[, ,"mbaf"] = baf2mbaf( genoset.ds[, , "baf"], hom.cutoff=0.9 ) # Put mbaf back into the BAFSet object as a new element
##' @export baf2mbaf
baf2mbaf <- function(baf, hom.cutoff=0.95, calls=NULL, call.pairs=NULL) {
  mbaf = abs(baf[,] - 0.5) + 0.5
  is.na(mbaf) <- mbaf > hom.cutoff

  if (!is.null(calls) && !is.null(call.pairs)) {

    # Use genotypes for/from samples specified by call.pairs to NA some HOMs
    if (! all(names(call.pairs) %in% colnames(baf)) ) {
      stop("call.pairs names and baf colnames mismatch\n")
    }
    if (! all(call.pairs %in% colnames(calls)) ) {
      stop("call.pairs values and calls colnames mismatch\n")
    }
    if ( ! identical( rownames(calls), rownames(baf) ) ) {
      stop("rownames mismatch between calls and baf.")
    }

    # Check row matching between baf and calls.
    # Some calls rows will be missing from mbaf because PennCNV threw those features out.
    # Some rows of mbaf will not be in calls because some arrays have copy-only probes without calls
    # Can't subset both to intersection. Can only subset calls to those in mbaf
    # All baf values HET in calls will be NA, so no false positives
    # False negatives very possible, percent row overlap output will warn at least

    # NA all mbaf data for which there is a genotype and it is not HET

    hom.genotypes = c("AA","CC","GG","TT","AA","BB")
    if (nlevels(calls) > 0) { hom.genotypes = which( levels(calls) %in% hom.genotypes ) }
    is.na(mbaf[rownames(calls),names(call.pairs)]) <- calls[,unlist(call.pairs)] %in% hom.genotypes
  }
  return(mbaf)
}

##' Calculate GC Percentage in windows
##'
##' Local GC content  can be used to remove GC artifacts from copynumber data
##' (see Diskin et al, Nucleic Acids Research, 2008, PMID: 18784189). This
##' function will calculate GC content fraction in expanded windows around
##' a set of ranges following example in
##' http://www.bioconductor.org/help/course-materials/2012/useR2012/Bioconductor-tutorial.pdf. Currently
##' all ranges are tabulated, later I may do letterFrequencyInSlidingWindow for big windows and then
##' match to the nearest.
##' @param object GenomicRanges or GenoSet
##' @param bsgenome BSgenome, like Hsapiens from BSgenome.Hsapiens.UCSC.hg19 or DNAStringSet.
##' @param expand scalar integer, amount to expand each range before calculating gc
##' @param bases character, alphabet to count, usually c("G", "C"), but "N" is useful too
##' @return numeric vector, fraction of nucleotides that are G or C in expanded ranges of \code{object}
##' @examples
##  \dontrun{ data(genoset) }
##' \dontrun{ library(BSgenome.Hsapiens.UCSC.hg19) }
##' \dontrun{ gc = calcGC(genoset.ds, Hsapiens) }
##' @export calcGC
calcGC <- function(object, bsgenome, expand=1e6, bases=c("G", "C")) {
  if (!requireNamespace("BSgenome",quietly=TRUE)) {
    stop("Failed to require BSgenome package.\n")
  }
  if (!requireNamespace("Biostrings",quietly=TRUE)) {
    stop("Failed to require Biostrings package.\n")
  }
  if (is(bsgenome, "GmapGenome")) { bsgenome = as(bsgenome, "DNAStringSet") }
  
  chr.ind = chrIndices(object)
  rownames(chr.ind) = gsub("^chr", "", rownames(chr.ind)) # always take off 'chr' prefix to having 'chrchr'
  if (grepl("^chr", seqlevels(bsgenome)[1])) {
      rownames(chr.ind) = paste0("chr", rownames(chr.ind)) 
  }
  start = start(object) - expand
  end = end(object) + expand
  gc.list = lapply(rownames(chr.ind), function(chr.name) {
    range = seq.int(chr.ind[chr.name, 1], chr.ind[chr.name, 2])
    seq = bsgenome[[chr.name]]
    v = suppressWarnings( Views(seq,  start=start[range], end=end[range]) )
    alf = alphabetFrequency(v, as.prob = TRUE)
    gc = rowSums(alf[, bases, drop=FALSE])
  })
  gc = do.call(c, gc.list)
  return(gc)
}

##' Calculate GC Percentage in sliding window
##'
##' Local GC content  can be used to remove GC artifacts from copynumber data
##' (see Diskin et al, Nucleic Acids Research, 2008, PMID: 18784189). This
##' function will calculate GC content fraction in expanded windows around
##' a set of ranges following example in
##' http://www.bioconductor.org/help/course-materials/2012/useR2012/Bioconductor-tutorial.pdf.
##' Values are as.integer( 1e4 * fraction ) for space reasons.
##' @param dna BSgenome or DNAStringSet
##' @param window scalar integer, calculate GC content in a sliding (by one base) window of this size.
##' @return SimpleRleList, integer 1e4 * GC fraction, chromosomes 1:22, X and Y
##' @examples
##' \dontrun{ library(BSgenome.Hsapiens.UCSC.hg19) }
##' \dontrun{ gc = calcGC2(Hsapiens) }
##' @export
calcGC2 <- function(dna) {
    dna = as(gmapGenome, "DNAStringSet")
    dna = dna[ c(1:22, "X", "Y") ]
    window = 1e6
    gc.list = RleList(
        lapply( dna,
               function(x) {
                   gc = rowSums( letterFrequencyInSlidingView(x, window, c("G", "C"), as.prob=TRUE), na.rm=TRUE )
                   gc = Rle(1e4 * gc)
               }), compress=FALSE)
    return(gc.list)
}

##' Center continuous data on mode
##'
##' Copynumber data distributions are generally multi-modal. It is often assumed that
##' the tallest peak represents "normal" and should therefore be centered on a
##' log2ratio of zero. This function uses the density function to find the mode of
##' the dominant peak and subtracts that value from the input data.
##'
##' @param ds numeric matrix
##' @return numeric matrix
##' @export modeCenter
##' @examples
##'   modeCenter( matrix( rnorm(150, mean=0), ncol=3 ))
modeCenter <- function(ds) {
  if (!requireNamespace("stats",quietly=TRUE)) {
    stop("Failed to require stats package.\n")
  }
  column.modes = apply(ds,2, function(x) {
    l2r.density = stats::density(x,na.rm=TRUE)
    density.max.index = which.max(l2r.density$y)
    return(l2r.density$x[density.max.index])
  })
  ds = sweep(ds, 2, column.modes)
  return(ds)
}

##' Load a GenoSet from a RData file
##'
##' Given a rds file or a rda file with one object (a GenoSet or related object), load it,
##' and return.
##' @param path character, path to rds or rda file
##' @return GenoSet or related object (only object in RData file)
##' @examples
##' \dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
##' \dontrun{ ds = readGenoSet("/path/to/genoset.rda") }
##' \dontrun{ ds = readGenoSet("/path/to/genoset.rds") }
##' @export readGenoSet
readGenoSet <- function(path) {
  header = readLines(path, 1)
  if (grepl("^RD", header)[1] == TRUE) {
    object = get(load(path)[1])
  } else {
    object = readRDS(path)
  }
  if (!is(object,"eSet")) { stop("Loaded object is not an eSet or derived class.") }
  return( object )
}

.simple_rbind_dataframe <- function(dflist, element.colname) {
  numrows = vapply(dflist, nrow, integer(1))
  if (!missing(element.colname)) {
    list.name.col = factor(rep(names(dflist), numrows), levels=names(dflist))
  }
  dflist = dflist[ numrows > 0 ] # ARGH, if some data.frames have zero rows, factors become integers
  myunlist = base::unlist
  mylapply = base::lapply
  cn = names(dflist[[1]])
  inds = structure(1L:length(cn), names=cn)
  big <- mylapply(inds,
                function(x) {
                  myunlist(
#                    mylapply(dflist, function(y) { y[[x]] }),
                    mylapply(dflist, function(y) { .subset2(y, x) }),
                    recursive=FALSE, use.names=FALSE)
                })
  if (!missing(element.colname)) {
    big[[element.colname]] = list.name.col
  }
  class(big) <- "data.frame"
  attr(big, "row.names") <- .set_row_names(length(big[[1]]))
  return(big)
}
