# Tests for functions utilizing boundingIndices
library(RUnit)

test_bounds2Rle <- function() {
  locs = GRanges(IRanges(start=1:20,width=1),seqnames=c(rep("1",5),rep("2",5),rep("3",5),rep("4",5)))

  bounds1 = matrix(c(3,5, 6,7, 9,9, 13,15, 16,19),ncol=2,byrow=TRUE)
  bounds2 = matrix(c(1,3, 4,5, 6,10, 11,15, 16,20), byrow=2, ncol=2)
  bounds3 = matrix(c(1,5, 6,7, 9,9, 13,15, 16,19),ncol=2,byrow=TRUE)
  bounds4 = matrix(c(1,5, 6,7, 7,10, 11,15, 16,20),ncol=2,byrow=TRUE)
  bounds5 = matrix(c(1,5, 6,7, 7,9, 11,15, 16,20),ncol=2,byrow=TRUE)
  bounds6 = matrix(c(1,2, 3,4, 5,5, 6,6, 7, 9), ncol=2, byrow=TRUE)

  rle1 = Rle( c(NA,"A","B",NA,"C",NA,"D","E",NA), c(2,3,2,1,1,3,3,4,1) )
  rle2 = Rle( LETTERS[1:5], c(3,2,5,5,5) )
  rle3 = Rle( c("A","B",NA,"C",NA,"D","E",NA), c(5,2,1,1,3,3,4,1) )

  values1 = as.vector(na.omit(runValue(rle1)))
  values2 = as.vector(na.omit(runValue(rle2)))
  values3 = as.vector(na.omit(runValue(rle3)))
  values4 = LETTERS[1:5]
  values5 = values4

  checkIdentical( rle1, bounds2Rle( bounds1, values1, length(locs) ), "Gaps at beginning, end, middle" )
  checkIdentical( rle2, bounds2Rle( bounds2, values2, length(locs) ), "No NA segments")
  checkIdentical( rle3, bounds2Rle( bounds3, values3, length(locs) ), "Gaps in middle, end" )
  checkException( bounds2Rle( bounds4, values4, length(locs) ), silent=TRUE, "Exception when Rle too long, no NA" )
  checkException( bounds2Rle( bounds5, values4, length(locs) ), silent=TRUE, "Exception when Rle too long, with NA" )
}

test_boundingIndices <- function() {

  # Test with exact matches
  gene.starts = seq( 0, 42, 2)
  gene.stops = gene.starts + 2
  probes = 1:40
  bounds = matrix( c(c(seq(0,40,2),40),c(seq(2,40,2),40,40)), ncol=2, dimnames=list(NULL, c("left", "right")))
  bounds[1] = 1
  bounds[21:22, 1] = 39

  checkEquals( boundingIndices(gene.starts, gene.stops, probes), bounds)

  # Test random order input with some exact matches
  gene.starts = c(6,2,9,1,14,7,50)
  gene.stops = gene.starts + 2
  probes = seq(3,39,3)
  bounds = matrix( c(c(2,1,3,1,4,2,12),c(3,2,4,2,6,3,13)), ncol=2, dimnames=list(NULL, c("left", "right")))

  checkEquals( boundingIndices(gene.starts, gene.stops, probes), bounds)

  # Test with some not matching exactly
  gene.starts = seq(1,16,4)
  gene.stops = seq(1,16,4) + 1
  probes = seq(2,14,3)
  bounds = matrix(c(1,2,3,4,2,3,4,5),ncol=2, dimnames=list(NULL, c("left", "right")))

  checkEquals( boundingIndices(gene.starts, gene.stops, probes), bounds)
}

test_rangeSampleMeans <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]

  subject = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  query.gr = GRanges( ranges=IRanges(start=c(2,3,7,8),width=2,names=c("joe","bob","fred","tom")), seqnames=factor(c("chr1","chr1","chrX","chrX"),levels=c("chr1","chrX")))

  means = matrix(c(32,42,52,33,43,53,37,47,57,38,48,58)+0.5,ncol=nrow(query.gr),nrow=ncol(subject),dimnames=list(sampleNames(subject),rownames(query.gr)))
  means = t(means)
  checkEquals( rangeSampleMeans( query.gr, subject, "cn" ), means)

  rle.genoset = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    cn=RleDataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  rle.means = matrix(c(2.5,3.5,7.5,8.5,12.5,13.5,17.5,18.5,22.5,23.5,27.5,28.5), nrow=nrow(query.gr), ncol=ncol(rle.genoset), dimnames=list(rownames(query.gr),sampleNames(rle.genoset)))
  checkEquals( rangeSampleMeans( query.gr, rle.genoset, "cn", na.rm=TRUE), rle.means, "RleDataFrame")
  checkEquals( rangeSampleMeans( query.gr, subject, "cn", na.rm=TRUE ), means)
}

test_rangeMeans <- function() {
  bounds = matrix(as.integer(c(2,3,3,5,7,8,9,10)),ncol=2,byrow=TRUE)
  x = matrix(31:60,nrow=10,ncol=3)
  means = matrix(c(32.5,34,37.5,39.5,42.5,44,47.5,49.5,52.5,54,57.5,59.5),nrow=nrow(bounds),ncol=ncol(x))
  checkEquals( rangeMeans( x, bounds), means, "Matrix without dimnames")
  checkEquals( rangeMeans( x[,1], bounds), means[,1], "Vector without dimnames")

  rownames(x) = letters[1:nrow(x)]
  colnames(x) = letters[1:ncol(x)]
  rownames(bounds) = LETTERS[1:nrow(bounds)]
  rownames(means) = rownames(bounds)
  colnames(means) = colnames(x)
  checkEquals( rangeMeans(x, bounds), means, "Matrix with dimnames")
  checkEquals( rangeMeans(x[, 1], bounds), means[,1], "Vector without dimnames")

  na.cells = matrix(c(3,1,4,2,8,1,8,3,2,3,3,3),ncol=2,byrow=TRUE)
  x.w.na = x
  x.w.na[ na.cells ] = NA
  x.w.na[ 8,3 ] = NaN
  means.w.na = matrix(c(32,34.5,37,39.5, 42.5,44,47.5,49.5, NA,54.5,57,59.5),nrow=nrow(bounds),ncol=ncol(x),dimnames=dimnames(means))
  means.w.any.na = matrix(c(NA,NA,NA,39.5, 42.5,NA,47.5,49.5, NA,NA,NA,59.5),nrow=nrow(bounds),ncol=ncol(x),dimnames=dimnames(means))
  checkEquals( rangeMeans( x.w.na, bounds, na.rm=TRUE), means.w.na, "Matrix with dimnames and NAs, na.rm=TRUE")
  checkEquals( rangeMeans( x.w.na, bounds, na.rm=FALSE), means.w.any.na, "Matrix with dimnames and NAs, na.rm=FALSE")
}

test_boundingIndicesByChr <- function() {
  subject= GRanges(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=2,names=as.character(1:12)),
    seqnames=c(rep("1",4),rep("2",4),rep("3",4)))
  query = GRanges(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1109,1139,1150),width=2,names=as.character(1:12)),seqnames=c(rep("1",4),rep("2",4),rep("3",4)))
  res = matrix(as.integer(c(1,1, 1,1, 3,4, 4,4, 5,5, 5,5, 7,8, 8,8, 9,9, 9,9, 11,12, 12,12)),byrow=TRUE,ncol=2,dimnames=list(rownames(query),c("left","right")))
  checkIdentical(res, boundingIndicesByChr(query,subject))

  subject2= GRanges(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=2,names=as.character(1:12)),
    seqnames=c(rep("1",4),rep("2",4),rep("5",4)))
  query2 = GRanges(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1109,1139,1150),width=2,names=as.character(1:12)),seqnames=c(rep("1",4),rep("3",4),rep("5",4)))
  query2 = query2[chr(query2) %in% c("1", "5"), ]
  res2 = res[c(1:4,9:12),]
  checkIdentical(res2, boundingIndicesByChr(query2,subject2))

  subject3 = GRanges(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=1,names=as.character(1:12)),
    seqnames=c(rep("1",4),rep("2",4),rep("3",4)))
  query3 = GRanges(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1110,1139,1150),width=c(rep(2,9),1,2,2),names=as.character(1:12)),seqnames=c(rep("1",4),rep("2",4),rep("3",4)))
  res3 = matrix(as.integer(c(1,1, 1,1, 3,4, 4,4, 5,5, 5,5, 7,8, 8,8, 9,9, 9,9, 11,12, 12,12)),byrow=TRUE,ncol=2,dimnames=list(rownames(query3),c("left","right")))
  checkIdentical(res3, boundingIndicesByChr(query3,subject3))

  subject4 = GRanges(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=1,names=as.character(1:12)),
    seqnames=c(rep("1",4),rep("2",4),rep("4",4)))
  query4 = GRanges(
    ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1110,1139,1150),width=c(rep(2,9),1,2,2),names=as.character(1:12)),
    seqnames=c(rep("1",4),rep("2",2),rep("3", 2), rep("4",4)))
  checkException(boundingIndicesByChr(query4,subject4), silent=TRUE, "Not OK to have extra chrs in query.")

  subject5 = GRanges(ranges=IRanges(start=c(1, 2, seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=1,names=c("A", "B", as.character(1:12))),
    seqnames=c(rep("0", 2), rep("1",4),rep("2",4),rep("3",4)))
  query5 = GRanges(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1110,1139,1150),width=c(rep(2,9),1,2,2),names=c(as.character(1:12))),seqnames=c(rep("1",4),rep("2",4),rep("3",4)))

  res5 = matrix(as.integer(c(1,1, 1,1, 3,4, 4,4, 5,5, 5,5, 7,8, 8,8, 9,9, 9,9, 11,12, 12,12) + 2),byrow=TRUE,ncol=2,dimnames=list(rownames(query5),c("left","right")))
  checkIdentical(res5, boundingIndicesByChr(query5,subject5), "OK to have extra chrs in subject.")

}
