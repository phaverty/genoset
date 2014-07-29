### TO DO
# tidy up genome order tests
# GRanges tests for subsetting

test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]

test_creation <- function() {

  pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  locs=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  bad.locs=GRanges(ranges=IRanges(start=c(5,6,10:7,1:4),width=1,names=probe.names[c(5,6,10:7,1:4)]),seqnames=c(rep("chr3",2),rep("chrX",4),rep("chr1",4)))  

  tom = GenoSet( locData=locs, cn=cn, pData=pData )

  gs.from.ad = GenoSet( locData=locs, assayData=assayDataNew(storage.mode="environment",cn=cn), pData=pData )

  rle.genoset = GenoSet(
    locData=locs,
    lrr=DataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    baf=DataFrame(K=Rle(31:40),L=Rle(41:50),M=Rle(51:60),row.names=probe.names),
    pData=pData)
  
  misordered.genoset = GenoSet( locData=locs, cn=cn[ rev(probe.names), ], foo=cn[ rev(probe.names),], pData=pData[rev(test.sample.names),] )

  bad.locData.genoset = GenoSet( locData=bad.locs, cn=cn, foo=cn, pData=pData )

  checkTrue(validObject(tom),"Regular GenoSet")
  checkTrue(validObject(gs.from.ad),"GenoSet with provided assayData")
  checkTrue(validObject(rle.genoset),"GenoSet with Rle data")
  checkTrue(validObject(misordered.genoset),"Starting with some sample name and feature name misordering")
  checkTrue( identical(misordered.genoset[,,"cn"],cn) && identical(misordered.genoset[,,"foo"],cn))
  checkIdentical( pData, pData(misordered.genoset), "Misordered pData gets fixed" )
  checkTrue(validObject(bad.locData.genoset), "Can fix locData not in strict genome order")
  checkIdentical( toGenomeOrder(locData(bad.locData.genoset),strict=TRUE), locData(bad.locData.genoset), "badly ordered locData gets fixed" )
}

test_featureNames <- function() {
  ds = GenoSet(
    locData=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))), 
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  gr = locData(ds)
  checkEquals( featureNames(ds), probe.names, "Get featureNames from GenoSet" )
  checkEquals( featureNames(ds), featureNames(featureData(ds)), "GenoSet FeatureNames match featureData featureNames" )
  checkEquals( featureNames(ds), rownames(fData(ds)), "GenoSet FeatureNames match fData rownames" )
  checkEquals( featureNames(ds), featureNames(gr), "Get featureNames from GRanges" )
  checkEquals( featureNames(ds), rownames(ds), "featureNames and rownames are the same thing for a GenoSet.")
  checkEquals( featureNames(gr), rownames(gr), "featureNames and rownames are the same thing for a GRanges.")
  checkEquals( names(gr), rownames(gr), "names and rownames are the same thing for a GRanges.")
  new.fnames = paste("f",featureNames(ds),sep="")
  featureNames(ds) = new.fnames
  featureNames(gr) = new.fnames
  checkEquals( featureNames(ds), new.fnames, "Set featureNames in GenoSet")
  checkEquals( featureNames(gr), new.fnames, "Set featureNames in GRanges")
  new.fnames = paste("g",featureNames(ds),sep="")
  rownames(ds) = new.fnames
  rownames(gr) = new.fnames
  checkEquals( rownames(ds), new.fnames, "Set rownames in GenoSet")
  checkEquals( rownames(gr), new.fnames, "Set rownames in GRanges")

}

test_sampleNames <- function() {
  ds = GenoSet(
    locData=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkIdentical( suppressWarnings(sampleNames(ds)), test.sample.names )
  checkIdentical( colnames(ds), test.sample.names )
  colnames(ds) = LETTERS[1:3]
  checkIdentical( colnames(ds), LETTERS[1:3] )
}

test_locData <- function() {
  # With GRanges
  locs.gr = GRanges(ranges=IRanges(start=1:10, width=1, names=probe.names), seqnames=factor(c(rep("chr1", 4), rep("chr3", 2), rep("chrX", 4)), levels=c("chr1", "chr3", "chrX")))
  locs.gr.new = GRanges(ranges=IRanges(start=1:10, width=1, names=probe.names), seqnames=factor(c(rep("chr1", 4), rep("chr3", 2), rep("chrX", 4)), levels=c("chr1", "chr3", "chrX")))
  locs.gr.bad = locs.gr.new
  names(locs.gr.bad)[1] = "FOO"
  ds = GenoSet(
    locData=locs.gr,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  ds.new = GenoSet(
    locData=locs.gr.new,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkEquals(locs.gr,locData(ds))
  checkEquals(ds,ds.new,check.attributes=FALSE)
  checkException( {locData(ds) = locs.gr.bad},silent=TRUE )
  checkTrue({locData(ds) = locs.gr; is.null(names(ds@locData))}, "Setting locData makes locData GRanges names null")
  checkEquals(rownames(locData(ds)), names(locs.gr), "However, getting locData back out resets the GRanges names")
}

test_getters.and.setters <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  pData = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  bob = GenoSet(
    locData=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))), 
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=pData
    )
  checkIdentical(c("baf", "lrr"), names(bob))
}

test_rd.gs.shared.api.and.getting.genome.info <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  point.locData = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.locData.gr = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.bad.chr.order.locData = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr5",4),rep("chrX",2),rep("chr3",4)))
  wide.locData =  GRanges(ranges=IRanges(start=seq(1,30,by=3),width=3,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  gs = GenoSet(
    locData=point.locData,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  gr = as(point.locData,"GRanges")

  checkEquals( start( point.locData ), start( gs ) )
  checkEquals( width( point.locData ), width( gs ) )
  checkEquals( end( point.locData ), end( gs ) )
  checkEquals( chr(point.locData), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  checkEquals( chr(gr), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  checkEquals( chr( point.locData ), chr( gs ) )
  checkEquals( chr( point.locData ), chr( gr ) )
  checkEquals( pos(point.locData), 1L:10L )
  checkEquals( pos(wide.locData), seq(from=2L, length=10, by=3L ) )
  checkEquals( pos( point.locData ), pos( gs ) )
  checkEquals( pos( point.locData ), pos( gr ) )
  checkEquals( chrNames( point.locData ), c("chr1","chr3","chrX") )
  checkEquals( chrNames( point.locData ), chrNames( gs ) )
  checkEquals( chrNames( gr[1:3,] ), c("chr1"), "chrNames on GRanges with empty levels should give just unique values" )
  point.locData2 = point.locData
  chrNames(point.locData2) = sub("chr","",chrNames(point.locData2))
  checkEquals( chrNames( point.locData2 ), c("1","3","X") )
  gs2 = gs
  locData(gs2) = point.locData2
  checkEquals( chrNames( gs2 ), c("1","3","X") )
  checkEquals( chrNames( point.locData ), c("chr1","chr3","chrX") )
  checkEquals( elementLengths( point.locData ), elementLengths( gs ) )
  checkEquals( elementLengths( point.locData ), elementLengths( point.locData.gr ) )
  checkEquals( chrInfo( point.locData ), chrInfo( gs ) )
  checkEquals( chrInfo( point.locData ), chrInfo( gr ) )
  checkEquals( chrInfo( point.locData ), matrix(c(1,5,11,4,10,20,0,4,10),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("start","stop","offset") ) ))
  checkEquals( chrIndices( point.locData, "chr3"), c(5,6) )
  checkException( chrIndices( point.locData, "chrFOO"), silent=TRUE )
  checkEquals( chrIndices( point.locData ), chrIndices( gs ) )
  checkEquals( chrIndices( point.locData ), matrix(c(1,5,7,4,6,10,0,4,6),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("first","last","offset") ) ))
  checkEquals( chrIndices( point.locData ), chrIndices(point.locData.gr) )
  checkEquals( chrIndices( point.locData[1:6,] ), chrIndices(point.locData.gr)[1:2,], "Empty levels ignored" )
  checkEquals( genoPos( point.locData ), genoPos( gs ) )
  checkEquals( genoPos( point.locData ), genoPos( gs ) )

  # Universe
  gr.uni = GRanges(IRanges(start=1:4,width=1),seqnames=c("chr1","chr2","chr3","chr4"))
  genome(gr.uni) = c("hg18","hg19","hg19","hg19")
  genome(gr.uni) = c("hg19")
  geno = rep("hg19", 3)
  locData(gs) = gr
  genome(gs) = geno
  checkEquals(geno, genome(gs), "Get and set genome of GenoSet", checkNames=FALSE)
}

test_subset <- function() {
  test.gr = GRanges(ranges=IRanges(start=8:14,width=1,names=letters[8:14]),seqnames=rep("chrX",7))
  test.pdata = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])),stringsAsFactors=FALSE)
  test.phenodata = phenoData=new("AnnotatedDataFrame",test.pdata)

  test.ds = new("GenoSet",
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    phenoData=test.phenodata
    )
  
  expected.ds = new("GenoSet",
    locData=GRanges(ranges=IRanges(start=8:10,width=1,names=probe.names[8:10]),seqnames="chrX"),
    lrr=matrix(c(8:10,18:20,28:30),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    baf=matrix(c(38:40,48:50,58:60),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    phenoData=test.phenodata
    )
  
  chr3.ds = new("GenoSet",
    locData=GRanges(ranges=IRanges(start=5:6,width=1,names=probe.names[5:6]),seqnames="chr3"),
    lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    phenoData=test.phenodata
    )

  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  
  ds = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=test.pdata,
    annotation="SNP6"
    )

  subset.rows.ds = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))[2:3,,drop=TRUE],
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    pData=test.pdata,
    annotation="SNP6"
    )
  
  subset.cols.ds = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    pData=test.pdata[2:3,],
    annotation="SNP6"
    )

  gene.gr = GRanges(ranges=IRanges(start=2:3,width=1),seqnames=c("chr1","chr1"))
  
  # Subsetting whole object
  checkEquals( ds[ ,2:3], subset.cols.ds, check.attributes=TRUE)
  checkEquals( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  checkEquals( ds[ 2:3, , "baf" ], subset.rows.ds[, , "baf"], check.attributes=TRUE)
  ds.two.rows = ds[ 2:3, ]
  checkEquals( featureNames(assayData(ds.two.rows)), featureNames(locData(ds.two.rows)), "featureNames from locData and assayData should be the same when rows subset.")
  checkEquals( ds[ gene.gr, ], subset.rows.ds, check.attributes=FALSE )
  
  # Subsetting assayData / extracting
  checkEquals( ds[ 5, 3, "baf"], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ 5, 3, 1], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ , , "lrr"], assayDataElement(ds,"lrr"), "Extract whole matrix" )
  checkException( ds[ , , "foo"], "Fail to extract assayDataElement with bad character k", silent=TRUE)
  checkException( ds[ , , 8], "Fail to extract assayDataElement with bad integer k", silent=TRUE)
  
  # Test subsetting by location
  checkEquals( test.ds[test.gr,], expected.ds, check.attributes=FALSE )
  checkEquals( test.ds[as(test.gr,"GRanges"),], expected.ds, check.attributes=FALSE )
  checkEquals( test.ds[8:10,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds , check.attributes=FALSE )

  # Replace
  ds = GenoSet(
    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  ds[,,"baf"] = ds[,,"lrr"]
  checkEquals(ds[,,"baf"],ds[,,"lrr"],"Replace whole element")
  bad.names.lrr = ds[,,"lrr"]
  rownames(bad.names.lrr)[1] = "FOO"
  colnames(bad.names.lrr)[1] = "FOO"
  checkException({ds[,,"baf"] = bad.names.lrr}, "Incoming ad element must have dimnames that matches genoset.", silent=TRUE)
  lrr.mat = ds[,,"lrr"]
  lrr.mat[1:2,1:2] = 5
  ds[1:2,1:2,"lrr"] = 5
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with integer indices")
  lrr.mat[6:8,2] = 3
  ds[locData(ds)[6:8,],2,"lrr"] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with GRanges subsetting of rows")
  ds[,3,"lrr"] = 3
  lrr.mat[,3] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace column")
  ds[, , "FOO"] = ds[, , "lrr"]
  checkEquals(ds[, , "FOO"], ds[, , "lrr"], "Adding a whole new matrix is OK.")
  checkException({ds[1, 1, "foo"] = 5}, "Fail to replace with bad character assayDataElement index k", silent=TRUE)
  checkException({ds[1, 1, 8] = 5}, "Fail to replace with bad integer assayDataElement index k", silent=TRUE)
}

test_subset_w_granges <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  test.gr = GRanges(ranges=IRanges(start=8:14,width=1),names=letters[8:14],seqnames=rep("chrX",7))
  lrr = matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf = matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  pData = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  locs = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  test.ds = GenoSet(
    locData=locs,
    lrr=lrr,
    baf=baf,
    pData=pData
    )
  
  expected.ds = GenoSet(
    locData=locs[8:10,],
    lrr=lrr[8:10,],
    baf=baf[8:10,],
    pData=pData
    )
  
  chr3.ds = GenoSet(
    locData=locs[5:6,],
    lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    pData=pData
    )
  
  ds = GenoSet(
    locData=locs,
    lrr=lrr,
    baf=baf,
    pData=pData
    )

  subset.rows.ds = GenoSet(
    locData=locs[2:3,],
    lrr=lrr[2:3,],
    baf=baf[2:3,],
    pData=pData
    )
  
  subset.cols.ds = GenoSet(
    locData=locs,
    lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    pData=pData[2:3,]
    )

  gene.gr = GRanges(ranges=IRanges(start=2:3,width=1),seqnames=c("chr1","chr1"))
  
  # Subsetting whole object
  checkEquals( ds[ ,2:3], subset.cols.ds, check.attributes=FALSE)
  checkEquals( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  checkEquals( ds[ gene.gr, ], subset.rows.ds, check.attributes=FALSE)
  
  # Subsetting assayData / extracting
  checkEquals( ds[ 5, 3, "baf"], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ 5, 3, 1], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ , , "lrr"], assayDataElement(ds,"lrr"), "Extract whole matrix" )
  
  # Test subsetting by location
  checkEquals( test.ds[test.gr,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[8:10,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds, check.attributes=FALSE)

  # Replace
  ds = GenoSet(
    locData=locs,
    lrr=lrr,
    baf=baf,
    pData=pData
    )

  ds[,,"baf"] = ds[,,"lrr"]
  checkEquals(ds[,,"baf"],ds[,,"lrr"],"Replace whole element")
  bad.names.lrr = ds[,,"lrr"]
  rownames(bad.names.lrr)[1] = "FOO"
  colnames(bad.names.lrr)[1] = "FOO"
  checkException({ds[,,"baf"] = bad.names.lrr}, "Incoming ad element must have dimnames that matches genoset.",silent=TRUE)
  lrr.mat = ds[,,"lrr"]
  lrr.mat[1:2,1:2] = 5
  ds[1:2,1:2,"lrr"] = 5
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with integer indices")
  lrr.mat[6:8,2] = 3
  ds[locData(ds)[6:8,],2,"lrr"] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with RangedData subsetting of rows")
  ds[,3,"lrr"] = 3
  lrr.mat[,3] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace column")
}

test_genomeOrder <- function() {
  chr.names = c(rep("chr1",3),rep("chr2",3),rep("chr10",4))

  ok.locs = GRanges( ranges = IRanges(start=1:10,width=1,names=paste("p",1:10,sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  checkTrue( isGenomeOrder(ok.locs), "Good locs" )

  bad.locs = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  bad.locs.bad.chr = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr2","chr1","chr10")))
  checkTrue( ! isGenomeOrder(bad.locs, strict=TRUE), "Bad within chr, OK chr levels, fail")

  good.ds = GenoSet(
    locData=ok.locs,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  bad.ds = good.ds[ c(6,5,4,3,2,1,10:7),]
  bad.ds.bad.chrs = GenoSet(
    locData=bad.locs.bad.chr,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names))[c(4,6,5,2,3,1,10:7),],
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkTrue(isGenomeOrder(good.ds))
  checkTrue(!isGenomeOrder(bad.ds))
  checkEquals( good.ds, toGenomeOrder(bad.ds,strict=TRUE), check.attributes=FALSE, "GenoSet disordered within chrs" )
  checkEquals( good.ds, toGenomeOrder(bad.ds.bad.chrs,strict=TRUE), check.attributes=FALSE, "GenoSet disordered within chrs, disordered chrs" )

  gr1 = GRanges(ranges=IRanges(start=c(9,1,5,4,6,2),width=1,names=LETTERS[c(9,1,5,4,6,2)]),seqnames=Rle(factor(c("A","B","C","C","B","A"),levels=c("A","C","B"))))
  gr2 = GRanges(ranges=IRanges(start=c(2,9,4,5,1,6),width=1,names=LETTERS[c(2,9,4,5,1,6)]),seqnames=Rle(factor(c("A","A","C","C","B","B"),levels=c("A","C","B"))))
  gr3 = GRanges(ranges=IRanges(start=c(2,9,1,6,4,5),width=1,names=LETTERS[c(2,9,1,6,4,5)]),seqnames=Rle(factor(c("A","A","B","B","C","C"),levels=c("A","B","C"))))
  checkIdentical(toGenomeOrder(gr1,strict=FALSE),gr2,"GRanges with mis-ordered chromosomes, without strict")
  checkIdentical(toGenomeOrder(gr1,strict=TRUE),gr3,"GRanges with mis-ordered chromosomes, with strict")
  checkTrue(isGenomeOrder(gr2,strict=FALSE))
  checkTrue(isGenomeOrder(gr3,strict=TRUE))
  checkTrue(!isGenomeOrder(gr2,strict=TRUE))
  checkTrue(!isGenomeOrder(gr1,strict=TRUE), "Not in blocks by chromsome, strict")
  checkTrue(!isGenomeOrder(gr1,strict=FALSE), "Not in blocks by chromsome, strict")
}
