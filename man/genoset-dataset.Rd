\name{genoset-datasets}
\docType{data}
\alias{genoset-datasets}
\alias{fake.lrr}
\alias{fake.baf}
\alias{fake.cn}
\alias{fake.pData}
\alias{genoset.ds}
\alias{locData.gr}
\title{Example GenoSet, BAFSet, and CNSet objects and the data to create them.}
\description{Fake LRR, BAF, pData and location data were generated and saved as
  fake.lrr, fake.cn, fake.baf, fake.pData and locData.rd. These were used to
  construct the objects genoset.ds, baf.ds, and cn.ds}
\usage{data(genoset)}
\format{
  \describe{
    \item{fake.lrr}{A matrix with some randomly generated LRR (log2ratio copynumber) data}
    \item{fake.cn}{A matrix with some randomly generated LRR (log2ratio copynumber) data}
    \item{fake.baf}{A matrix with some randomly generated BAF (B-Allele
      Frequency) data}
    \item{fake.pData}{A data.frame of sample annotation to go with fake.lrr
      and fake.baf}
    \item{locData.gr}{A GRanges object describing the genomic locations
      of the probes in fake.baf and fake.lrr}
    \item{genoset.ds}{A GenoSet object created with fake.lrr as the "foo"
      element, locData.rd as the locData, and fake.pData as the phenoData}
  }
}
\source{Fake data generated using rnorm and the like.}
\keyword{datasets}
