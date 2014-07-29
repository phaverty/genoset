\name{genome}
\alias{genome}
\title{Get and set the genome universe annotation.}
\arguments{
  \item{x}{GenoSet}
}
\value{
  character, e.g. hg19
}
\description{
  Genome version
}
\details{
  The genome positions of the features in locData. The UCSC
  notation (e.g. hg18, hg19, etc.) should be used.
}
\examples{
data(genoset)
  genome(genoset.ds)
  genome(genoset.ds) = "hg19"
}

