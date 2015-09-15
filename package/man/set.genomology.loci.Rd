\name{set.genomology.loci}
\alias{set.genomology.loci}
\title{Set genome locus for genomology package function}
\description{Set simple (non-overlapping) genome locus (genome) data to
be used by functions in genomology package}
\usage{
  set.genomology.loci(data)
}
\arguments{
  \item{data}{data frame of genome locus (genes) positions, in format of
              http://ref.castedo.com/dna/loci.csv}
}

\examples{
  set.genomology.loci(read.csv("http://ref.castedo.com/dna/loci.csv", comment.char="#"))
}

