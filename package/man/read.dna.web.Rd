\name{read.dna.web}
\alias{read.dna.web}
\alias{read.dna.zip}
\alias{read.ancestrydna}
\title{Read AncestryDNA or 23andme file}
\usage{
  read.dna.web(url)
  read.dna.zip(file.zip)
}
\arguments{
  \item{url}{URL to raw DNA zip file from AncestryDNA or 23andme}
  \item{file.zip}{file path to zip file from AncestryDNA or 23andme}
}

\examples{
  mydna <- read.dna.web("http://example.com/dna-data.zip")
}

