\name{read.23andme}
\alias{read.23andme.web}
\alias{read.23andme.zip}
\alias{read.23andme}
\title{Read 23andme file}
\usage{
  read.23andme.zip(file.zip)
  read.23andme.web(url)
  read.23andme(file)
}
\arguments{
  \item{file.zip}{file path to zip file from 23andme}
  \item{url}{URL web address to raw DNA zip file downloaded from AncestryDNA}
  \item{file}{file path to unzipped file from 23andme}
}

\examples{
  mydna <- read.23andme.web("http://example.com/dna-data.zip")
}

