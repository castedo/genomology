\name{read.ancestrydna.web}
\alias{read.ancestrydna.web}
\alias{read.ancestrydna.zip}
\alias{read.ancestrydna}
\title{Read AncestryDNA file}
\usage{
  read.ancestrydna.zip(file.zip)
  read.ancestrydna.web(url)
  read.ancestrydna(file)
}
\arguments{
  \item{file.zip}{file path to zip file from AncestryDNA}
  \item{url}{URL web address to raw DNA zip file downloaded from AncestryDNA}
  \item{file}{file path to unzipped file from AncestryDNA}
}

\examples{
  mydna <- read.ancestrydna.web("http://example.com/dna-data.zip")
}

