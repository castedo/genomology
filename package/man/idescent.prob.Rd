\name{idescent.prob}
\alias{idescent.prob}
\title{Probability of identical-by-descent}
\usage{
  idescent.prob(recomb.prob, like.idescent, like.no.idescent)
}
\arguments{
  \item{recomb.prob}{vector of probabilities for recombination out of
identical-by-descent from preceding identical-by-descent position}
  \item{like.idecent}{probability of observed evidence at positions assuming identical-by-descent at positions}
  \item{like.no.idecent}{probability of observed evidence at positions assuming NOT identical-by-descent at positions}
}
\section{Math}{ \ifelse{latex}{

The argument \strong{recomb.prob} is a vector of \eqn{\Pr(\overline{I_k} |
I_{k-1})} where \eqn{I_k} is the event that nucleotide at position \eqn{k} is
identical-by-decent (IBD).

You can use the function \strong{recombination.prob} R function to calculate
such recombination probabilities.

For convenience \eqn{\Pr(I_1 | I_0) = \alpha = \Pr(I_k)} for
all \eqn{k}, the unconditional probability that a nucleotide is
identical-by-descent (i.e. "chance identity-by-descent"). Position zero is
modeled as a position infinitely far away.

A useful formula is
\deqn{ \Pr(I_k | \overline{I_{k-1}}) = \Pr(\overline{I_k} | I_{k-1}) \, \alpha / (1 - \alpha) }


}{ See PDF (latex) version.
} } # section

