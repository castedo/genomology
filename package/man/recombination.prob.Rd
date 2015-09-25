\name{recombination.prob}
\alias{recombination.prob}
\title{Probability of recombination out of identical-by-descent}
\usage{
  recombination.prob(centimorgan, num.meioses, chance.decent)
}
\arguments{
  \item{centimorgan}{vector of distances in centiMorgans}
  \item{num.meioses}{number of meioses that can happen between two candidate
identical-by-descent segments (e.g. 1 for parent-child, 2 for siblings, 4 for
1st cousins, etc...}
  \item{chance.descent}{chance that two chromosomes from two relatives are
identical-by-descent at a nucleotide position (SNP)}
}
\section{Math}{ \ifelse{latex}{

Let \eqn{d} denote the distance between nucleotides in Morgans.

Let \eqn{m} denote number of meioses.

Let \eqn{\alpha = \Pr(I_i)} denote probability of nucleotides being
identical-by-descent (unconditional IBD SNP probability).

Let \eqn{q} denote the multi-meiosis recombination function.

\deqn{ q(d, m, \alpha)  = (1-\alpha) \left( 1- e^{-d \cdot m / (1-\alpha)} \right) }

Note appropriate properties:

\deqn{ q(0, m, \alpha) = 0 }

\deqn{ q(\infty, m, \alpha) = 1 - \alpha }

\deqn{ \frac{dq}{dd}(0, m, \alpha) = m }

}{ See PDF (latex) version.
} } # section

