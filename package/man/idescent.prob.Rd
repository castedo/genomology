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

\subsection{recomb.prob}{

Define \eqn{I_k} to be the event that nucleotide at position \eqn{k} is
identical-by-decent (IBD).

The argument \strong{recomb.prob} corresponds to vector \eqn{\rho} with definition
\deqn{ \rho_k = \Pr(\overline{I_k} | I_{k-1}) = \Pr(\overline{I_{k-1}} | I_k) }

You can use the function \strong{recombination.prob} R function to calculate
such recombination probabilities.

For convenience, position zero is modeled as a position infinitely far away such that
for all \eqn{k}
\deqn{ \rho_1 = \Pr(\overline{I_1} | I_0) = \Pr(\overline{I_k}) }
\deqn{ \Pr(I_k) = 1 - \rho_1 }
This is the unconditional probability that a nucleotide is
not identical-by-descent.

A useful formula is
\deqn{ \Pr(I_k | \overline{I_{k-1}}) = \rho_k (1 - \rho_1) / \rho_1 }

}

\subsection{like.idescent and like.no.idescent}{

Let \eqn{W_k} be the event of evidence observed for nucleotides at position \eqn{k}.

The argument \strong{like.idescent} corresponds to vector \eqn{\phi} with definition
\deqn{ \phi_k = \Pr(W_k | I_k) }
and argument \strong{like.no.idescent} corresponds to vector \eqn{\psi} with definition
\deqn{ \psi_k = \Pr(W_k | \overline{I_k}) }

The function \strong{haploids.likeli.idescent} can be used to generate
\strong{like.idescent} values.  For \strong{like.no.idescent}, the default of
just \strong{0.5} is an OK default for comparing one sequence (halploid) to one
other sequence (haploid).

}

\subsection{one-sided conditional identical-by-descent probability}{

We define the event of evidence observed for nucleotides \strong{after} position \eqn{k}:

\deqn{ A_k = \bigcap_{i > k} W_i }

Likewise, we define the event of evidence observed for nucleotides \strong{before} position \eqn{k}:

\deqn{ B_k = \bigcap_{i < k} W_i }

Note that \eqn{B_1} is trivially the event of all possibilities (i.e. there is no \eqn{W_0}).

The output of R function \strong{roll.idescent.prob} is vector \eqn{ \Pr(I_i | B_i) }.

We assume that events \eqn{W_i} and \eqn{B_i} are conditionally independent
given \eqn{J \in \{ I_i, \overline{I_i} \}}, that is

\deqn{ \Pr(W_i \cap B_i | J) = \Pr(W_i | J) \Pr(B_i | J) }

which implies

\deqn{ \Pr(W_i | J) 
     = \frac{\Pr(W_i \cap B_i | J)}{\Pr(B_i | J)}
     = \Pr(W_i | B_i \cap J) }

\deqn{ \Pr(W_i | J) \Pr(J | B_i)
     = \Pr(W_i | B_i \cap J) \Pr(J | B_i) 
     = \Pr(W_i \cap J | B_i) }

which can be applied to calculate

\deqn{ \Pr(W_i \cap I_i | B_i) + \Pr(W_i \cap \overline{I_i} | B_i)
     = \phi_i \Pr(I_i | B_i) + \psi_i \Pr(\overline{I_i} | B_i) }
\deqn{ \Pr(W_i | B_i)
     = \phi_i \Pr(I_i | B_i) + \psi_i (1 - \Pr(I_i | B_i)) }
\deqn{ \Pr(I_i | W_i \cap B_i)
     = \frac{ \Pr(W_i \cap I_i | B_i) }{ \Pr(W_i | B_i) } }
\deqn{ \Pr(I_i | B_{i+1})
     = \frac{\phi_i \Pr(I_i | B_i)}{\phi_i \Pr(I_i | B_i) + \psi_i (1 - \Pr(I_i | B_i))} }

We assume \eqn{B_{i+1}} is independent of \eqn{I_{i+1}} for any given
\eqn{J \in \{ I_i, \overline{I_i} \}}, that is

\deqn{ \Pr(I_{i+1} \cap B_{i+1} | J) = \Pr(I_{i+1} | J) \Pr(B_{i+1} | J) }

thus

\deqn{ \Pr(I_{i+1} | J \cap B_{i+1}) = \Pr(I_{i+1} | J) }

Applying this to the expression of \eqn{\Pr(I_{i+1} | B_{i+1})} we get

\deqn{ \Pr(I_{i+1} | B_{i+1})
     = \Pr(I_{i+1} \cap I_i | B_{i+1})
     + \Pr(I_{i+1} \cap \overline{I_i} | B_{i+1})
}

\deqn{ \Pr(I_{i+1} | B_{i+1})
     = \Pr(I_{i+1} | I_i \cap B_{i+1}) \Pr(I_i | B_{i+1})
     + \Pr(I_{i+1} | \overline{I_i} \cap B_{i+1}) \Pr(\overline{I_i} | B_{i+1})
}

\deqn{ \Pr(I_{i+1} | B_{i+1})
     = \Pr(I_{i+1} | I_i) \Pr(I_i | B_{i+1})
     + \Pr(I_{i+1} | \overline{I_i}) \Pr(\overline{I_i} | B_{i+1})
}

\deqn{ \Pr(I_{i+1} | B_{i+1})
     = (1 - \rho_{i+1}) \Pr(I_i | B_{i+1})
     + (\rho_{i+1} (1 - \rho_1) / \rho_1)(1 - \Pr(I_i | B_{i+1}))
}

}

\subsection{idenscent.prob}{

We assume that \eqn{A = A_k} and \eqn{B = B_k} are conditionally independent
given fixed \eqn{I = I_k}. That is, if the nucleotide at position k is
inherited identically, then the observed nucleotides before and after are
independent of each other.

\deqn{ \Pr(A \cap B | I) = \Pr(A | I) \Pr(B | I) }

therefore

\deqn{ \frac{ \Pr(A \cap B | I) \Pr(I) }{
              \Pr(A) \Pr(B) }
       = \frac{ \Pr(A | I) \Pr(B | I) \Pr(I)^2 }{
                \Pr(A) \Pr(B) \Pr(I) }
}

\deqn{ \frac{ \Pr(A \cap B \cap I) }{
              \Pr(A) \Pr(B) }
     = \frac{ \Pr(I | A) \Pr(I | B) }{ \Pr(I) }
 }

\deqn{ \frac{ \Pr(A \cap B) }{ \Pr(A) \Pr(B) }
     = \frac{ \Pr(I | A) \Pr(I | B) }{ \Pr(I) }
     + \frac{ \Pr(\bar{I} | A) \Pr(\bar{I} | B) }{ \Pr(\bar{I}) }
}

\deqn{ \Pr(I | A \cap B)
     = \frac{ \Pr(I | A) \Pr(I | B) / \Pr(I)
           }{ \Pr(I | A) \Pr(I | B) / \Pr(I) 
            + \Pr(\bar{I} | A) \Pr(\bar{I} | B) / \Pr(\bar{I}) }
}

We assume \eqn{W = W_k} and \eqn{A \cap B} are conditionally independent
given \eqn{J \in \{I_k, \overline{I_k}\}}, that is

\deqn{ \Pr(W \cap A \cap B | J) = \Pr(W | J) \Pr(A \cap B | J) }

and thus

\deqn{ \Pr(W | J) 
     = \frac{\Pr(W \cap A \cap B | J)}{\Pr(A \cap B | J)}
      = \Pr(W | J \cap A \cap B) }

\deqn{ \Pr(J | A \cap B) \Pr(W | J)
     = \Pr(J | A \cap B) \Pr(W | J \cap A \cap B)
     = \Pr(J \cap W | A \cap B) }

Finally,

\deqn{ \Pr(I \cap W | A \cap B) + \Pr(\bar{I} \cap W | A \cap B)
     = \Pr(I | A \cap B) \phi_i + \Pr(\bar{I} | A \cap B) \psi_i }

with which we can get

\deqn{ \Pr(I | W \cap A \cap B)
     = \frac{ \Pr(I \cap W | A \cap B) }{
              \Pr(W | A \cap B) }
     = \frac{\Pr(I | A \cap B) \phi_i}{\Pr(I | A \cap B) \phi_i + \Pr(\bar{I} | A \cap B) \psi_i} }
}

which is the probability of identity-by-descent given all available observed evidence.

}

}{See PDF (latex) version.} }

