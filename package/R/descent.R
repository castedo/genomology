
centimorgan.dist <- function(position, ave.centimorgan=1e6) {
  c(Inf, diff(position)/ave.centimorgan)
}

recombination.prob <- function(centimorgan, num.meioses=1, chance.idescent=0.5) {
  limit <- (1 - chance.idescent)
  morgan <- centimorgan / 100
  return(limit * (1 - exp(-morgan * num.meioses / limit)))
}

idescent.continue.prob <- function(prev.prob,
                                   recomb.prob,
                                   chance.idescent,
                                   likeli.idescent,
                                   likeli.not.idescent)
{
  a <- prev.prob * likeli.idescent
  p <- a / (a + (1 - prev.prob) * likeli.not.idescent)
  q <- recomb.prob * chance.idescent / (1 - chance.idescent)
  return(p * (1 - recomb.prob) + (1 - p) * q)
}

roll.idescent.prob <- function(recomb.prob,
                               likeli.idescent,
                               likeli.not.idescent)
{
  ret <- numeric(length(recomb.prob))
  ret[1] <- recomb.prob[1]
  for (i in 2:length(recomb.prob)) {
    ret[i] <- idescent.continue.prob(ret[i-1],
                                     recomb.prob[i],
                                     ret[1],
                                     likeli.idescent[i],
                                     likeli.not.idescent[i])
  }
  return(ret)
}

haploids.likeli.idescent <- function(hap1, hap2, p.error=0.005) {
  match <- (hap1 == hap2)
  return(match + (1 - 2 * match) * p.error)
}

