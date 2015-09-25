
centimorgan.dist <- function(position, ave.centimorgan=1e6) {
  c(Inf, diff(position)/ave.centimorgan)
}

recombination.prob <- function(centimorgan, num.meioses=1, chance.descent=0.5) {
  limit <- (1 - chance.descent)
  morgan <- centimorgan / 100
  return(limit * (1 - exp(-morgan * num.meioses / limit)))
}

