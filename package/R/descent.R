
centimorgan.dist <- function(position, ave.centimorgan=1e6) {
  c(Inf, diff(position)/ave.centimorgan)
}

recombination.prob <- function(centimorgan, num.meioses=1, chance.idescent=0.5) {
  limit <- (1 - chance.idescent)
  morgan <- centimorgan / 100
  return(limit * (1 - exp(-morgan * num.meioses / limit)))
}

step.idescent.prob <- function(prev.prob,
                               recomb.prob,
                               chance.idescent,
                               like.idescent,
                               like.no.idescent)
{
  yes <- prev.prob * like.idescent
  no <- (1 - prev.prob) * like.no.idescent
  p <- yes / (yes + no)
  ## recomb.prob = \Pr(\overline{I_k} | I_{k-1})
  ## anti.recomb.prob = \Pr(I_k | \overline{I_{k-1}})
  anti.recomb.prob <- recomb.prob * chance.idescent / (1 - chance.idescent)
  return(p * (1 - recomb.prob) + (1 - p) * anti.recomb.prob)
}

roll.idescent.prob <- function(recomb.prob,
                               like.idescent,
                               like.no.idescent)
{
  ret <- numeric(length(recomb.prob))
  ret[1] <- recomb.prob[1]
  for (i in 2:length(recomb.prob)) {
    ret[i] <- step.idescent.prob(ret[i-1],
                                 recomb.prob[i],
                                 1 - ret[1],
                                 like.idescent[i],
                                 like.no.idescent[i])
  }
  return(ret)
}

idescent.prob <- function(recomb.prob, like.idescent, like.no.idescent=0.5) {
  if (length(like.no.idescent) == 1) {
     like.no.idescent <- rep_len(like.no.idescent, length(like.idescent))
  }
  before <- roll.idescent.prob(recomb.prob, like.idescent, like.no.idescent)
  recomb.prob <- c(recomb.prob[1], rev(recomb.prob[-1]))
  like.idescent <- rev(like.idescent)
  like.no.idescent <- rev(like.no.idescent)
  after <- rev(roll.idescent.prob(recomb.prob, like.idescent, like.no.idescent))
  chance.idescent <- 1 - recomb.prob[1]
  yes <- after * before / chance.idescent
  no <- (1 - after) * (1 - before) / (1 - chance.idescent)
  prob <- yes / (yes + no)
  yes <- prob * like.idescent
  no <- (1 - prob) * (1 - like.no.idescent)
  return(yes / (yes + no))
}

haploids.likeli.idescent <- function(hap1, hap2, p.error=0.005) {
  match <- (hap1 == hap2)
  return(match + (1 - 2 * match) * p.error)
}

