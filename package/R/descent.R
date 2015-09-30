
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
                                 like.idescent[i-1],
                                 like.no.idescent[i-1])
  }
  return(ret)
}

idescent.prob <- function(recomb.prob, like.idescent, like.no.idescent=0.5) {
  if (length(like.no.idescent) == 1) {
     like.no.idescent <- rep_len(like.no.idescent, length(like.idescent))
  }
  chance.idescent <- 1 - recomb.prob[1]
  before <- roll.idescent.prob(recomb.prob, like.idescent, like.no.idescent)
  rev.recomb.prob <- c(recomb.prob[1], rev(recomb.prob[-1]))
  after <- rev(roll.idescent.prob(rev.recomb.prob,
                                  rev(like.idescent),
                                  rev(like.no.idescent)))
  yes <- after * before / chance.idescent
  no <- (1 - after) * (1 - before) / (1 - chance.idescent)
  prob <- yes / (yes + no)
  yes <- prob * like.idescent
  no <- (1 - prob) * like.no.idescent
  return(yes / (yes + no))
}

fuzzy <- function(x, p.error=0.01) {
  stopifnot(is.logical(x))
  ret <- x + (1 - 2 * x) * p.error
  ret[is.na(x)] <- 0.5
  return(ret)
}

assign.descent <- function(position, assignment) {
  stopifnot(!is.unsorted(position))
  stopifnot(all(assignment %in% c(-1:1,NA)))
  recomb <- recombination.prob(centimorgan.dist(position))
  like1 <- fuzzy(assignment >= 0)
  like2 <- fuzzy(assignment <= 0)
  return(idescent.prob(recomb, like1, like2))
}

