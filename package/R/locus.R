
dataEnv <- new.env()

set.genomology.loci <- function(df) {
  stopifnot(is(df, "data.frame"))
  stopifnot(all(c("locus", "chromosome", "start", "stop") %in% names(df)))
  stopifnot(all(is.numeric(df$start)))
  stopifnot(all(is.numeric(df$stop)))
  assign("loci", df, envir=dataEnv)
}

get.genomology.loci <- function() {
  if (!exists("loci", envir=dataEnv)) stop("must call set.genomology.loci")
  get("loci", envir=dataEnv)
}

prob.triage <- function(prob, cutoff=0.05) {
  as.integer(sign(prob - 0.5) * (prob < cutoff | prob > 1 - cutoff))
}

chromorder22 <- function(chr, pos) { which(chr <= 22) }

par1.end <- 2.699e6
par2.start <- 154.932e6

chromorderX <- function(chr, pos) {
  c(which(chr == 25 & pos < par1.end),
    which(chr == 23),
    which(chr == 25 & pos > par2.start))
}

chromorder23 <- function(chr, pos) {
  c(which(chr <= 22), chromorderX(chr,pos))
}

chromorderY <- function(chr, pos) {
  c(which(chr == 25 & pos < par1.end),
    which(chr == 24),
    which(chr == 25 & pos > par2.start))
}

# convert GRCh37.p13 X (PAR) chromosome positions to Y (PAR)
as.Yposition <- function(chromosome, position) {
  ret <- position %% 95896994 # X Y chromsome size diff
  par1 <- (chromosome == 25 & position < 1e7)
  ret[par1] <- ret[par1] - 50e3
  ret[chromosome < 24] <- NA
  return(ret)
}

chromosome.size <- c(
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
    159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392,  90354753,  81195210,  78077248,
     59128983,  63025520,  48129895,  51304566, 155270560,  59373566)

chromordinate <- function(chromo, position) {
  n <- as.integer(chromo)
  return(position + 2e9 * (pmin(23, n) + (n == 24)))
}

chromorange <- function(chromo, position) {
  n <- as.integer(chromo)
  stopifnot(all(n != 24))
  chromosome.end <- chromordinate(chromo, chromosome.size[pmin(23, n)])
  chromosome.begin <- chromordinate(chromo, 0)
  ord <- chromordinate(chromo, position)
  stopifnot(!is.unsorted(ord))
  mids <- 0.5 * (c(-Inf, ord) + c(ord, Inf))
  after <- pmin(mids[-1], chromosome.end)
  before <- pmax(mids[-length(mids)], chromosome.begin)
  return(after - before)
}

PAR2.startY <- 59034050

chromorangeY <- function(positionY) {
  lower <- c(-Inf, positionY)
  upper <- c(positionY, Inf)
  mids <- 0.5 * (lower + upper)
  i <- which(lower <= PAR2.startY & upper >= PAR2.startY)
  mids[i] <- PAR2.startY
  after <- pmin(mids[-1], chromosome.size[24])
  before <- pmax(mids[-length(mids)], 0)
  return(after - before)
}

chromozone.levels <- c(1:22, "Xonly", "Yonly", "PAR1", "PAR2")

chromozone <- function(chromosome, position) {
  n <- chromosome + (chromosome == 25 & position > 1e7)
  return(factor(n, 1:26, chromozone.levels))
}

genoposition <- function(chromosome, position) { chromosome * 2^28 + position }

locus.at <- function(chromosome, position) {
  stopifnot(length(chromosome) == length(position))
  stopifnot(is.numeric(chromosome))
  stopifnot(is.numeric(position))
  loci <- get.genomology.loci()
  goff <- genoposition(chromosome, position)
  start <- genoposition(loci$chromosome, loci$start)
  stop <- genoposition(loci$chromosome, loci$stop)
  gind <- findInterval(goff, start)
  past <- findInterval(goff, stop + 1)
  gind[gind != (1 + past)] <- NA
  return(loci$locus[gind])
}

homozygous.loci <- function(genotype, locus) {
  stopifnot(all.equal(levels(genotype), genotype.levels))
  stopifnot(length(genotype) == length(locus))
  fun <- function(gt) all(homozygous(gt))
  df <- aggregate(genotype, list(locus=locus), fun)
  is.homoz <- df$x
  return(df$locus[is.homoz])
}

