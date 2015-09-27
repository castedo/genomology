
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

chromorderX <- function(chr, pos) {
  c(which(chr == 25 & pos < 1e7), which(chr == 23), which(chr == 25 & pos > 1e7))
}

chromorder23 <- function(chr, pos) {
  c(which(chr <= 22), chromorderX(chr,pos))
}

chromorderY <- function(chr, pos) {
  c(which(chr == 25 & pos < 1e7), which(chr == 24), which(chr == 25 & pos > 1e7))
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

chromopos <- function(chromosome, position) {
  position + 2e9 * (pmin(23, chromosome) + (chromosome == 24))
}

chromorange <- function(chromosome, position) {
  idx <- pmin(23, chromosome) + (chromosome == 24)
  chromosome.end <- chromopos(chromosome, chromosome.size[idx])
  chromosome.begin <- chromopos(chromosome, 0)
  pos <- chromopos(chromosome, position)
  stopifnot(!is.unsorted(pos))
  mids <- 0.5 * (c(-Inf, pos) + c(pos, Inf))
  after <- pmin(mids[-1], chromosome.end)
  before <- pmax(mids[-length(mids)], chromosome.begin)
  return(after - before)
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

