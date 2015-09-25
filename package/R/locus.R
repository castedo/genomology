
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

