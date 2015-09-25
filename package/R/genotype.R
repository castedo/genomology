nucleotide.levels <- c('C','G','T','A')
genotype.levels <- c('CC','GG','TT','AA','GC','TG','AT','TC','AG','AC')

as.genotype <- function(allele1, allele2) {
  stopifnot(identical(levels(allele1), nucleotide.levels))
  stopifnot(identical(levels(allele2), nucleotide.levels))
  i1 <- as.integer(allele1)
  i2 <- as.integer(allele2)
  ## this is just a convenient arithemetic mapping to the range 1 to 10
  d <- abs(i1 - i2)
  nn <- pmin(i1, i2) + (d * (9 - d)) %/% 2
  factor(nn, 1:10, genotype.levels)
}

homozygous <- function(genotype) { as.integer(genotype) <= 4 }
heterozygous <- function(genotype) { as.integer(genotype) > 4 }

allele1 <- function(genotype) {
  n1 <- c(nucleotide.levels, 'G', 'T', 'A', 'T', 'A', 'A')[genotype]
  factor(n1, nucleotide.levels)
}

allele2 <- function(genotype) {
  n2 <- c(nucleotide.levels, 'C', 'G', 'T', 'C', 'G', 'C')[genotype]
  factor(n2, nucleotide.levels)
}

nucleotide.in.genotype <- function(nucleotide, genotype) {
  stopifnot(identical(levels(nucleotide), nucleotide.levels))
  stopifnot(identical(levels(genotype), genotype.levels))
  n <- 2^(as.integer(nucleotide) - 1)
  #        'CC','GG','TT','AA','GC','TG','AT','TC','AG','AC'
  mask <- c(0x1, 0x2, 0x4, 0x8, 0x3, 0x6, 0xC, 0x5, 0xA, 0x9)
  return(as.logical(bitwAnd(n, mask[as.integer(genotype)])))
}

read.ancestrydna.raw <- function(file) {
  cols <- c('character', 'integer', 'integer', 'character', 'character')
  df <- read.table(file, header=TRUE, colClasses=cols, as.is=TRUE)
  df$allele1 <- factor(df$allele1, nucleotide.levels)
  df$allele2 <- factor(df$allele2, nucleotide.levels)
  return(df)
}

read.ancestrydna <- function(file) {
  df <- read.ancestrydna.raw(file)
  df$genotype <- as.genotype(df$allele1, df$allele2)
  # remove rare highly suspect data (often bogus)
  bad <- (as.integer(df$allele1) < as.integer(df$allele2) |
          df$genotype %in% c('GC', 'AT'))
  df$genotype[bad] <- NA
  return(df)
}

read.ancestrydna.zip <- function(file) {
  read.ancestrydna(unz(file, "AncestryDNA.txt"))
}

read.ancestrydna.web <- function(url) {
  tmpfilename <- tempfile()
  on.exit(unlink(tmpfilename))
  download.file(url, tmpfilename)
  read.ancestrydna(unz(tmpfilename, "AncestryDNA.txt"))
}

read.23andme.raw <- function(file) {
  cols <- c('character', 'character', 'integer', 'character')
  df <- read.table(file, colClasses=cols, as.is=TRUE)
  names(df) <- c("rsid", "chromosome", "position", "genotype")
  return(df) 
}

read.23andme <- function(file) {
  df <- read.23andme.raw(file)
  del <- with(df, chromosome == 'MT'
                  | substr(rsid,1,2) != "rs"
                  | genotype %in% c("--", "D", "DD", "DI", "I", "II"))
  df <- df[!del,]
  df$chromosome[df$chromosome == 'X'] <- "23"
  df$chromosome[df$chromosome == 'Y'] <- "24"
  df$chromosome <- as.integer(df$chromosome)
  haplo <- df$genotype %in% nucleotide.levels
  df$chromosome[df$chromosome == 23 & !haplo] <- 25
  allele1 <- factor(substr(df$genotype, 1, 1), nucleotide.levels)
  allele2 <- factor(substr(df$genotype, 2, 2), nucleotide.levels)
  allele2[haplo] <- allele1[haplo]
  df$genotype <- as.genotype(allele1, allele2)
  return(df)
}

read.23andme.zip <- function(file) {
  zip.list <- unzip(file, list=TRUE)
  if (NROW(zip.list) > 1) warning("More than one file inside file zip ", file)
  read.23andme(unz(file, zip.list$Name[1]))
}

read.23andme.web <- function(url) {
  tmpfilename <- tempfile()
  on.exit(unlink(tmpfilename))
  download.file(url, tmpfilename)
  read.23andme.zip(tmpfilename)
}

