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

read.rawancestrydna <- function(file) {
  cols <- c('character', 'integer', 'integer', 'character', 'character')
  df <- read.table(file, header=TRUE, colClasses=cols, as.is=TRUE)
  df$allele1 <- factor(df$allele1, nucleotide.levels)
  df$allele2 <- factor(df$allele2, nucleotide.levels)
  return(df)
}

read.ancestrydna <- function(file) {
  df <- read.rawancestrydna(file)
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
  download.file(url, tmpfilename)
  read.ancestrydna(unz(tmpfilename, "AncestryDNA.txt"))
  unlink(tmpfilename)
}

