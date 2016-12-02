nucleotide.levels <- c('C','G','T','A')
genotype.levels <- c('CC','GG','TT','AA','GC','TG','AT','TC','AG','AC')

as.nucleotide <- function(x) {
  factor(x, nucleotide.levels)
}

is.nucleotide <- function(x) {
  identical(levels(x), nucleotide.levels)
}

transition.nucleotide <- function(nucleotide) {
  as.nucleotide(c('T','A','C','G')[nucleotide])
}

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

          #        'CC','GG','TT','AA','GC','TG','AT','TC','AG','AC'
genotype.mask <- c(0x1, 0x2, 0x4, 0x8, 0x3, 0x6, 0xC, 0x5, 0xA, 0x9)

nucleotide.in.genotype <- function(nucleotide, genotype) {
  stopifnot(identical(levels(nucleotide), nucleotide.levels))
  stopifnot(identical(levels(genotype), genotype.levels))
  n <- 2^(as.integer(nucleotide) - 1)
  return(as.logical(bitwAnd(n, genotype.mask[as.integer(genotype)])))
}

genotypes.overlap <- function(genotype1, genotype2) {
  stopifnot(identical(levels(genotype1), genotype.levels))
  stopifnot(identical(levels(genotype2), genotype.levels))
  mask1 <- genotype.mask[as.integer(genotype1)]
  mask2 <- genotype.mask[as.integer(genotype2)]
  return(as.logical(bitwAnd(mask1, mask2)))
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
  df$allele1 <- NULL
  df$allele2 <- NULL
  return(df)
}

read.23andme.raw <- function(file) {
  cols <- c('character', 'character', 'integer', 'character')
  df <- read.table(file, colClasses=cols, as.is=TRUE)
  names(df) <- c("rsid", "chromosome", "position", "genotype")
  return(df) 
}

read.23andme <- function(file) {
  df <- read.23andme.raw(file)
  del <- with(df, chromosome == 'MT' | substr(rsid,1,2) != "rs")
  df <- df[!del,]
  df$chromosome[df$chromosome == 'X'] <- "23"
  df$chromosome[df$chromosome == 'Y'] <- "24"
  df$chromosome <- as.integer(df$chromosome)
  haplo <- df$genotype %in% nucleotide.levels
  df$chromosome[df$chromosome == 23 & !haplo] <- 25
  allele1 <- as.nucleotide(substr(df$genotype, 1, 1))
  allele2 <- as.nucleotide(substr(df$genotype, 2, 2))
  allele2[haplo] <- allele1[haplo]
  df$genotype <- as.genotype(allele1, allele2)
  return(df)
}

read.familyfinder <- function(file) {
  df <- read.csv(file, colClasses="character")
  expect <- c("RSID", "CHROMOSOME", "POSITION", "RESULT")
  stopifnot(all.equal(names(df), expect))
  names(df) <- c("rsid", "chromosome", "position", "genotype")
  del <- with(df, substr(rsid,1,2) != "rs")
  df <- df[!del,]
  df$chromosome[df$chromosome == 'X'] <- "23"
  df$chromosome[df$chromosome == 'Y'] <- "24"
  df$chromosome <- as.integer(df$chromosome)
  df$position <- as.integer(df$position)
  allele1 <- as.nucleotide(substr(df$genotype, 1, 1))
  allele2 <- as.nucleotide(substr(df$genotype, 2, 2))
  df$genotype <- as.genotype(allele1, allele2)
  return(df) 
}

read.dna.text <- function(file, filename) {
  len <- nchar(filename)
  ancestrydna <- (substr(filename, len - 14, len) == "AncestryDNA.txt")
  if (ancestrydna) return(read.ancestrydna(file))
  extension <- substr(filename, nchar(filename)-3, nchar(filename))
  if (extension == ".csv") return(read.familyfinder(file))
  return(read.23andme(file))
}

read.dna.zip <- function(file) {
  zip.list <- unzip(file, list=TRUE)
  if (NROW(zip.list) > 1) warning("More than one file inside zip file ", file)
  filename <- zip.list$Name[1]
  return(read.dna.text(unz(file, filename), filename))
}

read.dna <- function(file) {
  parts <- strsplit(file, '/', fixed=TRUE)[[1]]
  filename <- parts[length(parts)]
  parts <- strsplit(filename, '.', fixed=TRUE)[[1]]
  extension <- parts[length(parts)]
  if (extension %in% c("txt", "csv")) {
    return(read.dna.text(file, filename))
  }
  if (extension == "gz") {
    filename <- substr(filename, 1, nchar(filename) - 3)
    return(read.dna.text(gzfile(file), filename))
  }
  return(read.dna.zip(file))
}

read.dna.web <- function(url) {
  tmpfilename <- tempfile()
  on.exit(unlink(tmpfilename))
  download.file(url, tmpfilename)
  read.dna.zip(tmpfilename)
}

