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

read.ancestrydna <- function(first.line, file) {
  if (!grepl("AncestryDNA", first.line, ignore.case=TRUE)) return(NULL)
  cols <- c('character', 'integer', 'integer', 'character', 'character')
  df <- read.table(file, header=TRUE, colClasses=cols, as.is=TRUE)
  return(read.ancestrydna.data(df))
}

read.ancestrydna.data <- function(df) {
  df$allele1 <- factor(df$allele1, nucleotide.levels)
  df$allele2 <- factor(df$allele2, nucleotide.levels)
  df$genotype <- as.genotype(df$allele1, df$allele2)
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
  allele1 <- as.nucleotide(substr(df$genotype, 1, 1))
  allele2 <- as.nucleotide(substr(df$genotype, 2, 2))
  allele2[haplo] <- allele1[haplo]
  df$genotype <- as.genotype(allele1, allele2)
  return(df)
}

read.myheritage <- function(first.line, file) {
  if (!grepl("MyHeritage", first.line, ignore.case=TRUE)) return(NULL)
  df <- read.csv(file, header=TRUE, colClasses="character", comment.char="#")
  return(read.myheritage.data(df))
}

read.familyfinder <- function(first.line, file) {
  if (grepl("famfinder", first.line, ignore.case=TRUE)) {
    col.names <- c("rsid", "chromosome", "position", "allele1", "allele2")
    cls <- c('character', 'character', 'integer', 'character', 'character')
    df <- read.csv(file, col.names=col.names, colClasses=cls, comment.char="#")
    df$chromosome[df$chromosome == 'X'] <- "23"
    df$chromosome[df$chromosome == 'Y'] <- "24"
    df$chromosome[df$chromosome == 'XY'] <- "25"
    df$chromosome <- as.integer(df$chromosome)
    ok <- with(df, substr(rsid,1,2) == "rs" & chromosome > 0)
    df <- df[ok,]
    # new 2019 family finder data rows cleaned same as AncestryDNA
    return(read.ancestrydna.data(df))
  }
  col.names <- c("RSID", "CHROMOSOME", "POSITION", "RESULT")
  if (first.line != paste(col.names, collapse=",")) return(NULL)
  df <- read.csv(file, col.names=col.names, colClasses="character")
  # old pre-2019 family finder data rows same as MyHeritage format
  return(read.myheritage.data(df))
}

read.myheritage.data <- function(df) {
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

read.geno2 <- function(first.line, file) {
  # first.line sometimes header, sometimes non-rsid data, but maybe rsid data
  # regardless it gets ignored here
  df <- read.csv(file, colClasses="character", comment.char="#")
  return(read.geno2.data(df))
}

read.geno2.data <- function(df) {
  names(df) <- c("rsid", "chromosome", "allele1", "allele2")
  del <- with(df, substr(rsid,1,2) != "rs")
  df <- df[!del,]
  df$chromosome[df$chromosome == 'X'] <- "23"
  df$chromosome[df$chromosome == 'Y'] <- "24"
  df$chromosome <- as.integer(df$chromosome)
  df$position <- NA
  df$allele1 <- factor(df$allele1, nucleotide.levels)
  df$allele2 <- factor(df$allele2, nucleotide.levels)
  df$genotype <- as.genotype(df$allele1, df$allele2)
  df$allele1 <- NULL
  df$allele2 <- NULL
  return(df)
}

read.dna.text <- function(file, filename=NULL) {
  if (is.character(file)) {
    file <- file(file, "rt")
    on.exit(close(file))
  }
  if (!isOpen(file, "rt")) {
    open(file, "rt")
    on.exit(close(file))
  }
  line <- readLines(file, n=1)
  ret <- read.familyfinder(line, file)
  if (!is.null(ret)) return(ret);
  ret <- read.myheritage(line, file)
  if (!is.null(ret)) return(ret);
  ret <- read.ancestrydna(line, file)
  if (!is.null(ret)) return(ret);
  if (!is.null(filename) && endsWith(filename, ".all.csv")) {
    return(read.geno2(line, file))
  }
  return(read.23andme(file))
}

read.dna.zip <- function(file) {
  zip.list <- unzip(file, list=TRUE)
  filenames = zip.list$Name
  if (length(filenames) > 1) {
    # check for National Geographic Geno 2.0 genetic.zip file
    geno2 <- endsWith(filenames, ".all.csv")
    if (any(geno2)) {
      filenames <- filenames[geno2]
    }
    if (length(filenames) > 1) warning("More than one file choice inside zip file ", file)
  }
  filename <- filenames[1]
  return(read.dna.text(unz(file, filename), filename))
}

read.dna <- function(file) {
  parts <- strsplit(file, '/', fixed=TRUE)[[1]]
  filename <- parts[length(parts)]
  parts <- strsplit(filename, '.', fixed=TRUE)[[1]]
  extension <- parts[length(parts)]
  if (tolower(extension) %in% c("txt", "csv")) {
    return(read.dna.text(file, filename))
  }
  if (extension == "gz") {
    return(read.dna.text(gzfile(file)))
  }
  return(read.dna.zip(file))
}

read.dna.web <- function(url) {
  tmpfilename <- tempfile()
  on.exit(unlink(tmpfilename))
  download.file(url, tmpfilename)
  read.dna.zip(tmpfilename)
}

