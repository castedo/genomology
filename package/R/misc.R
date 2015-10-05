
read.omim.alleles <- function(gene.id) {
  url <- paste0("http://www.omim.org/allelicVariant/", gene.id, "?format=tab")
  ret <- read.delim2(url, colClasses="character", skip=8)
  return(ret)
}

