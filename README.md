# Genomic genealogy library in R

R code for analysis of family genome-wide DNA data (e.g. AncestryDNA and
23andme raw DNA data).

Visit the [R homepage](https://www.r-project.org) to learn more about the R
statistical programming language. R is a popular system for doing professional
DNA data analysis with [Bioconductor](http://www.bioconductor.org) packages.

This R package, genomology, is for doing genetic geneology and is a bit
"hobby-ish" compared to what bioinformatics professionals would typically do
with R and Bioconductor.

An example usage of this packge could be determining what percentage of genes a
grand-child has received from grand-parents. This can be done with raw DNA data
from AncestryDNA or 23andme for a grandparent, a "middle" parent and a
grandchild.

Help
====

You can get package and function documentation from within R by doing:

```
library(genomology)
help(package="genomology")
```

Example
=======

In this example:
* three family members DNA data is read (grandma, dad, girl) from different DNA
  services
* DNA data is merged together with extra information like genes and basepair
  ranges
* genotypes for lactase persistance (lactose tolerance) regulating gene MCM6 is
  printed
* estimation of which chromosome regions of girl are descended from grandma
* bar plot by chromosome of how much of girls DNA is inherited from dad's mom
  (grandma) vs dad's dad

```
library(genomology)
set.genomology.loci(read.csv("http://ref.castedo.com/dna/loci.csv", comment.char="#"))

grandma <- read.ancestrydna.web("http://example.com/dna-data-grandma.zip")
dad <- read.23andme.web("http://example.com/dna-data-dad.zip")
girl <- read.ancestrydna.web("http://example.com/dna-data-kid.zip")

dna <- grandma[, 1:3]
dna$zone <- chromozone(dna$chromosome, dna$position)
dna$gene <- locus.at(dna$chromosome, dna$position)
dna$gma <- grandma$genotype
dna$dad <- dad$genotype[match(dna$rsid, dad$rsid)]
dna$girl <- girl$genotype

dna <- dna[chromorder23(dna$chromosome, dna$position),]
dna$range <- chromorange(dna$chromosome, dna$position)

subset(dna, gene == 'LCT')

attach(dna)

kidNuc <- allele.assigned(assign.parent(dad, girl), girl)
dadMa <- assign.parent(gma, dad)
kidGrandSign <- assign.grandparents(dadMa, assign.passage(dad, kidNuc))
kidGrandSign[chromosome == 23] <- 1 # dad passes grandma's X to girl
gma.prob <- assign.descent(chromordinate(chromosome, position), kidGrandSign)

tab <- tapply(range/sum(range), list(call=prob.triage(gma.prob), zone=zone), sum)
tab[is.na(tab)] <- 0

barplot(tab, horiz=TRUE, las=1, ylab="Chromosome")
```

MATH
====

The core mathematical functions used in this package are `idescent.prob` and `recombination.prob`.
Documentation is in the package help, but only the PDF (LaTeX) versions have
the mathematical derivations. For convenience, the PDF versions are also online:

* [idescent.prob](http://ref.castedo.com/dna/idescent.prob.pdf)
* [recombination.prob](http://ref.castedo.com/dna/recombination.prob.pdf)

INSTALL
=======

Clone this repository and from inside do

```
R CMD INSTALL package
```

Please email castedo@castedo.com if you install this package or want a built R
package online (which would be easier to install).

