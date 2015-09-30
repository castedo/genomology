# genomology
Genomic genealogy library in R

R code for analysis of family genome-wide DNA data (e.g. AncestryDNA raw DNA data).

Example
=======

See SNP sequence for lactase persistance (lactose tolerance) gene LCT:

```
library(genomology)
set.genomology.loci(read.csv("http://ref.castedo.com/dna/loci.csv", comment.char="#"))
mydna <- read.ancestrydna.web("http://example.com/dna-data-2010-01-01.zip")
mydna$gene <- locus.at(mydna$chromosome, mydna$position)
subset(mydna, gene == 'LCT')
```

INSTALL
=======

Clone this repository and from inside do

```
R CMD INSTALL package
```

Please email castedo@castedo.com if you install this package or want a built R
package online (which would be easier to install).

