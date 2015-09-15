# genomology
Genomic genealogy library in R

Synopsis
========

See SNP sequence for lactase persistance (lactose tollerance) gene LCT:

```
library(genomology)
set.genomology.loci(read.csv("http://ref.castedo.com/dna/loci.csv", comment.char="#"))
mydna <- read.ancestrydna.web("http://example.com/dna-data-2010-01-01.zip")
mydna$gene <- locus.at(mydna$chromosome, mydna$position)
subset(mydna, gene == 'LCT')
```
