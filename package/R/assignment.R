# genotype allele assignment signs
# 0 = assignment ambiguous or genotype is homozygous
# 1 = (allele1, allele2) must be assigned to (person1, person2) (no flip)
# -1 = (allele1, allele2) must be assigned to (person2, person1) (flip)
# NA = assignment impossible

allele.assigned <- function(assignment, genotype) {
  stopifnot(all(assignment %in% c(-1:1,NA)))
  stopifnot(identical(levels(genotype), genotype.levels))
  ret <- allele1(genotype)
  want.allele2 <- (assignment == -1) & !is.na(assignment)
  ret[want.allele2] <- allele2(genotype)[want.allele2]
  ambiguous <- heterozygous(genotype) & !assignment
  ret[ambiguous | is.na(assignment)] <- NA
  return(ret)
}

assign.parent <- function(parent1, kid) {
  stopifnot(identical(levels(parent1), genotype.levels))
  stopifnot(identical(levels(kid), genotype.levels))
  not.a1 <- !nucleotide.in.genotype(allele1(kid), parent1)
  not.a2 <- !nucleotide.in.genotype(allele2(kid), parent1)
  ret <- (not.a2 - not.a1)
  ret[not.a1 & not.a2] <- NA # not possible
  return(ret)
}

reconcile.signs <- function(sign1, sign2) {
  stopifnot(all(sign1 %in% c(-1:1,NA)))
  stopifnot(all(sign2 %in% c(-1:1,NA)))
  ret <- sign(sign1 + sign2)
  ret[sign1 * sign2 < 0] <- NA # conflicting signs
  return(ret)
}

assign.parents <- function(mom, pop, kid, from.mom.only=FALSE) {
  stopifnot(identical(levels(kid), genotype.levels))
  stopifnot(identical(levels(mom), genotype.levels))
  stopifnot(identical(levels(pop), genotype.levels))
  stopifnot(!any(from.mom.only & heterozygous(kid)))
  mom.sign <- assign.parent(mom, kid)
  ret <- reconcile.signs(mom.sign, -assign.parent(pop, kid))
  ret[from.mom.only] <- mom.sign[from.mom.only]
  return(ret)
}

# assign positive to parents's mom and negative to parent's pop
assign.grandparents <- function(parent.sign, passage.sign) {
  stopifnot(all(parent.sign %in% c(-1:1,NA)))
  stopifnot(all(passage.sign %in% c(-1:1,NA)))
  return(parent.sign * passage.sign)
}

assign.passage <- function(parent, kid) {
  stopifnot(identical(levels(parent), genotype.levels))
  stopifnot(identical(levels(kid), nucleotide.levels))
  not.a1 <- kid != allele1(parent)
  not.a2 <- kid != allele2(parent)
  ret <- (not.a2 - not.a1)
  ret[not.a1 & not.a2] <- NA # not possible
  return(ret)
}

