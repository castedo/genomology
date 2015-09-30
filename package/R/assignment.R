assignmentClass <- setClass("assignment", contains="integer")

# genotype allele assignment codes
# 0 = assignment ambiguous or genotype is homozygous
# 1 = (allele1, allele2) must be assigned to (person2, person1) (flip)
# 2 = (allele1, allele2) must be assigned to (person1, person2) (do not flip)
# 3 = assignment impossible

as.sign <- function(assignment) {
  stopifnot(is(assignment, assignmentClass))
  c(0, -1, 1, NA)[assignment + 1]
}

is.impossible.assignment <- function(assignment) {
  stopifnot(is(assignment, assignmentClass))
  return(assignment == 3)
}

allele.assigned <- function(person, assignment, genotype) {
  stopifnot(isTRUE(person %in% 1:2))
  stopifnot(is(assignment, assignmentClass))
  stopifnot(identical(levels(genotype), genotype.levels))
  ret <- allele1(genotype)
  want.allele2 <- (assignment == person) & !is.na(assignment)
  ret[want.allele2] <- allele2(genotype)[want.allele2]
  ambiguous <- heterozygous(genotype) & !assignment
  stomp <- ambiguous | is.impossible.assignment(assignment)
  ret[stomp] <- NA
  ret[is.na(assignment)] <- NA
  return(ret)
}

assign.parent <- function(parent1, kid) {
  stopifnot(identical(levels(parent1), genotype.levels))
  stopifnot(identical(levels(kid), genotype.levels))
  bit0 <- !nucleotide.in.genotype(allele1(kid), parent1)
  bit1 <- !nucleotide.in.genotype(allele2(kid), parent1)
  return(as(bit0 + bit1 * 2, "assignment"))
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
  mom.sign <- as.sign(assign.parent(mom, kid))
  ret <- reconcile.signs(mom.sign, -as.sign(assign.parent(pop, kid)))
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

