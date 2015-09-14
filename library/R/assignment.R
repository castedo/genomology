assignmentClass <- setClass("assignment", contains="integer")

# genotype assignment codes
# 0 = assignment ambiguous or genotype is homozygous
# 1 = (allele1, allele2) must be assigned to (person2, person1) (flip)
# 2 = (allele1, allele2) must be assigned to (person1, person2) (do not flip)
# 3 = assignment impossible

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

flip.assignment <- function(assignment) {
  stopifnot(is(assignment, assignmentClass))
  bit0 <- assignment %/% 2
  bit1 <- assignment %% 2
  return(as(bit0 + bit1 * 2, "assignment"))
}

reconcile.assignments <- function(assignment1, assignment2) {
  stopifnot(is(assignment1, assignmentClass))
  stopifnot(is(assignment2, assignmentClass))
  return(as(bitwOr(assignment1, assignment2), "assignment"))
}

assign.parents <- function(mom, pop, kid, from.mom.only=FALSE) {
  stopifnot(identical(levels(kid), genotype.levels))
  stopifnot(identical(levels(mom), genotype.levels))
  stopifnot(identical(levels(pop), genotype.levels))
  stopifnot(!any(from.mom.only & heterozygous(kid)))
  momass <- assign.parent(mom, kid)
  ass <- reconcile.assignments(momass,
                               flip.assignment(assign.parent(pop, kid)))
  ass[from.mom.only] <- momass[from.mom.only]
  return(ass)
}

# assign parents's mom as person1 and parent's pop as person2
assign.grandparents <- function(parental.assignment, passage.assignment) {
  stopifnot(is(parental.assignment, assignmentClass))
  stopifnot(is(passage.assignment, assignmentClass))
  code <- parental.assignment
  reverse <- (passage.assignment == 1) # allele2 was passed down
  code[reverse] <- flip.assignment(code)[reverse] # so reserve assignment
  code[passage.assignment == 0] <- 0 # ambiguous/homozygous allele passed down
  code[passage.assignment == 3] <- 3 # not possible
  return(as(code, "assignment"))
}

assign.passage <- function(parent, kid) {
  stopifnot(identical(levels(parent), genotype.levels))
  stopifnot(identical(levels(kid), nucleotide.levels))
  bit0 <- kid != allele1(parent)
  bit1 <- kid != allele2(parent)
  code <- bit0 + bit1 * 2
  code[is.na(kid)] <- 0
  return(as(code, "assignment"))
}

