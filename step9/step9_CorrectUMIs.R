
#!/usr/bin/env Rscript
# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
# R SCRIPT
# UMI CORRECTION
version <- '4.12.0'
library(tidyverse)
source('./VariantCalling_functions.R')

# Arguments: Mutations: tibble, PointMutations: tibble, ReadMutations: tibble
SingleCellMutations <- read_csv()
Mutations <- read_csv(args[1])
PointMutations <- read_tsv(args[2], col_types = "cdccccdcc", col_names = FALSE)
ReadMutations <- read_tsv(args[3], col_names = FALSE)

# test if there is at least one argument: if not, return an error; 
# check if the arguments are of the proper type
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (is_tibble(SingleCellMutations) == F) {
  stop("wrong object type: Make sure mutations are in a tibble")
}

# score mutations for umi support:
filtered.mutations <- score.mutations(mutations, point, read)

# Write scored mutations to drive:
write(filtered.mutations, file = 'ScoredMutations')
