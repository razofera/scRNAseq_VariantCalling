
#!/usr/bin/env Rscript
# parameter arguments:
args = commandArgs(trailingOnly=TRUE)
# R SCRIPT
# UMI CORRECTION
version <- '4.12.0'
library(tidyverse)

# Arguments: metadata: tibble, mutations: tibble
#sc_metadata <- read_csv(args[1])
sc_mutations <- read_csv(args[1])

# test if there is at least one argument: if not, return an error; 
# check if the arguments are of the proper type
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (is_tibble(sc_mutations) == F) {
  # error if incorrect type of input
  stop("wrong object type: Make sure mutations are in a tibble")
}

# Filter passing mutations:

filter.mutations <- all.mutations %>%
  filter(TLOD >= 5.3 & DP > 10 & ECNT > 1)

# Write arguments for next script:
args <- as.vector(rbind(filter.mutations$Chr,filter.mutations$POS,filter.mutations$bc,filter.mutations$POS))
dir.create('/TL')
write(args, file = '/TL/TL_cell')
dir.create('/umi_correction_out')
