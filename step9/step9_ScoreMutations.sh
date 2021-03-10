#!/bin/bash

echo "mutations_list: $1";
echo "mutations_Reads: $2";
echo "mutations_Metadata: $3";

# GENERATE MUTATION SCORES:

Rscript step9_ScoreMutations.R ${1} ${2} ${3}


