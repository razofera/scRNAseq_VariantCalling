#!/bin/bash

echo "cell: $1";
echo "metadata: $2";
echo "mutations: $3";

# INITIALIZE VARIABLES/ FUNCTIONS / DIRECTORY PATHS:

SAM2TSV=/dist

function extract_meta {
 	## use awk to find fields that match patterns
 	awk \'{ for (i=1; i<=NF; ++i) { if ($i ~ /J00/) { for (j=1; j<=NF; ++j) { if ($j ~ /CB/) { for (k=1; k<=NF; ++k) { if ($k ~ /UB/) {print $i"___"$j"___"$k} } } } } } }\'
}; export -f extract_meta

# GENERATE METADATA:

Rscript UMI_CORRECTION_4.12.0.R ${2} ${3}

# EXTRACT METADATA FOR CELL MUTATIONS:

samtools index ${1}

cat /TL/TL_cell \
| parallel --jobs=30 --max-args=4 samtools view -b -S -h ${1} {1}:{2}-{2} \
| java -jar ${SAM2TSV}/sam2tsv.jar \
| grep -w {4} \
>> reads.tsv

cat /TL/TL_cell \
| parallel --jobs=30 --max-args=4 samtools view ${1} {1}:{2}-{2} \\
| extract_meta \
>> meta.tsv

ls -l

