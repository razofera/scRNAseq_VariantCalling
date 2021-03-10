#!/bin/bash

# ARGUMENTS
Sample=${Sample:default_sample_num}
DataDirectory=${DataDirectory:default_data_dir}
Targets=${Targets:default_targets}

READ_EXOME=$(samtools view -c -F 4 -L ${Targets} ${DataDirectory}/${Sample})
READ_SAM=$(samtools view -c -F 4 ${DataDirectory}/${Sample})
UMI=$(samtools view  ${DataDirectory}/${Sample} | grep -o \'UB:............\' | grep -o \'..........$\' | uniq -c | wc -l)
COV_EXOME=$(samtools depth -b ${Targets} ${DataDirectory}/${Sample} | wc -l)
COV=$(samtools depth ${DataDirectory}/${Sample} | wc -l)
echo $Sample, $READ_EXOME, $READ_SAM, $UMI, $COV_EXOME, $COV >> mutations.csv
