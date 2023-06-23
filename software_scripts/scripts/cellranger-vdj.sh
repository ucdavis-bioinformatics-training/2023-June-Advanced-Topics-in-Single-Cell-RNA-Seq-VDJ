#!/bin/bash

# edit this section to reflect your cellranger installation and reference location
crBase=/share/workshop/vdj_workshop/Software/
crBin=${crBase}/cellranger-7.1.0/bin
crRef=${crBase}/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

# edit this section to reflect the location and organization of your project directory
basepath=/share/workshop/intro_scrnaseq/$USER/scrnaseq_example
fastqs=${basepath}/00-RawData
outdir=${basepath}/01-Cellranger
sample=`sed "$1q;d" ${basepath}/samples.txt`

export PATH=${crBin}:$PATH
# run cellranger
cellranger vdj \
    --id=$sample \
    --fastqs=$fastqs \
    --sample=$sample \
    --reference=$crRef
