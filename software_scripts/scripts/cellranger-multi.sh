#!/bin/bash

## Where cellranger executable is located
export PATH=/share/workshop/vdj_workshop/Software/cellranger-7.1.0/bin:$PATH

## Set the parameters for the run
basedir="/share/workshop/vdj_workshop/${USER}"
reference="/share/workshop/vdj_workshop/Software/refdata-gex-GRCh38-2020-A"
fastqs="${basedir}/00-RawData"
outdir="${basedir}/01-Cellranger"
configs="${fastqs}/multiconfig"

mkdir -p $outdir
cd $outdir

## loop over samples in sample sheet, running cellranger on each
## NOTE: this script expects a sample sheet with one line per sample to be present in the basedir
## NOTE: this script expects a config file for each sample to be located in fastqs/multiconfig
for sample in `cat ${basedir}/samples.txt`
do
  ## Create the call
  call="cellranger multi \
    --id=${sample} \
    --csv=${configs}/${sample}_multi.csv"

## Echo the call
  echo $call
  ## Evaluate the call
  eval $call
done

