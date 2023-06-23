#!/bin/bash
#SBATCH --time=0-1  # days-hours
#SBATCH --job-name=vdj  # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of cores
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --reservation=rnaworkshop
#SBATCH --account=workshop
#SBATCH --array=1
#SBATCH --output=slurmout/cellranger-vdj-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=cslurmout/cellranger-vdj-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=user@domain.com # Email to which notifications will be sent

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

crBase=/share/workshop/vdj_workshop/Software
crBin=${crBase}/cellranger-7.1.0/bin
crRef=${crBase}/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
## Where cellranger and bcl2fastq executables are located
## a) by loading a module
#module load cellranger/7.0.0

## b) or, by placing the location of the executables on the path (edit to your location)
export PATH=${crBin}:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basepath=/share/workshop/intro_scrnaseq/$USER/scrnaseq_example
fastqs=${basepath}/00-RawData
outdir=${basepath}/01-Cellranger

## provide the script the row # of the sample to be run
sample=`sed "$1q;d" ${basepath}/samples.txt`

## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
## Create the call
call="cellranger vdj \
    --id=$sample \
    --fastqs=$fastqs \
    --sample=$sample \
    --reference=$crRef"

## Echo the call
echo $call
## Evaluate the call
eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
