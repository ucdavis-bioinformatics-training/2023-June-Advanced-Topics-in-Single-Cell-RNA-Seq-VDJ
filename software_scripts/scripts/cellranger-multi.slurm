#!/bin/bash
#SBATCH --time=1-00 # days-hours
#SBATCH --job-name=vdj  # Job name
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --reservation=rnaworkshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/cellranger-vdj-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/cellranger-vdj-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=user@email.com # Email to which notifications will be sent

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

## Where cellranger executable is located
## a) by loading a module
#module load cellranger/7.0.0

## b) or, by placing the location of the executables on the path (edit to your location)
export PATH=/share/workshop/vdj_workshop/Software/cellranger-7.1.0/bin:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basedir="/share/workshop/vdj_workshop/$USER"
reference=/share/workshop/vdj_workshop/Software/refdata-gex-GRCh38-2020-A
fastqs="${basedir}/00-RawData"
outdir="${basedir}/01-Cellranger"
configs="${fastqs}/multiconfig"

## provide the script the row # of the sample to be run
## NOTE: this script expects a file called samples.txt with one sample name per line in the basedir
## NOTE: this script expects a config file for each sample located in fastqs/multiconfig
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${basedir}/samples.txt`

## Create the call
call="cellranger multi \
  --id=${sample} \
  --csv=${configs}/${sample}_multi.csv"

## Create output directory and change directory
mkdir -p $outdir
cd $outdir
## Echo the call
echo $call
## Evaluate the call
eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime

