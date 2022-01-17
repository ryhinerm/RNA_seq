#!/bin/bash

#SBATCH --job-name=QC
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Load modue
module add UHTS/Quality_control/fastqc/0.11.9

# Create and set directory
mkdir ../ex4_QC
cd ../ex4_QC

# Perform quality check
fastqc -o . -f fastq ../ex3_clip_trim/*_tr.fastq.gz
