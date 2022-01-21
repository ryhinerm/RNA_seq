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
mkdir ../ex1_QC
cd ../ex1_QC

# Perform quality check
fastqc -o . -f fastq ../ex1_download/*.fastq.gz
