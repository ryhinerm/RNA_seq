#!/bin/bash

#SBATCH --job-name=multiQC
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Create and set directory
mkdir ../exA_multiQC
cd ../exA_multiQC

# Copy files
cp ../ex1_QC/*.zip .
cp ../ex4_QC/*.zip .
cp ../ex4_map_non-mRNA/*_log.txt .
cp ../ex4_map_genome/*_log.txt .

# Load module
module add UHTS/Analysis/MultiQC/1.8

# Perform quality check
multiqc . 
