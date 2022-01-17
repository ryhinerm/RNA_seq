#!/bin/bash

#SBATCH --job-name=multiQC
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Set directory
cd ../exA_multiQC

# Load module
module add UHTS/Analysis/MultiQC/1.8

# Perform quality check
multiqc . 
