#!/bin/bash

#SBATCH --job-name=concatenate
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Set directory
cd ../ex2_Rnor_seq

# Concatenate sequences
cat *.txt > Rnor_r_sn_sno_t_RNA.fa
