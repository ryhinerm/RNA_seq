#!/bin/bash

#SBATCH --job-name=faTo2bit
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Create and set directory
mkdir ../ex6_2bit
cd ../ex6_2bit

# Load module
module add SequenceAnalysis/blat/36

# Convert fasta to 2bit
faToTwoBit ../ex2_Rnor_seq/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus.Rnor_6.0.dna.toplevel.2bit
