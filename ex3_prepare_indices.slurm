#!/bin/bash

#SBATCH --job-name=prepare_indices
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end

# Load module
module add UHTS/Aligner/bowtie/1.2.0

# Create and set directory
mkdir ../ex3_prepare_indices
cd ../ex3_prepare_indices

# Create and set subdirectories, prepare indices
mkdir Rnor_r_sn_sno_t_RNA
cd Rnor_r_sn_sno_t_RNA
bowtie-build ../ex2_Rnor_seq/Rnor_r_sn_sno_t_RNA.fa Rnor_r_sn_sno_t_RNA

mkdir ../Rattus_norvegicus.Rnor_6.0.dna.toplevel
cd ../Rattus_norvegicus.Rnor_6.0.dna.toplevel
bowtie-build ../ex2_Rnor_seq/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus.Rnor_6.0.dna.toplevel
 
