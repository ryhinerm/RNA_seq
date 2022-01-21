#!/bin/bash

#SBATCH --job-name=biotype_count
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Create and set directory
mkdir ../ex7_biotype_counts
cd ../ex7_biotype_counts

# Load module
module add UHTS/Analysis/subread/2.0.1

# Generate gene counts
featureCounts -T 4 -t exon -g gene_biotype -a ../ex2_Rnor_seq/Rattus_norvegicus.Rnor_6.0.104.gtf -o biotype_counts_rawfile.txt ../ex4_map_genome/*_sorted.bam

# Extract sample columns
cut -f 1,7-12 biotype_counts_rawfile.txt > biotype_counts_processed.txt
