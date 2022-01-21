#!/bin/bash

#SBATCH --job-name=map_to_genome
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=02:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Create and set directory
mkdir ../ex4_map_genome
cd ../ex4_map_genome

# Load modules
module add UHTS/Aligner/bowtie/1.2.0
module add UHTS/Analysis/samtools/1.10

# Map to genome, convert to BAM
for x in $(ls -d ../ex4_map_non-mRNA/*_RNA.fastq); \
do echo ${x}; \
bowtie \
-S \
-t \
-p 4 \
-v 1 \
-m 1 \
--best \
--strata \
../ex3_prepare_indices/Rattus_norvegicus.Rnor_6.0.dna.toplevel/Rattus_norvegicus.Rnor_6.0.dna.toplevel \
-q ${x} \
2> $(basename ${x} .fastq)_Rnor_6_log.txt | \
samtools view -h -F 4 -b > $(basename ${x} .fastq)_Rnor_6.bam; done

# Sort BAM
for x in $(ls -d *.bam); \
do echo ${x}; \
samtools sort \
-@ 4 \
${x} \
-o $(basename ${x} .bam)_sorted.bam; done
