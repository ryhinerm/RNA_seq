#!/bin/bash

#SBATCH --job-name=map_non-mRNA
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=02:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Create and set directory
mkdir ../ex4_map_non-mRNA
cd ../ex4_map_non-mRNA

# Load module
module add UHTS/Aligner/bowtie/1.2.0

# Unzip trimmed reads
gunzip ../ex3_clip_trim/*_tr.fastq.gz

# Map to non-mRNA
for x in $(ls -d ../ex3_clip_trim/*_tr.fastq); \
do echo ${x}; \
bowtie \
-S \
-t \
-p 4 \
../ex3_prepare_indices/Rnor_r_sn_sno_t_RNA/Rnor_r_sn_sno_t_RNA ${x} \
--un $(basename ${x} .fastq.gz)_no_r_sno_sn_t_RNA.fastq \
2> $(basename ${x} .fastq.gz)_no_r_sno_sn_t_RNA_log.txt > /dev/null; done
