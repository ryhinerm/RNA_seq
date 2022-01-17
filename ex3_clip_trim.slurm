#!/bin/bash

#SBATCH --job-name=clip_trim
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=08:00:00
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end


# Load module
module add UHTS/Quality_control/cutadapt/2.5

# Create and set directory
mkdir ../ex3_clip_trim
cd ../ex3_clip_trim

# Clip 3' adapter and 2 nt from the 5' end
for x in $(ls -d ../ex1_download/*.fastq.gz); do echo ${x}; \
cutadapt \
-j 12 \
-a AGATCGGAAGAGCACACGTCTGAA \
-q 25 \
--cut 2 \
--minimum-length 22 \
--discard-untrimmed \
--overlap 3 \
-e 0.2 \
-o $(basename ${x} .fastq.gz)_clpd.fastq.gz \
${x} 1> $(basename ${x} .fastq.gz)_clpd_log.txt; done


# Trim 10 nt (5 nt randomized bases + 5 nt Barcode) from the 3' end
for x in $(ls -d *_clpd.fastq.gz); do echo ${x}; \
cutadapt \
-j 12 \
-q 25 \
--cut -10 \
--minimum-length 22 \
-o $(basename ${x} .fastq.gz)_tr.fastq.gz \
${x} 1> $(basename ${x} .fastq.gz)_tr_log.txt; done
