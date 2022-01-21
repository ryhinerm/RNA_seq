#!/usr/bin/env bash

#SBATCH --job-name=prefetch
#SBATCH --partition=pall
#SBATCH --mail-user=marc.ryhiner@students.unibe.ch
#SBATCH --mail-type=begin,fail,end

# Load modules
module add UHTS/Analysis/sratoolkit/2.10.7

# Create and set directory, create empty text file
mkdir ../ex1_download
cd ../ex1_download
touch samples.txt

# Write samples, which are to be downloaded, to text file
echo "SRR9596295.sra" > samples.txt
echo "SRR9596296.sra" >> samples.txt
echo "SRR9596300.sra" >> samples.txt
echo "SRR9596310.sra" >> samples.txt
echo "SRR9596303.sra" >> samples.txt
echo "SRR9596304.sra" >> samples.txt

# Download .sra files
prefetch --option-file samples.txt

# Move .sra files one direcory layer upwards
mv SRR*/*.sra .
while read line
	do rm -r ${line::10}
	done < samples.txt
rm samples.txt

# Covert .sra to .fastq files
fastq-dump --gzip SRR*.sra

# Rename files
mv SRR9596295.fastq.gz Biever_Somata_Poly_1.fastq.gz
mv SRR9596296.fastq.gz Biever_Somata_Poly_2.fastq.gz
mv SRR9596300.fastq.gz Biever_Somata_Poly_3.fastq.gz
mv SRR9596310.fastq.gz Biever_Neuropil_Poly_1.fastq.gz
mv SRR9596303.fastq.gz Biever_Neuropil_Poly_2.fastq.gz
mv SRR9596304.fastq.gz Biever_Neuropil_Poly_3.fastq.gz
