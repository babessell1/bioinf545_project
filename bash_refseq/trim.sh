#!/bin/bash

#SBATCH --job-name=trim
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/trim.log

module load Bioinformatics
module load fastqc/0.11.9-p6ckgle
module load trimgalore/0.6.7-ztb2tpz

mkdir -p ../data/trimmed_reads

# for each srr, trim sample
while read srr; do
    trim_galore \
        --paired \
        --gzip \
        --fastqc \
        --output_dir ../data/trimmed_reads/ \
        "../data/untrimmed_reads/${srr}_pass_1.fastq.gz" \
        "../data/untrimmed_reads/${srr}_pass_2.fastq.gz"
done <../data/srr.txt