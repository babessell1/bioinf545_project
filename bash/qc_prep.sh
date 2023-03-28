#!/bin/bash

#SBATCH --job-name=qc
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/qc.log

module load Bioinformatics
module load samtools/1.13-fwwss5n
module load picard-tools/2.8.1

mkdir -p ../data/read_bam
mkdir -p ../reports/complexity
mkdir -p ../data/sorted_sheep_align

while read srr; do
    # convert fq to sam
    picard-tools FastqToSam \
        F1="../data/untrimmed_reads/${srr}_pass_1.fastq.gz" \
        O="../data/read_bam/${srr}_sorted_reads_1.bam" \
        SM="${srr}"

    # convert fq to sam
    picard-tools FastqToSam \
        F1="../data/untrimmed_reads/${srr}_pass_2.fastq.gz" \
        O="../data/read_bam/${srr}_sorted_reads_2.bam" \
        SM="${srr}"

    # use picard to estimate library complexity from sorted bams (from fq)
    picard-tools EstimateLibraryComplexity \
        I="../data/read_bam/${srr}_sorted_reads_1.bam" \
        O="../reports/complexity/${srr}_1_complexity.txt"

    # use picard to estimate library complexity from sorted bams (from fq)
    picard-tools EstimateLibraryComplexity \
        I="../data/read_bam/${srr}_sorted_reads_2.bam" \
        O="../reports/complexity/${srr}_2_complexity.txt"

    # sort alignment bam files for qorts
    samtools sort \ 
        "../data/genomes/sheep_star/${srr}_Aligned.sortedByCoord.out.bam" \
        -o "../data/sorted_sheep_align/${srr}_aligned.bam"

    # run qorts
    mkdir -p "../reports/qorts_sheep/${srr}"
done <../data/srr.txt


