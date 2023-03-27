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
module load samtools/1.13-fwwss5n
module load picard-tools/2.8.1

# download qorts if dne
if [ ! -f ../software/QoRTs-1.3.6/QoRTs.jar ]; then
    wget https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz \
        --directory-prefix="../software/"
    tar -xf ../software/v1.3.6.tar.gz
fi

mkdir -p ../data/trimmed_bam
mkdir -p ../reports/complexity
mkdir -p ../data/sorted_align

while read srr; do
    # convert fq to sam
    picard-tools FastqToSam \
        F1="../data/trimmed_reads/${srr}_pass_1_val_1.fq.gz" \
        F2="../data/trimmed_reads/${srr}_pass_2_val_2.fq.gz" \
        O="../data/trimmed_bam/${srr}_sorted_reads.bam" \
        SM="${srr}"

    # use picard to estimate library complexity from sorted bams (from fq)
    picard-tools EstimateLibraryComplexity \
        I="../data/trimmed_bam/${srr}_sorted_reads.bam" \
        O="../reports/complexity/${srr}_complexity.txt"

    # sort alignment bam files for qorts
    samtools sort \
        ****unsorted.bam sorted \
        -o "../data/sorted_align/${srr}_aligned.bam"

    # run qorts
    mkdir -p "../reports/qorts_sheep/${srr}"
    mkdir -p "../reports/qorts_cow/${srr}"
    java-Xmx4G-jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
        --generatePlots \
        --genomeFA ../data/genomes/sheep.fna \
        --rawfastq \
            ../data/trimmed_reads/${srr}_pass_1_val_1.fastq.gz.gz, ../data/trimmed_reads/${srr}_pass_2_val_2.fastq.gz \
        ****mybamfile.bam \
        ../data/genomes/sheep.gff \
        "../reports/qorts_sheep/${srr}"

    java-Xmx4G-jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
        --generatePlots \
        --genomeFA ../data/genomes/cow.fna \
        --rawfastq \
            ../data/trimmed_reads/${srr}_pass_1_val_1.fastq.gz.gz, ../data/trimmed_reads/${srr}_pass_2_val_2.fastq.gz \
        ****mybamfile.bam \
        ../data/genomes/cow.gff \
        "../reports/qorts_cow/${srr}"

done <../data/srr.txt


