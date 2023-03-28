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

# download qorts if dne
if [ ! -f ../software/QoRTs-1.3.6/QoRTs.jar ]; then
    wget https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz \
        --directory-prefix="../software/"
    tar -xf ../software/v1.3.6.tar.gz
fi

while read srr; do
    mkdir -p ../data/sorted_sheep_align

    # sort alignment bam files
    samtools sort "../data/star_sheep/${srr}_Aligned.sortedByCoord.out.bam" -o "../data/sorted_sheep_align/${srr}_sorted_aligned_sheep.bam"

    # use picard to estimate library complexity from sorted bams (from fq)
    picard-tools EstimateLibraryComplexity \
        I="../data/sorted_sheep_align/${srr}_sorted_aligned_sheep.bam" \
        O="../reports/complexity/${srr}_sheep_complexity.txt"

    # qorts qc reports
    mkdir -p "../reports/qorts_sheep/${srr}"
    java -Xmx4G -jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
        --generatePlots \
        --genomeFA ../data/genomes/sheep.fna \
        --rawfastq \
            "../data/trimmed_reads/${srr}_pass_1.fastq.gz.gz","../data/trimmed_reads/${srr}_pass_2.fastq.gz" \
        "../data/sorted_sheep_align/${srr}_sorted_aligned_sheep.bam" \
        ../data/genomes/sheep.gff \
        "../reports/qorts_sheep/${srr}/"
done <../data/srr.txt
