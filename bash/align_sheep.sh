#!/bin/bash

#SBATCH --job-name=build_sheep
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/build_sheep.log    
       
module load Bioinformatics
module load star/2.7.6a-aekjdpr

mkdir -p ../data/star_output

while read srr; do
    STAR \
        --runThreadN 4 \
        --readFilesCommand "gunzip" \
        -c \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --readFilesIn \
            "../data/untrimmed_reads/${srr}_pass_1.fastq.gz" \
            "../data/untrimmed_reads/${srr}_pass_2.fastq.gz" \
        --genomeDir ../data/genomes/sheep_star/ \
        --outFileNamePrefix "../data/genomes/sheep_star/${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
done <../data/srr.txt