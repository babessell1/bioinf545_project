#!/bin/bash

#SBATCH --job-name=align_sheep
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/align_sheep.log    
       
module load Bioinformatics
module load star/2.7.6a-aekjdpr

mkdir -p ../data/star_sheep

while read srr; do
    STAR \
        --runThreadN 4 \
        --readFilesCommand "gunzip" \
        -c \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 5 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --sjdbGTFfile ../data/genomes/sheep.gtf \
        --readFilesIn \
            "../data/trimmed_reads/${srr}_pass_1_val_1.fq.gz" \
            "../data/trimmed_reads/${srr}_pass_2_val_2.fq.gz" \
        --genomeDir ../data/genomes/sheep_star/ \
        --outFileNamePrefix "../data/star_sheep/${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
done <../data/srr.txt