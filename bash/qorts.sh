#!/bin/bash

#SBATCH --job-name=qorts
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/qorts.log

# download qorts if dne
if [ ! -f ../software/QoRTs-1.3.6/QoRTs.jar ]; then
    wget https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz \
        --directory-prefix="../software/"
    tar -xf ../software/v1.3.6.tar.gz
fi

while read srr; do
    java-Xmx4G-jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
        --generatePlots \
        --genomeFA ../data/genomes/sheep.fna \
        --rawfastq \
            "../data/trimmed_reads/${srr}_pass_1.fastq.gz.gz", "../data/trimmed_reads/${srr}_pass_2.fastq.gz" \
        "../data/sorted_sheep_align/${srr}_aligned.bam" \
        ../data/genomes/sheep.gff \
        "../reports/qorts_sheep/${srr}/"
done <../data/srr.txt
