#!/bin/bash

#SBATCH --job-name=download_srr
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/download.log

module load Bioinformatics
module load sratoolkit/2.10.9-udmejx7

mkdir -p ../data/genomes
mkdir -p ../data/untrimmed_reads

#get genomes
# sheep
curl -OJX \
GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000298735.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT&filename=GCF_000298735.2.zip" -H "Accept: application/zip"

unzip ./*.zip
mv ./ncbi_dataset/data/GCF_000298735.2/*.fna ../data/genomes/sheep.fna
mv ./ncbi_dataset/data/GCF_000298735.2/*.gff ../data/genomes/sheep.gff
rm -r ./ncbi_dataset
rm ./*.zip
rm ./*.md

gzip ./*.fna
gzip ./*.gff

# cow
curl -OJX \
GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000003055.6/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT&filename=GCF_000003055.6.zip" -H "Accept: application/zip"

unzip ./*.zip
mv ./ncbi_dataset/data/GCF_000003055.6/*.fna ../data/genomes/cow.fna
mv ./ncbi_dataset/data/GCF_000003055.6/*.gff ../data/genomes/cow.gff
rm -r ./ncbi_dataset
rm ./*.zip
rm ./*.md

gzip ./*.fna
gzip ./*.gff

# get rnaseq
while read srr; do
    fastq-dump \
        --gzip \
        --readids \
        --read-filter pass \
        --dumpbase \
        --split-files \
        --clip \
        --outdir ../data/untrimmed_reads \
        $srr
done <../data/srr.txt

