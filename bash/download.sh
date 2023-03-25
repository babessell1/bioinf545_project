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
curl -OJX \
--output-dir data/genomes/ \
GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000298735.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000298735.2.zip" -H "Accept: application/zip"

curl -OJX \
--output-dir data/genomes/ \
GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000003055.6/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT&filename=GCF_000003055.6.zip" -H "Accept: application/zip"

# get rnaseq
while read srr; do
    fasterq-dump \
    --gzip \
    --readids \
    --read-filter pass \
    --dumpbase \
    --split-files \
    --clip \
    --outdir ../data/untrimmed_reads \
    $srr
done <srr.txt

