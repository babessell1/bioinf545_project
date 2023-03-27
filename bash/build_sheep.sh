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

STAR --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ../data/genomes/sheep_star/ \
    --genomeFastaFiles ../data/genomes/sheep.fna \
    --sjdbGTFfile ../data/genomes/sheep.gff \
    --genomeSAindexNbases 13 \
    genomeSAsparseD 3 \
    --sjdbOverhang 100