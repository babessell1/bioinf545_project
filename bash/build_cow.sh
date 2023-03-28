#!/bin/bash

#SBATCH --job-name=build_cow
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/build_cow.log    
       
module load Bioinformatics
module load star/2.7.6a-aekjdpr

# convert gff to gtf (have to download gff_read first)
if [ ! -f ../data/genomes/cow.gtf ]; then
    ../software/gffread-0.12.7.Linux_x86_64/gffread \
        -T ../data/genomes/cow.gff \
        -o ../data/genomes/cow.gtf
fi

STAR --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ../data/genomes/cow_star/ \
    --genomeFastaFiles ../data/genomes/cow.fna \
    --sjdbGTFfile ../data/genomes/cow.gtf \
    --genomeSAindexNbases 13 \
    genomeSAsparseD 3 \
    --sjdbOverhang 100