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

 convert gff to gtf if not already done
../software/gffread-0.12.7.Linux_x86_64/gffread \
    -T ../data/genomes/sheep.gff \
    -o ../data/genomes/sheep.gtf
grep 'gene_id' ../data/genomes/sheep.gtf > ../data/genomes/sheep_qortsfix.gtf


# index genome
STAR --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ../data/genomes/sheep_star/ \
    --genomeFastaFiles ../data/genomes/sheep_ens.fna \
    --sjdbGTFfile ../data/genomes/sheep_ens.gtf \
    --genomeSAindexNbases 13 \
    --genomeSAsparseD 1 \
    --sjdbOverhang 99