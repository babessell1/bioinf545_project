#!/bin/bash

sbatch <<EOF 
#!/bin/bash
#SBATCH --job-name="align_${2}_${1}"
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output="../logs/align_${2}_${1}.log"    
       
module load Bioinformatics
module load star/2.7.6a-aekjdpr

mkdir -p "../data/star_${2}"

# align species=2, srr=1 combination
STAR \
    --runThreadN 4 \
    --readFilesCommand "gunzip" -c \
    --outFilterMultimapNmax 1 \
    --alignSJoverhangMin 5 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --sjdbGTFfile "../data/genomes/${2}_qortsfix.gtf" \
    --readFilesIn \
        "../data/trimmed_reads/${1}_pass_1_val_1.fq.gz" \
        "../data/trimmed_reads/${1}_pass_2_val_2.fq.gz" \
    --genomeDir "../data/genomes/${2}_star/" \
    --outFileNamePrefix "../data/star_${2}/${2}_${1}_" \
    --outSAMtype BAM SortedByCoordinate 
EOF