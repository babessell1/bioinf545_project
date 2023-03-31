#!/bin/bash
sbatch <<EOF
#!/bin/bash

#SBATCH --output="count_${2}_${1}"
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output="../logs/count_${2}_${1}.log"

module load Bioinformatics
module load samtools/1.13-fwwss5n
module load picard-tools/2.8.1

mkdir -p "../data/${2}_gene_cnt"
mkdir -p "../data/${2}_aligned_rmdups"
mkdir -p "../data/${2}_gene_cnt_rmdup"

# index star outputed bam file if not exist
if [ ! -f "../data/star_${2}/${2}_${1}_Aligned.sortedByCoord.out.bam.bai" ]; then
    samtools index \
        -b \
        "../data/star_${2}/${2}_${1}_Aligned.sortedByCoord.out.bam"
fi

# count from star outputed bam file if not exist
htseq-count \
    -f=bam \
    -s=no \
    -r=pos \
    -i=gene_name \
    "../data/star_${2}/${2}_${1}_Aligned.sortedByCoord.out.bam" \
    "../data/genomes/${2}_qortsfix.gtf" \
    > "../data/${2}_gene_cnt/${2}_${1}_gene_cnt.htseq.out"

# mark and remove duplicates in star outputted bam file
mkdir -p "../reports/${2}_dups/"
java -jar ../software/picard-3.0.0/picard.jar MarkDuplicates \
    I="../data/star_${2}/${2}_${1}_Aligned.sortedByCoord.out.bam" \
    O="../data/${2}_aligned_rmdups/${2}_${1}_aligned_rmdup.bam" \
    M="../reports/${2}_dups/${2}_${1}_dups.bam" \
    REMOVE_DUPLICATES=true

# index bam file with removed duplicates
if [ ! -f "../data/${2}_aligned_rmdups/${2}_${1}_aligned_rmdup.bam.bai" ]; then
    samtools index \
        -b \
        "../data/${2}_aligned_rmdups/${2}_${1}_aligned_rmdup.bam"
fi

# count from bam file with removed duplicates
htseq-count \
    -f=bam \
    -s=no \
    -r=pos \
    -i=gene_name \
    "../data/${2}_aligned_rmdups/${2}_${1}_aligned_rmdup.bam" \
    "../data/genomes/${2}_qortsfix.gtf" \
    > "../data/${2}_gene_cnt/${2}_${1}_gene_cnt_rmdup.htseq.out"
EOF