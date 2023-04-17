#!/bin/bash
sbatch <<EOF
#!/bin/bash

#SBATCH --output="qc_${2}_${1}"
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/qc_${2}_${1}.log

module load Bioinformatics
module load samtools/1.13-fwwss5n
module load picard-tools/2.8.1
module load R

mkdir -p "../data/sorted_${2}"

# sort alignment bam files
if [ ! -f "../data/sorted_${2}/${2}_${1}_sorted.bam" ]; then
    samtools sort "../data/star_${2}/${2}_${1}_Aligned.sortedByCoord.out.bam" -o "../data/sorted_${2}/${2}_${1}_sorted.bam"
fi

# qorts qc reports
mkdir -p "../reports/qorts_${2}/${2}_${1}"
java -Xmx4G -jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
    --generatePlots \
    --genomeFA "../data/genomes/${2}.fna" \
    --rawfastq \
        "../data/untrimmed_reads/${1}_pass_1.fastq.gz","../data/untrimmed_reads/${1}_pass_2.fastq.gz" \
    "../data/sorted_${2}/${2}_${1}_sorted.bam" \
    "../data/genomes/${2}_qortsfix.gtf" \
    "../reports/qorts_${2}/${2}_${1}/"

mkdir -p "../data/${2}_sorted_aligned_rmdups"
if [ ! -f "../data/${2}_sorted_aligned_rmdups/${2}_${1}_sorted_aligned_rmup.bam" ]; then
    samtools sort "../data/${2}_aligned_rmdups/${2}_${1}_aligned_rmdup.bam" -o "../data/${2}_sorted_aligned_rmdups/${2}_${1}_sorted_aligned_rmup.bam"
fi

mkdir -p "../reports/qorts_${2}/${2}_${1}_rmdups"
java -Xmx4G -jar ../software/QoRTs-1.3.6/QoRTs.jar QC \
    --generatePlots \
    --genomeFA "../data/genomes/${2}.fna" \
    --rawfastq \
        "../data/untrimmed_reads/${1}_pass_1.fastq.gz","../data/untrimmed_reads/${1}_pass_2.fastq.gz" \
    "../data/${2}_sorted_aligned_rmdups/${2}_${1}_sorted_aligned_rmup.bam" \
    "../data/genomes/${2}_qortsfix.gtf" \
    "../reports/qorts_${2}/${2}_${1}_rmdups/"
EOF