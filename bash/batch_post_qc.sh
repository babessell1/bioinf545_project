#!/bin/bash

#SBATCH --job-name=qc_batch
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/qc.log

module load Bioinformatics
module load samtools/1.13-fwwss5n
module load picard-tools/2.8.1

# download qorts if dne
if [ ! -f ../software/QoRTs-1.3.6/QoRTs.jar ]; then
    wget https://github.com/hartleys/QoRTs/archive/refs/tags/v1.3.6.tar.gz \
        --directory-prefix="../software/"
    tar -xf ../software/v1.3.6.tar.gz
fi
# for each species, srr combination run qorts post qc script
for spec in cow sheep; do
    while read srr; do
        bash ./post_qc.sh $srr $spec
    done <../data/srr.txt
done
