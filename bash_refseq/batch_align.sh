#!/bin/bash

#SBATCH --job-name=batch_align
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/batch_align.log    

# for each species, srr combination, run alignment script
for spec in cow sheep; do
    while read srr; do
        bash ./align.sh $srr $spec
    done <../data/srr.txt
done