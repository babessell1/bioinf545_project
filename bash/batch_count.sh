#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=5:00:00
#SBATCH --account=bioinf545w23_class
#SBATCH --partition=standard
#SBATCH --output=../logs/batch_count.log

# for each species, srr combo, run count genes alignments script
for spec in cow sheep; do
    while read srr; do
        bash ./count.sh $srr $spec
    done <../data/srr.txt
done