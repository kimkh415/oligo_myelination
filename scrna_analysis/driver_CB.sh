#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=4g
#$ -l h_rt=72:00:00
#$ -l os=RedHat7
#$ -pe smp 16
#$ -binding linear:16
#$ -R y

#$ -t 1-6

source /broad/software/scripts/useuse
export PATH=/stanley/levin_dr/kwanho/anaconda3/bin:$PATH:$HOME/bin
use UGER
use .r-3.6.0
use .samtools-1.8
use .hdf5-1.8.16
use .boost-1.70.0
use .java-jdk-1.8.0_181-x86-64
use .openssl-1.0.2g

source ~/kwanho/anaconda3/etc/profile.d/conda.sh
conda activate CellBender

SEEDFILE=data_list_CB.txt ##File with parameters

#Here, $SGE_TASK_ID is the line number
raw_h5=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}') ##gets the first column of the current line
outfile=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}') ##gets the second column of the current line

echo "executing command:"
echo "cellbender remove-background --input ${raw_h5} --output ${outfile} --epochs 150 --low-count-threshold 5"

echo BEGIN
cellbender remove-background --input ${raw_h5} --output ${outfile} --epochs 150 --low-count-threshold 5

echo DONE
