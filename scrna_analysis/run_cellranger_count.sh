#!/bin/sh

#SBATCH --job-name=count
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --mem=80G

source ~/.bashrc
#conda activate CellBender
module load cellranger/3.0.1-fasrc01


##Here, $SLURM_ARRAY_TASK_ID is the line number
SEEDFILE=count.data.txt ##File with parameters

sample=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE | awk '{print $1}') ##gets the first column of the current line
fastq_dir=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE | awk '{print $2}') ##gets the second column of the current line

##Here, $SLURM_ARRAY_TASK_ID is the line number

echo "Processing ${sample}"
echo "FASTQ directory: ${fastq_dir}"

cellranger count --id=${sample} --fastqs=${fastq_dir} --transcriptome=/n/boslfs02/LABS/arlotta_lab/Lab/Kwanho/ref/mm10_3.0.0/mm10 --sample=${sample} --localmem=64


echo "DONE"

