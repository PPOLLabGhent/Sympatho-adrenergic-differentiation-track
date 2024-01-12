#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH --array=0-4
#SBATCH --mem 64G


module load CellRanger/6.0.0

# Go to the RNA FASTQ file directory, each sample has it's own folder with it's own FASTQ files
cd GSE211664/geo_submission_hSAP/scRNAseq/RNA/sample${SLURM_ARRAY_TASK_ID}/

id="sample${SLURM_ARRAY_TASK_ID}"
fastqs="GSE211664/geo_submission_hSAP/scRNAseq/RNA/sample${SLURM_ARRAY_TASK_ID}/FASTQ"
transcriptome="refdata-gex-GRCh38-2020-A"
sample="$(cat samName.txt)"
expCells=25000

cellranger count --id=$id --fastqs=$fastqs --transcriptome=$transcriptome --sample=$sample --expect-cells=$expCells --localcores=36
