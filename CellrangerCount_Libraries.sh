#!/bin/bash

#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH --array=0-4
#SBATCH --mem 64G

module load CellRanger/6.0.0

# Go to the HTO (libraries) FASTQ file directory, each sample has it's own folder with it's own FASTQ files
cd GSE211664/geo_submission_hSAP/scRNAseq/HTO/sample${SLURM_ARRAY_TASK_ID}/

id="sample${SLURM_ARRAY_TASK_ID}"
libraries="libraries.csv"
transcriptome="refdata-gex-GRCh38-2020-A"
featRef="feature_ref.csv"
expCells=25000

cellranger count --id=$id --libraries=$libraries --transcriptome=$transcriptome --feature-ref=$featRef --expect-cells=$expCells --localcores=18
