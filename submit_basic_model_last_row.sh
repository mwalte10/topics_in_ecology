#!/bin/csh

#$ -N basic_
#$ -t 1

module load bio/R/3.4.0

Rscript basic_model_last_row.R $SGE_TASK_ID