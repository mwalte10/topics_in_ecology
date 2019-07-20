#!/bin/csh

#$ -N val_
#$ -t 1-5

module load bio/R/3.4.0

Rscript simp_model_7.7_model_val.R $SGE_TASK_ID