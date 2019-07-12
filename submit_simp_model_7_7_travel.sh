#!/bin/csh

#$ -N simp_
#$ -t 1-20

module load bio/R/3.4.0

Rscript simp_model_7_7_travel.R $SGE_TASK_ID