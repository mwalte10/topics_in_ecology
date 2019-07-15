#!/bin/csh

#$ -N simp_
#$ -t 1-2

module load bio/R/3.4.0

Rscript simp_model_7_7_skewed_h.R $SGE_TASK_ID