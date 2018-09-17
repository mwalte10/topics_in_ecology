#!/bin/csh

#$ -N test_
#$ -t 1-220

module load bio/R/3.3.1-gcc

Rscript model_parms_mat_high.R $SGE_TASK_ID