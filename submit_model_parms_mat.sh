#!/bin/csh

#$ -N test_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript model_parms_mat.R $SGE_TASK_ID