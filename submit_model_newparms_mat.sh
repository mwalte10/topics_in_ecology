#!/bin/csh

#$ -N test_new_
#$ -t 1-4

module load bio/R/3.3.1-gcc

Rscript model_newparms_mat.R $SGE_TASK_ID