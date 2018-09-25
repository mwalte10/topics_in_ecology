#!/bin/csh

#$ -N test_old_
#$ -t 1-220

module load bio/R/3.3.1-gcc

Rscript model_oldparms_mat.R $SGE_TASK_ID