#!/bin/csh

#$ -N test_
#$ -t 1-20

module load bio/R/3.3.1-gcc

Rscript model_y_init.R $SGE_TASK_ID