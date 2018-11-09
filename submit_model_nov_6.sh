#!/bin/csh

#$ -N test_new_
#$ -t 1-220

module load bio/R/3.3.1-gcc

Rscript model_nov_6.R $SGE_TASK_ID