#!/bin/csh

#$ -N test_new_
#$ -t 1-6

module load bio/R/3.3.1-gcc

Rscript model_nov_6_temp.R $SGE_TASK_ID