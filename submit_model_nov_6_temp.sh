#!/bin/csh

#$ -N test_temp_
#$ -t 1

module load bio/R/3.3.1-gcc

Rscript model_nov_6.R $SGE_TASK_ID