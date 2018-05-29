#!/bin/csh

#$ -N model_test_
#$ -t 10-16

module load bio/R/3.3.1-gcc

Rscript model_fixed.R $SGE_TASK_ID
