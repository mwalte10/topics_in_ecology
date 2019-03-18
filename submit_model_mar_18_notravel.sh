#!/bin/csh

#$ -N mar_18_test_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_mar_18_test.R $SGE_TASK_ID