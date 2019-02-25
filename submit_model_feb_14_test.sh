#!/bin/csh

#$ -N feb_14_test_
#$ -t 1-465

module load bio/R/3.4.0

Rscript model_feb_14_test.R $SGE_TASK_ID