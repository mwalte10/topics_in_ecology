#!/bin/csh

#$ -N feb_14_notest_
#$ -t 1-20

module load bio/R/3.4.0

Rscript model_feb_14_notest.R $SGE_TASK_ID